"""
This will set up a ddPCR, generic PCR, or Illumina Dual Indexing PCR using the Opentrons OT-2 robot.
Done in Gaorav Gupta's lab.
Author: Dennis Simpson
Address:    University of North Carolina at Chapel Hill
            Lineberger Comprehensive Cancer Center
            Chapel Hill, NC 27599
email: dennis@email.unc.edu
Copyright:  2025
"""

import datetime
import os
import csv
import platform
from csv import excel

import serial
import time
from types import SimpleNamespace
from contextlib import suppress
from collections import defaultdict
from opentrons import protocol_api
import math
import Tool_Box

# metadata
metadata = {
    'protocolName': 'PCR v4.1.3',
    'author': 'Dennis Simpson <dennis@email.unc.edu>',
    'description': 'Setup a ddPCR, Generic PCR, or Dual Indexing PCR'
    }

# requirements
requirements = {"robotType": "OT-2", "apiLevel": "2.20"}


def add_parameters(parameters: protocol_api.Parameters):

    """
    Parse the TSV file and fill in some parameter information.  This is duplicated from Utilities.  I don't know
    another method to pass the information when on the robot.
    @param parameters:
    """
    # TSV file location on OT-2
    tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)

    # If not on the OT-2, get temp TSV file location on Windows Computers for simulation
    if not os.path.isfile(tsv_file_path):
        tsv_file_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())

    line_num = 0
    options_dictionary = defaultdict(str)
    sample_dictionary = defaultdict(list)
    index_file = list(csv.reader(open(tsv_file_path), delimiter='\t'))

    for line in index_file:
        if line_num == 0:
            options_dictionary["Version"] = line[1]
            options_dictionary["Template"] = line[0].strip("#")
        line_num += 1
        col_count = len(line)
        tmp_line = []
        sample_key = ""
        if col_count > 0 and "#" not in line[0] and len(line[0].split("#")[0]) > 0:
            # Skip any lines that are blank or comments.
            for i in range(7):
                try:
                    line[i] = line[i].split("#")[0]  # Strip out end of line comments.
                except IndexError:
                    continue

                if i == 0 and "--" in line[0]:
                    key = line[0].strip('--')
                    key_value = line[1]
                    if "Target_" in key or "PositiveControl_" in key:
                        try:
                            key_value = (line[1], line[2], line[3])
                        except IndexError:
                            pass
                    options_dictionary[key] = key_value
                elif "--" not in line[0] and int(line[0]) < 12:
                    sample_key = line[0], line[1]
                    tmp_line.append(line[i])
            if sample_key:
                sample_dictionary[sample_key] = tmp_line

    args = SimpleNamespace(**options_dictionary)

    # We have limited space for the run_label.  To make sure the label is unique, I use a unix timestamp for the run_date.
    run_date = datetime.datetime.today().strftime("%f")
    if "Illumina_Dual_Indexing" in args.Template:
        run_type = "Dual_Indexing"
    else:
        run_type = args.Template.split(" ")[1]

    # This is used by the Opentrons app to make each run unique.
    parameters.add_str(
        variable_name="run_label",
        display_name="{} {} {}".format(run_type, args.User, run_date),
        choices=[
            {"display_name": "Run Label", "value": "{} {} {}".format(run_type, args.User, run_date)},],
        default="{} {} {}".format(run_type, args.User, run_date),
    )

    parameters.add_str(
        variable_name="left_pipette",
        display_name="Left Pipette".format(run_type, args.User, run_date),
        choices=[
            {"display_name": "P20 Single Gen2", "value": "p20_single_gen2"},
            {"display_name": "P300 Single Gen2", "value": "p300_single_gen2"},],
        default="p300_single_gen2",
    )

    parameters.add_str(
        variable_name="right_pipette",
        display_name="Right Pipette".format(run_type, args.User, run_date),
        choices=[
            {"display_name": "P20 Single Gen2", "value": "p20_single_gen2"},
            {"display_name": "P300 Single Gen2", "value": "p300_single_gen2"},],
        default="p20_single_gen2",
    )

    """   
    parameters.add_csv_file(
        variable_name="dragons_run",
        display_name="Dragon Hunting",
        description="Looking to initialize Opentrons csv commands."
    )
    """


def calculate_volumes(args, sample_concentration, template_in_rxn):
    """
    Calculates volumes for dilution and distribution of sample.
    Returns a list of tuples consisting of
    (uL of sample to dilute, uL of water for dilution), (uL of diluted sample in reaction, uL of water in reaction)

    :param args:
    :param sample_concentration:
    :param template_in_rxn:
    :return:
    """

    max_template_vol = round(float(args.PCR_Volume)-float(args.MasterMixPerRxn), ndigits=1)

    # If at least 2 uL of sample is needed then no dilution is necessary
    if template_in_rxn/sample_concentration >= 2:
        sample_vol = round(template_in_rxn/sample_concentration, ndigits=1)
        return sample_vol, 0, 0, max_template_vol-sample_vol, max_template_vol

    # This will test a series of dilutions up to a 1:200.
    for i in range(50):
        dilution = (i+1)*2
        diluted_dna_conc = sample_concentration/dilution

        # Want to pipette at least 2 uL of diluted sample per well
        if 2 <= template_in_rxn/diluted_dna_conc <= max_template_vol:
            diluted_sample_vol = round(template_in_rxn/diluted_dna_conc, ndigits=1)
            reaction_water_vol = max_template_vol-diluted_sample_vol

            return 1, dilution - 1, diluted_sample_vol, reaction_water_vol, max_template_vol


def sample_processing(args, sample_parameters, target_info_dict, utility):
    sample_data_dict = defaultdict(list)
    target_well_dict = defaultdict(list)
    water_well_dict = defaultdict(float)
    plate_layout_by_column, layout_data = utility.plate_layout(args.PCR_PlateSlot)
    used_wells = []
    dest_well_count = 0
    target_list = []

    for sample_key in sample_parameters:
        sample_name = sample_parameters[sample_key][2]
        sample_concentration = float(sample_parameters[sample_key][3])

        if "Illumina_Dual_Indexing" in args.Template:
            # Extract indices.  Always going to be two values.
            # sample_targets = sample_parameters[sample_key][4].split("+")
            sample_targets = [sample_parameters[sample_key][4]]
            replicates = 1
        else:
            sample_targets = sample_parameters[sample_key][4].split(",")
            replicates = int(sample_parameters[sample_key][5])

        # Generic PCR allows different amounts of DNA in each reaction.
        if "Generic PCR" in args.Template:
            template_in_rxn = float(sample_parameters[sample_key][6])
        else:
            template_in_rxn = float(args.DNA_in_Reaction)

        # sample_string = sample_name
        ucv = calculate_volumes(args, sample_concentration, template_in_rxn)
        sample_vol = round(ucv[0], ndigits=1)
        diluent_vol = round(ucv[1], ndigits=1)
        diluted_sample_vol = round(ucv[2], ndigits=1)
        reaction_water_vol = ucv[3]
        max_template_vol = ucv[4]

        sample_wells = []
        for target in sample_targets:
            target_list.append(target)
            if "Illumina_Dual_Indexing" in args.Template:
                target_name = target
            else:
                target_name = target_info_dict[int(target)][1]

            for i in range(replicates):
                well = plate_layout_by_column[dest_well_count]
                try:
                    row = well[0]
                except IndexError:
                    row = well
                column = int(well[1:])-1
                s_volume = diluted_sample_vol

                if diluent_vol == 0:
                    dilution = "Neat"
                    s_volume = sample_vol
                else:
                    dilution = "1:{}".format(int((sample_vol + diluent_vol) / sample_vol))

                    required_vol = round(
                        (diluted_sample_vol * len(sample_targets) * replicates) + (3 * diluted_sample_vol), ndigits=1)
                    dilution_factor = round(sample_vol / (sample_vol + diluent_vol), ndigits=4)
                    final_sample_vol = round(dilution_factor * required_vol, ndigits=1)

                    sample_vol = round(dilution_factor * required_vol, ndigits=1)
                    diluent_vol = round(required_vol - final_sample_vol, ndigits=1)

                layout_data[row][column] = "{}|{}|{}|{}"\
                    .format(sample_name, target_name, dilution, s_volume)

                water_well_dict[well] = reaction_water_vol
                target_well_dict[target].append(well)
                sample_wells.append(well)
                used_wells.append(well)
                dest_well_count += 1

        sample_data_dict[sample_key] = [sample_vol, diluent_vol, diluted_sample_vol, sample_wells]

    # Define our no template control wells for the targets.
    for target in target_well_dict:
        target_list.append(target)

    target_list = list(set(target_list))

    # There is never a no template control for the Illumina Dual Indexing PCR.
    if "Illumina_Dual_Indexing" not in args.Template:
        for target in target_list:
            control_name = "Water"
            target_name = target_info_dict[int(target)][1]
            well = plate_layout_by_column[dest_well_count]
            used_wells.append(well)
            row = well[0]
            column = int(well[1:])-1
            layout_data[row][column] = "{}|{}|NA|{}".format(control_name, target_name, max_template_vol)
            water_well_dict[well] = max_template_vol
            dest_well_count += 1

            target_well_dict[target].append(well)

    return sample_data_dict, water_well_dict, target_well_dict, used_wells, layout_data, max_template_vol


def run(protocol: protocol_api.ProtocolContext):
    protocol.comment(protocol.params.run_label)
    protocol.set_rail_lights(True)

    utility = Utilities(protocol)
    sample_parameters, args = utility.parse_sample_template()
    utility.labware_parsing()

    left_tipracks, right_tipracks = utility.tipracks
    labware, slot_dict = utility.deck_layout

    # Pipettes
    left_pipette = protocol.load_instrument(protocol.params.left_pipette, 'left', tip_racks=left_tipracks)
    right_pipette = protocol.load_instrument(protocol.params.right_pipette, 'right', tip_racks=right_tipracks)

    # Set the location of the first tip in box.
    with suppress(IndexError):
        left_pipette.starting_tip = left_tipracks[0].wells_by_name()[args.LeftPipetteFirstTip.upper()]
    with suppress(IndexError):
        right_pipette.starting_tip = right_tipracks[0].wells_by_name()[args.RightPipetteFirstTip.upper()]

    # Turn off rail lights for actual run.
    if not protocol.is_simulating():
        protocol.set_rail_lights(False)

    if args.UseTemperatureModule == True and not protocol.is_simulating():
        temp_mod = ColdPlateSlimDriver(protocol)
        temp_mod.quick_temp(float(args.Temperature))
        # print("Using Temperature Module: ", temp_mod.get_info())

    target_info_dict = defaultdict(list)

    # Read targeting parameters into dictionary if not running an Indexing PCR.
    if "Illumina_Dual_Indexing" not in args.Template:
        for i in range(10):
            target = getattr(args, "Target_{}".format(i + 1), "")
            if target:
                # target_info_dict[i + 1] = target.split("|")
                target_info_dict[i + 1] = target

            # if len(target[0]) > 1:
            else:
                """
                if not all('' == s or s.isspace() for s in target):
                    target_info_dict[i + 1] = target
                """
                target_info_dict[i + 1] = target

    sample_data_dict, water_well_dict, target_well_dict, used_wells, layout_data, max_template_vol = \
        sample_processing(args, sample_parameters, target_info_dict, utility)

    # This will output a plate layout file.  Only does it during the simulation from our GUI
    if protocol.is_simulating() and platform.system() == "Windows":
        try:
            run_date = datetime.datetime.today().strftime("%a %b %d %H:%M %Y")
            plate_layout_string = \
                "## {} Setup\n## Setup Date:\t{}\n## Template User:\t{}\n" \
                "# Format:\tTemplate | Target | Template Dilution | Template Volume in Reaction\n\n\t"\
                .format(args.Template, run_date, args.User)

            for i in range(12):
                plate_layout_string += "{}\t".format(i+1)

            # I have to import this here because I have been unable to get natsort on the robot.
            import natsort

            for well in natsort.natsorted(layout_data):
                well_string = "\t".join(layout_data[well])
                plate_layout_string += "\n{}\t{}\t".format(well, well_string)

            plate_layout_file = \
                open("C:{0}Users{0}{1}{0}Documents{0}{2}_PlateLayout.tsv".
                     format(os.sep, os.getlogin(), args.Template), 'w')

            plate_layout_file.write(plate_layout_string)
            plate_layout_file.close()

        except ModuleNotFoundError:
            pass

    # Now do the actual dispensing.
    water_aspirated = utility.dispense_water(water_well_dict, left_pipette, right_pipette)
    utility.dispense_reagent_mix(labware, target_well_dict, target_info_dict, left_pipette, right_pipette)

    if "Illumina_Dual_Indexing" in args.Template:
        dispense_indexing_primers(args, protocol, utility, left_pipette, right_pipette, labware, sample_parameters,
                                  sample_data_dict)

    water_aspirated = dispense_samples(args, labware, sample_data_dict, sample_parameters, left_pipette,
                                       right_pipette, water_aspirated, utility, protocol)
    if "ddPCR" in args.Template:
        fill_empty_wells(args, used_wells, water_aspirated, labware, left_pipette, right_pipette, utility)

    # If using Temperature Module, hold PCR plate at set temperature until user removes it and closes program.
    if args.UseTemperatureModule != "FALSE" and not protocol.is_simulating():
        protocol.set_rail_lights(True)
        protocol.comment("Program is complete.  Temperature is holding at {}. Click RESUME to exit.".format(args.Temperature))

        protocol.pause()
        temp_mod.deactivate()
        protocol.set_rail_lights(False)
    else:
        protocol.comment("Program Complete")

    if not protocol.is_simulating():
        os.remove(utility.parameter_file)


def dispense_indexing_primers(args, protocol, utility, left_pipette, right_pipette, labware, sample_parameters,
                              sample_data_dict):
    protocol.comment("\nDispensing Indexing Primers")

    # Extract Index Primer information
    index_primers = \
        ["D501", "D502", "D503", "D504", "D505", "D506", "D507", "D508", "D701", "D702", "D703", "D704", "D705",
         "D706", "D707", "D708", "D709", "D710a", "D711", "D712"]

    primer_wells = {}

    # Build a dictionary of the index primer wells
    for i in range(len(index_primers)):
        # primer_well = getattr(args, "{}".format(index_primers[i]))
        primer_wells[index_primers[i]] = getattr(args, "{}".format(index_primers[i]))

    # Determine primer volumes and dispense them.
    # 6.25 uM = 2 uL per 50 uL
    # 10 uM = 1.25 uL per 50 uL
    # primer_volume = (float(args.PCR_Volume)/50) * 1.25
    primer_volume = (float(args.PCR_Volume) / 50) * 2.0
    selected_pipette = utility.pipette_selection(left_pipette, right_pipette, primer_volume)

    # Step 3: Identify labware for primers and sample destinations
    primer_labware = labware[args.IndexPrimerSlot]
    sample_destination_labware = labware[args.PCR_PlateSlot]

    for sample_key in sample_parameters:
        destination_well = sample_data_dict[sample_key][3][0]
        sample_data = sample_parameters[sample_key]
        sample_well = sample_data[1]
        d500, d700 = sample_data[4].split("+")

        # Dispense D500 primer
        utility.pipette_reagents(selected_pipette, primer_labware[primer_wells[d500]].bottom(float(args.BottomOffset)),
                                 sample_destination_labware[destination_well], primer_volume, NewTip=True, MixReaction=False,
                                 touch=True)

        # Dispense D700 primer
        utility.pipette_reagents(selected_pipette, primer_labware[primer_wells[d700]].bottom(float(args.BottomOffset)),
                                 sample_destination_labware[destination_well], primer_volume, NewTip=True, MixReaction=False,
                                 touch=True)


def fill_empty_wells(args, used_wells, water_aspirated, labware_dict, left_pipette, right_pipette, utility):
    """
    This will fill the remaining wells in a column with water.  Needed to for the droplet generator.
    """

    last_used_well = used_wells[-1]
    row = last_used_well[0]
    column = int(last_used_well.split(row)[1])
    row_list = ["A", "B", "C", "D", "E", "F", "G", "H"]
    row_index = row_list.index(row)
    wells_remaining = len(row_list)-row_index-1

    if wells_remaining > 0:
        sample_destination_labware = labware_dict[args.PCR_PlateSlot]
        reagent_labware = labware_dict[args.ReagentSlot]
        water_res_well_dia = reagent_labware[args.WaterResWell].diameter
        # cone_vol = utility.labware_cone_volume(args.ReagentSlot)
        fill_pipette = \
            utility.pipette_selection(left_pipette, right_pipette, water_aspirated)
        water_tip_height = \
            utility.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia)

        for i in range(wells_remaining):
            blank_well = "{}{}".format(row_list[i+row_index+1], column)

            utility.pipette_reagents(fill_pipette, reagent_labware[args.WaterResWell].bottom(water_tip_height),
                                     sample_destination_labware[blank_well], float(args.PCR_Volume),
                                     NewTip=False, MixReaction=False
                                     )

            water_aspirated = water_aspirated+float(args.PCR_Volume)
            water_tip_height = utility.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia)

        fill_pipette.drop_tip()


def dispense_samples(args, labware_dict, sample_data_dict, sample_parameters, left_pipette, right_pipette,
                     water_aspirated, utility, protocol):
    """
    Dilute and dispense samples
    @param utility:
    @param args:
    @param labware_dict:
    @param sample_data_dict:
    @param sample_parameters:
    @param left_pipette:
    @param right_pipette:
    @param water_aspirated:
    """
    protocol.comment("\nDiluting and Dispensing Samples")
    try:
        dilution_labware = labware_dict[args.DilutionPlateSlot]
    except KeyError:
        dilution_labware = ""

    bottom_offset = float(args.BottomOffset)
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    reagent_labware = labware_dict[args.ReagentSlot]
    water_res_well_dia = labware_dict[args.ReagentSlot][args.WaterResWell].diameter
    water_tip_height = utility.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia)

    # If the user determines no dilutions are required they can leave that slot blank.  I don't like this approach,
    # users could leave the information out and dilutions might still be required.
    if dilution_labware:
        slot = args.DilutionPlateSlot
        # labware_name = slot_dict[args.DilutionPlateSlot]
    else:
        slot = args.PCR_PlateSlot
        # labware_name = slot_dict[args.PCR_PlateSlot]

    dilution_plate_layout, unused_layout = utility.plate_layout(slot)

    dilution_well_index = 0

    for sample_key in sample_parameters:
        sample_source_labware = labware_dict[sample_parameters[sample_key][0]]
        sample_source_well = sample_parameters[sample_key][1]
        sample_dest_wells = sample_data_dict[sample_key][3]
        sample_vol = sample_data_dict[sample_key][0]
        diluent_vol = sample_data_dict[sample_key][1]
        diluted_sample_vol = sample_data_dict[sample_key][2]
        mix_volume = None

        if float(args.PCR_Volume) > 20:
            mix_volume = 17

        # If no dilution is necessary, dispense sample and continue
        if diluted_sample_vol == 0:

            sample_pipette = \
                utility.pipette_selection(left_pipette, right_pipette, sample_vol)

            for well in sample_dest_wells:
                utility.pipette_reagents(sample_pipette, sample_source_labware[sample_source_well],
                                         sample_destination_labware[well], sample_vol,
                                         NewTip=True, MixReaction=True, touch=True, MixVolume=mix_volume
                                         )
        else:
            water_aspirated, dilution_well_index = sample_dilution(args, sample_source_labware, sample_source_well, sample_vol, diluent_vol,
                                              dilution_plate_layout, dilution_well_index, reagent_labware,
                                              water_tip_height, dilution_labware, diluted_sample_vol, water_res_well_dia,
                                              sample_dest_wells, sample_destination_labware, bottom_offset, left_pipette,
                                              right_pipette, water_aspirated, utility)

    utility.drop_any_tips([left_pipette, right_pipette])
    return water_aspirated

def sample_dilution(args, sample_source_labware, sample_source_well, sample_vol, diluent_vol, dilution_plate_layout,
                    dilution_well_index, reagent_labware, water_tip_height, dilution_labware, diluted_sample_vol,
                    water_res_well_dia, sample_dest_wells, sample_destination_labware, bottom_offset, left_pipette, right_pipette,
                    water_aspirated, utility):
    # Adjust volume of diluted sample to make sure there is enough
    # diluted_template_needed = round(diluted_sample_vol*(len(sample_dest_wells)+1.5), ndigits=1)
    # diluted_template_factor = round(diluted_template_needed/(sample_vol+diluent_vol), ndigits=1)
    # adjusted_sample_vol = round((sample_vol * diluted_template_factor), ndigits=1)
    # diluent_vol = round((diluent_vol*diluted_template_factor), ndigits=1)

    # Reset the pipettes for the new volumes
    diluent_pipette = utility.pipette_selection(left_pipette, right_pipette, diluent_vol)
    sample_pipette = utility.pipette_selection(left_pipette, right_pipette, sample_vol)

    # Make dilution, diluent first
    dilution_well = dilution_plate_layout[dilution_well_index]

    utility.pipette_reagents(diluent_pipette, reagent_labware[args.WaterResWell].bottom(water_tip_height),
                             dilution_labware[dilution_well], diluent_vol, NewTip=True, MixReaction=False,
                             touch=True)
    mix_volume = None
    if diluted_sample_vol < 20:
        mix_volume = 18

    utility.pipette_reagents(sample_pipette, sample_source_labware[sample_source_well],
                             dilution_labware[dilution_well], sample_vol, NewTip=True, MixReaction=True,
                             MixVolume=mix_volume)

    water_aspirated += diluent_vol
    dilution_well_index += 1
    water_tip_height = \
        utility.res_tip_height(float(args.WaterResVol) - water_aspirated, water_res_well_dia)

    # Add diluted sample to PCR plate
    for well in sample_dest_wells:
        sample_pipette = \
            utility.pipette_selection(left_pipette, right_pipette, diluted_sample_vol)
        mix_volume = None

        if diluted_sample_vol < 20:
            mix_volume = 18

        utility.pipette_reagents(sample_pipette, dilution_labware[dilution_well].bottom(bottom_offset),
                                 sample_destination_labware[well], diluted_sample_vol, NewTip=True,
                                 MixReaction=True, MixVolume=mix_volume
                                 )
    return water_aspirated, dilution_well_index

class ColdPlateSlimDriver:
    """
    (С) Parhelia Biosciences Corporation 2024-2025
    Class to control their temperature module.
    @todo Need to find out if the Opentrons built in module will work.
    """
    def __init__(self, protocol_context, temp_mode_number=0):
        self.serial_number = "29533"
        self.device_name = "/dev/ttyUSB" + str(temp_mode_number)
        self.baudrate = 9600
        self.bytesize = serial.EIGHTBITS
        self.parity = serial.PARITY_NONE
        self.stopbits = serial.STOPBITS_ONE
        self.read_timeout = 2
        self.write_timeout = 2
        self.height = 45
        self.target_temp = 20
        self.protocol = protocol_context

        # check context, skip if simulating Linux
        if protocol_context.is_simulating():
            # self.protocol.comment("Simulation detected. Initializing Temperature Module in the dummy mode\n")
            self.serial_object = None
        else:
            # self.protocol.comment("Execution mode detected.  Initializing Temperature Module in the run mode\n")
            self.serial_object = serial.Serial(
                port=self.device_name,
                baudrate=self.baudrate,
                bytesize=self.bytesize,
                parity=self.parity,
                stopbits=self.stopbits,
                timeout=self.read_timeout,
                write_timeout=self.write_timeout,
            )

    def _reset_buffers(self):
        """
        Worker function
        """
        if self.serial_object is None:
            return
        self.serial_object.reset_input_buffer()
        self.serial_object.reset_output_buffer()

    def _read_response(self):
        """
        Worker function
        """
        if self.serial_object is None:
            return "Program is simulating"

        output_lines = self.serial_object.readlines()
        output_string = "".join(l.decode("utf-8") for l in output_lines)

        return output_string

    def _send_command(self, my_command):
        """
        Worker function
        @param my_command:
        @return:
        """
        SERIAL_ACK = "\r\n"

        command = my_command
        command += SERIAL_ACK

        if self.serial_object is None:
            print("Simulating: " + my_command)

            return

        self.serial_object.write(command.encode())
        self.serial_object.flush()

        return self._read_response()

    def get_info(self):
        if self.serial_object is None:
            return "Simulating or no temperature module detected."

        return self._send_command("info")

    def get_temp(self):
        """
        Get the module temperature.
        @return:
        """
        if self.serial_object is None:
            return self.target_temp

        actual_temp = round(float(self._send_command("getTempActual")), 1)
        return actual_temp

    def set_temperature(self, my_temp):
        """
        Send temperature command to the module.
        @param my_temp:
        @return:
        """
        if self.serial_object is None and self.protocol.is_simulating():
            self.target_temp = my_temp
            return

        temp = float(my_temp) * 10
        temp = int(temp)
        self._send_command(f"setTempTarget{temp:03}")
        self._send_command("tempOn")

    def deactivate(self):
        """
        Shutdown the temperature module and close the serial connection.
        """
        if self.serial_object is None:
            self.target_temp = 25
        else:
            self._send_command("tempOff")
            self.serial_object.close()

    @staticmethod
    def time_to_reach_sample_temp(delta_temp):
        """
        Estimate the time in minutes required for module to reach target temperature.
        @param delta_temp:
        @return:
        """
        x = delta_temp

        if x > 0:
            time_min = 0.364 + 0.559 * x - 0.0315 * x ** 2 + 1.29E-03 * x ** 3 - 2.46E-05 * x ** 4 + 2.21E-07 * x ** 5 - 7.09E-10 * x ** 6
        else:
            time_min = -0.1 - 0.329 * x - 0.00413 * x ** 2 - 0.0000569 * x ** 3 + 0.0000000223 * x ** 4

        return round(time_min, 2)

    def quick_temp(self, temp_target, overshot=3, undershot=2):
        """
        Set the module temperature and apply a delay, if needed.
        @param temp_target:
        @param overshot:
        @param undershot:
        """
        start_temp = self.get_temp()
        delta_temp = temp_target - start_temp

        if delta_temp > 0:
            overshot_temp = min(temp_target + overshot, 99.9)
            undershot_temp = delta_temp - undershot
        else:
            overshot_temp = max(temp_target - overshot, -10)
            undershot_temp = delta_temp + undershot

        delay_min = self.time_to_reach_sample_temp(undershot_temp)
        delay_seconds = round((delay_min * 60), 1)

        """
        if not self.protocol.is_simulating:
            print("Temp Module is {}C on its way to {}C.".format(self.get_temp(), temp_target, delay_seconds))
        """

        if delay_min > 3:
            # Set temperature to rapidly cool or heat and delay program to allow temperature change.
            self.set_temperature(overshot_temp)
            """
            print("Delaying program for {} seconds to allow Temp Module to reach {}C."
                                  .format(delay_seconds, temp_target))
            """
            if not self.protocol.is_simulating():
                time.sleep(delay_seconds)

        self.set_temperature(temp_target)


class Utilities:
    def __init__(self, protocol):

        # TSV file location on OT-2
        tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)
        self.on_ot2 = True
        # If not on the OT-2, get temp TSV file location on Windows Computers for simulation
        if not os.path.isfile(tsv_file_path):
            self.on_ot2 = False
            tsv_file_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())

        self.parameter_file = tsv_file_path
        self.sample_dictionary = defaultdict(list)
        self.protocol = protocol
        self.args = None
        self.slot_list = \
            ["Slot1", "Slot2", "Slot3", "Slot4", "Slot5", "Slot6", "Slot7", "Slot8", "Slot9", "Slot10", "Slot11"]
        self.tipbox_dict = \
            {"p10_multi": "opentrons_96_tiprack_10ul", "p10_single": "opentrons_96_tiprack_10ul",
             "p20_single_gen2": ["opentrons_96_tiprack_20ul", "opentrons_96_filtertiprack_20ul"],
             "p300_single_gen2": ["opentrons_96_tiprack_300ul", "opentrons_96_filtertiprack_200ul"]
             }
        self._labware_dict = {}
        self._slot_dict = {}
        self._left_tiprack_list = []
        self._right_tiprack_list = []
        self.left_pipette = None
        self.right_pipette = None

    def dispense_reagent_mix(self, labware_dict, target_well_dict, target_info_dict, left_pipette, right_pipette):
        """
        This will dispense our master mixes into each well.

        @param target_info_dict:
        @param labware_dict:
        @param target_well_dict:
        @param left_pipette:
        @param right_pipette:
        @return:
        """

        sample_destination_labware = labware_dict[self.args.PCR_PlateSlot]

        # Dispense reagents into all wells
        reagent_aspirated = float(self.args.MasterMixPerRxn)

        for target in target_well_dict:
            reagent_slot = self.args.ReagentSlot
            if "Illumina_Dual_Indexing" in self.args.Template:
                reagent_source_well = self.args.PCR_ReagentWell
                reagent_well_vol = float(self.args.TotalReagentVolume)
            else:
                reagent_source_well = target_info_dict[int(target)][1]
                reagent_well_vol = float(target_info_dict[int(target)][2])

            target_well_list = target_well_dict[target]
            reagent_source_labware = labware_dict[reagent_slot]
            reagent_well_dia = reagent_source_labware[reagent_source_well].diameter

            reagent_pipette = \
                self.pipette_selection(left_pipette, right_pipette, float(self.args.MasterMixPerRxn))

            if "Illumina_Dual_Indexing" not in self.args.Template:
                self.protocol.comment("\nDispensing {} target with {}"
                                      .format(target_info_dict[int(target)][1], reagent_pipette))
            else:
                self.protocol.comment("\nDispensing Master Mix with {}".format(reagent_pipette))

            for well in target_well_list:
                reagent_tip_height = self.res_tip_height(reagent_well_vol - reagent_aspirated, reagent_well_dia)
                self.pipette_reagents(reagent_pipette,
                                         reagent_source_labware[reagent_source_well].bottom(reagent_tip_height),
                                         sample_destination_labware[well], float(self.args.MasterMixPerRxn),
                                         NewTip=False, MixReaction=False, touch=True, MixVolume=None
                                         )

                reagent_aspirated += float(self.args.MasterMixPerRxn)

            # Drop any tips the pipettes might have.
            if "Illumina_Dual_Indexing" not in self.args.Template:
                self.drop_any_tips([left_pipette, right_pipette])
        self.drop_any_tips([left_pipette, right_pipette])

    def drop_any_tips(self, pipettes):

        for pipette in pipettes:
            try:
                if pipette.has_tip:
                    pipette.drop_tip()
            except AttributeError:
                pass

    def pipette_reagents(self, pipette, source_location, destination_location, volume, NewTip, MixReaction,
                         touch=False, MixVolume=None):
        """
        Generic function to dispense material into designated well.
        @param MixVolume:
        @param pipette:
        @param source_location:
        @param destination_location:
        @param volume:
        @param NewTip:
        @param MixReaction:
        @param touch:
        @return:
        """

        def tip_touch():
            pipette.touch_tip(radius=0.79, v_offset=-2, speed=10)

        if NewTip:
            if pipette.has_tip:
                pipette.drop_tip()

        if not pipette.has_tip:
            pipette.pick_up_tip()

        pipette.aspirate(volume, source_location, rate=0.75)
        if touch:
            tip_touch()

        pipette.dispense(volume, destination_location, rate=0.75)

        if not MixReaction:
            pipette.blow_out()

        if MixReaction:
            v = float(self.args.PCR_Volume)
            if MixVolume:
                v = MixVolume
            vol = round(v * 0.65, ndigits=1)
            pipette.mix(repetitions=4, volume=vol, rate=2.0)
            pipette.blow_out()
            tip_touch()

        if touch:
            tip_touch()

        if NewTip:
            pipette.drop_tip()

        return pipette

    def res_tip_height(self, res_vol, well_dia):
        """
        Calculate the height of the liquid in a reservoir and return the value to set the pipette tip height.
        This works for both conical shapes and cylinders.
        @param res_vol:
        @param well_dia:
        @return:
        """
        cone_vol = self.labware_cone_volume(self.args.ReagentSlot)
        bottom_offset = float(self.args.BottomOffset)
        cone_height = (3 * cone_vol / (math.pi * ((well_dia / 2) ** 2)))

        if cone_vol > 500:
            offset = 11
            cone_offset = 5
        else:
            offset = 5
            cone_offset = 2

        if res_vol > cone_vol:
            height = ((res_vol - cone_vol) / (math.pi * ((well_dia / 2) ** 2))) - offset + cone_height
        else:
            height = (3 * res_vol / (math.pi * ((well_dia / 2) ** 2))) - cone_offset

        if height < bottom_offset:
            height = bottom_offset

        return round(height, ndigits=1)

    def labware_cone_volume(self, labware_slot):
        """
        Based on the labware and reservoir return the volume at which the cylinder shape transitions to the
        conical shape.

        @param labware_slot:
        @return:
        """
        cone_vol = 200
        labware = getattr(self.args, "Slot{}".format(labware_slot), "")

        if "_1.5ml" in labware:
            cone_vol = 415

        elif "_5000ul_" in labware:
            cone_vol = 1300

        return cone_vol

    @staticmethod
    def pipette_selection(left_pipette, right_pipette, volume):
        """
        Function to select pipette based on expected volumes.  Will also adjust volume is pipette needs to pick up >1x
        @param left_pipette:
        @param right_pipette:
        @param volume:
        @return:
        """
        # ToDo: This will not run on a FLEX and is error prone.  Need to allow more pipettes
        pipette = ""
        if volume > 20:
            if "P300 Single-Channel GEN2" in str(right_pipette):
                pipette = right_pipette
            else:
                pipette = left_pipette
        elif volume <= 20:
            if "P20 Single-Channel GEN2" in str(left_pipette):
                pipette = left_pipette
            else:
                pipette = right_pipette

        return pipette

    def dispense_water(self, water_well_dict, left_pipette, right_pipette):
        """
        This will dispense water into any wells that require it in the PCR destination plate/tubes only.
        Water for dilutions is dispensed as needed.

        @param water_well_dict:
        @param left_pipette:
        @param right_pipette:
        @return:
        """

        # reagent_labware = self._labware_dict[self.args.ReagentSlot]
        sample_destination_labware = self._labware_dict[self.args.PCR_PlateSlot]
        '''
        bottom_offset = float(args.BottomOffset)
        cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
        water_res_well_dia = reagent_labware[args.WaterResWell].diameter
        water_tip_height = Utilities.res_tip_height(float(args.WaterResVol), water_res_well_dia, cone_vol, bottom_offset)
        water_aspirated = 0
        '''

        destination_wells = []
        dispense_vol = []
        water_aspirated = 0

        for well in water_well_dict:
            destination_wells.append(sample_destination_labware[well])
            # destination_wells.append(well)
            dispense_vol.append(round(float(water_well_dict[well]), 2))
            water_aspirated += water_well_dict[well]

        water_aspirated = round(water_aspirated, 2)
        if min(dispense_vol) <= 9:
            volume = max(dispense_vol)
        else:
            volume = water_aspirated

        # Define the pipette for dispensing the water.
        water_pipette = self.pipette_selection(left_pipette, right_pipette, volume)
        self.protocol.comment("\nDistributing water with {} pipette".format(water_pipette))

        # Use custom distribute command to dispense water.
        self.distribute_reagents(water_pipette, destination_wells, dispense_vol)

        self.drop_any_tips([left_pipette, right_pipette])

        return water_aspirated

    def distribute_reagents(self, pipette, destination_wells, dispense_vol):
        """
        Dispense reagents using a custom distribute function.
        @param pipette:
        @param destination_wells:
        @param dispense_vol:
        """

        # ToDo: This needs work.
        p20_default_rate = 7.50
        p300_default_rate = 75.0
        #  p300_default_rate = 92.86
        # p20 = False
        source_well = self._labware_dict[self.args.ReagentSlot][self.args.WaterResWell]

        p20_tips = False
        p200_tips = False
        p300_tips = False

        if "p300_single_gen2" in self.protocol.params.left_pipette and "P300 Single-Channel GEN2" in str(pipette):
            if any("300" in str(s) for s in self._left_tiprack_list):
                p300_tips = True
            elif any("200" in str(s) for s in self._left_tiprack_list):
                p200_tips = True
        elif "p20_single_gen2" in self.protocol.params.left_pipette and "P20 Single-Channel GEN2" in str(pipette):
            p20_tips = True
        elif "p300_single_gen2" in self.protocol.params.right_pipette and "P300 Single-Channel GEN2" in str(pipette):
            if any("300" in str(s) for s in self._left_tiprack_list):
                p300_tips = True
            elif any("200" in str(s) for s in self._left_tiprack_list):
                p200_tips = True
        elif "p20_single_gen2" in self.protocol.params.right_pipette and "P20 Single-Channel GEN2" in str(pipette):
            p20_tips = True

        if p20_tips:
            max_tip_vol = 19.0
        elif p200_tips:
            max_tip_vol = 195.0
        elif p300_tips:
            max_tip_vol = 295.0

        if "P300 Single-Channel GEN2" in str(pipette):
            # print("Distributing Water With P300 Single-Channel GEN2")
            # touch = False
            # r = 0.25
            # s = 10
            default_rate = p300_default_rate
            pipette.flow_rate.aspirate = 30
            pipette.flow_rate.dispense = 10
            pipette.flow_rate.blow_out = 50
            disposal_vol = 30

        elif "P20 Single-Channel GEN2" in str(pipette):
            print("Distributing Water With P20 Single-Channel GEN2")
            # p20 = True
            # touch = False
            # r = 0.80
            # s = 10
            default_rate = p20_default_rate
            pipette.flow_rate.aspirate = 6.5
            pipette.flow_rate.dispense = 5.0
            pipette.flow_rate.blow_out = 7.0
            disposal_vol = 2.0

        total_vol = 0
        # p20_vol = 0.0
        # p20_dispense_list = []
        # p20_destination_wells = []
        water_res_vol = float(self.args.WaterResVol)

        tip_vol = 0.0
        dispense_list = []
        well_distribution = []

        if not pipette.has_tip:
            pipette.pick_up_tip()

        i = 0

        # Trying to keep the tip from being submerged in the source well liquid
        for volume, dest_well in zip(dispense_vol, destination_wells):
            total_vol += volume
            tip_vol += volume
            dispense_list.append(volume)
            well_distribution.append(dest_well)
            i += 1

            # For some reason this value is getting 1e-5 added to it occasionally.  Rounding corrects this.
            tip_vol = round(tip_vol, 1)

            # Need to keep the volume in the tips below their max vol while dynamically changing the tip height.
            #  My hack to get a dispense like function that will keep the same tip
            if i == len(dispense_vol) or dispense_vol[i] + tip_vol + disposal_vol >= max_tip_vol:
                water_res_vol = round(water_res_vol, 1)
                height = self.res_tip_height(water_res_vol, source_well.diameter)
                aspirated_vol = tip_vol + disposal_vol
                pipette.aspirate(volume=aspirated_vol, location=source_well.bottom(height))

                for destination_well, dispensed_vol in zip(well_distribution, dispense_list):
                    pipette.dispense(volume=dispensed_vol, location=destination_well)

                pipette.blow_out(source_well)
                water_res_vol -= tip_vol
                tip_vol = 0.0
                del dispense_list[:i]
                del well_distribution[:i]

        # Reset flow rates to default values
        pipette.flow_rate.aspirate = default_rate
        pipette.flow_rate.dispense = default_rate
        pipette.flow_rate.blow_out = default_rate

    @ property
    def tipracks(self):
        return self._left_tiprack_list, self._right_tiprack_list

    @ property
    def deck_layout(self):
        return self._labware_dict, self._slot_dict

    def labware_parsing(self):
        for i in range(11):
            labware = getattr(self.args, "{}".format(self.slot_list[i]))

            if labware:
                self._slot_dict[str(i + 1)] = labware
                self._labware_dict[str(i + 1)] = self.protocol.load_labware(labware, str(i + 1))

                if labware in self.tipbox_dict[self.protocol.params.left_pipette]:
                    self._left_tiprack_list.append(self._labware_dict[str(i + 1)])
                elif labware in self.tipbox_dict[self.protocol.params.right_pipette]:
                    self._right_tiprack_list.append(self._labware_dict[str(i + 1)])

    def plate_layout(self, slot):
        """
        Define the destination layout for the reactions.  Can be 384-well, 96-well plate or 8-well strip tubes
        @return:
        """

        labware = self._slot_dict[slot]

        column_count = 0
        row_count = 0
        if "384_" in labware:
            column_count = 32
            row_count = 12
        elif "96_" or "8_well" or "ddpcr_plate" in labware:
            column_count = 12
            row_count = 8

        layout_data = defaultdict(list)
        plate_layout_by_column = []

        # This is the index when using strip tubes.
        column_index = [1, 3, 5, 7, 9, 11, 12]
        row_labels = \
            ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
             "U", "V", "W", "X", "Y", "Z"]

        for c in range(column_count):

            if "8_well" in labware:
                # Strip tubes are in every other column.
                c = column_index[c]
            else:
                # Humans don't do well with 0-based labels.
                c += 1

            for r in range(row_count):
                layout_data[row_labels[r]] = [''] * column_count
                plate_layout_by_column.append("{}{}".format(row_labels[r], c))

        return plate_layout_by_column, layout_data

    def parse_sample_template(self):
        """
        Parse the TSV file and return data objects to run def.
        @return:
        """
        line_num = 0
        options_dictionary = defaultdict(str)
        index_file = list(csv.reader(open(self.parameter_file), delimiter='\t'))

        for line in index_file:
            if line_num == 0:
                options_dictionary["Version"] = line[1]
                options_dictionary["Template"] = line[0].strip("#")
            line_num += 1
            col_count = len(line)
            tmp_line = []
            sample_key = ""
            if col_count > 0 and "#" not in line[0] and len(line[0].split("#")[0]) > 0:
                # Skip any lines that are blank or comments.
                for i in range(7):
                    try:
                        line[i] = line[i].split("#")[0]  # Strip out end of line comments.
                    except IndexError:
                        continue

                    if i == 0 and "--" in line[0]:
                        key = line[0].strip('--')
                        key_value = line[1]
                        if options_dictionary["Template"] != " Illumina_Dual_Indexing" and "Target_" in key or "PositiveControl_" in key:
                            try:
                                key_value = (line[1], line[2], line[3])
                            except IndexError:
                                pass
                            options_dictionary[key] = key_value

                        options_dictionary[key] = key_value
                    elif "--" not in line[0] and int(line[0]) < 12:
                        sample_key = line[0], line[1]
                        tmp_line.append(line[i])
                if sample_key:
                    self.sample_dictionary[sample_key] = tmp_line

        self.args = SimpleNamespace(**options_dictionary)
        return self.sample_dictionary, self.args
