import datetime
import os
# import sys
import csv
import platform
import serial
import time
from types import SimpleNamespace
from contextlib import suppress
from collections import defaultdict
from opentrons import protocol_api
from opentrons.simulate import simulate, format_runlog
import math
# import Tool_Box

# metadata
metadata = {
    'protocolName': 'PCR v3.0.0b',
    'author': 'Dennis Simpson <dennis@email.unc.edu>',
    'description': 'Setup a ddPCR or Generic PCR'
    }

# requirements
requirements = {"robotType": "OT-2", "apiLevel": "2.20"}

def add_parameters(parameters: protocol_api.Parameters):

    """
    Parse the TSV file and fill in some parameter information.  This is duplicated from Utilities.  I don't know
    another method to pass the information.
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
                        key_value = (line[1], line[2], line[3])
                    options_dictionary[key] = key_value
                elif "--" not in line[0] and int(line[0]) < 12:
                    sample_key = line[0], line[1]
                    tmp_line.append(line[i])
            if sample_key:
                sample_dictionary[sample_key] = tmp_line

    args = SimpleNamespace(**options_dictionary)

    # We have limited space for the run_label.  To make sure the label is unique, I use a unix timestamp for the run_date.
    run_date = datetime.datetime.today().strftime("%f")
    run_type = args.Template.split(" ")[1]

    # This is used by the Opentrons app to make each run unique.
    parameters.add_str(
        variable_name="run_label",
        display_name="Begin {} for {} {}".format(run_type, args.User, run_date),
        choices=[
            {"display_name": "Run Label", "value": "Begin {} for {} {}".format(run_type, args.User, run_date)},],
        default="Begin {} for {} {}".format(run_type, args.User, run_date),
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


def movedparse_sample_template(input_file):
    """
    Parse the TSV file and return data objects to run def.
    @param input_file:
    @return:
    """
    line_num = 0
    options_dictionary = defaultdict(str)
    sample_dictionary = defaultdict(list)
    index_file = list(csv.reader(open(input_file), delimiter='\t'))
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
                        key_value = (line[1], line[2], line[3])
                    options_dictionary[key] = key_value
                elif "--" not in line[0] and int(line[0]) < 12:
                    sample_key = line[0], line[1]
                    tmp_line.append(line[i])
            if sample_key:
                sample_dictionary[sample_key] = tmp_line
    return sample_dictionary, SimpleNamespace(**options_dictionary)


def movedlabware_parsing(args, protocol):
    # Extract Slot information
    slot_list = ["Slot1", "Slot2", "Slot3", "Slot4", "Slot5", "Slot6", "Slot7", "Slot8", "Slot9", "Slot10", "Slot11"]
    labware_dict = {}
    slot_dict = {}
    tipbox_dict = \
        {"p10_multi": "opentrons_96_tiprack_10ul", "p10_single": "opentrons_96_tiprack_10ul",
         "p20_single_gen2": ["opentrons_96_tiprack_20ul", "opentrons_96_filtertiprack_20ul"],
         "p300_single_gen2": ["opentrons_96_tiprack_300ul", "opentrons_96_filtertiprack_300ul"]
         }

    # Pipette Tip Boxes
    left_tiprack_list = []
    right_tiprack_list = []
    for i in range(11):
        labware = getattr(args, "{}".format(slot_list[i]))

        if labware:
            slot_dict[str(i + 1)] = labware
            labware_dict[str(i + 1)] = protocol.load_labware(labware, str(i + 1))

            if labware in tipbox_dict[protocol.params.left_pipette]:
                left_tiprack_list.append(labware_dict[str(i + 1)])
            elif labware in tipbox_dict[protocol.params.right_pipette]:
                right_tiprack_list.append(labware_dict[str(i + 1)])

    return labware_dict, slot_dict, left_tiprack_list, right_tiprack_list


def plate_layout(labware):
    """
    Define the destination layout for the reactions.  Can be 384-well, 96-well plate or 8-well strip tubes
    @param labware:
    @return:
    """
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


def sample_processing(args, sample_parameters, target_info_dict, slot_dict):
    sample_data_dict = defaultdict(list)
    target_well_dict = defaultdict(list)
    water_well_dict = defaultdict(float)
    plate_layout_by_column, layout_data = plate_layout(slot_dict[args.PCR_PlateSlot])
    used_wells = []
    dest_well_count = 0
    target_list = []

    for sample_key in sample_parameters:
        # sample_source_slot = sample_parameters[sample_key][0]
        # sample_source_well = sample_parameters[sample_key][1]
        sample_name = sample_parameters[sample_key][2]
        sample_concentration = float(sample_parameters[sample_key][3])
        sample_targets = sample_parameters[sample_key][4].split(",")
        replicates = int(sample_parameters[sample_key][5])

        if "ddPCR" in args.Template:
            template_in_rxn = float(args.DNA_in_Reaction)
        elif "Generic PCR" in args.Template:
            template_in_rxn = float(sample_parameters[sample_key][6])
        else:
            raise SystemExit("There is an error in the Template name.")

        # sample_string = sample_name
        ucv = calculate_volumes(args, sample_concentration, template_in_rxn)
        sample_vol = ucv[0]
        diluent_vol = ucv[1]
        diluted_sample_vol = ucv[2]
        reaction_water_vol = ucv[3]
        max_template_vol = ucv[4]

        sample_wells = []
        for target in sample_targets:
            target_list.append(target)
            target_name = target_info_dict[int(target)][1]

            for i in range(replicates):
                well = plate_layout_by_column[dest_well_count]
                row = well[0]
                column = int(well[1:])-1
                s_volume = diluted_sample_vol

                if diluent_vol == 0:
                    dilution = "Neat"
                    s_volume = sample_vol
                else:
                    dilution = "1:{}".format(int((sample_vol + diluent_vol) / sample_vol))

                layout_data[row][column] = "{}|{}|{}|{}"\
                    .format(sample_name, target_name, dilution, s_volume)

                water_well_dict[well] = reaction_water_vol
                target_well_dict[target].append(well)
                sample_wells.append(well)
                used_wells.append(well)
                dest_well_count += 1
                
        sample_data_dict[sample_key] = \
            [round(sample_vol, ndigits=1), round(diluent_vol, ndigits=1), round(diluted_sample_vol, ndigits=1), sample_wells]

    # Define our no template control wells for the targets.
    for target in target_well_dict:
        target_list.append(target)

    target_list = list(set(target_list))

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
    # Turn off rail lights.
    protocol.set_rail_lights(False)

    utility = Utilities(protocol)
    sample_parameters, args = utility.parse_sample_template()
    utility.labware_parsing()
    protocol.comment(protocol.params.run_label)

    robot = False
    if os.path.isfile("{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)):
        robot = True

    if bool(args.UseTemperatureModule):
        temp_mod = ColdPlateSlimDriver(protocol, robot)
        temp_mod.set_temperature(args.Temperature)

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

    target_info_dict = defaultdict(list)

    # Read targeting parameters into dictionary
    for i in range(10):
        target = getattr(args, "Target_{}".format(i + 1))

        if len(target[0]) > 1:
            """
            if not all('' == s or s.isspace() for s in target):
                target_info_dict[i + 1] = target
            """
            target_info_dict[i + 1] = target

    sample_data_dict, water_well_dict, target_well_dict, used_wells, layout_data, max_template_vol = \
        sample_processing(args, sample_parameters, target_info_dict, slot_dict)

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
    water_aspirated = dispense_water(args, labware, water_well_dict, left_pipette, right_pipette)

    dispense_reagent_mix(args, labware, target_well_dict, target_info_dict, left_pipette, right_pipette)

    water_aspirated = dispense_samples(args, labware, sample_data_dict, slot_dict, sample_parameters, left_pipette,
                                       right_pipette, water_aspirated)
    if "ddPCR" in args.Template:
        fill_empty_wells(args, used_wells, water_aspirated, labware, left_pipette, right_pipette)

    if args.UseTemperatureModule == "True":
        protocol.set_rail_lights(True)
        protocol.comment("\nProgram is complete. Click RESUME to exit.  Temperature is holding at {}"
                         .format(temp_mod.get_temp()))
        protocol.pause()
        temp_mod.deactivate()
        protocol.set_rail_lights(False)
    else:
        protocol.comment("\nProgram Complete")

    if not protocol.is_simulating():
        os.remove(utility.parameter_file)


def fill_empty_wells(args, used_wells, water_aspirated, labware_dict, left_pipette, right_pipette):
    """
    This will fill the remaining wells in a column with water.  Needed to for the droplet generator.
    """
    bottom_offset = float(args.BottomOffset)
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
        # cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
        # cone_vol = Utilities.labware_cone_volume(args, args.ReagentSlot)
        cone_vol = labware_cone_volume(args, args.ReagentSlot)
        fill_pipette, fill_loop, fill_vol = \
            pipette_selection(left_pipette, right_pipette, float(args.PCR_Volume))
        water_tip_height = \
            res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                     bottom_offset)

        for i in range(wells_remaining):
            blank_well = "{}{}".format(row_list[i+row_index+1], column)

            dispensing_loop(args, fill_loop, fill_pipette,
                                      reagent_labware[args.WaterResWell].bottom(water_tip_height),
                                      sample_destination_labware[blank_well], fill_vol,
                                      NewTip=False, MixReaction=False)

            water_aspirated = water_aspirated+fill_vol
            water_tip_height = res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia,
                                                        cone_vol, bottom_offset)
        fill_pipette.drop_tip()


def dispense_reagent_mix(args, labware_dict, target_well_dict, target_info_dict, left_pipette, right_pipette):
    """
    This will dispense our master mixes into each well.

    @param target_info_dict:
    @param args:
    @param labware_dict:
    @param target_well_dict:
    @param left_pipette:
    @param right_pipette:
    @return:
    """

    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    bottom_offset = float(args.BottomOffset)

    # Dispense reagents into all wells
    for target in target_well_dict:
        reagent_slot = args.ReagentSlot
        reagent_source_well = target_info_dict[int(target)][0]
        reagent_source_labware = labware_dict[reagent_slot]
        target_well_list = target_well_dict[target]

        reagent_aspirated = float(args.MasterMixPerRxn)
        reagent_well_vol = float(target_info_dict[int(target)][2])
        reagent_well_dia = reagent_source_labware[reagent_source_well].diameter
        # reagent_cone_vol = Utilities.labware_cone_volume(args, reagent_source_labware)
        reagent_cone_vol = labware_cone_volume(args, reagent_source_labware)
        reagent_pipette, reagent_loop, reagent_volume = \
            pipette_selection(left_pipette, right_pipette, reagent_aspirated)
        '''
        Distribute does not work with Multiplex master mixes.
        # Use distribute command to dispense master mix.
        destination_wells = []
        for well in target_well_list:
            destination_wells.append(sample_destination_labware[well])

        distribute_reagents(right_pipette,
                            reagent_source_labware[reagent_source_well].bottom(z=bottom_offset),
                            destination_wells,
                            float(args.MasterMixPerRxn)
                            )
        
        # Drop any tips the pipettes might have.
        if left_pipette.has_tip:
            left_pipette.drop_tip()
        if right_pipette.has_tip:
            right_pipette.drop_tip()

        '''
        for well in target_well_list:
            reagent_tip_height = res_tip_height(reagent_well_vol-reagent_aspirated, reagent_well_dia,
                                                          reagent_cone_vol, bottom_offset)

            dispensing_loop(args, reagent_loop, reagent_pipette,
                                      reagent_source_labware[reagent_source_well].bottom(reagent_tip_height),
                                      sample_destination_labware[well], reagent_volume,
                                      NewTip=False, MixReaction=False, touch=True, MixVolume=None)

            reagent_aspirated += reagent_volume

        # Drop any tips the pipettes might have.
        if left_pipette.has_tip:
            left_pipette.drop_tip()
        if right_pipette.has_tip:
            right_pipette.drop_tip()

    return


def dispense_water(args, labware_dict, water_well_dict, left_pipette, right_pipette):
    """
    This will dispense water into any wells that require it in the PCR destination plate/tubes only.
    Water for dilutions is dispensed as needed.

    @param args:
    @param labware_dict:
    @param water_well_dict:
    @param left_pipette:
    @param right_pipette:
    @return:
    """
    reagent_labware = labware_dict[args.ReagentSlot]
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    '''
    bottom_offset = float(args.BottomOffset)
    cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
    water_res_well_dia = reagent_labware[args.WaterResWell].diameter
    water_tip_height = Utilities.res_tip_height(float(args.WaterResVol), water_res_well_dia, cone_vol, bottom_offset)
    water_aspirated = 0
    '''
    # left_pipette = utility.left_pipette
    # right_pipette = utility.right_pipette

    destination_wells = []
    dispense_vol = []
    water_aspirated = 0
    for well in water_well_dict:
        destination_wells.append(sample_destination_labware[well])
        dispense_vol.append(round(float(water_well_dict[well]), 2))
        water_aspirated += water_well_dict[well]

    # Define the pipette for dispensing the water.
    water_pipette, water_loop, water_volume = pipette_selection(left_pipette, right_pipette, water_aspirated)

    # Use distribute command to dispense water.
    distribute_reagents(water_pipette,
                        reagent_labware[args.WaterResWell].bottom(z=float(args.BottomOffset)),
                        destination_wells,
                        dispense_vol
                        )
    
    if left_pipette.has_tip:
        left_pipette.drop_tip()
    if right_pipette.has_tip:
        right_pipette.drop_tip()

    return water_aspirated


def dispense_samples(args, labware_dict, sample_data_dict, slot_dict, sample_parameters, left_pipette, right_pipette,
                     water_aspirated):
    """
    Dilute and dispense samples
    @param slot_dict:
    @param args:
    @param labware_dict:
    @param sample_data_dict:
    @param sample_parameters:
    @param left_pipette:
    @param right_pipette:
    @param water_aspirated:
    """

    try:
        dilution_labware = labware_dict[args.DilutionPlateSlot]
    except KeyError:
        dilution_labware = ""
    bottom_offset = float(args.BottomOffset)
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    reagent_labware = labware_dict[args.ReagentSlot]
    # water_res_well_dia = reagent_labware[args.WaterResWell].diameter
    water_res_well_dia = labware_dict[args.ReagentSlot][args.WaterResWell].diameter

    # cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
    # cone_vol = Utilities.labware_cone_volume(args, args.ReagentSlot)
    cone_vol = labware_cone_volume(args, args.ReagentSlot)
    water_tip_height = res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                                bottom_offset)
    # If the user determines no dilutions are required they can leave that slot blank.  I don't like this approach,
    # users could leave the information out and dilutions might still be required.
    if dilution_labware:
        labware_name = slot_dict[args.DilutionPlateSlot]
    else:
        labware_name = slot_dict[args.PCR_PlateSlot]

    dilution_plate_layout, unused_layout = plate_layout(labware_name)

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
            mix_volume = 18

        # If no dilution is necessary, dispense sample and continue
        if diluted_sample_vol == 0:
            sample_pipette, sample_loop, sample_vol = \
                pipette_selection(left_pipette, right_pipette, sample_vol)

            for well in sample_dest_wells:
                dispensing_loop(args, sample_loop, sample_pipette,
                                          sample_source_labware[sample_source_well],
                                          sample_destination_labware[well], sample_vol,
                                          NewTip=True, MixReaction=True, touch=True, MixVolume=mix_volume)
            continue

        # Adjust volume of diluted sample to make sure there is enough
        diluted_template_needed = round(diluted_sample_vol*(len(sample_dest_wells)+1.5), ndigits=1)
        diluted_template_factor = round(diluted_template_needed/(sample_vol+diluent_vol), ndigits=1)
        adjusted_sample_vol = round((sample_vol * diluted_template_factor), ndigits=1)
        diluent_vol = round((diluent_vol*diluted_template_factor), ndigits=1)

        # Reset the pipettes for the new volumes
        diluent_pipette, diluent_loop, diluent_vol = \
            pipette_selection(left_pipette, right_pipette, diluent_vol)
        sample_pipette, sample_loop, sample_vol = \
            pipette_selection(left_pipette, right_pipette, adjusted_sample_vol)

        # Make dilution, diluent first
        dilution_well = dilution_plate_layout[dilution_well_index]

        dispensing_loop(args, diluent_loop, diluent_pipette,
                                  reagent_labware[args.WaterResWell].bottom(water_tip_height),
                                  dilution_labware[dilution_well], diluent_vol, NewTip=True, MixReaction=False,
                                  touch=True)
        mix_volume = None
        if diluted_sample_vol < 20:
            mix_volume = 18
        dispensing_loop(args, sample_loop, sample_pipette, sample_source_labware[sample_source_well],
                                  dilution_labware[dilution_well], sample_vol, NewTip=True, MixReaction=True,
                                  MixVolume=mix_volume)

        water_aspirated += diluent_vol
        dilution_well_index += 1
        water_tip_height = \
            res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                     bottom_offset)

        # Add diluted sample to PCR plate
        for well in sample_dest_wells:
            sample_pipette, sample_loop, diluted_sample_vol = \
                pipette_selection(left_pipette, right_pipette, diluted_sample_vol)
            mix_volume = None
            if diluted_sample_vol < 20:
                mix_volume = 18
            dispensing_loop(args, sample_loop, sample_pipette,
                                      dilution_labware[dilution_well].bottom(bottom_offset),
                                      sample_destination_labware[well], diluted_sample_vol, NewTip=True,
                                      MixReaction=True, MixVolume=mix_volume)

        if sample_pipette.has_tip:
            sample_pipette.drop_tip()

    return water_aspirated


def distribute_reagents(pipette, source_well, destination_wells, dispense_vol):
    """
    Dispense reagents using the distribute function.
    @param pipette:
    @param source_well:
    @param destination_wells:
    @param dispense_vol:
    """

    # ToDo: This needs work.
    p20_default_rate = 7.50
    p300_default_rate = 75.0
   #  p300_default_rate = 92.86

    if "P300 Single-Channel GEN2" in str(pipette):
        default_rate = p300_default_rate
        pipette.flow_rate.aspirate = 30
        pipette.flow_rate.dispense = 10
        pipette.flow_rate.blow_out = 50
    elif "P20 Single-Channel GEN2" in str(pipette):
        default_rate = p20_default_rate
        pipette.flow_rate.aspirate = 7.0
        pipette.flow_rate.dispense = 5.0
        pipette.flow_rate.blow_out = 7.0

    pipette.distribute(volume=dispense_vol, source=source_well, dest=destination_wells,
                       touch_tip=True, radius=0.75, v_offset=-4, speed=40, blow_out=True, disposal_volume=10,
                       blowout_location='source well')

    pipette.flow_rate.aspirate = default_rate
    pipette.flow_rate.dispense = default_rate
    pipette.flow_rate.blow_out = default_rate

def labware_cone_volume(args, labware_slot):
        """
        Based on the labware and reservoir return the volume at which the cylinder shape transitions to the conical shape.
        @param args:
        @param labware_slot:
        @return:
        """
        cone_vol = 200

        # labware = getattr(args, "Slot{}".format(str(labware_name)[-1:]))
        # labware = getattr(args, "Slot{}".format(str(labware_slot)))
        labware = getattr(args, "Slot{}".format(labware_slot), "")

        if "e1.5ml" in labware:
            cone_vol = 450

        elif "e5ml" in labware:
            cone_vol = 1200

        return cone_vol


def res_tip_height(res_vol, well_dia, cone_vol, bottom_offset):
        """
        Calculate the height of the liquid in a reservoir and return the value to set the pipette tip height.
        This works for both conical shapes and cylinders.
        @param bottom_offset:
        @param res_vol:
        @param well_dia:
        @param cone_vol:
        @return:
        """
        if res_vol > cone_vol:
            cone_height = (3 * cone_vol / (math.pi * ((well_dia / 2) ** 2)))
            height = ((res_vol - cone_vol) / (math.pi * ((well_dia / 2) ** 2))) - 5 + cone_height
        else:
            height = (3 * res_vol / (math.pi * ((well_dia / 2) ** 2))) - 3

        if height < 3:
            height = bottom_offset

        return round(height, ndigits=1)


def pipette_selection(left_pipette, right_pipette, volume):
        """
        Function to select pipette based on expected volumes.  Will also adjust volume is pipette needs to pick up >1x
        @param left_pipette:
        @param right_pipette:
        @param volume:
        @return:
        """
        # ToDo: This will not run on a FLEX and is error prone.  Need to allow more pipettes
        loop = 1
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

        return pipette, loop, round(volume, ndigits=1)


def dispensing_loop(args, loop_count, pipette, source_location, destination_location, volume, NewTip, MixReaction,
                        touch=False, MixVolume=None):
        """
        Generic function to dispense material into designated well.
        @param MixVolume:
        @param args:
        @param loop_count:
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
            pipette.touch_tip(radius=0.75, v_offset=-8, speed=40)

        if NewTip:
            if pipette.has_tip:
                pipette.drop_tip()

        if not pipette.has_tip:
            pipette.pick_up_tip()

        while loop_count > 0:
            pipette.aspirate(volume, source_location, rate=0.75)

            if touch:
                tip_touch()

            pipette.dispense(volume, destination_location, rate=0.75)
            loop_count -= 1

            if not MixReaction:
                pipette.blow_out()
                if touch:
                    tip_touch()

        if MixReaction:
            v = float(args.PCR_Volume)
            if MixVolume:
                v = MixVolume
            vol = round(v*0.65, ndigits=1)
            pipette.mix(repetitions=4, volume=vol, rate=2.0)
            pipette.blow_out()
            tip_touch()

        if NewTip:
            pipette.drop_tip()

        return pipette


class ColdPlateSlimDriver:
    """
    From Parhelia.  Class to control their temperature module.
    @todo Need to find out if the Opentrons built in module will work.
    """
    def __init__(self, protocol_context, temp_mode_number=0, robot=None):
        self.serial_number = "29533"
        self.device_name = "/dev/ttyUSB" + str(temp_mode_number)
        self.baudrate = 9600
        self.bytesize = serial.EIGHTBITS
        self.parity = serial.PARITY_NONE
        self.stopbits = serial.STOPBITS_ONE
        self.read_timeout = 2
        self.write_timeout = 2
        self.height = 45
        # self.temp = 20
        self.temp = None
        self.protocol = protocol_context

        # check context, skip if simulating Linux
        if protocol_context.is_simulating() and not robot:
            self.protocol.comment("Simulation detected. Initializing Temperature Module in the dummy mode\n")
            self.serial_object = None
        else:
            self.protocol.comment("Execution mode detected.Initializing Temperature Module in the real deal mode\n")
            self.serial_object = serial.Serial(
                port=self.device_name,
                baudrate=self.baudrate,
                bytesize=self.bytesize,
                parity=self.parity,
                stopbits=self.stopbits,
                timeout=self.read_timeout,
                write_timeout=self.write_timeout,
            )

    @property
    def temperature(self):
        return self.get_temp()

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
            return "dummy response"

        output_lines = self.serial_object.readlines()
        output_string = "".join(l.decode("utf-8") for l in output_lines)
        return output_string

    def _send_command(self, my_command):
        """
        Worker function
        """
        SERIAL_ACK = "\r\n"

        command = my_command
        command += SERIAL_ACK

        if self.serial_object is None:
            print("sending dummy command: " + my_command)
            return

        self.serial_object.write(command.encode())
        self.serial_object.flush()
        return self._read_response()

    def get_info(self):
        if self.serial_object is None:
            return "dummy info"
        return self._send_command("info")

    def get_temp(self):
        if self.serial_object is None:
            return self.temp
        t = self._send_command("getTempActual")
        return float(t)

    def set_temperature(self, my_temp):
        if self.serial_object is None:
            self.temp = my_temp
            return
        # temp = float(my_temp) * 10
        # temp = int(temp)
        temp = int(my_temp)
        self._send_command(f"setTempTarget{temp:03}")
        self._send_command("tempOn")
        self._set_temp_andWait(my_temp)

    def _set_temp_andWait(self, target_temp, timeout_min=30, tolerance=1.0):
        interval_sec = 10
        seconds_in_min = 60
        time_elapsed = 0

        while abs(self.get_temp() - target_temp) > tolerance:
            self.protocol.comment("Temp Module is {}C on its way to {}C.  Waiting for {} seconds\n"
                                  .format(self.get_temp(), target_temp, time_elapsed))

            if not self.protocol.is_simulating():  # Skip delay during simulation
                time.sleep(interval_sec)
            time_elapsed += interval_sec
            if time_elapsed > timeout_min * seconds_in_min:
                raise Exception("Temperature timeout")

        # return target_temp

    def deactivate(self):
        if self.serial_object is None:
            self.temp = 25
        else:
            self._send_command("tempOff")
            self.serial_object.close()

    def time_to_reach_sample_temp(self, delta_temp):
        x = delta_temp
        if (x > 0):
            time_min = 0.364 + 0.559 * x - 0.0315 * x ** 2 + 1.29E-03 * x ** 3 - 2.46E-05 * x ** 4 + 2.21E-07 * x ** 5 - 7.09E-10 * x ** 6
        else:
            time_min = -0.1 - 0.329 * x - 0.00413 * x ** 2 - 0.0000569 * x ** 3 + 0.0000000223 * x ** 4
        return time_min

    def quick_temp(self, temp_target, overshot=10, undershot=2):
        if testmode:
            pass
        else:
            start_temp = self.get_temp()
            delta_temp = temp_target - start_temp

            if delta_temp > 0:
                overshot_temp = min(temp_target + overshot, 99.9)
                undershot_temp = delta_temp - undershot
            else:
                overshot_temp = max(temp_target - overshot, -10)
                undershot_temp = delta_temp + undershot

            delay_min = self.time_to_reach_sample_temp(undershot_temp)
            delay_seconds = delay_min * 60

            self.protocol.comment(
                f"Setting temperature Target temp: {temp_target}\n C. Temp transition time: {delay_min} minutes")
            self.set_temp(overshot_temp)
            if not self.protocol.is_simulating():
                time.sleep(delay_seconds)
            self.set_temp(temp_target)

        if testmode:
            self.protocol.comment(
                f"Setting temperature Target temp: {temp_target}\n C.Ramp skipped for testmode")
            self.set_temp(temp_target)

class Utilities:
    def __init__(self, protocol):

        # TSV file location on OT-2
        tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)

        # If not on the OT-2, get temp TSV file location on Windows Computers for simulation
        if not os.path.isfile(tsv_file_path):
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
             "p300_single_gen2": ["opentrons_96_tiprack_300ul", "opentrons_96_filtertiprack_300ul"]
             }
        self._labware_dict = {}
        self._slot_dict = {}
        self._left_tiprack_list = []
        self._right_tiprack_list = []
        self.left_pipette = None
        self.right_pipette = None
        # self.labware_parsing()

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
                        if "Target_" in key or "PositiveControl_" in key:
                            key_value = (line[1], line[2], line[3])
                        options_dictionary[key] = key_value
                    elif "--" not in line[0] and int(line[0]) < 12:
                        sample_key = line[0], line[1]
                        tmp_line.append(line[i])
                if sample_key:
                    self.sample_dictionary[sample_key] = tmp_line

        self.args = SimpleNamespace(**options_dictionary)
        return self.sample_dictionary, self.args


if __name__ == "__main__":
    protocol_file = open('PCR.py')
    labware_path = "{}{}custom_labware".format(os.getcwd(), os.sep)
    run_log, __bundle__ = simulate(protocol_file, custom_labware_paths=[labware_path])
    run_date = datetime.datetime.today().strftime("%a %b %d %H:%M %Y")
    i = 1
    t = format_runlog(run_log).split("\n")

    outstring = "Opentrons OT-2 Steps for {}.\nDate:  {}\nProgram File: PCR.py\n\nStep\tCommand\n" \
        .format(metadata['protocolName'], run_date)

    for l in t:
        outstring += "{}\t{}\n".format(i, l)
        i += 1
    if platform.system() == "Windows":
        outfile = open("C:{0}Users{0}{1}{0}Documents{0}ProgramFileSimulation.txt"
                       .format(os.sep, os.getlogin()), 'w', encoding="UTF-16")
        outfile.write(outstring)
        outfile.close()
    protocol_file.close()
