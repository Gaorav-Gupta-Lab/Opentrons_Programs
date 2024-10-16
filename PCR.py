import datetime
import os
import sys
import csv
import platform
from types import SimpleNamespace
from contextlib import suppress
from collections import defaultdict
from opentrons import protocol_api
from opentrons.simulate import simulate, format_runlog
import math
# import Tool_Box

# Check if we are on the OT-2, Robotron, or some other computer.
template_parser_path = "{0}var{0}lib{0}jupyter{0}notebooks".format(os.sep)
if not os.path.exists(template_parser_path):
    # This is Robotron, the OT-2 controller computer.
    template_parser_path = "C:{0}Opentrons_Programs".format(os.sep)

    # This is the development computer
    if not os.path.exists(template_parser_path):
        template_parser_path = \
            "C:/Users/dennis/OneDrive - University of North Carolina at Chapel Hill/Projects/Programs/Opentrons_Programs"
sys.path.insert(0, template_parser_path)

# metadata
metadata = {
    'protocolName': 'PCR v2.1.0',
    'author': 'Dennis Simpson <dennis@email.unc.edu>',
    'description': 'Setup a Generic PCR or a ddPCR'
    }

# requirements
requirements = {"robotType": "OT-2", "apiLevel": "2.20"}

def add_parameters(parameters: protocol_api.Parameters):
    # TSV file location on OT-2
    tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)

    # If not on the OT-2, get temp TSV file location on Windows Computers for simulation
    if not os.path.isfile(tsv_file_path):
        tsv_file_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())

    sample_parameters, args = parse_sample_template(tsv_file_path)

    # run_date = datetime.datetime.today().strftime("%a %b %d %H:%M %Y")
    # run_date = datetime.datetime.today().strftime("%d-%b-%Y")
    run_date = datetime.datetime.today().strftime("%f")
    run_type = args.Template.split(" ")[1]

    # This is used by the Opentrons app to make each run unique.
    parameters.add_str(
        variable_name="run_label",
        display_name="Run Label",
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
        parameters.add_bool(
        variable_name="{}1".format(run_type),
        display_name="{} {} {}".format(run_type, args.User, run_date),
        description="Loading run information.",
        default=True
    )
    
    parameters.add_csv_file(
        variable_name="dragons_run",
        display_name="Dragon Hunting",
        description="Looking to initialize Opentrons csv commands."
    )
    """

def parse_sample_template(input_file):
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


def labware_parsing(args, protocol):
    # Extract Slot information
    slot_list = ["Slot1", "Slot2", "Slot3", "Slot4", "Slot5", "Slot6", "Slot7", "Slot8", "Slot9", "Slot10", "Slot11"]
    labware_dict = {}
    slot_dict = {}
    tipbox_dict = \
        {"p10_multi": "opentrons_96_tiprack_10ul", "p10_single": "opentrons_96_tiprack_10ul",
         "p20_single_gen2": ["opentrons_96_tiprack_20ul", "opentrons_96_filtertiprack_20ul"],
         "p300_single_gen2": ["opentrons_96_tiprack_300ul", "opentrons_96_filtertiprack_300ul"]}
    # Pipette Tip Boxes
    left_tiprack_list = []
    right_tiprack_list = []
    for i in range(11):
        labware = getattr(args, "{}".format(slot_list[i]))

        if labware:
            slot_dict[str(i + 1)] = labware
            labware_dict[str(i + 1)] = protocol.load_labware(labware, str(i + 1))
            # if labware in tipbox_dict[args.LeftPipette]:
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
    elif "96_" or "8_well" in labware:
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

    # TSV file location on OT-2
    tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)

    # If not on the OT-2, get temp TSV file location on Windows Computers for simulation
    if not os.path.isfile(tsv_file_path):
        tsv_file_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())

    sample_parameters, args = parse_sample_template(tsv_file_path)
    protocol.comment(protocol.params.run_label)

    # Turn on rail lights and pause program so user can load robot deck.
    # protocol.set_rail_lights(True)
    # protocol.pause("Load Labware onto robot deck and click resume when ready to continue")
    # protocol.home()
    protocol.set_rail_lights(False)
    labware_dict, slot_dict, left_tiprack_list, right_tiprack_list = labware_parsing(args, protocol)

    # Pipettes
    left_pipette = protocol.load_instrument(protocol.params.left_pipette, 'left', tip_racks=left_tiprack_list)
    right_pipette = protocol.load_instrument(protocol.params.right_pipette, 'right', tip_racks=right_tiprack_list)

    # Set the location of the first tip in box.
    with suppress(IndexError):
        left_pipette.starting_tip = left_tiprack_list[0].wells_by_name()[args.LeftPipetteFirstTip.upper()]
    with suppress(IndexError):
        right_pipette.starting_tip = right_tiprack_list[0].wells_by_name()[args.RightPipetteFirstTip.upper()]

    target_info_dict = defaultdict(list)

    # Read targeting parameters into dictionary
    for i in range(10):
        target = getattr(args, "Target_{}".format(i + 1))
        if target:
            if not all('' == s or s.isspace() for s in target):
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
    water_aspirated = dispense_water(args, labware_dict, water_well_dict, left_pipette, right_pipette)

    dispense_reagent_mix(args, labware_dict, target_well_dict, target_info_dict, left_pipette, right_pipette)

    water_aspirated = dispense_samples(args, labware_dict, sample_data_dict, slot_dict, sample_parameters, left_pipette,
                                       right_pipette, water_aspirated)
    if "ddPCR" in args.Template:
        fill_empty_wells(args, used_wells, water_aspirated, labware_dict, left_pipette, right_pipette)

    protocol.comment("\nProgram Complete")

    if not protocol.is_simulating():
        os.remove(tsv_file_path)


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
        cone_vol = Utilities.labware_cone_volume(args, args.ReagentSlot)
        fill_pipette, fill_loop, fill_vol = \
            Utilities.pipette_selection(left_pipette, right_pipette, float(args.PCR_Volume))
        water_tip_height = \
            Utilities.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                     bottom_offset)

        for i in range(wells_remaining):
            blank_well = "{}{}".format(row_list[i+row_index+1], column)

            Utilities.dispensing_loop(args, fill_loop, fill_pipette,
                                      reagent_labware[args.WaterResWell].bottom(water_tip_height),
                                      sample_destination_labware[blank_well], fill_vol,
                                      NewTip=False, MixReaction=False)

            water_aspirated = water_aspirated+fill_vol
            water_tip_height = Utilities.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia,
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
        '''
        reagent_aspirated = float(args.MasterMixPerRxn)
        reagent_well_vol = float(target_info_dict[int(target)][2])
        reagent_well_dia = reagent_source_labware[reagent_source_well].diameter
        reagent_cone_vol = Utilities.labware_cone_volume(args, reagent_source_labware)
        reagent_pipette, reagent_loop, reagent_volume = \
            Utilities.pipette_selection(left_pipette, right_pipette, reagent_aspirated)
        '''
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
            reagent_tip_height = Utilities.res_tip_height(reagent_well_vol-reagent_aspirated, reagent_well_dia,
                                                          reagent_cone_vol, bottom_offset)

            Utilities.dispensing_loop(args, reagent_loop, reagent_pipette,
                                      reagent_source_labware[reagent_source_well].bottom(reagent_tip_height),
                                      sample_destination_labware[well], reagent_volume,
                                      NewTip=False, MixReaction=False, touch=True, MixVolume=None)

            reagent_aspirated += reagent_volume

        # Drop any tips the pipettes might have.
        if left_pipette.has_tip:
            left_pipette.drop_tip()
        if right_pipette.has_tip:
            right_pipette.drop_tip()
        '''
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

    destination_wells = []
    dispense_vol = []
    water_aspirated = 0
    for well in water_well_dict:
        destination_wells.append(sample_destination_labware[well])
        dispense_vol.append(float(water_well_dict[well]))
        water_aspirated += water_well_dict[well]

    # Define the pipette for dispensing the water.
    water_pipette, water_loop, water_volume = Utilities.pipette_selection(left_pipette, right_pipette, water_aspirated)

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

    """
    # This block is only required for stand alone simulations.  Our GUI does not need it.
    for well in water_well_dict:
        water_volume = water_well_dict[well]

        # There is a bug in the Opentrons software that results in a 0 uL aspiration defaulting to pipette max.
        if water_volume < 0.1:
            continue

        # Define the pipette for dispensing the water.
        water_pipette, water_loop, water_volume = \
            Utilities.pipette_selection(left_pipette, right_pipette, water_volume)

        Utilities.dispensing_loop(args, water_loop, water_pipette,
                                  reagent_labware[args.WaterResWell].bottom(water_tip_height),
                                  sample_destination_labware[well], water_volume, NewTip=False, MixReaction=False,
                                  touch=True)

        water_aspirated += water_volume
        water_tip_height = Utilities.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia,
                                                    cone_vol, bottom_offset)

    if left_pipette.has_tip:
        left_pipette.drop_tip()
    if right_pipette.has_tip:
        right_pipette.drop_tip()

    return water_aspirated

    """

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
    cone_vol = Utilities.labware_cone_volume(args, args.ReagentSlot)
    water_tip_height = Utilities.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
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
                Utilities.pipette_selection(left_pipette, right_pipette, sample_vol)

            for well in sample_dest_wells:
                Utilities.dispensing_loop(args, sample_loop, sample_pipette,
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
            Utilities.pipette_selection(left_pipette, right_pipette, diluent_vol)
        sample_pipette, sample_loop, sample_vol = \
            Utilities.pipette_selection(left_pipette, right_pipette, adjusted_sample_vol)

        # Make dilution, diluent first
        dilution_well = dilution_plate_layout[dilution_well_index]

        Utilities.dispensing_loop(args, diluent_loop, diluent_pipette,
                                  reagent_labware[args.WaterResWell].bottom(water_tip_height),
                                  dilution_labware[dilution_well], diluent_vol, NewTip=True, MixReaction=False,
                                  touch=True)
        mix_volume = None
        if diluted_sample_vol < 20:
            mix_volume = 18
        Utilities.dispensing_loop(args, sample_loop, sample_pipette, sample_source_labware[sample_source_well],
                                  dilution_labware[dilution_well], sample_vol, NewTip=True, MixReaction=True,
                                  MixVolume=mix_volume)

        water_aspirated += diluent_vol
        dilution_well_index += 1
        water_tip_height = \
            Utilities.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                     bottom_offset)

        # Add diluted sample to PCR plate
        for well in sample_dest_wells:
            sample_pipette, sample_loop, diluted_sample_vol = \
                Utilities.pipette_selection(left_pipette, right_pipette, diluted_sample_vol)
            mix_volume = None
            if diluted_sample_vol < 20:
                mix_volume = 18
            Utilities.dispensing_loop(args, sample_loop, sample_pipette,
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

class Utilities:

    def labware_cone_volume(args, labware_slot):
        """
        Based on the labware and reservoir return the volume at which the cylinder shape transitions to the conical shape.
        @param args:
        @param labware_slot:
        @return:
        """
        cone_vol = 200

        # labware = getattr(args, "Slot{}".format(str(labware_name)[-1:]))
        labware = getattr(args, "Slot{}".format(str(labware_slot)))

        if "e5ml_" in labware:
            cone_vol = 1200

        elif "e1.5ml_" in labware:
            cone_vol = 450

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

    def load_tipracks(protocol, tiprack_list, labware_dict):
        """
        Creates a list of the pipette tip labware.
        @param protocol:
        @param tiprack_list:
        @param labware_dict:
        @return:
        """
        tiprack_labware = []
        for slot in tiprack_list:
            if slot not in protocol.loaded_labwares:
                tiprack_labware.append(labware_dict[str(slot)])
        return tiprack_labware

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

    def build_labware_dict(protocol, sample_parameters, slot_dict):
        sample_reagent_labware_dict = {}
        for key in sample_parameters:
            sample_slot = sample_parameters[key][0]
            sample_dest_slot = sample_parameters[key][5]

            if sample_dest_slot not in sample_reagent_labware_dict:
                sample_reagent_labware_dict[sample_dest_slot] = \
                    protocol.load_labware(slot_dict[sample_dest_slot], sample_dest_slot)

            if sample_slot not in sample_reagent_labware_dict:
                sample_reagent_labware_dict[sample_slot] = protocol.load_labware(slot_dict[sample_slot], sample_slot)

        return sample_reagent_labware_dict

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

        max_template_vol = round(float(args.PCR_Volume) - float(args.MasterMixPerRxn), ndigits=1)

        # If at least 2 uL of sample is needed then no dilution is necessary
        if template_in_rxn / sample_concentration >= 2:
            sample_vol = round(template_in_rxn / sample_concentration, ndigits=1)
            return sample_vol, 0, 0, max_template_vol - sample_vol, max_template_vol

        # This will test a series of dilutions up to a 1:200.
        for i in range(50):
            dilution = (i + 1) * 2
            diluted_dna_conc = sample_concentration / dilution

            # Want to pipette at least 2 uL of diluted sample per well
            if 2 <= template_in_rxn / diluted_dna_conc <= max_template_vol:
                diluted_sample_vol = round(template_in_rxn / diluted_dna_conc, ndigits=1)
                reaction_water_vol = max_template_vol - diluted_sample_vol

                return 1, dilution - 1, diluted_sample_vol, reaction_water_vol, max_template_vol

    def distribute_reagents(pipette, source_well, destination_wells, dispense_vol):
        """
        Dispense reagents using the distribute function.
        @param pipette:
        @param source_well:
        @param destination_wells:
        @param dispense_vol:
        """
        # ToDo:  Once api 2.15 is out, simplify this.
        p20_default_rate = 7.56
        p300_default_rate = 92.86

        if "P300 Single-Channel GEN2" in str(pipette):
            default_rate = p300_default_rate
        elif "P20 Single-Channel GEN2" in str(pipette):
            default_rate = p20_default_rate

        pipette.flow_rate.aspirate = 30
        pipette.flow_rate.dispense = 10

        pipette.distribute(volume=dispense_vol, source=source_well, dest=destination_wells,
                           touch_tip=True, blow_out=True, radius=0.75, speed=40, disposal_volume=10,
                           blowout_location='source well')

        pipette.flow_rate.aspirate = default_rate
        pipette.flow_rate.dispense = default_rate

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
