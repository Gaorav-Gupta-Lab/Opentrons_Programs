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
# Check if we are on the OT-2, Robotron, or some other computer.
template_parser_path = "{0}var{0}lib{0}jupyter{0}notebooks".format(os.sep)
if not os.path.exists(template_parser_path):
    template_parser_path = "C:{0}Opentrons_Programs".format(os.sep)
    # this is the development computer
    if not os.path.exists(template_parser_path):
        template_parser_path = \
            "C:/Users/dennis/OneDrive - University of North Carolina at Chapel Hill/Projects/Programs/Opentrons_Programs"
sys.path.insert(0, template_parser_path)
import Utilities

# metadata
metadata = {
    'protocolName': 'PCR v1.0.2',
    'author': 'Dennis Simpson <dennis@email.unc.edu>',
    'description': 'Setup a Generic PCR or a ddPCR',
    'apiLevel': '2.13'
}


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


def labware_parsing(args, ctx):
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
    for i in range(len(slot_list)):
        labware = getattr(args, "{}".format(slot_list[i]))
        if labware:
            slot_dict[str(i + 1)] = labware
            labware_dict[str(i + 1)] = ctx.load_labware(labware, str(i + 1))
            if labware in tipbox_dict[args.LeftPipette]:
                left_tiprack_list.append(labware_dict[str(i + 1)])
            elif labware in tipbox_dict[args.RightPipette]:
                right_tiprack_list.append(labware_dict[str(i + 1)])

    return labware_dict, slot_dict, left_tiprack_list, right_tiprack_list


def plate_layout(labware):
    """
    Define the destination layout for the reactions.  Can be 96-well plate or 8-well strip tubes
    @param labware:
    @return:
    """

    layout_data = defaultdict(list)
    for k in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
        layout_data[k] = ['', '', '', '', '', '', '', '', '', '', '', '', ]

    if labware == "stacked_96_well" or labware == "8_well_strip_dilution_tubes":
        column_index = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

    elif labware == "8_well_strip_tubes_200ul":
        column_index = [1, 3, 5, 7, 9, 11, 12]

    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    plate_layout_by_column = []
    for i in column_index:
        for row in rows:
            plate_layout_by_column.append("{}{}".format(row, i))

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

    max_template_vol = round(float(args.PCR_Volume)-float(args.MasterMixPerRxn), 1)

    # If at least 2 uL of sample is needed then no dilution is necessary
    if template_in_rxn/sample_concentration >= 2:
        sample_vol = round(template_in_rxn/sample_concentration, 2)
        return sample_vol, 0, 0, max_template_vol-sample_vol, max_template_vol

    # This will test a series of dilutions up to a 1:200.
    for i in range(50):
        dilution = (i+1)*2
        diluted_dna_conc = sample_concentration/dilution

        # Want to pipette at least 2 uL of diluted sample per well
        if 2 <= template_in_rxn/diluted_dna_conc <= max_template_vol:
            diluted_sample_vol = round(template_in_rxn/diluted_dna_conc, 2)
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

        if args.Template == " Generic PCR" or args.Template == "Generic PCR":
            template_in_rxn = float(sample_parameters[sample_key][6])
        elif args.Template == " ddPCR" or args.Template == "ddPCR":
            template_in_rxn = float(args.DNA_in_Reaction)
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
            [round(sample_vol, 1), round(diluent_vol, 1), round(diluted_sample_vol, 1), sample_wells]

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


def run(ctx: protocol_api.ProtocolContext):

    # TSV file location on OT-2
    tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)

    # If not on the OT-2, get temp TSV file location on Windows Computers for simulation
    if not os.path.isfile(tsv_file_path):
        tsv_file_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())

    sample_parameters, args = parse_sample_template(tsv_file_path)

    ctx.comment("Begin {} {}".format(args.Template, args.Version))

    # Turn on rail lights and pause program so user can load robot deck.
    # ctx.set_rail_lights(True)
    # ctx.pause("Load Labware onto robot deck and click resume when ready to continue")
    # ctx.home()
    ctx.set_rail_lights(False)

    labware_dict, slot_dict, left_tiprack_list, right_tiprack_list = labware_parsing(args, ctx)

    # Pipettes
    left_pipette = ctx.load_instrument(args.LeftPipette, 'left', tip_racks=left_tiprack_list)
    right_pipette = ctx.load_instrument(args.RightPipette, 'right', tip_racks=right_tiprack_list)

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
    if ctx.is_simulating() and platform.system() == "Windows":
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

    # Now do the actual dispensing.
    water_aspirated = dispense_water(args, labware_dict, water_well_dict, left_pipette, right_pipette)

    dispense_reagent_mix(args, labware_dict, target_well_dict, target_info_dict, left_pipette, right_pipette)

    water_aspirated = dispense_samples(args, labware_dict, sample_data_dict, slot_dict, sample_parameters, left_pipette,
                                       right_pipette, water_aspirated)
    if args.Template == " ddPCR":
        fill_empty_wells(args, used_wells, water_aspirated, labware_dict, left_pipette, right_pipette)

    ctx.comment("\nProgram Complete")

    if not ctx.is_simulating():
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
        cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
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
        reagent_well_vol = float(target_info_dict[int(target)][2])
        reagent_aspirated = float(args.MasterMixPerRxn)
        reagent_source_labware = labware_dict[reagent_slot]
        target_well_list = target_well_dict[target]
        reagent_well_dia = reagent_source_labware[reagent_source_well].diameter
        reagent_cone_vol = Utilities.labware_cone_volume(args, reagent_source_labware)
        reagent_pipette, reagent_loop, reagent_volume = \
            Utilities.pipette_selection(left_pipette, right_pipette, reagent_aspirated)

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
    bottom_offset = float(args.BottomOffset)
    cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
    water_res_well_dia = reagent_labware[args.WaterResWell].diameter
    water_tip_height = Utilities.res_tip_height(float(args.WaterResVol), water_res_well_dia, cone_vol, bottom_offset)
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    water_aspirated = 0

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
    water_res_well_dia = reagent_labware[args.WaterResWell].diameter
    cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
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
        diluted_template_needed = round(diluted_sample_vol*(len(sample_dest_wells)+1.5), 2)
        diluted_template_factor = round(diluted_template_needed/(sample_vol+diluent_vol), 2)
        '''
        diluted_template_on_hand = sample_vol+diluent_vol
        diluted_template_factor = 1.0
        if diluted_template_needed <= diluted_template_on_hand:
            diluted_template_factor = diluted_template_needed/diluted_template_on_hand
            if diluted_template_factor <= 1.5 and (sample_vol * diluted_template_factor) < 10:
                diluted_template_factor = 2.0
        '''
        adjusted_sample_vol = round((sample_vol * diluted_template_factor),1)
        diluent_vol = round((diluent_vol*diluted_template_factor), 1)

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
        outfile = open("C:{0}Users{0}{1}{0}Documents{0}Simulation.txt"
                       .format(os.sep, os.getlogin()), 'w', encoding="UTF-16")
        outfile.write(outstring)
        outfile.close()
    protocol_file.close()
