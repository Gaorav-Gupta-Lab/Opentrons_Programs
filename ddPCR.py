import datetime
import os
import sys
import platform
from contextlib import suppress
from collections import defaultdict
from opentrons import protocol_api
from opentrons.simulate import simulate, format_runlog

# Check if we are on the OT-2, Robotron, or some other computer.
if not os.path.exists("C:{0}Opentrons_Programs".format(os.sep)):
    template_parser_path = "{0}var{0}lib{0}jupyter{0}notebooks".format(os.sep)
    if not os.path.exists(template_parser_path):
        template_parser_path = \
            "C:/Users/dennis/OneDrive - University of North Carolina at Chapel Hill/Projects/Programs/Opentrons_Programs"
    sys.path.insert(1, template_parser_path)

from Utilities import parse_sample_template, res_tip_height, labware_cone_volume, labware_parsing, pipette_selection, dispensing_loop

# metadata
metadata = {
    'protocolName': 'ddPCR v0.6.0',
    'author': 'Dennis Simpson',
    'description': 'Setup a ddPCR using either 2x or 4x SuperMix',
    'apiLevel': '2.9'
}


def plate_layout():
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    plate_layout_by_column = []
    for i in range(12):
        for row in rows:
            plate_layout_by_column.append("{}{}".format(row, i+1))
    return plate_layout_by_column


def calculate_volumes(args, sample_concentration):
    """
    Calculates volumes for dilution and distribution of sample.
    Returns a list of tuples consisting of
    (uL of sample to dilute, uL of water for dilution), (uL of diluted sample in reaction, uL of water in reaction)
    @param args:
    @param sample_concentration:
    @return:
    """

    template_in_reaction = float(args.DNA_in_Reaction)
    # Target Concentration is used to keep volumes > 1 uL
    target_concentration = template_in_reaction/2
    max_template_vol = \
        round(float(args.PCR_Volume)-((float(args.PCR_Volume) / int(args.PCR_MixConcentration.split("x")[0])) + float(args.TargetVolume)), 2)

    min_dna_in_reaction = template_in_reaction/max_template_vol

    # If template concentration per uL is less than desired template in reaction then no dilution is necessary.
    if sample_concentration <= target_concentration:
        sample_vol = template_in_reaction/sample_concentration
        return 0, 0, sample_vol, round(max_template_vol-sample_vol, 2), max_template_vol

    # This will test a series of dilutions up to a 1:200.
    for i in range(50):
        dilution = (i+1)*2
        diluted_dna_conc = sample_concentration/dilution

        if target_concentration >= diluted_dna_conc >= min_dna_in_reaction:
            dilution_data = (1, dilution - 1)
            diluted_sample_vol = round(template_in_reaction/diluted_dna_conc, 2)
            reaction_water_vol = max_template_vol-diluted_sample_vol

            return dilution_data[0], dilution_data[1], diluted_sample_vol, reaction_water_vol, max_template_vol


def sample_processing(args, sample_parameters):
    sample_data_dict = defaultdict(list)
    target_well_dict = defaultdict(list)
    water_well_dict = defaultdict(float)
    layout_data = defaultdict(list)

    # Builds the data frame for printing the plate layout file
    for k in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
        layout_data[k] = ['', '', '', '', '', '', '', '', '', '', '', '', ]

    plate_layout_by_column = plate_layout()
    used_wells = []
    dest_well_count = 0

    for sample_key in sample_parameters:
        sample_name = sample_parameters[sample_key][2]
        sample_string = sample_name
        concentration = float(sample_parameters[sample_key][3])
        targets = sample_parameters[sample_key][4].split(",")
        replicates = int(sample_parameters[sample_key][5])
        sample_vol, diluent_vol, diluted_sample_vol, reaction_water_vol, max_template_vol = \
            calculate_volumes(args, concentration)
        diluted_sample_vol = round(diluted_sample_vol, 2)
        sample_wells = []
        for target in targets:
            for i in range(replicates):
                well = plate_layout_by_column[dest_well_count]
                row = well[0]
                column = int(well[1:])-1

                try:
                    dilution = "1:{}".format(int((sample_vol+diluent_vol)/sample_vol))
                except ZeroDivisionError:
                    dilution = "Neat"

                layout_data[row][column] = "{}|{}|{}|{}"\
                    .format(sample_string, getattr(args, "Target_{}".format(target)).split(",")[2],
                            dilution, diluted_sample_vol)

                water_well_dict[well] = reaction_water_vol
                target_well_dict[target].append(well)
                sample_wells.append(well)
                used_wells.append(well)
                dest_well_count += 1
                
        sample_data_dict[sample_key] = [sample_vol, diluent_vol, diluted_sample_vol, sample_wells]

    # Define our positive control wells for the targets.
    for target in target_well_dict:
        loop_count = 2
        template = getattr(args, "PositiveControl_{}".format(target)).split(",")[2]
        target_name = getattr(args, "Target_{}".format(target)).split(",")[2]

        while loop_count > 0:
            well = plate_layout_by_column[dest_well_count]
            used_wells.append(well)
            target_well_dict[target].append(well)
            row = well[0]
            column = int(well[1:])-1
            layout_data[row][column] = "{}|{}|NA|{}".format(template, target_name, max_template_vol)

            if template == "Water":
                water_well_dict[well] = max_template_vol

            template = "Water"
            dest_well_count += 1
            loop_count -= 1

    return sample_data_dict, water_well_dict, target_well_dict, used_wells, layout_data, max_template_vol


def run(ctx: protocol_api.ProtocolContext):

    ctx.comment("Begin {}".format(metadata['protocolName']))

    # Turn on rail lights and pause program so user can load robot deck.
    ctx.set_rail_lights(True)
    ctx.pause("Load Labware onto robot deck and click resume when ready to continue")
    ctx.home()
    ctx.set_rail_lights(False)

    # TSV file location on OT-2
    tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)

    # If not on the OT-2, get temp TSV file location on Win10 Computers for simulation
    if not os.path.isfile(tsv_file_path):
        tsv_file_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())

    sample_parameters, args = parse_sample_template(tsv_file_path)
    labware_dict, left_tiprack_list, right_tiprack_list = labware_parsing(args, ctx)

    # Pipettes
    left_pipette = ctx.load_instrument(args.LeftPipette, 'left', tip_racks=left_tiprack_list)
    right_pipette = ctx.load_instrument(args.RightPipette, 'right', tip_racks=right_tiprack_list)

    # Set the location of the first tip in box.
    with suppress(IndexError):
        left_pipette.starting_tip = left_tiprack_list[0].wells_by_name()[args.LeftPipetteFirstTip]
    with suppress(IndexError):
        right_pipette.starting_tip = right_tiprack_list[0].wells_by_name()[args.RightPipetteFirstTip]

    sample_data_dict, water_well_dict, target_well_dict, used_wells, layout_data, max_template_vol = \
        sample_processing(args, sample_parameters)

    # This will output a plate layout file.  Only does it during the simulation from our GUI
    if ctx.is_simulating() and platform.system() == "Windows":
        # if not os.path.isdir("C:{0}Users{0}{1}{0}Documents".format(os.sep, os.getlogin())):
        run_date = datetime.datetime.today().strftime("%a %b %d %H:%M %Y")
        plate_layout_string = \
            "## ddPCR Setup\n## Setup Date:\t{}\n## Template User:\t{}\n" \
            "# Format:\tTemplate | Target | Template Dilution | Template Volume in Reaction\n\n\t"\
            .format(run_date, args.User)

        for i in range(12):
            plate_layout_string += "{}\t".format(i+1)

        # I have to import this here because I have been unable to get natsort on the robot.
        import natsort

        for well in natsort.natsorted(layout_data):
            well_string = "\t".join(layout_data[well])
            plate_layout_string += "\n{}\t{}\t".format(well, well_string)
        plate_layout_file = open("C:{0}Users{0}{1}{0}Documents{0}PlateLayout.tsv".format(os.sep, os.getlogin()), 'w')
        plate_layout_file.write(plate_layout_string)
        plate_layout_file.close()

    # Now do the actual dispensing.
    water_aspirated = dispense_water(args, labware_dict, water_well_dict, left_pipette, right_pipette)

    positive_control_dict = dispense_primer_mix(args, labware_dict, target_well_dict, left_pipette, right_pipette)

    water_aspirated = dispense_samples(args, labware_dict, sample_data_dict, sample_parameters, left_pipette,
                                       right_pipette, water_aspirated)

    dispense_positive_controls(args, positive_control_dict, labware_dict, left_pipette, right_pipette, max_template_vol)

    dispense_supermix(args, labware_dict, left_pipette, right_pipette, used_wells)

    fill_empty_wells(args, used_wells, water_aspirated, labware_dict, left_pipette, right_pipette)

    ctx.comment("\nProgram Complete")

    if not ctx.is_simulating():
        os.remove(tsv_file_path)


def fill_empty_wells(args, used_wells, water_aspirated, labware_dict, left_pipette, right_pipette):
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
        water_res_well_dia = reagent_labware[args.WaterWell].diameter
        cone_vol = labware_cone_volume(args, reagent_labware)
        fill_pipette, fill_loop, fill_vol = pipette_selection(left_pipette, right_pipette, float(args.PCR_Volume))
        water_tip_height = res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                          float(args.BottomOffset))

        for i in range(wells_remaining):
            blank_well = "{}{}".format(row_list[i+row_index+1], column)
            dispensing_loop(args, fill_loop, fill_pipette,
                            reagent_labware[args.WaterWell].bottom(water_tip_height),
                            sample_destination_labware[blank_well], fill_vol,
                            NewTip=False, MixReaction=False)
            water_aspirated = water_aspirated+fill_vol
            water_tip_height = res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                              float(args.BottomOffset))
        fill_pipette.drop_tip()


def dispense_primer_mix(args, labware_dict, target_well_dict, left_pipette, right_pipette):
    """

    @param args:
    @param labware_dict:
    @param target_well_dict:
    @param left_pipette:
    @param right_pipette:
    @return:
    """
    target_pipette, target_loop, target_volume = \
        pipette_selection(left_pipette, right_pipette, float(args.TargetVolume))
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    positive_control_dict = defaultdict(list)

    # Dispense target primers into all wells
    for target in target_well_dict:
        target_info = getattr(args, "Target_{}".format(target)).split(",")
        target_slot = target_info[0]
        target_source_well = target_info[1]
        target_source_labware = labware_dict[target_slot]
        target_well_list = target_well_dict[target]

        if not target_pipette.has_tip:
            target_pipette.pick_up_tip()

        # The penultimate well in the list is the positive control for the target primers
        pos_control_well = target_well_list[-2]
        pos_control_info = getattr(args, "PositiveControl_{}".format(target))
        positive_control_dict[pos_control_well] = pos_control_info

        for well in target_well_list:
            dispensing_loop(args, target_loop, target_pipette, target_source_labware[target_source_well],
                            sample_destination_labware[well], target_volume, NewTip=False, MixReaction=False,
                            touch=True)

        # Drop any tips the pipettes might have.
        if left_pipette.has_tip:
            left_pipette.drop_tip()
        if right_pipette.has_tip:
            right_pipette.drop_tip()

    return positive_control_dict


def dispense_water(args, labware_dict, water_well_dict, left_pipette, right_pipette):
    """

    @param args:
    @param labware_dict:
    @param water_well_dict:
    @param left_pipette:
    @param right_pipette:
    @return:
    """
    reagent_labware = labware_dict[args.ReagentSlot]
    cone_vol = labware_cone_volume(args, reagent_labware)
    water_res_well_dia = reagent_labware[args.WaterWell].diameter
    water_tip_height = res_tip_height(float(args.WaterResVol), water_res_well_dia, cone_vol, float(args.BottomOffset))
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    water_aspirated = 0

    for well in water_well_dict:
        water_volume = water_well_dict[well]

        # There is a bug in the Opentrons software that results in a 0 uL aspiration defaulting to pipette max.
        if water_volume < 0.2:
            continue

        # Define the pipette for dispensing the water.
        water_pipette, water_loop, water_volume = pipette_selection(left_pipette, right_pipette, water_volume)

        if not water_pipette.has_tip:
            water_pipette.pick_up_tip()

        dispensing_loop(args, water_loop, water_pipette, reagent_labware[args.WaterWell].bottom(water_tip_height),
                        sample_destination_labware[well], water_volume, NewTip=False, MixReaction=False)
        water_aspirated += water_volume
        water_tip_height = res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                          float(args.BottomOffset))

    if left_pipette.has_tip:
        left_pipette.drop_tip()
    if right_pipette.has_tip:
        right_pipette.drop_tip()

    return water_aspirated


def dispense_samples(args, labware_dict, sample_data_dict, sample_parameters, left_pipette, right_pipette,
                     water_aspirated):
    """
    Dilute and dispense samples
    @param args:
    @param labware_dict:
    @param sample_data_dict:
    @param sample_parameters:
    @param left_pipette:
    @param right_pipette:
    @param water_aspirated:
    """
    dilution_labware = labware_dict[args.DilutionPlateSlot]
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    reagent_labware = labware_dict[args.ReagentSlot]
    water_res_well_dia = reagent_labware[args.WaterWell].diameter
    cone_vol = labware_cone_volume(args, reagent_labware)
    water_tip_height = res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol)
    dilution_plate_layout = plate_layout()
    dilution_well_index = 0

    for sample_key in sample_parameters:
        sample_source_labware = labware_dict[sample_parameters[sample_key][0]]
        sample_source_well = sample_parameters[sample_key][1]
        sample_dest_wells = sample_data_dict[sample_key][3]
        undiluted_sample_vol = sample_data_dict[sample_key][0]
        diluent_vol = sample_data_dict[sample_key][1]
        diluted_sample_vol = sample_data_dict[sample_key][2]

        # If no dilution is necessary, dispense sample and continue
        if undiluted_sample_vol == 0:
            sample_pipette, sample_loop, undiluted_sample_vol = \
                pipette_selection(left_pipette, right_pipette, diluted_sample_vol)
            for well in sample_dest_wells:
                dispensing_loop(args, sample_loop, sample_pipette,
                                sample_source_labware[sample_source_well],
                                sample_destination_labware[well], undiluted_sample_vol,
                                NewTip=True, MixReaction=False)
            continue

        # Adjust volume of diluted sample to make sure there is enough
        volume_ratio = (diluted_sample_vol*len(sample_dest_wells)) / (undiluted_sample_vol+diluent_vol)
        if volume_ratio >= 0.66:
            for i in range(len(sample_dest_wells)+2):
                undiluted_sample_vol = undiluted_sample_vol*(i+1)
                diluent_vol = diluent_vol*(i+1)
                volume_ratio = (diluted_sample_vol * len(sample_dest_wells)) / (undiluted_sample_vol + diluent_vol)
                if volume_ratio < 0.66:
                    break

        # Reset the pipettes for the new volumes
        diluent_pipette, diluent_loop, diluent_vol = pipette_selection(left_pipette, right_pipette, diluent_vol)
        sample_pipette, sample_loop, undiluted_sample_vol = \
            pipette_selection(left_pipette, right_pipette, undiluted_sample_vol)

        # Make dilution, diluent first
        dilution_well = dilution_plate_layout[dilution_well_index]
        dispensing_loop(args, diluent_loop, diluent_pipette, reagent_labware[args.WaterWell].bottom(water_tip_height),
                        dilution_labware[dilution_well], diluent_vol, NewTip=True, MixReaction=False)
        dispensing_loop(args, sample_loop, sample_pipette, sample_source_labware[sample_source_well],
                        dilution_labware[dilution_well], undiluted_sample_vol, NewTip=True, MixReaction=True)
        water_aspirated += diluent_vol
        dilution_well_index += 1
        water_tip_height = res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                          float(args.BottomOffset))

        # Add diluted sample to PCR plate
        for well in sample_dest_wells:
            sample_pipette, sample_loop, diluted_sample_vol = \
                pipette_selection(left_pipette, right_pipette, diluted_sample_vol)

            dispensing_loop(args, sample_loop, sample_pipette, dilution_labware[dilution_well].bottom(0.1),
                            sample_destination_labware[well], diluted_sample_vol, NewTip=True, MixReaction=False)
        if sample_pipette.has_tip:
            sample_pipette.drop_tip()
    return water_aspirated


def dispense_positive_controls(args, positive_control_dict, labware_dict, left_pipette, right_pipette,
                               max_template_vol):
    """
    Dispense positive controls to appropriate wells
    @param args:
    @param positive_control_dict:
    @param labware_dict:
    @param left_pipette:
    @param right_pipette:
    @param max_template_vol:
    """

    sample_destination_labware = labware_dict[args.PCR_PlateSlot]

    for destination_well in positive_control_dict:
        positive_control_info = positive_control_dict[destination_well].split(",")
        control_sample_vol = max_template_vol

        # There is a bug in the Opentrons software that results in a 0 uL aspiration defaulting to pipette max.
        if control_sample_vol < 0.4:
            continue

        positive_control_slot = positive_control_info[0]
        positive_control_source_well = positive_control_info[1]
        positive_control_source_labware = labware_dict[positive_control_slot]
        control_pipette, control_loop, control_sample_vol = \
            pipette_selection(left_pipette, right_pipette, control_sample_vol)

        if not control_pipette.has_tip:
            control_pipette.pick_up_tip()

        dispensing_loop(args, control_loop, control_pipette,
                        positive_control_source_labware[positive_control_source_well],
                        sample_destination_labware[destination_well], control_sample_vol,
                        NewTip=True, MixReaction=False)


def dispense_supermix(args, labware_dict, left_pipette, right_pipette, used_wells):
    """
    Dispense the SuperMix into each well.
    @param args:
    @param labware_dict:
    @param left_pipette:
    @param right_pipette:
    @param used_wells:
    """

    reagent_labware = labware_dict[args.ReagentSlot]
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    supermix_vol = float(args.PCR_Volume) / int(args.PCR_MixConcentration.split("x")[0])
    supermix_pipette, supermix_loop, supermix_vol = pipette_selection(left_pipette, right_pipette, supermix_vol)
    supermix_source_well = args.PCR_MixWell
    supermix_well_dia = reagent_labware[supermix_source_well].diameter
    cone_vol = labware_cone_volume(args, reagent_labware)
    supermix_tip_height = res_tip_height(float(args.PCR_MixResVolume), supermix_well_dia, cone_vol,
                                         float(args.BottomOffset))
    supermix_aspirated = 0

    for well in used_wells:
        dispensing_loop(args, supermix_loop, supermix_pipette,
                        reagent_labware[supermix_source_well].bottom(supermix_tip_height),
                        sample_destination_labware[well], supermix_vol,
                        NewTip=True, MixReaction=True)

        supermix_aspirated += supermix_vol
        supermix_tip_height = res_tip_height(float(args.PCR_MixResVolume)-supermix_aspirated, supermix_well_dia,
                                             cone_vol, float(args.BottomOffset))


if __name__ == "__main__":
    protocol_file = open('ddPCR.py')
    labware_path = "{}{}custom_labware".format(os.getcwd(), os.sep)
    run_log, __bundle__ = simulate(protocol_file, custom_labware_paths=[labware_path])
    run_date = datetime.datetime.today().strftime("%a %b %d %H:%M %Y")
    i = 1
    t = format_runlog(run_log).split("\n")

    outstring = "Opentrons OT-2 Steps for {}.\nDate:  {}\nProgram File: ddPCR.py\n\nStep\tCommand\n" \
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
