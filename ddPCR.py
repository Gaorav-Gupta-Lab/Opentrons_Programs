import datetime
import os
import sys
import platform
from contextlib import suppress
from collections import defaultdict
from opentrons import protocol_api
from opentrons.simulate import simulate, format_runlog
# Check if we are on the OT-2, Robotron, or some other computer.

template_parser_path = "{0}var{0}lib{0}jupyter{0}notebooks".format(os.sep)
if not os.path.exists(template_parser_path):
    template_parser_path = "C:{0}Opentrons_Programs".format(os.sep)
    if not os.path.exists(template_parser_path):
        template_parser_path = \
            "C:/Users/dennis/OneDrive - University of North Carolina at Chapel Hill/Projects/Programs/Opentrons_Programs"
sys.path.insert(0, template_parser_path)
import Utilities

# metadata
metadata = {
    'protocolName': 'ddPCR v0.9.6',
    'author': 'Dennis Simpson',
    'description': 'Setup a ddPCR using either 2x or 4x SuperMix',
    'apiLevel': '2.11'
}


def sample_processing(args, sample_parameters):
    sample_data_dict = defaultdict(list)
    target_well_dict = defaultdict(list)
    water_well_dict = defaultdict(float)
    layout_data = defaultdict(list)

    # Builds the data frame for printing the plate layout file
    for k in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
        layout_data[k] = ['', '', '', '', '', '', '', '', '', '', '', '', ]

    plate_layout_by_column = Utilities.plate_layout()
    used_wells = []
    dest_well_count = 0
    target_list = []

    for sample_key in sample_parameters:
        sample_name = sample_parameters[sample_key][2]
        sample_string = sample_name
        sample_concentration = float(sample_parameters[sample_key][3])
        targets = sample_parameters[sample_key][4].split(",")
        replicates = int(sample_parameters[sample_key][5])

        sample_vol, diluent_vol, diluted_sample_vol, reaction_water_vol, max_template_vol = \
            Utilities.calculate_volumes(args, sample_concentration)

        sample_wells = []
        for target in targets:
            target_list.append(target)
            target_name = getattr(args, "Target_{}".format(target))[1]
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
                    .format(sample_string, target_name, dilution, s_volume)

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
        target_data = getattr(args, "Target_{}".format(target))
        target_name = target_data[1]

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

    ctx.comment("Begin {}".format(metadata['protocolName']))

    # Turn on rail lights and pause program so user can load robot deck.
    # ctx.set_rail_lights(True)
    # ctx.pause("Load Labware onto robot deck and click resume when ready to continue")
    # ctx.home()
    ctx.set_rail_lights(False)

    # TSV file location on OT-2
    tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)

    # If not on the OT-2, get temp TSV file location on Win10 Computers for simulation
    if not os.path.isfile(tsv_file_path):
        tsv_file_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())

    sample_parameters, args = Utilities.parse_sample_template(tsv_file_path)
    labware_dict, left_tiprack_list, right_tiprack_list = Utilities.labware_parsing(args, ctx)

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
        plate_layout_file = \
            open("C:{0}Users{0}{1}{0}Documents{0}ddPCR_PlateLayout.tsv".format(os.sep, os.getlogin()), 'w')
        plate_layout_file.write(plate_layout_string)
        plate_layout_file.close()

    # Now do the actual dispensing.
    water_aspirated = dispense_water(args, labware_dict, water_well_dict, left_pipette, right_pipette)

    dispense_reagent_mix(args, labware_dict, target_well_dict, left_pipette, right_pipette)

    water_aspirated = dispense_samples(args, labware_dict, sample_data_dict, sample_parameters, left_pipette,
                                       right_pipette, water_aspirated)

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
        water_res_well_dia = reagent_labware[args.WaterWell].diameter
        cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
        fill_pipette, fill_loop, fill_vol = \
            Utilities.pipette_selection(left_pipette, right_pipette, float(args.PCR_Volume))
        water_tip_height = \
            Utilities.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                     bottom_offset)

        for i in range(wells_remaining):
            blank_well = "{}{}".format(row_list[i+row_index+1], column)
            Utilities.dispensing_loop(args, fill_loop, fill_pipette,
                                      reagent_labware[args.WaterWell].bottom(water_tip_height),
                                      sample_destination_labware[blank_well], fill_vol,
                                      NewTip=False, MixReaction=False)
            water_aspirated = water_aspirated+fill_vol
            water_tip_height = Utilities.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia,
                                                        cone_vol, bottom_offset)
        fill_pipette.drop_tip()


def dispense_reagent_mix(args, labware_dict, target_well_dict, left_pipette, right_pipette):
    """
    This will dispense our primer+probe+BSA+Supermix reagent mixture into each well and mix it with the template.

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
        target_info = getattr(args, "Target_{}".format(target))
        reagent_slot = args.ReagentSlot
        reagent_source_well = target_info[0]
        reagent_well_vol = float(target_info[2])
        reagent_aspirated = float(args.ReagentVolume)
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
                                      NewTip=False, MixReaction=False, touch=True, MixVolume=None, speed=0.75)

            reagent_aspirated += reagent_volume

        # Drop any tips the pipettes might have.
        if left_pipette.has_tip:
            left_pipette.drop_tip()
        if right_pipette.has_tip:
            right_pipette.drop_tip()

    return


def dispense_water(args, labware_dict, water_well_dict, left_pipette, right_pipette):
    """
    This will dispense water into any wells that require it in the PCR destination plate only.  Water for dilutions
    is dispensed as needed.

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
    water_res_well_dia = reagent_labware[args.WaterWell].diameter
    water_tip_height = Utilities.res_tip_height(float(args.WaterResVol), water_res_well_dia, cone_vol, bottom_offset)
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    water_aspirated = 0

    for well in water_well_dict:
        water_volume = water_well_dict[well]

        # There is a bug in the Opentrons software that results in a 0 uL aspiration defaulting to pipette max.
        if water_volume < 0.2:
            continue

        # Define the pipette for dispensing the water.
        water_pipette, water_loop, water_volume = \
            Utilities.pipette_selection(left_pipette, right_pipette, water_volume)

        Utilities.dispensing_loop(args, water_loop, water_pipette,
                                  reagent_labware[args.WaterWell].bottom(water_tip_height),
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

    bottom_offset = float(args.BottomOffset)
    dilution_labware = labware_dict[args.DilutionPlateSlot]
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    reagent_labware = labware_dict[args.ReagentSlot]
    water_res_well_dia = reagent_labware[args.WaterWell].diameter
    cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
    water_tip_height = Utilities.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                                bottom_offset)
    dilution_plate_layout = Utilities.plate_layout()
    dilution_well_index = 0

    for sample_key in sample_parameters:
        sample_source_labware = labware_dict[sample_parameters[sample_key][0]]
        sample_source_well = sample_parameters[sample_key][1]
        sample_dest_wells = sample_data_dict[sample_key][3]
        sample_vol = sample_data_dict[sample_key][0]
        diluent_vol = sample_data_dict[sample_key][1]
        diluted_sample_vol = sample_data_dict[sample_key][2]

        # If no dilution is necessary, dispense sample and continue
        if diluted_sample_vol == 0:
            sample_pipette, sample_loop, sample_vol = \
                Utilities.pipette_selection(left_pipette, right_pipette, sample_vol)
            for well in sample_dest_wells:
                Utilities.dispensing_loop(args, sample_loop, sample_pipette,
                                          sample_source_labware[sample_source_well],
                                          sample_destination_labware[well], sample_vol,
                                          NewTip=True, MixReaction=True, touch=True)
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
                                  reagent_labware[args.WaterWell].bottom(water_tip_height),
                                  dilution_labware[dilution_well], diluent_vol, NewTip=True, MixReaction=False,
                                  touch=True)

        Utilities.dispensing_loop(args, sample_loop, sample_pipette, sample_source_labware[sample_source_well],
                                  dilution_labware[dilution_well], sample_vol, NewTip=True, MixReaction=True)
        water_aspirated += diluent_vol
        dilution_well_index += 1
        water_tip_height = \
            Utilities.res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol,
                                     bottom_offset)

        # Add diluted sample to PCR plate
        for well in sample_dest_wells:
            sample_pipette, sample_loop, diluted_sample_vol = \
                Utilities.pipette_selection(left_pipette, right_pipette, diluted_sample_vol)

            Utilities.dispensing_loop(args, sample_loop, sample_pipette,
                                      dilution_labware[dilution_well].bottom(bottom_offset),
                                      sample_destination_labware[well], diluted_sample_vol, NewTip=True,
                                      MixReaction=True)

        if sample_pipette.has_tip:
            sample_pipette.drop_tip()

    return water_aspirated


if __name__ == "__main__":
    protocol_file = open('ddPCR_v2.py')
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
