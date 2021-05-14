import datetime
import os
import sys
import platform
from collections import defaultdict
from contextlib import suppress
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
    'protocolName': 'Generic PCR v0.7.0',
    'author': 'Dennis Simpson',
    'description': 'Sets up a PCR from concentrated template',
    'apiLevel': '2.9'
}


def run(ctx: protocol_api.ProtocolContext):

    # Turn on rail lights and pause program so user can load robot deck.
    ctx.set_rail_lights(True)
    ctx.pause("Load Labware onto robot deck and click resume when ready to continue")
    ctx.home()
    ctx.set_rail_lights(False)

    # TSV file location on OT-2
    tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)
    if not os.path.isfile(tsv_file_path):
        # Temp TSV file location on Win10 Computers for simulation
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

    # Dispense Samples
    sample_data_dict, aspirated_water_vol = \
        dispense_samples(args, sample_parameters, labware_dict, left_pipette, right_pipette)

    # Dispense PCR Reagents.
    dispense_pcr_reagents(args, labware_dict, left_pipette, right_pipette, aspirated_water_vol, sample_data_dict)

    if not ctx.is_simulating():
        os.remove(tsv_file_path)


def dispense_pcr_reagents(args, labware_dict, left_pipette, right_pipette, aspirated_water_vol, sample_data_dict):
    """
    Dispense PCR reagents into each well and setup negative control well.
    @param args:
    @param labware_dict:
    @param left_pipette:
    @param right_pipette:
    @param aspirated_water_vol:
    @param sample_data_dict:
    """

    # Define PCR reagent pipette and reagent volume.
    reagent_pipette, reagent_loop, pcr_reagent_vol = \
        pipette_selection(left_pipette, right_pipette, float(args.PCR_Volume) * 0.5)

    reagent_labware = labware_dict[args.ReagentSlot]
    reagent_source_well = args.PCR_MixWell
    pcr_reagent_well_dia = reagent_labware[args.PCR_MixWell].diameter
    cone_vol = labware_cone_volume(args, reagent_labware)
    tip_height = res_tip_height(float(args.PCR_MixResVolume), pcr_reagent_well_dia, cone_vol, float(args.BottomOffset))
    reagent_aspirated = 0

    for sample_key in sample_data_dict:
        sample_data = sample_data_dict[sample_key][0]
        reagent_dest_labware = labware_dict[sample_data[0]]
        reagent_dest_wells = sample_data[1]
        for reagent_dest_well in reagent_dest_wells:
            dispensing_loop(args, reagent_loop, reagent_pipette,
                            reagent_labware[reagent_source_well].bottom(tip_height),
                            reagent_dest_labware[reagent_dest_well], pcr_reagent_vol, NewTip=True, MixReaction=True,
                            touch=True)

            reagent_aspirated += pcr_reagent_vol
            tip_height = \
                res_tip_height(float(args.PCR_MixResVolume) - reagent_aspirated, pcr_reagent_well_dia, cone_vol,
                               float(args.BottomOffset))

    # Setup the water control sample(s)
    if args.WaterControl:
        water_reservoir_dia = reagent_labware[args.WaterWell].diameter
        water_tip_height = \
            res_tip_height(float(args.WaterResVol) - aspirated_water_vol, water_reservoir_dia, cone_vol,
                           float(args.BottomOffset))

        slot = args.WaterControl.split(',')[0]
        well = args.WaterControl.split(',')[1]
        neg_control_labware = labware_dict[slot]
        reagent_pipette, reagent_loop, pcr_reagent_vol = \
            pipette_selection(left_pipette, right_pipette, float(args.PCR_Volume) * 0.5)
        reagent_pipette.pick_up_tip()
        temp_loop = reagent_loop

        while temp_loop > 0:
            # Dispense Water into negative control well
            reagent_pipette.aspirate(pcr_reagent_vol, reagent_labware[args.WaterWell].bottom(water_tip_height))
            reagent_pipette.dispense(pcr_reagent_vol, neg_control_labware[well])
            reagent_pipette.blow_out()
            temp_loop -= 1

        while reagent_loop > 0:
            # Dispense PCR reagents into negative control well
            reagent_pipette.aspirate(pcr_reagent_vol, reagent_labware[args.PCR_MixWell])
            reagent_pipette.dispense(pcr_reagent_vol, neg_control_labware[well])
            reagent_loop -= 1

        reagent_pipette.mix(repetitions=4, volume=float(args.PCR_Volume)*0.7, rate=5.0)
        reagent_pipette.blow_out()
        reagent_pipette.touch_tip(radius=0.75, v_offset=-8)
        reagent_pipette.drop_tip()


def dispense_samples(args, sample_parameters, labware_dict, left_pipette, right_pipette):
    """
    Add water and template to the destination wells for the PCR.
    @param args:
    @param sample_parameters:
    @param labware_dict:
    @param left_pipette:
    @param right_pipette:
    @return:
    """

    sample_mass = float(args.DNA_in_Reaction)
    reaction_vol = float(args.PCR_Volume)
    reagent_labware = labware_dict[args.ReagentSlot]
    water_reservoir_dia = reagent_labware[args.WaterWell].diameter
    cone_vol = labware_cone_volume(args, reagent_labware)
    tip_height = res_tip_height(float(args.WaterResVol), water_reservoir_dia, cone_vol, float(args.BottomOffset))
    aspirated_water_vol = 0

    # Figure out volumes and dispense water into each destination well.
    sample_data_dict = defaultdict(list)
    for sample_key in sample_parameters:
        concentration = float(sample_parameters[sample_key][3])
        sample_volume = round(sample_mass/concentration, 2)
        sample_dest_slot = sample_parameters[sample_key][4]
        sample_destination_labware = labware_dict[sample_dest_slot]
        sample_source_labware = labware_dict[sample_parameters[sample_key][0]]
        sample_source_well = sample_parameters[sample_key][1]

        # If user has replicates there can be more than one sample destination well.
        sample_dest_wells = sample_parameters[sample_key][5].split(",")

        water_volume = round((0.5*reaction_vol)-sample_volume, 2)
        sample_data_dict[sample_key].append((sample_dest_slot, sample_dest_wells, sample_volume))

        # Define the pipettes for dispensing the water.
        water_pipette, water_loop, water_volume = pipette_selection(left_pipette, right_pipette, water_volume)
        sample_pipette, sample_loop, sample_volume = pipette_selection(left_pipette, right_pipette, sample_volume)

        # Add water to all the destination wells for this sample.
        for sample_dest_well in sample_dest_wells:
            if not water_pipette.has_tip:
                water_pipette.pick_up_tip()

            dispensing_loop(args, water_loop, water_pipette, reagent_labware[args.WaterWell].bottom(tip_height),
                            sample_destination_labware[sample_dest_well], water_volume, NewTip=False, MixReaction=False)

            aspirated_water_vol += water_volume
            tip_height = res_tip_height(float(args.WaterResVol) - aspirated_water_vol, water_reservoir_dia, cone_vol,
                                        float(args.BottomOffset))

        # Add template to all the destination wells for this sample.
        for sample_dest_well in sample_dest_wells:

            if not sample_pipette.has_tip:
                sample_pipette.pick_up_tip()

            dispensing_loop(args, sample_loop, sample_pipette, sample_source_labware[sample_source_well],
                            sample_destination_labware[sample_dest_well], sample_volume, NewTip=False,
                            MixReaction=False)

        # Drop any tips the pipettes might have.
        if left_pipette.has_tip:
            left_pipette.drop_tip()
        if right_pipette.has_tip:
            right_pipette.drop_tip()

    return sample_data_dict, aspirated_water_vol


if __name__ == "__main__":
    protocol_file = open('Generic_PCR.py')
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
