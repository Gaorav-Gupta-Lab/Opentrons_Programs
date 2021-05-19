import datetime
import os
import platform
import sys
from collections import defaultdict
from contextlib import suppress
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
from Utilities import parse_sample_template, res_tip_height, labware_cone_volume, labware_parsing, pipette_selection, dispensing_loop

# metadata
metadata = {
    'protocolName': 'Illumina Dual Indexing v0.5.1',
    'author': 'Dennis Simpson',
    'description': 'Add Illumina dual indexing to library',
    'apiLevel': '2.9'
}


def add_pcr_mix(args, labware_dict, sample_dest_dict, left_pipette, right_pipette):
    pcr_reagent_vol = float(args.PCR_Volume) * 0.5
    pcr_reagent_labware = labware_dict[args.ReagentSlot]
    pcr_reagent_well = args.PCR_MixWell
    reagent_labware = labware_dict[args.ReagentSlot]
    pcr_reservoir_dia = reagent_labware[args.PCR_MixWell].diameter
    cone_vol = labware_cone_volume(args, reagent_labware)
    tip_height = res_tip_height(float(args.PCR_MixResVolume), pcr_reservoir_dia, cone_vol, float(args.BottomOffset))
    pcr_pipette, pcr_loop, pcr_reagent_vol = pipette_selection(left_pipette, right_pipette, pcr_reagent_vol)
    aspirated_vol = 0

    for sample_key in sample_dest_dict:
        sample_dest_slot = sample_dest_dict[sample_key][0][0]
        sample_dest_well = sample_dest_dict[sample_key][0][1]
        sample_destination_labware = labware_dict[sample_dest_slot]

        pcr_pipette = \
            dispensing_loop(args, pcr_loop, pcr_pipette, pcr_reagent_labware[pcr_reagent_well].bottom(tip_height),
                            sample_destination_labware[sample_dest_well], pcr_reagent_vol,
                            NewTip=True, MixReaction=True)

        aspirated_vol += pcr_reagent_vol

        tip_height = res_tip_height(float(args.PCR_MixResVolume)-aspirated_vol, pcr_reservoir_dia, cone_vol,
                                    float(args.BottomOffset))


def dispense_samples(args, sample_parameters, labware_dict, left_pipette, right_pipette):
    """
    Add water, template, and indexing primers to the destination wells for the PCR.
    @param args:
    @param sample_parameters:
    @param labware_dict:
    @param left_pipette:
    @param right_pipette:
    @return:
    """
    # Extract Index Primer information
    primer_list = \
        ["D501", "D502", "D503", "D504", "D505", "D506", "D507", "D508", "D701", "D702", "D703", "D704", "D705",
         "D706", "D707", "D708", "D709", "D710a", "D711", "D712"]

    primer_dict = {}

    for i in range(len(primer_list)):
        primer_well = getattr(args, "{}".format(primer_list[i]))
        primer_dict[primer_list[i]] = primer_well

    sample_mass = float(args.DNA_in_Reaction)
    reaction_vol = float(args.PCR_Volume)
    reagent_labware = labware_dict[args.ReagentSlot]
    water_reservoir_dia = reagent_labware[args.WaterWell].diameter
    cone_vol = labware_cone_volume(args, reagent_labware)
    tip_height = res_tip_height(float(args.WaterResVol), water_reservoir_dia, cone_vol, float(args.BottomOffset))
    aspirated_water_vol = 0
    bottom_offset = float(args.BottomOffset)

    # First pipette water into each well.
    sample_dest_dict = defaultdict(list)
    for sample_key in sample_parameters:
        sample_concentration = float(sample_parameters[sample_key][4])
        sample_volume_needed = sample_mass/sample_concentration
        sample_dest_slot = sample_parameters[sample_key][5]
        sample_dest_well = sample_parameters[sample_key][6]
        sample_destination_labware = labware_dict[sample_dest_slot]

        # Apply some sanity to the sample volumes
        sample_volume = round(sample_volume_needed, 1)
        water_volume = (0.5*reaction_vol)-sample_volume

        # Define the pipette for dispensing the water.
        water_pipette, water_loop, water_volume = pipette_selection(left_pipette, right_pipette, water_volume)

        sample_dest_dict[sample_key].append((sample_dest_slot, sample_dest_well, sample_volume))

        # Add water to all the destination wells for this sample.
        dispensing_loop(args, water_loop, water_pipette, reagent_labware[args.WaterWell].bottom(tip_height),
                        sample_destination_labware[sample_dest_well], water_volume, NewTip=False, MixReaction=False)

        aspirated_water_vol += water_volume
        tip_height = res_tip_height(float(args.WaterResVol) - aspirated_water_vol, water_reservoir_dia, cone_vol,
                                    float(args.BottomOffset))

    # Drop any tips the pipettes might have.
    if left_pipette.has_tip:
        left_pipette.drop_tip()
    if right_pipette.has_tip:
        right_pipette.drop_tip()

    # Change pipettes.
    primer_labware = labware_dict[args.IndexPrimersSlot]
    for sample_key in sample_dest_dict:
        d500, d700 = sample_parameters[sample_key][2].split("+")
        sample_slot = sample_parameters[sample_key][0]
        sample_source_well = sample_parameters[sample_key][1]
        sample_source_labware = labware_dict[sample_slot]
        sample_dest_slot = sample_dest_dict[sample_key][0][0]
        sample_dest_well = sample_dest_dict[sample_key][0][1]
        sample_destination_labware = labware_dict[sample_dest_slot]
        sample_volume = sample_dest_dict[sample_key][0][2]

        sample_pipette, sample_loop, sample_volume = \
            pipette_selection(left_pipette, right_pipette, sample_volume)

        # Add template to the destination well for this sample.
        dispensing_loop(args, sample_loop, sample_pipette,
                        sample_source_labware[sample_source_well].bottom(bottom_offset),
                        sample_destination_labware[sample_dest_well], sample_volume, NewTip=True, MixReaction=False,
                        touch=True)

        # Determine primer volumes and dispense them.
        primer_volume = (float(args.PCR_Volume)/50) * 1.25
        primer_pipette, primer_loop, primer_volume = \
            pipette_selection(left_pipette, right_pipette, primer_volume)

        # D500 primers
        dispensing_loop(args, primer_loop, primer_pipette, primer_labware[primer_dict[d500]].bottom(bottom_offset),
                        sample_destination_labware[sample_dest_well], primer_volume, NewTip=True, MixReaction=False,
                        touch=True)
        # D700 primers
        dispensing_loop(args, primer_loop, primer_pipette, primer_labware[primer_dict[d700]].bottom(bottom_offset),
                        sample_destination_labware[sample_dest_well], primer_volume, NewTip=True, MixReaction=False,
                        touch=True,)

    return sample_dest_dict


def run(ctx: protocol_api.ProtocolContext):
    ctx.comment("Begin {}".format(metadata['protocolName']))

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

    # Dispense Samples and primers
    sample_dest_dict = dispense_samples(args, sample_parameters, labware_dict, left_pipette, right_pipette)

    # Add PCR mix to each destination well.
    add_pcr_mix(args, labware_dict, sample_dest_dict, left_pipette, right_pipette)

    if not ctx.is_simulating():
        os.remove(tsv_file_path)

    ctx.comment("Program End")


if __name__ == "__main__":
    protocol_file = open('Illumina_Dual_Indexing.py')
    labware_path = "{}{}custom_labware".format(os.getcwd(), os.sep)
    run_log, __bundle__ = simulate(protocol_file, custom_labware_paths=[labware_path])
    run_date = datetime.datetime.today().strftime("%a %b %d %H:%M %Y")
    i = 1
    t = format_runlog(run_log).split("\n")

    outstring = "Opentrons OT-2 Steps for {}.\nDate:  {}\nProgram File: Illumina_Dual_Indexing.py\n\nStep\tCommand\n" \
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
