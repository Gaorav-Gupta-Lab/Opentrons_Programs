import datetime
import os
import platform
import sys
import csv
from types import SimpleNamespace
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

import Utilities

# metadata
metadata = {
    'protocolName': 'Illumina Dual Indexing v1.1.0',
    'author': 'Dennis Simpson',
    'description': 'Add Illumina dual indexing to NGS library',
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


def add_pcr_mix(args, labware_dict, sample_dest_dict, left_pipette, right_pipette):
    pcr_reagent_vol = float(args.PCR_Volume) * 0.5
    pcr_reagent_labware = labware_dict[args.ReagentSlot]
    pcr_reagent_well = args.PCR_ReagentWell.upper()
    reagent_labware = labware_dict[args.ReagentSlot]
    pcr_reservoir_dia = reagent_labware[args.PCR_ReagentWell.upper()].diameter
    cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
    tip_height = \
        Utilities.res_tip_height(float(args.TotalReagentVolume), pcr_reservoir_dia, cone_vol, float(args.BottomOffset))
    pcr_pipette, pcr_loop, pcr_reagent_vol = Utilities.pipette_selection(left_pipette, right_pipette, pcr_reagent_vol)
    aspirated_vol = 0

    for sample_key in sample_dest_dict:
        sample_dest_slot = sample_dest_dict[sample_key][0][0]
        sample_dest_well = sample_dest_dict[sample_key][0][1]
        sample_destination_labware = labware_dict[sample_dest_slot]

        pcr_pipette = \
            Utilities.dispensing_loop(args, pcr_loop, pcr_pipette,
                                      pcr_reagent_labware[pcr_reagent_well].bottom(tip_height),
                                      sample_destination_labware[sample_dest_well], pcr_reagent_vol,
                                      NewTip=True, MixReaction=True)

        aspirated_vol += pcr_reagent_vol

        tip_height = Utilities.res_tip_height(float(args.TotalReagentVolume)-aspirated_vol, pcr_reservoir_dia, cone_vol,
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
    water_reservoir_dia = reagent_labware[args.WaterResWell.upper()].diameter
    cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
    tip_height = \
        Utilities.res_tip_height(float(args.WaterResVol), water_reservoir_dia, cone_vol, float(args.BottomOffset))
    aspirated_water_vol = 0
    bottom_offset = float(args.BottomOffset)
    sample_dest_slot = args.PCR_PlateSlot

    # First pipette water into each well.
    sample_dest_dict = defaultdict(list)
    for sample_key in sample_parameters:
        sample_concentration = float(sample_parameters[sample_key][4])
        sample_volume_needed = sample_mass/sample_concentration
        sample_dest_well = sample_parameters[sample_key][5]
        sample_destination_labware = labware_dict[sample_dest_slot]

        # Apply some sanity to the sample volumes
        sample_volume = round(sample_volume_needed, 1)
        water_volume = (0.5*reaction_vol)-sample_volume

        # Define the pipette for dispensing the water.
        water_pipette, water_loop, water_volume = Utilities.pipette_selection(left_pipette, right_pipette, water_volume)

        sample_dest_dict[sample_key].append((sample_dest_slot, sample_dest_well, sample_volume))

        # Add water to all the destination wells for this sample.
        Utilities.dispensing_loop(args, water_loop, water_pipette, reagent_labware[args.WaterResWell.upper()].bottom(tip_height),
                                  sample_destination_labware[sample_dest_well], water_volume, NewTip=False,
                                  MixReaction=False)

        aspirated_water_vol += water_volume
        tip_height = \
            Utilities.res_tip_height(float(args.WaterResVol) - aspirated_water_vol, water_reservoir_dia, cone_vol,
                                     float(args.BottomOffset))

    # Drop any tips the pipettes might have.
    if left_pipette.has_tip:
        left_pipette.drop_tip()
    if right_pipette.has_tip:
        right_pipette.drop_tip()

    # Change pipettes.
    primer_labware = labware_dict[args.IndexPrimerSlot]
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
            Utilities.pipette_selection(left_pipette, right_pipette, sample_volume)

        # Add template to the destination well for this sample.
        Utilities.dispensing_loop(args, sample_loop, sample_pipette,
                                  sample_source_labware[sample_source_well].bottom(bottom_offset),
                                  sample_destination_labware[sample_dest_well], sample_volume, NewTip=True,
                                  MixReaction=False, touch=True)

        # Determine primer volumes and dispense them.
        # 6.25 uM = 2 uL per 50
        # primer_volume = (float(args.PCR_Volume)/50) * 1.25
        primer_volume = (float(args.PCR_Volume)/50) * 2.0
        primer_pipette, primer_loop, primer_volume = \
            Utilities.pipette_selection(left_pipette, right_pipette, primer_volume)

        # D500 primers
        Utilities.dispensing_loop(args, primer_loop, primer_pipette,
                                  primer_labware[primer_dict[d500]].bottom(bottom_offset),
                                  sample_destination_labware[sample_dest_well], primer_volume, NewTip=True,
                                  MixReaction=False, touch=True)
        # D700 primers
        Utilities.dispensing_loop(args, primer_loop, primer_pipette,
                                  primer_labware[primer_dict[d700]].bottom(bottom_offset),
                                  sample_destination_labware[sample_dest_well], primer_volume, NewTip=True,
                                  MixReaction=False, touch=True,)

    return sample_dest_dict


def run(ctx: protocol_api.ProtocolContext):

    # TSV file location on OT-2
    tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)
    if not os.path.isfile(tsv_file_path):
        # Temp TSV file location on Win10 Computers for simulation
        tsv_file_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())

    sample_parameters, args = parse_sample_template(tsv_file_path)

    ctx.comment("Begin {}".format(metadata['protocolName']))

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
