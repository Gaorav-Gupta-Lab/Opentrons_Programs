import csv
import os
import re
from collections import defaultdict
from types import SimpleNamespace
from opentrons.simulate import simulate, format_runlog
from opentrons import protocol_api
import opentrons

# metadata
metadata = {
    'protocolName': 'Generic PCR v0.1.0',
    'author': 'Dennis Simpson',
    'description': 'Sets up a PCR from concentrated template',
    'apiLevel': '2.8'
}


def parse_sample_file(input_file):
    """
    Parse the input file.
    :return:
    """

    line_num = 0
    options_dictionary = defaultdict(str)
    sample_dictionary = defaultdict(list)
    index_file = list(csv.reader(open(input_file), delimiter='\t'))
    for line in index_file:
        line_num += 1
        col_count = len(line)
        tmp_line = []
        sample_key = ""
        if col_count > 0 and "#" not in line[0] and len(line[0].split("#")[0]) > 0:
            # Skip any lines that are blank or comments.

            for i in range(6):
                try:
                    line[i] = line[i].split("#")[0]  # Strip out end of line comments and white space.
                except IndexError:
                    raise SystemExit("There is a syntax error in file {0} on line {1}, column {2} "
                                     .format(input_file, str(line_num), str(i)))

                line[i] = re.sub(",", '', line[i])  # Strip out any commas.

                if i == 0 and "--" in line[0]:
                    key = line[0].strip('--')
                    options_dictionary[key] = line[1]
                elif "--" not in line[0] and int(line[0]) < 12:
                    sample_key = line[0], line[1]
                    tmp_line.append(line[i])
            if sample_key:
                sample_dictionary[sample_key] = tmp_line

    return sample_dictionary, SimpleNamespace(**options_dictionary)


def run(protocol: protocol_api.ProtocolContext):
    # load sample data
    # opentrons.protocol_api.labware.LocationLabware()
    # Turn on rail lights and pause program so user can load robot deck.
    protocol.set_rail_lights(True)
    protocol.pause("Load Labware onto robot deck and hit resume when ready to continue")

    sample_parameters, args = parse_sample_file("{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep))
    protocol.set_rail_lights(False)
    sample_mass = float(args.DNA_in_Reaction)
    reaction_vol = float(args.PCR_Volume)

    # Extract Slot information
    slot_list = ["Slot1", "Slot2", "Slot3", "Slot4", "Slot5", "Slot6", "Slot7", "Slot8", "Slot9", "Slot10", "Slot11"]
    slot_dict = {}

    for i in range(len(slot_list)):
        slot_dict[str(i+1)] = getattr(args, "{}".format(slot_list[i]))

    # Labware
    reagent_source = protocol.load_labware(slot_dict[args.ReagentSlot], args.ReagentSlot)
    left_tiprack_list = ""
    right_tiprack_list = ""

    if args.LeftPipetteTipRackSlot:
        left_tiprack_list = args.LeftPipetteTipRackSlot.split(",")
    if args.RightPipetteTipRackSlot:
        right_tiprack_list = args.RightPipetteTipRackSlot.split(",")

    left_tipracks = []
    if left_tiprack_list:
        for slot in left_tiprack_list:
            if slot not in protocol.loaded_labwares:
                left_tipracks.append(protocol.load_labware(slot_dict[args.LeftPipetteTipRackSlot], str(slot)))

    right_tipracks = []
    if right_tiprack_list:
        for slot in right_tiprack_list:
            if slot not in protocol.loaded_labwares:
                right_tipracks.append((protocol.load_labware(slot_dict[args.RightPipetteTipRackSlot], str(slot))))

    # Pipettes
    left_pipette = protocol.load_instrument(args.LeftPipette, 'left', tip_racks=left_tipracks)
    right_pipette = protocol.load_instrument(args.RightPipette, 'right', tip_racks=right_tipracks)

    # Load sample data
    sample_dict = defaultdict(list)
    destination_dict = {}
    for sample_key in sample_parameters:
        concentration = float(sample_parameters[sample_key][3])
        volume = sample_mass/concentration
        sample_slot = sample_parameters[sample_key][0]
        sample_position = sample_parameters[sample_key][1]
        sample_name = sample_parameters[sample_key][2]
        sample_dest_slot = sample_parameters[sample_key][4]
        sample_dest_well = sample_parameters[sample_key][5]

        if volume > 20:
            dna_volume = round(volume, 1)
        else:
            dna_volume = round(volume, 2)

        water_volume = round((0.5*reaction_vol)-dna_volume, 2)

        if sample_dest_slot not in destination_dict:
            destination_dict[sample_dest_slot] = protocol.load_labware(slot_dict[sample_dest_slot], sample_dest_slot)

        sample_dict[sample_slot].append(
            (sample_position, sample_name, dna_volume, water_volume, sample_dest_slot, sample_dest_well))

    # Dilute DNA in PCR wells
    temp_well_list = []
    for sample_slot in sample_dict:
        sample_source = protocol.load_labware(slot_dict[sample_slot], sample_slot)
        for sample in sample_dict[sample_slot]:
            sample_pipette = left_pipette
            water_pipette = left_pipette
            sample_destination_slot = sample[4]
            sample_well = sample[0]
            sample_name = sample[1]
            water_volume = sample[3]

            # Let's deal with some potential errors.
            if (float(args.PCR_Volume) * 0.5) < water_volume or water_volume < 0:
                raise SystemExit("ERROR: Water Volume is {} uL on sample {}. PCR volume of {} uL using a 2x PCR reagent"
                                 " mix will not work with this sample."
                                 .format(water_volume, sample_name, args.PCR_Volume))

            sample_volume = sample[2]
            sample_dest_well = sample[5]
            sample_destination_labware = destination_dict[sample_destination_slot]
            temp_well_list.append(sample_well)

            # Define the pipettes for dispensing our samples and the water.
            sample_loop = 1
            if sample_volume >= 20:
                sample_pipette = right_pipette
            elif 10 <= sample_volume < 20:
                sample_pipette = left_pipette
                sample_volume = sample_volume * 0.5
                sample_loop = 2
            print(water_volume)
            water_loop = 1
            if water_volume >= 20:
                water_pipette = right_pipette
            elif 10 <= water_volume < 20:
                water_pipette = left_pipette
                water_volume = water_volume * 0.5
                water_loop = 2

            # Dispense the water.
            water_pipette.pick_up_tip()
            while water_loop > 0:
                water_pipette.aspirate(water_volume, reagent_source[args.Water])
                water_pipette.dispense(water_volume, sample_destination_labware[sample_dest_well])
                water_loop -= 1

            if not water_pipette == sample_pipette:
                water_pipette.drop_tip()
                sample_pipette.pick_up_tip()

            # Dispense the samples into the wells that have water.
            while sample_loop > 0:
                sample_pipette.aspirate(sample_volume, sample_source[sample_well])
                sample_pipette.dispense(sample_volume, sample_destination_labware[sample_dest_well])
                sample_loop -= 1

            sample_pipette.drop_tip()

    # Define PCR reagent pipette.
    pcr_reagent_vol = float(args.PCR_Volume) * 0.5
    pcr_pipette = left_pipette
    pcr_loop = 1
    neg_control_loop = 1
    if pcr_reagent_vol >= 20:
        pcr_pipette = right_pipette
    elif 10 < pcr_reagent_vol < 20:
        pcr_pipette = left_pipette
        pcr_reagent_vol = pcr_reagent_vol * 0.5
        pcr_loop = 2
        neg_control_loop = 2

    # Add PCR Reagents
    mix_vol = float(args.PCR_Volume) * 0.75
    for destination_slot in destination_dict:
        pcr_plate = destination_dict[destination_slot]
        for i in range(len(temp_well_list)):
            pcr_pipette.pick_up_tip()

            while pcr_loop > 0:
                pcr_pipette.aspirate(pcr_reagent_vol, reagent_source[args.PCR_Mix])
                pcr_pipette.dispense(pcr_reagent_vol, pcr_plate[temp_well_list[i]])
                pcr_loop -= 1

            pcr_pipette.mix(3, mix_vol, pcr_plate[temp_well_list[i]])
            pcr_pipette.drop_tip()

    # The negative control pipette will be the same as the pcr_pipette
    neg_control_pipette = pcr_pipette
    if args.WaterControl:
        slot = args.WaterControl.split(',')[0]
        well = args.WaterControl.split(',')[1]
        pcr_plate = destination_dict[slot]
        mix_vol = float(args.PCR_Volume) * 0.75
        neg_control_pipette.pick_up_tip()

        while neg_control_loop > 0:
            neg_control_pipette.aspirate(pcr_reagent_vol, reagent_source[args.Water])
            neg_control_pipette.dispense(pcr_reagent_vol, pcr_plate[well])
            neg_control_pipette.aspirate(pcr_reagent_vol, reagent_source[args.PCR_Mix])
            neg_control_pipette.dispense(pcr_reagent_vol, pcr_plate[well])
            neg_control_loop -= 1

        neg_control_pipette.mix(3, mix_vol, pcr_plate[well])
        neg_control_pipette.drop_tip()


if __name__ == "__main__":
    protocol_file = open('Generic_PCR_v0.1.0.py')
    run_log, __bundle__ = simulate(protocol_file)
    print(format_runlog(run_log))
    protocol_file.close()
