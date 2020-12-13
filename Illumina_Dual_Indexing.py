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
    'protocolName': 'Illumina Dual Indexing v0.1.0',
    'author': 'Dennis Simpson',
    'description': 'Add Illumina dual indexing to library',
    'apiLevel': '2.8'
}


def parse_sample_file(input_file, protocol):
    """
    Parse the input file.
    :return:
    """
    protocol.set_rail_lights(False)
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


def add_pcr_mix(reagent_source, pipette, args, sample_parameters, sample_reagent_labware_dict):
    pcr_reagent_vol = float(args.PCR_Volume) * 0.5
    for sample in sample_parameters:
        sample_dest_slot = sample_parameters[sample][5]
        sample_dest_well = sample_parameters[sample][6]
        sample_destination_labware = sample_reagent_labware_dict[sample_dest_slot]
        reagent_destination = sample_destination_labware[sample_dest_well]
        mix_vol = 20
        pipette.pick_up_tip()
        pipette.aspirate(pcr_reagent_vol, reagent_source[args.PCR_Mix])
        pipette.dispense(pcr_reagent_vol, reagent_destination)
        pipette.mix(3, mix_vol, reagent_destination)
        pipette.drop_tip()


def run(protocol: protocol_api.ProtocolContext):
    # load sample data
    # opentrons.protocol_api.labware.LocationLabware()
    # Turn on rail lights and pause program so user can load robot deck.
    protocol.set_rail_lights(True)
    protocol.pause("Load Labware onto robot deck and hit resume when ready to continue")

    sample_parameters, args = \
        parse_sample_file("{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep), protocol)

    sample_mass = float(args.DNA_in_Reaction)
    reaction_vol = float(args.PCR_Volume)

    # Extract Slot information
    slot_list = ["Slot1", "Slot2", "Slot3", "Slot4", "Slot5", "Slot6", "Slot7", "Slot8", "Slot9", "Slot10", "Slot11"]
    slot_dict = {}

    for i in range(len(slot_list)):
        slot_dict[str(i+1)] = getattr(args, "{}".format(slot_list[i]))

    # Labware
    reagent_source = protocol.load_labware(slot_dict[args.ReagentSlot], args.ReagentSlot)
    index_primer_source = protocol.load_labware(slot_dict[args.IndexPrimersSlot], args.IndexPrimersSlot)
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

    sample_reagent_labware_dict = build_labware_dict(protocol, sample_parameters, slot_dict)

    # Load sample data
    index_primer_vol = round(float(args.PCR_Volume)*0.025, 1)
    for sample_key in sample_parameters:
        sample_source_slot = sample_parameters[sample_key][0]
        sample_source_well = sample_parameters[sample_key][1]
        sample_index = sample_parameters[sample_key][2]
        # sample_name = sample_parameters[sample_key][3]
        concentration = float(sample_parameters[sample_key][4])
        sample_dest_slot = sample_parameters[sample_key][5]
        sample_dest_well = sample_parameters[sample_key][6]
        sample_volume = sample_mass / concentration

        if sample_volume > 20:
            dna_volume = round(sample_volume, 1)
        else:
            dna_volume = round(sample_volume, 2)

        water_volume = round((0.5*reaction_vol)-dna_volume, 2)

        sample_source_labware = sample_reagent_labware_dict[sample_source_slot]
        sample_destination_labware = sample_reagent_labware_dict[sample_dest_slot]

        # Define Pipettes
        sample_pipette = left_pipette
        if sample_volume > 20:
            sample_pipette = right_pipette
        water_pipette = left_pipette
        if water_volume > 20:
            water_pipette = right_pipette

        # Dispense the water.
        water_pipette.pick_up_tip()
        water_pipette.aspirate(water_volume, reagent_source[args.Water])
        water_pipette.dispense(water_volume, sample_destination_labware[sample_dest_well])

        if water_pipette == sample_pipette:
            sample_pipette = water_pipette
        else:
            water_pipette.drop_tip()
            sample_pipette.pick_up_tip()

        # Dispense the samples into the wells that have water.
        sample_pipette.aspirate(sample_volume, sample_source_labware[sample_source_well])
        sample_pipette.dispense(sample_volume, sample_destination_labware[sample_dest_well])
        sample_pipette.drop_tip()

        # Dispense Indexing Primers.  D500 first then D700
        primer_pipette = left_pipette
        D500_sample_index = sample_index.split("+")[0]
        D700_sample_index = sample_index.split("+")[1]
        index_primer_wells = \
            [getattr(args, "{}".format(D500_sample_index)), getattr(args, "{}".format(D700_sample_index))]

        for index_primer_well in index_primer_wells:
            primer_pipette.pick_up_tip()
            primer_pipette.aspirate(index_primer_vol, index_primer_source[index_primer_well])
            primer_pipette.dispense(index_primer_vol, sample_destination_labware[sample_dest_well])
            primer_pipette.drop_tip()

    # Add PCR reagents to each well
    pcr_reagent_vol = float(args.PCR_Volume) * 0.5
    pcr_pipette = left_pipette
    if pcr_reagent_vol > 20:
        pcr_pipette = right_pipette

    add_pcr_mix(reagent_source, pcr_pipette, args, sample_parameters, sample_reagent_labware_dict)


if __name__ == "__main__":
    protocol_file = open('Illumina_Dual_Indexing.py')
    run_log, __bundle__ = simulate(protocol_file)
    print(format_runlog(run_log))
    protocol_file.close()
