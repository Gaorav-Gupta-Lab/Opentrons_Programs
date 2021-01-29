import csv
import math
import os
from collections import defaultdict
from types import SimpleNamespace
from opentrons.simulate import simulate, format_runlog
# import opentrons

# metadata
metadata = {
    'protocolName': 'Illumina Dual Indexing v0.4.0',
    'author': 'Dennis Simpson',
    'description': 'Add Illumina dual indexing to library',
    'apiLevel': '2.8'
}


def parse_sample_file(input_file):
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
        line_num += 1
        col_count = len(line)
        tmp_line = []
        sample_key = ""
        if col_count > 0 and "#" not in line[0] and len(line[0].split("#")[0]) > 0:
            # Skip any lines that are blank or comments.
            for i in range(7):
                try:
                    line[i] = line[i].split("#")[0]  # Strip out end of line comments and white space.
                except IndexError:
                    continue

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


def add_pcr_mix(args, labware_dict, sample_dest_dict, left_pipette, right_pipette):
    pcr_reagent_vol = float(args.PCR_Volume) * 0.5
    pcr_reagent_labware = labware_dict[args.ReagentSlot]
    pcr_reagent_well = args.PCR_MixWell
    reagent_labware = labware_dict[args.ReagentSlot]
    pcr_reservoir_dia = reagent_labware[args.PCR_MixWell].diameter
    cone_vol = labware_cone_volume(args, reagent_labware)
    tip_height = res_tip_height(float(args.PCR_MixResVolume), pcr_reservoir_dia, cone_vol)
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

        tip_height = res_tip_height(float(args.PCR_MixResVolume)-aspirated_vol, pcr_reservoir_dia, cone_vol)



def pipette_selection(left_pipette, right_pipette, volume):
    """
    Function to select pipette based on expected volumes.  Will also adjust volume is pipette needs to pick up >1x
    @param left_pipette:
    @param right_pipette:
    @param volume:
    @return:
    """
    loop = 1
    pipette = ""
    if volume > 20 and "P300 Single-Channel GEN2" in str(right_pipette):
        pipette = right_pipette
    elif volume <= 20 and "P20 Single-Channel GEN2" in str(left_pipette):
        pipette = left_pipette
    elif volume < 10 and "P10 Single-Channel GEN1" in str(left_pipette):
        pipette = left_pipette
    elif volume < 10 and "P10 Single-Channel GEN1" in str(right_pipette):
        pipette = right_pipette
    elif 10 <= volume <= 20 and "P10 Single-Channel GEN1" in str(left_pipette):
        pipette = left_pipette
        volume = volume * 0.5
        loop = 2
    elif 10 <= volume <= 20 and "P10 Single-Channel GEN1" in str(right_pipette):
        pipette = right_pipette
        volume = volume * 0.5
        loop = 2

    return pipette, loop, volume


def res_tip_height(res_vol, well_dia, cone_vol):
    """
    Calculate the the height of the liquid in a reservoir and return the value to set the pipette tip height.
    This works for both conical shapes and cylinders.
    @param res_vol:
    @param well_dia:
    @param cone_vol:
    @return:
    """
    if res_vol > cone_vol:
        cone_height = (3*cone_vol / (math.pi * ((well_dia / 2) ** 2)))
        height = ((res_vol-cone_vol)/(math.pi*((well_dia/2)**2)))-9+cone_height
    else:
        height = (3*res_vol / (math.pi * ((well_dia / 2) ** 2))) - 4

    if height <= 8:
        height = 1

    return int(height)


def labware_cone_volume(args, labware_name):
    """
    Based on the labware and reservoir return the volume at which the cylinder shape transitions to the conical shape.
    @param args:
    @param labware_name:
    @return:
    """
    cone_vol = 200
    labware = getattr(args, "Slot{}".format(str(labware_name)[-1:]))

    if "e5ml_" in labware:
        cone_vol = 1200

    elif"1.5ml_24" in labware:
        cone_vol = 500

    return cone_vol


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
    tip_height = res_tip_height(float(args.WaterResVol), water_reservoir_dia, cone_vol)
    aspirated_water_vol = 0

    # First pipette water into each well.
    sample_dest_dict = defaultdict(list)
    for sample_key in sample_parameters:
        sample_concentration = float(sample_parameters[sample_key][4])
        sample_volume_needed = sample_mass/sample_concentration
        sample_dest_slot = sample_parameters[sample_key][5]
        sample_dest_well = sample_parameters[sample_key][6]
        sample_destination_labware = labware_dict[sample_dest_slot]

        # Apply some sanity to the sample volumes
        if sample_volume_needed > 20:
            sample_volume = round(sample_volume_needed, 1)
        else:
            sample_volume = round(sample_volume_needed, 2)

        water_volume = (0.5*reaction_vol)-sample_volume

        # Define the pipette for dispensing the water.
        water_pipette, water_loop, water_volume = pipette_selection(left_pipette, right_pipette, water_volume)

        sample_dest_dict[sample_key].append((sample_dest_slot, sample_dest_well, sample_volume))

        # Add water to all the destination wells for this sample.
        if not water_pipette.has_tip:
            water_pipette.pick_up_tip()
        dispensing_loop(args, water_loop, water_pipette, reagent_labware[args.WaterWell].bottom(tip_height),
                        sample_destination_labware[sample_dest_well], water_volume, NewTip=False, MixReaction=False)
        aspirated_water_vol += water_volume
        tip_height = res_tip_height(float(args.WaterResVol) - aspirated_water_vol, water_reservoir_dia, cone_vol)

    # Drop any tips the pipettes might have.
    left_pipette.drop_tip()
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
        dispensing_loop(args, sample_loop, sample_pipette, sample_source_labware[sample_source_well],
                        sample_destination_labware[sample_dest_well], sample_volume, NewTip=True, MixReaction=False)

        # Determine primer volumes and dispense them.
        primer_volume = (float(args.PCR_Volume)/50) * 1.25
        primer_pipette, primer_loop, primer_volume = \
            pipette_selection(left_pipette, right_pipette, primer_volume)

        # D500 primers
        dispensing_loop(args, primer_loop, primer_pipette, primer_labware[primer_dict[d500]],
                        sample_destination_labware[sample_dest_well], primer_volume, NewTip=True, MixReaction=False)
        # D700 primers
        dispensing_loop(args, primer_loop, primer_pipette, primer_labware[primer_dict[d700]],
                        sample_destination_labware[sample_dest_well], primer_volume, NewTip=True, MixReaction=False)

    return sample_dest_dict


def dispensing_loop(args, loop_count, pipette, source_location, destination_location, volume, NewTip, MixReaction):
    if NewTip:
        if pipette.has_tip:
            pipette.drop_tip()
        pipette.pick_up_tip()
    while loop_count > 0:
        pipette.aspirate(volume, source_location)
        pipette.dispense(volume, destination_location)
        loop_count -= 1
        if not MixReaction:
            pipette.blow_out()

    if MixReaction:
        pipette.mix(3, float(args.PCR_Volume)*0.7)
    if NewTip:
        pipette.drop_tip()

    return pipette


def run(ctx):

    # TSV file location on OT-2
    tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)
    if not os.path.isfile(tsv_file_path):
        # Temp TSV file location on Win10 Computers for simulation
        tsv_file_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())

    sample_parameters, args = parse_sample_file(tsv_file_path)

    # Turn on rail lights and pause program so user can load robot deck.
    ctx.set_rail_lights(True)
    ctx.pause("Load Labware onto robot deck and click resume when ready to continue")
    ctx.home()
    ctx.set_rail_lights(False)

    # Extract Slot information
    slot_list = ["Slot1", "Slot2", "Slot3", "Slot4", "Slot5", "Slot6", "Slot7", "Slot8", "Slot9", "Slot10", "Slot11"]

    labware_dict = {}
    for i in range(len(slot_list)):
        labware = getattr(args, "{}".format(slot_list[i]))
        if labware:
            labware_dict[str(i+1)] = ctx.load_labware(labware, str(i+1))

    # Pipette Tip Boxes
    left_tipracks = ""
    right_tipracks = ""
    if args.LeftPipetteTipRackSlot:
        left_tipracks = load_tipracks(ctx, args.LeftPipetteTipRackSlot.split(","), labware_dict)

    if args.RightPipetteTipRackSlot:
        right_tipracks = load_tipracks(ctx, args.RightPipetteTipRackSlot.split(","), labware_dict)

    # Pipettes
    left_pipette = ctx.load_instrument(args.LeftPipette, 'left', tip_racks=left_tipracks)
    right_pipette = ctx.load_instrument(args.RightPipette, 'right', tip_racks=right_tipracks)

    # Set the location of the first tip in box.
    left_pipette.starting_tip = left_tipracks[0].wells_by_name()[args.LeftPipetteFirstTip]
    right_pipette.starting_tip = right_tipracks[0].wells_by_name()[args.RightPipetteFirstTip]

    # Dispense Samples and primers
    sample_dest_dict = dispense_samples(args, sample_parameters, labware_dict, left_pipette, right_pipette)

    # Add PCR mix to each destination well.
    add_pcr_mix(args, labware_dict, sample_dest_dict, left_pipette, right_pipette)

    if not ctx.is_simulating():
        os.remove(tsv_file_path)


if __name__ == "__main__":
    protocol_file = open('Illumina_Dual_Indexing.py')
    labware_path = "{}{}custom_labware".format(os.getcwd(), os.sep)
    run_log, __bundle__ = simulate(protocol_file, custom_labware_paths=[labware_path])
    print(format_runlog(run_log))
    protocol_file.close()
