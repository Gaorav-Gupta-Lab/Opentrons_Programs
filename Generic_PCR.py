import csv
import os
import math
from collections import defaultdict
from types import SimpleNamespace
from opentrons.simulate import simulate, format_runlog
# from opentrons import protocol_api
# import opentrons

# metadata
metadata = {
    'protocolName': 'Generic PCR v0.3.4',
    'author': 'Dennis Simpson',
    'description': 'Sets up a PCR from concentrated template',
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
            for i in range(6):
                try:
                    line[i] = line[i].split("#")[0]  # Strip out end of line comments and white space.
                except IndexError:
                    raise SystemExit("There is a syntax error in file {0} on line {1}, column {2} "
                                     .format(input_file, str(line_num), str(i)))

                if i == 0 and "--" in line[0]:
                    key = line[0].strip('--')
                    options_dictionary[key] = line[1]
                elif "--" not in line[0] and int(line[0]) < 12:
                    sample_key = line[0], line[1]
                    tmp_line.append(line[i])
            if sample_key:
                sample_dictionary[sample_key] = tmp_line

    return sample_dictionary, SimpleNamespace(**options_dictionary)


def sample_pipette_selection(left_pipette, right_pipette, volume):
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


def res_tip_height(res_vol, well_dia, cone_vol):
    """
    Calculate the the height of the liquid in a conical tube and return the value to set the pipette tip height.
    @param res_vol:
    @param well_dia:
    @param cone_vol:
    @return:
    """
    if res_vol > cone_vol:
        cone_height = (3*cone_vol / (math.pi * ((well_dia / 2) ** 2)))
        height = ((res_vol-cone_vol)/(math.pi*((well_dia/2)**2)))-9+cone_height
    else:
        height = (3*res_vol / (math.pi * ((well_dia / 2) ** 2))) - 3

    if height <= 10:
        height = 1

    return int(height)


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

    # Dispense Samples
    sample_dest_dict, remaining_water_vol = \
        dispense_samples(args, sample_parameters, labware_dict, left_pipette, right_pipette)

    # Dispense PCR Reagents.
    dispense_pcr_reagents(args, labware_dict, left_pipette, right_pipette, remaining_water_vol, sample_dest_dict)

    if not ctx.is_simulating():
        os.remove(tsv_file_path)


def dispense_pcr_reagents(args, labware_dict, left_pipette, right_pipette, remaining_water_vol, sample_dest_dict):
    """
    Dispense PCR reagents into each well and setup negative control well.
    @param args:
    @param labware_dict:
    @param left_pipette:
    @param right_pipette:
    @param remaining_water_vol:
    @param sample_dest_dict:
    """

    # Define PCR reagent pipette and reagent volume.
    reagent_pipette, pcr_loop, pcr_reagent_vol = \
        sample_pipette_selection(left_pipette, right_pipette, float(args.PCR_Volume) * 0.5)

    reagent_labware = labware_dict[args.ReagentSlot]

    pcr_reagent_well_dia = reagent_labware[args.PCR_MixWell].diameter
    trigger_vol, cone_vol = height_adjust_trigger(args, reagent_labware)
    tip_height = res_tip_height(float(args.PCR_MixResVolume), pcr_reagent_well_dia, cone_vol)
    incremental_aspirated_vol = 0
    reagent_aspirated = 0

    # Add PCR Reagents
    mix_vol = float(args.PCR_Volume) * 0.80

    for destination_slot in sample_dest_dict:
        destination_labware = labware_dict[destination_slot]
        destination_well_list = sample_dest_dict[destination_slot]

        for destination_well in destination_well_list:
            reagent_pipette.pick_up_tip()
            temp_pcr_loop = pcr_loop

            while temp_pcr_loop > 0:
                reagent_pipette.aspirate(pcr_reagent_vol, reagent_labware[args.PCR_MixWell].bottom(tip_height))
                reagent_pipette.dispense(pcr_reagent_vol, destination_labware[destination_well])
                incremental_aspirated_vol += pcr_reagent_vol
                reagent_aspirated += pcr_reagent_vol
                tip_height = \
                    res_tip_height(float(args.PCR_MixResVolume) - reagent_aspirated, pcr_reagent_well_dia, cone_vol)

                temp_pcr_loop -= 1

            reagent_pipette.mix(3, mix_vol, destination_labware[destination_well])
            reagent_pipette.drop_tip()

    # Setup the water control sample(s)
    if args.WaterControl:
        water_reservoir_dia = reagent_labware[args.WaterWell].diameter
        water_tip_height = \
            res_tip_height(float(args.PCR_MixResVolume) - remaining_water_vol, water_reservoir_dia, cone_vol)
        slot = args.WaterControl.split(',')[0]
        well = args.WaterControl.split(',')[1]
        neg_control_labware = labware_dict[slot]
        temp_pcr_loop = pcr_loop

        reagent_pipette.pick_up_tip()

        while temp_pcr_loop > 0:
            # Dispense Water into negative control well
            reagent_pipette.aspirate(pcr_reagent_vol, reagent_labware[args.WaterWell].bottom(water_tip_height))
            reagent_pipette.dispense(pcr_reagent_vol, neg_control_labware[well])

            # Dispense PCR reagents into negative control well
            reagent_pipette.aspirate(pcr_reagent_vol, reagent_labware[args.PCR_MixWell])
            reagent_pipette.dispense(pcr_reagent_vol, neg_control_labware[well])
            temp_pcr_loop -= 1

        reagent_pipette.mix(3, mix_vol, neg_control_labware[well])
        reagent_pipette.drop_tip()


def height_adjust_trigger(args, labware_name):
    trigger_vol = 200
    cone_vol = 200
    labware = getattr(args, "Slot{}".format(str(labware_name)[-1:]))

    if "5ml_" in labware:
        trigger_vol = 500
        cone_vol = 1200

    elif"1.5ml" in labware:
        trigger_vol = 225
        cone_vol = 500

    return trigger_vol, cone_vol


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
    trigger_vol, cone_vol = height_adjust_trigger(args, reagent_labware)
    tip_height = res_tip_height(float(args.WaterResVol), water_reservoir_dia, cone_vol)
    incremental_aspirated_vol = 0
    aspirated_vol = 0

    # This is used to dispense the PCR reagents.
    sample_dest_dict = defaultdict(list)

    for sample_key in sample_parameters:
        concentration = float(sample_parameters[sample_key][3])
        sample_volume_needed = sample_mass/concentration
        sample_slot = sample_parameters[sample_key][0]
        sample_well = sample_parameters[sample_key][1]
        # sample_name = sample_parameters[sample_key][2]
        sample_dest_slot = sample_parameters[sample_key][4]
        sample_destination_labware = labware_dict[sample_dest_slot]
        sample_source_labware = labware_dict[sample_slot]

        # If user has replicates there can be more than one sample destination well.
        sample_dest_wells = sample_parameters[sample_key][5].split(",")

        # Apply some sanity to the sample volumes
        if sample_volume_needed > 20:
            sample_volume = round(sample_volume_needed, 1)
        else:
            sample_volume = round(sample_volume_needed, 2)

        water_volume = round((0.5*reaction_vol)-sample_volume, 2)

        # Define the pipettes for dispensing our samples and the water.
        sample_pipette, sample_loop, sample_volume = sample_pipette_selection(left_pipette, right_pipette, sample_volume)
        water_pipette, water_loop, water_volume = sample_pipette_selection(left_pipette, right_pipette, water_volume)

        # Add water to all the destination wells for this sample.
        water_pipette.pick_up_tip()
        for destination_well in sample_dest_wells:
            temp_water_loop = water_loop
            sample_dest_dict[sample_dest_slot].append(destination_well)

            # Dispense the water.
            while temp_water_loop > 0:
                # water_pipette.aspirate(water_volume, reagent_labware[args.WaterWell])
                water_pipette.aspirate(water_volume, reagent_labware[args.WaterWell].bottom(tip_height))
                water_pipette.dispense(water_volume, sample_destination_labware[destination_well])
                water_pipette.blow_out()
                aspirated_vol += water_volume
                incremental_aspirated_vol += water_volume

                # Adjust tip height from tube bottom.
                tip_height = res_tip_height(float(args.WaterResVol)-aspirated_vol, water_reservoir_dia, cone_vol)

                temp_water_loop -= 1

        # Change pipettes and get new tips if necessary.
        if not water_pipette == sample_pipette:
            water_pipette.drop_tip()
            sample_pipette.pick_up_tip()
        else:
            sample_pipette = water_pipette

        # Add template to all the destination wells for this sample.
        for destination_well in sample_dest_wells:
            temp_sample_loop = sample_loop
            while temp_sample_loop > 0:
                sample_pipette.aspirate(sample_volume, sample_source_labware[sample_well])
                sample_pipette.dispense(sample_volume, sample_destination_labware[destination_well])
                sample_pipette.blow_out()
                temp_sample_loop -= 1

        sample_pipette.drop_tip()

    return sample_dest_dict, aspirated_vol


if __name__ == "__main__":
    protocol_file = open('Generic_PCR.py')
    labware_path = "{}{}custom_labware".format(os.getcwd(), os.sep)
    run_log, __bundle__ = simulate(protocol_file, custom_labware_paths=[labware_path])
    print(format_runlog(run_log))
    protocol_file.close()
