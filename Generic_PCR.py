import csv
import datetime
import os
import math
import platform
from collections import defaultdict
from contextlib import suppress
from types import SimpleNamespace
from opentrons.simulate import simulate, format_runlog

# metadata
metadata = {
    'protocolName': 'Generic PCR v0.6.0',
    'author': 'Dennis Simpson',
    'description': 'Sets up a PCR from concentrated template',
    'apiLevel': '2.9'
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
    Calculate the the height of the liquid in a reservoir and return the value to set the pipette tip height.
    This works for both conical shapes and cylinders.
    @param res_vol:
    @param well_dia:
    @param cone_vol:
    @return:
    """
    if res_vol > cone_vol:
        cone_height = (3*cone_vol/(math.pi*((well_dia/2)**2)))
        height = ((res_vol-cone_vol)/(math.pi*((well_dia/2)**2)))-5+cone_height
    else:
        height = (3*res_vol/(math.pi*((well_dia/2)**2)))-3

    if height < 1:
        height = 0

    return int(height)


def labware_parsing(args, ctx):
    # Extract Slot information
    slot_list = ["Slot1", "Slot2", "Slot3", "Slot4", "Slot5", "Slot6", "Slot7", "Slot8", "Slot9", "Slot10", "Slot11"]
    labware_dict = {}
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
            labware_dict[str(i + 1)] = ctx.load_labware(labware, str(i + 1))
            if labware in tipbox_dict[args.LeftPipette]:
                left_tiprack_list.append(labware_dict[str(i + 1)])
            elif labware in tipbox_dict[args.RightPipette]:
                right_tiprack_list.append(labware_dict[str(i + 1)])

    return labware_dict, left_tiprack_list, right_tiprack_list


def run(ctx):

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

    sample_parameters, args = parse_sample_file(tsv_file_path)
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
    tip_height = res_tip_height(float(args.PCR_MixResVolume), pcr_reagent_well_dia, cone_vol)
    reagent_aspirated = 0

    for sample_key in sample_data_dict:
        sample_data = sample_data_dict[sample_key][0]
        reagent_dest_labware = labware_dict[sample_data[0]]
        reagent_dest_wells = sample_data[1]
        for reagent_dest_well in reagent_dest_wells:
            if not reagent_pipette.has_tip:
                reagent_pipette.pick_up_tip()

            dispensing_loop(args, reagent_loop, reagent_pipette,
                            reagent_labware[reagent_source_well].bottom(tip_height),
                            reagent_dest_labware[reagent_dest_well], pcr_reagent_vol, NewTip=True, MixReaction=True,
                            touch=True)

            reagent_aspirated += pcr_reagent_vol
            tip_height = \
                res_tip_height(float(args.PCR_MixResVolume) - reagent_aspirated, pcr_reagent_well_dia, cone_vol)

    # Setup the water control sample(s)
    if args.WaterControl:
        water_reservoir_dia = reagent_labware[args.WaterWell].diameter
        water_tip_height = \
            res_tip_height(float(args.WaterResVol) - aspirated_water_vol, water_reservoir_dia, cone_vol)

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


def labware_cone_volume(args, labware_name):

    cone_vol = 200
    labware = getattr(args, "Slot{}".format(str(labware_name)[-1:]))

    if "e5ml_" in labware:

        cone_vol = 1200

    elif"1.5ml" in labware:

        cone_vol = 450

    return cone_vol


def dispensing_loop(args, loop_count, pipette, source_location, destination_location, volume, NewTip, MixReaction, touch=False):
    """
    Generic function to dispense material into designated well.
    @param args:
    @param loop_count:
    @param pipette:
    @param source_location:
    @param destination_location:
    @param volume:
    @param NewTip:
    @param MixReaction:
    @param touch:
    @return:
    """
    def tip_touch():
        pipette.touch_tip(radius=0.75, v_offset=-8)

    if NewTip:
        if pipette.has_tip:
            pipette.drop_tip()
        pipette.pick_up_tip()

    while loop_count > 0:
        pipette.aspirate(volume, source_location)
        tip_touch()
        pipette.dispense(volume, destination_location)

        loop_count -= 1
        if not MixReaction:
            pipette.blow_out()
            if touch:
                tip_touch()

    if MixReaction:
        pipette.mix(repetitions=4, volume=float(args.PCR_Volume)*0.7, rate=5.0)
        tip_touch()

    if NewTip:
        pipette.drop_tip()

    return pipette


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
    tip_height = res_tip_height(float(args.WaterResVol), water_reservoir_dia, cone_vol)
    aspirated_water_vol = 0

    # Figure out volumes and dispense water into each destination well.
    sample_data_dict = defaultdict(list)
    for sample_key in sample_parameters:
        concentration = float(sample_parameters[sample_key][3])
        sample_volume_needed = sample_mass/concentration
        sample_dest_slot = sample_parameters[sample_key][4]
        sample_destination_labware = labware_dict[sample_dest_slot]

        # If user has replicates there can be more than one sample destination well.
        sample_dest_wells = sample_parameters[sample_key][5].split(",")

        # Apply some sanity to the sample volumes
        if sample_volume_needed > 20:
            sample_volume = round(sample_volume_needed, 1)
        else:
            sample_volume = round(sample_volume_needed, 2)

        water_volume = round((0.5*reaction_vol)-sample_volume, 2)
        sample_data_dict[sample_key].append((sample_dest_slot, sample_dest_wells, sample_volume))

        # Define the pipettes for dispensing the water.
        water_pipette, water_loop, water_volume = pipette_selection(left_pipette, right_pipette, water_volume)

        # Add water to all the destination wells for this sample.
        for sample_dest_well in sample_dest_wells:
            if not water_pipette.has_tip:
                water_pipette.pick_up_tip()

            dispensing_loop(args, water_loop, water_pipette, reagent_labware[args.WaterWell].bottom(tip_height),
                            sample_destination_labware[sample_dest_well], water_volume, NewTip=False, MixReaction=False)

            aspirated_water_vol += water_volume
            tip_height = res_tip_height(float(args.WaterResVol) - aspirated_water_vol, water_reservoir_dia, cone_vol)

    # Drop any tips the pipettes might have.
    left_pipette.drop_tip()
    right_pipette.drop_tip()

    # Now dispense samples into PCR wells.
    for sample_key in sample_data_dict:
        sample_source_well = sample_parameters[sample_key][1]
        sample_source_labware = labware_dict[sample_parameters[sample_key][0]]
        sample_data = sample_data_dict[sample_key][0]
        sample_volume = sample_data[2]
        sample_dest_wells = sample_data[1]
        sample_dest_labware = labware_dict[sample_data[0]]

        # Define the pipettes for dispensing the samples.
        sample_pipette, sample_loop, sample_volume = pipette_selection(left_pipette, right_pipette, sample_volume)

        # Add template to all the destination wells for this sample.
        for sample_dest_well in sample_dest_wells:

            if not sample_pipette.has_tip:
                sample_pipette.pick_up_tip()

            dispensing_loop(args, sample_loop, sample_pipette, sample_source_labware[sample_source_well],
                            sample_dest_labware[sample_dest_well], sample_volume, NewTip=True, MixReaction=False)

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
