import csv
import datetime
import math
import os
import platform
from collections import defaultdict
from types import SimpleNamespace
from opentrons.simulate import simulate, format_runlog
# import opentrons

# metadata
metadata = {
    'protocolName': 'ddPCR v0.2.0',
    'author': 'Dennis Simpson',
    'description': 'Setup a ddPCR using either 2x or 4x SuperMix',
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
        cone_height = (3*cone_vol / (math.pi * ((well_dia / 2) ** 2)))
        height = ((res_vol-cone_vol)/(math.pi*((well_dia/2)**2)))-6+cone_height
    else:
        height = (3*res_vol / (math.pi * ((well_dia / 2) ** 2))) - 4

    if height <= 6:
        height = 1

    return int(height)


def labware_cone_volume(args, labware_name):

    cone_vol = 200
    labware = getattr(args, "Slot{}".format(str(labware_name)[-1:]))

    if "5ml_" in labware:

        cone_vol = 1200

    elif"1.5ml" in labware:

        cone_vol = 500

    return cone_vol


def plate_layout():
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    plate_layout_by_column = []
    for i in range(12):
        for row in rows:
            plate_layout_by_column.append("{}{}".format(row, i+1))
    return plate_layout_by_column


def dispensing_loop(args, loop_count, pipette, source_location, destination_location, volume, NewTip, MixReaction):
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
    @return:
    """
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
    target_concentration = template_in_reaction/2
    max_template_vol = \
        (float(args.PCR_Volume) / int(args.PCR_MixConcentration.split("x")[0])) - float(args.TargetVolume)

    # There is the potential to need the max template volume for the lower dilutions.  This defines those.
    dilution_template = [(4, 4), (2, 6), (2, 10), (1, 7), (1, 9), (1, 11)]
    if args.PCR_MixConcentration == "4x":
        dilution_template = [(7, 7), (4, 12), (3, 15), (2, 14), (2, 18), (2, 22)]

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
            if i < 6:
                dilution_data = dilution_template[i]
            else:
                dilution_data = (1, dilution-1)
            diluted_sample_vol = round(template_in_reaction/diluted_dna_conc, 2)
            reaction_water_vol = max_template_vol-diluted_sample_vol

            return dilution_data[0], dilution_data[1], diluted_sample_vol, reaction_water_vol, max_template_vol


def sample_processing(args, sample_parameters):
    sample_data_dict = defaultdict(list)
    target_well_dict = defaultdict(list)
    water_well_dict = defaultdict(float)
    layout_data = defaultdict(list)
    all_data_dict = defaultdict(list)

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

        # Add no primer control well for sample and increment well count
        no_primer_well = plate_layout_by_column[dest_well_count]
        used_wells.append(no_primer_well)
        sample_wells.append(no_primer_well)
        water_vol = float(args.TargetVolume)+(max_template_vol-diluted_sample_vol)
        water_well_dict[no_primer_well] = water_vol
        sample_data_dict[sample_key] = [sample_vol, diluent_vol, diluted_sample_vol, sample_wells]
        water_row = no_primer_well[0]
        water_column = int(no_primer_well[1:])-1
        layout_data[water_row][water_column] = "{}|No Primers|{}|{}".format(sample_string, dilution, diluted_sample_vol)
        dest_well_count += 1

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

    return sample_data_dict, water_well_dict, target_well_dict, used_wells, layout_data


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

    sample_data_dict, water_well_dict, target_well_dict, used_wells, layout_data = \
        sample_processing(args, sample_parameters)

    # This will output a plate layout file.  Only does it during the simulation from our GUI
    if ctx.is_simulating():
        # if not os.path.isdir("C:{0}Users{0}{1}{0}Documents".format(os.sep, os.getlogin())):
        if platform.system() == "Windows":
            run_date = datetime.datetime.today().strftime("%a %b %d %H:%M %Y")

            header = "## ddPCR Setup\n## Setup Date:\t{}\n## Template User:\t{}\n" \
                     "# Format:\tTemplate | Target | Template Dilution | Template Volume in Reaction\n\n\t"\
                .format(run_date, args.User)

            for i in range(12):
                header += "{}\t".format(i+1)
            import natsort
            outstring = ""
            for well in natsort.natsorted(layout_data):
                well_string = "\t".join(layout_data[well])
                outstring += "\n{}\t{}\t".format(well, well_string)
            outfile = open("C:{0}Users{0}{1}{0}Documents{0}PlateLayout.csv".format(os.sep, os.getlogin()), 'w')
            outfile.write(header+outstring)
            outfile.close()

    water_aspirated = dispense_water(args, labware_dict, water_well_dict, left_pipette, right_pipette)

    positive_control_dict = dispense_primer_mix(args, labware_dict, target_well_dict, left_pipette, right_pipette)
    dispense_samples(args, labware_dict, sample_data_dict, sample_parameters, left_pipette, right_pipette, water_aspirated)
    dispense_positive_controls(args, positive_control_dict, labware_dict, left_pipette, right_pipette)
    dispense_supermix(args, labware_dict, left_pipette, right_pipette, used_wells)

    if not ctx.is_simulating():
        os.remove(tsv_file_path)


def dispense_primer_mix(args, labware_dict, target_well_dict, left_pipette, right_pipette):
    """

    @param args:
    @param labware_dict:
    @param target_well_dict:
    @param left_pipette:
    @param right_pipette:
    @return:
    """
    target_pipette, target_loop, target_volume = pipette_selection(left_pipette, right_pipette, float(args.TargetVolume))
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

        # The last well in the list is the positive control for the target primers
        pos_control_well = target_well_list[-1]
        pos_control_info = getattr(args, "PositiveControl_{}".format(target))
        positive_control_dict[pos_control_well] = pos_control_info

        for well in target_well_list:

            dispensing_loop(args, target_loop, target_pipette, target_source_labware[target_source_well],
                            sample_destination_labware[well], target_volume, NewTip=False, MixReaction=False)

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
    water_tip_height = res_tip_height(float(args.WaterResVol), water_res_well_dia, cone_vol)
    sample_destination_labware = labware_dict[args.PCR_PlateSlot]
    water_aspirated = 0

    for well in water_well_dict:
        water_volume = water_well_dict[well]

        # There is a bug in the Opentrons software that results in a 0 uL aspiration defaulting to pipette max.
        if water_volume < 0.4:
            continue

        # Define the pipette for dispensing the water.
        water_pipette, water_loop, water_volume = pipette_selection(left_pipette, right_pipette, water_volume)

        if not water_pipette.has_tip:
            water_pipette.pick_up_tip()

        dispensing_loop(args, water_loop, water_pipette, reagent_labware[args.WaterWell].bottom(water_tip_height),
                        sample_destination_labware[well], water_volume, NewTip=False, MixReaction=False)
        water_aspirated += water_volume
        water_tip_height = res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol)

    if left_pipette.has_tip:
        left_pipette.drop_tip()
    if right_pipette.has_tip:
        right_pipette.drop_tip()

    return water_aspirated


def dispense_samples(args, labware_dict, sample_data_dict, sample_parameters, left_pipette, right_pipette, water_aspirated):
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
    water_tip_height = res_tip_height(float(args.WaterResVol), water_res_well_dia, cone_vol)
    dilution_plate_layout = plate_layout()
    dilution_well_index = 0

    for sample_key in sample_parameters:
        sample_source_labware = labware_dict[sample_parameters[sample_key][0]]
        sample_source_well = sample_parameters[sample_key][1]
        sample_dest_wells = sample_data_dict[sample_key][3]

        # Remove no template control well
        sample_dest_wells = sample_dest_wells[:-1]
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
        volume_ratio = (undiluted_sample_vol + diluent_vol) / (len(sample_dest_wells) * diluted_sample_vol)
        if volume_ratio < 1:
            undiluted_sample_vol = undiluted_sample_vol / volume_ratio
            diluent_vol = diluent_vol / volume_ratio

        # Reset the pipettes for the new volumes
        diluent_pipette, diluent_loop, diluent_vol = pipette_selection(left_pipette, right_pipette, diluent_vol)
        sample_pipette, sample_loop, undiluted_sample_vol = \
            pipette_selection(left_pipette, right_pipette, undiluted_sample_vol)

        # Make dilution, diluent first
        dilution_well = dilution_plate_layout[dilution_well_index]
        dispensing_loop(args, diluent_loop, diluent_pipette, reagent_labware[args.WaterWell].bottom(water_tip_height),
                        dilution_labware[dilution_well], diluent_vol, NewTip=True, MixReaction=False)
        dispensing_loop(args, sample_loop, sample_pipette, sample_source_labware[sample_source_well],
                        dilution_labware[dilution_well], diluent_vol, NewTip=True, MixReaction=True)
        water_aspirated += diluent_vol
        water_tip_height = res_tip_height(float(args.WaterResVol)-water_aspirated, water_res_well_dia, cone_vol)

        # Add diluted sample to PCR plate
        for well in sample_dest_wells:
            sample_pipette, sample_loop, diluted_sample_vol = \
                pipette_selection(left_pipette, right_pipette, diluted_sample_vol)

            dispensing_loop(args, sample_loop, sample_pipette, dilution_labware[dilution_well],
                            sample_destination_labware[well], diluted_sample_vol, NewTip=True, MixReaction=False)
        if sample_pipette.has_tip:
            sample_pipette.drop_tip()


def dispense_positive_controls(args, positive_control_dict, labware_dict, left_pipette, right_pipette):
    """
    Dispense positive controls to appropriate wells
    @param args:
    @param positive_control_dict:
    @param labware_dict:
    @param left_pipette:
    @param right_pipette:
    """

    sample_destination_labware = labware_dict[args.PCR_PlateSlot]

    for destination_well in positive_control_dict:
        positive_control_info = positive_control_dict[destination_well].split(",")
        control_sample_vol = float(args.TargetVolume) + (
                float(args.PCR_Volume) / int(args.PCR_MixConcentration.split("x")[0]))

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
    supermix_tip_height = res_tip_height(float(args.PCR_MixResVolume), supermix_well_dia, cone_vol)
    supermix_aspirated = 0

    for well in used_wells:
        dispensing_loop(args, supermix_loop, supermix_pipette,
                        reagent_labware[supermix_source_well].bottom(supermix_tip_height),
                        sample_destination_labware[well], supermix_vol,
                        NewTip=True, MixReaction=True)

        supermix_aspirated += supermix_vol
        supermix_tip_height = res_tip_height(float(args.PCR_MixResVolume)-supermix_aspirated, supermix_well_dia, cone_vol)


if __name__ == "__main__":
    protocol_file = open('ddPCR.py')
    labware_path = "{}{}custom_labware".format(os.getcwd(), os.sep)
    run_log, __bundle__ = simulate(protocol_file, custom_labware_paths=[labware_path])
    print(format_runlog(run_log))
    protocol_file.close()
