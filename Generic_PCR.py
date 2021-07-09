import datetime
import os
import sys
import platform
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
    'protocolName': 'Generic PCR v0.8.1',
    'author': 'Dennis Simpson <dennis@email.unc.edu>',
    'description': 'Sets up a PCR from concentrated template',
    'apiLevel': '2.10'
}


def run(ctx: protocol_api.ProtocolContext):

    # Turn on rail lights and pause program so user can load robot deck.
    # ctx.set_rail_lights(True)
    # ctx.pause("Load Labware onto robot deck and click resume when ready to continue")
    # ctx.home()
    ctx.set_rail_lights(False)

    # TSV file location on OT-2
    tsv_file_path = "{0}var{0}lib{0}jupyter{0}notebooks{0}ProcedureFile.tsv".format(os.sep)
    if not os.path.isfile(tsv_file_path):
        # Temp TSV file location on Win10 Computers for simulation
        tsv_file_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())

    sample_parameters, args = Utilities.parse_sample_template(tsv_file_path)
    labware_dict, left_tiprack_list, right_tiprack_list = Utilities.labware_parsing(args, ctx)

    # Pipettes
    left_pipette = ctx.load_instrument(args.LeftPipette, 'left', tip_racks=left_tiprack_list)
    right_pipette = ctx.load_instrument(args.RightPipette, 'right', tip_racks=right_tiprack_list)

    # Set the location of the first tip in box.
    left_pipette.starting_tip = left_tiprack_list[0].wells_by_name()[args.LeftPipetteFirstTip]
    right_pipette.starting_tip = right_tiprack_list[0].wells_by_name()[args.RightPipetteFirstTip]

    # Make sample dilutions.  Calculate sample and water volumes.
    sample_data_dict, aspirated_water_vol = \
        process_samples(args, ctx, sample_parameters, labware_dict, left_pipette, right_pipette)

    # Dispense Water
    aspirated_water_vol = \
        dispense_water(args, sample_data_dict, labware_dict, left_pipette, right_pipette, aspirated_water_vol)

    # Dispense Samples
    aspirated_water_vol = \
        dispense_samples(args, sample_data_dict, left_pipette, right_pipette, aspirated_water_vol)

    # Dispense PCR Reagents.
    dispense_pcr_reagents(args, labware_dict, left_pipette, right_pipette, aspirated_water_vol, sample_data_dict)

    if not ctx.is_simulating():
        os.remove(tsv_file_path)


def process_samples(args, ctx, sample_parameters, labware_dict, left_pipette, right_pipette):
    """

    @param args:
    @param ctx:
    @param sample_parameters:
    @param labware_dict:
    @param left_pipette:
    @param right_pipette:
    @return:
    """
    sample_mass = float(args.DNA_in_Reaction)
    reaction_vol = float(args.PCR_Volume)
    reagent_labware = labware_dict[args.ReagentSlot]
    dilution_slot = getattr(args, "DilutionPlateSlot")
    if dilution_slot:
        dilution_labware = labware_dict[args.DilutionPlateSlot]
    water_res_well_dia = reagent_labware[args.WaterResWell].diameter
    cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
    water_tip_height = Utilities.res_tip_height(float(args.WaterResVol), water_res_well_dia, cone_vol, float(args.BottomOffset))
    aspirated_water_vol = 0
    dilution_well_index = 0
    water_aspirated = 0
    dilution_plate_layout = Utilities.plate_layout()
    output_string = "Well\tSample\tDilution\tDispensed Volume\n"

    # Figure out volumes and dispense water into each destination well.
    sample_data_dict = defaultdict(list)
    for sample_key in sample_parameters:
        sample_concentration = float(sample_parameters[sample_key][3])
        sample_name = sample_parameters[sample_key][2]
        sample_volume = round(sample_mass/sample_concentration, 2)
        sample_dest_slot = sample_parameters[sample_key][4]
        sample_destination_labware = labware_dict[sample_dest_slot]
        sample_source_labware = labware_dict[sample_parameters[sample_key][0]]
        sample_source_well = sample_parameters[sample_key][1]
        sample_source = sample_source_labware[sample_source_well]
        dilution = "Neat"

        # Calculate any dilutions required.
        sample_vol, diluent_vol, diluted_sample_vol, reaction_water_vol, max_template_vol = \
            Utilities.calculate_volumes(args, sample_concentration)

        dispensed_sample = sample_vol
        diluted_sample_vol = round(diluted_sample_vol, 2)

        # If user has replicates there can be more than one sample destination well.
        sample_dest_wells = sample_parameters[sample_key][5].split(",")

        # Dilute samples, if required, and collect information for dispensing.
        if diluent_vol > 0:
            dispensed_sample = diluted_sample_vol

            # Adjust volume of diluted sample to make sure there is enough
            diluted_template_needed = diluted_sample_vol * (len(sample_dest_wells) + 1)
            diluted_template_on_hand = sample_vol + diluent_vol
            diluted_template_factor = 1.0
            if diluted_template_needed <= (diluted_template_on_hand+5):
                diluted_template_factor = diluted_template_on_hand/diluted_template_needed
                if diluted_template_factor <= 1.5 and (sample_vol * diluted_template_factor) < 10:
                    diluted_template_factor = 2.2

            adjusted_sample_vol = sample_vol * diluted_template_factor
            diluent_vol = diluent_vol * diluted_template_factor
            dilution = "1:{}".format(int((adjusted_sample_vol + diluent_vol)/adjusted_sample_vol))
            # Reset the pipettes for the new volumes
            diluent_pipette, diluent_loop, diluent_vol = \
                Utilities.pipette_selection(left_pipette, right_pipette, diluent_vol)
            sample_pipette, sample_loop, sample_vol = \
                Utilities.pipette_selection(left_pipette, right_pipette, adjusted_sample_vol)

            # Make dilution, diluent first
            if dilution_slot:
                dilution_well = dilution_plate_layout[dilution_well_index]
                Utilities.dispensing_loop(args, diluent_loop, diluent_pipette,
                                          reagent_labware[args.WaterResWell].bottom(water_tip_height),
                                          dilution_labware[dilution_well], diluent_vol, NewTip=False, MixReaction=False,
                                          touch=True)

                mvol = sample_vol+diluent_vol
                Utilities.dispensing_loop(args, sample_loop, sample_pipette, sample_source_labware[sample_source_well],
                                          dilution_labware[dilution_well], sample_vol, NewTip=False, MixReaction=True,
                                          touch=True, MixVolume=mvol)

                water_aspirated += diluent_vol
                dilution_well_index += 1
                sample_source = dilution_labware[dilution_well]
                water_tip_height = \
                    Utilities.res_tip_height(float(args.WaterResVol) - water_aspirated, water_res_well_dia,
                                             cone_vol, float(args.BottomOffset))

            if diluent_pipette.has_tip:
                diluent_pipette.drop_tip()
            if sample_pipette.has_tip:
                sample_pipette.drop_tip()

        for well in sample_dest_wells:
            dispensed_sample = round(dispensed_sample, 1)
            output_string += "{}\t{}\t{}\t{}\n".format(well, sample_name, dilution, dispensed_sample)
            sample_data_dict[well] = [sample_source, dispensed_sample, reaction_water_vol, sample_destination_labware]

    # This will output a plate layout file.  Only does it during the simulation from our GUI
    if ctx.is_simulating() and platform.system() == "Windows":
        run_date = datetime.datetime.today().strftime("%a %b %d %H:%M %Y")
        header = "Sample Information for Generic PCR\nSimulation Date:\t{}\tTemplate User:\t{}\n\n{}"\
            .format(run_date, args.User, output_string)
        plate_layout_file = \
            open("C:{0}Users{0}{1}{0}Documents{0}GenericPCR_SampleData.tsv".format(os.sep, os.getlogin()), 'w')
        plate_layout_file.write(header)
        plate_layout_file.close()

    return sample_data_dict, water_aspirated


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
        Utilities.pipette_selection(left_pipette, right_pipette, float(args.PCR_Volume) * 0.5)

    reagent_labware = labware_dict[args.ReagentSlot]
    reagent_source_well = args.PCR_MixWell
    pcr_reagent_well_dia = reagent_labware[args.PCR_MixWell].diameter
    cone_vol = Utilities.labware_cone_volume(args, reagent_labware)
    tip_height = \
        Utilities.res_tip_height(float(args.PCR_MixResVolume), pcr_reagent_well_dia, cone_vol, float(args.BottomOffset))
    reagent_aspirated = 0

    for well in sample_data_dict:
        reagent_dest_labware = sample_data_dict[well][3]
        Utilities.dispensing_loop(args, reagent_loop, reagent_pipette,
                                  reagent_labware[reagent_source_well].bottom(tip_height),
                                  reagent_dest_labware[well], pcr_reagent_vol, NewTip=True,
                                  MixReaction=True, touch=True)
        reagent_aspirated += pcr_reagent_vol
        tip_height = \
            Utilities.res_tip_height(float(args.PCR_MixResVolume) - reagent_aspirated, pcr_reagent_well_dia, cone_vol,
                                     float(args.BottomOffset))

    # Setup the water control sample(s)
    if args.WaterControl:
        water_reservoir_dia = reagent_labware[args.WaterResWell].diameter
        water_tip_height = \
            Utilities.res_tip_height(float(args.WaterResVol) - aspirated_water_vol, water_reservoir_dia, cone_vol,
                                     float(args.BottomOffset))

        slot = args.WaterControl.split(',')[0]
        well = args.WaterControl.split(',')[1]
        neg_control_labware = labware_dict[slot]
        reagent_pipette, reagent_loop, pcr_reagent_vol = \
            Utilities.pipette_selection(left_pipette, right_pipette, float(args.PCR_Volume) * 0.5)
        reagent_pipette.pick_up_tip()
        temp_loop = reagent_loop

        while temp_loop > 0:
            # Dispense Water into negative control well
            reagent_pipette.aspirate(pcr_reagent_vol, reagent_labware[args.WaterResWell].bottom(water_tip_height))
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


def dispense_water(args, sample_data_dict, labware_dict, left_pipette, right_pipette, aspirated_water_vol):
    """
    Dispense water into wells that need it using a single tip.

    @param args:
    @param sample_data_dict:
    @param labware_dict:
    @param left_pipette:
    @param right_pipette:
    @param aspirated_water_vol:
    @return:
    """
    reagent_labware = labware_dict[args.ReagentSlot]
    water_res_well_dia = reagent_labware[args.WaterResWell].diameter
    cone_vol = Utilities.labware_cone_volume(args, reagent_labware)

    for well in sample_data_dict:
        water_destination = sample_data_dict[well][3]
        water_aspirated = sample_data_dict[well][2]
        water_pipette, water_loop, water_dispensed = \
            Utilities.pipette_selection(left_pipette, right_pipette, water_aspirated)

        water_tip_height = \
            Utilities.res_tip_height(float(args.WaterResVol) - aspirated_water_vol, water_res_well_dia, cone_vol,
                                     float(args.BottomOffset))

        Utilities.dispensing_loop(args, water_loop, water_pipette,
                                  reagent_labware[args.WaterResWell].bottom(water_tip_height),
                                  water_destination[well], water_dispensed, NewTip=False, MixReaction=False, touch=True)

        aspirated_water_vol += water_dispensed

    # Drop any tips the pipettes might have.
    if left_pipette.has_tip:
        left_pipette.drop_tip()
    if right_pipette.has_tip:
        right_pipette.drop_tip()

    return aspirated_water_vol


def dispense_samples(args, sample_data_dict, left_pipette, right_pipette, aspirated_water_vol):
    """
    Add water and template to the destination wells for the PCR.
    @param args:
    @param sample_data_dict:
    @param left_pipette:
    @param right_pipette:
    @param aspirated_water_vol:
    @return:
    """

    for well in sample_data_dict:
        sample_destination = sample_data_dict[well][3]
        sample_source = sample_data_dict[well][0]
        dispensed_sample = sample_data_dict[well][1]

        sample_pipette, sample_loop, sample_dispensed = \
            Utilities.pipette_selection(left_pipette, right_pipette, dispensed_sample)

        Utilities.dispensing_loop(args, sample_loop, sample_pipette, sample_source, sample_destination[well],
                                  sample_dispensed, NewTip=True, MixReaction=False, touch=True)

    # Drop any tips the pipettes might have.
    if left_pipette.has_tip:
        left_pipette.drop_tip()
    if right_pipette.has_tip:
        right_pipette.drop_tip()

    return aspirated_water_vol


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
