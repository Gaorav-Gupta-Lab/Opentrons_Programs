"""
Functions common for all robot programs.

Dennis Simpson
University of North Carolina at Chapel Hill
Chapel Hill NC, 27599

@copyright 2022

"""
import math

__version__ = "1.0.1"


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
        cone_vol = 450

    return cone_vol


def res_tip_height(res_vol, well_dia, cone_vol, bottom_offset):
    """
    Calculate the the height of the liquid in a reservoir and return the value to set the pipette tip height.
    This works for both conical shapes and cylinders.
    @param bottom_offset:
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

    if height < 3:
        height = bottom_offset

    return round(height, 1)


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

    return pipette, loop, round(volume, 1)


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


def calculate_volumes(args, sample_concentration, template_in_rxn):
    """
    Calculates volumes for dilution and distribution of sample.
    Returns a list of tuples consisting of
    (uL of sample to dilute, uL of water for dilution), (uL of diluted sample in reaction, uL of water in reaction)

    :param args:
    :param sample_concentration:
    :param template_in_rxn:
    :return:

    """

    max_template_vol = round(float(args.PCR_Volume)-float(args.MasterMixPerRxn), 1)

    # If at least 2 uL of sample is needed then no dilution is necessary
    if template_in_rxn/sample_concentration >= 2:
        sample_vol = round(template_in_rxn/sample_concentration, 2)
        return sample_vol, 0, 0, max_template_vol-sample_vol, max_template_vol

    # This will test a series of dilutions up to a 1:200.
    for i in range(50):
        dilution = (i+1)*2
        diluted_dna_conc = sample_concentration/dilution

        # Want to pipette at least 2 uL of diluted sample per well
        if 2 <= template_in_rxn/diluted_dna_conc <= max_template_vol:
            diluted_sample_vol = round(template_in_rxn/diluted_dna_conc, 2)
            reaction_water_vol = max_template_vol-diluted_sample_vol

            return 1, dilution - 1, diluted_sample_vol, reaction_water_vol, max_template_vol


def dispensing_loop(args, loop_count, pipette, source_location, destination_location, volume, NewTip, MixReaction,
                    touch=False, MixVolume=None):
    """
    Generic function to dispense material into designated well.
    @param MixVolume:
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

    if not pipette.has_tip:
        pipette.pick_up_tip()

    while loop_count > 0:
        pipette.aspirate(volume, source_location, rate=0.75)

        if touch:
            tip_touch()

        pipette.dispense(volume, destination_location, rate=0.75)
        loop_count -= 1

        if not MixReaction:
            pipette.blow_out()
            if touch:
                tip_touch()

    if MixReaction:
        v = float(args.PCR_Volume)
        if MixVolume:
            v = MixVolume
        pipette.mix(repetitions=4, volume=v*0.65, rate=2.0)
        pipette.blow_out()
        tip_touch()

    if NewTip:
        pipette.drop_tip()

    return pipette
