import csv
import datetime
import math
import os
import platform
from contextlib import suppress
from collections import defaultdict
from types import SimpleNamespace
from opentrons.simulate import simulate, format_runlog

# metadata
metadata = {
    'protocolName': 'ddPCR v0.4.6',
    'author': 'Dennis Simpson',
    'description': 'Setup a ddPCR using either 2x or 4x SuperMix',
    'apiLevel': '2.8'
}


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
        height = ((res_vol-cone_vol)/(math.pi*((well_dia/2)**2)))-5+cone_height
    else:
        height = (3*res_vol / (math.pi * ((well_dia / 2) ** 2))) - 5

    if height <= 6:
        height = 1

    return int(height)


def dispensing_loop(loop_count, pipette, source_location, destination_location, volume, NewTip, MixReaction,
                    touch=False):
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
        pipette.mix(repetitions=4, volume=25*0.7, rate=5.0)
        pipette.blow_out()
        tip_touch()

    if NewTip:
        pipette.drop_tip()

    return pipette


def run(ctx):
    ctx.comment("Begin {}".format(metadata['protocolName']))

    # Turn on rail lights and pause program so user can load robot deck.
    ctx.set_rail_lights(True)
    right_tipracks = [ctx.load_labware('opentrons_96_tiprack_300ul', '3')]
    water_pipette = ctx.load_instrument('p300_single_gen2', 'right', tip_racks=right_tipracks)
    print(water_pipette)
    reagent_labware = ctx.load_labware('vwrmicrocentrifugetube1.5ml_24_tuberack_1500ul', '2')
    cone_vol = 500
    initial_water_volume = 1500
    dispensed_vol = 25
    source_well = "D1"
    destination_well = "D3"
    water_res_well_dia = reagent_labware["A1"].diameter
    water_tip_height = res_tip_height(initial_water_volume, water_res_well_dia, cone_vol)
    water_aspirated = 0
    count = 0

    while count < 11:
        ctx.comment(str(water_tip_height))
        if not water_pipette.has_tip:
            water_pipette.pick_up_tip()

        dispensing_loop(1, water_pipette, reagent_labware[source_well].bottom(water_tip_height),
                        reagent_labware[destination_well], dispensed_vol, NewTip=False, MixReaction=False)

        water_aspirated += dispensed_vol
        water_tip_height = res_tip_height(initial_water_volume-water_aspirated, water_res_well_dia, cone_vol)

        count += 1

    if water_pipette.has_tip:
        water_pipette.drop_tip()


if __name__ == "__main__":
    protocol_file = open('Dispensing_Test.py')
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
