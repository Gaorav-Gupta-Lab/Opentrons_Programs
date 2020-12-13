import csv
import os
import re
from collections import namedtuple, defaultdict
from types import SimpleNamespace
from opentrons import protocol_api
import opentrons
import Tool_Box


def get_values(*names):
    import json
    _all_values = json.loads("""{"transfer_csv":"Source Labware,Source Slot,Source Well,Source Aspiration Height Above Bottom (in mm),Dest Labware,Dest Slot,Dest Well,Volume (in ul)\\nagilent_1_reservoir_290ml,1,A1,1,nest_96_wellplate_100ul_pcr_full_skirt,4,A11,1\\nnest_12_reservoir_15ml,2,A1,1,nest_96_wellplate_2ml_deep,5,A5,3\\nnest_1_reservoir_195ml,3,A1,1,nest_96_wellplate_2ml_deep,5,H12,7","pipette_type":"p10_single","pipette_mount":"right","tip_type":"standard","tip_reuse":"always"}""")
    return [_all_values[n] for n in names]


metadata = {
    'protocolName': 'Cherrypicking',
    'author': 'Nick <protocols@opentrons.com>',
    'source': 'Custom Protocol Request',
    'apiLevel': '2.3'
}


def parse_sample_file(input_file):
    """
    Parse the input file and return a list of lists.
    :return:
    """

    if not os.path.isfile(input_file):
        raise SystemExit("{} Not Found.  Check File Name and Path.".format(input_file))

    index_list = []
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


def test(protocol: protocol_api.ProtocolContext):
    tsv = "D:{0}Working{0}MyTest.tsv".format(os.sep)
    [pipette_type, pipette_mount, tip_type,
     tip_reuse, transfer_csv] = get_values(  # noqa: F821
        "pipette_type", "pipette_mount", "tip_type", "tip_reuse",
        "transfer_csv")
    transfer_info = [[val.strip().lower() for val in line.split(',')]
                     for line in transfer_csv.splitlines()
                     if line.split(',')[0].strip()][1:]

    # load sample data
    sample_parameters, arg = parse_sample_file(tsv)
    sample_mass = float(arg.DNA_in_Reaction)
    reaction_vol = float(arg.PCR_Volume)
    left_pipette = arg.LeftPipette
    for sample_key in sample_parameters:
        print(sample_parameters[sample_key])
        concentration = float(sample_parameters[sample_key][3])
        volume = sample_mass/concentration
        sample_slot = sample_parameters[sample_key][0]
        if volume > 20:
            dna_volume = round(volume, 1)
        else:
            dna_volume = round(volume, 2)
        water_volume = (0.5*reaction_vol)-dna_volume
        print(dna_volume, water_volume)

    for line in transfer_info:

        s_lw, s_slot, d_lw, d_slot = line[:2] + line[4:6]
        print(line, s_lw, s_slot, d_lw, d_slot)
        for slot, lw in zip([s_slot, d_slot], [s_lw, d_lw]):
            print(slot, lw)
            if not int(slot) in protocol.loaded_labwares:
                protocol.load_labware(lw.lower(), slot)


def run(ctx):

    [pipette_type, pipette_mount, tip_type,
     tip_reuse, transfer_csv] = get_values(  # noqa: F821
        "pipette_type", "pipette_mount", "tip_type", "tip_reuse",
        "transfer_csv")

    tiprack_map = {
        'p10_single': {
            'standard': 'opentrons_96_tiprack_10ul',
            'filter': 'opentrons_96_filtertiprack_20ul'
        },
        'p50_single': {
            'standard': 'opentrons_96_tiprack_300ul',
            'filter': 'opentrons_96_filtertiprack_200ul'
        },
        'p300_single_gen1': {
            'standard': 'opentrons_96_tiprack_300ul',
            'filter': 'opentrons_96_filtertiprack_200ul'
        },
        'p1000_single_gen1': {
            'standard': 'opentrons_96_tiprack_1000ul',
            'filter': 'opentrons_96_filtertiprack_1000ul'
        },
        'p20_single_gen2': {
            'standard': 'opentrons_96_tiprack_20ul',
            'filter': 'opentrons_96_filtertiprack_20ul'
        },
        'p300_single_gen2': {
            'standard': 'opentrons_96_tiprack_300ul',
            'filter': 'opentrons_96_filtertiprack_200ul'
        },
        'p1000_single_gen2': {
            'standard': 'opentrons_96_tiprack_1000ul',
            'filter': 'opentrons_96_filtertiprack_1000ul'
        }
    }

    # load labware
    transfer_info = [[val.strip().lower() for val in line.split(',')]
                     for line in transfer_csv.splitlines()
                     if line.split(',')[0].strip()][1:]
    for line in transfer_info:
        s_lw, s_slot, d_lw, d_slot = line[:2] + line[4:6]
        for slot, lw in zip([s_slot, d_slot], [s_lw, d_lw]):
            if not int(slot) in ctx.loaded_labwares:
                ctx.load_labware(lw.lower(), slot)

    # load tipracks in remaining slots
    tiprack_type = tiprack_map[pipette_type][tip_type]
    tipracks = []
    for slot in range(1, 13):
        if slot not in ctx.loaded_labwares:
            tipracks.append(ctx.load_labware(tiprack_type, str(slot)))

    # load pipette
    pip = ctx.load_instrument(pipette_type, pipette_mount, tip_racks=tipracks)

    tip_count = 0
    tip_max = len(tipracks*96)

    def pick_up():
        nonlocal tip_count
        if tip_count == tip_max:
            ctx.pause('Please refill tipracks before resuming.')
            pip.reset_tipracks()
            tip_count = 0
        pip.pick_up_tip()
        tip_count += 1

    def parse_well(well):
        letter = well[0]
        number = well[1:]
        return letter.upper() + str(int(number))

    if tip_reuse == 'never':
        pick_up()
    for line in transfer_info:
        _, s_slot, s_well, h, _, d_slot, d_well, vol = line[:8]
        source = ctx.loaded_labwares[
            int(s_slot)].wells_by_name()[parse_well(s_well)].bottom(float(h))
        dest = ctx.loaded_labwares[
            int(d_slot)].wells_by_name()[parse_well(d_well)]
        if tip_reuse == 'always':
            pick_up()
        pip.transfer(float(vol), source, dest, new_tip='never')
        if tip_reuse == 'always':
            pip.drop_tip()
    if pip.hw_pipette['has_tip']:
        pip.drop_tip()

if __name__ == "__main__":
    test(opentrons.config)