"""
This is to run a simulation of PCR.py
"""
import datetime
import os
import platform
from opentrons.simulate import simulate, format_runlog

# metadata
metadata = {
    'protocolName': 'Simulate_PCR v1.0.0',
    'author': 'Dennis Simpson <dennis@email.unc.edu>',
    'description': 'Setup a ddPCR or Generic PCR'
    }

if __name__ == "__main__":
    protocol_file = open('PCR.py')
    labware_path = "{}{}custom_labware".format(os.getcwd(), os.sep)
    run_log, __bundle__ = simulate(protocol_file, custom_labware_paths=[labware_path])
    run_date = datetime.datetime.today().strftime("%a %b %d %H:%M %Y")
    i = 1
    t = format_runlog(run_log).split("\n")

    outstring = "Opentrons OT-2 Steps for {}.\nDate:  {}\nProgram File: PCR.py\n\nStep\tCommand\n" \
        .format(metadata['protocolName'], run_date)

    for l in t:
        outstring += "{}\t{}\n".format(i, l)
        i += 1
    if platform.system() == "Windows":
        outfile = open("C:{0}Users{0}{1}{0}Documents{0}ProgramFileSimulation.txt"
                       .format(os.sep, os.getlogin()), 'w', encoding="UTF-16")
        outfile.write(outstring)
        outfile.close()
    protocol_file.close()
