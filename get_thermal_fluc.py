import sys
import numpy as np
from Dump import Dump
from DumpSeries import DumpSeries

def get_thermal_fluc(filename_0, filename_end,dump_step):
    ds = DumpSeries(filename_0,filename_end,dump_step)
    print("Tracking {} planes".format(len(ds.planes)))
    print("RMS for Z coord: {}".format(round(ds.RMS,4)))
if len(sys.argv) != 3:
    print('Incorrect arguments: [0 filename, end filename]')
else:
    filename_0 = sys.argv[1]
    filename_end = sys.argv[2]
    dump_step = 100000 # hard coded for now
    get_thermal_fluc(filename_0,filename_end,dump_step)