#!/usr/bin/env python3


import csv
import sys
# from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
# import os
from PlotHeat import graph
from PlotHeat import graph22


cmdargs = list(map(str,sys.argv))

filename1 =  cmdargs[1]
filename2 =  cmdargs[2]
filename3 =  cmdargs[3]
filename4 =  cmdargs[4]

graph22(filename1,filename2,filename3,filename4)
