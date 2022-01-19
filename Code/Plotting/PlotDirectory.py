import csv
import sys
# from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
# import os
from PlotHeat import graph


cmdargs = list(map(str,sys.argv))

filename =  cmdargs[1]


graph(filename)