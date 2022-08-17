#!/usr/bin/env python3
import os
import sys
import csv
import numpy
from os.path import exists


def readincsv(filename):
    #mat = np.loadtxt(open(filename, "rb"), delimiter=",", skiprows=1)
    reader = csv.reader(open(filename, "r"), delimiter=",")
    x = list(reader)
    mat = numpy.array(x).astype("float")
    myfunc = lambda x: 1 if (x > 0.4) else 0
    myfunc_vec = numpy.vectorize(myfunc)
    mat=myfunc_vec(mat)
    return mat
