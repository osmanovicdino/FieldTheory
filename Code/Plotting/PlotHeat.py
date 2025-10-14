#!/usr/bin/env python3
import os
import sys
import csv
import numpy
import matplotlib.pyplot as plt
import matplotlib
from os.path import exists
matplotlib.use('Agg')
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import




# plt.imshow(np.random.random((50, 50)))
# plt.colorbar()
# plt.show()


# def graph(filename):
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     outputfilename = os.path.splitext(filename)[0]+'.jpg'
#     with open(filename, newline='') as csvfile:
#         spamreader = csv.reader(csvfile, delimiter=',',
#                                 quoting=csv.QUOTE_NONNUMERIC)
#         for eachLine in spamreader:
#             xpos = eachLine[0]
#             ypos = eachLine[1]
#             zpos = eachLine[2]

#             ax.scatter(xpos, ypos, zpos, color='blue')

#             ax.set_xlabel('x')
#             ax.set_ylabel('y')
#             ax.set_zlabel('z')

#     plt.savefig(outputfilename, format='jpg')
#     plt.close(fig)


def graph(filename):
    #mat = np.loadtxt(open(filename, "rb"), delimiter=",", skiprows=1)
    outputfilename = os.path.splitext(filename)[0]+'.png'
    file_exists = exists(outputfilename)
    if file_exists == 0:
        reader = csv.reader(open(filename, "r"), delimiter=",")
        x = list(reader)
        mat = numpy.array(x).astype("float")
        fig = plt.figure()
        plt.imshow(mat)
        #plt.clim(0., 10.)
        #plt.colorbar()
        plt.savefig(outputfilename, format='png')
        plt.close(fig)


def graph22(filename1,filename2,filename3,filename4) :
    outputfilename = os.path.splitext(filename1)[0]+'.png'
    file_exists = exists(outputfilename)
    if file_exists == 0:    
        reader1 = csv.reader(open(filename1, "r"), delimiter=",")
        x1 = list(reader1)
        mat1 = numpy.array(x1).astype("float")
        reader2 = csv.reader(open(filename2, "r"), delimiter=",")
        x2 = list(reader2)
        mat2 = numpy.array(x2).astype("float")
        reader3 = csv.reader(open(filename3, "r"), delimiter=",")
        x3 = list(reader3)
        mat3 = numpy.array(x3).astype("float")
        reader4 = csv.reader(open(filename4, "r"), delimiter=",")
        x4 = list(reader4)
        mat4 = numpy.array(x4).astype("float")
        fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
        axes = axes.ravel()
        mats = [mat1, mat2, mat3, mat4]
        files = [filename1, filename2, filename3, filename4]

        # Optional: common color scale across all four
        vmin = vmax = None


        fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
        axes = axes.ravel()

        ims = []
        for ax, mat, fname in zip(axes, mats, files):
            # Use shared vmin/vmax if requested; otherwise let matplotlib autoscale
            im = ax.imshow(mat)
            ax.set_title(os.path.basename(fname))
            ax.set_xticks([])
            ax.set_yticks([])
            ims.append(im)

        # One shared colorbar if using shared scale; otherwise skip to stay minimal
        plt.savefig(outputfilename, format="png", dpi=200)
        plt.close(fig)

