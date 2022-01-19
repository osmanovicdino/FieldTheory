#!/usr/bin/env python3
import os
import sys
import csv
import numpy
import matplotlib.pyplot as plt
import matplotlib
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
    reader = csv.reader(open(filename, "r"), delimiter=",")
    x = list(reader)
    mat = numpy.array(x).astype("float")
    fig = plt.figure()
    outputfilename = os.path.splitext(filename)[0]+'.png'
    plt.imshow(mat)
    plt.colorbar()
    plt.savefig(outputfilename, format='png')
    plt.close(fig)



