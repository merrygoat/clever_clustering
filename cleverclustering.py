#!/usr/bin/env python

# cleverclustering.py - a python script to read in XYZ data and perform hierarchical aglomerative clustering
# Peter Crowther - November 2014

from time import time
from math import sqrt
from sys import argv, stdout, exit
from scipy.cluster.hierarchy import linkage
import numpy as np


def printxyzoutput(maxclusterlocation, linkagearray, coordarray):
    numparticles = len(coordarray)
    numAparticles = 0
    numBparticles = 0
    clusterparticles = []
    i = 0
    xyzoutputfile = open("clusteroutput.xyz", 'a')

    xyzoutputfile.write(str(numparticles) + "\n")  # Write numparticles to first line of output
    xyzoutputfile.write("Comment line\n")

    # retrieve 2 sub trees from cluster head
    clusterparticles.append(int(linkagearray[maxclusterlocation, 0]))
    clusterparticles.append(int(linkagearray[maxclusterlocation, 1]))

    # clusters start at zero. Clusters 0 to (numparticles-1) are individual particles. Clusters numparticles +
    # are compund clusters

    while i < len(clusterparticles):
        if clusterparticles[i] >= numparticles:
            clusterline = clusterparticles[i] - numparticles
            clusterparticles.append(int(linkagearray[clusterline, 0]))
            clusterparticles.append(int(linkagearray[clusterline, 1]))
        i += 1

    for i in range(0, numparticles):
        if i in clusterparticles:
            xyzoutputfile.write("A ")  # write A for particles in cluster
            numAparticles += 1
        else:
            xyzoutputfile.write("B ")  # write B for particles not in cluster
            numBparticles += 1
        xyzoutputfile.write(str(coordarray[i][0]) + " " + str(coordarray[i][1]) + " " + str(coordarray[i][2]))
    xyzoutputfile.close()


def read_box_size(box_file):
    box_size = []

    with open(box_file, 'r') as box_size_file:
        box_size_file.readline()  # Read in comment line to flush it
        for line_num, line in enumerate(box_size_file):
            split_line = line.split()
            if len(split_line) != 4:
                print("Error in box file on line %d.\n %s", line_num, line)
                exit()
            else:
                try:
                    split_line[1] = float(split_line[1])
                    split_line[2] = float(split_line[2])
                    split_line[3] = float(split_line[3])
                except ValueError as e:
                    print("Error on reading box file, line %d, cannot convert value to float.", line_num)
                    return 1
                box_size.append(split_line[1:])
    return box_size


def plot_linkage_array(cutoff, linkage_array, numparticles):
    import matplotlib.pyplot as plt
    plt.ion()
    plt.figure()
    plt.scatter(linkage_array[:, 2], linkage_array[:, 3], marker="o")
    plt.plot([cutoff, cutoff], [0, numparticles])
    plt.show()
    plt.draw()


def get_max_cluster_size(linkagearray, numparticles, cutoff):

    maxclustersize = [0, 0]  # max cluster size and location of max cluster

    for i in range(0, numparticles):
        if linkagearray[i, 3] > maxclustersize[0]:
            maxclustersize[0] = linkagearray[i, 3]
            maxclustersize[1] = i
        if float(linkagearray[i, 2]) > cutoff:
            break

    return maxclustersize


def clever_clustering(data_file, box_file, cutoff=2.2):
    start = time()
    # open and close output file to delete old copies.
    xyzoutputfile = open("clusteroutput.xyz", 'w')
    xyzoutputfile.close()

    xyzinput = open(data_file, 'r')

    frame_number = 0
    line = xyzinput.readline()

    box_size = read_box_size(box_file)

    while line != "":  # read until EOF
        coordarray = []
        distancearray = []
        numparticles = int(line)  # read number of particles from first line
        xyzinput.readline()  # read to clear out the comment line

        xlen = box_size[frame_number][0]
        ylen = box_size[frame_number][1]
        zlen = box_size[frame_number][2]

        # read in every particle coordinate in the current xyz frame into a list
        for i in range(1, numparticles + 1):
            line = xyzinput.readline()
            coordarray.append(line[1:].lstrip().split("\t"))

        # build a reduced distance array, taking into account PBCs
        for i in range(0, numparticles):
            for j in range(i + 1, numparticles):
                distx = abs(float(coordarray[i][0]) - float(coordarray[j][0]))
                disty = abs(float(coordarray[i][1]) - float(coordarray[j][1]))
                distz = abs(float(coordarray[i][2]) - float(coordarray[j][2]))
                if distx > xlen / 2:
                    distx = xlen - distx
                if disty > ylen / 2:
                    disty = ylen - disty
                if distz > zlen / 2:
                    distz = zlen - distz
                distxyz = sqrt(distx ** 2 + disty ** 2 + distz ** 2)
                distancearray.append(distxyz)

        linkagearray = linkage(distancearray)  # linkage outputs a numpy array

        #plot_linkage_array(cutoff, linkagearray, numparticles)

        maxclustersize = get_max_cluster_size(linkagearray, numparticles, cutoff)

        # Write to output files
        with open("clustersize.txt", 'a') as numberoutputfile:
            numberoutputfile.write(str(maxclustersize[0]) + "\n")
        printxyzoutput(maxclustersize[1], linkagearray, coordarray)

        # Update user on progress
        frame_number += 1
        if frame_number % 10 == 0:
            stdout.write('\rNumber of frames processed: ' + str(frame_number))
            stdout.flush()

        line = xyzinput.readline()

    xyzinput.close()

    print("\nTime taken: {:.2f} seconds".format(time() - start))


def main():

    if len(argv) == 3:
        print("Processing file " + argv[1] + "\n")
    else:
        print("Syntax: ./cleverclustering.py filetoanalyse.xyz boxsize.txt")
        exit()

    data_file = argv[1]
    box_file = argv[2]

    clever_clustering(data_file, box_file)


if __name__ == "__main__":
    main()
