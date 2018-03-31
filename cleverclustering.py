#!/usr/bin/env python

# cleverclustering.py - a python script to read in XYZ data and perform hierarchical aglomerative clustering

from time import time
from sys import argv, stdout, exit
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
import numpy as np


def printxyzoutput(maxclusterlocation, linkagearray, coordarray):
    """ Write an xyz file where A particles are in the in the largest cluster and
     B particles are not in the largest cluster """

    numparticles = len(coordarray)
    num_a_particles = 0
    num_b_particles = 0
    clusterparticles = []
    i = 0
    xyzoutputfile = open("clusteroutput.xyz", 'a')

    xyzoutputfile.write(str(numparticles) + "\n")  # Write numparticles to first line of output
    xyzoutputfile.write("Comment line\n")

    # retrieve 2 sub trees from cluster head
    clusterparticles.append(int(linkagearray[maxclusterlocation, 0]))
    clusterparticles.append(int(linkagearray[maxclusterlocation, 1]))

    # clusters start at zero. Clusters 0 to (numparticles-1) are individual particles. Clusters numparticles +
    # are compound clusters

    while i < len(clusterparticles):
        if clusterparticles[i] >= numparticles:
            clusterline = clusterparticles[i] - numparticles
            clusterparticles.append(int(linkagearray[clusterline, 0]))
            clusterparticles.append(int(linkagearray[clusterline, 1]))
        i += 1

    for i in range(0, numparticles):
        if i in clusterparticles:
            xyzoutputfile.write("A ")  # write A for particles in cluster
            num_a_particles += 1
        else:
            xyzoutputfile.write("B ")  # write B for particles not in cluster
            num_b_particles += 1
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
                except ValueError:
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


def get_max_cluster_size(linkagearray, cutoff):
    # Take a linkage array, and find the largest cluster under a certain distance cutoff.
    cutoff_location = np.searchsorted(linkagearray[:, 2], cutoff)
    largest_location = np.argmax(linkagearray[:cutoff_location, 3])
    largest_cluster = linkagearray[largest_location, 3]

    return [largest_cluster, largest_location]


def write_output_files(coordarray, linkagearray, maxclustersize, printxyz):
    # Write to output files
    with open("clustersize.txt", 'a') as numberoutputfile:
        numberoutputfile.write(str(maxclustersize[0]) + "\n")
    if printxyz == 1:
        printxyzoutput(maxclustersize[1], linkagearray, coordarray)


def build_distance_array(coordarray, xlen, ylen, zlen):
    # build a reduced distance array, taking into account PBCs
    distx = pdist(np.swapaxes(np.atleast_2d(coordarray[:, 0]), 0, 1))
    disty = pdist(np.swapaxes(np.atleast_2d(coordarray[:, 1]), 0, 1))
    distz = pdist(np.swapaxes(np.atleast_2d(coordarray[:, 2]), 0, 1))
    distx[distx > (xlen / 2)] = xlen - distx[distx > (xlen / 2)]
    disty[disty > (ylen / 2)] = ylen - disty[disty > (ylen / 2)]
    distz[distz > (zlen / 2)] = zlen - distz[distz > (zlen / 2)]
    distxyz = distx ** 2 + disty ** 2 + distz ** 2
    return distxyz


def read_particles_from_xyz(numparticles, xyzinput):
    coordarray = []
    for i in range(numparticles):
        line = xyzinput.readline()
        coordarray.append(line[1:].lstrip().split("\t"))
    return np.array(coordarray)


def clever_clustering(data_file, box_file, cutoff=2.2, printxyz=0):
    start = time()
    # open and close output file to delete old copies.
    xyzoutputfile = open("clusteroutput.xyz", 'w')
    xyzoutputfile.close()
    sizefile = open("clusteroutput.xyz", 'w')
    sizefile.close()

    frame_number = 0
    box_size = read_box_size(box_file)

    with open(data_file, 'r') as xyzinput:
        line = xyzinput.readline()
        while line != "":  # read until EOF
            numparticles = int(line)  # read number of particles from first line
            xyzinput.readline()  # read to clear out the comment line

            xlen = box_size[frame_number][0]
            ylen = box_size[frame_number][1]
            zlen = box_size[frame_number][2]

            # read in every particle coordinate in the current xyz frame into a list
            coordarray = read_particles_from_xyz(numparticles, xyzinput)

            distxyz = build_distance_array(coordarray, xlen, ylen, zlen)

            linkagearray = linkage(distxyz)  # linkage outputs a numpy array

            # plot_linkage_array(cutoff, linkagearray, numparticles)

            maxclustersize = get_max_cluster_size(linkagearray, cutoff*cutoff)

            write_output_files(coordarray, linkagearray, maxclustersize, printxyz)

            # Update user on progress
            frame_number += 1
            if frame_number % 10 == 0:
                stdout.write('\rNumber of frames processed: ' + str(frame_number))
                stdout.flush()

            line = xyzinput.readline()

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
