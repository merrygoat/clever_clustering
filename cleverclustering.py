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
    # retrieve 2 sub trees from cluster head
    clusterparticles = [int(linkagearray[maxclusterlocation, 0]),
                        int(linkagearray[maxclusterlocation, 1])]

    # Clusters are numbered starting at zero.
    # Clusters numbered 0 to (numparticles-1) are individual particles.
    # Clusters numbered > numparticles are compound clusters

    particle_counter = 0
    while particle_counter < len(clusterparticles):
        if clusterparticles[particle_counter] >= numparticles:
            clusterline = clusterparticles[particle_counter] - numparticles
            clusterparticles.append(int(linkagearray[clusterline, 0]))
            clusterparticles.append(int(linkagearray[clusterline, 1]))
        particle_counter += 1

    with open("clusteroutput.xyz", 'a') as xyzoutputfile:
        xyzoutputfile.write("{}\nComment line\n".format(str(numparticles)))

        for particle_id in range(numparticles):
            if particle_id in clusterparticles:
                xyzoutputfile.write("A ")  # write A for particles in cluster
            else:
                xyzoutputfile.write("B ")  # write B for particles not in cluster
            xyzoutputfile.write("{} {} {}\n".format(str(coordarray[particle_id][0]),
                                                    str(coordarray[particle_id][1]),
                                                    str(coordarray[particle_id][2])))


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
    plt.scatter(linkage_array[:, 2], linkage_array[:, 3], marker="o")
    plt.plot([cutoff, cutoff], [0, numparticles])
    plt.show()



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


def build_distance_array(coordarray, box_size):
    sqdist_wrapped = []
    for dimension in range(3):
        side_len = box_size[dimension]
        sqdist = pdist(np.swapaxes(np.atleast_2d(coordarray[:, dimension]), 0, 1), metric="sqeuclidean")
        dist_mask = sqdist > (side_len ** 2 / 2)
        sqdist[dist_mask] = side_len - sqdist[dist_mask]
        sqdist_wrapped.append(sqdist)

    distxyz = sqdist_wrapped[0] + sqdist_wrapped[1] + sqdist_wrapped[2]
    return distxyz


def read_particles_from_xyz(numparticles, xyzinput):
    coordarray = []
    for i in range(numparticles):
        line = xyzinput.readline()
        line = line[1:].lstrip().split()
        line = [float(x) for x in line]
        coordarray.append(line)
    return np.array(coordarray)


def clever_clustering(data_file, box_file, cutoff, printxyz=0):
    start = time()
    # open and close output file to delete old copies.
    open("clusteroutput.xyz", 'w').close()
    open("clustersize.txt", 'w').close()


    frame_number = 0
    box_size = read_box_size(box_file)

    with open(data_file, 'r') as xyzinput:
        line = xyzinput.readline()
        while line != "":  # read until EOF
            numparticles = int(line)  # read number of particles from first line
            xyzinput.readline()  # read to clear out the comment line

            # read in every particle coordinate in the current xyz frame into a list
            coordarray = read_particles_from_xyz(numparticles, xyzinput)

            if len(box_size) == 1:
                distxyz = build_distance_array(coordarray, box_size[0])
            else:
                distxyz = build_distance_array(coordarray, box_size[frame_number])

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

    if len(argv) == 4:
        print("Processing file " + argv[1] + "\n")
    else:
        print("Syntax: ./cleverclustering.py filetoanalyse.xyz boxsize.txt cutoff")
        exit()

    data_file = argv[1]
    box_file = argv[2]
    cut_off = float(argv[3])

    clever_clustering(data_file, box_file, cut_off)


if __name__ == "__main__":
    main()
