#!/usr/bin/env python

# cleverclustering.py - a python script to read in XYZ data and perform hierarchical aglomerative clustering
# Peter Crowther - November 2014

from time import time
from math import sqrt
from sys import argv, stdout, exit
from scipy.cluster.hierarchy import linkage
# import matplotlib.pyplot as plt


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


def clever_clustering(data_file, box_file):
    start = time()
    # open and close output file to delete old copies.
    xyzoutputfile = open("clusteroutput.xyz", 'w')
    xyzoutputfile.close()

    xyzinput = open(data_file, 'r')
    boxsizeinput = open(box_file, 'r')
    boxsizeinput.readline()  # Read in comment line to flush it
    framecounter = 0
    line = xyzinput.readline()

    while line != "":  # read until EOF
        coordarray = []
        distancearray = []
        numparticles = int(line)  # read number of particles from first line
        xyzinput.readline()  # read to clear out the comment line

        # process box coordinates to determine boundary conditions
        boxcoords = boxsizeinput.readline()  # read box coordinates from second line
        boxcoords = boxcoords.split()
        xlen = float(boxcoords[1])
        ylen = float(boxcoords[2])
        zlen = float(boxcoords[3])

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

        cutoff = 2.2

        # plt.ion()
        # plt.figure()
        # plt.scatter(linkagearray[:,2], linkagearray[:,3], marker="o")
        # plt.plot([cutoff, cutoff], [0,numparticles])
        # plt.show()
        # plt.draw()

        # np.savetxt("yyy.link", linkagearray, fmt="%1i, %1i, %10.5f, %1i")

        maxclustersize = [0, 0]  # max cluster size and location of max cluster

        for i in range(0, numparticles):
            if linkagearray[i, 3] > maxclustersize[0]:
                maxclustersize[0] = linkagearray[i, 3]
                maxclustersize[1] = i
            if float(linkagearray[i, 2]) > cutoff:
                break

        numberoutputfile = open("clustersize.txt", 'a')
        numberoutputfile.write(str(maxclustersize[0]))
        numberoutputfile.write("\n")
        numberoutputfile.close()
        printxyzoutput(maxclustersize[1], linkagearray, coordarray)

        framecounter += 1
        if framecounter % 10 == 0:
            stdout.write('\rNumber of frames processed: ' + str(framecounter))
            stdout.flush()

        line = xyzinput.readline()

    xyzinput.close()
    boxsizeinput.close()

    print("\nTime taken: " + str(time() - start))


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
