import cleverclustering as cc
from os import remove
import numpy as np


def test_run_clustering():
    try:
        cc.clever_clustering("test/test.xyz", "test/test_box.txt")
        return 0
    except Exception as e:
        print("Failed test_run_clustering")
        return 1


def test_read_box():
    try:
        cc.read_box_size("test/test_box.txt")
        return 0
    except Exception as e:
        print("Failed test_red_box")
        return 1


def test_cluster_output():
    measured_clusters_file = "clusteroutput.xyz"
    known_clusters_file = "test/sample_clusteroutput.xyz"
    with open(measured_clusters_file) as tmp:
        measured_particles = int(tmp.readline())
    measured_clusters = np.genfromtxt(measured_clusters_file, skip_header=2, max_rows=measured_particles, dtype=None)

    with open(known_clusters_file) as tmp:
        known_particles = int(tmp.readline())
    known_clusters = np.genfromtxt(known_clusters_file, skip_header=2, max_rows=known_particles, dtype=None)

    if np.all(measured_clusters == known_clusters):
        return_val = 0
    else:
        return_val = 1
    remove(measured_clusters_file)
    return return_val


def test_cluster_size():
    measured_cluster_size_file = "clustersize.txt"
    measured_cluster_size = np.loadtxt(measured_cluster_size_file)
    known_cluster_size = np.loadtxt("test/sample_clustersize.txt")
    if np.all(measured_cluster_size == known_cluster_size):
        return_val = 0
    else:
        return_val = 1
    remove(measured_cluster_size_file)
    return return_val


def test_clustering():
    assert test_read_box() == 0
    assert test_run_clustering() == 0
    assert test_cluster_output() == 0
    assert test_cluster_size() == 0


test_clustering()
