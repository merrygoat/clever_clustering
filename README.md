# Clever Clustering
A python script to read in XYZ data and perform hierarchical aglomerative clustering

Required files:
* An XYZ file containing the points to be clustered
* A box file containing describing the size of the simulation box at each step in order to implement periodic boundaries

## XYZ specification
An XYZ file has the format:

Number of particles in frame 1 \
Comment line \
Particle_1_identifier x_coordinate y_coordinate z_coordinate\
Particle_2_identifier x_coordinate y_coordinate z_coordinate\
...\
Number of particles in frame 2\
Comment line\
Particle_1_identifier x_coordinate y_coordinate z_coordinate\
Particle_2_identifier x_coordinate y_coordinate z_coordinate\
...

For this code, the comment line and particle identifier are ignored.

## Box-file specification
The code reads in the size of the simulation box at each timestep from a box file. It uses this to calculate the periodic boundaries. The format of this file is designed to be the same as the TCC in order that the box file can be reused. The format is:

Comment line\
Frame_1 x_length y_length z_length\
Frame_2 x_length y_length z_length\

The comment line and frame numbers are ignored by teh code but are still required.

## Testing

Before running the code, the internal tests should be run to check the code is working correctly. Testing requires the pytest library. To test run:

"py.test test.py"
 
from the root directory