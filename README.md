# particleDynamics1D
## simulation code for simulating particle transport, interactions, fusions, on locally 1D domains
### Read doc/README.pdf for a detailed breakdown of keywords.

1. Compile files in the folder `source` using gfortran or similar compiler to generate an executable (partdynamics1D.exe). The default compilation requires navigating to the `source` folder in a terminal, and typing `make`
2. Navigate to the directory in which the executable is located, and type `./partdynamics.exe example` to run an example simulation of particle dynamics. Parameters for the example simulation are located in `paramfiles/param.example`
3. The executable generates a snapshot file which is stored in `results/snapfiles/example.snap.out`. The MATLAB script `plotSimDensities.m` reads the snapshot file, and saves the data to `results/example_processed.mat`. The MATLAB script also plots the distribution of different particle types in the system.