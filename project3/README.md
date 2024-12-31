# Architecture of the repository

In `lib`you will find the executable.
In `src`you will find the source files for the project.

# Sources

The main source is `dynalics.c`. It contains the main function and private functions to compute the trajectory of the atoms in a molecule with Verlet algorithm. A documentation for those functions is available at the beginning of the file.

# How to run

To compile and run, you just need to launch the `compile_and_run.sh` file. To use it, you need to call it with three arguments:
`compile_and_run.sh <file_path> <delta_time> <max_steps>` where `file_path` is the input file containing the initial state of the molecule, `delta_time` is the time step used for the Verlet Algorithm and `max_steps` is the maximum number of iterations of the Verlet Algorithm.

# Output

The output of the algorithm is a file `trajectory.txt` stored in the same directory as the `compile_and_run.sh` file.