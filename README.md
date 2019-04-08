<div align="center">
  <a><img src="https://fenix.tecnico.ulisboa.pt/api/bennu-portal/configuration/logo"></a><br><br>
</div>

# This is an effin particle Simulator

## Files:
- `simpar.h`: Header file for the serial solution, includes function headers, data structures and other library headers.
- `simpar.c`: Source file for the serial solution wich implements main and functions declared in `simpar.h`.
- `simpar-omp.h`: Header file for the parallel solution, includes function headers, data structures and other library headers.
- `simpar-omp.c`: Source file for the sparallel solution wich implements main and functions declared in `simpar-omp.h`.


## Objective

The Goal of this project is to develop a Particle simulator which calculates at each instant the updated position of each particle in space, due to the effect of gravitational pull from other particles.
A first instance is a serial implementation of the simulator, then we will produce a parallelized version of the solution using openMP, and finally the last instance will be implemented using MPI.

## Command Line Parameters

| Parameter     | data type  | Description                                         |
|:--------------|:----------:|:----------------------------------------------------|
| seed          | long       | Seed for random number generator.                   |
| ncside        | long       | Size of the 2D grid - number of cells on the side.  |
| n_part        | long long  | Number of particles to simulate.                    |
| n_steps       | long       | Number of time-steps to simulate.                   |


## Compiling and Running

To run this project solutions GCC is advised.

### Compile
- Serial solution:    `gcc simpar.c -lm -o simpar`
- Parallel solution:  `gcc simpar-omp.c -lm -fpopenmp -o simpar-omp`

### Run
- serial solution:    `./simpar <seed> <ncside> <n_part> <n_steps>`
- parallel solution:  `./simpar-omp <seed> <ncside> <n_part> <n_steps>`
