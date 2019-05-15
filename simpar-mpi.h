/**
 * @file simpar.h
 * @authors: Filipe Marques, Luís Fonseca
 * @date 16 Mar 2019
 * @brief Header File containing the particle simulation data structures and functions headers.
 *
 * Both the structs and relevant functions are here defined for usage in
 * the n-body gravitational simulation. Some physical constants are defined
 * as variables. For a detailed overview of the project visit the link below.
 */

#ifndef simpar_mpi_h
#define simpar_mpi_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define RND0_1 ((double) random() / ((long long)1<<31)) /*Random number generator*/
#define G 6.67408e-11 /*Gravitation*/
#define EPSLON 0.0005 /*distance threshold*/

typedef double vect_dbl[2];
typedef long   vect_lng[2];

int id, p;
MPI_Datatype vect_mpi_dbl;
MPI_Datatype vect_mpi_lng;

/**
 * @brief Particle.
 *
 * A particle containing information about its position,
 * velocity, mass and grid position at a given time-step
 * and its constant mass.
 **/
typedef struct particle_t{
    double x, y;    /**<position (x, y) in 2D space. */
    double vx, vy;  /**<velocity (vx, vy) in 2D space. */
    double m;       /**<mass. */
    long cx, cy;    /**<grid index position (cx, cy). */
}particle_t;

/**
 * @brief Grid cell in 2D space.
 *
 * A grid cell containing a point representing the center of mass and the total
 * mass of all the particles inside this cell..
 **/
typedef struct cell_t{
    double x, y; /**< center of mass of the cell. */
    double M;    /**< total mass of the cell. */
}cell_t;

double t_mass = 0.0;            /* global variable to hold the total mass of the grid */
double t_cx = 0.0, t_cy = 0.0;  /* global variables to hold  the position of the total center of mass*/

/**
 * @brief Output usage command error to console.
 *
 */
void usg_err(int rank);

/**
 * @brief Validate a char to parse the corresponding positive integer.
 *
 * Parses the value pointed by a character pointer to extarct a positive long
 * numerical value if this can't be accomplished or the value is negative then
 * the function returns 0 with an error message.
 *
 * @param arg char ptr from which's pointed value to parse a numeric.
 * @return long integer parsed from arg on success, or 0 on fail.
 */
long long val_l(const char* arg, int rank);

/**
 * @brief Initialize grid  by allocating memory for each cell.
 *
 * Initialization of a square 2D array of type cell_t, using a
 * given ncside to allow for cell_t memory allocation. Each cell
 * center of mass positions and mass are initialized to 0.
 *
 * @param ncside Size of the grid (number of cells on the side).
 * @return pointer to cell_t 2D array (cell_t**)
 */
cell_t** init_grid(const long ncside);

/**
 * @brief Frees up the space previously allocated for the grid;
 *
 * Iterates through input grid, freeing each cell line by line,
 *  and then finally frees space allocated for the grid.
 *
 * @param g double pointer to grid to be freed.
 * @param ncside Size of the grid (number of cells on the side).
 */
void free_grid(cell_t** g, long ncside);

/**
 * @brief Initialize particle system with random physical properties.
 *
 * Initialization of an array of type particle_t, using a
 * given random seed to allow for result replication. Each particle
 * position, velocity and mass are initialized randomly.
 *
 * @param seed Seed for the random number generator.
 * @param ncside Size of the grid (number of cells on the side).
 * @param n_part Number of particles.
 * @param p Pointer to an array of particles.
 */
void init_particles(long seed, long ncside, long long n_part, particle_t *par);

/**
 * @brief Initialize Environment, i.e, grid and set particles grid positions.
 *
 * Initialization of each particle's initial grid indices, and from there
 *
 *
 * @param g Double pointer to grid to be freed.
 * @param ncside Size of the grid (number of cells on the side).
 * @param p Pointer to an array of particles.
 * @param n_part Number of particles.
 */
void init_env(cell_t** g, long ncside, particle_t* p, long long n_par);

/**
 * @brief Compute the gravitational accelleration applied to a particle, by a given cell's center of mass.
 *
 * The accelleration relative to a given center of mass applied on a particle
 * is calculated only from knoing that F=m_p*a and F=G*(M*m_a)/d^2 which holds
 * a=G*M/d^2, if the distance (d) is bellow EPSLON then it is counted as a zero
 * accelleration contributionnction.
 *
 * @param ax Pointer to a variable (initialized to 0.0) to which to add the x component of acceleration due to c's center of mass.
 * @param ax Pointer to a variable (initialized to 0.0) to which to add the y component of acceleration due to c's center of mass.
 * @param c Pointer to struct of type cell_t.
 * @param m Mass of the particle being considered
 * @param x Position, in axis x, of the particle.
 * @param y Position, in axis y, of the particle.
 */
void accellerate_p(double* ax, double* ay, const cell_t* c, double m, double x, double y);

/**
 * @brief Calculate the new velocity and then the new position of each particle.
 *
 * Update each particle in the system state at a given time-step of the simulation.
 * The function time complexity is O(n) where n is the number of particles.
 * Update the centers of mass in the grid.
 * On the last iteration the total center of mass is also calculated.
 *
 * @param x Pointer to a 2D array of cells.
 * @param ncside Size of the grid (number of cells on the side).
 * @param p Pointer to an array of particles.
 * @param n_part Number of particles.
 * @param n_step Number of timesteps.
 * @param step Current timestep.
 */
void update_particles(cell_t** x, long ncside, particle_t* par, long long n_par, long n_step, long step);

#endif /* simpar_mpi_h */
