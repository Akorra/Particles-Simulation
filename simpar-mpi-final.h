/**
 * @file simpar-mpi.h
 * @authors: Filipe Marques
 * @date 25 Apr 2019
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
#define G 6.67408e-11                                   /*Gravitation*/
#define EPSLON 0.0005                                   /*distance threshold*/

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
    double m;
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

int rank, comm_sz;
double t_mass = 0.0;            /* global variable to hold the total mass of the grid */
double t_cx = 0.0, t_cy = 0.0;  /* global variables to hold  the position of the total center of mass*/

double* masses;
particle_t* par;
cell_t* grid;

MPI_Datatype MPI_PARTICLE_T;

/**
 * @brief Output usage command error to console.
 **/
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

#endif
