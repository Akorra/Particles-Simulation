#include "simpar-mpi-v.h"

void usg_err(int rank){
    if(!rank){
      printf("\t[-] usage : mpirun -np <np> simpar <seed> <ncside> <n_par> <n_step>\n");
      printf("\t\t[-] int <np> : number of processes to be used.\n");
      printf("\t\t[-] int <seed> : seed for random number generation.\n");
      printf("\t\t[-] int <ncside> :  size of the grid (number of cells on the side.\n");
      printf("\t\t[-] int <n_par> :  number of particles\n");
      printf("\t\t[-] int <n_par> :  number of time-steps\n");
    }
    exit(1);
}

long long val_l(const char* arg, int rank){
    char *endptr;
    long long x = strtol(arg, &endptr, 10);
    if (endptr == arg) {
        if(!rank) printf("[-] ERROR: Invalid number: %s\n", arg);
        return 0;
    } else if (*endptr) {
        if(!rank) printf("[-] ERROR: Trailing characters after number: %s\n", arg);
        return 0;
    } else if (x <= 0) {
        if(!rank) printf("[-] ERROR: Number must be positive: %llu\n", x);
        return 0;
    }
    return x;
}

void init_particles(long seed, long ncside, long long n_part, long long loc_n, particle_t *par, particle_t* loc_par, cell_t** grid, double* masses){
    long long i;
    long cx, cy;

    srandom(seed);
    for(i=0; i < n_part; i++){
        par[i][X] = RND0_1;
        par[i][Y] = RND0_1;
        par[i][VX] = RND0_1 / ncside / 10.0;
        par[i][VY] = RND0_1 / ncside / 10.0;
        masses[i] = RND0_1 * ncside / (G * 1e6 * n_part);

        /* init_env*/
        cx = (long) par[i][X] * ncside;
        cy = (long) par[i][Y] * ncside;
        grid[cx][cy][M] += masses[i];
        grid[cx][cy][X] += masses[i] * par[i][X];
        grid[cx][cy][Y] += masses[i] * par[i][X];
    }
    MPI_Bcast(masses, n_part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid, ncside*ncside, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(par, loc_n, particle_mpi_t, loc_par, loc_n, particle_mpi_t, 0, MPI_COMM_WORLD);
}

cell_t** init_grid(const long ncside){
    cell_t** grid = (cell_t**)calloc(ncside, sizeof(cell_t*));
    for(long c=0; c<ncside; c++){
        grid[c] = (cell_t*)calloc(ncside, sizeof(cell_t));
        if(grid[c]==NULL) exit(0);
    }
    return grid;
}

void free_grid(cell_t** g, long ncside){
    for(long c=0; c<ncside; c++){
        free(g[c]);
    }
    free(g);
}

int main(int argc, char* argv[]){
  int n;                      /* Total number of particles  */
  int loc_n;                  /* Number of my particles     */
  int n_steps;                /* Number of timesteps        */
  int step;                   /* Current step               */
  int loc_part;               /* Current local particle     */
  double start_t, end_t;      /* Time Tags                  */
  double* par_m;               /* Particle Masses dont change */
  particle_t* par;            /* Positions of my particles  */
  particle_t* loc_par;            /* Positions of my particles  */
  cell_t** grid;              /* Velocities of my particles */

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(argc!=5){
    MPI_Finalize();
    if(!rank) printf("[-] ERROR: Invalid number of arguments... Expected 4 but got %d\n", argc-1);
    usg_err(rank);
  }

  const long seed = (long) val_l(argv[1], rank);
  const long ncside = (long) val_l(argv[2], rank);
  const long long n_par = val_l(argv[3], rank);
  const long n_step = (long) val_l(argv[4], rank);
  if(!(seed*ncside*n_par*n_step)){
    MPI_Finalize();
    usg_err(rank);
  }

  loc_n = n/comm_sz; /*should be evenly divisible by comm_sz*/
  loc_par = par + rank*loc_n;

  grid = init_grid(ncside);
  par = (particle_t*)calloc(n_par, sizeof(particle_t));
  loc_par = (particle_t*)calloc(loc_n, sizeof(particle_t));
  par_m = (double*)calloc(n_par, sizeof(double));

  MPI_Type_contiguous(4, MPI_DOUBLE, &particle_mpi_t);
  MPI_Type_commit(&particle_mpi_t);
  MPI_Type_contiguous(3, MPI_DOUBLE, &cell_mpi_t);
  MPI_Type_commit(&cell_mpi_t);

  init_particles(seed, ncside, n_par, loc_n, par, loc_par, grid, par_m);

  MPI_Gather(loc_par, loc_n, particle_mpi_t, par, loc_n, particle_mpi_t, 0, MPI_COMM_WORLD);
  if(rank){
    printf("Hey I'm %d and I got stuff %f\n", rank, loc_par[0][X]);
  }

  MPI_Finalize();
  return(0);
}
