#include "simpar-mpi.h"

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

void init_particles(long seed, long ncside, long long n_part, double* m, vect_dbl* pos, vect_dbl* loc_vel, int loc_n){
    long long i;
    srandom(seed);
    for(i=0; i < n_part; i++){
        pos[i][0] = RND0_1;
        pos[i][1] = RND0_1;
        vel[i][0] = RND0_1 / ncside / 10.0;
        vel[i][1] = RND0_1 / ncside / 10.0;
        m[i] = RND0_1 * ncside / (G * 1e6 * n_part);
    }
    MPI_Bcast(m, n_part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(pos, n_part, vect_mpi_dbl, 0, MPI_COMM_WORLD);
    MPI_Scatter(vel, loc_n, vect_mpi_dbl, loc_vel, loc_n, vect_mpi_dbl, 0, MPI_COMM_WORLD);
}

int main(int argc, char *argv[]){
  long long loc_n; /* Number of particles local to process */
  double* m;             /* All the masses             */
  vect_dbl* loc_pos;            /* Positions of my particles  */
  vect_lng* loc_idx;            /* Positions of my particles  */
  vect_dbl* pos;                /* Positions of all particles */
  vect_lng* idx;
  vect_dbl* vel;
  vect_dbl* loc_vel;            /* Velocities of my particles */
  vect_dbl* loc_a;         /* Forces on my particles     */

  double start, finish;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  if(argc!=5){
    MPI_Finalize();
    if(!id) printf("[-] ERROR: Invalid number of arguments... Expected 4 but got %d\n", argc-1);
    usg_err(id);
  }

  const long seed = (long) val_l(argv[1], id);
  const long ncside = (long) val_l(argv[2], id);
  const long long n_par = val_l(argv[3], id);
  const long n_step = (long) val_l(argv[4], id);
  if(!(seed*ncside*n_par*n_step)){
    MPI_Finalize();
    usg_err(id);
  }

  loc_n = n_par/p;
  m = malloc(n_par * sizeof(double));
  pos = malloc(n_par * sizeof(vect_dbl));
  loc_a = malloc(loc_n * sizeof(vect_dbl));
  loc_pos = pos + id*loc_n;
  loc_vel = malloc(loc_n * sizeof(vect_dbl));
  if(!id) vel = malloc(n_par * sizeof(vect_dbl));

  MPI_Type_contiguous(2, MPI_DOUBLE, &vect_mpi_dbl);
  MPI_Type_commit(&vect_mpi_dbl);
  MPI_Type_contiguous(2, MPI_LONG, &vect_mpi_lng);
  MPI_Type_commit(&vect_mpi_lng);

  init_particles(seed, ncside, n_par, m, pos, loc_vel, loc_n);

  start = MPI_Wtime();
  /*YOUR CODE HERE*/

  finish = MPI_Wtime();

  MPI_Finalize();
  exit(0);
}
