#include "simpar-mpi-final.h"

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

void init_particles(long seed, long ncside, long long n_part, particle_t *par, cell_t* grid){
    long long i;
    long idx;

    srandom(seed);
    for(i=0; i < n_part; i++){
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;
        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);

        /* init_env*/
        idx = ((par[i].x*ncside)+par[i].y)*ncside;
        grid[idx].M += par[i].m;
        grid[idx].x += par[i].m * par[i].x;
        grid[idx].y += par[i].m * par[i].y;
    }
}

void scatter_particles(long long n_par, long long loc_n, particle_t *loc_par){
  int i, sum=0;
  int counts[comm_sz], disps[comm_sz];
  double rem = n_par%comm_sz;

  for(i = 0; i<comm_sz; i++){
  	counts[i] = n_par/comm_sz;
  	if(rem>0){
  		counts[i]++;
  		rem--;
  	}
  	disps[i] = sum;
  	sum += counts[i];
  }
  MPI_Scatterv(par, counts, disps, MPI_PARTICLE_T, loc_par, loc_n, MPI_PARTICLE_T, 0, MPI_COMM_WORLD);
}

void init_grid(long ncside){
  for(long idx=0; idx<ncside*ncside; idx++){
    MPI_Allreduce(&(grid[idx].M), &(grid[idx].M), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&(grid[idx].x), &(grid[idx].x), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&(grid[idx].y), &(grid[idx].y), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
}

void exit_routine(particle_t* loc){
  if(!rank && par){
    free(par);
  }
  if(grid) free(grid);
  if(loc) free(loc);

  MPI_Finalize();
  exit(1);
}

void accellerate_p(double* ax, double* ay, const cell_t* c, double m, double x, double y){
    if((c->M) == 0.0) return ;

    double mag;
    double dx = ((c->x)/(c->M)) - x;
    double dy = ((c->y)/(c->M)) - y;

    double d_2 = (dx*dx)+(dy*dy);
    if(sqrt(d_2) < EPSLON){
        return;
    }

    mag = (((c->M)*G)/d_2);
    *ax += dx * mag;
    *ay += dy * mag;
}

void update_particles(long ncside, particle_t* par, long long n_par, long n_step){
    double m, px, py, ax, ay;
    long cx, cy, nx, ny, ux, uy, lx, ly, o_idx, n_idx;
    cell_t* x = (cell_t*)calloc(ncside*ncside, sizeof(cell_t));

    for(long step = 0; step < n_step; step++){
      for(long long i=0; i<n_par; i++){
        m = par[i].m;
        px = par[i].x;
        py = par[i].y;
        cx = (long) px * ncside, nx;
        cy = (long) py * ncside, ny;
        ux = cx+1; uy = cy+1; lx = cx-1; ly = cy-1;
        if(ux >= ncside) ux = 0;
        else if(lx < 0) lx = ncside-1;
        if(uy >= ncside) uy = 0;
        else if(ly < 0) ly = ncside-1;

        ax = 0.0;
        ay = 0.0;

        accellerate_p(&ax, &ay, &(grid[cx*ncside+cy]), m, px, py); // current cell
        accellerate_p(&ax, &ay, &(grid[ux*ncside+cy]), m, px, py); // right cell
        accellerate_p(&ax, &ay, &(grid[lx*ncside+cy]), m, px, py); // left cell
        //upper adjacents
        accellerate_p(&ax, &ay, &(grid[cx*ncside+uy]), m, px, py); // upper cell
        accellerate_p(&ax, &ay, &(grid[lx*ncside+uy]), m, px, py); // upper left cell
        accellerate_p(&ax, &ay, &(grid[ux*ncside+uy]), m, px, py); // upper right cell
        //lower adjacents
        accellerate_p(&ax, &ay, &(grid[cx*ncside+ly]), m, px, py); // lower cell
        accellerate_p(&ax, &ay, &(grid[lx*ncside+ly]), m, px, py); // lower left cell
        accellerate_p(&ax, &ay, &(grid[ux*ncside+ly]), m, px, py); // lower right cell

        //update velocity
        par[i].vx += ax;
        par[i].vy += ay;

        //update position
        par[i].x += par[i].vx + ax*0.5;
        while(par[i].x >= 1.0) par[i].x -= 1.0;
        while(par[i].x < 0.0) par[i].x += 1.0;

        par[i].y += par[i].vy + ay*0.5;
        while(par[i].y >= 1.0) par[i].y -= 1.0;
        while(par[i].y < 0.0) par[i].y += 1.0;

        //update cells if cell changed maybe outside loop?
        nx = (long) par[i].x*ncside;
        ny = (long) par[i].y*ncside;
        if(cx-nx || cy-ny){
          o_idx = cx*ncside+cy;
          n_idx = nx*ncside+ny;

          x[o_idx].M -= m;
          x[o_idx].x -= m * px;
          x[o_idx].y -= m * py;

          x[n_idx].M += m;
          x[n_idx].x += m * par[i].x;
          x[n_idx].y += m * par[i].y;
        }

        if(n_step-1-step == 0){
          t_mass += m;
          t_cx += m * par[i].x;
          t_cy += m * par[i].y;
        }
      }

      for(long c = 0; c<ncside*ncside; c++){
    		 MPI_Allreduce(&(x[c].M), &(grid[c].M), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(&(x[c].x), &(grid[c].x), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(&(x[c].y), &(grid[c].y), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         x[c].M = 0.0;
         x[c].x = 0.0;
         x[c].y = 0.0;
       }
       MPI_Barrier(MPI_COMM_WORLD);
     }

     MPI_Allreduce(&(t_mass), &(t_mass), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     MPI_Allreduce(&(t_cx), &(t_cx), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     MPI_Allreduce(&(t_cy), &(t_cy), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

int main(int argc, char *argv[]){
  double strt_t, end_t;
  long long loc_n;
  particle_t* loc_par;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //declare particle datatype
  MPI_Type_contiguous(5, MPI_DOUBLE, &MPI_PARTICLE_T);
  MPI_Type_commit(&MPI_PARTICLE_T);

  //check parameters
  if(argc!=5){
    if(!rank) printf("[-] ERROR: Invalid number of arguments... Expected 4 but got %d\n", argc-1);
    usg_err(rank);
    exit_routine(NULL);
  }

  const long seed = (long) val_l(argv[1], rank);
  const long ncside = (long) val_l(argv[2], rank);
  const long long n_par = val_l(argv[3], rank);
  const long n_step = (long) val_l(argv[4], rank);
  if(!(seed*ncside*n_par*n_step)){
    usg_err(rank);
    exit_routine(NULL);
  }

  //get current process number of particles
  loc_n = (n_par%comm_sz && (n_par%comm_sz) > rank) ? (n_par/comm_sz)+1 : (n_par/comm_sz);

  //initialize arrays for gird, particles and particles local to current process
  grid = (cell_t*)calloc(ncside*ncside, sizeof(cell_t));
  loc_par = (particle_t*)calloc(loc_n, sizeof(particle_t));
  if(!rank){
    par = (particle_t*)calloc(n_par, sizeof(particle_t));
    init_particles(seed, ncside, n_par, par, grid);
  }
  
  scatter_particles(n_par, loc_n, loc_par);
  init_grid(ncside);

  MPI_Barrier(MPI_COMM_WORLD);

  update_particles(ncside, loc_par, loc_n, n_step);

  MPI_Barrier(MPI_COMM_WORLD);

  if(!rank){
    t_cx /= t_mass;
    t_cy /= t_mass;
    printf("%.2f %.2f\n", loc_par[0].x, loc_par[0].y);
    printf("%.2f %.2f\n", t_cx, t_cy);
    //printf("%.2f\n", finish-start);
  }

  exit_routine(loc_par);
  return 0;
}
