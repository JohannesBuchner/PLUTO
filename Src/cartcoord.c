#include

GetRankCoordinates (periods, nprocs, coords);

/* ********************************************************************* */
void GetRankCoordinates (int periods[3], int nprocs[3], int coords[3])
/*
 *
 *
 *********************************************************************** */
{
   MPI_Comm cartcomm;
   MPI_Status status;
   int nprocs[3], periods[3], coords[3];
   int rank;

   MPI_Comm_rank(MPI_COMM_WORLD, &prank);

   coords[0]  = coords[1]  = coords[2]  = 0;
   periods[0] = periods[1] = periods[2] = 0;
   nprocs[0]  = nprocs[1]  = nprocs[2]  = 0;

   AL_Get_cart_comm(SZ, &cartcomm);
   MPI_Cart_get(cartcomm, 3, nprocs, periods, coords);
   MPI_Cart_rank(cartcomm, coords, &rank);
}


/* ********************************************************************* */
int CoordOffRank(Grid *grid, int dir, int poff, MPI_comm *cartcomm)
/*
 * PURPOSE
 *
 *  Obtain the rank associated to the processor with 
 *  Cartesian grid coordinates 
 *  [ix, iy, iz] + poff*[dir==IDIR, dir==JDIR, dir==KDIR] 
 *
 *********************************************************************** */
{
  int coords[3], rnk;

  coords[IDIR] = grid[IDIR].rank_coord;
  coords[JDIR] = grid[JDIR].rank_coord;
  coords[KDIR] = grid[KDIR].rank_coord;

  coords[dir] += poff;

  MPI_Cart_rank(cartcomm, coords, &rnk);
  return rnk;
}