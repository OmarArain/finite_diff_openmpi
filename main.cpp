#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "Matrix.h"
#include "mpi.h"
#define DEBUG
#define ALPHA .00001
#define XMAX 100
#define PROCS_PER_DIM 2
#define IMAX 3 //1000
#define DT 1
#define TMAX 2
#define TOUTPUT 10
#define NDIMS 3
#define SENDTAG 999
#define RECVTAG 1000
#define DEBUGRANK 4
using namespace std;

inline void initialize_matrix_test(Matrix3d<double>& M, int xstart, int ystart, int zstart)
{
  int xmax = M.xsize();
  int ymax = M.ysize();
  int zmax = M.zsize();
  int data_ix_l = 0;
  for(int x=0; x<xmax; ++x)
  {
    for(int y=0; y<ymax; ++y)
    {
      for(int z=0; z<zmax; ++z)
      {
        M(x,y,z) = data_ix_l++;
      }
    }
  }
}

void exchange_ghost_cells(int * mpi_nbrs[6], int mpi_rank_l, Matrix3d<double>& data_l, 
                          MPI_Request mpi_sreqs[6], MPI_Request mpi_rreqs[6])
{
  for (int i = 0; i<6; i++)
    {
      if(*(mpi_nbrs[i]) >= 0 )
      {
        #ifdef DEBUG
        if(mpi_rank_l==DEBUGRANK) cout<<"Sending to *(mpi_nbrs[i]:"<<*(mpi_nbrs[i])<<endl;
        #endif
        MPI_Isend(&(data_l[i]), 1, MPI_DOUBLE, *(mpi_nbrs[i]),
                  SENDTAG, MPI_COMM_WORLD, &(mpi_sreqs[i]));
        //++mpi_msg_send_pending;
      }
      else { mpi_sreqs[i] = MPI_REQUEST_NULL; }
    }

  for (int i = 0; i<6; i++)
    {
      if(*(mpi_nbrs[i]) >= 0 )
      {
        #ifdef DEBUG
        if(mpi_rank_l==DEBUGRANK) cout<<"Recving to *(mpi_nbrs[i]:"<<*(mpi_nbrs[i])<<endl;
        #endif  
        MPI_Irecv(&(data_l[i+6]), 1, MPI_DOUBLE, *(mpi_nbrs[i]),
                  SENDTAG, MPI_COMM_WORLD, &(mpi_rreqs[i]));
        //++mpi_msg_recv_pending;
      }
      else { mpi_rreqs[i] = MPI_REQUEST_NULL; }
    }
    // #ifdef DEBUG
    // if(mpi_rank_l == DEBUGRANK) cout<<mpi_rank_l<<" Waiting for sends..."<<endl; 
    // #endif
    MPI_Waitall(6, mpi_sreqs, MPI_STATUSES_IGNORE);
    // #ifdef DEBUG
    // if(mpi_rank_l == DEBUGRANK) cout<<mpi_rank_l<<" Waiting for recvs..."<<endl;  
    // #endif
    MPI_Waitall(6, mpi_rreqs, MPI_STATUSES_IGNORE);
    #ifdef DEBUG
     cout<<mpi_rank_l<<" DONE"<<endl;  
    #endif
    MPI_Barrier(MPI_COMM_WORLD); 
}


inline void calculate_matrix_test(Matrix3d<double> &M)
{
  int xmax = M.xsize();
  int ymax = M.ysize();
  int zmax = M.zsize();
  int data_ix_l = 0;
  for(int x=0; x<xmax; ++x)
    for(int y=0; y<ymax; ++y)
      for(int z=0; z<zmax; ++z)
        M(x,y,z) = 100+data_ix_l++;
}
int main(int argc, char **argv)
{
  /* calculation related variables*/
  int imax      = IMAX;                // num gridpoints per side per proc
  int dt        = DT;                  // delta time 
  int tmax      = TMAX;                // total num of timesteps
  int toutput   = TOUTPUT;             //timesteps to output results
  int procs_per_dim   = PROCS_PER_DIM; //num processes per side 

  // Get command line parameters
  if (argc==6)
  {
    imax      = atoi(argv[1]);         // num gridpoints per side per proc
    dt        = atoi(argv[2]);         // delta time 
    tmax      = atoi(argv[3]);         // total num of timesteps
    toutput   = atoi(argv[4]);         //timesteps to output results
    procs_per_dim   = atoi(argv[5]);   //num processes per side
  }

  /* more calculation realated vars*/
  int ndims = NDIMS;
  double alpha  = ALPHA;
  double xmax, ymax, zmax;
  xmax = ymax = zmax = XMAX;
  double boundaryT = 0;
  double dx = xmax / ((double) imax * procs_per_dim);
  Matrix3d<double> data_l(imax+2, imax+2, imax+2); //wrapper over std::vector, contiguous, ROW MAJOR
  Matrix3d<double> ghost_cells[6]; //array of ghost cell arrays
  int data_ix_l = 0;

  // Check stability
  if (((alpha*dt)/(dx*dx*dx)) >= (.125))
  {
    cout<<"Stability condition failed: "<<((alpha*dt)/(dx*dx*dx));
    cout<<"< .125"<<endl;
    return EXIT_FAILURE;
  }
  /* MPI related vars */
  int mpi_rank_l;
  int mpi_size;
  MPI_Comm mpi_comm;
  int mpi_cart_dim[NDIMS] = { procs_per_dim, 
                          procs_per_dim, 
                          procs_per_dim};
  int mpi_cart_period[NDIMS] = {0, 0, 0};
  int mpi_coords_l[NDIMS];

  //rank of neighbors in each direction
  int mpi_nbr_xup_l, mpi_nbr_xdn_l, 
      mpi_nbr_yup_l, mpi_nbr_ydn_l, 
      mpi_nbr_zup_l, mpi_nbr_zdn_l;
  int * mpi_nbrs[6] = { &mpi_nbr_xup_l, &mpi_nbr_xdn_l, 
                        &mpi_nbr_yup_l, &mpi_nbr_ydn_l, 
                        &mpi_nbr_zup_l, &mpi_nbr_zdn_l };

  // send requests
  MPI_Request mpi_sreqs[6] = {MPI_REQUEST_NULL};
  //recv requests
  MPI_Request mpi_rreqs[6] = {MPI_REQUEST_NULL};

  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_l);


  // Create border datatypes 
  MPI_Datatype border_xup, border_xdn,
               border_yup, border_ydn,
               border_zup, border_zdn;
  int mpi_subarray_bigsizes[3] = {imax, imax, imax};
  //xupborder
  {
    int mpi_subarray_subsizes[3] =  {1, imax, imax};
    int mpi_subarray_starts[3]   =  {imax-1, 0, 0};
    MPI_Type_create_subarray(ndims, mpi_subarray_bigsizes, 
                          mpi_subarray_subsizes, mpi_subarray_starts, 
                          MPI_ORDER_C, MPI_DOUBLE, &border_xup);
    MPI_Type_commit(&border_xup);
  }
  // xdnborder
  {
    int mpi_subarray_subsizes[3] = {1, imax, imax};
    int mpi_subarray_starts[3]   = {0, 0, 0};
    MPI_Type_create_subarray(ndims, mpi_subarray_bigsizes, 
                          mpi_subarray_subsizes, mpi_subarray_starts, 
                          MPI_ORDER_C, MPI_DOUBLE, &border_xdn);
    MPI_Type_commit(&border_xdn);
  }
  //yupborder
  {
    int mpi_subarray_subsizes[3] = {imax, 1, imax};
    int mpi_subarray_starts[3]   = {0, imax-1, 0 };
    MPI_Type_create_subarray(ndims, mpi_subarray_bigsizes, 
                          mpi_subarray_subsizes, mpi_subarray_starts, 
                          MPI_ORDER_C, MPI_DOUBLE, &border_yup);
    MPI_Type_commit(&border_yup);
  }
  // ydnborder
  {
    int mpi_subarray_subsizes[3] = {imax, 1, imax};
    int mpi_subarray_starts[3]   = {0, 0, 0};
    MPI_Type_create_subarray(ndims, mpi_subarray_bigsizes, 
                          mpi_subarray_subsizes, mpi_subarray_starts, 
                          MPI_ORDER_C, MPI_DOUBLE, &border_ydn);
    MPI_Type_commit(&border_ydn);
  }
  // zupborder
  {
    int mpi_subarray_subsizes[3] = {imax, imax, 1};
    int mpi_subarray_starts[3]   = {0, 0, imax-1};
    MPI_Type_create_subarray(ndims, mpi_subarray_bigsizes, 
                          mpi_subarray_subsizes, mpi_subarray_starts, 
                          MPI_ORDER_C, MPI_DOUBLE, &border_zup);
    MPI_Type_commit(&border_zup);
  }
  // zdnborder
  {
    int mpi_subarray_subsizes[3] = {imax, imax, 1};
    int mpi_subarray_starts[3]   = {0, 0, 0};
    MPI_Type_create_subarray(ndims, mpi_subarray_bigsizes, 
                          mpi_subarray_subsizes, mpi_subarray_starts, 
                          MPI_ORDER_C, MPI_DOUBLE, &border_zdn);
    MPI_Type_commit(&border_zdn);
  }

  // check that num processors make sense
  if (mpi_size!=(procs_per_dim*procs_per_dim*procs_per_dim))
  {
    printf("MPI size: %d does not match procs_per_dim:%d^3 \n", 
            mpi_size, procs_per_dim);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

// #ifdef DEBUG  
//   if(mpi_rank_l == DEBUGRANK)
//   {
//     cout<<"imax: "<<imax<<endl;
//     cout<<"dt: "<<dt<<endl;
//     cout<<"tmax: "<<tmax<<endl;
//     cout<<"toutput: "<<toutput<<endl;
//     cout<<"num_proc: "<<procs_per_dim<<endl;
//     cout<<"alpha: "<<alpha<<endl;
//     cout<<"xmax,ymax,zmax: "<<xmax<<endl;
//     cout<<"boundaryT: "<<boundaryT<<endl;
//   }
// #endif

  // Create MPI Cartesian topology
  MPI_Cart_create(MPI_COMM_WORLD, ndims, mpi_cart_dim, 
                  mpi_cart_period, true, &mpi_comm);
  MPI_Cart_coords(mpi_comm, mpi_rank_l, ndims, mpi_coords_l);
  // get neighbors
  MPI_Cart_shift(mpi_comm, 0, 1, &mpi_nbr_xdn_l, &mpi_nbr_xup_l);
  MPI_Cart_shift(mpi_comm, 1, 1, &mpi_nbr_ydn_l, &mpi_nbr_yup_l);
  MPI_Cart_shift(mpi_comm, 2, 1, &mpi_nbr_zdn_l, &mpi_nbr_zup_l);

  // #ifdef DEBUG
  // printf("mpi_size: %d\n", mpi_size);
  // printf("Rank %d coordinates are %d %d %d\n", 
  //         mpi_rank_l, mpi_coords_l[0],
  //         mpi_coords_l[1], mpi_coords_l[2]);
  // printf("neighbors: xup: %d xdn: %d yup: %d ydn: %d zup: %d zdn: %d\n",
  //       mpi_nbr_xup_l, mpi_nbr_xdn_l, 
  //       mpi_nbr_yup_l, mpi_nbr_ydn_l,
  //       mpi_nbr_zup_l, mpi_nbr_zdn_l);
  // fflush(stdout);
  // #endif
  //initialize matrix
  initialize_matrix_test(data_l, 1, 1, 1);

  /* CALC LOOP */
  for(int t = 0; t<tmax; t+=dt)
  {
    exchange_ghost_cells(mpi_nbrs, mpi_rank_l, data_l, mpi_sreqs, mpi_rreqs);
    calculate_matrix_test(data_l);
  }
    // if t is a multiple of touput, write to file (not the ghost cells!)
    //    write results

    MPI_Finalize();
    return EXIT_SUCCESS;
 }