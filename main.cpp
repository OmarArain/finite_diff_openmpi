#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "Matrix.h"
#include "mpi.h"
#define DEBUG
#define ALPHA .00001
#define XMAX 100
#define PROCS_PER_DIM 2
#define IMAX 2 //1000
#define DT 1
#define TMAX 2
#define TOUTPUT 10
#define NDIMS 3
#define SENDTAG 999
using namespace std;


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

  //send requests, declared like this for clarity
  MPI_Request mpi_sreqs[6] = {MPI_REQUEST_NULL};
  // MPI_Request &mpi_sreq_xup = mpi_sreqs[0], &mpi_sreq_xdn = mpi_sreqs[1],
  //             &mpi_sreq_yup = mpi_sreqs[2], &mpi_sreq_ydn = mpi_sreqs[3],
  //             &mpi_sreq_zup = mpi_sreqs[4], &mpi_sreq_zdn = mpi_sreqs[5];
  //recv requests
  MPI_Request mpi_rreqs[6] = {MPI_REQUEST_NULL};
  // MPI_Request &mpi_rreq_xup = mpi_rreqs[0], &mpi_rreq_xdn = mpi_rreqs[1],
  //             &mpi_rreq_yup = mpi_rreqs[2], &mpi_rreq_ydn = mpi_rreqs[3],
  //             &mpi_rreq_zup = mpi_rreqs[4], &mpi_rreq_zdn = mpi_rreqs[5];

  // int mpi_test_send = mpi_rank_l;
  // int mpi_msg_send_pending = 0;
  // int mpi_msg_recv_pending = 0;

  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_l);

  // check that num processors make sense
  if (mpi_size!=(procs_per_dim*procs_per_dim*procs_per_dim))
  {
    printf("MPI size: %d does not match procs_per_dim:%d^3 \n", 
            mpi_size, procs_per_dim);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

#ifdef DEBUG  
  if(mpi_rank_l == 0)
  {
    cout<<"imax: "<<imax<<endl;
    cout<<"dt: "<<dt<<endl;
    cout<<"tmax: "<<tmax<<endl;
    cout<<"toutput: "<<toutput<<endl;
    cout<<"num_proc: "<<procs_per_dim<<endl;
    cout<<"alpha: "<<alpha<<endl;
    cout<<"xmax,ymax,zmax: "<<xmax<<endl;
    cout<<"boundaryT: "<<boundaryT<<endl;
  }
#endif

  // Create MPI Cartesian topology
  MPI_Cart_create(MPI_COMM_WORLD, ndims, mpi_cart_dim, 
                  mpi_cart_period, true, &mpi_comm);
  MPI_Cart_coords(mpi_comm, mpi_rank_l, ndims, mpi_coords_l);
  // get neighbors
  MPI_Cart_shift(mpi_comm, 0, 1, &mpi_nbr_xdn_l, &mpi_nbr_xup_l);
  MPI_Cart_shift(mpi_comm, 1, 1, &mpi_nbr_ydn_l, &mpi_nbr_yup_l);
  MPI_Cart_shift(mpi_comm, 2, 1, &mpi_nbr_zdn_l, &mpi_nbr_zup_l);

#ifdef DEBUG
  printf("mpi_size: %d\n", mpi_size);
  printf("Rank %d coordinates are %d %d %d\n", 
          mpi_rank_l, mpi_coords_l[0],
          mpi_coords_l[1], mpi_coords_l[2]);
  printf("neighbors: xup: %d xdn: %d yup: %d ydn: %d zup: %d zdn: %d\n",
        mpi_nbr_xup_l, mpi_nbr_xdn_l, 
        mpi_nbr_yup_l, mpi_nbr_ydn_l,
        mpi_nbr_zup_l, mpi_nbr_zdn_l);
  fflush(stdout);
#endif
  //  make sub_matrix (a 3d Matrix with x+2 size, where x is the interior size)
  //  TODO: assign gaussian data as a function of processor rank/coords
  // TODO: make array of ghost cell MPI_TYPES for each data_l (i.e. each surface of the cube)
  data_ix_l = 0;
  for(int x=0; x<(imax+2); ++x)
  {
    for(int y=0; y<(imax+2); ++y)
    {
      for(int z=0; z<(imax+2); ++z)
      {
        data_l[data_ix_l] = mpi_rank_l;
        ++data_ix_l;
      }
    }
  }

  /* CALC LOOP */
  for(int t = 0; t<tmax; t+=dt)
  {
    // TODO: do calculation
    for(int x=0; x<(imax+2); ++x)
    {
      for(int y=0; y<(imax+2); ++y)
      {
        for(int z=0; z<(imax+2); ++z)
        {
          //data_l[data_ix_l] += 100;
          data_l(x,y,z) += 100;
        }
      }
    }
    // If you have neighbor in certain direction do nonblocking send
    for (int i = 0; i<6; i++)
    {
      if(*(mpi_nbrs[i]) > 0 )
      {
        MPI_Isend(&(data_l[0]), 1, MPI_DOUBLE, *(mpi_nbrs[i]),
                  SENDTAG, MPI_COMM_WORLD, &(mpi_sreqs[i]));
        //++mpi_msg_send_pending;
      }
      else
      {
        mpi_sreqs[i] = MPI_REQUEST_NULL;
      }
    }
    for (int i = 0; i<6; i++)
    {
      if(*(mpi_nbrs[i]) > 0 )
      {  
        MPI_Irecv(&(data_l[i+6]), 1, MPI_DOUBLE, *(mpi_nbrs[i]),
                  SENDTAG, MPI_COMM_WORLD, &(mpi_rreqs[i]));
        //++mpi_msg_recv_pending;
      }
      else
      {
        mpi_sreqs[i] = MPI_REQUEST_NULL;
      }
    } 

    // if t is a multiple of touput, write to file (not the ghost cells!)
    //    write results
#ifdef DEBUG
    if(mpi_rank_l == 0)
    {
      for(int x=0; x<(imax+2); ++x)
      {
        for(int y=0; y<(imax+2); ++y)
        {
          for(int z=0; z<(imax+2); ++z)
          {
            //data_l[data_ix_l] += 100;
            cout<<data_l(x,y,z)<<",";
          }
          cout<<";"<<endl;;
        }
        cout<<";";
      }
      cout<<endl;
    }

#endif
    //    while messages not yet recieved
    //        check messages
#ifdef DEBUG
    if(mpi_rank_l == 0) cout<<"Waiting for sends..."<<endl; 
#endif
    MPI_Waitall(6, mpi_sreqs, MPI_STATUSES_IGNORE);
#ifdef DEBUG
    if(mpi_rank_l == 0) cout<<"Waiting for recvs..."<<endl; 
#endif
    MPI_Waitall(6, mpi_rreqs, MPI_STATUSES_IGNORE);
    //    MPI barrier?
  }
    //    

    MPI_Finalize();
    return EXIT_SUCCESS;
 }
/*
(a) A cubic domain (xmax = ymax = zmax) (b) A uniform grid (∆x = ∆y = ∆z)
(c) The same maximum bounds for each dimension (xmax = ymax = zmax) (d) A constant T = 0 at all boundaries
(e) An initial Gaussian heat distribution (see above)
2. Your program must take the following as runtime parameters:
(a) Either the grid spacing (∆x) or the number of gridpoints (imax) (b) The size of each timestep (∆t)
(c) The total number of timesteps (tmax)
(d) A number of timesteps, toutput, to output temperature data (see below)
(e) The number of MPI processes 
*/
