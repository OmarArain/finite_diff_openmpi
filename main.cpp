#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include "Matrix.h"
#include "mpi.h"
#include <math.h>
//#define DEBUG
#define ALPHA .0001
#define XMAX 2.
#define PROCS_PER_DIM 2
#define IMAX 100 //1000
#define DT .01
#define TMAX 1000
#define TOUTPUT 1
#define NDIMS 3
#define SENDTAG 999
#define DEBUGRANK 0
#define BOUNDARY_T 0
#define MEAN_T 1.
#define VAR_T .2
using namespace std;

void print_output(Matrix3d<double> &M, int mpi_rank_l,
                   MPI_Datatype * io_subarray, MPI_Datatype * interior,
                   int timestep, MPI_Datatype *io_sub_xslice, MPI_Datatype *xslice,
                   int offset, int mpi_coords_l[], MPI_Comm * mpi_io_comm)
{

  MPI_File mpi_file, mpi_file_io;
  MPI_Status mpi_status;
  // double __test_buf[1000];
  // for(int i=0;i<1000;i++) __test_buf[i]= mpi_rank_l;
  char filename_2d[80], filename_3d[80];
  snprintf(filename_2d, sizeof(filename_2d), "output/heat_output_2d_%d.bin", timestep);
  snprintf(filename_3d, sizeof(filename_3d), "output/heat_output_3d_%d.bin", timestep);
// do stuff
    #ifdef DEBUG
   cout<<mpi_rank_l<<" printing  at time "<<timestep<<endl;  
  #endif
  MPI_File_open(MPI_COMM_WORLD, filename_3d,
                MPI_MODE_CREATE|MPI_MODE_WRONLY,
                MPI_INFO_NULL, &mpi_file);
  MPI_File_set_view(mpi_file, 0, MPI_DOUBLE, *io_subarray, 
                     "native", MPI_INFO_NULL);

  
  MPI_File_write_all(mpi_file, &(M[0]), 1, *interior, &mpi_status);
  MPI_File_close(&mpi_file);
    #ifdef DEBUG
   cout<<mpi_rank_l<<" finished  printing at time "<<timestep<<endl;  
  #endif

  if (mpi_coords_l[0] == 0 )
  {
    #ifdef DEBUG
   cout<<mpi_rank_l<<" printing xslice at time "<<timestep<<endl;  
  #endif
    MPI_File_open(*mpi_io_comm, filename_2d,
                MPI_MODE_CREATE|MPI_MODE_WRONLY,
                MPI_INFO_NULL, &mpi_file_io);
 
    MPI_File_set_view(mpi_file_io, 0, MPI_DOUBLE, *io_sub_xslice, 
                      "native", MPI_INFO_NULL);
    MPI_File_write_all(mpi_file_io, &(M[offset]), 1, *xslice, &mpi_status);
    // MPI_File_write_all(mpi_file_io, &(__test_buf[0]), 1, *xslice, &mpi_status);
    MPI_File_close(&mpi_file_io);
    #ifdef DEBUG
   cout<<mpi_rank_l<<" finished xslice printing at time "<<timestep<<endl;  
  #endif
  }


  //MPI_Barrier(MPI_COMM_WORLD);



}

void exchange_ghost_cells_and_print(
                          int * mpi_nbrs[], int mpi_rank_l, 
                          Matrix3d<double>& data_l, MPI_Request mpi_sreqs[], 
                          MPI_Request mpi_rreqs[], int send_offsets[], int recv_offsets[],
                          MPI_Datatype *mpi_types[], int timestep, int toutput,
                          int mpi_coords_l[], MPI_Comm * mpi_io_comm)
{
  MPI_Datatype *interior = mpi_types[6];
  MPI_Datatype *io_subarray = mpi_types[7];
  MPI_Datatype *io_sub_xslice = mpi_types[8];
  MPI_Datatype *xslice = mpi_types[0];
  #ifdef DEBUG
   cout<<mpi_rank_l<<" sending at time "<<timestep<<endl;  
  #endif
  for (int i = 0; i<6; i++)
    {
      if(*(mpi_nbrs[i]) >= 0 )
      {
        int offset = send_offsets[i];
        MPI_Datatype *mpi_type = mpi_types[i];
        MPI_Isend(&(data_l[offset]), 1, *mpi_type, *(mpi_nbrs[i]),
                  SENDTAG, MPI_COMM_WORLD, &(mpi_sreqs[i]));
      }
      else { mpi_sreqs[i] = MPI_REQUEST_NULL; }
    }

  if ((timestep%toutput) == 0)
      print_output(data_l, mpi_rank_l, io_subarray, interior, 
                    timestep, io_sub_xslice, xslice, 
                    send_offsets[0], mpi_coords_l, mpi_io_comm);
  
  #ifdef DEBUG
   cout<<mpi_rank_l<<" recving at time "<<timestep<<endl;  
  #endif
  for (int i = 0; i<6; i++)
    {
      if(*(mpi_nbrs[i]) >= 0 )
      {
        int offset = recv_offsets[i];
        MPI_Datatype *mpi_type = mpi_types[i];
        MPI_Irecv(&(data_l[offset]), 1, *mpi_type, *(mpi_nbrs[i]),
                  SENDTAG, MPI_COMM_WORLD, &(mpi_rreqs[i]));
      }
      else { mpi_rreqs[i] = MPI_REQUEST_NULL; }
    }
  #ifdef DEBUG
   cout<<mpi_rank_l<<" finished recving at time "<<timestep<<endl;
   // if (mpi_rank_l ==DEBUGRANK)
   // {
   //  for (int i =0; i<6; i++)
   //    cout<<mpi_rank_l<<"sreqs "<<i<<": "<<mpi_sreqs[i]<<endl;
   // }
  #endif
    MPI_Waitall(6, mpi_sreqs, MPI_STATUSES_IGNORE);
  #ifdef DEBUG
   cout<<mpi_rank_l<<" done waiting for send at time "<<timestep<<endl;  
  #endif
    MPI_Waitall(6, mpi_rreqs, MPI_STATUSES_IGNORE);
    #ifdef DEBUG
   cout<<mpi_rank_l<<" done waiting for recv at time "<<timestep<<endl;  
  #endif
    #ifdef DEBUG
     cout<<mpi_rank_l<<" DONE with step "<<timestep<<endl; 
    #endif
    MPI_Barrier(MPI_COMM_WORLD); 
}


inline void calculate_matrix_test(Matrix3d<double> &M, int mpi_rank_l)
{
  int xmax = M.xsize();
  int ymax = M.ysize();
  int zmax = M.zsize();
  for(int x=0; x<xmax; ++x)
    for(int y=0; y<ymax; ++y)
      for(int z=0; z<zmax; ++z)
        M(x,y,z) += mpi_rank_l;
}
int main(int argc, char **argv)
{
  /* calculation related variables*/
  int imax      = IMAX;                // num gridpoints per side
  double dt     = DT;                  // delta time 
  int tmax      = TMAX;                // total num of timesteps
  int toutput   = TOUTPUT;             //timesteps to output results
  int procs_per_dim   = PROCS_PER_DIM; //num processes per side
  double mean_t= MEAN_T;
  double var_t = VAR_T;


  // Get command line parameters
  if (argc==6)
  {
    // cout<<"ghi"<<endl;
    imax      = atoi(argv[1]);         // num gridpoints per side per proc
    dt        = atof(argv[2]);         // delta time 
    tmax      = atoi(argv[3]);         // total num of timesteps
    toutput   = atoi(argv[4]);         //timesteps to output results
    procs_per_dim   = atoi(argv[5]);   //num processes per side
  }

  /* more calculation realated vars*/
  int ndims = NDIMS;
  int imax_l = imax/procs_per_dim;
  double alpha  = ALPHA;
  double xmax, ymax, zmax;
  xmax = ymax = zmax = XMAX;
  double boundary_t = BOUNDARY_T;
  double dx = xmax / ((double) imax );
  double start_time, end_time;
 
  //wrapper over std::vector, contiguous, ROW MAJOR
  Matrix3d<double> data_l(imax_l+2, imax_l+2, imax_l+2); 

  // Check stability
  if ( ((alpha*dt)/(dx*dx*dx)) >= (.125))
  {
    cout<<"Stability condition failed: "<<((alpha*dt)/(dx*dx*dx));
    cout<<"> .125"<<endl;
    return EXIT_FAILURE;
  }


  /* Initialize MPI */
  if (argc==6) MPI_Init(&argc, &argv);
  else MPI_Init(NULL, NULL);
  /* MPI related vars */
  int mpi_rank_l;
  int mpi_size;
  MPI_Comm mpi_comm, mpi_io_comm;
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
  // datatypes
  MPI_Aint lb, sz_dbl;
  MPI_Datatype xslice, yslice, zslice, 
                interior, io_subarray,
                 zchunk, io_sub_xslice;
  MPI_Datatype * mpi_types[9] = { &xslice, &xslice,
                                  &yslice, &yslice,
                                  &zslice, &zslice,
                                  &interior, &io_subarray,
                                  &io_sub_xslice};
  //offsets of ACTUAL BORDERS OF ARRAY
  int send_offsets[6] = { (imax_l+2) * (imax_l+2) * (imax_l+2-2)+1+imax_l, // ydim*zdim*xrow
                          (imax_l+2) * (imax_l+2) * (1) + 1+imax_l,          // ydim*zdim*xrow
                          (imax_l+2) * (imax_l+2-2) + (imax_l+2) * (imax_l+2) + 1,              //zdim*yrow
                          (imax_l+2) * (1)+ (imax_l+2) * (imax_l+2) + 1,                       //zdim*yrow
                          (imax_l+2-2) + (imax_l+2)*(imax_l+2)+(imax_l+2),                           //zrow
                          (1) + (imax_l+2)*(imax_l+2)+(imax_l+2)};                                   //zrow+1
  //offsets of GHOST CELLS OF ARRAY
  int recv_offsets[6]  = { (imax_l+2) * (imax_l+2) * (imax_l+2-1)+1+imax_l, // ydim*zdim*xrow
                          (imax_l+2) * (imax_l+2) * (0)+1+imax_l,          // ydim*zdim*xrow
                          (imax_l+2) * (imax_l+2-1) + (imax_l+2) * (imax_l+2) + 1,              //zdim*yrow
                          (imax_l+2) * (0)+ (imax_l+2) * (imax_l+2) + 1,                       //zdim*yrow
                          (imax_l+2-1) + (imax_l+2)*(imax_l+2)+(imax_l+2),                           //zrow
                          (0)+ (imax_l+2)*(imax_l+2)+(imax_l+2)};                                   //zrow
  // send requests
  MPI_Request mpi_sreqs[6] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                              MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                              MPI_REQUEST_NULL, MPI_REQUEST_NULL};
  //recv requests
  MPI_Request mpi_rreqs[6] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                              MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                              MPI_REQUEST_NULL, MPI_REQUEST_NULL};

  start_time = MPI_Wtime();
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_l);
  MPI_Type_get_extent (MPI_DOUBLE, &lb, &sz_dbl);

  #ifdef DEBUG
   cout<<mpi_rank_l<<" INITIALIZED"<<endl;  
  #endif
  // Create MPI Cartesian topology
  MPI_Cart_create(MPI_COMM_WORLD, ndims, mpi_cart_dim, 
                  mpi_cart_period, true, &mpi_comm);
  MPI_Cart_coords(mpi_comm, mpi_rank_l, ndims, mpi_coords_l);
  //color CANNOT BE NEGATIVE
  MPI_Comm_split(MPI_COMM_WORLD, mpi_coords_l[0]+1, mpi_rank_l, &mpi_io_comm);
  // get neighbors
  MPI_Cart_shift(mpi_comm, 0, 1, &mpi_nbr_xdn_l, &mpi_nbr_xup_l);
  MPI_Cart_shift(mpi_comm, 1, 1, &mpi_nbr_ydn_l, &mpi_nbr_yup_l);
  MPI_Cart_shift(mpi_comm, 2, 1, &mpi_nbr_zdn_l, &mpi_nbr_zup_l);

  // Create slice datatypes, interior, and io_subarray
  // and io_sub_xslice for filewriting

  MPI_Type_vector(imax_l, imax_l, (imax_l+2), 
                  MPI_DOUBLE, &xslice);
  MPI_Type_vector(imax_l, imax_l, (imax_l+2)*(imax_l+2),
                  MPI_DOUBLE, &yslice);
  MPI_Type_vector(imax_l, 1, imax_l+2, MPI_DOUBLE, &zchunk);
  MPI_Type_hvector(imax_l, 1, (imax_l+2)*(imax_l+2)*sz_dbl, 
                  zchunk, &zslice);

  int global_sizes_i[3] = {imax_l+2, imax_l+2, imax_l+2};
  int local_sizes_i[3] = {imax_l, imax_l, imax_l};
  int starts_i[3] = {1, 1, 1};
  MPI_Type_create_subarray(3, global_sizes_i, local_sizes_i, 
                          starts_i, MPI_ORDER_C, 
                          MPI_DOUBLE, &interior);

  int global_sizes_f[3] = { imax_l*procs_per_dim, 
                            imax_l*procs_per_dim, 
                            imax_l*procs_per_dim };
  int local_sizes_f[3] = {imax_l, imax_l, imax_l};
  int starts_f[3] = { mpi_coords_l[0]*imax_l, 
                      mpi_coords_l[1]*imax_l, 
                      mpi_coords_l[2]*imax_l};
  MPI_Type_create_subarray(3, global_sizes_f, local_sizes_f, 
                          starts_f, MPI_ORDER_C, 
                          MPI_DOUBLE, &io_subarray);
  int global_sizes_xs[2] = { imax_l*procs_per_dim, 
                            imax_l*procs_per_dim};
  int local_sizes_xs[2] = {imax_l, imax_l};
  int starts_xs[2] = { mpi_coords_l[2]*imax_l, 
                       mpi_coords_l[1]*imax_l}; 
  MPI_Type_create_subarray(2, global_sizes_xs, local_sizes_xs, 
                          starts_xs, MPI_ORDER_C, 
                          MPI_DOUBLE, &io_sub_xslice);

  MPI_Type_commit(&xslice);
  MPI_Type_commit(&yslice);
  MPI_Type_commit(&zslice);
  MPI_Type_commit(&interior);
  MPI_Type_commit(&io_subarray);
  MPI_Type_commit(&io_sub_xslice);


//FILE IO TEST

  // check that num processors make sense
  if (mpi_size!=(procs_per_dim*procs_per_dim*procs_per_dim))
  {
    printf("MPI size: %d does not match procs_per_dim:%d^3 \n", 
            mpi_size, procs_per_dim);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  //initialize matrix
  int _i_xstart, _i_ystart, _i_zstart;
  _i_xstart = mpi_coords_l[0]*imax_l;
  _i_ystart = mpi_coords_l[1]*imax_l;
  _i_zstart = mpi_coords_l[2]*imax_l;
  data_l.set_gaussian_interior(mean_t, var_t, dx, 
                                _i_xstart, _i_ystart, _i_zstart);
  data_l.reset_boundaries(boundary_t);
  #ifdef DEBUG
   cout<<mpi_rank_l<<" MATRIX INITIALIZED"<<endl;  
  #endif
  /* CALC LOOP */
  for(int t = 0; t<tmax; t++)
  {
    #ifdef DEBUG
    cout<<mpi_rank_l<<" started calc "<<t<<endl;  
    #endif
    exchange_ghost_cells_and_print(
                          mpi_nbrs, mpi_rank_l, 
                          data_l, mpi_sreqs, 
                          mpi_rreqs, send_offsets, recv_offsets,
                          mpi_types, t, toutput, mpi_coords_l,
                          &mpi_io_comm);
    data_l.calc_heat_equation(dx, dt, alpha);

  }
  end_time = MPI_Wtime();
  if(mpi_rank_l==0)
  printf("mat size: %d, timesteps: %d, num_processors: %d, total_runtime(s): %.2f\n", 
          imax, tmax, mpi_size, (double)(end_time-start_time));
  MPI_Finalize();
  return EXIT_SUCCESS;
 }