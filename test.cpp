#include "Matrix.h"
#include <iostream>
#include <math.h>
using namespace std;

int main()
{
	Matrix3d<double> M(20,20,20);
	M(0,0,0) = 0;
	M(0,0,1) = 1;
	M(0,1,0) = 2;
	M(0,1,1) = 3;
	M(1,0,0) = 4;
	M(1,0,1) = 5;
	M(1,1,0) = 6;
	M(1,1,1) = 7;
	for(int i=0; i<8; i++) cout<<M[i]<<", "<<endl;

	M.set_gaussian_interior(10., .2, 1., 5,5,5);
	M.reset_boundaries(-1);
	print_xslice(M, 0);
}


/*
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
*/