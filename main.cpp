#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "mpi.h"
#define DEBUG
#define ALPHA 1.0
#define XMAX 100
#define PROCS_PER_DIM 2
#define IMAX 1000
#define DT 1
#define TMAX 100
#define TOUTPUT 10
#define NDIMS 3
using namespace std;

int main(int argc, char **argv)
{
	// get command line parameters

	/* calculation related variables*/
	int imax 			= IMAX; 	// num gridpoints per side per proc
	int dt 				= DT; 	// delta time 
	int tmax 			= TMAX; 	// total num of timesteps
	int toutput 	= TOUTPUT; //timesteps to output results
	int procs_per_dim 	= PROCS_PER_DIM; //num processes per side	

	if (argc==6)
	{
		int imax 			= atoi(argv[1]); 	// num gridpoints per side per proc
		int dt 				= atoi(argv[2]); 	// delta time 
		int tmax 			= atoi(argv[3]); 	// total num of timesteps
		int toutput 	= atoi(argv[4]); //timesteps to output results
		int procs_per_dim 	= atoi(argv[5]); //num processes per side
	}

	int ndims = NDIMS;
	double alpha 	= ALPHA;
	double xmax, ymax, zmax;
	xmax = ymax = zmax = XMAX;
	double boundaryT = 0;
	double dx = xmax / ((double) imax * procs_per_dim);

	/* MPI related Vars */
	int mpi_rank_l;
	int mpi_size;
	int mpi_cart_dim[NDIMS] = {	procs_per_dim, 
													procs_per_dim, 
													procs_per_dim};
	int mpi_cart_period[NDIMS] = {0, 0, 0};
	int mpi_coords_l[NDIMS];

	MPI_Comm mpi_comm;
	// check that parameters make sense
#ifdef DEBUG	
	cout<<"imax: "<<imax<<endl;
	cout<<"dt: "<<dt<<endl;
	cout<<"tmax: "<<tmax<<endl;
	cout<<"toutput: "<<toutput<<endl;
	cout<<"num_proc: "<<procs_per_dim<<endl;
	cout<<"alpha: "<<alpha<<endl;
	cout<<"xmax,ymax,zmax: "<<xmax<<endl;
	cout<<"boundaryT: "<<boundaryT<<endl;
#endif
	// check stability
	if (((alpha*dt)/(dx*dx*dx)) >= (.125))
	{
		cout<<"stability condition failed."<<endl;
		return EXIT_FAILURE;
	}

	//TODO: check that num processors make sense
	// TODO:initialize matrix with gaussian distribution value
	
	// create MPI Cart topology
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_l);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, mpi_cart_dim, 
									mpi_cart_period, true, &mpi_comm);
	MPI_Cart_coords(mpi_comm, mpi_rank_l, ndims, mpi_coords_l);
#ifdef DEBUG
printf("Rank %d coordinates are %d %d %d\n", 
				mpi_rank_l, mpi_coords_l[0], mpi_coords_l[1], mpi_coords_l[2]);
fflush(stdout);
#endif
	MPI_Finalize();
	return EXIT_SUCCESS;
}
// for each processor:
//    make sub_matrix (a 3d Matrix with x+1 size, where x is the interior size)
//		assign data as a function of processor rank
// 		get neighbors of cell

//	for each timestep
//    do calculation
//		
//    if you have neighbor in certain direction, nonblocking send
//			define left gc type/ right gc type
//    	define up gc type down gc type
//			define top gc type bottom gc type
//			non blocking send to each neighbor
// 
//		for each neighbor:
//				nonblocking recieve
//		while messages not yet recieved
//				check messages
//		MPI barrier
//		write results
//		next timestep
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