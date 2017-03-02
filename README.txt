MPCS 51087 Homework 3
=====================
Omar Arain

### The 3D Heat Equation
- See hw3.pdf

### MPI Topologies
MPI offers several library functions for decomposing a spatial domain and exchanging ghost cells. Section 4.2 in Using MPI by Gropp, Lusk, and Skjellum offers a very detailed example for using these subroutines. Example programs are also available on the book’s website ( link ); note, in particular, the “Jacobi with two-dimensional decomposition” examples.
The most important subroutines are MPI Cart create and MPI Cart shift.



//define MPI_TYPE

typedef struct CalcGrid3d
{
	double ** grid_cells;

};


rank to x_val, y_val, z_val - use MPI_Cart_coords to get coords

x_coord_mpi 

matrix of size x_dim, y_dim, z_dim;

x_coord = rank % y_dim;


// get command line parameters
// check that parameters make sense
// check stability
// check that num processors make sense
// initialize matrix with gaussian distribution value
// create MPI Cart topology
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

)