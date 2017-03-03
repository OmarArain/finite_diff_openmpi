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

4.1 Program Requirements
1. The program may assume the following:
(a) A cubic domain (xmax = ymax = zmax) (b) A uniform grid (∆x = ∆y = ∆z)
(c) The same maximum bounds for each dimension (xmax = ymax = zmax) (d) A constant T = 0 at all boundaries
(e) An initial Gaussian heat distribution (see above)
2. Your program must take the following as runtime parameters:
(a) Either the grid spacing (∆x) or the number of gridpoints (imax) (b) The size of each timestep (∆t)
(c) The total number of timesteps (tmax)
(d) A number of timesteps, toutput, to output temperature data (see below)
(e) The number of MPI processes
3. Your program must write the following output to file(s) every (toutput) timesteps:
(a) A 2D slice of the temperature data
(b) The mean and standard deviation of the temperature in the entire domain
4. Your program must meet these additional requirements
(a) No process may allocate memory for more than its subdomain
(b) The simulation should not proceed if the stability condition is not met
4.2 Analysis
Along with your assignment, you must the submit the following:
1. A strong scaling study for a large problem, including a graph and explanation
2. A weak scaling study for a large problem, including a graph and explanation
3. A graph of mean temperature with standard deviation vs. time for a large problem.
4. A 2D animation of temperature vs. time for a problem run in parallel. This can be significantly smaller than the problems simulated in the other analyses. You can either submit a graphic file (an MPEG, a GIF, etc) or a script that runs an animation from some output that you also submit.