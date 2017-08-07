3D Heat Equation Simulator
=====================
Omar Arain

Please see main.cpp and Matrix.h for my heat diffusion
simulation.

To build and run:

make heat 
mpirun -n <num_nodes> ./heat <gridpoints_per_dim> <dt>
					<max_timesteps> <toutput> <procs_per_dim>
e.g.
mpirun -n 8 ./heat 100 .01 25 5 2

would run the simulation with 100 gridpoints per side (100x100x100),
with dt = .01 and 100 total timesteps (i.e. t ranges from 0 to 1),
and would output/write to file every 5 timesteps.  Additionally,
it would have 2 processors per dimension (2*2*2 procs total).

image.py is used  to create images from the output files, 
and puts them in the images directory.

mean_var.py is used to generate lists of means and 
variances from the output files.

Both of these python scripts are "fragile" and have parameters
that need to be changed depending on what files are in the 
output directory.

If the processors per dim and the mpi num_nodes parameters are
not consistent, the program will exit.

The real range of the matrix is 0 to 2 for all dimensions,
and the alpha value is .000001.  These are coded as macro
#defines in main.cpp (XMAX and ALPHA).

The program outputs a 2d slice of the matrix at every toutput
timestep as a seperate .bin file in the output directory. 
These are named "heat_output_2d_T.bin", where T is the timestep.
Additioanlly, it also creates heat_output_mean.bin and heat_output_var.bin

The program can also output the full global matrix at every toutput,
if it is compiled with -DOUTPUT3D.  However this is quite memory and
i/o intensive for even modest size matrices (200x200x200).

My 3d matrices output nice looking values, however my 2d slices 
contain some artifacts along some of the borders of some of the
submatrices, which I have yet to resolve.
