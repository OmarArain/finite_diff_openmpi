MPCS 51087 Homework 3
=====================
Omar Arain

### The 3D Heat Equation Simulator

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


Analysis:

I ran this analysis using:

sinteractive --partition=broadwl --nodes=64 --ntasks-per-node=16

1. Please see plots/weak_scaling_plot.png.

The weak scaling study shows that the processing time goes up as the number of
processes go up, given that the problem size grows proportionately.  However,
I notice spikes at 27 and 125 processes.  I think this is because 8 processes
can exactly fit on one node, and 64 can fit exactly on 4 nodes, and thus are
able to minimize network traffic. However, 27 and 125 processes are less able
to benefit from that optimization, since it is not clear how they would be
able to efficiently "fit" to minimize the network usage.

Additionally, the I think line is trending up because of the i/o operations
involved - each process has to coordinate with every other process for
collective io - additionally, the larger the problem size, the more data to be
written, and the i/o hardware does not scale up with the number of processes.

2. Please see plots/strong_scaling_plot.png

The strong scaling study shows that the processing time decreases quickly as
the number of processes goes up, but start to slowly increas again around 512
processes.  The trendline initially decreases quickly because of the normal
benefits of parrellism.  I think it starts to increase again because of the
added overhead of i/o - each process does i/o, and at some point the i/o
overhead outweighs the benefits.

3. Please see plots/temp_vs_time_plot.png

Since I ran this for a large problem size, I had to use a small timestep,
which makes changes in the plot hard to see.  The mean temperature decreases
with time, but the standard deviation is increasing slightly, which goes
against my intuition.

4. Please see plots/heat.gif.



WEAK SCALING DATA 
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 8 ./heat 200 .001 100 5 2
mat size: 200, timesteps: 100, num_processors: 8, total_runtime(s): 9.49
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 27 ./heat 300 .001 100 5 3
mat size: 300, timesteps: 100, num_processors: 27, total_runtime(s): 15.22
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 64 ./heat 400 .001 100 5 4
mat size: 400, timesteps: 100, num_processors: 64, total_runtime(s): 10.40
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 125 ./heat 500 .001 100 5 5
mat size: 500, timesteps: 100, num_processors: 125, total_runtime(s): 17.61
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 216 ./heat 600 .001 100 5 6
mat size: 600, timesteps: 100, num_processors: 216, total_runtime(s): 14.18
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 343 ./heat 700 .001 100 5 7
mat size: 700, timesteps: 100, num_processors: 343, total_runtime(s): 15.99
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 512 ./heat 800 .001 100 5 8
mat size: 800, timesteps: 100, num_processors: 512, total_runtime(s): 16.54
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 729 ./heat 900 .001 100 5 9
mat size: 900, timesteps: 100, num_processors: 729, total_runtime(s): 21.68
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 1000 ./heat 900 .001 100 5 10
mat size: 900, timesteps: 100, num_processors: 1000, total_runtime(s): 24.08
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 1000 ./heat 1000 .001 100 5 10
mat size: 1000, timesteps: 100, num_processors: 1000, total_runtime(s): 27.98


STRONG SCALING DATA
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 8 ./heat 400 .001 100 5 2
mat size: 400, timesteps: 100, num_processors: 8, total_runtime(s): 77.48
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 27 ./heat 400 .001 100 5 3
mat size: 400, timesteps: 100, num_processors: 27, total_runtime(s): 23.00
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 64 ./heat 400 .001 100 5 4
mat size: 400, timesteps: 100, num_processors: 64, total_runtime(s): 10.20
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 125 ./heat 400 .001 100 5 5
mat size: 400, timesteps: 100, num_processors: 125, total_runtime(s): 10.18
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 216 ./heat 400 .001 100 5 6
mat size: 400, timesteps: 100, num_processors: 216, total_runtime(s): 7.65
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 343 ./heat 400 .001 100 5 7
mat size: 400, timesteps: 100, num_processors: 343, total_runtime(s): 7.11
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 512 ./heat 400 .001 100 5 8
mat size: 400, timesteps: 100, num_processors: 512, total_runtime(s): 7.43
^[[A^[[A[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 729 ./heat 400 .001 100 5 9
mat size: 400, timesteps: 100, num_processors: 729, total_runtime(s): 11.87
[oarain@midway2-0044 mpcs51087_hw3_oarain]$ mpirun -n 1000 ./heat 400 .001 100 5 10
mat size: 400, timesteps: 100, num_processors: 1000, total_runtime(s): 13.15