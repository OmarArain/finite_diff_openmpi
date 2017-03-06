MPCS 51087 Homework 3
=====================
Omar Arain

### The 3D Heat Equation

on is not met
4.2 Analysis
Along with your assignment, you must the submit the following:
1. A strong scaling study for a large problem, including a graph and explanation


2. A weak scaling study for a large problem, including a graph and explanation
3. A graph of mean temperature with standard deviation vs. time for a large problem.
4. A 2D animation of temperature vs. time for a problem run in parallel. This can be significantly smaller than the problems simulated in the other analyses. You can either submit a graphic file (an MPEG, a GIF, etc) or a script that runs an animation from some output that you also submit.




Analysis:
1. 
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


2.

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