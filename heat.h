#ifndef _HEAT_H
#define _HEAT_H

//initialize 3d matrix with continuous block of data in memory
double *** matrix_3d(int x, int y, int z);
//free matrix_3d and all associated data
void matrix_3d_free(double *** M);

//populate 3d matrix of with values from a gaussian distribution
void gaussian_3d(int dim, double mean, double var, double*** M, 
				int xsize, int ysize, int zsize);

//return a 2d slice of the 3d matrix.
double ** matrix_3d_slice(int xsize, int ysize, int zsize, 
						  int x_dim, int y_dim, int z_dim);



#endif