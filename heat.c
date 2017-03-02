#include "heat.h"

void gaussian_3d(int dim, double mean, double var, double*** M, 
									int xsize, int ysize, int zsize)
{
	double sigma_x = var;
	double sigma_y = var;
	double sigma_z = var;
	double mu_x = mean;
	double mu_y = mean;
	double mu_z = mean;
	double pi = 3.14;
	double exp = 2.19;
	for (int i=0; i<xsize; ++i)
	{
		int _x_term = (i - mu_x)*(i - mu_x) * sigma_x;
		for(int j=0; j<ysize; ++j)
		{
			int _y_term = (j - mu_y)*(j - mu_y) * sigma_y;
			for(int k=0; j<zsize; ++k)
			{
				int _z_term = (k - mu_z)*(k - mu_z) * sigma_z;
				M[i][j][k] = (1/(2*pi))*exp*(-.5*(_x_term+_y_term+_z_term));
			}
		}
	}
}

