#ifndef _MATRIX3D_H
#define _MATRIX3D_H
#include <vector>
#include <iostream>
#include <math.h>
using namespace std;

template<typename T>
class Matrix3d: public vector<T>
{
public:
	//constructors
	Matrix3d(size_t xsize, size_t ysize, size_t zsize) :
		vector<T>(xsize*ysize*zsize),
		_xsize(xsize), _ysize(ysize), _zsize(zsize)
		{};
	Matrix3d():
		vector<T>(0),
		_xsize(0), _ysize(0), _zsize(0)
		{};
	//() operator overloads
	T operator()(size_t i, size_t j, size_t k) const 
  { 
     return vector<T>::operator[](i*_ysize*_zsize + j*_zsize + k); 
  }

  T& operator()(size_t i, size_t j, size_t k) 
  { 
     return vector<double>::operator[](i*_ysize*_zsize + j*_zsize + k); 
  }

  void set_gaussian_interior(double mean, double var, double gridsize,
  													int xstart_ix, int ystart_ix, int zstart_ix)
  {
		double sigma_x = var;
		double sigma_y = var;
		double sigma_z = var;
		double mu_x = mean;
		double mu_y = mean;
		double mu_z = mean;
	  int xmax = _xsize;
	  int ymax = _ysize;
	  int zmax = _zsize;
		double pi   = 3.14159;
		for (int i=0; i<xmax; ++i)
		{
			double _x = (i+xstart_ix)*gridsize;
			double _x_term = (_x - mu_x)*(_x - mu_x) * sigma_x;
			for(int j=0; j<ymax; ++j)
			{	
				double _y = (j+ystart_ix)*gridsize;
				double _y_term = (_y - mu_y)*(_y - mu_y) * sigma_y;
				for(int k=0; k<zmax; ++k)
				{	
					double _z = (k+zstart_ix)*gridsize;
					double _z_term = (_z - mu_z)*(_z - mu_z) * sigma_z;
					(*this)(i,j,k) = ( 1/(pow(2*pi*sigma_x, 1.5)) 
											* exp(-.5*(_x_term+_y_term+_z_term)) );
				}
			}
	}
  }

  void set_test_interior(double mean, double var, double gridsize,
  													int xstart_ix, int ystart_ix, int zstart_ix)
  {
		double sigma_x = var;
		double sigma_y = var;
		double sigma_z = var;
		double mu_x = mean;
		double mu_y = mean;
		double mu_z = mean;
	  int xmax = _xsize;
	  int ymax = _ysize;
	  int zmax = _zsize;
		double pi   = 3.14159;
		cout<<"xmax: "<<xmax;
		for (int i=0; i<xmax; i++)
		{
			double _x = (i+xstart_ix)*gridsize;
			double _x_term = (_x - mu_x)*(_x - mu_x) * sigma_x;
			for(int j=0; j<ymax; j++)
			{	
				double _y = (j+ystart_ix)*gridsize;
				double _y_term = (_y - mu_y)*(_y - mu_y) * sigma_y;
				for(int k=0; k<zmax; k++)
				{	
					double _z = (k+zstart_ix)*gridsize;
					double _z_term = (_z - mu_z)*(_z - mu_z) * sigma_z;
					(*this)(i,j,k) = (i)*100+(j)*10+(k);
				}
			}
	}
  }

  void reset_boundaries(double boundary_value=0)
  {
  	int xmax = _xsize;
	  int ymax = _ysize;
	  int zmax = _zsize;
	  //set xup and xdn borders
		for(int j=0; j<ymax; ++j)
		{	
			for(int k=0; k<zmax; ++k)
			{	
				(*this)(0,j,k) = boundary_value;
				(*this)(xmax-1, j, k) = boundary_value;
			}
		}

		//set yup and ydn borders
		for(int i=0; i<xmax; ++i)
		{	
			for(int k=0; k<zmax; ++k)
			{	
				(*this)(i,0,k) = boundary_value;
				(*this)(i, ymax-1, k) = boundary_value;
			}
		}

		//set zup and zdn borders
		for(int i=0; i<xmax; ++i)
		{	
			for(int j=0; j<ymax; ++j)
			{	
				(*this)(i,j,0) = boundary_value;
				(*this)(i,j, zmax-1) = boundary_value;
			}
		}
  }


  void init_mean_var()
  {
  	 double sum=0; double varsum=0;
  	 double mean, var;
  	 double t;
  	 int n = 0;
  	 for (int x=1; x<_xsize-1; x++)
  		{
  			for (int y=1; y<_ysize-1; y++)
  			{
  				for (int z=1; z<_zsize-1; z++)
  				{
						sum += (*this)(x,y,z);
						n++;	
  				}
  			}
  		}
  	 mean = sum / (double) n;
  	 for (int x=1; x<_xsize-1; x++)
  		{
  			for (int y=1; y<_ysize-1; y++)
  			{
  				for (int z=1; z<_zsize-1; z++)
  				{
						t = (*this)(x,y,z);
  					varsum += t;	
  				}
  			}
  		}
 
  		var = varsum / (double) n;
  		this->_var = var;
  		this-> _mean = mean;
  }
  void calc_heat_equation(double dx, double dt, double alpha)
  {
  		double xup, xdn, yup, ydn, zup, zdn, tn0, tn1;
  		Matrix3d <T> M_new(_xsize, _ysize, _zsize);
  		double sum=0,  varsum=0;
  		int n=0;
  		double mean, var;

  		for (int x=1; x<_xsize-1; x++)
  		{
  			for (int y=1; y<_ysize-1; y++)
  			{
  				for (int z=1; z<_zsize-1; z++)
  				{

  					xup = (*this)(x+1, y, z);
  					xdn = (*this)(x-1, y, z);
  					yup = (*this)(x, y+1, z);
  					ydn = (*this)(x, y-1, z);
  					zup = (*this)(x, y, z+1);
  					zdn = (*this)(x, y, z-1);
  					tn0 = (*this)(x,y,z);
  					tn1 = tn0+ ((alpha*dt)/(dx*dx)) 
  								* (xup+xdn+yup+ydn+zup+zdn-6*tn0);
  					M_new(x,y,z) =tn1; 
  					n++;
  					sum+=tn1;
  				}
  			}
  		}

  		mean = sum / (double) n;

  	 for (int x=1; x<_xsize-1; x++)
  		{
  			for (int y=1; y<_ysize-1; y++)
  			{
  				for (int z=1; z<_zsize-1; z++)
  				{
  					tn1 = M_new(x,y,z);
  					varsum += (tn1-mean)*(tn1-mean);
						(*this)(x,y,z) = M_new(x,y,z);
  					
  				}
  			}
  		}
  		var = varsum / (double) n;
  		this->_var = var;
  		this-> _mean = mean;
  }
  size_t xsize(){return _xsize;};
  size_t ysize(){return _ysize;};
  size_t zsize(){return _zsize;};
  double _mean, _var;
private:
	size_t _xsize, _ysize, _zsize;
};


//prints cross section of matrix to screen, for debugging
void print_xslice(Matrix3d<double> &M, int xcoord)
{
	for (int y=0; y<M.ysize(); y++)
	{
		for (int z=0; z<M.zsize(); z++)
		{
			cout<<M(xcoord, y, z)<<',';
		}
		cout<<endl;
	}
}

#endif