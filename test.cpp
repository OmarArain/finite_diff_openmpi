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