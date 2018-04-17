//////////////////////////////////////////////////////
//  solve 2D reaction diffusion equation using CN and D'yakonov method
//	with zero-flux boundary
//	arguement: r (multiplication factor), D (ratio of diffusion coefficients), (optional: starting from the existing configuation at a given time t)
//////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <fenv.h>
using namespace std;

#include "cn.h"

int main(int argc, char** argv)
{
	if(argc!=3 && argc!=4)
	{
		cout << "Usage: " << argv[0] << " r D existing-t(optional)\n";
		exit(1);
	}
	double dt=0.1;
	double dx=1.4;
	int L = 283;
	int N = 8;
	double Dc = 1.;
	double r = atof(argv[1]);
	double Dd = atof(argv[2]);
	double death = 1.2; 
	double birth = 1;
	double cost = 1;
	double tmax = 10000;
	int Nt = tmax/dt+1;
	int record = 100./dt+0.5;

	Eco sys;
	init(sys, N, L, dt, dx, Dc, Dd, r, death, birth, cost);
	set_zeroflux_matrix(sys);//matrix form for updating x axis and y asix are the same in our case.
	init_disk(sys);//or you can us "init_rand(sys)" for randomized initial condition. 
	//if there is an existing configuration at time t
	if(argc==4)
	{
		double t = atof(argv[3]);
		char fname[200];
		sprintf(fname, "data/r%f_D%f_L%d_t%f.d", r, Dd/Dc, L, t);
		init_read(sys, t, fname);
	}

	for(int i=0; i<Nt; i++)
	{
		if(i%record==0)
		{
			write(sys);
		}
		cn(sys);
	}


	return 0;
}
