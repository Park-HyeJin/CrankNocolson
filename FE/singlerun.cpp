#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "diffusion.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <fenv.h>
using namespace std;

int main (int argc, char** argv)
{
	if(argc!=4 && argc!=5)
	{
		cout << "Usage: " << argv[0] << " r, sigma, max (maximum of the density dependent function), (optioanl)exsting-t\n";
		exit(1);
	}
	Eco sys;
	int N = 8; int L = 283;
	long double d = 1.2; long double b = 1.; long double c = 1.; 
	long double r = atof(argv[1]);
	long double Dc = 1.; long double Dd; long double sigma = atof(argv[2]);
	long double t=0;
	long double dt=0.005; long double dx=1.4;
	long double tmax = 100;
	long double max = atof(argv[3]);

	int record = (int)(10./dt+0.5);
	int Nt = tmax/dt+1;
	char str[100];
	char fname[200];

	//check the Von Neumann stability analysis
	if( dt > dx*dx/(1.+sigma)/4. )
	{
		cout << "dt has to be more smaller\n";
		exit(1);
	}

	//initialize 
	init(sys, N, L, d, b, c, Dc, sigma, max, r, dt, dx);

	//initialize the initial distribution of x and y
	if(argc==4)
	{
		init_disk(sys);
	}
	else if(argc==5)
	{
		t = atof(argv[4]);
		sprintf(fname, "data/r%Lf_sigma%Lf_L%d_t%Lf.d", sys.r, sys.sigma, sys.L, t);
		init_read(sys, t, fname);
		Nt = (int)(tmax-sys.t)/dt+1;
	}
	
	//run
	for(int idx=0; idx<Nt; idx++)
	{
		if(idx%record==0)
		{
			write(sys);
		}
		calc(sys);
	}


	return 0;
}
