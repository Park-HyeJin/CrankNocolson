#ifndef diffusion_feuler_h
#define diffusion_feuler_h
#include"twist.h"
#include <iostream>
#include <fenv.h>
#include <limits>
using namespace std;
 
struct Eco
{
	//setting
    int N;//maximum group size
    int L;//linear system size, we will treat L*L 2D square lattice
	int L2;//L*L

	//for checking the numerical underflow or overflow
	long double Max;
	long double Min;
	
	//parameters
	long double t;//time
	long double dt;//time increasement
	long double dx;//space increasment
	long double sigma;//sensitivity to density functions
	long double max;//maximum of the density dependent function
	long double Dc;//diffusion coefficient for cooperator
	long double r; //multiplication factor
	long double death; //death rate
	long double birth; //baseline birth rate
	long double cost;  //cost for cooperating 
	
	//densities; one dimensional array u[i] = u[y*L+x]
	long double* u;
	long double* v;
	long double* z;
	long double* du1;
	long double* dv1;
};

void init(Eco &sys, int N, int L, long double death, long double birth, long double cost, long double Dc, long double sigma, long double max, long double r, long double dt, long double dx)
{
	sys.max = max;
    sys.N = N;
	sys.L = L;
	sys.L2 = sys.L*sys.L;
    sys.death = death;
    sys.birth = birth;
    sys.cost = cost;
    sys.Dc = Dc;
    sys.sigma = sigma;
	sys.r = r;
	sys.dt = dt; 
	sys.dx = dx; 
	sys.t = 0;
	sys.Max = numeric_limits<long double>::max();
	sys.Min = numeric_limits<long double>::min();

	sys.u = (long double*)calloc(sys.L2, sizeof(long double));
	sys.v = (long double*)calloc(sys.L2, sizeof(long double));
	sys.z = (long double*)calloc(sys.L2, sizeof(long double));
	sys.du1 = (long double*)calloc(sys.L2, sizeof(long double));
	sys.dv1 = (long double*)calloc(sys.L2, sizeof(long double));
}

void init_disk(Eco &sys)
{
	int L = sys.L;
	long double mid = (long double)(L-1)/2;
	long double th = (long double)L/10;

	int i;
	for(int y=0; y<L; y++)
	{
		for(int x=0; x<L; x++)
		{
			i = y*L+x;
			sys.u[i] = 0;
			sys.v[i] = 0;
			if( ((double)x-mid)*((double)x-mid) + ((double)y-mid)*((double)y-mid) < (double)th*th )
			{
				sys.u[i] = 0.1;
				sys.v[i] = 0.1;
			}
			sys.z[i] = 1- sys.u[i] - sys.v[i];
		}
	}
}
 

void init_read(Eco &sys, long double t, char* fname)
{
	FILE *ofp = fopen(fname, "r");
	int x,y,i;
	long double u,v;

	if(ofp == NULL)
	{
		cout << fname << " cannot open\n";
		init_disk(sys);
	}
	else
	{
		sys.t = t;
		for(int l=0; l<sys.L*sys.L; l++)
		{
			fscanf(ofp, "%d %d %Lg %Lg", &x, &y, &u, &v);
			i = y*sys.L+x;
			sys.u[i] = u;
			sys.v[i] = v;
			sys.z[i] = 1. - sys.u[i] - sys.v[i];
		}
	}
	fclose(ofp);
}


//laplacian term
long double lap(Eco &sys, int x, int y, long double* arr)
{
	//Nueman boundary condition: u[-1] = u[0] (first-order) is used
	int L = sys.L;
	int i = y*L + x;
	long double sum = 0;
	//right
	if( x != L-1 )
	{
		sum += (arr[i+1] - arr[i]);  
	}
	//left
	long double temp;
	if( x != 0 )
	{
		sum += (arr[i-1] - arr[i]);  
	}
	//up
	if( y != L-1 )
	{
		sum += (arr[i+L] - arr[i]);  
	}
	//down
	if( y != 0 )
	{
		sum += (arr[i-L] - arr[i]);  
	}

	return sum/(sys.dx)/(sys.dx);
}

//finding index = y*L+x at a given x and y
int pos(int L, int x, int y)
{
	int i;
	if(x==-1) 
	{
		x = 0;
	}
	else if(x==L)
	{
		x = L-1;
	}
	if(y==-1)
	{
		y = 0;
	}
	else if(y==L)
	{
		y = L-1;
	}
	i = y*L + x;
	return i;
}

//Diffusion coefficient function for defectors
long double Dd_ftn(Eco &sys, long double sigma, long double max, int x, int y)
{
	long double* u = sys.u;
	long double* v = sys.v;
	int L = sys.L;
	int i = pos(L, x, y);
	long double rho = u[i] + v[i];
	if(rho)
	{
		long double f = u[i]/rho;
		return sys.Dc*(1.+ sigma*(1.-rho)*rho*f*(1.-f)/max);
	}
	return sys.Dc;
}


//calculate diffusion term: D laplacian(u) + gradienet(D) product gradient(u)
//for the second term, I used average the left and right derivatives for mass conservation and isotropic symmetry
long double NonlinearDiffusion(Eco &sys, long double sigma, int x, int y, long double* arr)
{
	long double Dc = sys.Dc;
	long double temp=0.;
	long double max = sys.max;
	int L = sys.L;
	int fp_exception = fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);

	temp = Dd_ftn(sys, sigma, max, x, y)*lap(sys, x, y, arr);
	//if (x==141 and y==141) 	cout << x<< " " << y << " " << Dd_ftn(sys, sigma, max, x, y) << " " << lap(sys, x, y, arr) <<  endl;
	temp += (Dd_ftn(sys, sigma, max, x+1, y)-Dd_ftn(sys, sigma, max, x, y))*(arr[pos(L, x+1, y)]-arr[pos(L, x, y)])/sys.dx/sys.dx/2;
	temp += (Dd_ftn(sys, sigma, max, x, y)-Dd_ftn(sys, sigma, max, x-1, y))*(arr[pos(L, x, y)]-arr[pos(L, x-1, y)])/sys.dx/sys.dx/2;
	temp += (Dd_ftn(sys, sigma, max, x, y+1)-Dd_ftn(sys, sigma, max, x, y))*(arr[pos(L, x, y+1)]-arr[pos(L, x, y)])/sys.dx/sys.dx/2;
	temp += (Dd_ftn(sys, sigma, max, x, y)-Dd_ftn(sys, sigma, max, x, y-1))*(arr[pos(L, x, y)]-arr[pos(L, x, y-1)])/sys.dx/sys.dx/2;
	if(fp_exception & FE_DIVBYZERO) {cout << "FE_DIVBYZERO\n"; exit(1);}
	if(fp_exception & FE_UNDERFLOW) {cout << "underflow\n"; exit(1);}
	if(fp_exception & FE_OVERFLOW) {cout << "overflow\n"; exit(1);}
	return temp;
}

//calculating time derivative: diffusion + reaction
void v(Eco &sys, int x, int y, long double* iu, long double* iv, long double* du, long double* dv)
{
	long double u, v, z;
	long double sumz, temp;
	long double fc, fd;
	long double Dc = sys.Dc;
	long double sigma = sys.sigma;
	int L = sys.L, L2 = sys.L2, N = sys.N, i;
	long double birth = sys.birth;
	long double death = sys.death;
	long double r = sys.r;
	int fp_exception = fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);

	i = y*L+x;
	u = iu[i];
	v = iv[i];
	z = 1.-u-v;
	fc=0;fd=0;
	if (u>sys.Min)
	{
		sumz=0., temp=1.;
		for(int l=0; l<N; l++)
		{
			sumz += temp;
			temp *= z;
		}
		fd = r*u*( N - sumz )/N/(u+v);
		fc = fd - 1. - (r-1.)*pow(z, N-1) + r*sumz/N;
		if(fp_exception & FE_DIVBYZERO) {cout << "FE_DIVBYZERO\n"; exit(1);}
		if(fp_exception & FE_UNDERFLOW) {cout << "underflow\n"; exit(1);}
		if(fp_exception & FE_OVERFLOW) {cout << "overflow\n"; exit(1);}
	}
	//cout << sys.t << " " << i << " x=" << x << " y=" << y <<  " fd=" << fd << " fc=" << fc << endl;
	du[i] = Dc*lap(sys, x, y, sys.u) + u*(z*(fc+birth)-death);
	dv[i] = NonlinearDiffusion(sys, sigma, x, y, sys.v) + v*(z*(fd+birth)-death);
	if(fp_exception & FE_DIVBYZERO) {cout << "FE_DIVBYZERO\n"; exit(1);}
	if(fp_exception & FE_UNDERFLOW) {cout << "underflow\n"; exit(1);}
	if(fp_exception & FE_OVERFLOW) {cout << "overflow\n"; exit(1);}
}


//time evolution of densities. Forward Euler method is used.
void calc(Eco &sys)
{
	long double* du1 = sys.du1;
	long double* dv1 = sys.dv1;
	long double dt = sys.dt;
	long double x, y;
	int L = sys.L, L2 = sys.L2, N = sys.N;
	int fp_exception = fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);
	int i;

	//calculate f'(t)
	for(int y=0; y<L; y++)
	{
		for(int x=0; x<L; x++)
		{
			i = y*L + x;
			v(sys, x, y, sys.u, sys.v, du1, dv1);
		}
	}

	//update
	for(int i=0; i<L2; i++)
	{
		sys.u[i] += dt*du1[i];
		sys.v[i] += dt*dv1[i];
		sys.z[i] = 1. - sys.u[i] - sys.v[i];
	}
	if(fp_exception & FE_DIVBYZERO) {cout << "FE_DIVBYZERO\n"; exit(1);}
	if(fp_exception & FE_UNDERFLOW) {cout << "underflow\n"; exit(1);}
	if(fp_exception & FE_OVERFLOW) {cout << "overflow\n"; exit(1);}

	sys.t += sys.dt;
}

void write(Eco &sys)
{
    char fname[200];
	sprintf(fname, "data/r%Lf_sigma%Lf_L%d_t%Lf.d", sys.r, sys.sigma,  sys.L, sys.t);
	
    //write the data : x y u(x,y) v(x,y) 
	FILE *ofp = fopen(fname, "w");
 
	int L = sys.L;
	int i;
    for(int y=0; y<L; y++)
    {
	    for(int x=0; x<L; x++)
		{
			i = y*L + x;
			fprintf(ofp, "%d %d %Lf %Lf\n", x, y, sys.u[i], sys.v[i]);
		}
    }
    fclose(ofp);
}


#endif
 
