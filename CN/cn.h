//2017-07-24 HJP
#ifndef cn_h
#define cn_h
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
	long double F;//=dt/(dx)^2/2
	long double Fc;//Dc*dt/(dx)^2/2
	long double Fd;//Dd*dt/(dx)^2/2	
	long double Dc;//diffusion coefficient for cooperator
	long double Dd;//diffusion coefficient for defector
	long double r; //multiplication factor
	long double death; //death rate
	long double birth; //baseline birth rate
	long double cost;  //cost for cooperating 

	//densities
	long double** u;//cooperator density at each site
	long double** v;//defector density at each site
	long double** z;//1-u-v; empty density at each site
	long double** ustar;//intermediate variable
	long double** vstar;//intermediate variable
	long double** uT;//density at each site; transpose matrix of u
	long double** ustarT;//intermediate variable; transpose matrix of ustar
	long double** vT;//density at each site; transpose matrix of v
	long double** vstarT;//intermediate variable; transpose matrix of vstar
	
	//calculate tridaigonal matrix
	//for tridiagonal matrix
	long double* au;
	long double* bu;
	long double* cu;
	long double* av;
	long double* bv;
	long double* cv;
	//for d_i
	long double** du;
	long double** dv;
};

void init(Eco &sys, int N, int L, long double dt, long double dx, long double Dc, long double Dd, long double r, long double death, long double birth, long double cost)
{
	sys.N = N;
	sys.L = L;
	sys.L2 = L*L;
	sys.Max = numeric_limits<long double>::max();
	sys.Min = numeric_limits<long double>::min();
	

	sys.t=0;
	sys.dt = dt;
	sys.dx = dx;
	sys.Dc = Dc;
	sys.Dd = Dd;
	sys.F = dt/dx/dx/2.;
	sys.Fc = Dc*dt/dx/dx/2.;
	sys.Fd = Dd*dt/dx/dx/2.;
	sys.r = r;
	sys.death = death;
	sys.birth = birth;
	sys.cost = cost;


	sys.au= (long double*)calloc(L, sizeof(long double));
	sys.bu= (long double*)calloc(L, sizeof(long double));
	sys.cu= (long double*)calloc(L, sizeof(long double));
	sys.av= (long double*)calloc(L, sizeof(long double));
	sys.bv= (long double*)calloc(L, sizeof(long double));
	sys.cv= (long double*)calloc(L, sizeof(long double));
	sys.u= (long double**)calloc(L, sizeof(long double*));
	sys.ustar= (long double**)calloc(L, sizeof(long double*));
	sys.uT= (long double**)calloc(L, sizeof(long double*));
	sys.ustarT= (long double**)calloc(L, sizeof(long double*));
	sys.v= (long double**)calloc(L, sizeof(long double*));
	sys.vstar= (long double**)calloc(L, sizeof(long double*));
	sys.vT= (long double**)calloc(L, sizeof(long double*));
	sys.z= (long double**)calloc(L, sizeof(long double*));
	sys.vstarT= (long double**)calloc(L, sizeof(long double*));
	sys.du= (long double**)calloc(L, sizeof(long double*));
	sys.dv= (long double**)calloc(L, sizeof(long double*));
	for(int y=0; y<L; y++)
	{
		sys.u[y] = (long double*)calloc(L, sizeof(long double));
		sys.ustar[y] = (long double*)calloc(L, sizeof(long double));
		sys.uT[y] = (long double*)calloc(L, sizeof(long double));//we treat square lattice 
		sys.ustarT[y] = (long double*)calloc(L, sizeof(long double));
		sys.du[y] = (long double*)calloc(L, sizeof(long double));
		sys.v[y] = (long double*)calloc(L, sizeof(long double));
		sys.vstar[y] = (long double*)calloc(L, sizeof(long double));
		sys.vT[y] = (long double*)calloc(L, sizeof(long double));
		sys.vstarT[y] = (long double*)calloc(L, sizeof(long double));
		sys.dv[y] = (long double*)calloc(L, sizeof(long double));
		sys.z[y] = (long double*)calloc(L, sizeof(long double));
	}
}

void init_disk(Eco &sys)
{
	int L = sys.L;
	long double mid = (long double)(L-1)/2;
	long double th = (long double)L/10;

	for(int y=0; y<L; y++)
	{
		for(int x=0; x<L; x++)
		{
			sys.u[y][x] = 0;
			sys.v[y][x] = 0;
			if( ((long double)x-mid)*((long double)x-mid) + ((long double)y-mid)*((long double)y-mid) < (long double)th*th )
			{
				sys.u[y][x] = 0.1;
				sys.v[y][x] = 0.1;
			}
			sys.z[y][x] = 1. - sys.u[y][x] - sys.v[y][x];
		}
	}
}

void init_rand(Eco &sys)
{
	int L = sys.L;

	for(int y=0; y<L; y++)
	{
		for(int x=0; x<L; x++)
		{
			sys.u[y][x] = 0.1*drnd();
			sys.v[y][x] = 0.1*drnd();
			sys.z[y][x] = 1. - sys.u[y][x] - sys.v[y][x];
		}
	}
}
 
void init_read(Eco &sys, long double t, char* fname)
{
	FILE *fp = fopen(fname, "r");
	int L = sys.L;
	long double u,v;
	if(fp == NULL)
	{
		cout << fname << " cannot open\n";
		init_disk(sys);
	}
	else
	{
		sys.t = t;
		for(int y=0; y<L; y++)
		{
			for(int x=0; x<L; x++)
			{
				fscanf(fp, "%d %d %Lg %Lg", &x, &y, &u, &v);
				sys.u[y][x] = u;
				sys.v[y][x] = v;
				sys.z[y][x] = 1.- sys.u[y][x] - sys.v[y][x];
			}
		}
	}
	fclose(fp);
}

void transpose(long double** in, long double** out, int L)//out = in.transpose
{
	for(int y=0; y<L; y++)
	{
		for(int x=0; x<L; x++)
		{
			out[x][y] = in[y][x];
		}
	}
}

 
void set_zeroflux_matrix(Eco &sys) //set tridaigonal elements {a_i}, {b_i}, {c_i} as a 1D zero-flux boundary Crank-Nicolson method
{
	int L = sys.L;
	long double Fc = sys.Fc;
	long double Fd = sys.Fd;
	long double* au = sys.au;
	long double* bu = sys.bu;
	long double* cu = sys.cu;
	long double* av = sys.av;
	long double* bv = sys.bv;
	long double* cv = sys.cv;
	
	//Neuman boundary condition: u_{-1}=u_{1} (second order correction) is used
	//for u
	for(int i=0; i<L; i++)
	{
		au[i] = -Fc;
		bu[i] = 1.+2.*Fc;
		cu[i] = -Fc;
	}
	au[0] = 0;
	au[L-1] = -2.*Fc;
	cu[0] = -2.*Fc; 
	cu[L-1] = 0;
	
	//for v
	for(int i=0; i<L; i++)
	{
		av[i] = -Fd;
		bv[i] = 1.+2.*Fd;
		cv[i] = -Fd;
	}
	av[0] = 0;
	av[L-1] = -2.*Fd;
	cv[0] = -2.*Fd; 
	cv[L-1] = 0;
}



void thomas(long double* a, long double* b, long double* c, long double* u, long double* d, int L)//thomas algorithm (solving tridiagonal matrix: calculate \vec{u} from A\vec{u}=\vec{d})
{
	//a[i]u[i-1] + b[i]u[i] + c[i]u[i+1] = d[i] => u[i] + alpha[i]u[i+1] = beta[i]
	int fp_exception = fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);
	long double* alpha = (long double*)calloc(L, sizeof(long double));
	long double* beta = (long double*)calloc(L, sizeof(long double));
	alpha[0] = c[0]/b[0]; 
	beta[0] = d[0]/b[0];
	for(int i=1; i<L; i++)
	{
		alpha[i] = c[i]/(b[i]-a[i]*alpha[i-1]);
		beta[i] = (d[i]-a[i]*beta[i-1])/(b[i]-a[i]*alpha[i-1]);
	}

	u[L-1] = beta[L-1];
	for(int i=L-2; i>=0; i-=1)
	{
		u[i] = beta[i] - alpha[i]*u[i+1];
	}

	if(fp_exception & FE_DIVBYZERO) {cout << "FE_DIVBYZERO\n"; exit(1);}
	if(fp_exception & FE_UNDERFLOW) {cout << "underflow\n"; exit(1);}
	if(fp_exception & FE_OVERFLOW) {cout << "overflow\n"; exit(1);}
	free(alpha);
	free(beta);
}

//calculating d[y][x] to get \vec{d} at a given x or y
//functions named calc_dl* are for diffusion term, reaction terms are calculated in the 
//Neuman boundary condition: u_{-1}=u_{1} (second order correction) is used
void calc_dl(long double** in, long double** out, int L, long double F)
{
	for(int y=1; y<L-1; y++)
	{
		out[y][0] = F*( in[y+1][0] + in[y-1][0] + 2.*in[y][1] - 4.*in[y][0] ) + in[y][0]; 
		for(int x=1; x<L-1; x++)
		{
			out[y][x] = F*( in[y+1][x] + in[y-1][x] + in[y][x-1] + in[y][x+1] - 4.*in[y][x] ) + in[y][x]; 
		}
		out[y][L-1] = F*( in[y+1][L-1] + in[y-1][L-1] + 2.*in[y][L-2] - 4.*in[y][L-1] ) + in[y][L-1]; 
	}
}

void calc_dl0(long double** in, long double** out, int L, long double F)
{
	out[0][0] = F*( 2.*in[1][0] + 2.*in[0][1] - 4.*in[0][0] ) + in[0][0]; 
	for(int x=1; x<L-1; x++)
	{
		out[0][x] = F*( 2.*in[1][x] + in[0][x-1] + in[0][x+1] - 4.*in[0][x] ) + in[0][x]; 
	}
	out[0][L-1] = F*( 2.*in[1][L-1] + 2.*in[0][L-2] - 4.*in[0][L-1] ) + in[0][L-1]; 
}

void calc_dlL(long double** in, long double** out, int L, long double F)
{
	out[L-1][0] = F*( 2.*in[L-2][0] + 2.*in[L-1][1] - 4.*in[L-1][0] ) + in[L-1][0]; 
	for(int x=1; x<L-1; x++)
	{
		out[L-1][x] = F*( 2*in[L-2][x] + in[L-1][x-1] + in[L-1][x+1] - 4.*in[L-1][x] ) + in[L-1][x]; 
	}
	out[L-1][L-1] = F*( 2.*in[L-2][L-1] + 2.*in[L-1][L-2] - 4.*in[L-1][L-1] ) + in[L-1][L-1]; 
}

void calc_d(Eco &sys)
{
	int L= sys.L;
	int N= sys.N;
	long double dt = sys.dt;
	long double Fc = sys.Fc;
	long double Fd = sys.Fd;
	long double r = sys.r;
	long double birth = sys.birth;
	long double death = sys.death;
	long double** u = sys.u;
	long double** v = sys.v;
	long double** z = sys.z;//1-u-v
	long double** du = sys.du;
	long double** dv = sys.dv;
	long double fc, fd;
	long double sumz, sumw, sumab, temp;
	long double uxy, vxy, zxy;


	//for l=1, ..., L-1 (except boundary)
	calc_dl(u, du, L, Fc);
	calc_dl(v, dv, L, Fd);

	//l=0
	calc_dl0(u, du, L, Fc);
	calc_dl0(v, dv, L, Fd);

	//l=L-1
	calc_dlL(u, du, L, Fc);
	calc_dlL(v, dv, L, Fd);

	int fp_exception = fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);
	//calcuating f at for all patches
	for(int x=0; x<L; x++)
	{
		for(int y=0; y<L; y++)
		{
			uxy = u[y][x];
			vxy = v[y][x];
			zxy = z[y][x];
			fc=0;fd=0;
			if(uxy>sys.Min)
			{
				sumz=0.; temp=1.;	
				for(int i=0; i<N; i++)
				{
					sumz += temp;
					temp *= zxy;
				}
				fd = r*uxy*( N - sumz )/N/(uxy+vxy);
				fc = fd - 1. - (r-1.)*pow(zxy, N-1) + r*sumz/N;
			}
			du[y][x] += dt*uxy*(zxy*(fc+birth)-death);
			dv[y][x] += dt*vxy*(zxy*(fd+birth)-death);
			if(fp_exception & FE_DIVBYZERO) {cout << "FE_DIVBYZERO\n"; exit(1);}
			if(fp_exception & FE_UNDERFLOW) {cout << "underflow\n"; exit(1);}
			if(fp_exception & FE_OVERFLOW) {cout << "overflow\n"; exit(1);}
		}
	}
}

void cn(Eco &sys)
{
	int L = sys.L;
	
	//calculate r.h.s of step (1)
	calc_d(sys);
	
	//=========== for u =========================
	//step1: {d} -> {ustar} (update x direction)
	for(int y=0; y<L; y++)
	{
		thomas(sys.au, sys.bu, sys.cu, sys.ustar[y], sys.du[y], L);
	}

	//step2: {ustar} -> {u^(n+1)} (update y direction)
	transpose(sys.ustar, sys.ustarT, L);
	for(int x=0; x<L; x++)
	{
		thomas(sys.au, sys.bu, sys.cu, sys.uT[x], sys.ustarT[x], L);
	}
	transpose(sys.uT, sys.u, L);
	
	//=========== for v =========================
	//step1: {d} -> {vstar} 
	for(int y=0; y<L; y++)
	{
		thomas(sys.av, sys.bv, sys.cv, sys.vstar[y], sys.dv[y], L);
	}

	//step2: {vstar} -> {v^(n+1)}
	transpose(sys.vstar, sys.vstarT, L);
	for(int x=0; x<L; x++)
	{
		thomas(sys.av, sys.bv, sys.cv, sys.vT[x], sys.vstarT[x], L);
	}
	transpose(sys.vT, sys.v, L);

	//update z = 1-u-v
	for(int x=0; x<L; x++)
	{
		for(int y=0; y<L; y++)
		{
			sys.z[y][x] = 1.-sys.u[y][x]-sys.v[y][x];
		}
	}
	sys.t += sys.dt;
}





void write(Eco &sys)
{
    char fname[200];
	sprintf(fname, "data/r%Lf_D%Lf_L%d_t%Lf.d", sys.r, sys.Dd/sys.Dc , sys.L, sys.t);
	
    //write the data : x, y, u(x, y), v(x, y)
	FILE *ofp = fopen(fname, "w");
 
	int L = sys.L;
    for(int y=0; y<L; y++)
    {
	    for(int x=0; x<L; x++)
		{
			fprintf(ofp, "%d %d %Lf %Lf\n", x, y, sys.u[y][x], sys.v[y][x]);
		}
    }
    fclose(ofp);
}



void show(Eco &sys)
{
	int L = sys.L;
	cout << "============== u ==============" << endl;
	for (int y=0; y<L; y++)
	{
		for(int x=0; x<L; x++)
		{
			cout << sys.u[y][x] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << "============== v ==============" << endl;
	for (int y=0; y<L; y++)
	{
		for(int x=0; x<L; x++)
		{
			cout << sys.v[y][x] << " ";
		}
		cout << endl;
	}
}

#endif
