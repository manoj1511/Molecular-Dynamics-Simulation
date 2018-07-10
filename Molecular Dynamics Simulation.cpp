//============================================================================
// Name        : Molecular.cpp
// Author      : Manoj
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

int grid_dim = 10;
int npart = grid_dim*grid_dim*grid_dim;
float temp = 1;
float dt = 0.001;
double en;

double randf();
void init(int x[], double v[], double xm[]);
void force(double f[], int x[]);
void integrate(double f[], int x[], double xm[]);

int main() {
	srand(time(NULL));
	double *v,*xm,*f;
	int *x;

	x = (int*) malloc (npart*sizeof(int));
	v = (double*) malloc (npart*sizeof(double));
	xm = (double*) malloc (npart*sizeof(double));
	f = (double*) malloc (npart*sizeof(double));


	//int t;
	init(x,v,xm); 				//Initialize the position and velocities of particles.
//		cout<< x[1]<<" "<<x[100]<<"\n";
//		cout<< v[1]<<" "<<v[100]<<"\n";
//		cout<< xm[1]<<" "<<xm[100]<<"\n";
//		cout<<pow(2,3);

	//	t = 0;
//	while(t<tmax)
//	{
		force(f,x);
		integrate(f,x,xm);
//		t=t+delt;
//		sample();
//	}
	free(x);
	free(v);
	free(xm);
	free(f);
	return 0;
}

void init(int x[], double v[], double xm[])
{
	double sumv=0;
	double sumv2=0;
	for(int i=1;i<=npart;i++)
	{
		x[i] = i;						//index for particles									///ask if this is ok.
		v[i] = ( randf() - 0.5);		//randomly initializing velocity between [-0.5,0.5] 	////Ask Y velocity scaled from -0.5 to 0.5
		sumv = sumv + v[i];
		sumv2 = sumv2 + v[i] * v[i] ;

	}

	sumv = sumv/npart;
	sumv2 = sumv2/npart;
	double fs= (3*temp/sumv2) * 0.5;	//scaling factor
	for(int i=1; i <= npart; i++)
	{
		v[i]=(v[i]-sumv)*fs;			//To make velocity of the overall system as 0
		xm[i]=x[i]-v[i]*dt;				//previous     											////ask to explain output of this
	}
	double sum = 0;
	for(int i=1; i <= npart; i++)
		{
			sum+=v[i];

		}
	cout<< sum <<"\n";

}

double randf()
{
	double num;
	num =  rand();
	num /= RAND_MAX;
	return num;
}

void force(double f[], int x[])
{
	en = 0;
	int rc = 2;						//cut off distance
	int box = npart;
	double ecut = 4 * ((1/pow(rc,12)) - (1/pow(rc,12)));

	for(int i=1; i<=npart; i++)
		f[i] = 0;
	for(int i = 1; i<npart; i++)
		for(int j=i+1; j<=npart; j++ )
		{
			double xr = x[i] - x[j];
			xr = xr - box*round(xr/box);
			double r2 = xr * xr;
			if(r2<(rc*rc))					//rc2 cutoff = 4
			{
				double r2i = 1/r2;
				double r6i = r2i* r2i* r2i;
				double ff = 48*r2i*r6i*(r6i-0.5);
				f[i] = f[i] + ff*xr;
				f[j] = f[j] - ff*xr;
				en = en + 4* r6i* (r6i - 1) - ecut;
			}
		}
}

void integrate(double f[], int x[], double xm[])
{
	double sumv = 0;
	double sumv2 = 0;
	for(int i = 1; i<=npart; i++)
	{
		double xx = 2*x[i] -xm[i] + dt*dt*f[i];
		double vi=(xx-xm[i])/(2*dt);
		sumv=sumv+vi;
		sumv2=sumv2+pow(vi,2);
		xm[i]=x[i];
		x[i]=xx;
	}
	cout<<"sumv "<<sumv<<endl;
	temp=sumv2/(3*npart);
	double etot=(en+0.5*sumv2)/npart;
	cout<<"temp :"<<temp<<endl;
	cout<<"etot :"<<etot<<endl;
}


