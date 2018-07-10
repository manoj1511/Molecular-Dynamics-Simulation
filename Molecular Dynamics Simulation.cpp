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

using namespace std;

int grid_dim = 10;
int npart = grid_dim*grid_dim*grid_dim;
int temp = 1;
float dt = 0.001;

double randf();
void init(int x[], double v[], double xm[]);

int main() {
	double *v,*xm;
	int *x;

	x = (int*) malloc (npart*sizeof(int));
	v = (double*) malloc (npart*sizeof(double));
	xm = (double*) malloc (npart*sizeof(double));


	//int t;
	init(x,v,xm); 				//Initialize the position and velocities of particles.
		cout<< x[1]<<" "<<x[100]<<"\n";
		cout<< v[1]<<" "<<v[100]<<"\n";
		cout<< xm[1]<<" "<<xm[100]<<"\n";


	//	t = 0;
//	while(t<tmax)
//	{
//		force(f,en);
//		integrate(f,en);
//		t=t+delt;
//		sample();
//	}
	free(x);
	free(v);
	free(xm);
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

	cout<<x[1]<<" "<<v[1]<<endl;
	sumv = sumv/npart;
	sumv2 = sumv2/npart;
	double fs= (3*temp/sumv2) * 0.5;	//scaling factor
	cout<< "Fs" << fs<<endl;
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
	num =  (rand() % 10) + 1;
	num /= 10;
	return num;
}
