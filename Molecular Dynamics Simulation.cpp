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
#include <string>
#include <fstream>

using namespace std;

int grid_dim = 10;
int npart = grid_dim*grid_dim*grid_dim;
float temp = 1;
float dt = 0.001;
double en;


double *x,*y,*z,*vx,*vy,*vz,*xm,*ym,*zm,*fx,*fy,*fz;

double randf();
void init();
void force();
void integrate();
void lattice_pos(int i);

int main() {
	srand(time(NULL));

	x = new double[npart];
	y = new double[npart];
	z = new double[npart];
	vx = new double[npart];
	vy = new double[npart];
	vz = new double[npart];
	xm = new double[npart];
	ym = new double[npart];
	zm = new double[npart];
	fx = new double[npart];
	fy = new double[npart];
	fz = new double[npart];



	//int t;
	init(); 				//Initialize the position and velocities of particles.
		cout<< x[1]<<" "<<x[100]<<"\n";
		cout<< y[1]<<" "<<y[100]<<"\n";
		cout<< vx[1]<<" "<<vx[100]<<"\n";
		cout<< vy[1]<<" "<<vy[100]<<"\n";
		cout<< xm[1]<<" "<<xm[100]<<"\n";
		cout<< ym[1]<<" "<<ym[100]<<"\n";
		cout<< zm[1]<<" "<<zm[100]<<"\n";
//		cout<<pow(2,3);

//	t = 0;
//	while(t<tmax)
//	{
		force();
		integrate();
		string filename = "test.xyz";
		ofstream file(filename);
		file << npart<<endl;
		file <<"Filename "<< filename<<endl;

		for(int index = 0; index<npart; index++)
		{
			file<<"Atom"<<index<<"    "<<x[index]<<"    "<<y[index]<<"    "<<z[index]<<endl;
		}
		file.close();



//		t=t+delt;
//		sample();
//	}
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] vx;
	delete[] vy;
	delete[] vz;
	delete[] xm;
	delete[] ym;
	delete[] zm;
	delete[] fx;
	delete[] fy;
	delete[] fz;
	return 0;
}

void init()
{
	double sumvx=0;
	double sumvy=0;
	double sumvz=0;
	double sumvx2=0;
	double sumvy2=0;
	double sumvz2=0;
	for(int i=0;i<npart;i++)
	{

		lattice_pos(i);						//index for particles

		vx[i] = ( randf() - 0.5);
		vy[i] = ( randf() - 0.5);
		vz[i] = ( randf() - 0.5);			//randomly initializing velocity between [-0.5,0.5]

		sumvx += vx[i];
		sumvy += vy[i];
		sumvz += vz[i];

		sumvx2 += vx[i] * vx[i];
		sumvy2 += vy[i] * vy[i];
		sumvz2 += vz[i] * vz[i];
	}

	sumvx /= npart;
	sumvy /= npart;
	sumvz /= npart;

	sumvx2 /= npart;
	sumvy2 /= npart;
	sumvz2 /= npart;

	double fsx= (3*temp/sumvx2) * 0.5;
	double fsy= (3*temp/sumvy2) * 0.5;
	double fsz= (3*temp/sumvz2) * 0.5;//scaling factor

	for(int i=0; i < npart; i++)
	{
		vx[i]=(vx[i]-sumvx)*fsx;
		vy[i]=(vy[i]-sumvy)*fsy;
		vz[i]=(vz[i]-sumvz)*fsz;//To make velocity of the overall system as 0

		xm[i]=x[i]-vx[i]*dt;
		ym[i]=y[i]-vy[i]*dt;
		zm[i]=z[i]-vz[i]*dt;
	}

	double sum=0;
	for(int i =0; i<npart; i++)
	{
		sum+=vx[i];
	}
	cout<<sum<<endl;
}

double randf()
{
	double num;
	num =  rand();
	num /= RAND_MAX;
	return num;
}

void force()
{
	en = 0;

	int rc = 2;						//cut off distance
	int box = npart;
	double ecut = 4 * ((1/pow(rc,12)) - (1/pow(rc,12)));

	for(int i=0; i<npart; i++)
	{
		fx[i] = 0;
		fy[i] = 0;
		fz[i] = 0;
	}
	for(int i = 0; i < npart-1; i++)
		for(int j = i+1; j < npart; j++)
		{
			double xr = x[i] - x[j];
			double yr = y[i] - y[j];
			double zr = z[i] - z[j];

			xr -= box*round(xr/box);
			yr -= box*round(yr/box);
			zr -= box*round(zr/box);

			double r2 = pow(xr, 2) + pow(yr, 2) + pow(zr, 2);

			if(r2<(rc*rc))					//rc2 cutoff = 4
			{
				double r2i = 1/r2;
				double r6i = r2i* r2i* r2i;
				double ff = 48*r2i*r6i*(r6i-0.5);

				fx[i] += ff*xr;
				fx[j] -= ff*xr;

				fy[i] += ff*yr;
				fy[j] -= ff*yr;

				fz[i] += ff*zr;
				fz[j] -= ff*zr;

				en += 4* r6i* (r6i - 1) - ecut;
			}
		}
}



void integrate()
{
	double sumvx=0;
	double sumvy=0;
	double sumvz=0;

	double sumv2=0;

	for(int i = 0; i<npart; i++)
	{
		double xx = 2*x[i] -xm[i] + dt*dt*fx[i];
		double yy = 2*y[i] -ym[i] + dt*dt*fy[i];
		double zz = 2*z[i] -zm[i] + dt*dt*fz[i];

		double vix=(xx-xm[i])/(2*dt);
		double viy=(yy-ym[i])/(2*dt);
		double viz=(zz-zm[i])/(2*dt);

		sumvx += vix;
		sumvy += viy;
		sumvz += viz;

		sumv2 += pow(vix,2) + pow(viy,2) + pow(viz,2);

		xm[i]=x[i];
		x[i]=xx;

		ym[i]=y[i];
		y[i]=yy;

		zm[i]=z[i];
		z[i]=zz;
	}
	cout<<"sumvx "<<sumvx<<endl;
	cout<<"sumvy "<<sumvy<<endl;
	cout<<"sumvz "<<sumvz<<endl;

	temp=sumv2/(3*npart);							//doubt should temp have 3 values

	double etot=(en+0.5*sumv2)/npart;

	cout<<"temp :"<<temp<<endl;
	cout<<"etot :"<<etot<<endl;
}


void lattice_pos(int i)
{
	x[i] = i%10 ;
	y[i] = i%100;
	y[i] = floor(y[i]/10);
	z[i] = i%1000;
	z[i] = floor(z[i]/100);
}

