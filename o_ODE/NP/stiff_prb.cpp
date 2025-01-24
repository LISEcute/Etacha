//# include <stdlib.h>
# include <iostream.h>
//# include <iomanip.h>

extern double* vector(int start, int dimension);
extern int* ivector(int start, int dimension);
extern double** matrix(int startx, int dimx, int starty, int dimy);
extern void free_vector(double *vv, int start, int dimension);
extern void free_ivector(int *vv, int /*start*/, int /*dimension*/);
extern void free_matrix(double **aa, int startx, int dimx, int /*starty*/, int /*dimy*/);

#define NRANSI
//#include "o_phys\NP\nrutil.h"

extern void stiff(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));

extern void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));

extern void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []));

//------------------------------7-------------------------------
void STIFs_test(int i);
void STIFF_test();
void STIFBS_test();
void simpr_test();

void jacobn(double x, double y[], double dfdx[], double **dfdy, int n);
void derivs(double x, double y[], double dydx[]);


//---------------------------------
int kmaxO, kountO;
double *xpO,**ypO, dxsav;

                                                   /// *dydx    == *yp
                                                   /// *y       == *y  
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void STIFs_test(int option)    // 0 - STIFF, 1- STIFBS
{
double *y, *dydx, *yscal;
//double  *y0;
double hnext, x, htry, eps, hdid;
int n = 3;
double thickmax=50;

y     = vector(1,n);
//y0     = vector(1,n);
dydx  = vector(1,n);
yscal = vector(1,n);

y[1]=1; 	y[2]=0; 	y[3]=0;
yscal[1]=1; yscal[2]=1; yscal[3]=1;
//y0[1]=1;    y0[2]=1;    y0[3]=1;

eps = 1e-4;
x=0;
htry=2.9e-4;
hdid=0;
hnext=1;
int i=0;

//  cout << "\n";
//  cout << "i       x       htry         hnext      Y[1]    Y[2]     Y[3] \n";
//  cout << "\n";



//kmaxO=10000;
//xpO  = vector(1,kmaxO);
//ypO  = matrix(1,n,1,kmaxO);
//dxsav=0;



  for(; ; )
	{
	derivs(x,y,dydx);

	if(option == 0) 	stiff (y, dydx, n, &x, htry, eps, yscal, &hdid, &hnext, derivs);
	else 			stifbs(y, dydx, n, &x, htry, eps, yscal, &hdid, &hnext, derivs);
      
      if(x>=thickmax) break;
      if(x+hnext >= thickmax) htry = thickmax-x;
      else 				htry = hnext;
      i++;
	cout << i <<  x  << htry  <<   hnext  << y[1]  <<  y[2]  << y[3] <<  "\n";
      }

/*
int nok, nbad;

derivs(x,y,dydx);
odeint(y0, n, 0, 50, eps, htry,
	0.1, &nok, &nbad, derivs);
*/


free_vector(y,	1,n);
//free_vector(y0,	1,n);
free_vector(dydx, 1,n);
free_vector(yscal,1,n);

//free_vector(xpO,	1,n);
//free_matrix(ypO,	1,n,1,kmaxO);

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void jacobn(double /*x*/, double y[], double dfdx[], double **dfdy, int n)
{
	int i;

	for (i=1;i<=n;i++) dfdx[i]=0.0;

	dfdy[1][1] = -0.013-1000.0*y[3];
	dfdy[1][2]=0.0;
	dfdy[1][3] = -1000.0*y[1];
	dfdy[2][1]=0.0;
	dfdy[2][2] = -2500.0*y[3];
	dfdy[2][3] = -2500.0*y[2];
	dfdy[3][1] = -0.013-1000.0*y[3];
	dfdy[3][2] = -2500.0*y[3];
	dfdy[3][3] = -1000.0*y[1]-2500.0*y[2];
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void derivs(double /*x*/, double y[], double dydx[])
{
	dydx[1] = -0.013*y[1]-1000.0*y[1]*y[3];
	dydx[2] = -2500.0*y[2]*y[3];
	dydx[3] = -0.013*y[1]-1000.0*y[1]*y[3]-2500.0*y[2]*y[3];
/*
	dydx[1] = -0.013*y[1];
	dydx[2] =  0.00025*y[2]*y[3];
	dydx[3] =  0.013*y[1]+0.00025*y[2]*y[3];
*/      
}
/* (C) Copr. 1986-92 Numerical Recipes Software $=)J"!:V. */

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void STIFF_test() {STIFs_test(0);}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void STIFBS_test(){STIFs_test(1);}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void simpr_test()
{
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

