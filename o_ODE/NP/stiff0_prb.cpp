//# include <stdlib.h>
# include <iostream.h>
//# include <iomanip.h>


#define NRANSI
//#include "o_phys\NP\nrutil.h"
extern double* vector0(int dimension);
extern int* ivector0(int dimension);
extern double** matrix0(int dimx, int dimy);
extern void free_vector0(double *vv);
extern void free_ivector0(int *vv);
extern void free_matrix0(double **aa, int dimx, int /*dimy*/);

extern void nrerror(char error_text[]);


extern void stiff0(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs0)(double, double [], double []));

extern void stifbs0(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs0)(double, double [], double []));

extern void odeint0(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs0)(double, double [], double []));

//------------------------------7-------------------------------
void STIF0s_test(int i);
void STIFF0_test();
void STIFBS0_test();
void simpr0_test();

void jacobn0_old(double x, double y[], double dfdx[], double **dfdy, int n);
void derivs0(double x, double y[], double dydx[]);


//---------------------------------
extern int kmaxO, kountO;
extern double *xpO,**ypO, dxsav;


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void STIF0s_test(int option)
{
double *y, *dydx, *yscal, *y0;
double hnext, x, htry, eps, hdid;
int n = 3;
double thickmax=50;

y     = vector0(n);
y0     = vector0(n);
dydx  = vector0(n);
yscal = vector0(n);

y[0]=1; 	y[1]=0; 	y[2]=0;
yscal[0]=1; yscal[1]=1; yscal[2]=1;
y0[0]=1;    y0[1]=1;    y0[2]=1;

eps = 1e-4;
x=0;
htry=2.9e-4;
hdid=0;
hnext=1;
int i=0;

  cout << "\n";
  cout << "i       x       htry         hnext      Y[0]    Y[1]     Y[2] \n";
  cout << "\n";



kmaxO=10000;
xpO  = vector0(kmaxO);
ypO  = matrix0(n,kmaxO);
dxsav=0;



  for(; ; )
	{
	derivs0(x,y,dydx);
	if(option==0)  stiff0 (y, dydx, n, &x, htry, eps, yscal, &hdid, &hnext, derivs0);
      else		   stifbs0(y, dydx, n, &x, htry, eps, yscal, &hdid, &hnext, derivs0);

      if(x>=thickmax) break;
      if(x+hnext >= thickmax) htry = thickmax-x;
      else 				htry = hnext;
      i++;
	cout << i <<  x  << htry  <<   hnext  << y[1]  <<  y[2]  << y[3] <<  "\n";
      }



/*
int nok, nbad;
derivs0(x,y,dydx);
odeint0(y0, n, 0, 50, eps, htry,
	0.1, &nok, &nbad, derivs0);
*/


free_vector0(y);
free_vector0(y0);
free_vector0(dydx);
free_vector0(yscal);

free_vector0(xpO);
free_matrix0(ypO,n,kmaxO);

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void jacobn0_old(double /*x*/, double y[], double dfdx[], double **dfdy, int n)
{
	int i;

	for (i=0;i<n;i++) dfdx[i]=0.0;

	dfdy[0][0] = -0.013-1000.0*y[2];
	dfdy[0][1]=0.0;
	dfdy[0][2] = -1000.0*y[0];
	dfdy[1][0]=0.0;
	dfdy[1][1] = -2500.0*y[2];
	dfdy[1][2] = -2500.0*y[1];
	dfdy[2][0] = -0.013-1000.0*y[2];
	dfdy[2][1] = -2500.0*y[2];
	dfdy[2][2] = -1000.0*y[0]-2500.0*y[1];
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void derivs0(double /*x*/, double y[], double dydx[])     // yp==dydx
{
	dydx[0] = -0.013*y[0]-1000.0*y[0]*y[2];
	dydx[1] = -2500.0*y[1]*y[2];
	dydx[2] = -0.013*y[0]-1000.0*y[0]*y[2]-2500.0*y[1]*y[2];
/*
	dydx[0] = -0.013*y[0];
	dydx[1] =  0.00025*y[1]*y[2];
	dydx[2] =  0.013*y[0]+0.00025*y[1]*y[2];
  */
}
/* (C) Copr. 1986-92 Numerical Recipes Software $=)J"!:V. */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void STIFF0_test() {STIF0s_test(0);}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void STIFBS0_test(){STIF0s_test(1);}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void simpr0_test(){}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

