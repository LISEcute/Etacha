#include <math.h>

#define NRANSI
//#include "nrutil.h"
#define EPS 1.0e-4

extern double* vector0(int dimension);
extern void free_vector0(double *vv);
//------------------------------------------------------------
// Forward Jacobian
//  df_i/dx =0

void fdjac(int n, double y[], double yp[], double **dfdy,
	void (*vecfunc)(int , double [], double []))
{
	int i,j;
	double h,temp,*yp1;

	yp1=vector0(n);

	for (j=0;j<n;j++)
      	{
		temp=y[j];
		h=EPS*fabs(temp);
		if (h == 0.0) h = EPS;
		y[j]=temp+h;
		h=y[j]-temp;
		(*vecfunc)(n,y,yp1);
		y[j]=temp;

		for (i=0;i<n;i++)
            		dfdy[i][j]=(yp1[i]-yp[i])/h;
		}

	free_vector0(yp1);
}

//------------------------------------------------------------

void FwdJacobian(int n, double x, double y[], double yp[], double **dfdy,
	void (*derivs)(double x, double [], double []))
{
	int i,j;
	double h,temp,*yp1;

	yp1=vector0(n);

	for (j=0;j<n;j++)
      	{
		temp=y[j];
		h=EPS*fabs(temp);
		if (h == 0.0) h = EPS;
		y[j]=temp+h;
		h=y[j]-temp;
		(*derivs)(x,y,yp1);
		y[j]=temp;

		for (i=0;i<n;i++)
            		dfdy[i][j]=(yp1[i]-yp[i])/h;
		}

	free_vector0(yp1);
}

#undef EPS
#undef NRANSI
