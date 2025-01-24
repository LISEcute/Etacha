#define NRANSI
//#include "o_phys\NP\nrutil.h"
#define max(a,b)  (a>b ? a : b)
#define min(a,b)  (a<b ? a : b)
#define SQR(a)   (a == 0.0 ? 0.0 : a*a)

extern void lubksb0(double **a, int n, int *indx, double b[]);
extern void ludcmp0(double **a, int n, int *indx, double *d);

extern double* vector0(int dimension);
extern int* ivector0(int dimension);
extern double** matrix0(int dimx, int dimy);
extern void free_vector0(double *vv);
extern void free_ivector0(int *vv);
extern void free_matrix0(double **aa, int dimx, int /*dimy*/);

extern void nrerror(char error_text[]);

void simpr0(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[],
	void (*derivs)(double, double [], double []));

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void simpr0(double y[], double dydx[], double dfdx[], double **dfdy, int n,                   // page 747
	double xs, double htot, int nstep, double yout[],
	void (*derivs)(double, double [], double []))
{
	int i,j,nn,*indx;
	double d,h,x,**a,*del,*ytemp;

	a    =matrix0(n,n);
	indx =ivector0(n);
	del  =vector0(n);
	ytemp=vector0(n);
	h=htot/nstep;

	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) a[i][j] = -h*dfdy[i][j];      // !!!
		++a[i][i];
		}

	ludcmp0(a,n,indx,&d);
	for (i=0;i<n;i++) yout[i]=h*(dydx[i]+h*dfdx[i]);

	lubksb0(a,n,indx,yout);
	for (i=0;i<n;i++) ytemp[i]=y[i]+(del[i]=yout[i]);

	x=xs+h;
	(*derivs)(x,ytemp,yout);

	for (nn=2;nn<=nstep;nn++)
      	{
		for (i=0;i<n;i++)  yout[i]=h*yout[i]-del[i];

		lubksb0(a,n,indx,yout);
		for (i=0;i<n;i++) ytemp[i] += (del[i] += 2.0*yout[i]);

		x += h;
		(*derivs)(x,ytemp,yout);
		}

	for (i=0;i<n;i++) yout[i]=h*yout[i]-del[i];

	lubksb0(a,n,indx,yout);
	for (i=0;i<n;i++) yout[i] += ytemp[i];

	free_vector0(ytemp);
	free_vector0(del);
	free_matrix0(a,n,n);
	free_ivector0(indx);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software $=)J"!:V. */
