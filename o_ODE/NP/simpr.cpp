#define NRANSI
//#include "o_phys\NP\nrutil.h"
#define max(a,b)  (a>b ? a : b)
#define min(a,b)  (a<b ? a : b)
#define SQR(a)   (a == 0.0 ? 0.0 : a*a)

extern void lubksb(double **a, int n, int *indx, double b[]);
extern void ludcmp(double **a, int n, int *indx, double *d);

extern double* vector(int start, int dimension);
extern int* ivector(int start, int dimension);
extern double** matrix(int startx, int dimx, int starty, int dimy);
extern void free_vector(double *vv, int start, int dimension);
extern void free_ivector(int *vv, int /*start*/, int /*dimension*/);
extern void free_matrix(double **aa, int startx, int dimx, int /*starty*/, int /*dimy*/);

extern void nrerror(char error_text[]);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[],
	void (*derivs)(double, double [], double []))
{
	int i,j,nn,*indx;
	double d,h,x,**a,*del,*ytemp;

	indx=ivector(1,n);
	a=matrix(1,n,1,n);
	del=vector(1,n);
	ytemp=vector(1,n);
	h=htot/nstep;
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) a[i][j] = -h*dfdy[i][j];
		++a[i][i];
	}
	ludcmp(a,n,indx,&d);
	for (i=1;i<=n;i++)
		yout[i]=h*(dydx[i]+h*dfdx[i]);
	lubksb(a,n,indx,yout);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+(del[i]=yout[i]);
	x=xs+h;
	(*derivs)(x,ytemp,yout);
	for (nn=2;nn<=nstep;nn++) {
		for (i=1;i<=n;i++)
			yout[i]=h*yout[i]-del[i];
		lubksb(a,n,indx,yout);
		for (i=1;i<=n;i++)
			ytemp[i] += (del[i] += 2.0*yout[i]);
		x += h;
		(*derivs)(x,ytemp,yout);
	}
	for (i=1;i<=n;i++)
		yout[i]=h*yout[i]-del[i];
	lubksb(a,n,indx,yout);
	for (i=1;i<=n;i++)
		yout[i] += ytemp[i];
	free_vector(ytemp,1,n);
	free_vector(del,1,n);
	free_matrix(a,1,n,1,n);
	free_ivector(indx,1,n);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software $=)J"!:V. */
