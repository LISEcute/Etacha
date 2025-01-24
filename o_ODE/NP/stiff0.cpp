#include <math.h>
#define NRANSI
//#include "o_phys\NP\nrutil.h"
#define max(a,b)  (a>b ? a : b)
#define min(a,b)  (a<b ? a : b)


#define SAFETY 0.9
#define GROW 1.5
#define PGROW -0.25
#define SHRNK 0.5
#define PSHRNK (-1.0/3.0)
#define ERRCON 0.1296
#define MAXTRY 40
#define GAM (1.0/2.0)
#define A21 2.0
#define A31 (48.0/25.0)
#define A32 (6.0/25.0)
#define C21 -8.0
#define C31 (372.0/25.0)
#define C32 (12.0/5.0)
#define C41 (-112.0/125.0)
#define C42 (-54.0/125.0)
#define C43 (-2.0/5.0)
#define B1 (19.0/9.0)
#define B2 (1.0/2.0)
#define B3 (25.0/108.0)
#define B4 (125.0/108.0)
#define E1 (17.0/54.0)
#define E2 (7.0/36.0)
#define E3 0.0
#define E4 (125.0/108.0)
#define C1X (1.0/2.0)
#define C2X (-3.0/2.0)
#define C3X (121.0/50.0)
#define C4X (29.0/250.0)
#define A2X 1.0
#define A3X (3.0/5.0)

extern double* vector0(int dimension);
extern int* ivector0(int dimension);
extern double** matrix0(int dimx, int dimy);
extern void free_vector0(double *vv);
extern void free_ivector0(int *vv);
extern void free_matrix0(double **aa, int dimx, int /*dimy*/);

extern void nrerror(char error_text[]);

//extern void jacobn0(double x, double y[], double dfdx[], double **dfdy, int n);
extern void lubksb0(double **a, int n, int *indx, double b[]);
extern void ludcmp0(double **a, int n, int *indx, double *d);

void stiff0(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs0)(double, double [], double []));

void stiff0_OT(double *y, double *dydx, int n, double *x, double htry, double eps,
	double *yscal, double *hdid, double *hnext, double **dfdy,
	void (*derivs0)(double x, double *y, double *yp),
      void (*jacobnL)(double x, double *y, double *dfdx, double **dfdy, int n));


void FwdJacobian(int n, double x, double y[], double yp[], double **dfdy,
	void (*derivs)(double x, double y[], double yp[]));

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void stiff0(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs0)(double, double [], double []))
{
	int i,j,jtry,*indx;
	double d,errmax,h,xsav,**a,*dfdx,**dfdy,*dysav,*err;
	double *g1,*g2,*g3,*g4,*ysav;

	indx =ivector0(n);     	a    = matrix0(n,n);
	dfdx =vector0(n);     	dfdy = matrix0(n,n);

	dysav=vector0(n);
      err  =vector0(n);
	g1   =vector0(n);
	g2   =vector0(n);
	g3   =vector0(n);
	g4   =vector0(n);
	ysav =vector0(n);
	xsav =(*x);

	for (i=0;i<n;i++)
      	{
		ysav[i]=y[i];
		dysav[i]=dydx[i];
		dfdx[i] = 0; // oleg
		}

//	jacobn0(xsav,ysav,dfdx,dfdy,n);

	FwdJacobian(n,xsav,ysav,dydx,dfdy,derivs0);   // Oleg

	h=htry;
//------------------------------------------------
	for (jtry=0; jtry<MAXTRY; jtry++)
      	{
		for (i=0; i<n; i++) {
			for (j=0;j<n;j++) a[i][j] = -dfdy[i][j];
			a[i][i] += 1.0/(GAM*h);
			}

		ludcmp0(a,n,indx,&d);
		for (i=0;i<n;i++)   g1[i]=dysav[i]+h*C1X*dfdx[i];

		lubksb0(a,n,indx,g1);
		for (i=0;i<n;i++)  y[i]=ysav[i]+A21*g1[i];

		*x=xsav+A2X*h;
		(*derivs0)(*x,y,dydx);
		for (i=0;i<n;i++)  g2[i]=dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h;

		lubksb0(a,n,indx,g2);
		for (i=0;i<n;i++)  y[i]=ysav[i]+A31*g1[i]+A32*g2[i];

		*x=xsav+A3X*h;
		(*derivs0)(*x,y,dydx);
		for (i=0;i<n;i++)  g3[i]=dydx[i]+h*C3X*dfdx[i]+(C31*g1[i]+C32*g2[i])/h;

		lubksb0(a,n,indx,g3);
		for (i=0;i<n;i++) g4[i]=dydx[i]+h*C4X*dfdx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;

		lubksb0(a,n,indx,g4);
		for (i=0;i<n;i++)
            	{
			y[i]=ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i];
			err[i]=E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i];
			}

		*x=xsav+h;
		if (*x == xsav) nrerror("stepsize not significant in stiff");
		errmax=0.0;

		for (i=0;i<n;i++) errmax=max(errmax,fabs(err[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0)
            	{
			*hdid=h;
			*hnext=(errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);
			free_vector0(ysav);
			free_vector0(g4);
			free_vector0(g3);
			free_vector0(g2);
			free_vector0(g1);
			free_vector0(err);
			free_vector0(dysav);
			free_vector0(dfdx);
			free_matrix0(dfdy,n,n);
			free_matrix0(a,n,n);
			free_ivector0(indx);
			return;
			}
            else  {
			*hnext=SAFETY*h*pow(errmax,PSHRNK);
			h=(h >= 0.0 ? max(*hnext,SHRNK*h) : min(*hnext,SHRNK*h));
			}
		}

nrerror("exceeded MAXTRY in stiff");
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void stiff0_OT(double *y, double *dydx, int n, double *x, double htry, double eps,
	double *yscal, double *hdid, double *hnext, double **dfdy,
	void (*derivs0)(double x, double *y, double *yp),
      void (*jacobnL)(double x, double *y, double *dfdx, double **dfdy, int n)
      )
{
	int i,j,jtry,*indx;
	double d,errmax,h,xsav,**a,*dfdx,*dysav,*err;
	double *g1,*g2,*g3,*g4,*ysav;

	indx =ivector0(n);     	a    = matrix0(n,n);
	dfdx =vector0(n);

	dysav=vector0(n);
      err  =vector0(n);
	g1   =vector0(n);
	g2   =vector0(n);
	g3   =vector0(n);
	g4   =vector0(n);
	ysav =vector0(n);
	xsav =(*x);

	for (i=0;i<n;i++)
      	{
		ysav[i]=y[i];
		dysav[i]=dydx[i];
		}

	jacobnL(xsav,ysav,dfdx,dfdy,n);
	h=htry;
//------------------------------------------------
	for (jtry=0; jtry<MAXTRY; jtry++)
      	{
		for (i=0; i<n; i++) {
			for (j=0;j<n;j++) a[i][j] = -dfdy[i][j];
			a[i][i] += 1.0/(GAM*h);
			}

		ludcmp0(a,n,indx,&d);
		for (i=0;i<n;i++)   g1[i]=dysav[i]+h*C1X*dfdx[i];

		lubksb0(a,n,indx,g1);
		for (i=0;i<n;i++)  y[i]=ysav[i]+A21*g1[i];

		*x=xsav+A2X*h;
		(*derivs0)(*x,y,dydx);
		for (i=0;i<n;i++)  g2[i]=dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h;

		lubksb0(a,n,indx,g2);
		for (i=0;i<n;i++)  y[i]=ysav[i]+A31*g1[i]+A32*g2[i];

		*x=xsav+A3X*h;
		(*derivs0)(*x,y,dydx);
		for (i=0;i<n;i++)  g3[i]=dydx[i]+h*C3X*dfdx[i]+(C31*g1[i]+C32*g2[i])/h;

		lubksb0(a,n,indx,g3);
		for (i=0;i<n;i++) g4[i]=dydx[i]+h*C4X*dfdx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;

		lubksb0(a,n,indx,g4);
		for (i=0;i<n;i++)
            	{
			y[i]=ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i];
			err[i]=E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i];
			}

		*x=xsav+h;
		if (*x == xsav) nrerror("stepsize not significant in stiff");
		errmax=0.0;

		for (i=0;i<n;i++) errmax=max(errmax,fabs(err[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0)
            	{
			*hdid=h;
			*hnext=(errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);
			free_vector0(ysav);
			free_vector0(g4);
			free_vector0(g3);
			free_vector0(g2);
			free_vector0(g1);
			free_vector0(err);
			free_vector0(dysav);
			free_vector0(dfdx);
			free_matrix0(a,n,n);
			free_ivector0(indx);
			return;
			}
            else  {
			*hnext=SAFETY*h*pow(errmax,PSHRNK);
			h=(h >= 0.0 ? max(*hnext,SHRNK*h) : min(*hnext,SHRNK*h));
			}
		}

nrerror("exceeded MAXTRY in stiff");
}
#undef SAFETY
#undef GROW
#undef PGROW
#undef SHRNK
#undef PSHRNK
#undef ERRCON
#undef MAXTRY
#undef GAM
#undef A21
#undef A31
#undef A32
#undef C21
#undef C31
#undef C32
#undef C41
#undef C42
#undef C43
#undef B1
#undef B2
#undef B3
#undef B4
#undef E1
#undef E2
#undef E3
#undef E4
#undef C1X
#undef C2X
#undef C3X
#undef C4X
#undef A2X
#undef A3X
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software $=)J"!:V. */
