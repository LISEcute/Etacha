#include <math.h>
#define NRANSI
//#include "nrutil.h"
#define FMAX(a,b)  (a>b ? a : b)
#define FMIN(a,b)  (a<b ? a : b)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

extern double* vector(int start, int dimension);
extern int* ivector(int start, int dimension);
extern double** matrix(int startx, int dimx, int starty, int dimy);
extern void free_vector(double *vv, int start, int dimension);
extern void free_ivector(int *vv, int /*start*/, int /*dimension*/);
extern void free_matrix(double **aa, int startx, int dimx, int /*starty*/, int /*dimy*/);

extern void nrerror(char error_text[]);

//#define NRANSI
//#include "nrutil.h"
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double []));


void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));

void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []));

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double []))
{
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.00/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

	ak2=vector(1,n);
	ak3=vector(1,n);
	ak4=vector(1,n);
	ak5=vector(1,n);
	ak6=vector(1,n);
	ytemp=vector(1,n);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,ytemp,ak2);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)(x+a3*h,ytemp,ak3);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)(x+a4*h,ytemp,ak4);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)(x+a5*h,ytemp,ak5);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs)(x+a6*h,ytemp,ak6);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=1;i<=n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
	free_vector(ytemp,1,n);
	free_vector(ak6,1,n);
	free_vector(ak5,1,n);
	free_vector(ak4,1,n);
	free_vector(ak3,1,n);
	free_vector(ak2,1,n);
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []))
{
	int i;
	double errmax,h,htemp,xnew,*yerr,*ytemp;

	yerr=vector(1,n);
	ytemp=vector(1,n);
	h=htry;
	for (;;) {
		rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
		errmax=0.0;
		for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax > 1.0) {
			htemp=SAFETY*h*pow(errmax,PSHRNK);
			h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
			xnew=(*x)+h;
			if (xnew == *x) nrerror("stepsize underflow in rkqs");
			continue;
		} else {
			if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
			else *hnext=5.0*h;
			*x += (*hdid=h);
			for (i=1;i<=n;i++) y[i]=ytemp[i];
			break;
		}
	}
	free_vector(ytemp,1,n);
	free_vector(yerr,1,n);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software $=)J"!:V. */


#define MAXSTP 50010
#define TINY 1.0e-30

extern int kmaxO, kountO;
extern double *xpO,**ypO,  dxsav;

/*void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,

	void (*derivs)(double, double [], double []),

	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])))
*/

void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []))
{
	int nstp,i;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;

	yscal=vector(1,nvar);
	y=vector(1,nvar);
	dydx=vector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kountO = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmaxO > 0) xsav=x-dxsav*2.0;

	for (nstp=1;nstp<=MAXSTP;nstp++)
      	{
		(*derivs)(x,y,dydx);

		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;

		if (kmaxO > 0 && kountO < kmaxO-1 && fabs(x-xsav) > fabs(dxsav)) {
			xpO[++kountO]=x;
			for (i=1;i<=nvar;i++) ypO[i][kountO]=y[i];
			xsav=x;
			}

		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;

		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);

		if (hdid == h) ++(*nok); else ++(*nbad);

		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmaxO) {
				xpO[++kountO]=x;
				for (i=1;i<=nvar;i++) ypO[i][kountO]=y[i];
				}
			free_vector(dydx,1,nvar);
			free_vector(y,1,nvar);
			free_vector(yscal,1,nvar);
			return;
			}

		if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
		h=hnext;
		}
            
nrerror("Too many steps in routine odeint");
}
#undef MAXSTP
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software $=)J"!:V. */
