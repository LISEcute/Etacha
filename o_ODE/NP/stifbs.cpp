#include <math.h>
#define NRANSI
//#include "o_phys\NP\nrutil.h"
#define KMAXX 7
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

#define max(a,b)  (a>b ? a : b)
#define min(a,b)  (a<b ? a : b)
#define SQR(a)   (a == 0.0 ? 0.0 : a*a)

extern double* vector(int start, int dimension);
extern int* ivector(int start, int dimension);
extern double** matrix(int startx, int dimx, int starty, int dimy);
extern void free_vector(double *vv, int start, int dimension);
extern void free_ivector(int *vv, int /*start*/, int /*dimension*/);
extern void free_matrix(double **aa, int startx, int dimx, int /*starty*/, int /*dimy*/);

extern void nrerror(char error_text[]);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

extern void jacobn(double x, double y[], double dfdx[], double **dfdy, int n);
extern void simpr(double y[], double dydx[], double dfdx[], double **dfdy,
		int n, double xs, double htot, int nstep, double yout[],
		void (*derivs)(double, double [], double []));

//-------------------------------------------------------------------------------------------
void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
		int nv, double *x, double **d);

//-------------------------------------------------------------------------------------------
//extern double **d,*x;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []))
{
	int i,iq,k,kk,km;
	static int first=1,kmax,kopt,nvold = -1;
	static double epsold = -1.0,xnew;
	double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
	double *dfdx,**dfdy,*err,*yerr,*ysav,*yseq;
	static double a[IMAXX+1];
	static double alf[KMAXX+1][KMAXX+1];
	static int nseq[IMAXX+1]={0,2,6,10,14,22,34,50,70};
	int reduct,exitflag=0;

      double **d,*x;

	d=matrix(1,nv,1,KMAXX);
      dfdy=matrix(1,nv,1,nv);

      dfdx=vector(1,nv);
	err=vector(1,KMAXX);
	x=vector(1,KMAXX);
	yerr=vector(1,nv);
	ysav=vector(1,nv);
	yseq=vector(1,nv);


	if(eps != epsold || nv != nvold) {
		*hnext = xnew = -1.0e29;
		eps1=SAFE1*eps;
		a[1]=nseq[1]+1;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (iq=2;iq<=KMAXX;iq++) {
			for (k=1;k<iq;k++)
				alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/
					((a[iq+1]-a[1]+1.0)*(2*k+1))));
		}
		epsold=eps;
		nvold=nv;
		a[1] += nv;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (kopt=2;kopt<KMAXX;kopt++)
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
	}
	h=htry;
	for (i=1;i<=nv;i++) ysav[i]=y[i];
	jacobn(*xx,y,dfdx,dfdy,nv);
	if (*xx != xnew || h != (*hnext)) {
		first=1;
		kopt=kmax;
	}
	reduct=0;
	for (;;) {
		for (k=1;k<=kmax;k++) {
			xnew=(*xx)+h;
			if (xnew == (*xx)) nrerror("step size underflow in stifbs");
			simpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq,derivs);
			xest=SQR(h/nseq[k]);
			pzextr(k,xest,yseq,y,yerr,nv,x,d);
			if (k != 1) {
				errmax=TINY;
				for (i=1;i<=nv;i++) errmax=max(errmax,fabs(yerr[i]/yscal[i]));
				errmax /= eps;
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
			}
			if (k != 1 && (k >= kopt-1 || first)) {
				if (errmax < 1.0) {
					exitflag=1;
					break;
				}
				if (k == kmax || k == kopt+1) {
					red=SAFE2/err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
						red=1.0/err[km];
						break;
					}
				else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
						red=alf[km][kmax-1]*SAFE2/err[km];
						break;
					}
				else if (alf[km][kopt] < err[km]) {
					red=alf[km][kopt-1]/err[km];
					break;
				}
			}
		}
		if (exitflag) break;
		red=min(red,REDMIN);
		red=max(red,REDMAX);
		h *= red;
		reduct=1;
	}
	*xx=xnew;
	*hdid=h;
	first=0;
	wrkmin=1.0e35;
	for (kk=1;kk<=km;kk++) {
		fact=max(err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin) {
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
		}
	}
	*hnext=h/scale;
	if (kopt >= k && kopt != kmax && !reduct) {
		fact=max(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			*hnext=h/fact;
			kopt++;
		}
	}
	free_vector(yseq,1,nv);
	free_vector(ysav,1,nv);
	free_vector(yerr,1,nv);
	free_vector(x,1,KMAXX);
	free_vector(err,1,KMAXX);
	free_matrix(dfdy,1,nv,1,nv);
	free_vector(dfdx,1,nv);
	free_matrix(d,1,nv,1,KMAXX);
}



//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv, double *x, double **d)
{
	int k1,j;
	double q,f2,f1,delta,*c;

	c=vector(1,nv);
	x[iest]=xest;
	for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
	if (iest == 1) {
		for (j=1;j<=nv;j++) d[j][1]=yest[j];
	} else {
		for (j=1;j<=nv;j++) c[j]=yest[j];
		for (k1=1;k1<iest;k1++) {
			delta=1.0/(x[iest-k1]-xest);
			f1=xest*delta;
			f2=x[iest-k1]*delta;
			for (j=1;j<=nv;j++) {
				q=d[j][k1];
				d[j][k1]=dy[j];
				delta=c[j]-q;
				dy[j]=f1*delta;
				c[j]=f2*delta;
				yz[j] += dy[j];
			}
		}
		for (j=1;j<=nv;j++) d[j][iest]=dy[j];
	}
	free_vector(c,1,nv);
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software $=)J"!:V. */
