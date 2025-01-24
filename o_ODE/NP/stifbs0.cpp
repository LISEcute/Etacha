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

extern double* vector0(int dimension);
extern int* ivector0(int dimension);
extern double** matrix0(int dimx, int dimy);
extern void free_vector0(double *vv);
extern void free_ivector0(int *vv);
extern void free_matrix0(double **aa, int dimx, int /*dimy*/);

extern void nrerror(char error_text[]);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

//extern void jacobn0(double x, double y[], double dfdx[], double **dfdy, int n);
extern void FwdJacobian(int n, double x, double y[], double yp[], double **dfdy,
	void (*derivs)(double x, double y[], double yp[]));

extern void simpr0(double y[], double dydx[], double dfdx[], double **dfdy,
		int n, double xs, double htot, int nstep, double yout[],
		void (*derivs)(double, double [], double []));

//-------------------------------------------------------------------------------------------

void stifbs0(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs0)(double, double [], double []));

/*void stifbs0_OT(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext, double **dfdy,
	void (*derivs0)(double, double [], double []),
      void (*jacobnL)(double x, double *y, double *dfdx, double **dfdy, int n));      */

void stifbs0_OT(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext, double **dfdy,
	void (*derivs0)(double, double [], double []),
      void (*FwdJacobian)(int n, double x, double y[], double yp[], double **dfdy, void *));

void pzextr0(int iest, double xest, double yest[], double yz[], double dy[],
		int nv, double *x, double **d);
//-------------------------------------------------------------------------------------------
//extern double **d,*x;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void stifbs0(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs0)(double, double [], double []))
{
	int i,iq,k,kk,km;
	static int first=1,kmax,kopt,nvold = -1;
	static double epsold = -1.0,xnew;
	double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
	double *dfdx,**dfdy,*err,*yerr,*ysav,*yseq;
	static double a[IMAXX];
	static double alf[KMAXX][KMAXX];
	static int nseq[IMAXX]={2,6,10,14,22,34,50,70};
	int reduct,exitflag=0;

      double **d,*x;

	d    = matrix0(nv,KMAXX);
      dfdy = matrix0(nv,nv);

      dfdx = vector0(nv);
	err  = vector0(KMAXX);
	x    = vector0(KMAXX);
	yerr = vector0(nv);
	ysav = vector0(nv);
	yseq = vector0(nv);

//------------------
	if(eps != epsold || nv != nvold)
      	{
		*hnext = xnew = -1.0e29;
		eps1=SAFE1*eps;
		a[0]=nseq[0]+1;
		for (k=0;k<KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (iq=1;iq<KMAXX;iq++)
            	{
			for (k=0;k<iq;k++)
				alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/
					((a[iq+1]-a[0]+1.0)*(2*k+3))));              	// !!!
			}
		epsold=eps;
		nvold=nv;
		a[0] += nv;										// !!!
		for (k=0;k<KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (kopt=1;kopt<KMAXX-1;kopt++)                                  // !!
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
		}

//------------------
	h=htry;
	for (i=0;i<nv;i++)
      			{
                        ysav[i]=y[i];
                        dfdx[i]=0;  // Oleg
                        }

//	jacobn0(*xx,y,dfdx,dfdy,nv);
     	FwdJacobian(nv,*xx,ysav,dydx,dfdy,derivs0);   // Oleg


	if (*xx != xnew || h != (*hnext)) {
		first=1;
		kopt=kmax;
		}

	reduct=0;
	for (;;)
      	{
		for (k=0;k<=kmax;k++)
            	{
			xnew=(*xx)+h;
			if (xnew == (*xx)) nrerror("step size underflow in stifbs");

			simpr0(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq,derivs0);
			xest=SQR(h/nseq[k]);
			pzextr0(k,xest,yseq,y,yerr,nv,x,d);

			if (k != 0)
                  	{
				errmax=TINY;
				for (i=0;i<nv;i++) errmax=max(errmax,fabs(yerr[i]/yscal[i]));
				errmax /= eps;
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+3));
				}
			if (k != 0 && (k >= kopt-1 || first))
                  	{
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
	for (kk=0;kk<=km;kk++)
      	{
		fact=max(err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin)
            	{
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
			}
		}

	*hnext=h/scale;
	if (kopt >= k && kopt != kmax && !reduct)
      	{
		fact=max(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			*hnext=h/fact;
			kopt++;
			}
		}
	free_vector0(yseq);
	free_vector0(ysav);
	free_vector0(yerr);
	free_vector0(x);
	free_vector0(err);
	free_vector0(dfdx);
	free_matrix0(dfdy,nv,nv);
	free_matrix0(d,nv,KMAXX);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*void stifbs0_OT(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext, double **dfdy,
	void (*derivs0)(double, double [], double []),
	void (*jacobn0)(double x, double y[], double dfdx[], double **dfdy, int n)
      )
*/
void stifbs0_OT(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext, double **dfdy,
	void (*derivs0)(double, double [], double []),
//	void (*jacobn0)(double x, double y[], double dfdx[], double **dfdy, int n)
      void (*FwdJacobian)(int n, double x, double y[], double yp[], double **dfdy, void *))

{
	int i,iq,k,kk,km;
	static int first=1,kmax,kopt,nvold = -1;
	static double epsold = -1.0,xnew;
	double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
	double *dfdx,*err,*yerr,*ysav,*yseq;
	static double a[IMAXX];
	static double alf[KMAXX][KMAXX];
	static int nseq[IMAXX]={2,6,10,14,22,34,50,70};
	int reduct,exitflag=0;

      double **d,*x;

	d    = matrix0(nv,KMAXX);

      dfdx = vector0(nv);
	err  = vector0(KMAXX);
	x    = vector0(KMAXX);
	yerr = vector0(nv);
	ysav = vector0(nv);
	yseq = vector0(nv);

//------------------
	if(eps != epsold || nv != nvold)
      	{
		*hnext = xnew = -1.0e29;
		eps1=SAFE1*eps;
		a[0]=nseq[0]+1;
		for (k=0;k<KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (iq=1;iq<KMAXX;iq++)
            	{
			for (k=0;k<iq;k++)
				alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/
					((a[iq+1]-a[0]+1.0)*(2*k+3))));              	// !!!
			}
		epsold=eps;
		nvold=nv;
		a[0] += nv;										// !!!
		for (k=0;k<KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (kopt=1;kopt<KMAXX-1;kopt++)                                  // !!
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
		}

//------------------
	h=htry;
	for (i=0;i<nv;i++) ysav[i]=y[i];

//	jacobn0(*xx,y,dfdx,dfdy,nv);
      FwdJacobian(nv, *xx, ysav, dydx, dfdy, derivs0); // Oleg

	if (*xx != xnew || h != (*hnext)) {
		first=1;
		kopt=kmax;
		}

	reduct=0;
	for (;;)
      	{
		for (k=0;k<=kmax;k++)
            	{
			xnew=(*xx)+h;
			if (xnew == (*xx))
                  		 nrerror("step size underflow in stifbs");
			simpr0(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq,derivs0);
			xest=SQR(h/nseq[k]);
			pzextr0(k,xest,yseq,y,yerr,nv,x,d);

			if (k != 0)
                  	{
				errmax=TINY;
				for (i=0;i<nv;i++) errmax=max(errmax,fabs(yerr[i]/yscal[i]));
				errmax /= eps;
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+3));
				}
			if (k != 0 && (k >= kopt-1 || first))
                  	{
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
	for (kk=0;kk<=km;kk++)
      	{
		fact=max(err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin)
            	{
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
			}
		}

	*hnext=h/scale;
	if (kopt >= k && kopt != kmax && !reduct)
      	{
		fact=max(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			*hnext=h/fact;
			kopt++;
			}
		}
	free_vector0(yseq);
	free_vector0(ysav);
	free_vector0(yerr);
	free_vector0(x);
	free_vector0(err);
	free_vector0(dfdx);
	free_matrix0(d,nv,KMAXX);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void pzextr0(int iest, double xest, double yest[], double yz[], double dy[], int nv, double *x, double **d)
{
	int k1,j;                                                                       //page 735
	double q,f2,f1,delta,*c;

	c=vector0(nv);
	x[iest]=xest;
	for (j=0;j<nv;j++) dy[j]=yz[j]=yest[j];

	if (iest == 0)
		for (j=0;j<nv;j++) d[j][0]=yest[j];
	else
      	{

		for (j=0;j<nv;j++) c[j]=yest[j];

		for (k1=0;k1<iest;k1++)
            	{
			delta=1.0/(x[iest-k1-1]-xest);                  // !!!
			f1=xest*delta;
			f2=x[iest-k1-1]*delta;  				// !!!
			for (j=0;j<nv;j++)
                  	{
				q=d[j][k1];
				d[j][k1]=dy[j];
				delta=c[j]-q;
				dy[j]=f1*delta;
				c[j]=f2*delta;
				yz[j] += dy[j];
				}
			}
		for (j=0;j<nv;j++) d[j][iest]=dy[j];
		}

free_vector0(c);
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
