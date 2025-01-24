//#pragma hdrstop

#include <math.h>
#include <stdio.h>
#include "../win/e_myextern.h"
#include "../win/e_Constant.h"

extern double Velocity_au(double E);
extern double E_to_Beta(double E);
extern double pow2(double par);
extern double pow2I(double par);
extern double pow6(double par);
extern double pow8(double par);

extern void INTG( double (*funct)(double x), double BOUND, int INF,
                 double EPSABS, double EPSREL, double &RESULT, double &ABSERR,
                 int &NEVAL, int &IER, int klim, int klenw, int &LAST,
                 int *IWORK, double *WORK);


void sexi(double Ecurr, double &e2s,double &e2p,double &e3s,double &e3p,double &e3d,double &es4,double &es5);
double FEX(double x);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void sexi(double Ecurr, double &e2s,double &e2p,double &e3s,double &e3p,double &e3d,double &es4,double &es5)
{
//c         Version du 02/03/94 simplifiee pour eta3
//c           dernieres corrections : 31/08/2000 (COA)
//c    version avec makefile utilisant la routine d'integration INTG
//c*********************************************************************
//c        Calcule les sections efficaces d'excitation 1s-2s, 1s-2p et
//c    1s-n dans l'approximation PWBA avec:
//c        -facteurs de forme analytiques d'Anholt
//c        -facteurs d'ecrantage et d'antiecrantage dans le modele
//c         de Thomas-Fermi
//c    Les facteurs de forme sont calcules en fonction de la variable
//c    reduite k =  q/Zp (q en u.a.)
//c
//c    Les resultats sont donnes en unites de 10-20 cm2
//c*********************************************************************

      double BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK[etaLENWc];
      double sig1,qmin;
      int INF,NEVAL,IER,LAST,IWORK[etaLIMITc];

      for (int i=0; i<etaLENWc;  i++)  WORK[i]=0;
      for (int i=0; i<etaLIMITc; i++) IWORK[i]=0;

      bool Flag34 = EtachaVersion >= etacha_v34;
      int nmin=2;
      int nmax=Flag34 ? 5 : 4;

      //------------------------------------------
      double beta = E_to_Beta(Ecurr);
      double vu = Velocity_au(Ecurr);

      double fs1=1.;

            if (Zt == 1.)  fs1=1.23;
       else if (Zt == 2.)  fs1=1.526;
       else if (Zt == 6.)  fs1=0.78;
       else if (Zt == 7.)  fs1=0.85;
       else if (Zt == 10.) fs1=1.04;
       else if (Zt == 13.) fs1=0.58;
       else if (Zt == 14.) fs1=0.59;
       else if (Zt == 18.) fs1=0.68;
       else if (Zt == 29.) fs1=0.672;
       else if (Zt == 36.) fs1=0.61;
       else if (Zt == 54.) fs1=0.535;

      n_scf = double(fs1*1.13*pow(Zt,1./3.));   // commom local
      double som1=0.;
      double so51=0.;

//c****************************
//c******integration***********
for(int n=nmin; n<=nmax; n++)  /*L20*/
   {
            if (n == 2) n_ll=1;
   else     if (n == 3) n_ll=4;
   else     if (n > 3)  n_ll=7;

L21:    n_rn=n;

        //c**********************************
        //c correction pour effet energie liaison:
        //c**********************************
    double epsi = 1;
    if(ibin == 1)
              {
              double y=2.*vu/(zk*tetak);
              double g=(1.+5.*y+7.14*pow(y,2.)+4.27*pow(y,3.)+0.947*pow(y,4.))/pow((1.+y),5.);
              epsi=pow((1+Zt*g/(zk*tetak)),2.);
              }

//c****************************
                // commented original     //     pmin=deltaE/2v => qmin=deltaE/2v/Zp

  qmin = epsi*(zk*tetak/(vu*2.))*(1.-1./(n_rn*n_rn));
        //     calcul du coefficient d'antiscreening (08/2000)
  double vs=1.75*epsi*zk*tetak*qmin;

  if (vu > vs) n_coa=1.-pow((vs/vu),2);
  else         n_coa=0.;

  BOUND = qmin;
  INF = 1;
  EPSABS = 0.0;
  EPSREL = 1.E-4;
//*************************************************


INTG(FEX,BOUND,INF, EPSABS, EPSREL, RESULT, ABSERR, NEVAL, IER, etaLIMIT, etaLENW, LAST, IWORK, WORK);

               sig1=RESULT;
               if (IER > 0)
                     {
                      FILE *f = fopen("messerr1.fil","wt");
                      fprintf(f," n=%d  ier=%d \n",n,IER);
                      fprintf(f," result=%g   abserr=%g\n",RESULT,ABSERR);
                      fclose(f);;
                      }

//*************************************************
        double sec1=sig1*(1.874)*Zt*Zt/(zk*zk*beta*beta);

        if(n_ll == 3 || n_ll == 7)  som1 += sec1;
        else if (n_ll < 3)
          {
                if(n_ll == 1) e2s=sec1;
           else if(n_ll == 2) e2p=sec1;

           n_ll=n_ll+1;
           goto L21;
           }
       else if (n_ll > 3 && n_ll < 7)
           {
                if(n_ll == 4) e3s=sec1;
           else if(n_ll == 5) e3p=sec1;
           else if(n_ll == 6) e3d=sec1;

           n_ll=n_ll+1;
           goto L21;
           }


                if (n == 4)   es4 = sec1;
           else if (n > 4 ) so51 += sec1;
      //-----------------------------------------------
       } //continue; L20:


// commented original     so51=so51+4.52*sec1
      es5=so51;

}

//c*************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double FEX(double x)
//void FEX(double x)
{

      double q   = x*zk;
      double stf = n_scf/q;
      double atf = q/n_scf;
      double corr  = pow2I(1.+pow2(stf)) +
                     n_coa*(  (1.-pow2I(1.+pow2(atf)))/Zt );


//	corr=dble(1./(1.+stf**2.)**2.)
//	corr=corr+coa*dble( (1.-(1./(1.+atf**2.)**2.))/zt  )

//c*********************************************************************
      double bp=x;
      double bp2=x*x;
      double fact;


      if (n_ll == 1)    fact = 64.     *bp /pow6(bp2+9./4.);
 else if (n_ll == 2)    fact = 144.    /bp /pow6(bp2+9./4.);
 else if (n_ll == 4)    fact = 1119744.*bp*pow2((16.+27.*bp2)) /pow8((9.*bp2+16.));
 else if (n_ll == 5)    fact = 2985984./bp*pow2((16.+27.*bp2)) /pow8((9.*bp2+16.));
 else if (n_ll == 6)    fact = 573308928.*bp                   /pow8((9.*bp2+16.));
 else                   fact = (512./(3.*n_rn*n_rn*n_rn))
                               *   (3.*bp2+1.-1./(n_rn*n_rn))
                                 * pow((bp2+pow2(1.-1./n_rn)),n_rn-3.)
                              /bp/ pow((bp2+pow2(1.+1./n_rn)),n_rn+3.);


      double FEX=fact*corr;
      return FEX;
}
//c*********************************************************************

