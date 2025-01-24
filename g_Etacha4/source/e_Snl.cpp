#include "../win/e_myextern.h"
#include "../win/e_Constant.h"
#include <QtMath>
#include <stdio.h>

extern double Velocity_au(double E);
extern double E_to_Beta(double E);

extern double pow1(double par);
extern double pow2(double par);
extern double pow3(double par);
extern double pow4(double par);
extern double pow5(double par);
extern double pow6(double par);
extern double pow8(double par);
extern double pow10(double par);
extern double pow2I(double par);

extern void INTG( double (*funct)(double x), double BOUND, int INF,
                 double EPSABS, double EPSREL, double &RESULT, double &ABSERR,
                 int &NEVAL, int &IER, int klim, int klenw, int &LAST,
                 int *IWORK, double *WORK);


void snl(double Ecurr, double &esp,double &esp3,double &epd3);
double FNL(double xx);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void snl(double Ecurr, double &esp,double &esp3,double &epd3)
{
//c*********************************************************************
//c        Calcule les sections efficaces d'excitation ns-np  dans
//c    l'approximation PWBA avec:
//c        -facteurs de forme analytiques utilisant la forme recurrente
//c        etablie a partir des calculs n=2 a n=6
//c        -des facteurs d'ecrantage et d'antiecrantage Thomas-Fermi
//c    Les facteurs de forme sont calcules en fonction de la variable
//c    reduite k =  q/Zp (q en u.a.)
//c
//c    Les resultats sont donnes ici en unites de 10-20 cm2
//c*********************************************************************
      double secnl[4];


      double BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK[etaLENWc];
      int INF,NEVAL,IER,LAST,IWORK[etaLIMITc];

      for (int i=0; i<etaLENWc;  i++)  WORK[i]=0;
      for (int i=0; i<etaLIMITc; i++) IWORK[i]=0;

      double beta = E_to_Beta(Ecurr);


//c******integration***********

      for(int n=1; n<=3; n++)
        {
        nsp_n = n;

        //*************************************************
        BOUND = 1.E-5/Zp;
        INF = 1;
        EPSABS = 0.0;
        EPSREL = 1.E-4;

        INTG(FNL,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
                  IER,etaLIMIT,etaLENW,LAST,IWORK,WORK);
        double sig1=RESULT;

               if (IER > 0)
                     {
                      FILE *f = fopen("messerr4L.fil","wt");
                      fprintf(f," n=%d  ier=%d \n",n,IER);
                      fprintf(f," result=%g   abserr=%g\n",RESULT,ABSERR);
                      fclose(f);;
                      }
        //*************************************************
         secnl[n]=sig1*(1.874)*Zt*Zt/(Zp*Zp*beta*beta);
         } //continue;


      esp =secnl[1];
      esp3=secnl[2];
      epd3=secnl[3];
}
//c*************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double FNL(double xx)
{

      double      q = xx*Zp;
      double    scf = 1.13*pow(Zt,1./3.);
      double    stf = scf / q;
      double    atf = q / scf;
      double   corr1 = pow2I(1.+pow2(stf));
      double   corr2 = (1.-pow2I(1.+pow2(atf)))/Zt ;
      double   corr = corr1 + corr2;

//c*********************************************************************
      double bp  = xx;
      double bp2 = xx*xx;
      double fac=1,t;



      if (nsp_n == 1) {                                         // 2s-2p
      t = 1.-bp2;
      fac=18.*pow2(t)/(bp*pow8(1.+bp2));
                      }
 else if (nsp_n == 2) {                                         // 3s-3p
      t=128.-1392.*bp2+3888.*pow4(bp)-2187.*pow6(bp);
      fac=110592.*pow2(t)/(bp*pow((4.+9.*bp2),12.));
      }
 else if (nsp_n == 3) {                                         // 3p-3d
              t=256.-3840.*bp2+30816.*pow4(bp)-81648.*pow6(bp)+76545.*pow8(bp);
              fac=2949120.*t/(bp*pow((4.+9.*bp2),12.));
      }

double fnl=fac*corr;

return fnl;
}
//c*********************************************************************
