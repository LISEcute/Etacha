#include <math.h>
#include <stdio.h>

#include "../win/e_myextern.h"
#include "../win/e_Constant.h"

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


void sex2(double Ecurr, double &s3s,double &s3p,double &s3d,double &p3s,double &p3p,double &p3d,
     double &e2s4,double &e2p4,double &e2s5,double &e2p5);
void sex3(double Ecurr, double &e3s4,double &e3p4,double &e3d4,double &e3s5,double &e3p5,double &e3d5,
double &e4s5,double &e4p5,double &e45,double &e456);

double FEX2(double x);
double FEX3(double x);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void sex2(double Ecurr, double &s3s,double &s3p,double &s3d,double &p3s,double &p3p,double &p3d,
     double &e2s4,double &e2p4,double &e2s5,double &e2p5)
{
//c         Version du 21/06/2000 simplifiee pour eta4
//c           dernieres corrections : 31/08/2000 (COA et certains facteurs)
//c    version utilisant la routine d'integration INTG
//c*********************************************************************
//c        Calcule les sections efficaces d'excitation 2l-3l' et 2l-4l'
//c        dans l'approximation PWBA avec facteurs d'ecrantage et
//c        d'antiecrantage dans le modele de Thomas-Fermi
//c    Les facteurs de forme sont calcules en fonction de la variable
//c    reduite k = q/Zp (q en u.a.)
//c
//c    Les resultats sont donnes en unites de 10-20 cm2
//c*********************************************************************
      double sec[11];
      for(int i=0; i<11; i++) sec[i]=0;
//**********************************************************************
      double BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK[etaLENWc];
      double sig1,qmin;
      int INF,NEVAL,IER,LAST,IWORK[etaLIMITc];

      for (int i=0; i<etaLENWc;  i++)  WORK[i]=0;
      for (int i=0; i<etaLIMITc; i++) IWORK[i]=0;

      bool Flag34 = EtachaVersion >= etacha_v34;
      int nmax=Flag34 ? 10 : 8;
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

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWww
//c****************************
//c******integration***********
for(int n=1; n<=nmax; n++) // L20
        {
        n_rn = n;
        //c**********************************
        //c correction pour effet energie liaison:
        //c**********************************
         double   epsi=1.;

         if (ibin == 1) {
              double y = 4.*vu/(zl2*tetal2);
              double g=1.+9.*y+34.7*pow2(y)+81.2*pow3(y)+112.*pow4(y)
                         +93.5*pow(y,5.)+46.6*pow(y,6.)+12.9*pow(y,7.)+1.549*pow(y,8.);
              g /= pow((1.+y),9.);

              epsi = pow2(1+Zt*g/(zl2*tetal2));
              }

//c****************************
        double temp = epsi*(zl2*tetal2/(vu*2.));

              if(n <= 6)                      //      2 -> 3
                      qmin=temp*(5./36.);
         else if(n <= 8)                       //      2 -> 4
                      qmin=temp*(3./16.);
         else                                 //      2 -> 5
                      qmin=temp*(21./100.);


                //     calcul du coefficient d'antiscreening (08/2000)
        double vs=1.75*epsi*zl2*tetal2*qmin;

        if (vu > vs)   n_coa=1.-pow2(vs/vu);
        else           n_coa=0.;

      BOUND = qmin;
      INF = 1;
      EPSABS = 0.0;
      EPSREL = 1.E-4;
//*************************************************
      INTG(FEX2,BOUND,INF,
                EPSABS,EPSREL,RESULT,ABSERR,
                NEVAL,IER,etaLIMIT,etaLENW,LAST,
                IWORK,WORK);

/*      CALL INTG[FEX2,BOUND,INF,
                EPSABS,EPSREL,RESULT,ABSERR,
                NEVAL,IER,LIMIT,LENW,LAST,
                IWORK,WORK];*/


               sig1=RESULT;
               if (IER > 0)
                     {
                      FILE *f = fopen("messerr2.fil","wt");
                      fprintf(f," n=%d  ier=%d \n",n,IER);
                      fprintf(f," result=%g   abserr=%g\n",RESULT,ABSERR);
                      fclose(f);;
                      }
        //*************************************************
        sec[n] = sig1*(1.874)*Zt*Zt/(zl2*zl2*beta*beta);

         } //L20
//*************************************************
//*************************************************



      s3s=sec[1];
      s3p=sec[2];
      s3d=sec[3];
      p3s=sec[4];
      p3p=sec[5];
      p3d=sec[6];
      e2s4=sec[7];
      e2p4=sec[8];
      e2s5=sec[9];
      e2p5=sec[10];
      }
//c*************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double FEX2(double x)
{

      double      q = x*zl2;
      double    stf = n_scf / q;
      double    atf = q / n_scf;
      double   corr1 = pow2I(1.+pow2(stf));
      double   corr2 = n_coa*( (1.-pow2I(1.+pow2(atf)))/Zt   );
      double   corr = corr1 + corr2;

//c*********************************************************************
      double bp  = x;
      double bp2 = x*x;
      double bp4 = bp2*bp2;
      double bp6 = bp2*bp4;
      double fact;

       if (n_rn == 1) {
//c 2s->3s
//  pow(2.,18.)*pow(3.,7.) ==  573308928
      fact = 573308928.*bp / pow10(25.+36.*bp2)
           *pow2(2875.-6984.*bp2+3888.*bp4);
      }
 else if (n_rn == 2) {
//c 2s->3p
//    PP 2.*pow2(27648.)=1528823808
      fact=1528823808. /(bp*pow10(25.+36.*bp2))
             *pow2((-625.+5040.*bp2-3888.*bp4));
      }
 else if (n_rn == 3) {
//c 2s->3d
//  2.*pow2(65536.)*pow(3.,7.) == 18786186952704
      fact=18786186952704.*bp/pow10(25.+36.*bp2)
           *pow2((-25.+18.*bp2));
      }
 else if (n_rn == 4) {
//c 2p->3s
//c corrige le 19/08/2009 : 3.**7 remplace par 3.**6
//         pow(2.,16.)*pow(3.,6.) ==  47775744
      fact= 47775744. /(bp*pow10(25.+36.*bp2)) *pow2(625.-14040.*bp2+11664.*bp4);
      }

 else if (n_rn == 5) {
//c 2p->3p
//   pow(2.,25.)*pow(3.,9.) ==  6.60452E+11
      fact=6.60451885E+11*bp/pow10(25.+36.*bp2)
          *(6875.-23400.*bp2+34992.*bp4);
      }
 else if (n_rn == 6) {
//c 2p->3d
// pow(2.,24.)*pow(3.,6.)*pow(5.,2.) ==  3.057647616E+11
      fact=3.057647616E+11/(bp*pow10(25.+36.*bp2))
      *(3125.-5400.*bp2+27216.*bp4);
      }
 else if (n_rn == 7) {
//c 2s->4
      fact=405.+5120.*bp2+98816.*bp4
        -163840.*bp6+65536.*pow8(bp);
      fact=fact*pow(2.,20.)/(bp*pow((9.+16.*bp2),9.));
      }
 else if (n_rn == 8) {
//c 2p->4
      fact=4194304.*(123. + 3728.*bp2- 3840.*bp4+ 12288.*bp6)
                /(bp*pow((9.+16.*bp2),9));
      }
 else if (n_rn == 9) {
//c 2s->5

      fact=(20480000000.*(1555848.+29336825.*bp2+376570000.*bp4+
                3689500000.*bp6-6300000000.*pow(bp,8)+2500000000.*pow10(bp)))
                / (bp*pow10(49.+100.*bp2));
      }
 else {
//c 2p->5
      fact=(128000000000.*(1643158523433.+74683227547040.*bp2+
                993723878526400.*bp4+4820995175360000.*pow6(bp)+
                14083087256000000.*pow8(bp)+45107440000000000.*pow10(bp)+
                81423200000000000.*pow(bp,12)-79200000000000000.*pow(bp,14)+
                90000000000000000.*pow(bp,16)))
                /(bp*pow(49.+100.*bp2,14));
      }

      double FEX2a=fact*corr;
      return FEX2a;
      }
//c*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void sex3(double Ecurr, double &e3s4,double &e3p4,double &e3d4,double &e3s5,double &e3p5,double &e3d5,
          double &e4s5,double &e4p5,double &e45,double &e456)
 {
//c         Version du 25/01/95 simplifiee pour eta4
//c           dernieres corrections : 31/08/2000 (COA)
//c    version utilisant la routine d'integration INTG
//c*********************************************************************
//c      Calcule les sections efficaces d'excitation 3l-4l'
//c        dans l'approximation PWBA avec facteurs d'ecrantage et
//c        d'antiecrantage dans le modäle de Thomas-Fermi
//c    Les facteurs de forme sont calcules en fonction de la variable
//c    reduite k = q/Zp (q en u.a.)
//c
//c    Les resultats sont donnes en unites de 10-20 cm2
//c*********************************************************************

      double sec[11];
      for(int i=0; i<11; i++) sec[i]=0;
//**********************************************************************
      double BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK[etaLENWc];
      double sig1,qmin;
      int INF,NEVAL,IER,LAST,IWORK[etaLIMITc];

      for (int i=0; i<etaLENWc;  i++)  WORK[i]=0;
      for (int i=0; i<etaLIMITc; i++) IWORK[i]=0;

      bool Flag34 = EtachaVersion >= etacha_v34;
      int nmax=Flag34 ? 10 : 3;
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

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWww
//c******integration***********
for(int n=1; n<=nmax; n++)   // L20
        {
        n_rn=n;
                //c**********************************
                //c correction pour effet energie liaison:
                //c**********************************
    double epsi = 1;
           if(ibin == 1)
                {
                double y = 18.*vu/(zm2*tetam2);
                double g= (1.+5.*y+7.14*pow2(y)+4.27*pow(y,3.)+0.947*pow(y,4.))/pow((1.+y),5.);
                epsi = pow2(1.+Zt*g/(zm2*tetam2));
                }
//c****************************
        double vs;

              if (n <= 3) {                                            //     3 vers 4
                      qmin=epsi*(zm2*tetam2/(vu*2.))*(7./144.);
                      vs = 1.75*epsi*zm2*tetam2*qmin;
              }
         else if (n <= 6) {                                             //     3 vers 5
              qmin= epsi*(zm2*tetam2/(vu*2.))*(16./225.);
              vs=1.75*epsi*zm2*tetam2*qmin ;
              }
         else if (n <= 9) {                                             //     4 vers 5
              qmin= (zn*tetan/(vu*2.))*(9./400.);
              vs=1.75*epsi*zn*tetan*qmin ;
              }
         else {                                                         //     4 vers 6
              qmin= (zn*tetan/(vu*2.))*(5./144.);
              vs=1.75*epsi*zn*tetan*qmin;
              }


//     calcul du coefficient d'antiscreening (08/2000)
       if (vu > vs)  n_coa = 1.-pow2(vs/vu);
       else          n_coa = 0.;

      BOUND = qmin;
      INF = 1;
      EPSABS = 0.0;
      EPSREL = 1.E-4;
//*************************************************
        INTG(FEX3,BOUND,INF, EPSABS, EPSREL, RESULT, ABSERR, NEVAL, IER, etaLIMIT, etaLENW, LAST, IWORK, WORK);

        sig1=RESULT;

        if (IER > 0)
                     {
                      FILE *f = fopen("messerr3.fil","wt");
                      fprintf(f," n=%d  ier=%d \n",n,IER);
                      fprintf(f," result=%g   abserr=%g\n",RESULT,ABSERR);
                      fclose(f);;
                      }
  //*************************************************

       if(n >= 7) sec[n]=sig1*(1.874)*Zt*Zt/(zn* zn* beta*beta);
       else       sec[n]=sig1*(1.874)*Zt*Zt/(zm2*zm2*beta*beta);

       } //continue;   L20
//------------------------------------

      e3s4=sec[1];
      e3p4=sec[2];
      e3d4=sec[3];
//     (facteur theorique en 1/n**3 = 3.049358 pour n=5)
//     les facteurs sont dans CSEC(EF) et valent 2.5 (13/09/01)
      e3s5=sec[4];
      e3p5=sec[5];
      e3d5=sec[6];
      e4s5=sec[7];
      e4p5=sec[8];
      e45=sec[9];
//     on suppose que la loi en 1/n**3 n'est pas etablie des n=6
//     (facteur theorique en 1/n**3 = 3.5412910805 pour n=6)
//     e456=sec(9)+2.5*sec(10)
      e456=sec[9]+sec[10];
      }
//c*************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double FEX3(double x)
{
double q;

if (n_rn <= 6.)       q = x*zm2;
else                  q = x*zn;

      double stf = n_scf/q;
      double atf = q/n_scf;
      double   corr1 = pow2I(1.+pow2(stf));
      double   corr2 = n_coa*( (1.-pow2I(1.+pow2(atf)))/Zt   );
      double   corr = corr1 + corr2;

//c*********************************************************************
      double bp=x;
      double bp2=x*x;
      double bp4 = bp2*bp2;
      double bp6 = bp2*bp4;
      double fact=0;


if (n_rn == 1) {
//c 3s->4
      fact=4250070125.+110475272992.*bp2;
      fact=fact+251535456000.*bp4;
      fact=fact-8858287816704.*bp6;
      fact=fact+30668348915712.*pow8(bp);
      fact=fact-29926726041600.*pow10(bp);
      fact=fact+8916100448256.*pow(bp,12.);
      fact=fact*18345885696./(bp*pow((49.+144.*bp2),11.));
      }
 else if (n_rn == 2) {
//c 3p->4
// pp pow(2.,30)*pow(3.,5) = 2.609192632E+11
      fact=401065441.+5869103184.*bp2-57396708864.*bp4;
      fact=fact+313343188992.*bp6-651422269440.*pow8(bp);
      fact=fact+557256278016.*pow(bp,10);
      fact=fact*2.609192632E+11/(bp*pow((49.+144.*bp2),11));
      }
 else if (n_rn == 3) {
//c 3d->4
// pp pow(2.,35)*pow(3.,7) ==7.514474781E+13
      fact=11008585.-92829520.*bp2+552524544.*bp4;
      fact=fact+417042432.*bp6+1528823808.*pow8(bp);
      fact=fact*7.514474781E+13/(5.*bp*pow((49.+144.*bp2),11));
      }
 else if (n_rn == 4) {
//c 3s->5
      fact=87480000000.*(183744069632.+7768401510400.*bp2+
                227789475840000.*bp4+3915483300000000.*bp6-
                41894536500000000.*pow8(bp)+115844706562500000.*pow10(bp)-
                103451080078125000.*pow(bp,12)+29192926025390625.*pow(bp,14))
                / (bp*pow(64.+225*bp2,12));
      }
 else if (n_rn == 5) {
//c 3p->5
      fact=1944.e9*(10020192256.+4232282112.e2*bp2+
                1934093376.e4*bp4-1692872865.e5*bp6+
                79597231875.e4*pow8(bp)-152235703125.e4*pow10(bp)+
                1167717041015625.*pow(bp,12))
                /(bp*pow((64.+225.*bp2),12));
      }
 else if (n_rn == 6) {
//c 3d->5
      fact=69984.e11*(3014656+382644224*bp2-33788808.e2*bp4
            +1613793375.e1*bp6+10308515625.*pow8(bp)+512578125.e2*pow10(bp))
             / (bp*pow(64.+225.*bp2,12));
      }
 else if (n_rn == 7) {
//c 4s->5
      fact=(2414107905106248.+122375509655954625.*bp2+
                141571857480192.e4*bp4-42952190122096.e6*bp6-
                3750709452544.e8*pow(bp,8)+898924011392.e10*pow(bp,10)-
                446969856.e14*pow(bp,12)+865128448.e14*pow(bp,14)-
                6488064.e16*pow(bp,16)+16384.e18*pow(bp,18));
      fact=fact*pow(2.,27)*pow(5.,7)/(bp*pow(81.+400.*bp2,14));
      }
 else if (n_rn == 8) {
//c 4p->5
// pp   pow(2.,23)*pow(5.,10) == 8.19200E+13
      fact=376041246141627.e0+143031288139776.e2*bp2;
      fact=fact-8084779967808.e4*bp4-4296764470272.e6*bp6;
      fact=fact+870387186176.e8*pow(bp,8)-63160877056.e10*pow(bp,10);
      fact=fact+22278144.e14*pow(bp,12)-33816576.e14*pow(bp,14);
      fact=fact+196608.e16*pow(bp,16);
      fact=fact*8.19200E+13/(bp*pow(81. + 400.*bp2,14)) ;
      }
 else if (n_rn == 9) {
//c 4->5
// pp pow(2.,19)*pow(5.,7) = 40960000000
      fact=2217735398973.e0-546176812856.e2*bp2;
      fact=fact+7931167704.e5*bp4-37782745600.e5*bp6;
      fact=fact+97729280000.e5*pow(bp,8)-71884800000.e5*pow(bp,10);
      fact=fact+40960000000.e5*pow(bp,12);
      fact=fact*40960000000 /(bp*pow(81.+400.*bp2,11));
      }
 else if (n_rn == 10) {
//c 4->6
      fact=12881310395.e0+2722585301712.e0*bp2;
      fact=fact-66947504285952.e0*bp4+87094878234624.e1*bp6;
      fact=fact-396354332491776.e1*pow(bp,8)+960059690975232.e1*pow(bp,10);
      fact=fact-7016971052777472.*pow(bp,12)+3851755393646592.*pow(bp,14);
      fact=fact*47775744./(bp*pow(25.+144.*bp2,12));
      }

double FEX3a=fact*corr;
return FEX3a;
}
//c*********************************************************************


