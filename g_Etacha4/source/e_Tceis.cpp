#include "e_Etacha4.h"
#include <complex>
using namespace std;

#include <math.h>
#include <stdio.h>
#include "../win/e_myextern.h"
#include "../win/e_Constant.h"
//#include "e_myextern.h"
#include "e_GaussDataBlock.h"
#define gama 0.577215664901532860606512

//#include <vcl.h>

extern double pow2(double par);
extern double pow3(double par);
extern double pow2I(double par);
extern complex<double> pow2(complex<double> par);

extern double Velocity_au(double E);
extern char BufRtf[1000];


double dcoufa(double x);
double GaussLaguerre(int n, double (*fun)(double x), double alpha, double qmin);
double AVINT(double *x, double *y,int n, double xlo, double xup);
void POCH(int n, complex<double> cx, complex<double>* cpoch);

double GQUAD(double (*Func)(double x),double AX, double BX, int NX, int &N);
double FRZT(double x);
double SRZT(double x);

void hetar  (double x, double al, double ar, double ag, double ad, complex<double> &crcei);
void hy1star(double x, double al, double ar, double ag, double ad, complex<double> &crcei);
void hy2star(double x, double al, double ar, double ag, double ad, complex<double> &crcei);
void hy2ptar(double x,double eta,double al,double ar,double ag,double ad,
                complex<double> &crx,complex<double> &cry,complex<double> &crz);


void hypcei(complex<double> ca, complex<double> cb, complex<double> cc, double x, complex<double> &cf);
void hypere (complex<double> ca, complex<double> cb, complex<double> cc, double x, complex<double> &cf);
complex<double> LogGammaFunc(complex<double> cz);
complex<double> cGamLn(complex<double> cz);
complex<double> diser (complex<double> z);

void ca1536 (complex<double> ca,complex<double> cb, complex<double> cc, double x, complex<double> &cf);
void ca1537 (complex<double> ca,complex<double> cb, complex<double> cc, double x, complex<double> &cf);
void ca1538 (complex<double> ca,complex<double> cb, complex<double> cc, double x, complex<double> &cf);
void ca15310(complex<double> ca,complex<double> cb,                     double x, complex<double> &cf);
void ca15313(complex<double> ca,                    complex<double> cc, double x, complex<double> &cf);

void chypser(complex<double> ca,complex<double> cb, complex<double> cc, double x, complex<double> &cf);

void digi1(int nz, double xz, double yz, complex<double>&  diser);
void digi2(complex<double> z,            complex<double>&  diser);
void digam(complex<double> z,            complex<double>&  psi);

void POCH(int n, complex<double> cx, complex<double>& cpoch);


void error(int kod);


//********************
//*Pour ETACHA4
//*  Total cross sections for single ionization of Hydrogenic or Helium
//*  targets by bare ion impact using the CDW-EIS (RHF) approximation.
//*  (see Fainstein \etal, JPB 24 3091 (1991) and references therein)
//*
//*  Uses: AVINT : Integration Program - 2nd order polinomial fit
//*        HYPER : Gauss Hypergeometric function
//*        GaussLaguerre : Integration Program - Gauss-Laguerre
//*        GQUAD : Integration Program - Gauss
//*
//*  Input Files  : gcdw.dat -->  NG   : no. of points for Gauss-Laguerre
//*                                      integration
//*                               ALPHA: scaling factor for Gauss-Laguerre
//*                                      integration
//*                               ND   : no. of points for Gauss integration
//*                               NA   : no. of points for angular integration
//*                               DTE  : electron angle (deg)
//*
//*
//*  Author: PabloF  (pablof@cab.cnea.gov.ar)     Last Version: 14/06/2000
//*          Hydrogenic 2s & 2p initial states by Mariel
//********************

// common/pa4/ca0,ca1,cb0,cb1;
const complex<double> c0(0,0), c1(1,0),c2(2,0), ci(0,1);
complex<double> ca0,ca1, cb0,cb1;


double zef, rnu, sc, pv, ek, de, fi, qmin, coa, epLoc, ze, fip, alpha_scale, te;
double zp, zt;

//      common/pa3/ninis,ng;
int ninis, ng;

//      common/pa5/icor;
int icor;


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ETACHA::Tceis(double Ecurr, double &tot, int netat)
{
extern double Zp,Zt;

        zt=Zp;   // INVERSE!!!!
        zp=Zt;   // INVERSE!!!!


#define a0 5.2917706e-9
#define au 27.2
#define enstar  -2.0
#define enstep  0.125


#define nang  50
#define nangc 51
#define nele  100
#define nelec 101
#define na    18
#define nac   19


      double eand[nangc],eanr[nangc],eene[nelec],table[nangc];
      double ddsia[nangc],ddsie[nelec];
      double sdene[nelec],sdang[nangc],dd100[nelec+1][nangc+1];

      int npedo;
      double epamu=1000.*Ecurr;
      ninis=netat;

      int jconv=0;
      ng=20;
      alpha_scale=5.;
      int nd=15;

      double tableLoc[nac] = {0,
          0.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,
           60.0,70.0,80.0,90.0,100.0,110.0,130.0,150.0,170.0,180.0 };

      for(int i=0; i<=na; i++) { table[i]=eand[i]=tableLoc[i];}


     double ga2=(1.+epamu/931.5e3)*(1.+epamu/931.5e3);
     pv=1.3703599976e2*sqrt(1.-1./ga2);
// temporary  Oleg      pv = Velocity_au(Ecurr);

//  commented originally    pv=dsqrt(epamu/24.98d0)
      rnu=zp/pv;
      double ei=1, dd0=0;
      double dte, d100;

      zef = zt;

      if(ninis == 10)
         {
          ei=-zt*zt/2.;
          dd0=exp(-2.*PI*rnu)*dcoufa(rnu)*pow2((zp/pv/PI))/2.;
          dd0=pow(zt,5)*dd0;
          //printf("Hydrogenic Target - 1s Initial State");
        }
      else if(ninis == 20) {
          ei=-zt*zt/8.;

          dd0=exp(-2.*PI*rnu)*dcoufa(rnu)*pow2((zp/pv/PI/4.));
          dd0=pow(zt,5)*dd0;
         // printf("Hydrogenic Target - 2s Initial State");
        }
      else if(ninis == 21) {
          ei=-zt*zt/8.;

          dd0=exp(-2.*PI*rnu)*dcoufa(rnu)*pow2((zt*zp/pv/PI))/6.;
          dd0=pow(zt/2.,5)*dd0;
        //  printf("Hydrogenic Target - 2p Initial State");
        }


      ca0=ci*rnu;
      ca1=ci*rnu+c1;

      double scf= 1.13*pow(zp,(1./3.));
      double fs1=1.0;

            if (zp == 1.) fs1=1.23;
       else if (zp == 2.) fs1=1.526;
       else if (zp == 6.) fs1=0.78;
       else if (zp == 7.) fs1=0.85;
       else if (zp == 10.) fs1=1.04;
       else if (zp == 13.) fs1=0.58;
       else if (zp == 14.) fs1=0.59;
       else if (zp == 18.) fs1=0.68;
       else if (zp == 29.) fs1=0.672;
       else if (zp == 36.) fs1=0.61;
       else if (zp == 54.) fs1=0.535;

       sc = scf*fs1;
       icor=1;


//*****Comienza el calculo de la DDCS.*******************************************
      int i=1;
      int nf;
      double enelog = enstar;
     const double sdene_precision = 1.e-08;
//     const double sdene_precision = 1.e-03;

//L310:
   //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
      sprintf(BufRtf," i &nbsp;  EEV[eV]  &nbsp;  SDCS[cm<sup>2</sup>/eV]"); emit appendShell(0,0);

while(1)
      {
      eene[i]=pow(10.,enelog);

      ek   = sqrt(2.*eene[i]/au);
      de   = 0.5*ek*ek-ei;
      fi   = zef/ek;
      qmin = de/pv;
      //  is not used ==> double qmm  = -ei/pv;

      if (pv > 1.75*qmin) coa=1.0-pow2(1.75*qmin/pv);
      else                coa=0.0;

      //--------------------------------------
      for(int j=1; j<=na; j++)
        {
         dte = eand[j];

              if(dte == 0. ) te=0.;
         else if(dte == 180) te=PI;
         else                te=dte*PI/180.;

        eanr[j] = te;

        epLoc=sqrt(ek*ek-2.*ek*pv*cos(te)+pv*pv);
        ze=zp/epLoc;
        cb0=ci*ze;
        cb1=ci*ze+c1;

         if(dte == 0. || dte == 180.)
                {
                fip=0.;
                d100=2.*PI*GaussLaguerre(ng,FRZT,alpha_scale,qmin);       //  it passed for ninis==21 !
                }
        else    {
                d100=2.*GQUAD(SRZT,0.,PI, nd, npedo);
                }

        double dd1=dcoufa(fi)*dcoufa(ze);
        dd100[i][j]=ek*d100*dd0*dd1*pow2(a0)/au ;
        ddsia[j]=sin(eanr[j])*dd100[i][j];
        }


      nf=i;
      dd100[i][na+1]=2.*PI*AVINT(eanr,ddsia,na,eanr[1],eanr[na]);
      sdene[i] = dd100[i][na+1];
      sprintf(BufRtf,"%02d  &nbsp; %10.3e &nbsp;  %10.3e",i,eene[i],sdene[i]); emit appendShell(0,0);

      if(sdene[i] < sdene_precision*sdene[1]) break; //goto L350;
// test      if(i==2) break; //goto L350;

      enelog+=enstep;
      i++;

      if(i > nele) {error(1); error(2); return;}
      }
   //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
   //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

          //goto L310;  changed for while
//      jconv=1;         // unreachable code

      //----------------------------------------------

//L350:
      for(int j=1; j<=na; j++)
        {
        for(int i=1; i<=nf; i++) ddsie[i]=dd100[i][j];

        dd100[nf+1][j] = AVINT(eene,ddsie,nf,eene[1],eene[nf]);
        sdang[j]       = sin(eanr[j])*dd100[nf+1][j];
        }

//*****Comienza el calculo de la TCS.********************************************
      double tcene =       AVINT(eene,sdene,nf,eene[1],eene[nf]);
      double tcang = 2.*PI*AVINT(eanr,sdang,na,eanr[1],eanr[na]);
      double tcs   = 0.5*(tcene+tcang);
//     tcs=tcene
      if(zp < 1.) tcs /=  pow2(zp);


      if(icor == 1) strcpy(BufRtf,"SC and ASC corrected<BR> CDW-EIS cross sections");
      else          strcpy(BufRtf,"uncorrected CDW-EIS cross sections");

      emit appendShell(0,0);

      sprintf(BufRtf,"tcene (DTE &rarr; EEV)   &nbsp;    %10.3e<br>"
                     "tcang (EEV &rarr; DTE)   &nbsp;    %10.3e<br>"
                     "tcs=0.5*(tcene+tcang) &nbsp; %10.3e cm<sup>2</sup><br>"
                     "EPAMU=%9.3e keV/amu &nbsp; Zp=%4.1f Zt=%4.1f<br>", tcene,tcang,tcs,epamu,zt,zp);
      emit appendText(0,0);


if(jconv != 0) { error(1); error(3);}

tot=tcs;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double SRZT(double x)
{
      fip=x;
      double srzt=GaussLaguerre(ng,FRZT,alpha_scale,qmin);
      return srzt;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

double FRZT(double x)
{

     complex<double> crcei,crx,cry,crz;

     double FRZT=0;

      double eta = sqrt (x*x-qmin*qmin);
      double al =  0.5*x*x;
      double ar = -ek*sin(te)*eta*cos(fip)-ek*cos(te)*qmin;
      double ag =  al+ar+de;
      double ad = -de+pv*(ek*cos(te)-pv-epLoc);

//commented originanlly       nt=dint(zt)
///commented originanlly      if(zt/nt.eq.1.) then

         if(ninis == 10) {
                  hy1star(x,al,ar,ag,ad,crcei);
                  FRZT=x*pow2(abs(crcei));                      //  Oleg  check result!!!!
                  }
    else if(ninis == 20) {
                  hy2star(x,al,ar,ag,ad,crcei);
                  FRZT=x*pow2(abs(crcei));                      //  Oleg  check result!!!!
                  }
     else if(ninis == 21) {
                  hy2ptar(x,eta,al,ar,ag,ad,crx,cry,crz);
                  FRZT= x*(pow2(abs(crx)) + pow2(abs(cry))+ pow2(abs(crz)));
                  }
//      else      {
//                hetar(x,al,ar,ag,ad,crcei);                     // Oleg check it for another  ninis
//                FRZT=x*pow2(abs(crcei));
//                }

//     SC and ASC corrections
       if (icor == 1) {
              double sq=pow2I(1.+pow2(sc/x))+coa*(1.-pow2I(1.+pow2(x/sc)))/zp;
//     write (6,*) 'sq=', sq
      FRZT *=sq;
      }

return FRZT;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void hy1star(double /*x*/, double al, double ar, double ag, double ad, complex<double> &crcei)
{
//
// R(\eta) for Hydrogenic Target 1s Initial State

complex<double> cf1,cf2,cgg,cjj,cbc,ccc,cgc,crbn1,cauxs;

      double zc=1.+al*ad/de/ag;
      hypcei(ca0,cb0,c1,zc,cf1);
      hypcei(ca1,cb1,c2,zc,cf2);

      cgg = (ar+ek*ek*(c1+ci*fi))/ag;
      cjj = c1-cgg;
      cbc = 2.*al+ar*(c1+ci*fi);
      ccc = pv/epLoc*ag*cgg - (1.0+pv/epLoc)*(-de+pv*ek*cos(te)*(c1+ci*fi));

      cgc = -al*(ad*cbc+ag*ccc)/de/ag/cbc;

      crbn1 = exp(-ci*fi*log(cjj)) * cbc/cjj/al/pow3(ag);   //  Oleg check it for complex functions
      cauxs = cf1-ci*rnu*cgc*cf2;
      crcei = crbn1*cauxs;
      }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void hy2star(double x, double al, double ar, double ag, double ad, complex<double> &crcei)
{

// R(\eta) for Hydrogenic Target 2s Initial State
//          Author: Mariel                     Last Version: 12/06/2000;

      complex<double> cf1, cf2;

      double zc = 1.0+al*ad/de/ag;
      hypcei(ca0,cb0,c1,zc,cf1);
      hypcei(ca1,cb1,c2,zc,cf2);

      double bet   =  zt/2.;
      double fibet = bet/ek;
      double gam   = 0.5*(bet*bet+x*x+2.*ar+ek*ek);

      complex<double>  cgg = (ar+ek*ek*(c1+ci*fibet))/gam;
      complex<double>  cjc =c1-cgg;
      complex<double>  cbc =x*x+ar*(c1+ci*fibet);
      complex<double>  ccc = pv/epLoc*gam*cgg  - (1.+pv/epLoc)*(-de+pv*ek*cos(te)*(c1+ci*fibet));

      complex<double> cgc   = -al*(ad*cbc+ag*ccc)/de/ag/cbc;
      complex<double> cr1   = cbc*exp(-ci*fi*log(cjc))/cjc/pow2(gam);
      complex<double> crh1s = (cr1*cf1-ci*rnu*cr1*cgc*cf2)/(al*ag);

//***********************Estado de Slater 2s************************

      complex<double> cor   = -2.*bet*cjc +(c1+ci*fi)*(ci*ek-bet*cgg);
      complex<double> caux1 = cbc*cor+ci*cjc*gam*ar/ek;
      complex<double> cfac  = -(al/de/ag)*(ad*ar/ek+ag*pv/epLoc*(ek-pv*cos(te)-epLoc*cos(te)));
      complex<double> caux2 = cbc*cor*cgc+ci*cjc*gam*cfac;
      complex<double> cr2   = exp(-ci*fi*log(cjc))/pow2(cjc)/pow3(gam);
      complex<double> crh2s = cr2*(caux1*cf1-ci*rnu*caux2*cf2)/al/ag;

//******************************************************************

      crcei=crh1s+bet*crh2s;
      }

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void hy2ptar(double x,double eta,double al,double ar,double ag,double ad,
                complex<double> &crx, complex<double> &cry, complex<double> &crz)
                {

// R(\eta) for Hydrogenic Target 2p Initial State
//   Author: Mariel                     Last Version: 12/06/2000

      double mek[4],mpv[4],meta[4],mx[4],mep[4];
      complex<double> cf1, cf2, cr2p[4];

      double zc = 1. + al*ad/de/ag;
      hypcei(ca0,cb0,c1,zc,cf1);
      hypcei(ca1,cb1,c2,zc,cf2);

      mek[1]=ek*sin(te);
      mek[2]=0.;
      mek[3]=ek*cos(te);
      mpv[1]=0.;
      mpv[2]=0.;
      mpv[3]=pv;
      meta[1]=eta*cos(fip);
      meta[2]=eta*sin(fip);
      meta[3]=0.;
      mep[1]=mek[1]-mpv[1];
      mep[2]=mek[2]-mpv[2];
      mep[3]=mek[3]-mpv[3];

      double bet   = zt/2.;
      double fibet = bet/ek;
      double gam   = 0.5*(bet*bet+x*x+2.*ar+ek*ek);

      complex<double> cgg   = (ar+ek*ek*(c1+ci*fibet))/gam;
      complex<double> cjc   = c1-cgg;
      complex<double> cbc   = x*x+ar*(c1+ci*fibet);
      complex<double> cm2p  = ar+ek*ek*(c1+ci*fibet) +(1.+epLoc/pv)*(de-pv*ek*cos(te)*(c1+ci*fibet));
      complex<double> crsl1 = exp(-ci*fi*log(cjc))/pow2(cjc)/pow(gam,3);

      for(int j=1; j<=3; j++)
        {
        mx[j] = -meta[j]-de*mpv[j]/(pv*pv);
        complex<double> cor2p =(c1-ci*fi)*(mx[j]+mek[j])*cjc+(c1+ci*fi)*mx[j];
        complex<double> cu2p  =-cbc*cor2p+gam*cjc*mx[j];
        complex<double> cl2p  =-cm2p*cor2p+gam*cjc*(mep[j]-epLoc/pv*mpv[j]);
        complex<double> caux1 =cu2p*(cf1+ci*rnu*(al*ad/de/ag)*cf2);
        complex<double> caux2 =ci*rnu*(al*pv/de/epLoc)*cl2p*cf2;

        cr2p[j] = -crsl1*(caux1+caux2)/(al*ag);
        }

      crx=cr2p[1];
      cry=cr2p[2];
      crz=cr2p[3];

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void hetar(double /*x*/, double al, double ar, double ag, double ad, complex<double> &crcei)
{
//
// R(\eta) for Hydrogenic Target 1s Initial State
//

double d100[6] = {0,
               2.5925,1.6377,0.75254,-0.3315,0.103};

double b100[6]= {0,
               1.41714,2.37682,4.39628,6.52699,7.94252};
      complex<double> cf1, cf2, cs1, cs2, cgg, cjc, cbc, ccc, cgc, cr1, cr2;

      double zc=1.+al*ad/de/ag;
      hypcei(ca0,cb0,c1,zc,cf1);
      hypcei(ca1,cb1,c2,zc,cf2);

      cs1=c0;
      cs2=c0;

      for(int j=1; j<=5; j++)
        {
        double gam = 0.5*(pow2(b100[j])+2.*al+ek*ek+2.*ar);
        double fij = b100[j]/ek;

        cgg=(ar+ek*ek*(c1+ci*fij))/gam;
        cjc=c1-cgg;
        cbc=2.*al+ar*(c1+ci*fij);
        ccc=pv/epLoc*gam*cgg  -(1.+pv/epLoc)*(-de+pv*ek*cos(te)*(c1+ci*fij));
        cgc=-al*(ad*cbc+ag*ccc)/de/ag/cbc;

        cr1=d100[j]*cbc*exp(-ci*fi*log(cjc))/cjc/pow2(gam);
        cr2=cr1*cgc;

        cs1=cs1+cr1;
        cs2=cs2+cr2;
        }

  crcei=(cs1*cf1-ci*rnu*cs2*cf2)/(al*ag);

}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double dcoufa(double x)
{
double t = 2.*PI*x;
double  dcoufa =  t/(1.-exp(-t));
return dcoufa;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void error(int kod)
{
      if(kod == 1) printf("I am sorry but this calculation did not converge\n"
                          "Please address your complain to pablof@cab.cnea.gov.ar\n");
 else if(kod == 2) printf("Error Code Number 1: vector dimension\n");
 else if(kod == 3) printf("Error Code Number 2: convergence\n");
 else if(kod == 4) printf("This initial state is not yet available\n");
 else if(kod == 5) printf("l > n-1 , please read a book on Quantum Mechanics !\n");
}
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

double GQUAD(double (*Func)(double x), double AX, double BX, int NX, int &N)
{
//     N-POINT GAUSSIAN QUADRATURE OF FUNCTION F OVER INTERVAL (AX,BX).

//-----TEST N
      N=NX;
      double alpha = 0.5*(BX+AX);
      double beta  = 0.5*(BX-AX);
      N = max(1, N);
      N = min(96,N);

      if(N == 1) {return (BX-AX)*Func(alpha); }

       if(N > 16) N = 4*(N/4);
       if(N > 24) N = 8*(N/8);
       if(N > 48) N =16*(N/16);

//----- SET K EQUAL TO INITIAL SUBSCRIPT AND INTEGRATE

      int K=gKTAB[N];
      int M=N/2;
      double SUM=0.;
      int JMAX=K-1+M;

      for(int J=K; J<=JMAX; J++)
              {
              double delta = beta*gaussX[J];
              SUM=SUM+gaussA[J]*(Func(alpha+delta)+Func(alpha-delta));
              }

      if(N - 2*M != 0)
                SUM+=gaussA[K+M]*Func(alpha);


      double GQUAD=beta*SUM;
      return GQUAD;
      }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void hypcei(complex<double> ca, complex<double> cb, complex<double> cc, double x, complex<double> &cf)
{
//  Program: hypere.f
//  Version: 1.6 (see readme_hyper) 24/03/1997
//  Author: Pablof (pablof@cab.cnea.gov.ar)
//  Hypergeometric function for CDW-EIS calculations.
//  Double precision, see comment for V 1.4
//  We assume that x is REAL

       if(ca == c0 || cb == c0 )cf=c1;
 else  if(x == 0)               cf=c1;
 else  if(x == 1.)
                                cf= exp( LogGammaFunc(cc)
                                        -LogGammaFunc(cc-ca-cb)
                                        -LogGammaFunc(cc-ca)
                                        -LogGammaFunc(cc-cb)
                                        );

 else               hypere(ca,cb,cc,x,cf);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void hypere(complex<double> ca, complex<double> cb, complex<double> cc, double x, complex<double> &cf)
{
#define step1 -200.
#define step2 -0.618
#define step3  0.5
#define step4  1.1
double z;

if(x < step1)
        {
        if(ca == cb)
                {
                z = x/(x-1.);
                ca15310(cc-ca,cb,z,cf);
                cf=exp(-cb*log(1.-x))*cf;
                }
       else             ca1538(ca,cb,cc,x,cf);
      }
 else if(x < step2)
        {
        if(ca == cb)
                {
                z=x/(x-1.);
                chypser(cc-ca,cb,cc,z,cf);
                cf=exp(-cb*log(1.-x))*cf;
                }
        else            ca1538 (ca,cb,cc,x,cf);
       }
 else if(x < step3)     chypser(ca,cb,cc,x,cf);
 else if(x < step4)
        {
        if(cc == ca+cb) ca15310(ca,cb,   x,cf);
        else            ca1536 (ca,cb,cc,x,cf);
        }
 else   {
        if(ca == cb)    ca15313(ca,   cc,x,cf);
        else            ca1537 (ca,cb,cc,x,cf);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ca1536 (complex<double> ca,complex<double> cb, complex<double> cc, double x, complex<double> &cf)
{
// formula 15.3.6 from M. Abramowitz and I. Stegun

double z;
complex<double> cz,cca,ccb,ccab,cgam1,cf1,cgam2,cf2;

cca  = cc-ca;
ccb  = cc-cb;
ccab = cc-ca-cb;

if(ccab == c0) printf("There is a pole cc=ca+cb");

z=1.-x;
cz=complex<double>(z,0.);

cgam1 = exp( LogGammaFunc(cc) + LogGammaFunc( ccab) -LogGammaFunc(cca)-LogGammaFunc(ccb));
cgam2 = exp( LogGammaFunc(cc) + LogGammaFunc(-ccab) -LogGammaFunc( ca)-LogGammaFunc(cb ))*exp(ccab*log(cz));

chypser(ca, cb, c1-ccab,z,cf1);
chypser(cca,ccb,c1+ccab,z,cf2);

cf=cgam1*cf1+cgam2*cf2;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ca1537(complex<double> ca,complex<double> cb, complex<double> cc, double x, complex<double> &cf)
{
//          s formula 15.3.7 from M. Abramowitz and I. Stegun

double z;
complex<double> cca,ccb,cba,cx,cgam1,cf1,cgam2,cf2;

      cca=cc-ca;
      ccb=cc-cb;
      cba=cb-ca;
      if(cba == c0) printf("There is a pole ca=cb");

      cx=complex<double>(x,0.);
      cgam1=exp(LogGammaFunc(cc)+LogGammaFunc( cba)-LogGammaFunc(cb)-LogGammaFunc(cca))* exp(-ca*log(-cx));
      cgam2=exp(LogGammaFunc(cc)+LogGammaFunc(-cba)-LogGammaFunc(ca)-LogGammaFunc(ccb))* exp(-cb*log(-cx));
      z=1./x;
      chypser(ca,c1-cca,c1-cba,z,cf1);
      chypser(cb,c1-ccb,c1+cba,z,cf2);
      cf=cgam1*cf1+cgam2*cf2;
      }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ca1538(
             complex<double> ca,
             complex<double> cb,
             complex<double> cc,
             double x,
             complex<double> &cf)
{
//formula 15.3.8 from M. Abramowitz and I. Stegun

double z;
complex<double> cca,ccb,cba,cgam1,cf1,cgam2,cf2;



      cca=cc-ca;
      ccb=cc-cb;
      cba=cb-ca;

      if(cba == c0) printf("There is a pole ca=cb");

      cgam1=exp(LogGammaFunc(cc)+LogGammaFunc( cba)-LogGammaFunc(cb)-LogGammaFunc(cca))*  exp(-ca*log(1.-x));
      cgam2=exp(LogGammaFunc(cc)+LogGammaFunc(-cba)-LogGammaFunc(ca)-LogGammaFunc(ccb))*  exp(-cb*log(1.-x));
      z=1./(1.-x);
      chypser(ca,ccb,c1-cba,z,cf1);
      chypser(cb,cca,c1+cba,z,cf2);
      cf=cgam1*cf1+cgam2*cf2;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ca15310(complex<double> ca,
             complex<double> cb,
             double x,
             complex<double> &cf)
{

//  Uses formula 15.3.10 from M. Abramowitz and I. Stegun
complex<double> caa,cbb,cc,psi1,psi2,psi3,coef,csum,cfac,ctemp;

      int n;

      n=1;
      caa=ca;
      cbb=cb;
      cc=caa+cbb;

      digam(c1,psi1);
      digam(caa,psi2);
      digam(cbb,psi3);

      complex<double> cmx(1.-x,0.);
      coef =-log(cmx);
      csum =coef+2.*psi1-psi2-psi3;
      cfac =c1;
      ctemp =cfac*csum;

      psi1=psi1+1.;
      psi2=psi2+1./caa;
      psi3=psi3+1./cbb;

L100:
        cfac=cfac*caa*cbb*cmx/double(n*n);
        csum=coef+2.*psi1-psi2-psi3;
        cf=ctemp+cfac*csum;
         if(cf == ctemp) goto L200;
        ctemp=cf;
        n=n+1;
        caa=caa+c1;
        cbb=cbb+c1;
        psi1=psi1+1./double(n);
        psi2=psi2+1./caa;
        psi3=psi3+1./cbb;
      goto L100;

L200:
      cf=cf*exp(LogGammaFunc(cc)-LogGammaFunc(ca)-LogGammaFunc(cb));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ca15313(complex<double> ca,
             complex<double> cc,
             double x,
             complex<double> &cf)
{
// formula 15.3.13 from M. Abramowitz and I. Stegun

complex<double> caa,cbb,psi1,psi2,psi3,cmx,coef,csum,cfac,ctemp;
int n;

      n=1;
      caa=ca;
      cbb=c1+ca-cc;
      digam(c1,psi1);
      digam(caa,psi2);
      digam(cbb,psi3);

      cmx  = complex<double>(-x,0.);
      coef = log(cmx)-PI*cos(PI*cbb)/sin(PI*cbb);
      csum = coef+2.*psi1-psi2-psi3;
      cfac = c1;
      ctemp= cfac*csum;

      psi1=psi1+1.;
      psi2=psi2+1./caa;
      psi3=psi3+1./cbb;

L100:
        cfac = cfac*caa*cbb/(n*n*x);
        csum = coef+2.*psi1-psi2-psi3;
        cf   = ctemp+cfac*csum;
        if(cf == ctemp) goto L200;
        ctemp=cf;
        n++;
        caa=caa+c1;
        cbb=cbb+c1;
        psi1=psi1+1./double(n);
        psi2=psi2+1./caa;
        psi3=psi3+1./cbb;
        goto L100;
L200:

      cf = cf * exp(LogGammaFunc(cc)-LogGammaFunc(ca)-LogGammaFunc(cc-ca)) * exp(-ca*log(cmx));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void chypser(complex<double> ca,
             complex<double> cb,
             complex<double> cc,
             double x,
             complex<double> &cf)
{
//  Hypergeometric series calculated to machine accuracy
complex<double> cfac,ctemp,caa,cbb,ccc;
int n;

       if(ca == c0 || cb == c0) {  cf=c1;  return; }

      n=1;
      caa=ca;
      cbb=cb;
      ccc=cc;
      cfac=c1;
      ctemp=cfac;

while(1)
        {
        cfac=((caa*cbb)/ccc)*cfac;
        cfac=cfac*x/double(n);
        cf=ctemp+cfac;
        if(cf == ctemp) break;

        ctemp=cf;
        n=n+1;
        caa=caa+c1;
        cbb=cbb+c1;
        ccc=ccc+c1;
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void digam(complex<double> z, complex<double>& psi)
{
/*  Subroutine to calculate the DIGAMMA function of a complex
*  argument using formulas from M. Abramowitz and I. Stegun.
*
*  If z = integer --> uses A6.3.2  --> gives very good results
*  If z = complex --> uses A6.3.16 --> not so good ...
*
*  If Real(z).lt.0.d0 it uses the reflection formula A6.3.7
*  to avoid the poles in the series expansion.
*
*  Author: Pablof    V 0.8  15/12/1994 */


if(z == c1)  {  psi = -c1*gama;  return;}


if(z.real() < 0)  psi = diser(c1-z) - PI*cos(PI*z)/sin(PI*z);
else              psi = diser(z);
}



/*
 *
subroutine digamma_complex(z, result)
    ! Compute the digamma (psi) function of a complex argument z
    ! Inputs:
    !   z      - Complex argument (input)
    ! Outputs:
    !   result - Complex digamma function value (output)

    implicit none
    complex(kind=8), intent(in) :: z      ! Input complex number
    complex(kind=8), intent(out) :: result ! Output digamma value

    ! Constants
    complex(kind=8) :: z0, psi_series, z_reduced, term
    integer :: n, max_iter
    real(kind=8), parameter :: epsilon = 1.0d-12, pi = 3.141592653589793d0

    ! Initialize variables
    z0 = z
    result = (0.0d0, 0.0d0)
    psi_series = (0.0d0, 0.0d0)
    max_iter = 100

    ! Shift z into the right half-plane using the recurrence relation
    ! psi(z) = psi(z + 1) - 1/z
    do while (real(z0) < 6.0d0)
        result = result - 1.0d0 / z0
        z0 = z0 + 1.0d0
    end do

    ! Use the asymptotic expansion for large |z|:
    ! psi(z) ~ ln(z) - 1/(2z) - sum_{n=1}^(inf) B_{2n} / (2n z^(2n))
    z_reduced = z0
    result = result + log(z_reduced) - 0.5d0 / z_reduced
    do n = 1, max_iter
        term = (-1.0d0)**n * bernoulli(2*n) / (2.0d0 * n * z_reduced**(2*n))
        if (abs(term) < epsilon) exit
        result = result + term
    end do

    ! Handle the imaginary part for the digamma of a complex number
    if (aimag(z) /= 0.0d0) then
        result = result - (0.0d0, pi) * (z - dble(int(real(z))))
    end if
end subroutine digamma_complex

! Function to calculate Bernoulli numbers B_{2n}
real(kind=8) function bernoulli(n)
    integer, intent(in) :: n
    integer :: i
    real(kind=8) :: B
    real(kind=8), dimension(20) :: B_table

    ! Precomputed Bernoulli numbers (only B_{2n} are non-zero)
    B_table = (/ 1.0d0, 1.0d0/6.0d0, -1.0d0/30.0d0, 1.0d0/42.0d0, -1.0d0/30.0d0, &
                5.0d0/66.0d0, -691.0d0/2730.0d0, 7.0d0/6.0d0, -3617.0d0/510.0d0, &
                43867.0d0/798.0d0, -174611.0d0/330.0d0, 854513.0d0/138.0d0, &
                -236364091.0d0/2730.0d0, 8553103.0d0/6.0d0 /)
    if (n >= 1 .and. n <= 20) then
        B = B_table(n)
    else
        B = 0.0d0
    end if
    bernoulli = B
end function bernoulli


********************************************************************************************
#include <complex>
#include <cmath>
#include <iostream>

// Constants
const double PI = 3.141592653589793;
const double EPSILON = 1e-12;
const int MAX_ITER = 100;

// Precomputed Bernoulli numbers B_{2n}
const double BERN[] = {1.0, -1.0 / 6.0, 1.0 / 30.0, -1.0 / 42.0, 1.0 / 30.0, -5.0 / 66.0,
                       691.0 / 2730.0, -7.0 / 6.0, 3617.0 / 510.0, -43867.0 / 798.0};

std::complex<double> bernoulli(int n) {
    if (n < 1 || n > 10) return 0.0;
    return BERN[n - 1];
}

// Digamma function for a complex argument
std::complex<double> digamma(const std::complex<double>& z) {
    std::complex<double> result = 0.0;
    std::complex<double> z_shifted = z;

    // Shift z into the right half-plane if necessary
    while (real(z_shifted) < 6.0) {
        result -= 1.0 / z_shifted;
        z_shifted += 1.0;
    }

    // Use the asymptotic expansion for large |z|
    std::complex<double> term;
    result += std::log(z_shifted) - 0.5 / z_shifted;
    for (int n = 1; n < MAX_ITER; ++n) {
        term = bernoulli(2 * n) / (2.0 * n * std::pow(z_shifted, 2 * n));
        if (std::abs(term) < EPSILON) break;
        result -= term;
    }

    // Handle the imaginary part for non-real arguments
    if (imag(z) != 0.0) {
        result -= std::complex<double>(0.0, PI) * std::floor(real(z));
    }

    return result;
}

// Example usage
int main() {
    std::complex<double> z(3.0, 2.0); // Example complex number
    std::complex<double> result = digamma(z);
    std::cout << "Digamma(" << z << ") = " << result << std::endl;
    return 0;
}

 */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
complex<double> diser(complex<double> z)
{
complex<double> cdiser;

      double xz,yz,divx;
      int nz;

      xz=z.real();
      yz=z.imag();
      nz=xz;

       if(xz != 0.)  divx=nz/xz;

       if(xz == 0. || divx == 1.)
                {
                if(yz == 0.)
                        {
                        cdiser=-gama;
                        for(int k=1; k<=nz-1; k++)
                                        cdiser = cdiser+1./k;
                        }
                 else   {
                        digi1(nz,xz,yz,cdiser);
                        }
                }
        else digi2(z,cdiser);

return cdiser;
}

/*
 * The function ∑k=1N1k∑k=1N​k1​ is the NN-th harmonic number, denoted by HNHN​. An asymptotic approximation for HNHN​ as N→∞N→∞ is:
HN∼ln⁡(N)+γ+12N−112N2+O(1N4),
HN​∼ln(N)+γ+2N1​−12N21​+O(N41​),

where:

    ln⁡(N)ln(N) is the natural logarithm of NN,
    γγ is the Euler-Mascheroni constant (γ≈0.577γ≈0.577).

For a simpler leading-order asymptotic, you can use:
HN∼ln⁡(N)+γ.
HN​∼ln(N)+γ.

This approximation is sufficient for most applications. The higher-order terms become significant only for small NN.
*/

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void digi1(int nz, double xz, double yz,complex<double>&  diser)
{
      complex<double> z;
      double fac,fac1,fac2,temp1,temp2,rpsi,cpsi;


      fac=1./(1.+yz*yz);
      temp1 = fac;
      temp2 = fac;

      for(int n=2; n<=10000; n++)
              {
              fac2=1./(n*n+yz*yz);
              fac1=fac2/double(n);
              rpsi=temp1+fac1;
              cpsi=temp2+fac2;
              if(rpsi == temp1 && cpsi == temp2) break;
              temp1=rpsi;
              temp2=cpsi;
              }

//     pause 'convergence failure in digi1'

      rpsi = -gama+yz*yz*rpsi;
      cpsi = yz*cpsi;

      if(xz == 0.)  { diser=complex<double>(rpsi,cpsi+1./yz);}
 else if(xz == 1.) {  diser=complex<double>(rpsi,cpsi);      }
 else {
        diser=complex<double>(rpsi,cpsi);

        for(int n=1; n<=nz-1; n++)
                        {
                        z=complex<double>(double(n),yz);
                        diser=diser+1./z;
                        }
      }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void digi2(complex<double> z, complex<double>&  diser)
{
complex<double> cfac,ctemp;
double rn;

cfac  = -gama-1./z + PI*PI/6.0;
ctemp =  cfac;

for(int n=1; n<=10000; n++)
        {
        rn    = n;
        cfac  = z/(rn*rn+rn*z)-1./(rn*rn);
        diser = ctemp+cfac;
        if(diser == ctemp) break;
        ctemp=diser;
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//**************************************
//  Subroutine tu calculate the logarithm of the Gamma function
//  It uses the program gammln from Numerical Recipes (2nd Ed.)
//  which is valid only when Real(cz).gt.0. Thus, the present
//  subroutine introduces appropiate formulas to calculate for
//  any value of cz.
//
//  Author: Pablof    V 1.0   9/12/1994
//**************************************
complex<double> LogGammaFunc(complex<double> cz)
{
complex<double> cz1,canswer;
const char *form = "  z = %e14.7 %e14.7i is considered as a pole";

double xz = cz.real();
double yz = cz.imag();

if(cz == c0) { printf(form, xz,yz);  return c0;}

if(xz == 0.)
        {
        cz1=c1+cz;
        canswer=cGamLn(cz1)-log(cz);
        }
 else if(xz < 0)
        {
        if(int(xz)/fabs(xz) == -1. && yz == 0.)
                {
                 printf(form,xz,yz);
                 return c0;
                }
        cz1=c1-cz;
        canswer=log(PI/sin(PI*cz))-cGamLn(cz1);
       }
 else  canswer=cGamLn(cz);

return canswer;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
complex<double>  cGamLn(complex<double> cz)
{
      double rcof[7] = {0,
                        76.18009172947146e0,-86.50532032941677e0,
                        24.01409824083091e0,-1.231739572450155e0,
                        0.1208650973866179e-2,-0.5395239384953e-5};

      double stp=2.5066282746310005;

      if(cz.real() <= 0.) return c0; //// 'real[z] <= 0';

      complex<double> cx(cz);
      complex<double> cy(cx);
      complex<double> ctmp(cx+5.5);
      complex<double> cser(1.000000000190015,0);

      ctmp=(cx+0.5)*log(ctmp)-ctmp;

      for(int j=1; j<=6; j++)
                {
                cy   = cy+1.;
                cser = cser+rcof[j]/cy;
                }

complex<double> canswer=ctmp+log(stp*cser/cx);

return canswer;
}
/*
c**********************************************************************
*
*    PabloF's library
*    Contents:           AVINT.f -----> Simpson
*                        galag.f -----> Gauss-Laguerre
*                        POCH.f ------> Pochammer symbol
*
*-----------------------------------------------------------------------------*
*                                                                             *
*    This program computes the integral of y(x) in XLO < x < XUP .            *
*                                                                             *
*    X and Y are N dimensional arrays containing the abscissas x and the      *
*    ordinates y respectively.The x are usually unequally spaced and must     *
*    be given in ascending order.The y are usually data.XLO must be less      *
*    than XUP,but otherwise there are no restrictions on XLO and XUP except,  *
*    of course that of accuracy which requires that if XLO is less than x1,   *
*    it be close to x1,and if XUP is greater than xN,it be close to xN.       *
*                                                                             *
*    AVINT is adapted from P.E.Hennion,Algorithm 77,CACM 5,1962,p. 96.        *
*    (CACM = Communications of the Association for Computing Machinery)       *
*    AVINT was copied from P.J.Davis and P.Rabinowitz,NUMERICAL INTEGRATION,  *
*    Blaisdell Publishing Company,1967.                                       *
*                                                                             *
*    Some parts coded by PabloF (20/07/1999) (pablof@cab.cnea.gov.ar)         *
*-----------------------------------------------------------------------------*
*/
double AVINT(double *x, double *y, int n, double xlo, double xup)
{
      int i,j,ib,jm;
      double sum,syl,syu,x1,x2,x3,term1,term2,term3,a,b,c,ca=0,cb=0,cc=0;

      syl=xlo;

      ib=2;
      for( i=1; i<=n; i++)
         if(x[i] >= xlo) break;
         else            ib++;


      j=n;
      for( i=1; i<=n; i++)
         if(xup >= x[j]) break;
         else            j--;


      j=j-1;
      sum=0.;
      for( jm=ib; jm<=j; jm++)
                {
                x1 = x[jm-1];
                x2 = x[jm];
                x3 = x[jm+1];

                term1 = y[jm-1]/((x1-x2)*(x1-x3));
                term2 = y[jm]  /((x2-x1)*(x2-x3));
                term3 = y[jm+1]/((x3-x1)*(x3-x2));

                a =  term1+term2+term3;
                b = -(x2+x3)*term1-(x1+x3)*term2-(x1+x2)*term3;
                c = x2*x3*term1+x1*x3*term2+x1*x2*term3;

                if(jm == ib) {ca=a; cb=b; cc=c;  }
                else            {
                                ca=0.5*(a+ca);
                                cb=0.5*(b+cb);
                                cc=0.5*(c+cc);
                                }
                syu=x[jm];
                sum += ca*(pow3(syu)-pow3(syl))/3. + cb*0.5*(pow2(syu)-pow2(syl)) + cc*(syu-syl);

                ca=a;
                cb=b;
                cc=c;
                syl=syu;
                }

double AVINTl = sum + ca*(pow3(xup)-pow3(syl))/3. + cb*0.5*(pow2(xup)-pow2(syl)) + cc*(xup-syl);

return AVINTl;
}
//c**********************************************************************
//*
//*    PabloF's library
//*    Contents:           AVINT.f -----> Simpson
//*                        GaussLaguerre.f -----> Gauss-Laguerre
//*                        POCH.f ------> Pochammer symbol
//*
//*-----------------------------------------------------------------------------*
//*                                                                             *
//*    This program computes the integral of y(x) in XLO < x < XUP .            *
//*                                                                             *
//*    X and Y are N dimensional arrays containing the abscissas x and the      *
//*    ordinates y respectively.The x are usually unequally spaced and must     *
//*    be given in ascending order.The y are usually data.XLO must be less      *
//*    than XUP,but otherwise there are no restrictions on XLO and XUP except,  *
//*    of course that of accuracy which requires that if XLO is less than x1,   *
//*    it be close to x1,and if XUP is greater than xN,it be close to xN.       *
//*                                                                             *
//*    AVINT is adapted from P.E.Hennion,Algorithm 77,CACM 5,1962,p. 96.        *
//*    (CACM = Communications of the Association for Computing Machinery)       *
//*    AVINT was copied from P.J.Davis and P.Rabinowitz,NUMERICAL INTEGRATION,  *
//*    Blaisdell Publishing Company,1967.                                       *
//*                                                                             *
//*    Some parts coded by PabloF (20/07/1999) (pablof@cab.cnea.gov.ar)         *
//*-----------------------------------------------------------------------------*
double GaussLaguerre(int n, double (*fun)(double x), double alpha, double qmin)
{

double sum=0.;

       if(n <= 15)
                {
                for(int i=1; i<=15; i++)
                      sum += Gw15[i]*fun(Gx15[i]/alpha+qmin);
                }
 else if(n <= 20)
                {
                for(int i=1; i<=20; i++)
                    sum += Gw20[i]*fun(Gx20[i]/alpha+qmin);
                }
 else if(n <= 24)
                {
                for(int i=1; i<=24; i++)
                    sum += Gw24[i]*fun(Gx24[i]/alpha+qmin);
                }
 else {
                for(int i=1; i<=48; i++)
                    sum += Gw48[i]*fun(Gx48[i]/alpha+qmin);
      }

double  GALAGs=sum/alpha;

return GALAGs;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//  Pochhammer's symbol

void POCH(int n, complex<double> cx, complex<double>& cpoch)
{

//
//  Routine  to calculate  the POCHAMMER symbol
// (A)n = A(A+1)(A+2)...(A+n-1) , (A)o = 1
//  Author: PabloF 3/2/94    Version 1.1
//  New version extended to complex arguments
//

#define nmax 100

      complex<double> cvpo[nmax+1];   //cvpo[0:nmax];
      complex<double> cnil(0.,1.);
      complex<double> cone(1.,1.);

      n = min(nmax,n);  // Oleg correction

        if(n == 0)                 { cpoch=cone;    }
         else      {
                   if(cx == cnil) {  cpoch=cnil;    }
                   else {
                        cvpo[0]=cone;
                        for(int i=0; i<=n-1; i++)  cvpo[i+1] = cvpo[i]*( cx + complex<double>(i,0.) );

                        cpoch=cvpo[n];
                        }
                   }
}

/*
To optimize this routine for calculating the Pochhammer symbol (A)n=A(A+1)(A+2)…(A+n−1)(A)n​=A(A+1)(A+2)…(A+n−1), we can incorporate an asymptotic expansion for large nn. The key idea is to avoid iterating through the entire product for large nn, leveraging the properties of the Gamma function and its logarithm (related to the digamma function).

The Pochhammer symbol can be expressed using the Gamma function:
(A)n=Γ(A+n)Γ(A).
(A)n​=Γ(A)Γ(A+n)​.
#include <complex>
#include <cmath>
#include <algorithm> // for std::min
#include <iostream>  // for debugging

void POCH(int n, std::complex<double> cx, std::complex<double>& cpoch) {
    const int nmax = 100; // Maximum n to calculate directly
    const double epsilon = 1e-12; // Precision threshold for asymptotic approximation

    n = std::min(nmax, n); // Limit n to nmax
    std::complex<double> cone(1.0, 0.0); // Represent 1 + 0i

    if (n == 0) {
        cpoch = cone; // (A)0 = 1
        return;
    }

    if (std::abs(cx) < epsilon) {
        // If A is very close to zero, the Pochhammer symbol collapses
        cpoch = std::complex<double>(0.0, 0.0);
        return;
    }

    // For small n, use direct computation
    if (n <= 10) {
        std::complex<double> product = cone;
        for (int i = 0; i < n; ++i) {
            product *= (cx + std::complex<double>(i, 0.0));
        }
        cpoch = product;
        return;
    }

    // For larger n, use the asymptotic expansion
    std::complex<double> z = cx + std::complex<double>(n - 1, 0.0);
    std::complex<double> log_gamma_ratio = std::log(std::tgamma(z + cone)) - std::log(std::tgamma(cx));
    cpoch = std::exp(log_gamma_ratio);
}
*/
