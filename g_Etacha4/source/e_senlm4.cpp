#include "e_Etacha4.h"
#include <complex>
using namespace std;

#include <QtMath>
#include "../win/e_myextern.h"
#include "../win/e_Constant.h"
#define sqrtPI 1.77245385090552

extern double pow2(double par);
extern double pow3(double par);
extern double pow2I(double par);
//extern double pow_int(double par, int power);
extern int pow_m1(int power);
extern complex<double> pow2(complex<double> par);
extern complex<double> LogGammaFunc(complex<double> cz);

extern double SeSE[67],StSE[13];
extern int EtachaVersion;
extern char BufRtf[1000];

//void SEnlm(double zp8, double zt8, double E8, TForm *tobject);
void gamac(complex<double>za, int ny, double test1, double test2, int& nt,complex<double> &logam);
double Factorial(int n);
void INTDEF(double PRE, double Ymin, double Ymax, double (*cfun)(double eta), double& cint);
double rns  (double eta);
double rnp0 (double eta);
double rnp1 (double eta);
double rnp1m(double eta);

double PlnmLegendr(int L, int M, double X);
void CF21D(complex<double> ca,complex<double> cb,complex<double> cc,double x, complex<double>& cf);
void ARMONIC(double s,int l,int m, double tk, double fik, complex<double>& cylm);

double ClebschGordan(int l1, int l2, int l3, int m1, int m2, int m3);
void LogFac(void);
double LogCle(int l1,int l2,int l3,int m1,int m2,int m3);
double HyperGeom_function(double A,double B,double C, double X);
double GAMMA_function(double X);
double PSI_function(double X);
double cotan(double x);

const complex<double>       ci    (0.,1.);
const complex<double>       cuno  (1.,0.);
const complex<double>       cdos  (2.,0.);
const complex<double>       czero (0.,0.);

//struct {int n0,l0,m0,nf,lf,mf;} m5;
int _n0,_l0,_m0,_nf,_lf,_mf;

double epsi;
double scoa;
complex<double> ca2,ca3,ca4;

complex<double> cinu,con0,calfa;

extern double sc;
extern int icor;
extern double zp,zt;
extern double pv;
extern complex<double> ca1;

#define _idf 140
#define _idfL 141
double fact_fac[_idfL];
bool  LogFacInit = false;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ETACHA::SEnlm(double zp8, double zt8, double E8)
{
//********1*********1*********1*********1*********1*********1********1**
// Pour ETACHA4
//.....calcula la secc. eficaz total de excitacion desde los
//.....estados ns, np1 y np-1 al estado final nlm, por S.Eikonal
//.....Autor: Cesar A. Ramirez            Fecha: 10-12-2012
//         zt=carga del nucleo blanco
//         zp=carga del nucleo proyectil
//         ei=energia del estado ligado inicial
//         ef=energia del estado ligado final
//         n,l,m = numeros cuanticos del estado final
//         eta=componente normal del momento transferido
//**********************************************************************
//----- Program to calculte total cross sections of
//----- monoelectronic atoms excitation, for colision with nude projectiles
//----- in the SE (Symmetric-Eikonal) aproximation.
//--------------------------------------------------
//......changes in the present version (JPR 02/2012)
//----- zp: Projectile charge (to be excited)
//----- zt: Target charge
//......SC and ASC corrections
//-------------------------------------------------
//----- (n0,l0,m0) quantum numbers of initial state
//----- (n,l,m) quantum numbers of final state
//----- SUBROUTINES:
//----- * rns, rnp0, rnp1, rnp1m: excitation from ns, np0, np1 np-1 to n'lm
//-----                    respectively with n'lm arbitraries (n'>n)
//----- * intdef: for definite integrals from 0 to infinity
//----- * Factorial: factorial function
//----- * ClebschGordan: Clebsch-Gordan coeficients
//----- * ARMONIC: spherical hARMONIC function
//----- * plgndr: Legendre polinomials
//----- * gamac: gama function
//----- * CF21D, chyper, hypgfx: hypergeometric functions F21(ca,cb,cc,cx)
//**********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

//Oleg       double sexnlm[251];
//      double signlm[5][3][4][6][6][7],signl[5][3][6][6] /*,signt[5][3][6]*/;
      double signlm[5][2][3][6][5][6],signl[5][2][6][5] /*,signt[5][3][6]*/;


//      int n0,l0,m0,n,l,m;

      bool Flag34 = EtachaVersion >= etacha_v34;

      int n0_max = Flag34 ? 4 : 3;
      int nf_max = Flag34 ? 5 : 4;

      for(int  i=0; i<67; i++) SeSE[i]=0;
      for(int  i=0; i<13; i++) StSE[i]=0;

      double *t = &signlm[0][0][0][0][0][0];
      for(int  i=0; i<5400; i++) *t++ = 0;;

      t = &signl[0][0][0][0];
      for(int  i=0; i<300; i++) *t++ = 0;;
      //-----------------------------------------------------------------
      zp=zp8;
      zt=zt8;
      double E=E8;
      double ga2=pow2I(1.+E/931.5);
      pv=1.3703599976e2*sqrt(1.-ga2);
      // temporary  Oleg      pv = Velocity_au(Ecurr);
//----------------------------------------------------------------------
      double scf=1.13*pow(zt,(1./3.));
      double fs1=1.;
            if (zt == 1.) fs1=1.23;
       else if (zt == 2.) fs1=1.526;
       else if (zt == 6.) fs1=0.78;
       else if (zt == 7.) fs1=0.85;
       else if (zt == 10.) fs1=1.04;
       else if (zt == 13.) fs1=0.58;
       else if (zt == 14.) fs1=0.59;
       else if (zt == 18.) fs1=0.68;
       else if (zt == 29.) fs1=0.672;
       else if (zt == 36.) fs1=0.61;
       else if (zt == 54.) fs1=0.535;

      sc=scf*fs1;
      icor=1;
//----------------------------------------------------------------------
      double xnu = zt / pv;

      cinu = ci*xnu;
      ca1 = 1.   +cinu;
      ca2 = 1.   -cinu;
      ca3 = 1.+2.*cinu;
      ca4 = 1.-2.*cinu;
      double senhy = sinh(PI*xnu);

      double  test1=1.e-6;
      double  test2=1.e-6;

      complex<double> cloga,cgama;
      int nt;
      //void gamac(complex<double>za, int ny, double test1, double test2, int nt,complex<double> &logam);
      gamac(cinu,1,test1,test2,nt,cloga);
      cgama = exp(cloga);
//-----------------------------------------------------------------------
//     open(40,file='tnlm4.dat',status='unknown')
//      write(40,200)zp,zt,e
//     if (icor.eq.1) then
//     write(40,*)' with SC+ASC corrections'
//     else
//     write(40,*)' no   SC+ASC corrections'
//     endif
//200   format(' zp=',d15.5,' zt=',d15.5,' E=',e15.5)
//-----------------------------------------------------------------------
    int index=1;
    for(_n0=1; _n0<=n0_max; _n0++)
        {
        int l0m=min(_n0-1,1);
        for(_l0=0; _l0<=l0m; _l0++)
           {
           for(_m0=0; _m0<=_l0; _m0++)
                {
                int nfm=_n0+1;
                for( _nf=nfm; _nf<=nf_max; _nf++) {   // L110
                   for( _lf=0; _lf<=_nf-1; _lf++) { //L120
                        for(_mf=0; _mf<=_lf; _mf++)  //L130
                              {
                              int n10=_n0-_l0-1;
                              int n20=_n0+_l0;
                              double fac10 = Factorial(n10);
                              double fac20 = Factorial(n20);
                              double t15 = pow(zp,1.5);
                              double xn0   = 2./pow2(_n0)* t15 * sqrt(fac20*fac10)*pow(2.*zp/_n0,_l0);
                              int n1=_nf-_lf-1;
                              int n2=_nf+_lf;
                              double fac1  = Factorial(n1);
                              double fac2  = Factorial(n2);
                              double xnl   = 2./pow2(_nf )* t15 * sqrt(fac2 *fac1 )*pow(2.*zp/_nf ,_lf );
                               //----------------------------------------------------------------------
                              double ei = -pow2(zp/_n0)/2.;
                              double ef = -pow2(zp/_nf )/2.;

                              epsi =(ef-ei)/pv;
                              if(epsi == 0.)epsi=1.e-6;
                                //.............................
                               if (pv > 1.75*epsi)  scoa=1.-pow2((1.75*epsi/pv));
                               else                 scoa=0.0;
                                //.............................
                               calfa = 2.*ci*pv*epsi;
                               double beta  = pow2((2.*pv*epsi));
                                //----------------------------
                               con0 = 4.*PI*pv*pow2(PI*xnu/cgama/senhy)/beta*xn0*xnl;
                                //  original commented   *        /cdexp(cinu*dlog(beta))         cte. de módulo=1
                                //----------------------------
                              double ymin=0.;
                              double ymax=3.;
                              double a02=2.8e-17;
                              double pre=1.e-3;
                                //---------------------------------------------
                              double xint=0, tnlm=0, tp1=0, tp1m=0;

                               if(_l0 == 0)
                                      {
                                      INTDEF(pre,ymin,ymax,rns,xint);
                                      tnlm=2*PI*xint*a02;
                                      }
                                 else {
                                       if(_m0 == 0)
                                                {
                                                INTDEF(pre,ymin,ymax,rnp0,xint);
                                                tnlm=2*PI*xint*a02;
                                                }
                                         else   {
                                                INTDEF(pre,ymin,ymax,rnp1,xint);
                                                tp1=2*PI*xint*a02;
                                                if(_mf != 0)tp1 *= 2;

                                                INTDEF(pre,ymin,ymax,rnp1m,xint);
                                                tp1m=2*PI*xint*a02;
                                                if(_mf != 0)tp1m *=2;
                                                tnlm=tp1+tp1m;
                                                }
                                      }
                        //--------------------------------------------------------------------
                        //---- Se multiplica por 2 la sección eficaz correspondiente a
                        //---- la transición desde un estado inicial con m=0
                        //---- a un estado final con m distinto a cero, para tener en cuenta los
                        //---- dos casos simétricos +m y -m.
                        //--------------------------------------------------------------------

                       if(_m0 == 0 && _mf != 0)tnlm *=2 ;

                                //      write(40,300)_n0,_l0,_m0,_nf,_lf,m,tnlm,index,sc,scoa
                                //L300: format(' inicial=',3i3,'   final=',3i3,'  tot = ',d15.5,
                                //   '   index=',i3,' coef ASC=',2d15.5);

//Oleg                       sexnlm[index]=tnlm;       // not used in calculation
                      signlm[_n0][_l0][_m0][_nf][_lf][_mf]=tnlm;
                      index ++;
                           //----------------------------------------------------------------------
                      sprintf(BufRtf,"<b>i</b>: %d.%d.%d  &rarr;  <b>f</b>: %d.%d.%d &nbsp; tnlm: %9.3e", _n0,_l0,_m0,_nf,_lf,_mf,tnlm); emit appendShell(0,0);

                      } // L130:
                 } //  L120:
               }  //  L110:
       }}}       //  L100:
//-----------------------------------
//     somme sur les _m0 et m
//     write(40,*)' '
      index=1;

      for( _n0=1; _n0<=n0_max; _n0++) //L140
        {
        int l0m=min(_n0-1,1);
        int nfm=_n0+1;
           for( _l0=0; _l0<=l0m; _l0++)  //L140
                for( _nf=nfm; _nf<=nf_max; _nf++)     //L140
                    for( _lf=0; _lf<=_nf-1; _lf++)     //L140
                      {
                      double sigt=0.0;
                      for( _m0=0; _m0<=l0m; _m0++)
                          for( _mf=0; _mf<=_lf; _mf++)
                              sigt += signlm[_n0][_l0][_m0][_nf][_lf][_mf];

                      if(_l0 == 1) sigt=sigt/3.;
                      signl[_n0][_l0][_nf][_lf]=sigt;
                      SeSE[index]=sigt;
                      index ++;
                                                        //      write(40,310)_n0,_l0,_nf,_lf,sigt
         }} // L140:
                                                        //     write(40,*)' '
                                                        //     somme sur les _lf pour _nf=4 et 5
      index=1;
      for( _nf=4; _nf<=nf_max; _nf++)    // L160
        {
        int nim=_nf-1;
        for(_n0=1; _n0<=nim; _n0++)   // L160
                {
                int l0m = min(_n0-1,1);
                for(_l0=0; _l0<=l0m; _l0++)   // L160
                      {
                      double sigtn=0.;
                      for(_lf=0; _lf<=_nf-1; _lf++)
                              sigtn += signl[_n0][_l0][_nf][_lf];

 // not used - Oleg                     signt[_n0][_l0][_nf]=sigtn;
                      StSE[index]=sigtn;
                      index++;
                                                        //      write(40,311)_n0,_l0,_nf,sigtn
        }}} //L160:
                                                        //      close(unit=40)

                                                        //L310: format(' inicial=',2i3,'      final=',2i3,'     tot = ',d15.5);
                                                        //L311: format(' inicial=',2i3,'      final=',i3,'        tot = ',d15.5);
}
//*************************************** fin  programa principal ******

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//.....function=eta * /R(eta)/ ^ 2      en la aprox. SE
//.....desde ns a nlm
//***********************************************************************
double rns(double eta)
{
      complex<double> chyp1,chyp2;
      double gama=pow2(eta)+pow2(epsi);
      double xk=sqrt(gama);
      double znz=zp/_nf+zp/_n0;
      double vari=-pow2(eta/epsi);
      CF21D(cinu,cinu,cdos,vari,chyp1);
      CF21D(ca1 ,ca1, cdos,vari,chyp2);

      chyp1=(ca2*chyp1-cinu*(1.-vari)*chyp2)/ca4;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      complex<double> csum1=czero;
      complex<double> csum2=czero;
      complex<double> csum3=czero;

      int minla=max(0,_lf-1);
      int maxla=_lf+1;

//-----

for(int la=minla; la<=maxla; la++)
      {
      double ga=sqrtPI;
      for(int  j=0; j<=la; j++)
                         ga *= (0.5+j);


      double xla=sqrtPI*pow(xk,la)/pow(2.,(la+1))/ga;
      complex<double> sumq1(czero);
      complex<double> sumq (czero);
      int maxiq=_nf-_lf-1;

      for(int iq=0; iq<=maxiq; iq++)
              {
              int n4=_nf-_lf-1-iq;
              int n5=iq;
              int n6=2*_lf+1+iq;
              double fac4=Factorial(n4);
              double fac5=Factorial(n5);
              double fac6=Factorial(n6);
              double xq=pow((-2.*zp/_nf),iq)/fac4/fac5/fac6;
              complex<double> sump1 (czero);
              complex<double> sump  (czero);
              int maxip=_n0-1;

              //-------------------------------------------------------
              for(int ip=0; ip<=maxip; ip++) {
                      int n40=_n0-1-ip;
                      int n50=ip;
                      int n60=1+ip;
                      double fac40=Factorial(n40);
                      double fac50=Factorial(n50);
                      double fac60=Factorial(n60);
                      double xp = pow((-2.*zp/_n0),ip)/fac40/fac50/fac60;
                        //-------------------------------------------------
                      double xmu=2.5+_lf+iq+ip;
                      double xnu=la+0.5;
                      double xa=(xmu+xnu)/2.;
                      double xc=xnu+1.;
                      double xz=-gama/pow(znz,2);
                      double xc1=2.*xa;
                      double xc2=2.*xa-xc+1.;
// OLEG                      complex<double> yz=(1.-sqrt(1.-xz))/(1.+sqrt(1.-xz))*cuno;
                      double yz=(1.-sqrt(1.-xz))/(1.+sqrt(1.-xz));
                      int n3 = 3+_lf+la+iq+ip;
                      double h1=pow((2./(1.+sqrt(1.-xz))),n3 );
                      double hyp31 = HyperGeom_function(xc1,xc2,xc,yz);
                      int nn2 = 2+_lf+la+iq+ip;
                      double fac3 = Factorial(nn2);
                      double xintx1 = xla*fac3/pow(znz,n3)*h1*hyp31;

                      if(la == _lf)      sump1 = sump1 + xp*xintx1;
                      else {
                              xc1=2.*xa-1. ;
                              xc2=2.*xa-xc ;
                              double h2 = pow((2./(1.+sqrt(1.-xz))),(n3-1));
                              double hyp32  = HyperGeom_function(xc1,xc2,xc,yz);
                              double xintx2 = xla*fac3/pow(znz,n3)*h2*hyp32;
                              sump += xp*(-zp/_n0*xintx1 + znz*ip/(2.*xa-1.)*xintx2);
                              }
                        }
              //-------------------------------------------------------

              if(la == _lf) sumq1  += sump1*xq;
              else          sumq   += sump *xq;

      }  // end do
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//              s=-1      indica que es un ARMONICo conjugado
      double s   = -1.;
      double xkz = -epsi;
      double tk  = acos(xkz/xk);
      complex<double> calga = calfa/gama;
      double xl=sqrt((2.*la+1.)/(2.*_lf+1.)/4./PI);

      complex<double> cylam,cylam1,cylam2,cylam3;

      if(la == _lf)
                {
                ARMONIC(s,_lf,_mf, tk, 0., cylam);
                csum1=pow( ci,la)*xl*sumq1*cylam*calfa*(ca2/cinu*chyp1+chyp2);
                }
         else   {

                double xc3 = ClebschGordan(la,1,_lf,    0, 0,  0);
                double xc4 = ClebschGordan(la,1,_lf,_mf+1,-1,_mf);
                double xc5 = ClebschGordan(la,1,_lf,_mf-1, 1,_mf);
                double xc6 = ClebschGordan(la,1,_lf,_mf  , 0,_mf);

                ARMONIC(s,la,_mf+1,tk,0.,cylam1);
                ARMONIC(s,la,_mf-1,tk,0.,cylam2);
                ARMONIC(s,la,_mf  ,tk,0.,cylam3);

                csum2 += pow(ci,la)*xl*sumq*xc3*(xc4*cylam1-xc5*cylam2)/sqrt(2.) *ci*eta*calga*ca2/cinu*chyp1 ;
                csum3 += pow(ci,la)*xl*sumq*xc3*xc6*cylam3  *(2.*pv*chyp2-ci*epsi*calga*ca2/cinu*chyp1) ;
              }
    }

//-----------------------------------------------------------------------
      complex<double>  csum = 4.*PI* (csum1 - 2. * (csum2+csum3));
//-----------------------------------------------------------------------
      complex<double> cf = ci/(2.*PI*pv) * con0 * csum ;
//commented original      *         * cdexp(2*cinu*dlog(gama))    cte de modulo=1

      double xf=eta*pow2(abs(cf));
//...............................
//     SC and ASC corrections
       if (icor == 1)
                {
                double sq=pow2I(1.+pow2(sc/xk))+scoa*(1.-pow2I(1.+pow2(xk/sc)))/zt;
                xf *=sq;
                }
//...............................
return xf;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//***********************************************************************
//.....function=eta * /R(eta)/^2      en la aprox. SE
//.....desde np0 a nlm
//***********************************************************************
double rnp0(double eta)
 {
      complex<double> chyp1,chyp2;
      double gama=pow2(eta)+pow2(epsi);
      double xk=sqrt(gama);
      double znz=zp/_nf+zp/_n0;
//-----
      double vari=-pow2(eta/epsi);

      CF21D(cinu,cinu,cdos,vari,chyp1);
      CF21D(ca1 ,ca1, cdos,vari,chyp2);
      chyp1=(ca2*chyp1-cinu*(1.-vari)*chyp2)/ca4;

      complex<double> csum1(czero);
      complex<double> csum2(czero);
      complex<double> csum3(czero);
      complex<double> csum4(czero);
      complex<double> calga(czero);

int minla2=abs(_lf-2);
int minla1=abs(_lf-1);
int minla=min(minla1,minla2);

if(_lf == 0)minla=0;

int maxla=_lf+2;

//-----------------------------------------------------
for(int la=minla; la<=maxla; la++)
      {
      complex<double> sumq2(czero);
      complex<double> sumq1(czero);
      complex<double> sumq (czero);
      int  maxiq=_nf-_lf-1;

      for(int iq=0; iq<=maxiq; iq++)
         {
         int n4 = _nf-_lf-1-iq;
         int n5 = iq;
         int n6 = 2*_lf+1+iq;
         double fac4=Factorial(n4);
         double fac5=Factorial(n5);
         double fac6=Factorial(n6);
         double xq =pow( (-2.*zp/_nf),iq)/fac4/fac5/fac6;
         complex<double> sump1(czero);
         complex<double> sump2(czero);
         complex<double> sump (czero);
         int maxip=_n0-_l0-1;

          for(int ip=0; ip<=maxip; ip++)
                {
                 int n40=_n0-_l0-1-ip;
                 int n50=ip;
                 int n60=2*_l0+1+ip;
                 double fac40=Factorial(n40);
                 double fac50=Factorial(n50);
                 double fac60=Factorial(n60);
                 double xp = pow( (-2*zp/_n0),ip)/fac40/fac50/fac60;
                //--------------------------------------------------
                 int nn3 = 3+_lf+_l0+la+iq+ip;
                 int nn2 = 2+_lf+_l0+la+iq+ip;
                 double fn2 = Factorial(nn2);
                 double xmu=  2.5+_lf+_l0+iq+ip;
                 double xnu=  la+0.5;
                 double xc=la+1.5;
                 double xxa=(xnu+xmu)/2;
                 double xxb=(xnu-xmu)/2+0.5;
                 double xxz=gama/(gama+pow2(znz));
                 double hyp31 = HyperGeom_function(xxa,xxb,xc,xxz);
                        hyp31 /=pow(sqrt(gama+pow2(znz)),nn3);

                 if(la == _lf+1 || la == abs(_lf-1))   sump1  += xp*fn2*hyp31;
                 else {
                       xxa=xxa-0.5;
                       xxb=xxb+0.5 ;
                       double hyp32 = HyperGeom_function(xxa,xxb,xc,xxz);
                       hyp32 /= pow(sqrt(gama+pow2(znz)),nn2);
                       sump  += xp*fn2*( ip*hyp32/nn2 - zp/_n0*hyp31 );
                       if(la == _lf)
                              sump2 += xp*fn2*( (ip+3)*hyp32/nn2 - zp/_n0*hyp31 );
                      }
                 }  // end do

            if(la == _lf+1 || la == abs(_lf-1))       sumq1 += xq*sump1;
            else                                      sumq  += xq*sump;

            if(la == _lf)      sumq2 += xq*sump2;

            }  // end do
                //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                //                           s=-1 (indica el conjugado)
     double  s=-1.;
     double  ga=sqrtPI;

     for(int  j=0; j<=la; j++) ga *= (0.5+j);

     double  xkz=-epsi;
     double  tk=acos(xkz/xk);

     calga=calfa/gama;

     double  xla=sqrtPI*pow(xk,la)/pow(2.,(la+1))/ga * sqrt((2*la+1.)/(2*_lf+1.)/4./PI);

     complex<double> cylam,cylam1,cylam2,cylam3;

     if(la == _lf+1 || la == abs(_lf-1))
                {
                double xc1 = ClebschGordan(la, 1, _lf,  0, 0,  0);
                double xc2 = ClebschGordan(la, 1, _lf,_mf, 0,_mf);
                ARMONIC(s,la,_mf,tk,0.,cylam1);
                csum1 += pow( ci,la)*xla*sqrt(3.)*sumq1*xc1*xc2*cylam1;
                }
     else       {
                double xc3 = ClebschGordan(la, 2, _lf,   0,    0,    0);
                double xc4 = ClebschGordan(la, 2, _lf, _mf,    0,  _mf);
                double xc5 = ClebschGordan(la, 2, _lf, _mf+1, -1,  _mf);
                double xc6 = ClebschGordan(la, 2, _lf, _mf-1,  1,  _mf);

                ARMONIC(s,la,_mf,  tk,0.,cylam );
                ARMONIC(s,la,_mf+1,tk,0.,cylam2);
                ARMONIC(s,la,_mf-1,tk,0.,cylam3);

                csum2 += pow( ci,la)*xla*sumq*xc3* (xc5*cylam2-xc6*cylam3) /sqrt(2.);
                csum3 += pow( ci,la)*xla*sumq*xc3*  xc4*cylam*2./sqrt(3.);

                if(la == _lf)
                        csum4 = pow(ci,_lf)*xla*sumq2/sqrt(3.)*cylam; // Attention!!      csum4 is not SUM
                }
      }  // end do


      csum1*=calfa*(ca2/cinu*chyp1+chyp2);
      csum2*=ci*eta*calga*ca2/cinu*chyp1;   // Attention!! calga is taken from cycle

      complex<double> ct= (2*pv*chyp2-ci*epsi*calga*ca2/cinu*chyp1);
      csum3*= ct;
      csum4*= ct;
//-----------------------------------------------------------------------
      complex<double> csum = 4*PI*(  csum1 - 2.*(csum2+csum3+csum4) );
//-----------------------------------------------------------------------
      complex<double> cf = ci/(2.*PI*pv) * con0 * csum ;
//     *          * cdexp(2.*cinu*dlog(gama))      cte de modulo=1
      double xf = eta * pow2(abs(cf));
//...............................
//     SC and ASC corrections
       if (icor == 1) {
                      double sq=pow2I(1.+pow2(sc/xk)) +  scoa*(1.-pow2I(1.+pow2(xk/sc)))/zt;
                      xf *=sq;
                      }
//...............................
return xf;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//***********************************************************************
//.....function=eta * /R(eta)/2      en la aprox. SE
//.....desde np1 a nlm
//***********************************************************************
double rnp1(double eta)
{

      complex<double> chyp1,chyp2;
      double gama=pow2(eta)+pow2(epsi);
      double xk=sqrt(gama);
      double znz=zp/_nf+zp/_n0;
//-----
      double vari=-pow2(eta/epsi);

      CF21D(cinu,cinu,cdos,vari,chyp1);
      CF21D(ca1 ,ca1, cdos,vari,chyp2);
      chyp1=(ca2*chyp1-cinu*(1.-vari)*chyp2)/ca4;

      complex<double> csum1(czero);
      complex<double> csum2(czero);
      complex<double> csum3(czero);
      complex<double> csum4(czero);
      complex<double> calga(czero);

int minla2=abs(_lf-2);
int minla1=abs(_lf-1);
int minla=min(minla1,minla2);

if(_lf == 0)minla=0;

int maxla=_lf+2;

//-----
for(int la=minla; la<=maxla; la++)
      {
      complex<double> sumq2(czero);
      complex<double> sumq1(czero);
      complex<double> sumq (czero);
      int  maxiq=_nf-_lf-1;

      for(int iq=0; iq<=maxiq; iq++)
         {
         int n4 = _nf-_lf-1-iq;
         int n5 = iq;
         int n6 = 2*_lf+1+iq;

         double fac4=Factorial(n4);
         double fac5=Factorial(n5);
         double fac6=Factorial(n6);
         double xq =pow( (-2.*zp/_nf),iq)/fac4/fac5/fac6;
         complex<double> sump1(czero);
         complex<double> sump2(czero);
         complex<double> sump (czero);
         int maxip=_n0-_l0-1;

         for(int ip=0; ip<=maxip; ip++)
                {
                int n40 = _n0-_l0-1-ip;
                int n50 = ip;
                int n60=2*_l0+1+ip;
                double fac40=Factorial(n40);
                double fac50=Factorial(n50);
                double fac60=Factorial(n60);
                double xp = pow( (-2*zp/_n0),ip)/fac40/fac50/fac60;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                int nn3 = 3+_lf+_l0+la+iq+ip;
                int nn2 = 2+_lf+_l0+la+iq+ip;
                double fn2 = Factorial(nn2);
                double xmu= 2.5+_lf+_l0+iq+ip;
                double xnu= la+0.5;
                double xa= (xnu+xmu)/2;                  // differ with rnp0
                double xb= xa+0.5;                       // differ with rnp0
                double xc= xnu+1.;                       // differ with rnp0
                double xxb= xc-xb;                       // differ with rnp0
                double xxz= gama/(gama+pow2(znz));
                double hyp31 = HyperGeom_function(xa,xxb,xc,xxz);
                       hyp31 /= pow(sqrt(gama+pow2(znz)),nn3);

               if(la == _lf+1 || la == abs(_lf-1))       sump1 += xp*fn2*hyp31;
               else {
                     double xxa = xa-0.5;                   // differ with rnp0
                            xxb = xc-xa ;                   // differ with rnp0
                     double hyp32 = HyperGeom_function(xxa,xxb,xc,xxz);
                            hyp32 /= pow(sqrt(gama+pow2(znz)),nn2);
                            sump += xp*fn2*( ip*hyp32/nn2 - zp/_n0*hyp31 );

                     if(la == _lf)     sump2 += xp*fn2*((ip+3)*hyp32/nn2 - zp/_n0*hyp31 );
                     }
               }  // end do

       if(la == _lf+1 || la == abs(_lf-1)) sumq1 += xq*sump1;
       else {
             sumq  += xq*sump;
             if(la == _lf)      sumq2 += xq*sump2;
             }
      }  // end do
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//                           s=-1 (indica el conjugado)
     double s=-1.;
     double ga=sqrtPI;

     for(int j=0; j<=la; j++)  ga *= (0.5+j);

      double xkz = -epsi;
      double tk  = acos(xkz/xk);
      calga=calfa/gama;

      double xla=sqrtPI*pow(xk,la)/pow(2.,(la+1))/ga  *sqrt((2*la+1.)/(2*_lf+1.)/4./PI);

      complex<double> cylam,cylam1,cylam2,cylam3;

      if(la == _lf+1 || la == abs(_lf-1))
                {
                double xc1 = ClebschGordan(la,1,_lf,    0, 0,   0);
                double xc2 = ClebschGordan(la,1,_lf,_mf-1, 1, _mf);
                ARMONIC(s,la,_mf-1,tk,0.,cylam1);
                csum1=csum1 +pow( ci,la)*xla*sqrt(3.)*sumq1*xc1*xc2*cylam1;
                }
      else      {
                double xc3 = ClebschGordan(la,2,_lf,    0, 0,  0);
                double xc4 = ClebschGordan(la,2,_lf,_mf  , 0, _mf);
                double xc5 = ClebschGordan(la,2,_lf,_mf-2, 2, _mf);
                double xc6 = ClebschGordan(la,2,_lf,_mf-1, 1, _mf);

                ARMONIC(s,la,_mf,  tk,0.,cylam );
                ARMONIC(s,la,_mf-2,tk,0.,cylam2);
                ARMONIC(s,la,_mf-1,tk,0.,cylam3);

                csum2 += pow( ci,la)*xla*sumq*xc3 /sqrt(6.)*(xc4*cylam-sqrt(6.)*xc5*cylam2);
                csum3 += pow( ci,la)*xla*sumq*xc3 *xc6*cylam3;

               if(la == _lf)  csum4=-pow(ci,_lf)*xla*sumq2/sqrt(6.)*cylam;    // Attention!!      csum4 is not SUM
            }

     }  // end do

      csum1 *= calfa*(ca2/cinu*chyp1+chyp2);
      csum2 *= ci*eta*calga*ca2/cinu*chyp1;
      csum3 *= (2*pv*chyp2-ci*epsi*calga*ca2/cinu*chyp1);
      csum4 *= ci*eta*calga*ca2/cinu*chyp1;
//-----------------------------------------------------------------------
      complex<double> csum = 4*PI* ( csum1 - 2. *(csum2+csum3+csum4));
//-----------------------------------------------------------------------
      complex<double>cf= ci/(2.*PI*pv) * con0 * csum;
//     *           * cdexp(2.*cinu*dlog(gama))   cte. de modulo=1
      double xf=eta*pow2(abs(cf));
//...............................
//     SC and ASC corrections
       if (icor == 1) {
                      double sq=pow2I(1.+pow2(sc/xk))+scoa*(1.-pow2I(1.+pow2(xk/sc)))/zt;
                      xf *= sq;
                      }
//...............................
return xf;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double rnp1m(double eta)
//***********************************************************************
//.....function=eta * /R(eta)/2      en la aprox. SE
//.....desde np-1 a nlm
//***********************************************************************
{
      complex<double> chyp1,chyp2;
      double gama=pow2(eta)+pow2(epsi);
      double xk=sqrt(gama);
      double znz=zp/_nf+zp/_n0;
//-----
      double vari=-pow2(eta/epsi);

      CF21D(cinu,cinu,cdos,vari,chyp1);
      CF21D(ca1 ,ca1, cdos,vari,chyp2);
      chyp1=(ca2*chyp1-cinu*(1.-vari)*chyp2)/ca4;

      complex<double> csum1(czero);
      complex<double> csum2(czero);
      complex<double> csum3(czero);
      complex<double> csum4(czero);
      complex<double> calga(czero);

int minla2=abs(_lf-2);
int minla1=abs(_lf-1);
int minla=min(minla1,minla2);

if(_lf == 0)minla=0;

int maxla=_lf+2;

//-----------------------------------------------------
for(int la=minla; la<=maxla; la++)
      {
      complex<double> sumq2(czero);
      complex<double> sumq1(czero);
      complex<double> sumq (czero);
      int  maxiq=_nf-_lf-1;

      for(int iq=0; iq<=maxiq; iq++)
         {
         int n4 = _nf-_lf-1-iq;
         int n5 = iq;
         int n6 = 2*_lf+1+iq;
         double fac4=Factorial(n4);
         double fac5=Factorial(n5);
         double fac6=Factorial(n6);
         double xq =pow( (-2.*zp/_nf),iq)/fac4/fac5/fac6;
         complex<double> sump1(czero);
         complex<double> sump2(czero);
         complex<double> sump (czero);
         int maxip=_n0-_l0-1;

          for(int ip=0; ip<=maxip; ip++)
                {
                 int n40=_n0-_l0-1-ip;
                 int n50=ip;
                 int n60=2*_l0+1+ip;
                 double fac40=Factorial(n40);
                 double fac50=Factorial(n50);
                 double fac60=Factorial(n60);
                 double xp = pow( (-2*zp/_n0),ip)/fac40/fac50/fac60;
                //--------------------------------------------------

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                //--------------------------------------------------
                 int nn3 = 3+_lf+_l0+la+iq+ip;
                 int nn2 = 2+_lf+_l0+la+iq+ip;
                 double fn2 = Factorial(nn2);
                 double xmu =  2.5+_lf+_l0+iq+ip;
                 double xnu =  la+0.5;
                 double xa = (xnu+xmu)/2;
                 double xb= xa+0.5;
                 double xc= xnu+1.;
                 double xxb= xc-xb;
                 double xxz= gama/(gama+pow2(znz));
                 double hyp31 =  HyperGeom_function(xa,xxb,xc,xxz);
                        hyp31 /= pow(sqrt(gama+pow2(znz)),nn3);

       if(la == _lf+1 || la == abs(_lf-1)) sump1 += xp*fn2*hyp31;
       else {
                double xxa=xa-0.5;
                double xxb=xc-xa ;
                double hyp32 = HyperGeom_function(xxa,xxb,xc,xxz);
                       hyp32 /= pow(sqrt(gama+pow2(znz)),nn2);
                sump += xp*fn2*( ip*hyp32/nn2 - zp/_n0*hyp31 );

                if(la == _lf)
                      sump2 += xp*fn2*((ip+3)*hyp32/nn2 - zp/_n0*hyp31 );
                }
      }  // end do

   if(la == _lf+1 || la == abs(_lf-1))    sumq1 += xq*sump1;
   else         {                         sumq  += xq*sump;
                if(la == _lf)             sumq2 += xq*sump2;
                }
   }  // end do
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//                           s=-1 (indica el conjugado)
      double s=-1.;
      double ga=sqrtPI;

      for(int j=0; j<=la; j++) ga *= (0.5+j);

      double xkz=-epsi;
      double tk=acos(xkz/xk);
      calga=calfa/gama;

      double xla=sqrtPI*pow(xk,la)/pow(2.,(la+1))/ga  *sqrt((2*la+1.)/(2*_lf+1.)/4./PI);

      complex<double> cylam,cylam1,cylam2,cylam3;

      if(la == _lf+1 || la == abs(_lf-1))
         {
         double xc1 = ClebschGordan(la,1,_lf,   0,  0,  0);
         double xc2 = ClebschGordan(la,1,_lf,_mf+1,-1,_mf);
         ARMONIC(s,la,_mf+1,tk,0.,cylam1);
         csum1 += pow( ci,la)*xla*sqrt(3.)*sumq1*xc1*xc2*cylam1;
      }
 else {
         double xc3 = ClebschGordan(la,2,_lf,  0,     0,   0);
         double xc4 = ClebschGordan(la,2,_lf, _mf,    0, _mf);
         double xc5 = ClebschGordan(la,2,_lf, _mf+2, -2, _mf);
         double xc6 = ClebschGordan(la,2,_lf, _mf+1, -1, _mf);
         ARMONIC(s,la,_mf,   tk, 0., cylam );
         ARMONIC(s,la,_mf+2, tk, 0., cylam2);
         ARMONIC(s,la,_mf+1, tk, 0., cylam3);

         csum2+= pow( ci,la)*xla*sumq*xc3   *(-xc4*cylam/sqrt(6.)+xc5*cylam2);
         csum3+= pow( ci,la)*xla*sumq*xc3   *xc6*cylam3;
          if(la == _lf)
                      csum4=pow(ci,_lf)*xla*sumq2/sqrt(6.)*cylam;
         }
 }  // end do


      csum1*= calfa*(ca2/cinu*chyp1+chyp2);
      complex<double> ct = ci*eta*calga*ca2/cinu*chyp1;
      csum2*= ct;
      csum3*= (2*pv*chyp2-ci*epsi*calga*ca2/cinu*chyp1);
      csum4*= ct;
//-----------------------------------------------------------------------
      complex<double> csum= 4*PI*( csum1 - 2. *(csum2+csum3+csum4)  );
//-----------------------------------------------------------------------
      complex<double> cf = ci/(2.*PI*pv) * con0 * csum;
//     *   * cdexp(2.*cinu*dlog(gama))  factor de módulo=1
      double xf = eta*pow2(abs(cf));
//...............................
//     SC and ASC corrections
       if (icor == 1) {
                      double sq=pow2I(1.+pow2(sc/xk))+scoa*(1.-pow2I(1.+pow2(xk/sc)))/zt;
                      xf *= sq;
                      }
//...............................
return xf;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//***********************************************************************
//     Subrutinas.for
//***********************************************************************
//     Armonicos esfericos
//     Autor: Cesar Ramirez          fecha: 10-4-2000
//     l,/m/:   numeros enteros
//     tk     angulo teta
//     s  +1: real  -1: imaginario
//     ma  valor absoluto de m
//-----------------------------------------------------------------------
void ARMONIC(double s, int l, int m, double tk, double fik, complex<double>& cylm)
{
int ma=abs(m);
cylm = complex<double>(0.,0.);

if(ma > l)return;

int n1=l-ma;
int n2=l+ma;

double fac11=Factorial(n1);
double fac12=Factorial(n2);

double x1 = sqrt((2.*l+1.)/(4.*PI)*fac11/fac12);
double x  = cos(tk);

double xlm = PlnmLegendr(l,ma,x);

complex<double> t = s*ci*double(ma)*fik;

if(m >= 0)      cylm=x1*exp(t)*xlm;
else            cylm=pow_m1(ma)*x1*exp(-t)*xlm;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//***********************************************************************
//     Legendre polynomials (Numerical Recipes Software: www.nr.com)
//     x=cos(t)
//     ------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double PlnmLegendr(int L, int M, double X)
{
double plgndr=0;

if(M < 0 || M > L || fabs(X) > 1.)
      {
      //printf(*,*)' m<0,   or   m>l,   or   fabs[x]>1';
      //pause 'bad arguments';
      return BadValue;
      }

double PMM=1.;
if(M > 0)
        {
        double SOMX2=sqrt((1.-X)*(1.+X));
        double FACT=1.;
        for(int I=1; I<=M; I++)
              {
              PMM = -PMM*FACT*SOMX2;
              FACT +=2.;
              }
        }

if(L == M) plgndr=PMM;
else    {
         double PMMP1=X*(2*M+1)*PMM;

         if(L == M+1)  plgndr=PMMP1;
         else {
                 double PLL=0;
                 for(int LL=M+2; LL<=L; LL++)
                      {
                      PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M);
                      PMM=PMMP1;
                      PMMP1=PLL;
                      }
                  plgndr=PLL;
                }
         }
return plgndr;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void gamac(complex<double>za, int ny, double test1, double test2, int& nt, complex<double> &logam)
{
//     ******************************************************************
//          auteur: r. gayet              date inconnue
//          version 2    (modif a. salin)      30/10/81
//     log de gamma(z) ecrit si ny=2
//     test1 est la difference minimum admise entre reel(z) et un pole
//     pour distinguer reel(z) du pole
//     test2 est la valeur minimum admise pour imag(z) quand la
//     difference entre reel(z) et un pole est inferieure a test1
//     nt=2 si z est un pole de gamma:dans ce cas logam=0 arbitrairement-
//     sinon nt=1

      complex<double> z,z2,c;
      double a[3],d, test10;
      int irz, ng;


#define PIss  0.918938533204673

//      equivalence[z][a[1]] ;
      z = za ;
      a[1] = z.real();
      a[2] = z.imag();

      logam  = complex<double>(0.,0.) ;
      irz    = a[1];
      test10 = test1*2. ;

       if(a[1]-test10 <= 0)
          {
          d=fabs(a[1] - double(irz)) ;
         if(d <= test1 )
               {
               if(fabs(a[2]) <= test2)
                        {
//                        printf("real(z)=%9.2g + %22.15g   imag(z)=,1pd22\n"
//                               "z est considere comme un pole de la fonction gamma\n",d,z.real(), z.imag());
                        nt=2 ;
                        goto L100 ;
                        }
                }
          }
//L1:
      nt = 1;
      ng = 10-irz ;
      if(ng > 0)
        {
        complex<double> c(ng,0.) ;
        z=z+c ;
        }
//L4:
      z2=z*z ;
      d=abs(z) ;

      if(d < 100)
            {
            logam=1./156./z2 - 691./360360. ;
            logam=logam/z2+1./1188. ;
            logam=logam/z2-1./1680. ;
            }
//L6:
   if(d < 1.e+04)
              {
              logam = logam/z2+1./1260. ;
              logam = logam/z2-1./360. ;
              }
//L7:
   if(d < 1.e+07)
              logam = (logam/z2+1./12.)/z ;

//L8:
      logam = logam+PIss-z+(z-0.5)*log(z) ;

      if(ng > 0)
        for(int i=1; i<=ng; i++) {
                z=z-cuno ;
                logam=logam-log(z) ;
                }
L100:

    if(ny == 2)
             printf("log de gamma %22.15e %22.15e = %22.15e = %22.15e\n",
             z.real(),z.imag(),logam.real(), logam.imag());

}

//***********************************************************************
//***********************************************************************
//  Rutina para calcular el factorial de n
//  n! = n(n-1)(n-2)...2 ,
//  0! = 1
//----------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double Factorial(int n)
{
double xn=n;
if(n == 0 || n == 1) return 1.;

if(n > 90) return sqrt(PI*(2*xn+1/3.))*pow(xn,xn)/exp(xn);

double xfac=0.;
for(int j=2; j<=n; j++)
        xfac += log(double(j));

return exp(xfac);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//*********************************************************************
//**************************** l o g f a c ******************************
//     calculates ln((n-1)!).
//     beware: fact_fac(n) = ln((n-1)!) = ln(gamma(n))
//-----------------------------------------------------------------------
void LogFac()
{
if(LogFacInit) return;
else LogFacInit = true;

fact_fac[1]=0. ;

for(int i=2; i<=_idf; i++)
        fact_fac[i]=fact_fac[i-1]+log(double(i-1));
}
//*********************************************************************
//.......... ClebschGordan Calcula los coeficientes de Clebsch-Gordan
//.......... Debe llamarse antes al LogFac.
//.......... logcle:  simbolos 3j
//.......... xcg:     coef de C-G
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double ClebschGordan(int l1, int l2, int l3, int m1, int m2, int m3)
{
LogFac();
double x3j = LogCle(l1,l2,l3,m1,m2,-m3);
double z=pow_m1(l1-l2-m3) * sqrt(2.*l3+1.)*x3j;
return z;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//**************************** l o g c l e ******************************
//         author:  a.salin     version 2  23/2/77 - 5/2/96

//     calculation of WIGNER's 3j  (definition of MESSIAH)
//     restriction: INTEGER MOMENTS ONLY.

//     before the first run of logcle, a to LogFac should be
//     performed. LogFac calculates ln((n-1)!) for n going from 1 to idf
//     and stores the results in the common fact.
//     for large values of the l's, increase idf.
//        input
//                l1, l2, l3, m1, m2, m3: momenta and their components.
//        output
//                z: 3j coefficient.
//***********************************************************************
double LogCle(int l1,int l2,int l3,int m1,int m2,int m3)
{
double ac[_idfL];

if(!LogFacInit) {printf("LogFac not called before logcle"); return BadValue;}

int i4=l1+l2+l3+2;
if(i4 > _idf)   {printf("LogCle: value of idf too small");   return BadValue;}

if(m1+m2+m3 != 0) return 0;

int  izmax=min(min(l1+l2-l3, l1-m1   ),l2+m2   ) +1 ;
int  izmin=max(max(       0, l2-l3-m1),l1+m2-l3) +1 ;

if(izmax < izmin) return 0;

int i1=l1+l2-l3+1 ;
int i2=l1-m1+1 ;
int i3=l2+m2+1 ;

double  abra=0.5*(fact_fac[i1]      +fact_fac[l3+l1-l2+1]  +fact_fac[l3+l2-l1+1]  -fact_fac[i4]
                 +fact_fac[l1+m1+1] +fact_fac[i2]          +fact_fac[i3]          +fact_fac[l2-m2+1]
                 +fact_fac[l3+m3+1] +fact_fac[l3-m3+1]);

int k1=l3-l2+m1+1;
int k2=l3-l1-m2+1 ;

double   gros=250. ;

for(int ii=izmin; ii<=izmax; ii++)
        {
        int i=ii-1 ;
        ac[ii]= fact_fac[i+1]   +fact_fac[i1-i]  +fact_fac[i2-i]
               +fact_fac[i3-i]  +fact_fac[k1+i]  +fact_fac[k2+i];
        if(ac[ii] < gros)  gros=ac[ii];
        }

double accu=0.;

double sig=pow_m1(izmin);

for(int ii=izmin; ii<=izmax; ii++)
        {
        sig = -sig ;
        ac[ii] -= gros ;
        accu   += sig*exp(-ac[ii]);
        }

double z=pow_m1(l1-l2-m3)*exp(abra-gros)*accu ;
return z;
}
//***********************************************************************

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void CF21D(complex<double> ca, complex<double> cb, complex<double> cc, double x, complex<double>& cf)
{
// Program: hypere.f
// Version: 1.6 (see readme_hyper) 24/03/1997
// Author: Pablof (pablof@cab.cnea.edu.ar)
// Hypergeometric function for CDW-EIS calculations.
// Double precision, see comment for V 1.4
// We assume that x is REAL

if(ca == czero || cb == czero) {  cf=cuno;  return; }

       if(x == 0.)             {  cf=cuno;  return; }
 else  if(x == 1.) {
                cf=exp(LogGammaFunc(cc)   -LogGammaFunc(cc-ca-cb)
                      -LogGammaFunc(cc-ca)-LogGammaFunc(cc-cb   ));
                      return;}
else    {
        extern void hypere(complex<double> ca, complex<double> cb, complex<double> cc, double x, complex<double> &cf);
        hypere(ca,cb,cc,x,cf);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double HyperGeom_function(double A, double B, double C, double X)  /// Oleg check format!
{
//       ====================================================
//       Purpose: Compute hypergeometric function F(a,b,c,x)
//       Input :  a --- Parameter
//                b --- Parameter
//                c --- Parameter, c <> 0,-1,-2,...
//                x --- Argument   ( x < 1 )
//       Output:  HF --- F(a,b,c,x)
//       Routines called:
//            (1) GAMMA for computing gamma function
//            (2) PSI for computing psi function
//       ====================================================

const double EL=0.5772156649015329;

bool L0= C == int(C)      && C     <  0.;
bool L1= 1.0-X < 1.0e-15  && C-A-B <= 0.;
bool L2= A == int(A)      && A     <  0.;
bool L3= B == int(B)      && B     <  0.;
bool L4= C-A == int(C-A)  && C-A   <= 0.;
bool L5= C-B == int(C-B)  && C-B   <= 0.;

if (L0 || L1) { printf("The hypergeometric series is divergent");   return BadValue;}

double EPS=1.0e-15;
double HF=1, R;
int NM, K=0;

if (X > 0.95) EPS=1.0e-8;

        if (X == 0.0 || A == 0.0 || B == 0.0)  {return HF;}
 else   if (1.-X == EPS && C-A-B > 0.)
           {
           double GC   = GAMMA_function(C);
           double GCAB = GAMMA_function(C-A-B);
           double GCA  = GAMMA_function(C-A);
           double GCB  = GAMMA_function(C-B);
           HF=GC*GCAB/(GCA*GCB);
           return HF;
           }
 else if (1.+X <= EPS && fabs(C-A+B-1.) <= EPS)
           {
           double G0=sqrtPI*pow(2.,-A);
           double G1 = GAMMA_function(C);
           double G2 = GAMMA_function(1.0 + A/2.-B);
           double G3 = GAMMA_function(0.5+  A/2.);
           HF=G0*G1/(G2*G3);
           return HF;
           }
 else if (L2 || L3)
           {
           if (L2) NM=int(fabs(A));
           if (L3) NM=int(fabs(B));

           HF = 1.;
           R  = 1.;

           for(K=1; K<=NM; K++)
                {
                R *=  (A+K-1.)*(B+K-1.)/(K*(C+K-1.))*X;
                HF += R;
                }
           return HF;
           }
 else if (L4 || L5)
           {
           if (L4) NM=int(fabs(C-A));
           if (L5) NM=int(fabs(C-B));

           HF = 1.0;
           R  = 1.0;
           for(K=1; K<=NM; K++)
                {
                R *=  (C-A+K-1.)*(C-B+K-1.)/(K*(C+K-1.))*X;
                HF += R;
                }

           HF = pow(1.-X, C-A-B) * HF;
           return HF;
           }

double AA=A;
double BB=B;
double X1=X;

if (X < 0.)
          {
          X = X/(X-1.);

          if (C > A && B < A && B > 0.0) { A=BB; B=AA;}
          B=C-B;
          }

if (X >= 0.75)
         {       //-----------------------------------------------------
         double GM=0.;
         if (fabs(C-A-B-int(C-A-B)) < 1.0e-15)
                {
                int M = int(C-A-B);
                double GA = GAMMA_function(A);
                double GB = GAMMA_function(B);
                double GC = GAMMA_function(C);
                double GAM = GAMMA_function(A+M);
                double GBM = GAMMA_function(B+M);
                double PA = PSI_function(A);
                double PB = PSI_function(B);

              if (M != 0) GM=1.;
              for(int  J=1; J<=abs(M)-1; J++) GM *= J;

              double RM=1.;

              for(int J=1;  J<=abs(M);   J++) RM *= J;

              double F0=1.;
              double R0=1.;
              double R1=1.;
              double SP0=0.;
              double SP=0.;

              if (M >= 0)
                        {
                        double C0 = GM*GC/(GAM*GBM);
                        double C1 =-GC*pow((X-1.),M)/(GA*GB*RM);

                        for(K=1; K<=M-1; K++)
                            {
                            R0 *=  (A+K-1.)*(B+K-1.)/(K*(K-M))*(1.-X);
                            F0 +=  R0;
                            }

                        for(K=1; K<=M; K++)
                                SP0+= (1./(A+K-1.)+1./(B+K-1.)-1./K);

                        double F1=PA+PB+SP0+2.*EL+log(1.-X);
                        double HW = 0;

                        for(K=1; K<=250; K++)   //L55
                                {
                                   SP=SP+(1.-A)/(K*(A+K-1.))+(1.-B)/(K*(B+K-1.));

                                   double SM=0.;
                                   for(int J=1; J<=M; J++)
                                                  SM += ( (1.-A)/((J+K)*(A+J+K-1.))+1./(B+J+K-1.));

                                double RP =  PA+PB+2.*EL+SP+SM+log(1.-X);
                                R1 *= (A+M+K-1.)*(B+M+K-1.0)/(K*(M+K))*(1.0-X);
                                F1 += F1+R1*RP;
                                if (fabs(F1-HW) < fabs(F1)*EPS) break;
                                HW=F1;
                                }       //L55

                        HF=F0*C0+F1*C1;
                        }
              else if (M < 0)
                        {
                         M = -M;
                         double C0 = GM*GC/(GA*GB*pow((1.-X),M));
                         double C1 = -pow_m1(M)*GC/(GAM*GBM*RM);

                         for(int K=1; K<=M-1; K++) {
                            R0 *= (A-M+K-1.)*(B-M+K-1.)/(K*(K-M))*(1.-X);
                            F0 += R0;
                            }

                         for(int K=1; K<=M; K++)  SP0 += 1./K;

                         double F1=PA+PB-SP0+2.*EL+logl(1.-X);
                         double HW = 0;

                         for(K=1; K<=250; K++)
                                {
                                SP += ((1.-A)/(K*(A+K-1.))+(1.-B)/(K*(B+K-1.)));
                                double SM=0.;

                                for(int J=1; J<=M; J++) SM += 1./(J+K);

                                double RP = PA+PB+2.*EL+SP-SM+logl(1.-X);
                                R1 *= (A+K-1.)*(B+K-1.0)/(K*(M+K))*(1.0-X);
                                F1 += R1*RP;
                                if (fabs(F1-HW) < fabs(F1)*EPS) break;
                                HW=F1;
                                }
                         HF=F0*C0+F1*C1;
                         }
         }       //-----------------------------------------------------
    else {
              double GA   = GAMMA_function(A);
              double GB   = GAMMA_function(B);
              double GC   = GAMMA_function(C);
              double GCA  = GAMMA_function(C-A);
              double GCB  = GAMMA_function(C-B);
              double GCAB = GAMMA_function(C-A-B);
              double GABC = GAMMA_function(A+B-C);

              double C0=GC*GCAB/(GCA*GCB);
              double C1=GC*GABC/(GA*GB)*pow((1.-X),(C-A-B));

              HF=0.;
              double R0=C0;
              double R1=C1;
              double HW=0;

              for(K=1; K<=250; K++)
                  {
                  R0 *= (A+K-1.)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X);
                  R1 *= (C-A+K-1.)*(C-B+K-1.0)/(K*(C-A-B+K))*(1.0-X);
                  HF += R0+R1;
                  if( fabs(HF-HW) < fabs(HF)*EPS ) break;
                  HW=HF;
                  }

             HF += C0+C1;
        }
       }
 else  {
        double A0=1.;
        if (C > A && C < 2.*A && C > B && C < 2.*B)
                {
                A0=pow((1.-X),(C-A-B));
                A=C-A;
                B=C-B;
                }
        HF=1.;
        R=1.;
        double HW = 0;
        for(K=1; K<=250; K++)
              {
              R *= (A+K-1.)*(B+K-1.)/(K*(C+K-1.))*X;
              HF += R;
              if (fabs(HF-HW) <= fabs(HF)*EPS) break;
              HW=HF;
              }

        HF *= A0;
        }

if (X1 < 0.)
          {
          X=X1;
          double C0 = 1./pow((1.-X),AA);
          HF *= C0;
          }

//A=AA; FORTRAN issue
//B=BB; FORTRAN issue

if (K > 120) printf("Warning! You should check the accuracy\n");

return HF;
}


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//       Purpose: Compute gamma function â(x)
//       Input :  x  --- Argument of â(x)
//                       ( x is not equal to 0,-1,-2,úúú)
//       Output:  GA --- â(x)
//       ==================================================
//       ==================================================
double GAMMA_function(double X)
{
const double G[27] = { 0.,
                    1.0e0,0.5772156649015329e0,
                    -0.6558780715202538e0, -0.420026350340952e-1,
                    0.1665386113822915e0, -0.421977345555443e-1,
                    -0.96219715278770e-2, 0.72189432466630e-2,
                    -0.11651675918591e-2, -0.2152416741149e-3,
                    0.1280502823882e-3, -0.201348547807e-4,
                    -0.12504934821e-5, 0.11330272320e-5,
                    -0.2056338417e-6, 0.61160950e-8,
                    0.50020075e-8, -0.11812746e-8,
                    0.1043427e-9, 0.77823e-11,
                    -0.36968e-11, 0.51e-12,
                    -0.206e-13, -0.54e-14, 0.14e-14, 1e-15};

double GA,Z;
double R=1.;

if (X == int(X))
         {
         if (X > 0.)
                  {
                  GA=1.;
                  int M1=X-1;
                  for(int K=2; K<=M1; K++)
                                                GA *= K;
                  }
         else     GA=1.0e+300; /// check nmuber!
        }
 else   {
        if(fabs(X) > 1.)
              {
              Z=fabs(X);
              int    M=int(Z);
              for(int K=1; K<=M; K++)
                                 R *= (Z-K);

              Z -= M;
              }
        else   Z=X;

        double GR=G[26];
        for(int  K=25; K<=1; K--)
                         GR=GR*Z+G[K];

        GA = 1./(GR*Z);
        if(fabs(X) > 1.)
                {
                GA *= R;
                if (X < 0.) GA = -PI/ (X*GA*sin(PI*X));
                }
        }

return GA;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//       ======================================
//       Purpose: Compute Psi function
//       Input :  x  --- Argument of psi(x)
//       Output:  PS --- psi(x)                           ???  the logarithmic derivative of the gamma function. ???
//       ======================================
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double PSI_function(double X)
{
double  XA=fabs(X);
double  EL=0.5772156649015329;
double  S=0.;
double PSI;
int N;

      if (X  == int(X) && X <= 0.0)               return 1.0e+300;
 else if (XA == int(XA))
                {
                N=XA;
                for(int K=1; K<=N-1; K++) S += 1./K;

                PSI = -EL + S;
                }
 else if (XA+0.5 == int(XA+0.5))
                {
                N=XA-0.5;
                for(int K=1; K<=N;   K++) S += 1./(2.*K-1.);

                PSI = -EL + 2.*S  -1.386294361119891;
                }
 else           {
                if (XA < 10.0)
                        {
                        N=10-int(XA);
                        for(int K=0; K<=N-1; K++) S += 1./(XA+K);

                        XA += N;
                        }

                 double X2=1.0e0/(XA*XA);
                 double A1=-0.8333333333333e-01;
                 double A2= 0.83333333333333333e-02;
                 double A3=-0.39682539682539683e-02;
                 double A4= 0.41666666666666667e-02;
                 double A5=-0.75757575757575758e-02;
                 double A6= 0.21092796092796093e-01;
                 double A7=-0.83333333333333333e-01;
                 double A8= 0.4432598039215686e0;
                 PSI = log(XA)-0.5/XA+X2*(((((((A8*X2+A7)*X2+A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1);
                 PSI -= S;
                }

if (X < 0) PSI -=  (PI*cotan(PI*X)  + 1./X);    // cot(PI*X)= cos(PI*X)/sin(PI*X)

return PSI;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//************************************************************
//     Integrador.for
//*************************************** fin  programa principal ******
//     INTEGRADOR POR SIMPSO en [ymin, infinito)
//     AJUSTANDO LA PRESICION EN CADA INTERVALO Y LUEGO
//     entre el ultimo intervalo Y LA SUMA TOTAL
//     Para funciones Reales
//      (10-8-94)                        autor:  Cesar Ramirez
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void INTDEF(double PRE, double Ymin, double Ymax, double (*cfun)(double a), double& cint)
{

      double preci=PRE;
      double x1=Ymin;
      double x2=Ymax;
      double csum1=0.;
      cint=0.;

L39:  int n=1;
      double csum=0.;
      double h=x2-x1;

      double cf1 = cfun(x1);
      double cf2 = cfun(x2);
      double cf12=cf1+cf2;

      double cfmc=0.;

L40:  csum1=csum;
      double d=h/n;
      csum=(cf12+2.*cfmc)*d/6.;

      for(int i=1; i<=n; i++) {
              double xm = h/n*(i-0.5)+x1;
              double cfm = cfun(xm);
              csum += 4.*cfm*d/6. ;
              cfmc += cfm;
              }

 if(fabs(csum-csum1) > preci*fabs(csum1))
        {
        n *= 2.;
        goto L40;
        }
 else   {
        x1=x2;
        x2=x1+h;
        cint += csum;
        }

 if(fabs(csum) > preci*fabs(cint))
                                 goto L39;
}
//***********************************************************************
double cotan(double x)
{
if(x==0) return BadValue;
return cos(x)/sin(x);
}
