#include "../win/e_myextern.h"
#include "../win/e_Constant.h"
#include <QtMath>

extern double pow1(double par);
extern double pow2(double par);
extern double pow3(double par);
extern double pow4(double par);
extern double pow5(double par);
extern double pow2I(double par);

extern double E_to_Beta(double E);
extern double E_to_Gamma(double E);
extern double Velocity_au(double E);

void SEIK(double Ecurr, double &ck, double &cl, double &cm, double &cn, double &ct, double &cs);
void REC(double Ecurr, double &srk,double &srls,double &srlp);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEIK(double Ecurr, double &ck, double &cl, double &cm, double &cn, double &ct, double &cs)
{
//c*********************************************************************
//c                      Version du 20/01/98, revisitee en 2012, avec :
//c            - calcul du nombre d'electrons modifie
//c            - possibilite d'avoir Qp et Zp differents
//c*********************************************************************
//c        Calcule les sections efficaces de capture non radiative
//c         dans l'approximation eikonale
//c    (formules 8 p 3295  Meyerhof et al. Phys. Rev A 32 (1985) 3291 )
//c*********************************************************************
      double oc[4];
      double SIG[31];
      double ceik[31][4];


      bool Flag34 = EtachaVersion >= etacha_v34;

      int N_jp = Flag34 ?  30 : 3;

      for(int j=0; j<31; j++) SIG[j]=0;  //   Oleg 04/21/19

      for(int j=0; j<=3; j++)
                {
                oc[j]=0;
                for(int i=0; i<=30; i++) ceik[i][j]=0;
                }


//     yo is the mean number of electrons in n=5
//     not yet used, but may be used one day...

      double s0=2.8e3;
      double yO=0.;
      double Qm=Zp-(y1s+yL+yM+yN+yO);
      double Qp=Qm;
//     .........
/*   good  solution start

      double alph= FineStructConst;
      double alph2 = alph*alph;


      double beta = E_to_Beta(Ecurr);
      double gam  = E_to_Gamma(Ecurr);

   //good  solution stop
   */

// old solution to get the same values  start

      double alph=1./137.036;
      double alph2 = alph*alph;
      double beta=sqrt(1.-pow2I(1.+Ecurr/931.5));
      double gam=1.+Ecurr/931.5;


// old solution to get the same values  stop

      double d2 =(gam-1.)/(gam+1.);
      double d1 = sqrt(d2);
      double gag=gam/(gam+1.);
      double eta = alph / beta;

      double y1=qMax(0.,y1s-1.);
      double y2=qMax(0.,yL-1.);
      double y3=qMax(0.,yM-1.);
      double y4=qMax(0.,yN-1.);
      double y5=qMax(0.,yO-1.);

//c****************************
//c calcul du nombre d'electrons par couche de la cible (version modifiee)
      int jm=3;
      double tk=1.;
      double tm=0, tl=0;
      oc[1]=2.;

        if(Zt < 29.) tm=Zt-10.;
        else         tm=18.;

        if(Zt < 11.)
              {
              tm=0.;
              tl=Zt-2.;
              jm=2;
              }
        else  tl=8.;

       if(Zt < 3)
              {
              tk=Zt-1.;
              tl=0.;
              jm=1;
              }
       else   tk=1.;


       if(Zt == 1) oc[1]=1.;
                               //commented in original:     tk=0.8

//..............................
      oc[2]=tl;
      oc[3]=tm;

 for(int jp=1; jp<=N_jp; jp++)    // L1
        {
        double zps, zts, Z2PFAC, z2q;

             if(jp == 1) zps =  Zp-0.15*y1;
        else if(jp == 2) zps = (Zp-0.85*y1s-0.35*y2)/2.;
        else if(jp == 3) zps = (Zp-y1s-0.85*yL-0.35*y3)/3.;
        else if(jp == 4) zps = (Zp-y1s-yL-0.85*yM-0.35*y4)/4.;
        else if(jp == 5) zps = (Zp-y1s-yL-yM-0.85*y4)/5.;
        else if(jp == 6) zps = (Zp-y1s-yL-yM-y4-0.85*y5)/6.;
        else if(jp >= 7) zps = Qm/jp;

              for(int jt=1; jt<=jm; jt++)   //L2
                        {
                            if(jt == 1) zts=Zt-0.3*tk;
                       else if(jt == 2) zts=(Zt-1.7-0.35*tl)/2.;
                       else if(jt == 3) zts=(Zt-8.8-0.35*tm)/3.;

                       if(jt == 1 && jp == 1) Z2PFAC=1.;
                       else                   Z2PFAC=1.16;

                       if(zps <= zts)
                              {
                              z1=zps;
                              z2=zts;
                              z2q=zts*Z2PFAC;
                              }
                         else {
                              z1=zts;
                              z2=zps;
                              z2q=zps*Z2PFAC;
                              }

                      double Ei = sqrt(1.-alph2*z2*z2);
                      double Ef = sqrt(1.-alph2*z1*z1);
                      double PM = eta*(Ef/gam-Ei)/alph2;
                      double ezp  = eta*z2q;
                      double ezp2 = ezp*ezp;
                      double pez  = PI*ezp;
                      double expez  = exp(pez);
                      double exzt   = -2.*ezp*atan(-PM/z2);
                      double zfac   = z1*z2/(z2*z2+PM*PM);
                      double obkfac = s0*128*PI*eta*eta*pow5(zfac)/(5.*gag*gam);
                      double eikfac = 2.*pez*exp(exzt)/(expez-1./expez);
                      double eik = 1.+5.*ezp*PM/(4.*z2)+5.*ezp2*PM*PM/(12.*z2*z2)+ezp2/6.;
                      double mag = -d2+5.*d2*d2/16.+5.*d2*gag*z2q/(8.*z2)+d2*ezp2/4.;
                      mag += 5.*d2*d2*ezp2/48.;
                      double orb  = 5.*PI*d1*alph*(z1+z2)*(1.-d2/2.)/18.;
                             orb -= 5.*d1*alph*z2*ezp*(1-d2/2.)/8.;
                             orb -= 5.*PI*d1*gag*alph*z1*z2q/(18.*z2);
                             orb += 5.*PI*d1*gag*gag*alph*z1*z2q*z2q/(28.*z2*z2);
                             orb -= 5.*PI*d1*gag*alph*(z1+z2-d2*z1)*z2q/(28*z2);

                      ceik[jp][jt]=jp*jp*oc[jt]*obkfac*eikfac*(eik+mag+orb);
                      }

        SIG[jp]=ceik[jp][1]+ceik[jp][2]+ceik[jp][3];
        }


//c*****************************************************************
      ck=SIG[1];
      cl=SIG[2];
      cm=SIG[3];
      cn=SIG[4];
// not used
// Oleg     double co=SIG[5];
// Oleg     double cp=SIG[6];
// ........
      ct=0.;
      for(int i=5; i<=N_jp; i++)  ct += SIG[i];

//*******formule de schlachter****************************
      double E1=1000.*Ecurr/((pow(Zt,1.25))*(pow(Qp,0.7)));
      double sig1  = (1.1e-8/(pow(E1,4.8)))*(1-exp(-0.037*(pow(E1,2.2))));
             sig1 *=(1-exp(-2.44e-5*(pow(E1,2.6))));
      double sigs=(sig1*(pow(Qp,0.5)))/(pow(Zt,1.8)) ;
//c*****************************************************************
      cs=sigs;
      }
//c*************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//void REC(double &srk, double &srls, double &srlp)
  void REC(double Ecurr, double &sk,  double &sl1,  double &slp)
{
//c                      Version du 3/11/93
//c*********************************************************************
//c        Calcule les sections efficaces de REC 1s,2s et 2p
//c         dans l'approximation type Bethe et Salpeter
//c     (formules 71.7 p 304 et 71.14 71.15 p 306 ; voir aussi
//c      formule 75.6 p 322)
//c*********************************************************************
//Oleg -- not used       double QM=Zp-(y1s+yL+yM+yN);

      double Tel=Ecurr * 511003.4 /931.5;

      double beta = E_to_Beta(Ecurr);
      double gam  = E_to_Gamma(Ecurr);

      double azk =(Zp-0.15*y1s)*FineStructConst;
      double azl =(Zp-0.85*y1s-0.35*yL)*FineStructConst;
      double azm =(Zp-y1s-0.85*yL-0.35*yM)*FineStructConst;
         //  BK  = -511003.4*((1.+(az/(1.-az**2.)**0.5)**2.)**(-0.5)-1.)
      double az  = azk;   double BK  = -511003.4*(1./sqrt(1.+pow2(az/    sqrt(1.-pow2(az)))) -1.);
             az  = azl;   double BL1 = -511003.4*(1./sqrt(1.+pow2(az/(1.+sqrt(1.-pow2(az)))))-1.);
             az  = azm;   double BL3 = -511003.4*(1./sqrt(1.+pow2(az/    sqrt(4.-pow2(az)))) -1.);

//Oleg -- not used       int IBK =BK+0.5;
//Oleg -- not used       int IBL1=BL1+0.5;
//Oleg -- not used       int IBL3=BL3+0.5;
//c****************************
       double gmabet511 = gam*beta*511003.4;

      double EK  = BK +Tel;
      double EL1 = BL1+Tel;
      double EL3 = BL3+Tel;
//Oleg -- not used      int   IEK  = EK+0.5;
//Oleg -- not used      int   IEL1 = EL1+0.5;
//Oleg -- not used      int   IEL3 = EL3+0.5;

      double zk  =sqrt(BK /Tel);
      double zl3 =sqrt(BL3/Tel);
      double zl1 =sqrt(BL1/Tel);


      double sphk = 256.*PI*3.5e3/3.*pow3(BK )/pow4(EK)                   *exp(-4.*zk* atan(1./zk ))/(1.-exp(-2.*PI*zk));

      double sphl1=2048.*PI*3.5e3/3.*pow3(BL1)/pow4(EL1)*(1.+3.*BL1/EL1) * exp(-8.*zl1*atan(1./zl1))/(1.-exp(-4.*PI*zl1));
      double sphl2=2048.*PI*3.5e3/9.*pow4(BL1)/pow5(EL1)*(3.+8.*BL1/EL1) * exp(-8.*zl1*atan(1./zl1))/(1.-exp(-4.*PI*zl1));
      double sphl3=4096.*PI*3.5e3/9.*pow4(BL3)/pow5(EL3)*(3.+8.*BL3/EL3) * exp(-8.*zl3*atan(1./zl3))/(1.-exp(-4.*PI*zl3));

             sk =Zt*sphk *pow2(EK /gmabet511);     // argument
             sl1=Zt*sphl1*pow2(EL1/gmabet511);     // argument
      double sl2=Zt*sphl2*pow2(EL1/gmabet511);
      double sl3=Zt*sphl3*pow2(EL3/gmabet511);

      slp = sl2 + sl3;                     // argument
//Oleg -- not used       double sl  = sl1 + slp;
      }
//c*************************************************************



