#include <QtMath>
#include "../win/e_myextern.h"
#include "../win/e_Constant.h"


extern double pow1(double par);
extern double pow2(double par);
extern double pow3(double par);
extern double pow4(double par);
extern double pow5(double par);
extern double pow2I(double par);
extern double pow_int(double par, int power);
extern double E_to_Beta(double E);
extern double ANINT_O (double x);

double BTot(double z1, double rk, double rl1, double rl2, double rm1, double rm2, double rn);
double enel(int i);
void SEKE (double *WK, int IMAX, double *DCKE);
void SEL1E(double *WL1,int IMAX, double *DCL1E);
void SEL2E(double *WL2,int IMAX, double *DCL2E);
void SEM1E(double *WM1,int IMAX, double *DCM1E);
void SEM2E(double *WM2,int IMAX, double *DCM2E);
void SENE (double *WN, int IMAX, double *DCNE);

double F1S(double Q);
double F2S(double Q);
double F2P(double Q);
double FM1(double Q);
double FM2(double Q);
double FN (double Q);

double QUAD1(double *SOM, double *T, int IMAX);
double QUAD2(double *SOM2, double *Q);

const double Tol0=-1e-3;
const double Tol1= 1e-3;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PION(double EF)
{
//------------------------------------------------PION.FOR-----------------
// SECTIONS EFFICACES D'IONISATION PWBA                Version du 21/04/94
// d'un projectile  rk,rl1,rl2,rm1 et rm2          -----------------------
// electrons K,2s,2p,3s,p et 3d                  DEBUT DU PROGRAMME PRINCIPAL.
//-------------------------------------------------------------------------

      double DCKE[81],DCL1E[81],DCL2E[81],DCM1E[81], DCM2E[81],DCNE[81],
             WK[81],WL1[81],WL2[81],WM1[81],WM2[81],WN[81],
             sk[81],sl1[81],sl2[81],sm1[81],sm2[81],sn[81],
             T[81];

bool Flag34 = EtachaVersion >= etacha_v34;
bool Flag45 = EtachaVersion >= etacha_v45;

        for(int i=0; i<81; i++)
               { DCKE[i]=DCL1E[i]=DCL2E[i]=DCM1E[i]=DCM2E[i]=DCNE[i]=0;
                 WK[i]=WL1[i]=WL2[i]=WM1[i]=WM2[i]=WN[i]=sk[i]=sl1[i]=sl2[i]=sm1[i]=sm2[i]=sn[i]=0;}


      z1=Zp;
      z2=Zt;
      epo=EF;
      double rk  = y1s;
      double rl1 = y2s;
      double rl2 = y2p;
      double rm1 = y3s+y3p;
      double rm2 = y3d;
      double rn  = yN;
      double rk0 =qMax(rk,1.);
      double rl10=qMax(rl1,1.);
      double rl20=qMax(rl2,1.);
      double rm10=qMax(rm1,1.);
      double rm20=qMax(rm2,1.);
      double rn0 =qMax(rn,1.);
      double rkm  = rk0-1.;
      double rl1m = rl10-1.;
      double rl2m = rl20-1.;
      double rm1m = rm10-1.;
      double rm2m = rm20-1.;
      double rnm  = rn0-1.;


      o_Bk=  BTot(z1,rk0,rl1,rl2,rm1,rm2,rn)
           -BTot(z1,rkm,rl1,rl2,rm1,rm2,rn);

      o_Bl1=  BTot(z1,rk,rl10,rl2,rm1,rm2,rn)
             -BTot(z1,rk,rl1m,rl2,rm1,rm2,rn);

      o_Bl2=  BTot(z1,rk,rl1,rl20,rm1,rm2,rn)
             -BTot(z1,rk,rl1,rl2m,rm1,rm2,rn);

      o_Bm1=  BTot(z1,rk,rl1,rl2,rm10,rm2,rn)
             -BTot(z1,rk,rl1,rl2,rm1m,rm2,rn);

      o_Bm2=  BTot(z1,rk,rl1,rl2,rm1,rm20,rn)
             -BTot(z1,rk,rl1,rl2,rm1,rm2m,rn);

      o_Bn=   BTot(z1,rk,rl1,rl2,rm1,rm2,rn0)
             -BTot(z1,rk,rl1,rl2,rm1,rm2,rnm);


       if (o_Bl1 <  3.4 ) o_Bl1=3.4;
       if (o_Bl2 <  3.4 ) o_Bl2=3.4;
       if (o_Bm1 <= 1.51) o_Bm1=1.51;
       if (o_Bm2 <= 1.51) o_Bm2=1.51;
       if (o_Bn  <= 0.85) o_Bn =0.85;


// ************** calcul des constantes *******************************
//      betal=(1.-(1./pow((1.+epo/931.5),2.)))* pow((137.036),2);
      double beta0 = E_to_Beta(epo) * FineStructConstInverse;
      double betal = beta0 * beta0;               //  beta_au ^2
      double cste = 3.519e4* z2*z2;
//-------- charges effectives ecrantage Slater --------
    if(EtachaVersion >= etacha_v34)                 // Oleg's modification
        {
        if(rk > 1.)   zk=z1-0.3;                 // data from common/par4/
        else          zk=z1;
        }
    else
        {
        if(rk == 2.)   zk=z1-0.3;                 // data from common/par4/
        else           zk=z1;
        }


      zl1=z1-( 0.8*rk+                      0.30*qMax(rl1-1.,0.));
      zl2=z1-(     rk+0.75*rl1             +0.35*qMax(rl2-1.,0.));
      zm1=z1-(     rk+0.85*(rl1+rl2)       +0.35*qMax(rm1-1.,0.));
      zm2=z1-(     rk+      rl1+rl2+rm1    +0.35*qMax(rm2-1.,0.));
      zn =z1-(     rk+      rl1+rl2+rm1+rm2+0.35*qMax(rn -1.,0.));
       if(zl1 <= 1.) zl1=1.;
       if(zl2 <= 1.) zl2=1.;
       if(zm1 <= 1.) zm1=1.;
       if(zm2 <= 1.) zm2=1.;
       if(zn  <= 1.) zn=1. ;

// -------------------------------------------------------
//     factors for screening and antiscreening corrections
//........................................................
//      vu=137.036*pow((1.-pow((1.+epo/931.5),(-2.))),(0.5));
      double vu = beta0;

      double zp1 =0.9354*zk;
      double zp21=0.9354*zl1/2.;
      double zp22=0.9354*zl2/2.;
      double zp31=0.9354*zm1/3.;
      double zp32=0.9354*zm2/3.;
      double zp4 =0.9534*zn /4.;
                                         // data from common/par3/

       if (vu > zp1)  coak=1.-pow4(zp1/vu);
       else           coak=0;

       if (vu > zp21) coal1=1.-pow4(zp21/vu);
       else           coal1=0;

       if (vu > zp22) coal2=1.-pow4(zp22/vu);
       else           coal2=0;

       if (vu > zp31) coam1=1.-pow4(zp31/vu);
       else           coam1=0;

       if (vu > zp32) coam2=1.-pow4(zp32/vu);
       else           coam2=0;

       if (vu > zp4)  coan =1.-pow4(zp4 /vu);
       else           coan=0;

// ----------------------------------------------------
//     l'ionisation en n=4 est calculee comme etant
//     l'ionisation 2p de z/2
// ----------------------------------------------------
      zn=zn/2.;                          // data from common/par1/ and /par2/

      double ipow2_zk  = 1./pow2(zk);
      double ipow2_zl1 = 1./pow2(zl1);
      double ipow2_zl2 = 1./pow2(zl2);
      double ipow2_zm1 = 1./pow2(zm1);
      double ipow2_zm2 = 1./pow2(zm2);
      double ipow2_zn  = 1./pow2(zn);
      double ipow2_2zn = 1./pow2(2.*zn);

      etak =betal * ipow2_zk;
      etal1=betal * ipow2_zl1;
      etal2=betal * ipow2_zl2;
      etam1=betal * ipow2_zm1;
      etam2=betal * ipow2_zm2;
      etan =betal * ipow2_zn;
                                        //     Bk from common/bind/

      tetak = 1.*o_Bk /(13.6058)*ipow2_zk;
      tetal1= 4.*o_Bl1/(13.6058)*ipow2_zl1;
      tetal2= 4.*o_Bl2/(13.6058)*ipow2_zl2;
      tetam1= 9.*o_Bm1/(13.6058)*ipow2_zm1;
      tetam2= 9.*o_Bm2/(13.6058)*ipow2_zm2;
      tetan =16.*o_Bn /(13.6058)*ipow2_2zn;

      double cstk  = cste * pow2(ipow2_zk);
      double cstl1 = cste * pow2(ipow2_zl1);
      double cstl2 = cste * pow2(ipow2_zl2);
      double cstm1 = cste * pow2(ipow2_zm1);
      double cstm2 = cste * pow2(ipow2_zm2);
      double cstn  = cste * pow2(ipow2_zn);

      double fs1=1.;

       if (Zt == 1. ) fs1=1.23;
       if (Zt == 2. ) fs1=1.526;
       if (Zt == 6. ) fs1=0.78;
       if (Zt == 7. ) fs1=0.85;
       if (Zt == 10.) fs1=1.04;
       if (Zt == 13.) fs1=0.58;
       if (Zt == 14.) fs1=0.59;
       if (Zt == 18.) fs1=0.68;
       if (Zt == 29.) fs1=0.672;
       if (Zt == 36.) fs1=0.61;
       if (Zt == 54.) fs1=0.535;

      double help1 = pow2(fs1)*1.277*pow(z2,(2./3.));
                                                             // common/par2/
      scfk  = help1 * ipow2_zk;
      scfl1 = help1 * ipow2_zl1;
      scfl2 = help1 * ipow2_zl2;
      scfm1 = help1 * ipow2_zm1;
      scfm2 = help1 * ipow2_zm2;
      scfn =  help1 * ipow2_zn;
// ...................................................................
//      correction energie de liaison
// ....................................................................
double ek=1.;
double el1=1.;
double el2=1.;
double em1=1.;
double em2=1.;
double en=1.;

if (ibin == 1)
     {
      double xk =2.*vu/(zk*tetak);
      double gk =(1.+5.*xk+7.14*pow2(xk)+4.27*pow3(xk)+0.947*pow4(xk))/pow((1.+xk),5.);

      ek = pow2(1.+Zt*gk/(zk*tetak));

      double xl1=4.*vu/(zl1*tetal1);

      double gl1=1.+9.*xl1+30.2*pow2(xl1)+66.8*pow3(xl1)+100.*pow4(xl1)+
           +94.1*pow(xl1,5.)+51.3*pow(xl1,6.)+15.2*pow(xl1,7.)+1.891*pow(xl1,8.);
      gl1=gl1/pow((1.+xl1),9.);

      el1 = pow2(1.+Zt*gl1/(zl1*tetal1));

      double xl2 =4.* vu /(zl2*tetal2);


      double gl2=1.+9.*xl2+34.7*pow2(xl2)+81.2*pow3(xl2)+112.*pow4(xl2)
                 +93.5*pow5(xl2)+46.6*pow(xl2,6.)+12.9*pow(xl2,7.)+1.549*pow(xl2,8.);

      gl2 = gl2/pow((1.+xl2),9.);
      el2 = pow2(1.+Zt*gl2/(zl2*tetal2));

      double xm1=18.*vu/(zm1*tetam1);

      double gm1=(1.+5.*xm1+7.14*pow2(xm1)+4.27*pow3(xm1)+0.947*pow4(xm1));
      gm1=gm1/pow5((1.+xm1));

      em1=pow2(1.+Zt*gm1/(zm1*tetam1));

      double xm2=18.*vu/(zm2*tetam2);
      double gm2=(1.+5.*xm2+7.14*pow(xm2,2.)+4.27*pow(xm2,3.)+0.947*pow(xm2,4.));
      gm2=gm2/pow5(1.+xm2);
      em2=pow2(1.+Zt*gm2/(zm2*tetam2));


    if(Flag45) {
        double xn=18.*vu/(zn*tetan);
        double gn=(1.+5.*xn+7.14*pow2(xn)+4.27*pow3(xn) +0.947*pow4(xn));
        gn=gn/pow5(1.+xn);
        en=pow2(1.+Zt*gn/(zn*tetan));
        }
      }


// --------------------------------------------------------------------
//   Table des energies cinetiques de l'electron
// ---------------------------------------------
      double TMAX1=2194.32*epo;
      double TMAX=qMax(TMAX1,4.* o_Bk);
      int IMAX;
// ---------------------------------------------
      if (TMAX <= 1000.)      IMAX=39;
 else if (TMAX <= 2000.)      IMAX=41;
 else if (TMAX <= 5000.)      IMAX=46;
 else if (TMAX <= 10000.)     IMAX=51;
 else if (TMAX <= 20000.)     IMAX=53;
 else if (TMAX <= 50000.)     IMAX=56;
 else if (TMAX <= 100000.)    IMAX=58;
 else if (TMAX <= 200000.)    IMAX=60;
 else if (TMAX <= 500000.)    IMAX=63;
 else if (TMAX <= 1000000.)   IMAX=65;
 else if (TMAX <= 2000000.)   IMAX=67;
 else if (TMAX <= 5000000.)   IMAX=70;
 else                         IMAX=73;

 for(int i=1; i<=IMAX; i++) T[i]=enel(i);

//*********** Energies transferees en unites reduites ******
for(int I=1; I<=IMAX; I++)  //L106
      {
      WK[I]=  ek*tetak      +T[I]/(13.6058)*ipow2_zk;
      WL1[I]=el1*tetal1/4.  +T[I]/(13.6058)*ipow2_zl1;
      WL2[I]=el2*tetal2/4.  +T[I]/(13.6058)*ipow2_zl2;
      WM1[I]=em1*tetam1/9.  +T[I]/(13.6058)*ipow2_zm1;
      WM2[I]=em2*tetam2/9.  +T[I]/(13.6058)*ipow2_zm2;

    if(Flag34)
      WN[I] =en *tetan/4.   +T[I]/(13.6058)*ipow2_zn;
      }


//      ****************************************************
//      *****   Calcul des DCS et des pertes d'energie *****
//      ****************************************************

                      SEKE(WK, IMAX, DCKE);
                      SEL1E(WL1,IMAX,DCL1E);
                      SEL2E(WL2,IMAX,DCL2E);
                      SEM1E(WM1,IMAX,DCM1E);
                      SEM2E(WM2,IMAX,DCM2E);

if(EtachaVersion >= etacha_v34) SENE (WN, IMAX,DCNE );

// --------traitement des donnees du calcul : ---------------------------
// ----------------------------------------------------debut de la boucle
      for(int i=1; i<=IMAX; i++)
        {
                                          // -----------------------------------
                                          // Sections efficaces differentielles:
                                          // -----------------------------------
        sk[i] = cstk *DCKE [i]/(etak *13.6058*pow2(zk));
        sl1[i]= cstl1*DCL1E[i]/(etal1*13.6058*pow2(zl1));
        sl2[i]= cstl2*DCL2E[i]/(etal2*13.6058*pow2(zl2));
        sm1[i]= cstm1*DCM1E[i]/(etam1*13.6058*pow2(zm1));
        sm2[i]= cstm2*DCM2E[i]/(etam2*13.6058*pow2(zm2));

    if(Flag34)
        sn[i] = cstn* DCNE [i]/(etan *13.6058*pow2(zn));
        }

//----------------------------------------------------------------
// Calcul de la section efficace totale:
// --------------------------------------
      o_TCk  = QUAD1(sk, T,IMAX);
      o_TCl1 = QUAD1(sl1,T,IMAX);
      o_TCl2 = QUAD1(sl2,T,IMAX);
      o_TCm1 = QUAD1(sm1,T,IMAX);
      o_TCm2 = QUAD1(sm2,T,IMAX);
    if(Flag34)
      o_TCn = QUAD1 (sn, T,IMAX);
    else o_TCn = 0;

//     write(16,*)'Bk Bl1 BM1 Bn =', Bk,Bl1,BM1,Bn
//     write(16,*) Tck,TCL1,TCm1,TCn
// ------------------------------------
// on retabli la valeur de zn
      zn=2.*zn;
//--------------------------------------
}
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//                                   FIN de BORN.FOR
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// -----------------------------------
// fonction d'integration sur les T(I):
// -----------------------------------
double QUAD1(double *SOM, double *T, int IMAX)
{
double quad1=0.;

      for(int I=1; I<=IMAX-1; I++)
              quad1 += (SOM[I+1]+SOM[I])*(T[I+1]-T[I])/2;

return quad1;
}
//-------------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double enel(int i)
       {
const double el[74] = {0.0,
                0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,
                6.0,7.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,
                70.0,80.0,90.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,
                600.0,700.0,800.0,900.0,1000.0,1500.0,2000.0,2500.0,3000.0,
                3500.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0,
                15000.0,20000.0,30000.0,40000.0,50000.0,70000.0,100000.0,
                150000.0,200000.0,300000.0,400000.0,500000.0,700000.0,1000000.0,
                1500000.0,2000000.0,3000000.0,4000000.0,5000000.0,7000000.0,
                10000000.0,20000000.0};

      return el[i];
      }

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// -----------------------------------
// fonction d'integration sur les Q(I):
// -----------------------------------
double QUAD2(double *SOM2, double *Q)
      {
      double quad2=0.;

      for(int I=1; I<=57; I++)
        quad2+= (SOM2[I+1]+SOM2[I])*(Q[I+1]-Q[I])/2.;

return quad2;
}

//-------------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double qred(int i)
{
const double red[59] = { 0,
                0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,
                6.0,7.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,
                70.0,80.0,90.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,
                600.0,700.0,800.0,900.0,1000.0,1500.0,2000.0,2500.0,3000.0,
                3500.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0,
                15000.0,20000.0,30000.0,40000.0,50000.0,70000.0,100000.0 };

return red[i];
}
// -----------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

double BTot(double z1, double rk, double rl1, double rl2, double rm1, double rm2, double rn)
{
      double rks =qMax(rk -1.,0.);
      double rl1s=qMax(rl1-1.,0.);
      double rl2s=qMax(rl2-1.,0.);
      double rm1s=qMax(rm1-1.,0.);
      double rm2s=qMax(rm2-1.,0.);
      double rns =qMax(rn -1.,0.);

      double B0=rk*pow2(z1-0.3125*rks);
      double B1=0.25*rl1*pow2(z1-0.8*rk-0.3*rl1s);
      double B2=0.25*rl2*pow2(z1-rk-0.75*rl1-0.35*rl2s);
      double B3=1./9.*rm1*pow2(z1-rk-0.85*(rl1+rl2)-0.35*rm1s);
      double B4=1./9.*rm2*pow2(z1-(rk+rl1+rl2+rm1)-0.35*rm2s);
      double B5=1./16.*rn*pow2(z1-(rk+rl1+rl2+rm1+rm2)-0.35*rns);
      double Btot=13.6058*(B0+B1+B2+B3+B4+B5);
      return Btot;
}
// -----------------------------------------
//                       Debut des routines de sections efficaces.
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSk(T)/dT                    pour les electrons 1S1/2
// -----------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEKE(double *WK, int IMAX, double *DCKE)
{
//      double WK[1],DCKE[1];
      double QP[59],FF[59];
      double Q, Qm, F0, Fn;

//      printf(6,*) '1s ionization';

      for(int i=1; i<59; i++) { QP[i]=FF[i]=0.;}

      for(int i=1; i<=IMAX; i++)  ///*L300
        {
        parf1s_W=WK[i];
        Qm = pow2(parf1s_W)/(4.*etak);

                //------ recherche du qmax de convergence -> pas d'integration------
        F0= F1S(Qm);
        Q=Qm;

L299:     Q=2.*Q;
          Fn=ANINT_O (100.*F1S(Q)/F0);
          if (Fn >= 1.) goto L299;

        double PQ=Q/1000.;
        //------------------------------------------------------------------
        for(int j=1; j<=39; j++)  //L301
                  {
                  QP[j]=qred(j)*PQ;
                  Q=Qm+QP[j];
                  FF[j]=F1S(Q);
                  } //L301:

        DCKE[i]=QUAD2(FF,QP);
        }

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// ~~~~~~~~~~~~~~~~~~~~~~~
// FONCTION F1S FACTEUR DE FORME POUR ELECTRON 1S
// ~~~~~~~~~~~~~~~~~~~~~~~
double F1S(double Q)
      {
      // Pi=4.*atan(1.);   from Constant.h
                                                        // common/par1/ epo,z1,z2,etak,etal1,etal2,etan;
                                                        // common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
                                                        // common/par3/coak,coal1,coal2,coam1,coam2,coan;
      double sq1 = 0;      if(Q > 1e-10   ) sq1 = pow2I(1.+scfk/Q);
      double sq2 = 0;      if(scfk > 1e-10) sq2 = coak*(1.-pow2I(1+Q/scfk))/z2;

      double sq = sq1 + sq2;
      double AK2=parf1s_W-1;

      double AS=Q+(AK2)/3.+(1./3.);
      double CS = 1;
      double ACC, BS;

if(AK2 < Tol0)
        {
        ACC=sqrt(fabs(AK2));
        BS=pow( (Q+pow2(1.-ACC)) / (Q+pow2(1.+ACC)) ,
                 1./ACC);
        }
else if(AK2 < 0)
        {
        ACC=sqrt(fabs(AK2));
        BS=exp(-4./(Q+pow2(1.+ACC)));
        }
 else if(AK2 < Tol1)
        {
        BS=exp(-4./(Q+1.-AK2));
        }
 else   {
        ACC=sqrt(AK2);
                                //ATAN2(Y,X) computes the arctangent of the complex number X + i Y.
        double ata = atan2(2.*ACC, Q-AK2+1.);
        BS=exp((-2./ACC)*ata);
        CS=(1.-exp((-2.*PI)/ACC));
        }

double DS=pow3(pow2(Q-AK2+1)+4.*AK2);
double ES=128;  //pow(2.,7.);
double F1S=sq*(ES*AS*BS)/(CS*DS*Q);
return F1S;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSl1(T)/dT                   pour les electrons 2S
// -----------------------------------------------------------------------
void SEL1E(double *WL1,int IMAX, double *DCL1E)
{
//      double WL1[1],DCL1E[1];
      double QP[59],FF[59];
      double Q, Qm, F0, Fn;

     // printf(6,*) '2s ionization';

     for(int i=1; i<59; i++) { QP[i]=FF[i]=0.;}

     for(int  i=1; i<=IMAX; i++)   ///*L400
        {
        parf2s_W1=WL1[i];
        Qm=pow2(parf2s_W1)/(4.*etal1);

                        //------ recherche du qmax de convergence -> pas d'integration------
        F0=F2S(Qm);
        Q=Qm;

L399:     Q=2.*Q;
          Fn=ANINT_O(100.*F2S(Q)/F0);
          if (Fn >= 1.) goto L399;


        double PQ=Q/1000.;
//------------------------------------------------------------------
        for(int j=1; j<=39; j++)   //L1401
                {
                QP[j]=qred(j)*PQ;
                Q=Qm+QP[j];
                FF[j]=F2S(Q);
                }

        DCL1E[i]=QUAD2(FF,QP);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// ~~~~~~~~~~~~~~~~~~~~~~~
// FONCTION F2S FACTEUR DE FORME POUR ELECTRON 2S
// ~~~~~~~~~~~~~~~~~~~~~~~
double F2S(double Q)
{

      double sq=pow2I(1.+scfl1/Q)+coal1*(1.-pow2I(1+Q/scfl1)) /z2;
      double AK2=parf2s_W1-0.25;

      double CS = 1;
      double ACC, BS;

if(AK2 < Tol0)
        {
        ACC=sqrt(fabs(AK2));
        BS=pow( (Q+pow2(0.5-ACC))/(Q+pow2(0.5+ACC)),
                 1./ACC);
        }
else if(AK2 < 0)
        {
        ACC=sqrt(fabs(AK2));
        BS=exp(-2./(Q+pow((0.5+ACC),2)));
        }
else if(AK2 < Tol1)
        {
        BS=exp(-2./(Q+0.25-AK2));
        }
else {
      ACC=sqrt(AK2);
      double ata = atan2(ACC,Q-AK2+0.25);
      BS=exp((-2./ACC)*ata);
      CS=(1.-exp((-2.*PI)/ACC));
      }

double AS=pow2(Q-AK2+0.25)+AK2;
double ES=16;
double AL1=BS*ES/(CS*pow5(AS)*Q);

double p2_AK2 = pow2(AK2);
double p3_AK2 = pow3(AK2);
double p4_AK2 = pow4(AK2);
double p5_AK2 = pow5(AK2);

double A5=                                                                                                        pow5(Q);
double A4=-(8./3.       +11./3. *AK2                                                                            )*pow4(Q);
double A3= (41./24.     +6.     *AK2      +14./3.*p2_AK2                                                        )*pow3(Q);
double A2= (5./48.      -31./24.*AK2      -10./3.*p2_AK2      -2. *p3_AK2                                       )*pow2(Q);
double A1= (47./3840.                   -41./120.*p2_AK2   -2./3. *p3_AK2  -1./3.*p4_AK2                        )*pow1(Q);
double A0=  1./768.     +17./768.*AK2   +7./48.  *p2_AK2  +11./24.*p3_AK2  +2./3.*p4_AK2 +1./3.*p5_AK2;

double F2S=sq*AL1*(A5+A4+A3+A2+A1+A0);

return F2S;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSl2(T)/dT                   pour les electrons 2P
// -----------------------------------------------------------------------
void SEL2E(double *WL2,int IMAX, double *DCL2E)
{
      double QP[59],FF[59];
      double Q, Qm, F0, Fn;

      //      printf(6,*) '2p ionization';

     for(int i=1; i<59; i++) { QP[i]=FF[i]=0.;}

      for(int i=1; i<=IMAX; i++)  //L500
        {
        parf2p_W2=WL2[i];
        Qm=pow2(parf2p_W2)/(4.*etal2);
                        //------ recherche du qmax de convergence -> pas d'integration------
        F0=F2P(Qm);
        Q=Qm;

L499:   Q=2.*Q;
        Fn=ANINT_O(100.*F2P(Q)/F0);
        if (Fn >= 1.) goto L499;

        double PQ=Q/1000.;
//------------------------------------------------------------------
        for(int j=1; j<=39; j++)   //L501
                {
                  QP[j]=qred(j)*PQ;
                  Q=Qm+QP[j];
                  FF[j]=F2P(Q);
                 }

        DCL2E[i]=QUAD2(FF,QP);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// ~~~~~~~~~~~~~~~~~~~~~~~~~
//    FONCTION F2P FACTEUR DE FORME POUR ELECTRON 2P
// ~~~~~~~~~~~~~~~~~~~~~~~~~
double F2P(double Q)
{
      double CS = 1;
      double ACC, BS;

      double sq=pow2I(1.+scfl2/Q)+coal2*(1.-pow2I(1+Q/scfl2))/z2;
      double AK2=parf2p_W2-0.25;

if(AK2 < Tol0) {
        ACC=sqrt(fabs(AK2));
        BS=pow( (Q+pow2(0.5-ACC))/(Q+pow2(0.5+ACC)),
                 1./ACC);
        }
 else if(AK2 < 0) {
        ACC=sqrt(fabs(AK2));
        BS=exp(-2./(Q+pow2(0.5+ACC)));
        }
 else if(AK2 < Tol1) {
      BS=exp(-2./(Q+0.25-AK2));
      CS=1.;
        }
 else   {
        ACC=sqrt(AK2);
        double ata = atan2(ACC,Q-AK2+0.25);
        BS=exp((-2./ACC)*(ata));
        CS=(1.-exp((-2.*PI)/ACC));
        }

double      AS=pow2(Q-AK2+0.25)+AK2;
double      ES=16;
double      AL2=BS*ES/(CS*pow5(AS)*Q);

double p2_AK2 = pow2(AK2);
double p3_AK2 = pow3(AK2);
double p4_AK2 = pow4(AK2);

double      A4=                                                            9./4.*pow4(Q);
double      A3=-(0.75          +3.*AK2                                         )*pow3(Q);
double      A2= (19./32.     -0.75*AK2      -0.5*p2_AK2                        )*pow2(Q);
double      A1= (107./960.+41./48.*AK2 +113./60.*p2_AK2 +       p3_AK2         )*pow1(Q);
double      A0=  11./3072  +3./64.*AK2 +  7./32.*p2_AK2+ 5./12.*p3_AK2 +1./4.*p4_AK2;

double      F2P=sq*AL2*(A4+A3+A2+A1+A0)*(1./3.);
return F2P;
}

//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSm1(T)/dT                   pour les electrons 3s,p
// -----------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEM1E(double *WM1,int IMAX,double *DCM1E)
{
      double QP[59],FF[59];
      double Q, Qm, F0, Fn;
      for(int i=1; i<59; i++) { QP[i]=FF[i]=0.;}

     //      printf(6,*) '3s,p ionization';

      for(int i=1; i<=IMAX; i++)  //L600
        {
        parfm1_W1M=WM1[i];
        Qm=pow2(parfm1_W1M)/(4.*etam1);
                        //------ recherche du qmax de convergence -> pas d'integration------
        F0=FM1(Qm);
        Q=Qm;

L599:     Q=2.*Q;
          Fn=ANINT_O(100.*FM1(Q)/F0);
          if (Fn >= 1.) goto L599;

        double PQ=Q/1000.;
//------------------------------------------------------------------
        for(int j=1; j<=39; j++)
                {
                QP[j]=qred(j)*PQ;
                Q=Qm+QP[j];
                FF[j]=FM1(Q);
                }

        DCM1E[i]=QUAD2(FF,QP);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// ~~~~~~~~~~~~~~~~~~~~~~~~~
// FONCTION Fm1 FACTEUR DE FORME POUR ELECTRON 3s,p
// ~~~~~~~~~~~~~~~~~~~~~~~~~
double FM1(double Q)
{
      double CS = 1;
      double ACC, BS;

      double sq=pow2I(1.+scfm1/Q) + coam1*(1.-pow2I(1+Q/scfm1))/z2;
      double AK2=parfm1_W1M-1/9.;

if(AK2 < Tol0)
        {
        ACC=sqrt(fabs(AK2));
        BS=pow(    (Q+pow2(1./3.-ACC)) / (Q+pow2(1./3.+ACC)),
               1./ACC);
        }
 else if(AK2 < 0)
        {
        ACC=sqrt(fabs(AK2));
        BS=exp(-4./3./(Q+pow2(1./3.+ACC)));
        }
 else if(AK2 < Tol1)
        {
        BS=exp(-4./3./(Q+1./9.-AK2));
        }
 else   {
        ACC=sqrt(AK2);
        double ata = atan2(2.*ACC/3.,Q-AK2+1./9.);
        BS=exp((-2./ACC)*(ata));
        CS=(1.-exp((-2.*PI)/ACC));
        }

double       AS=pow((Q-AK2+1./9.),2.)+4.*AK2/9.;
double       ES=4.74074074074074;                       //pow(2.,7.)/27.;
double       AL3=BS*ES/(Q*CS*pow5(AS));
double       A5=pow(Q,5.);

double p2_AK2 = pow2(AK2);
double p3_AK2 = pow3(AK2);
double p4_AK2 = pow4(AK2);
double p5_AK2 = pow5(AK2);

double       A4= -(43./27.       +11./3.*AK2                                                        )*pow4(Q);
double       A3=  (518./243.   +412./81.*AK2     +14./3.*p2_AK2                                     )*pow3(Q);
double       A2= -(442./729.   +310./81.*AK2   +122./27.*p2_AK2      +2.*p3_AK2                     )*pow2(Q);
double       A1=  (0.0987146  + 0.667581*AK2  +290./243.*p2_AK2  +4./27.*p3_AK2   -1./3.*p4_AK2     )*pow1(Q);
double       A0=   0.0021282  +0.0411692*AK2  +0.2728243*p2_AK2 +62./81.*p3_AK2 +71./81.*p4_AK2+ 1./3.*p5_AK2;

double       Fm1=sq*AL3*(A5+A4+A3+A2+A1+A0)*(1./9.);
return Fm1;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSm2(T)/dT                   pour les electrons 3d
// -----------------------------------------------------------------------
void SEM2E(double *WM2, int IMAX, double *DCM2E)
{
      double QP[59],FF[59];
      double Q, Qm, F0, Fn;
      for(int i=1; i<59; i++) { QP[i]=FF[i]=0.;}

//      printf(6,*) '3d ionization';

      for(int i=1; i<=IMAX; i++)
        {
        parfm2_W2M=WM2[i];
        Qm=pow2(parfm2_W2M)/(4.*etam2);
                        //------ recherche du qmax de convergence -> pas d'integration------
        F0=FM2(Qm);
        Q=Qm;

       if(F0 > 0.)
                {
L599:           Q=2.*Q;
                Fn=ANINT_O(100.*FM2(Q)/F0);
                if (Fn >= 1.) goto L599;
                }

        double PQ=Q/1000.;
//------------------------------------------------------------------
        for(int j=1; j<=39; j++)
                {
                QP[j]=qred(j)*PQ;
                Q=Qm+QP[j];
                FF[j]=FM2(Q);
                }

        DCM2E[i]=QUAD2(FF,QP);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// ~~~~~~~~~~~~~~~~~~~~~~~~~
// FONCTION Fm2 FACTEUR DE FORME POUR ELECTRON 3d
// ~~~~~~~~~~~~~~~~~~~~~~~~~
double FM2(double Q)
{
      double CS = 1;
      double ACC, BS;

      double sq=pow2I(1.+scfm2/Q) + coam2*(1.-pow2I(1+Q/scfm2))/z2;
      double AK2=parfm2_W2M-1/9.;

if(AK2 < Tol0)
        {
        ACC=sqrt(fabs(AK2));
        BS=pow(
                (Q+pow2(1./3.-ACC))/(Q+pow2(1./3.+ACC)),
                1./ACC);
        }
 else if(AK2 < 0)
        {
        ACC=sqrt(fabs(AK2));
        BS=exp(-4./3./(Q+pow2(1./3.+ACC)));
        }
 else if(AK2 < Tol1)
        {
        BS=exp(-4./3./(Q+1./9.-AK2));
        }
 else   {
        ACC=sqrt(AK2);
        double ata = atan2(2.*ACC/3.,Q-AK2+1./9.);
        BS=exp((-2./ACC)*(ata));
        CS=(1.-exp((-2.*PI)/ACC));
        }


double      AS=pow((Q-AK2+1./9.),2.)+4.*AK2/9.;
double      ES=4.74074074;
double      AL3=BS*ES/(Q*CS*pow5(AS));

double p2_AK2 = pow2(AK2);
double p3_AK2 = pow3(AK2);
double p4_AK2 = pow4(AK2);
double p5_AK2 = pow5(AK2);

double      A5=                                                                                   pow5(Q);
double      A4=  -(43./27.      +11./3.  *AK2                                                   )*pow4(Q);
double      A3=   (518./243.    +412./81.*AK2 +14./3.   *p2_AK2                                 )*pow3(Q);
double      A2=  -(442./729.    +310./81.*AK2 +122./27. *p2_AK2 +2.     *p3_AK2                 )*pow2(Q);
double      A1=   (0.09871463   +0.667581*AK2 +290./243.*p2_AK2 +4./27. *p3_AK2 -1./3.  *p4_AK2 )*pow1(Q);
double      A0=    0.002128176  +0.041169*AK2 +0.2728243*p2_AK2 +62./81.*p3_AK2 +71./81.*p4_AK2 + 1./3.*p5_AK2;

double      Fm2=sq*AL3*(A5+A4+A3+A2+A1+A0)*(1./9.);
return Fm2;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSN(T)/dT                    pour les electrons N -> "2P"
// -----------------------------------------------------------------------
void SENE(double *WN, int IMAX, double *DCNE)
{
      double QP[59],FF[59];
      double Q, Qm, F0, F4;
      for(int i=1; i<59; i++) { QP[i]=FF[i]=0.;}


//      printf(6,*) 'n=4 ionization';

      for(int i=1; i<=IMAX; i++)
        {
        parfn_W4 = WN[i];
        Qm=pow2(parfn_W4)/(4.*etan);
                        //------ recherche du qmax de convergence -> pas d'integration------
        F0=FN(Qm);
        Q=Qm;

L499:     Q=2.*Q;
          F4 = ANINT_O(100.*FN(Q)/F0);
          if (F4 >= 1.) goto L499;

        double PQ=Q/1000.;
//------------------------------------------------------------------
        for(int j=1; j<=39; j++)
                {
                QP[j]=qred(j)*PQ;
                Q=Qm+QP[j];
                FF[j]=FN(Q);
                }

        DCNE[i]=QUAD2(FF,QP);
        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
// ~~~~~~~~~~~~~~~~~~~~~~~~~
//    FONCTION FN FACTEUR DE FORME POUR ELECTRON 2P -> n=4
// ~~~~~~~~~~~~~~~~~~~~~~~~~
double FN(double Q)
{
      double CS = 1;
      double ACC, BS;

      double sq=pow2I(1.+scfn/Q) + coan*(1.-pow2I(1+Q/scfn))/z2;
      double AK2=parfn_W4-0.25;

if(AK2 < Tol0)
        {
        ACC=sqrt(fabs(AK2));
        BS=pow( (Q+pow2(0.5-ACC))/(Q+pow2(0.5+ACC)),
                 1./ACC);
          }
 else if(AK2 < 0)
        {
        ACC=sqrt(fabs(AK2));
        BS=exp(-2./(Q+pow2(0.5+ACC)));
        }
 else if(AK2 < Tol1)
        {
        BS=exp(-2./(Q+0.25-AK2));
        }
 else {
        ACC=sqrt(AK2);
        double ata = atan2(ACC,Q-AK2+0.25);
        BS=exp((-2./ACC)* ata);
        CS=(1.-exp((-2.*PI)/ACC));
        }

double      AS=pow2(Q-AK2+0.25)+AK2;
double      ES=16;
double      AL2=BS*ES/(CS*pow5(AS)*Q);

double p2_AK2 = pow2(AK2);
double p3_AK2 = pow3(AK2);
double p4_AK2 = pow4(AK2);

double      A4=                                                             9./4.*pow4(Q);
double      A3= -(0.75          +3.     *AK2                                    )*pow3(Q);
double      A2=  (19./32.       -0.75   *AK2   -0.5     *p2_AK2                 )*pow2(Q);
double      A1=  (107./960.     +41./48.*AK2   +113./60.*p2_AK2       +  p3_AK2 )*pow1(Q);
double      A0=   11./3072.     +3./64. *AK2   +7./32.  *p2_AK2  +5./12.*p3_AK2 +1./4.*p4_AK2;

double      FN=sq*AL2*(A4+A3+A2+A1+A0)*(1./3.);
return FN;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//                                  FIN DES ROUTINES DES SECTIONS EFFICACES.
//                           ----------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//           FIN DU PROGRAMME
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


