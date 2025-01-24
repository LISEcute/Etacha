//   SUB STOP96 (Z1,M1,Z2,EMAX1000,SE(1),ENUMB,SOLIDGAS)
//|=========================================================
//|      Variable definitions for STOP95                    |
//|==>   Z1 = ION ATOMIC NUMBER                             |
//|==>   MM1= ION ATOMIC MASS                               |
//|==>   M1 = ION ATOMIC WEIGHT (AMU)                       |
//|==>   Z2 = TARGET ATOMIC NUMBER                          |
//|==>   M2 = TARGET ATOMIC WEIGHT (AMU)                    |
//|==>   EMAX1000 = FINAL ION ENERGY FOR 1000 VALUES (KEV)  |
//|==>   RHODENSITY = TARGET DENSITY (G/CM3)                |
//|==>   ATDENSITY=TARGET DENSITY (ATOMS/CM3)               |
//|==>   Vfermi = (fermi VELOCITY OF SOLID) / V0            |
//|==>   SE = CALCULATED ELECTRONIC STOPPING (EV-A2)        |
//=========================================================

//====> The below arrays are defined and filled from the calling program.
//====> They are filled using the subroutine STOPCOEF.BAS
//#pragma  warn -pch
//#include <owl/window.h>
//#include <math.h>
//#define pow(x,y)  exp((y)*log(x))
#include "L_Loss/dll_eloss.h"
#include "L_Loss/dll_SCOEF.h"
#include "L_Loss/scoef95.h"

#include <QtMath>

#define RowFermi 131

double Z_HI(     int Z1,  double E, int Z2, int SOLIDGAS, int GasNumb);
double P_StopGas(int Z1,  double E, int Z2, int GasNumb);
double P_Stop(            double E, int Z2);
double Z_helium(          double E, int Z2, int SOLIDGAS, int GasNumb);
double FermiCorrections(  double E, int Z2);
double NuclearStopping(int zp, double Mp,  int zt, double Mt, double e);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

//double far _export ZieglerStopping(int Z1, double E, int Z2, int SOLIDGAS) {
double  ZieglerStopping(int Z1, double E, int Z2, int SOLIDGAS)
{

  int i, GasNumb=0;
  double SP=0;

  E*=1000.;
  //========= Determine if(GAS request is valid =================================

  if(SOLIDGAS!=0)
    {
      for( i=0; i<8; i++)  	if(Z2==GAS[i]) GasNumb=i+1;
      if(GasNumb==0) SOLIDGAS=0;     // Reset for( i=Solid target.
    }
  //==============================================================================
  //======================================== LOOP  Stopping values ============
  if(Z1==1) 				 //== Proton Stopping  Powers
    {
      if(SOLIDGAS==1)SP=P_StopGas(Z1,E,Z2,GasNumb);
      else		   SP=P_Stop(E,Z2);
    }
  else   if(Z1==2) SP=Z_helium(E,Z2,SOLIDGAS, GasNumb);
  else   if(Z1>2)  SP=Z_HI(Z1,E,Z2,SOLIDGAS,GasNumb);

  //	  SP*=10.;		     //This converts Stopping in eV/(1E15-cm2) to eV-A2

  return SP;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWW========== PROTON ELECTRONIC STOPPING POWERS ====================WWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double  P_Stop(double E, int Z2 )
{

  double PE, PE0, X, Ppower, SP, SL, SH;
  int Z2c = qMin (130,Z2);

  if(E < 1E4) {        //**** VELOCITY PROPORTIONAL STOPPING BELOW VELOCITY ** PE0 **.
      PE=PE0=10.;
      if(PE < E) PE=E;

      SL=(SCOEF[Z2c][9]*pow(PE,SCOEF[Z2c][10]))+SCOEF[Z2c][11]*pow(PE,SCOEF[Z2c][12]);
      SH=SCOEF[Z2c][13]/pow(PE,SCOEF[Z2c][14])*log((SCOEF[Z2c][15]/PE)+SCOEF[Z2c][16]*PE);
      SP=SL*SH/(SL+SH);

      if(E<PE0) { //**** Ppower IS THE POWER OF VELOCITY STOPPING BELOW PE0.
          Ppower=0.45;                        // Low Energy Stopping: S(E)**0.45
          if(Z2<7) Ppower-=0.1;		     // Z2=3-6 has low S(E)**0.35
          SP*=pow((E/PE0),Ppower);
        }
    }
  else  {
      X=log(E)/E;             //********* High Energy Stopping (6/87)
      SP=SCOEF[Z2c][17]+(SCOEF[Z2c][18]*X)+(SCOEF[Z2c][19]*X*X)+(SCOEF[Z2c][20]/X);
    }

  return SP;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//========== Proton Stopping Powers in GASES   ====================
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double P_StopGas(int Zt, double E, int Zp,  int GasNumb)
{

  double SP, PE, PE0, Ppower, X, SL, SH;
  int i;

  int Z2c = qMin (130,Zp);


  //if(Z1>2) goto P_StopGas1;
  //---- Use regular P_Stop if(Z1>2
  if(E<1.e4)  {
      i=GasNumb;                       //---- GasNumb is from 1 to 8
      if(Zt>1) i=GasNumb+8;
      if(Zt>2) i=GasNumb+16;
      //**** VELOCITY PROPORTIONAL STOPPING BELOW VELOCITY ** PE0 **.
      PE=PE0=10;
      if(PE < E) PE=E;

      SL= SCOEFGAS[i][4]*pow(PE,SCOEFGAS[i][5]) +
          SCOEFGAS[i][6]*pow(PE,SCOEFGAS[i][7]);

      SH=SCOEFGAS[i][8]/pow(PE,SCOEFGAS[i][9]) *
          log( (SCOEFGAS[i][10]/PE) + SCOEFGAS[i][11]*PE);

      SP=SL*SH/(SL+SH);
      if(E<PE0) { //**** Ppower IS THE POWER OF VELOCITY STOPPING BELOW PE0.
          Ppower=0.45;                        // Low Energy Stopping: S(E)**0.45
          if(Zp<7) Ppower=0.35;	            // Z2=3-6 has low S(E)**0.35
          SP*=pow( E/PE0, Ppower);
        }
    }
  else  {
      X=log(E)/E;             //********* Regular High Energy Stopping
      SP=SCOEF[Z2c][17]+(SCOEF[Z2c][18]*X)+(SCOEF[Z2c][19]*X*X)+(SCOEF[Z2c][20]/X);
    }

  return SP;
}
//=================================================================
//=======    HELIUM Stopping Powers    ============================
//---- VELOCITY PROPORTIONAL STOPPING BELOW KEV/AMU  >> HE0 <<.

double Z_helium(double E, int Z2, int SOLIDGAS, int GasNumb)
{

  double HE, A, B, HEH, HE0, SP;

  HE0=1.;                             //Units = keV/amu
  if(HE0 > E) HE=HE0;
  else 	 HE=E;

  B=log(HE);
  A=.2865+ B*( .1266    +
               B*(-.001429  +
                  B*( .02402   +
                      B*(-.01135   +
                         B*  .001475))));

  if(A>30) A=30;
  HEH=1.-exp(-A);
  //**** ADD Z1^3 EFFECT TO HE/H STOPPING POWER RATIO** HEH **.
  if(HE<1) HE=1.;

  B=7.6-log(HE);
  A=(1.+(.007+.00005*Z2)*exp(-B*B));
  HEH*=A*A;

  if(SOLIDGAS==1) SP=P_StopGas(2,HE,Z2,GasNumb);
  else		SP=P_Stop(HE, Z2);


  SP*=HEH*4.;
  if(E <=HE0 ) SP*=sqrt(E/HE0);		//*** CALC. HE VELOCITY PROPORTIONAL STOPPING

  return SP;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//========== HEAVY ION ELECTRONIC STOPPING POWERS ===================
//**** USE VELOCITY STOPPING for( i=(YRmin=VR/Z1**.67) <= 0.13
//**** OR for( i=VR <= 1.0
double Z_HI(int Z1, double E, int Z2,  int SOLIDGAS, int GasNumb)
{

  int Z2c = qMin (130,Z2);
  int Z1c = qMin (130,Z1);

  double  VR, YR, Q, A, LAMBDA0, LAMBDA1, L, VMIN, EEE, HIPOWER, V2, ZETA0, ZETA, SP;
  int J;
  double  YRmin=0.13;
  double  VRmin=1.0;

  double Vfermi= SCOEF[Z2c][7]; 		             // Vfermi of Target
  double  V=sqrt(E/25.)/Vfermi;                    // Relative Velocity
  double hv1;

  V2=V*V;
  if(V<1) VR=(3.*Vfermi/4.)*(1.+ V2*2./3. - V2*V2/15.);
  else    VR=V*Vfermi*(1.+1./(5.*V2));

  //**** SET YR = MAXIMUM OF (VR/Z1^.67),(VRmin/Z1^.67) OR YRmin.

  hv1 =  1./ pow(Z1,.6667);
  YR=VR    * hv1;      YR=qMax(YR, YRmin);
  A =VRmin * hv1;      YR=qMax(YR, A);

  A = -.803*pow(YR,0.3) + 1.3167*pow(YR,0.6) +.38157*YR + .008983*YR*YR;

  if(A>50) A=50;                                //Prevents Underflow
  Q=1-exp(-A);

  if(Q<0) Q=0;
  else if(Q>1) Q=1;
  //**** Q = IONIZATION LEVEL OF THE ION AT VELOCITY * YR *.
  //**** NOW WE CONVERT IONIZATION LEVEL TO EFFECTIVE CHARGE.
  //====================================================================

  //---------------------- Screening Distance of Ion (Lambda in B.& K.)
  // Lambda is in SCOEF[Z2,22-39), Interpolation values in SCOEF[RowFermi,22-39)
  for( J=22; J<=39; J++)                      // Find Q Interpolation
    if(Q<=SCOEF[RowFermi][J]) break;		      //   in SCOEF[RowFermi,22-39)


  J--;

  if(J<22)      J=22;
  else if(J>38) J=38;    // Boilerplate

  LAMBDA0=SCOEF[Z1c][J];
  LAMBDA1=(Q-SCOEF[RowFermi][J])*(SCOEF[Z1c][J+1]-SCOEF[Z1c][J])/(SCOEF[RowFermi][J+1]-SCOEF[RowFermi][J]);
  L=(LAMBDA0+LAMBDA1)/ pow(Z1,.33333);

  //====================================================================
  double temp=4.*L*Vfermi/1.919;
  ZETA0=Q+(1./(2.*Vfermi*Vfermi))*(1.-Q)*log(1. + temp*temp );
  // old      ZETA0=Q+(1./(2.*Vfermi*Vfermi))*(1.-Q)*log(1+pow(4*L*Vfermi/1.919,2) );

  //**** ADD Z1^3 EFFECT AS SHOWN IN REF. 779.

  A=log(E); if(A<0) A=0;

  ZETA=ZETA0*(1.+(1./Z1/Z1)*(.08+.0015*Z2)*exp(-(7.6-A)*(7.6-A)));
  ZETA *= Z1;        // ZETA=(ZETA*Z1)*(ZETA*Z1)
  ZETA *= ZETA;

  A=VRmin*hv1;   if(A<YRmin) A=YRmin;

  if( YR > A)
    {

      if(SOLIDGAS==1) SP=P_StopGas(Z1,E,Z2,GasNumb);
      else 		    SP=P_Stop(E,Z2);

      SP*=FermiCorrections(E,Z2);
    }
  else       //**** CALCULATE VELOCITY STOPPING FOR YR LESS THAN YRmin.
    {
      VRmin=qMax(VRmin,YRmin/hv1);
      A=qMax(VRmin*VRmin-0.8*Vfermi*Vfermi, 0.);

      VMIN=.5*(VRmin+sqrt(A));
      EEE=25*VMIN*VMIN;

      if(SOLIDGAS==1) SP=P_StopGas(Z1,EEE,Z2,GasNumb);
      else 		      SP=P_Stop(EEE,Z2);

      SP*=FermiCorrections(EEE,Z2);

      //=======================================================================
      //   Following corrects for( i=low-energy stopping, where little data exists.
      //   Traditionally, this is velocity-proportional stopping, however for
      //   light ions, light targets and semiconductors, a large variation exists.
      //
      //======= Note: HIPOWER down = Se up ==== LAMBDA down = Se down =========
      HIPOWER=.47;                                       //TRIM-88 used 0.50
      if(Z1==3) HIPOWER= 0.55;
      else if(Z2<7 ) HIPOWER= 0.375;
      else                 //======> Following compensates for( i=semiconductor band-gap
        if((Z1<18) && (Z2==14 || Z2==32)) HIPOWER=0.375;

      SP*= pow(E/EEE,HIPOWER);
    }

  SP*=ZETA;
  return SP;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double FermiCorrections(double E, int Z2)
{
  int Z2c = qMin (130,Z2);

  int J;	          //=========== Add fermi Velocity Correction - 1995 ========================
  // VFCORR is in SCOEF[Z2,41-54), Interpolation values in SCOEF[RowFermi,41-54)
  double EION=qMin(E,9999.);
  // Not valid >1E4 keV/amu
  for( J=41; J<=53; J++)                 	// Find E Interpolation
    if(EION<SCOEF[RowFermi][J]) break;		      //   in SCOEF[RowFermi,41-54)

  J--;
  if(J<41)J=41;
  else if(J>53)  J=53;    // Boilerplate

  double VFCORR0=SCOEF[Z2c][J];
  double VFCORR1=(EION-SCOEF[RowFermi][J])*(SCOEF[Z2c][J+1]-SCOEF[Z2c][J]) /
      (SCOEF[RowFermi][J+1]-SCOEF[RowFermi][J]);

  return VFCORR0+VFCORR1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//double far _export NuclearStopping(int zp, double Mp,  int zt, double Mt, double e)  //TOI
double  NuclearStopping(int zp, double Mp,  int zt, double Mt, double e)  //TOI
{
  double Sne, SN, Zt, Zp, sq, ReducedIonEnergy, znam;

  Zp=zp; 	Zt=zt;

  double Ekev= e * Mp * 1000.;   // ion energy in KeV

  sq=  (Mp + Mt) * ( pow(Zp, 0.23) + pow(Zt, 0.23));

  ReducedIonEnergy = 32.53 * Mt * Ekev / (Zp * Zt) / sq ;

  if(ReducedIonEnergy <= 30.) {  //  TOI
      znam= ReducedIonEnergy +
          0.01321* pow(ReducedIonEnergy, 0.21226) +
          0.19593* pow(ReducedIonEnergy, 0.5);

      Sne=log(1.+1.1383*ReducedIonEnergy) / 2. / znam;
    }
  else   Sne=log(ReducedIonEnergy) / 2./ ReducedIonEnergy;


  SN = 8.462*Zp*Zt*Mp*Sne / sq;


  SN *=0.6022 / Mt;      // MeV/ (mg/cm2)

  return SN;  // MeV/ (mg/cm2)
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//double far _export Get_SCOEF(int i, int k)
double Get_SCOEF(int i, int k)
{
  if(i>130)i=130;
  return SCOEF[i][k];
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//double far _export Get_SCOEFGAS(int i, int k)
double Get_SCOEFGAS(int i, int k)
{
  return SCOEFGAS[i][k];
}


