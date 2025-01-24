#include "../win/e_myextern.h"
#include "../win/e_Constant.h"
#include <QtMath>
//#include <stdio.h>

//extern double Velocity_au(double E);
//extern double E_to_Beta(double E);

//extern double pow1(double par);
extern double pow2(double par);


double SNRM2(int N, double *SX, int INCX);




/*
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int function ISAMAX[N][SX][INCX] {
      double SX[*],SMAX,XMAG;
      ISAMAX = 0;
       if(N <= 0) return;
      ISAMAX = 1;
       if(N <= 1)return;
       if(INCX == 1)goto L20;

      SMAX = fabs(SX[1]);
      NS = N*INCX;
      II = 1;
          for(int I=1; I<=NS; I+INCX) { //L10
          XMAG = fabs(SX[I]);
           if(XMAG <= SMAX) goto L5;
          ISAMAX = II;
          SMAX = XMAG;
L5:       II = II + 1;
L10:      } //continue;
      return;

L20:  SMAX = fabs(SX[1]);
      for(int I=2; I<=N; I++) { //L30
         XMAG = fabs(SX[I]);
          if(XMAG <= SMAX) goto L30;
         ISAMAX = I ;
         SMAX = XMAG;
L30:  } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function SASUM[N][SX][INCX] {

      double SX[*];
      SASUM = 0.0E0 ;
       if(N <= 0)return;
       if(INCX == 1)goto L20;

      NS = N*INCX;
          for(int I=1; I<=NS; I+INCX) {/// L10
          SASUM = SASUM + fabs(SX[I]);
L10:      } //continue;
      return;

L20:  M = ((N)%(6));
       if( M  ==  0 ) goto L40 ;
      for(int I=1; I<=M; I++) {/// L30
        SASUM = SASUM + fabs(SX[I]);
L30:  } //continue;
       if( N  <  6 ) return;
L40:  MP1 = M + 1;
      for(int I=MP1; I<=N; I+6) {// L50
        SASUM = SASUM + fabs(SX[I]) + fabs(SX[I + 1]) + ABS[SX[I + 2]]
            + fabs(SX[I + 3]) + fabs(SX[I + 4]) + ABS[SX[I + 5]];
L50:  } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SAXPY(N,SA,SX,INCX,SY,INCY) {

      double SX[*],SY[*],SA;
       if(N <= 0 || SA == 0.E0) return;
       if(INCX == INCY)  if(INCX-1) 5,20,60;
L5:   } //continue;

      IX = 1;
      IY = 1;
       if(INCX < 0)IX = (-N+1)*INCX + 1 ;
       if(INCY < 0)IY = (-N+1)*INCY + 1 ;
      for(int I=1; I<=N; I++) {// L10
        SY[IY] = SY[IY] + SA*SX[IX];
        IX = IX + INCX;
        IY = IY + INCY;
L10:  } //continue;
      return;

L20:  M = ((N)%(4));
       if( M  ==  0 ) goto L40 ;
      for(int I=1; I<=M; I++) {// L30
        SY[I] = SY[I] + SA*SX[I];
L30:  } //continue;
       if( N  <  4 ) return;
L40:  MP1 = M + 1;
      for(int I=MP1; I<=N; I+4) {///L50
        SY[I] = SY[I] + SA*SX[I];
        SY[I + 1] = SY[I + 1] + SA*SX[I + 1];
        SY[I + 2] = SY[I + 2] + SA*SX[I + 2];
        SY[I + 3] = SY[I + 3] + SA*SX[I + 3];
L50:  } //continue;
      return;

L60:  } //continue;
      NS = N*INCX;
          for(int I=1; I<=NS; I+INCX) {//L70
          SY[I] = SA*SX[I] + SY[I];
L70:      } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SCOPY(N,SX,INCX,SY,INCY) {

      double SX[*],SY[*];
       if(N <= 0)return;
       if(INCX == INCY)  if(INCX-1) 5,20,60;
L5:   } //continue;

      IX = 1;
      IY = 1;
       if(INCX < 0)IX = (-N+1)*INCX + 1 ;
       if(INCY < 0)IY = (-N+1)*INCY + 1 ;
      for( I=1; I<=N; I++) {/// L10
        SY[IY] = SX[IX];
        IX = IX + INCX;
        IY = IY + INCY;
L10:  } //continue;
      return;

L20:  M = ((N)%(7));
       if( M  ==  0 ) goto L40 ;
      for(int I=1; I<=M; I++) {// L30
        SY[I] = SX[I];
L30:  } //continue;
       if( N  <  7 ) return;
L40:  MP1 = M + 1;
      for(int I=MP1; I<=N; I+7) {//L50
        SY[I] = SX[I];
        SY[I + 1] = SX[I + 1] ;
        SY[I + 2] = SX[I + 2] ;
        SY[I + 3] = SX[I + 3] ;
        SY[I + 4] = SX[I + 4] ;
        SY[I + 5] = SX[I + 5] ;
        SY[I + 6] = SX[I + 6] ;
L50:  } //continue;
      return;

L60:  } //continue;
      NS = N*INCX;
          for(int I=1; I<=NS; I+INCX) {// L70
          SY[I] = SX[I];
L70:      } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function SDOT[N][SX][INCX][SY][INCY] {

      double SX[*],SY[*];
      SDOT = 0.0E0;
       if(N <= 0)return;
       if(INCX == INCY) IF[INCX-1]5,20,60;
L5:   } //continue;

      IX = 1;
      IY = 1;
       if(INCX < 0)IX = (-N+1)*INCX + 1 ;
       if(INCY < 0)IY = (-N+1)*INCY + 1 ;
      for(int I=1; I<=N; I++) {//L10
        SDOT = SDOT + SX[IX]*SY[IY];
        IX = IX + INCX;
        IY = IY + INCY;
L10:  } //continue;
      return;

L20:  M = ((N)%(5));
       if( M  ==  0 ) goto L40 ;
      for(int I=1; I<=M; I++) {//L30
        SDOT = SDOT + SX[I]*SY[I];
L30:  } //continue;
       if( N  <  5 ) return;
L40:  MP1 = M + 1;
      for(int I=MP1; I<=N; I+5) {  //L50
        SDOT = SDOT + SX[I]*SY[I] + SX[I + 1]*SY[I + 1] +
             SX[I + 2]*SY[I + 2] + SX[I + 3]*SY[I + 3] + SX[I + 4]*SY[I + 4];
L50:  } //continue;
      return;

L60:  } //continue;
      NS=N*INCX;
      for(int I=1; I<=NS; I+INCX) {// L70
        SDOT = SDOT + SX[I]*SY[I];
L70:    } //continue;
      return;
      }
//c*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double SNRM2(int N, double *SX, int INCX)           // POLNYJ MARAZM!!!
{
      int       NEXT, NN, I, J=0;
      double    HITEST, SUM, XMAX = 0;
//const double   ZERO = 0;
//const double   ONE= 1;
double SNRM2=0;

const double  CUTLO = 4.441E-16;
const double  CUTHI = 1.304E19;

if(N  <=  0) return 0;


      NEXT = 30;
      SUM = 0;
      NN = N * INCX ;
      I = 1;

do {

       switch(NEXT)
               {
//               case 30  : goto L30;
               case 50  : goto L50;
               case 70  : goto L70;
               case 110 : goto L110;
               }

//        L30:

      if( fabs(SX[I])  >  CUTLO) goto L85;

      NEXT = 50;
      XMAX = 0;


        L50:

       if( SX[I]  ==  0)          goto L200;
       if( fabs(SX[I])  >  CUTLO) goto L85;

        NEXT = 70;

      goto L105;

        L100:
        I = J;
        NEXT = 110;

      SUM /=  pow2(SX[I]);

      L105: XMAX = fabs(SX[I]);
      SUM  += pow2((SX[I]/XMAX));
      goto L200;

        L70:

        if( fabs(SX[I])  <= CUTLO ) // goto L75;
                {

                L110:  if( fabs(SX[I])  <=  XMAX ) {
                                  SUM  += pow2((SX[I]/XMAX));
                                }
                else    {
                        SUM = 1. + SUM *pow2(XMAX/SX[I]);
                        XMAX = fabs(SX[I]);
                        }
                }
        else    {
                                        //                L75:
                       SUM *= pow2(XMAX);

                L85:  HITEST = CUTHI/double(N);

                for(J=I; J<=NN; J+=INCX)                                 // bug found by Toshi Sumikama 09/28/2021
                        {
                        if(fabs(SX[J])  >=  HITEST) goto L100;
                        SUM += pow2(SX[J]);
                        }

                SNRM2 = sqrt( SUM );
                return SNRM2;
                }

      L200:
      I += INCX;
      }
      while ( I  <=  NN );

//----------------------------------------
SNRM2 = XMAX * sqrt(SUM);

return SNRM2;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*void SSCAL(N,SA,SX,INCX) {

      double SA,SX[*] ;
       if(N <= 0)return;
       if(INCX == 1)goto L20;

      NS = N*INCX;
          for(int I=1; I<=NS; I+INCX) {// L10
          SX[I] = SA*SX[I];
L10:      } //continue;
      return;

L20:  M = ((N)%(5));
       if( M  ==  0 ) goto L40 ;
      for(int I=1; I<=M; I++) {//L30
        SX[I] = SA*SX[I];
L30:  } //continue;
       if( N  <  5 ) return;
L40:  MP1 = M + 1;
      for(int I=MP1; I<=N; I+5) { // L50
        SX[I] = SA*SX[I];
        SX[I + 1] = SA*SX[I + 1];
        SX[I + 2] = SA*SX[I + 2];
        SX[I + 3] = SA*SX[I + 3];
        SX[I + 4] = SA*SX[I + 4];
L50:  } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SSWAP(N,SX,INCX,SY,INCY) {

      double SX[*],SY[*],STEMP1,STEMP2,STEMP3;
       if(N <= 0)return;
       if(INCX == INCY)  if(INCX-1) 5,20,60;
L5:   } //continue;

      IX = 1;
      IY = 1;
       if(INCX < 0)IX = (-N+1)*INCX + 1 ;
       if(INCY < 0)IY = (-N+1)*INCY + 1 ;
      for(int I=1; I<=N; I++) {// L10
        STEMP1 = SX[IX];
        SX[IX] = SY[IY];
        SY[IY] = STEMP1;
        IX = IX + INCX;
        IY = IY + INCY;
L10:  } //continue;
      return;

L20:  M = ((N)%(3));
       if( M  ==  0 ) goto L40 ;
      for(int I=1; I<=M; I++) {//L30
        STEMP1 = SX[I];
        SX[I] = SY[I];
        SY[I] = STEMP1;
L30:  } //continue;
       if( N  <  3 ) return;
L40:  MP1 = M + 1;
      for(int I=MP1; I<=N; I+3) {// L50
        STEMP1 = SX[I];
        STEMP2 = SX[I+1];
        STEMP3 = SX[I+2];
        SX[I] = SY[I];
        SX[I+1] = SY[I+1];
        SX[I+2] = SY[I+2];
        SY[I] = STEMP1;
        SY[I+1] = STEMP2;
        SY[I+2] = STEMP3;
L50:  } //continue;
      return;
L60:  } //continue;

      NS = N*INCX;
        for(int I=1; I<=NS; I+INCX) { //L70
        STEMP1 = SX[I];
        SX[I] = SY[I];
        SY[I] = STEMP1;
L70:    } //continue;
      return;
      }
//*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function UNI[] {
      PARAMETER[
              CSAVE=362436./16777216.  ][
              CD=7654321./16777216.][
              CM=16777213./16777216.  ];
//                            2**24=16777216
      double U[17],S,T,USTART,C,UNIB;
      int I,J,II,JJ,K,KK,I1,J1,K1,L1,M1,ISEED;
// 
      SAVE U,I,J,K,C;

      double  U/
          0.8668672834288,  0.3697986366357,  0.8008968294805,
          0.4173889774680,  0.8254561579836,  0.9640965269077,
          0.4508667414265,  0.6451309529668,  0.1645456024730,
          0.2787901807898,  0.06761531340295, 0.9663226330820,
          0.01963343943798, 0.02947398211399, 0.1636231515294,
          0.3976343250467,  0.2631008574685/;
      double  I,J,K,C/17,5,24,CSAVE/;
//
      UNI = U[I]-U[J];
       if(UNI < 0.0)UNI = UNI+1.0;
      U[I] = UNI;
      I = I-1;
       if(I == 0)I = 17;
      J = J-1;
       if(J == 0)J = 17;

      C = C-CD;
       if(C < 0.0) C=C+CM;

      UNI = UNI-C;
       if(UNI < 0.0)UNI = UNI+1.0;
      return;

      ENTRY USTART[ISEED];

        I1 = ((fabs(ISEED))%(177))+1;
        J1 = ((fabs(ISEED))%(167))+1;
        K1 = ((fabs(ISEED))%(157))+1;
        L1 = ((fabs(ISEED))%(147))+1;
//
        for(int II=1; II<=17; II++) {//L2
          S = 0.0;
          T = 0.5;

          for( JJ=1; JJ<=K; JJ++) { //L3
                  M1 = ((mod(I1*J1,179)*K1)%(179));
                  I1 = J1;
                  J1 = K1;
                  K1 = M1;
                  L1 = ((53*L1+1)%(169));
                   if(((L1*M1)%(64)) >= 32)S=S+T;
L3:               T = .5*T;
L2:     U[II] = S;
        USTART = FLOAT[ISEED];
        return;

      ENTRY UNIB[KK];
         if(KK <= 24){
             K=24;
        }
 else {
             K=KK;
        }
        UNIB=FLOAT[K];
      }
//c*********************************************************************
*/
