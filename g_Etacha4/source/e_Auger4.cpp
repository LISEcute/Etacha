#include "../win/e_myextern.h"
#include "../win/e_Constant.h"


int f_numPP(int NCO);

double Auger(double *Y,double Zp,double *PR);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double Auger(double *Y,double Zp,double *PR)
{
//c..........................................
//c    Calculates final charge state distribution for each target thickness
//c    of the output files;
//c    takes into account autoionisation at the exit of target
//c..........................................
double YA[1284],P1234[30][34];
                /*      double P1234[29][33];
                      common /etats/ICO1[63],ICO2[294],ICO3[1284],NUM1[733],NUM2[1911], NUM3[3329];
                      common/secKLM/C12[209],D12[209],RAD2[209],AKL2[209],PA2[209],AKM3[209],RAD3[209],ALM3[209],C3[209],D3[209],E3[209],DE3[209], PA3[209],PA23[209];
                      common/sec1234/P14[957],C13[957],D13[957],C4M[957],D4M[957],E4M[957],DE4M[957],PA4[957],AKLM4[957],PA123[957];
                      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
                      common/aug/AKLL,AKLM,ALMM,AM4;*/

      double YTOT=0.;
      for(int N=0;   N< 327;  N++)  { YA[N]=0;                  }
      for(int N=327; N<=1283; N++)  { YA[N]=Y[N]; YTOT += Y[N]; }

if(EtachaVersion < etacha_v4)  return YTOT;


                                //     Calculation of correlated states (n123,n4) probabilities after autoionisation
   for(int N1=0; N1<=29; N1++)
        for(int N2=0; N2<=33; N2++)
                         P1234[N1][N2]=0.;
//.....................

  for(int N1=1; N1<=29; N1++)
     for(int N4=33; N4>=1; N4--)
       {
       int I=N1-1;
       int M=N4-1;
       int INM=100*M+I;
       int N=f_numPP(INM);
                //     N stands for indices of correlated states (between 327 and 1283)
                //     effect of KLL,KLM and LMM (PA123) Auger effects on these states
                //     (gives rise to very little change in most cases)
       int L=N-326;

       if(I > 2)
              {
              P1234[I+1][M+1] += (1.-PA123[L])*YA[N];
              P1234[I  ][M+1] +=     PA123[L] *YA[N];
              }
         else P1234[I+1][M+1] += YA[N];
         }
//....................
//     effet of (KLM)NN (PA4)) Augers effect on P1234 fractions (as calculated above)
//     fractions with N-1 e- KLM and M-1 e- N initially
//     fractions with 28 e- KLM (N=29) or no N shell e- (M=1) are unchanged
  for(int N=1; N<=28; N++)
      for(int M=33; M>=2; M--)                 //     (states with at least one hole in K or L or M shells)
        {
        int INM=100*(M-1)+N-1;
        int NM=f_numPP(INM)-326;

        if (M > 2) {
                    P1234[N+1][M-1] += (1.-PA4[NM])*P1234[N][M];   //     radiative decay
                    P1234[N+1][M-2] +=     PA4[NM] *P1234[N][M];   //     Auger decay

                    //     NO Auger decay (as a test)
                    //     P1234(N+1,M-1)=P1234(N+1,M-1)+PA4(NM)*P1234(N,M)

                    P1234[N][M]=0.;
                    }
        else if (M == 2)                 //     only one N shell electron :
                        {
                        P1234[N+1][M-1] += P1234[N][M];
                        P1234[N][M]=0.;
                        }
        }

//..................
//     conbine configurations and calculates final probabilities

 for(int N=0; N<=62; N++)   PR[N]=0.;

 for(int N=1; N<=29; N++)
      for(int M=1; M<=33; M++)
          {
          int Ne1 = N+M-1;
          int Zp1 = int(Zp)+1;
          if  (Ne1 <= Zp1) PR[Ne1] += P1234[N][M];
          else             PR[Zp1] += P1234[N][M];
          }

  for(int N=1; N<=61; N++) if(PR[N] < 0.0000001) PR[N]=0;

return YTOT;
}

