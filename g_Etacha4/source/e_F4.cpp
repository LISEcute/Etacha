
#include "../win/e_myextern.h"
#include "../win/e_Constant.h"


extern void CalcOlegSum(double *v);
//void F(int NEQ, double X, double *U,double *UP);
//void F(double X, double *U, double *UP);
void INI_F4();

int f_num(int NCO);
int f_numP(int NCO);
int f_numPP(int NCO);
int f_II(int num);
int f_JJ(int I, int num);
int f_KK(int I,int J,int num);
int f_IM(int num);
int f_IKL(int M,int num);
int f_IN(int num);
int f_IKM(int N, int num);

int nequ = 0;  /// Oleg

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//void F(int /*NEQ*/, double /*X*/, double *U, double *UP)
//void F(double X, double *U, double *UP)

void F(double X, double *V, double *VP)
{
//c    *********************************************************
//c    creation des equations differentielles (fonction F)
//c    *********************************************************

//      common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911), NUM3(3329);
//c--------------------------------------------------------------
//c    ICO1=100*I+10*J+K
//c    ICO2=100*n3+n12
//c    ICO3=100*n4+n123
//c
//c    I= nb elec 2p , J= nb elec 2s , K= nb elec 1s
//c    n12=nb elec 1+2 , n3=nb elec 3, n123=nb elec 1+2+3 , n4=nb elec 4
//c
//c    y1s= nb moyen e- 1s , y2s= ...etc
//c
//c    This version takes into account in "partially" correlated way
//c    M shell lectrons, and similarly N shell electrons :
//c
//c    - for n=1 et 2 , 63 correlated states of the type {I,J,K}
//c    (fractions of ions with I e- 2p AND J e- 2s AND K e- 1s)
//c
//c    - for n=3 , 3+7+11 = 21 independent states
//c    (fractions of ions with 0,1,2 (,..) e- in 3s,3p and 3d)
//c    - for n=1+2+3 , 11*19=209 correlated states of the type {n12,n3}
//c
//c    - for n=4 , 33 independent states
//c    - for n=1+2+3+4 , 29*33=957 correlated states of the type {n123,n4}
//c
//c    corresponding to 63+21+209+33+957=1283 states and equations...
//c
//c--------------------------------------------------------------
//     common/seceff/sec(34),cor(48),secs(48);
//      common/secKLM/C12(209),D12(209),RAD2(209),AKL2(209),PA2(209),
//           AKM3(209),RAD3(209),ALM3(209),C3(209),D3(209),E3(209),DE3(209),
//           PA3(209),PA23(209);
//      common/sec1234/P14(957),C13(957),D13(957),C4M(957),D4M(957),
//                E4M(957),DE4M(957),PA4(957),AKLM4(957),PA123(957);
//      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
//      common/aug/AKLL,AKLM,ALMM,AM4;
//      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
//      common/proj/zp;

//c--------------------------------------------------------------
//c    Cross sections ...
//c--------------------------------------------------------------
//UP[0]=0;
//U [0]=0;

double *U  = V-1;
double *UP = VP-1;

//double *U  = V;
//double *UP = VP;

      nequ = nequ+1;
//not used      dum=x;
//not used      int Ndum = NEQ;
      int nzp  = Zp;
//------------------------------------------
      double C1s=GSEC[1];       // capture
      double C2s=GSEC[2];
      double C2p=GSEC[3];
      double C3s=GSEC[4];
      double C3p=GSEC[5];
      double C3d=GSEC[6];
      double C4 =GSEC[7];

      double D1s=GSEC[8];       // ionization
      double D2s=GSEC[9];
      double D2p=GSEC[10];
      double D3s=GSEC[11];
      double D3p=GSEC[12];
      double D3d=GSEC[13];
      double D4 =GSEC[14];

      double e1s2s =GSEC[15];
      double e1s2p =GSEC[16];

      double e1s3s =GSEC[17];  // e3s
      double e1s3p =GSEC[18];
      double e1s3d =GSEC[19];
      double e1s4  =GSEC[20];
      double e2s3s =GSEC[21];
      double e2s3p =GSEC[22];
      double e2s3d =GSEC[23];  // s3d
      double e2s4  =GSEC[24];
      double e2p3s =GSEC[25];
      double e2p3p =GSEC[26];   // p3p
      double e2p3d =GSEC[27];

      double e2p4  =GSEC[28];
      double e3s4  =GSEC[29];
      double e3p4  =GSEC[30];
      double e3d4  =GSEC[31];

      double e2s2p = GSEC[32];
      double e3s3p = GSEC[33];
      double e3p3d = GSEC[34];

//------------------------------------------
//c--------------------------------------------------------------
//c    actual cross sections for 1s,2s,2p states
//c    excitation included in rad3p1, rad3p2, rad 3s & rad 3d
//c    exchanges with n=4 included (wheighted with yN or 32-yN)
//c--------------------------------------------------------------
                                                // from Etacha4.cpp
                                                //      Rad3s   += GSEC[25];    e2p3s
                                                //      Rad3p1  += GSEC[18];    e1s3p
                                                //      Rad3p2  += GSEC[22];    e2s3p
                                                //      Rad3d   += GSEC[27];    e2p3d

// capture

      double C1se = C1s + y3s*e1s3s + y3p*Rad3p1 + y3d*e1s3d             +yN*e1s4;
      double C2se = C2s + y3s*e2s3s + y3p*Rad3p2 + y3d*e2s3d + ALMM*yMm  +yN*e2s4;
      double C2pe = C2p + y3s*Rad3s + y3p*e2p3p +  y3d*Rad3d + ALMM*yMm  +yN*e2p4;

// ionization

      double D1se = D1s + (2.-y3s)*e1s3s + (6.-y3p)*e1s3p + (10.-y3d)*e1s3d  + (32.-yN)*e1s4;
      double D2se = D2s + (2.-y3s)*e2s3s + (6.-y3p)*e2s3p + (10.-y3d)*e2s3d  + (32.-yN)*e2s4;
      double D2pe = D2p + (2.-y3s)*e2p3s + (6.-y3p)*e2p3p + (10.-y3d)*e2p3d  + (32.-yN)*e2p4;
//c--------------------------------------------------------------
//c--------------------------------------------------------------
//c    Loop for calculation of differential equations
//c    for the 63 correlated states {1s,2s,2p}
//c--------------------------------------------------------------
int M;

for(int  N=1; N<=63; N++)  //L200
      {
      int I=f_II(N);
      int J=f_JJ(I,N);
      int K=f_KK(I,J,N);
      int L=100*I+10*J+K;
      int nelec=I+J+K;
        //c------------------------------------------------------------------
        //c electron loss processes (exchanges with n=3 & 4 included)
        //c------------------------------------------------------------------
      UP[N] = -U[N]*(I*D2pe+J*D2se+K*D1se)
              -U[N]*((2.-K)*I*Rad+2*(K+J-K*J)*e1s2s)
              -U[N]*(2.*(3.*K+I-K*I)*e1s2p)            // bug found by Toshi Sumikama 09/28/2021
              -U[N]*(2.*(3.*J+I-I*J)*e2s2p);

      if (nelec < nzp)   UP[N] -= U[N]*((6.-I)*C2pe+(2.-J)*C2se+(2.-K)*C1se);

        //c.......... KLL & KLM Auger ....................................
        //c    (electron loss)
        //c.................................................................
       if (K < 2)
                {
                if (I >= 2) UP[N] -= U[N]*AKLL*I*(I-1)*(2-K);
                if (I >= 1) {
                                          UP[N] -=  U[N]*AKLM*I*yM*(2-K);
                             if (J >= 1)  UP[N] -=  U[N]*AKLL*2.*I*J*(2-K);
                             }
                 if (J >= 1)    UP[N] -= U[N]*AKLM*J*yM*(2-K);
                 if (J >= 2)    UP[N] -= U[N]*AKLL*J*(J-1)*(2-K);
                }
        //c.................................................................
        //c    electron gain processses
        //c.................................................................
        if (K != 0)
                {
                        //c........ Augers KLL ........................
                if (I <= 4)
                        {
                        M=f_num(L+199);
                        UP[N] += (3-K)*(I+2)*(I+1)*AKLL*U[M];
                        }
                if ((I <= 5) && (J <= 1))
                        {
                        M=f_num(L+109);
                        UP[N] += (3-K)*2*(I+1)*(J+1)*AKLL*U[M];
                        }
                if (J == 0)
                        {
                        M=f_num(L+19);
                        UP[N] += (3-K)*2.*AKLL*U[M];
                        }
                        //c........ Augers KLM ........................
                if (I <= 5)
                        {
                        M=f_num(L+99);
                        UP[N] += (3-K)*(I+1)*yM*AKLM*U[M];
                        }
                 if (J <= 1)
                        {
                        M=f_num(L+9);
                        UP[N] += (3-K)*(J+1)*yM*AKLM*U[M];
                        }
                }
                //c--------------------------------------------------------------
        //c States are populated by electron gain if the considered states possess
        //c electrons,
        //c and intershell processes if states possess holes
        //c--------------------------------------------------------------
        //c    test on unpopulated states if nelec > zp
     if(nelec <= nzp)
                {

                        //c..............................................................
                        //c              2p state
                        //c..............................................................
                if (I != 0) {
                      M=f_num(L-100);
                      UP[N] += (7-I)*C2pe*U[M];
                      if (K != 2)
                                {
                                M=f_num(L-99);
                                UP[N] += (7-I)*(K+1)*e1s2p*U[M];
                                }
                      if (J != 2)
                                {
                                M=f_num(L-90);
                                UP[N] += (7-I)*(J+1)*e2s2p*U[M];
                                }
                }
        //c..............................................................
        //c              2s state
        //c..............................................................
       if (J != 0) {
                M=f_num(L-10);
                UP[N] += (3-J)*C2se*U[M];
                if (K != 2) {
                      M=f_num(L-9);
                      UP[N] += (3-J)*(K+1)*e1s2s*U[M];
                      }
               if (I != 6) {
                      M=f_num(L+90);
                      UP[N] += (3-J)*(I+1)*e2s2p*U[M];
                      }
               }
        //c..............................................................
        //c              1s state
        //c..............................................................
       if (K != 0) {
                M=f_num(L-1);
                UP[N] += (3-K)*C1se*U[M];
                if (I != 6) {
                       M=f_num(L+99);
                       UP[N] += (3-K)*(I+1)*(Rad+e1s2p)*U[M];
                       }
                if (J != 2) {
                       M=f_num(L+9);
                       UP[N] += (3-K)*(J+1)*e1s2s*U[M];
                       }
                }
                //c--------------------------------------------------------------
                //c    end of test on unpopulated states if nelec > zp
             }
                //c--------------------------------------------------------------
                //c Sates are populated by ionization if if the considered states possess
                //c holes
                //c--------------------------------------------------------------
                //c

                //c              2p state
       if (I != 6) {
              M=f_num(L+100);
              UP[N] += (1+I)*D2pe*U[M];
              }
                //c              2s state
       if (J != 2) {
              M=f_num(L+10);
              UP[N] += (1+J)*D2se*U[M];
              }
                //c              1s state
       if (K != 2) {
              M=f_num(L+1);
              UP[N] += (K+1)*D1se*U[M];
              }
   } //continue;L200:



//c--------------------------------------------------------------
//c    actual cross sections for 3s,3p,3d states
//c    exchanges with n=4 included (wheighted with yN or 32-yN)
//c--------------------------------------------------------------
double  SLS = D3s + (6.-y2p)*Rad3s + (2.-y2s)*e2s3s  + (2.-y1s)*e1s3s  + (6.-y3p)*e3s3p                 + (2.-y1s)*yL*AKLM + (32.-yN)*e3s4;
double  SLP = D3p + (6.-y2p)*e2p3p + (2.-y2s)*Rad3p2 + (2.-y1s)*Rad3p1 + (2.-y3s)*e3s3p+(10.-y3d)*e3p3d + (2.-y1s)*yL*AKLM + (32.-yN)*e3p4;
double  SLD = D3d + (6.-y2p)*Rad3d + (2.-y2s)*e2s3d  + (2.-y1s)*e1s3d  + (6.-y3p)*e3p3d                 + (2.-y1s)*yL*AKLM + (32.-yN)*e3d4;

double  SGS  = C3s + y2p*e2p3s + y2s*e2s3s + y1s*e1s3s + y3p*e3s3p +             yN*e3s4 + AM4*yNn + Rad4*yN;
double  SGP  = C3p + y2p*e2p3p + y2s*e2s3p + y1s*e1s3p + y3s*e3s3p + y3d*e3p3d + yN*e3p4 + AM4*yNn + Rad4*yN;
double  SGD  = C3d + y2p*e2p3d + y2s*e2s3d + y1s*e1s3d + y3p*e3p3d +             yN*e3d4 + AM4*yNn + Rad4*yN;
//c--------------------------------------------------------------
//c    differential equations for 3s,3p et 3d populations :
//c    (no loop)
//c                64 to 66 -> 3s states with 0,1 or 2 e-
//c                67 to 73 -> 3p states with 0,1 ... 6 e-
//c                74 to 84 -> 3d states with 0,1 ... 10 e-
//c--------------------------------------------------------------
      UP[64] = -U[64]*2.*SGS + U[65]*(SLS);
      UP[65] = +U[64]*2.*SGS - U[65]*(SGS+SLS) +U[66]*2.*SLS;
      UP[66] =                +U[65]*(SGS)     -U[66]*2.*SLS;

      UP[67] = -U[67]*(6.*SGP       )              +U[68]*1.*SLP;
      UP[68] = -U[68]*(5.*SGP+1.*SLP)+U[67]*6.*SGP +U[69]*2.*SLP;
      UP[69] = -U[69]*(4.*SGP+2.*SLP)+U[68]*5.*SGP +U[70]*3.*SLP;
      UP[70] = -U[70]*(3.*SGP+3.*SLP)+U[69]*4.*SGP +U[71]*4.*SLP;
      UP[71] = -U[71]*(2.*SGP+4.*SLP)+U[70]*3.*SGP +U[72]*5.*SLP;
      UP[72] = -U[72]*(1.*SGP+5.*SLP)+U[71]*2.*SGP +U[73]*6.*SLP;
      UP[73] = -U[73]*(6.*SLP       )+U[72]*1.*SGP              ;
//c
      UP[74] = -U[74]*(10.*SGD       )                 + U[75]*1.*SLD;
      UP[75] = -U[75]*( 9.*SGD+1.*SLD) + U[74]*10.*SGD + U[76]*2.*SLD;
      UP[76] = -U[76]*( 8.*SGD+2.*SLD) + U[75]*9.*SGD  + U[77]*3.*SLD;
      UP[77] = -U[77]*( 7.*SGD+3.*SLD) + U[76]*8.*SGD  + U[78]*4.*SLD;
      UP[78] = -U[78]*( 6.*SGD+4.*SLD) + U[77]*7.*SGD  + U[79]*5.*SLD;
      UP[79] = -U[79]*( 5.*SGD+5.*SLD) + U[78]*6.*SGD  + U[80]*6.*SLD;
      UP[80] = -U[80]*( 4.*SGD+6.*SLD) + U[79]*5.*SGD  + U[81]*7.*SLD;
      UP[81] = -U[81]*( 3.*SGD+7.*SLD) + U[80]*4.*SGD  + U[82]*8.*SLD;
      UP[82] = -U[82]*( 2.*SGD+8.*SLD) + U[81]*3.*SGD  + U[83]*9.*SLD;
      UP[83] = -U[83]*( 1.*SGD+9.*SLD) + U[82]*2.*SGD  + U[84]*10.*SLD;
      UP[84] = -U[84]*(       10.*SLD) + U[83]*1.*SGD                 ;
//c
//c........... AUGERS LMM .......................................
      double AL =  ALMM*(8.-yL);
      double ysp = 2.*(y3s+y3p);
      double ysd = 2.*(y3s+y3d);
      double ypd = 2.*(y3p+y3d);
        //c.................................................

      UP[64] += AL*(+U[65]*1.*ypd  + U[66]*2.         );
      UP[65] += AL*(-U[65]*1.*ypd  + U[66]*2.* ypd    );
      UP[66] += AL*(               - U[66]*2.*(ypd+1.));

      UP[67] += AL*(+U[68]*1.*ysd                     + U[69]*2. );
      UP[68] += AL*(-U[68]*1.*(ysd+0.) + U[69]*2.*ysd + U[70]*6. );
      UP[69] += AL*(-U[69]*2.*(ysd+1.) + U[70]*3.*ysd + U[71]*12.);
      UP[70] += AL*(-U[70]*3.*(ysd+2.) + U[71]*4.*ysd + U[72]*20.);
      UP[71] += AL*(-U[71]*4.*(ysd+3.) + U[72]*5.*ysd + U[73]*30.);
      UP[72] += AL*(-U[72]*5.*(ysd+4.) + U[73]*6.*ysd);
      UP[73] += AL*(-U[73]*6.*(ysd+5.));

      UP[74] += AL*(+U[75]*1.* ysp                     + U[76]*2. );
      UP[75] += AL*(-U[75]*1.* (ysp+0.) + U[76]*2.*ysp + U[77]*6. );
      UP[76] += AL*(-U[76]*2.* (ysp+1.) + U[77]*3.*ysp + U[78]*12.);
      UP[77] += AL*(-U[77]*3.* (ysp+2.) + U[78]*4.*ysp + U[79]*20.);
      UP[78] += AL*(-U[78]*4.* (ysp+3.) + U[79]*5.*ysp + U[80]*30.);
      UP[79] += AL*(-U[79]*5.* (ysp+4.) + U[80]*6.*ysp + U[81]*42.);
      UP[80] += AL*(-U[80]*6.* (ysp+5.) + U[81]*7.*ysp + U[82]*56.);
      UP[81] += AL*(-U[81]*7.* (ysp+6.) + U[82]*8.*ysp + U[83]*72.);
      UP[82] += AL*(-U[82]*8.* (ysp+7.) + U[83]*9.*ysp + U[84]*90.);
      UP[83] += AL*(-U[83]*9.* (ysp+8.) + U[84]*10.*ysp);
      UP[84] += AL*(-U[84]*10.*(ysp+9.));
//c    **********************************************************
//c                  correlated states
//c--------------------------------------------------------------
//c    Loop for calculation of differential equations
//c    for the 209 correlated states {(1+2)*N,3l*M}
//c    Actual cross sections as a function of state code are
//c    calculated in secmean.for
//c--------------------------------------------------------------
if(EtachaVersion == etacha_v23)  return;

      for(int  N=85; N<=293; N++) //L300
           {
           int M=f_IM(N);
           int K=f_IKL(M,N);
           int L=100*M+K;
           int nelec=M+K;
           int NN=N-84;
                //c------------------------------------------------------------------
                //c    electron loss processes (exchanges with n=3 & 4 included)
                //c    E3=2->3,DE3=3->2
                //c    !! RAD3 is included in DE3, and RAD2 do not change populations !!
                //c------------------------------------------------------------------

           UP[N] = -U[N]*(E3[NN]+DE3[NN]+AKM3[NN]+ALM3[NN]+D12[NN]+AKL2[NN]+D3[NN]);

           if(nelec < nzp)    UP[N] = UP[N]-U[N]*(C12[NN]+C3[NN]);

        //c.................................................................
        //c    electron gain processes
        //c.................................................................
       if(nelec <= nzp)
                {                //............................
                if (K != 0)
                        {
                        int IJ = f_numP(L-1);
                        int KL = IJ-84;
                        UP[N] += C12[KL]*U[IJ];
                        if (M != 18)
                                {
                                IJ= f_numP(L+99);
                                KL=IJ-84;
                                UP[N] += DE3[KL]*U[IJ];
                                }
                        }

               if (K != 10) {
                        int IJ = f_numP(L+1);
                        int KL = IJ-84;
                        UP[N] += D12[KL]*U[IJ];
                         if (M != 0) {
                                  IJ= f_numP(L-99);
                                  KL=IJ-84;
                                  UP[N] += E3[KL]*U[IJ];
                                  }
                        }

               if (M != 0) {
                        int IJ = f_numP(L-100);
                        int KL = IJ-84;
                        UP[N] += C3[KL]*U[IJ];
                        }

               if (M != 18) {
                        int IJ = f_numP(L+100);
                        int KL = IJ-84;
                        UP[N] += D3[KL]*U[IJ];
                        }

                //c.......... Augers   ....................................
               if (K <= 8) {
                        int IJ= f_numP(L+1);
                        int KL=IJ-84;
                        UP[N] += AKL2[KL]*U[IJ];
                        }

               if((K != 0) && (M <= 17)) {
                         int IJ= f_numP(L+100);
                         int KL=IJ-84;
                         UP[N] += AKM3[KL]*U[IJ];
                         }

               if((K != 0) && (M <= 16)) {
                         int IJ= f_numP(L+199);
                         int KL=IJ-84;
                         UP[N] += ALM3[KL]*U[IJ];
                         }
              }

 } //continue;L300:
//c--------------------------------------------------------------
//c    actual cross sections for n=4 states
//c--------------------------------------------------------------
if(EtachaVersion <= etacha_v3)  return;


      double SG4 = C4 + y1s*e1s4 + y2s*e2s4 + y2p*e2p4 + y3s*e3s4 + y3p*e3p4 + y3d*e3d4;
      double SL4 = D4 +(2.-y1s)*e1s4 + (2.-y2s)*e2s4 +  (6.-y2p)*e2p4 +
                       (2.-y3s)*e3s4 +( 6.-y3p)*e3p4 + (10.-y3d)*e3d4 +
                       ((2.-y3s)+(6.-y3p)+(10.-y3d))*Rad4;

      double AM = AM4*(18.-yM);
//c--------------------------------------------------------------
//c    differential equations for n=4 population :
//c                294 to 326 -> 33 n=4 states with 0,1 .....or 32 e-
//c    (corrected 12/12)
//c--------------------------------------------------------------
 for(int i=1; i<=33; i++)  //L400
        {
        double ri=i;
        double UP1=0.;
        double UP2=0.;
        double UP3=0.;
        double UP4=0.;
        double UP5=0.;
        double UP6=0.;

               if (i != 1 )  UP1 = (34.-ri)*SG4       *U[i+292];
               if (i != 33)  UP2 = (ri)*SL4           *U[i+294];

               if (i != 33)  UP3 = -(33.-ri)*SG4      *U[i+293];
               if (i != 1 )  UP4 = -(ri -1.)*SL4      *U[i+293];

               if (i <  32)  UP5 =  AM* ri    *(ri+1.)*U[i+295];
               if (i >  2 )  UP6 = -AM*(ri-1.)*(ri-2.)*U[i+293];

        UP[i+293]=UP1+UP2+UP3+UP4+UP5+UP6;
        }
//c--------------------------------------------------------------
//c    Loop for calculation of differential equations
//c    for the 957 correlated {(1+2+3)*N1,4l*N4} states
//c--------------------------------------------------------------
if(EtachaVersion <= etacha_v34)  return;

 for(int  NN=1; NN<=957; NN++)  //L500
      {
      int N=NN+326;
      int NEN=f_IN(N);
      int K=f_IKM(NEN,N);
      int L=100*NEN+K;
      int nelec=NEN+K;
        //c------------------------------------------------------------------
        //c electron loss processes (exchanges with 4 included)
        //c
        //c
        //c------------------------------------------------------------------

      UP[N] = -U[N]*(D13[NN]+D4M[NN]+E4M[NN]+DE4M[NN]+AKLM4[NN]);
      if(nelec < nzp)    UP[N] -= U[N]*(C13[NN]+C4M[NN]);

        //c.................................................................
        //c    electron gain processes
        //c.................................................................
       if(nelec <= nzp)
            {
                                //............................
            if (K != 0) {
                        //     123 capture and 4->123 desexcitation
                        int IJ = f_numPP(L-1);
                        int KL = IJ-326;
                        UP[N] += C13[KL]*U[IJ];
                        if (NEN != 32) {
                              int IJ = f_numPP(L+99);
                              int KL = IJ-326;
                              UP[N] += DE4M[KL]*U[IJ];
                              }
                        }

           if (K != 28) {
                        //     123 ionization and 123->4 excitation
                       int IJ = f_numPP(L+1);
                       int KL=IJ-326;
                       UP[N] += D13[KL]*U[IJ];
                       if (NEN != 0) {
                              int IJ = f_numPP(L-99);
                              int KL=IJ-326;
                              UP[N] += E4M[KL]*U[IJ];
                              }
                      }

           if (NEN != 0) {
                                //     n=4 capture
                      int IJ = f_numPP(L-100);
                      int KL=IJ-326;
                      UP[N] += C4M[KL]*U[IJ];
                      }

            if (NEN != 32) {
                        //     n=4 ionization
                      int IJ = f_numPP(L+100);
                      int KL = IJ-326;
                      UP[N] += D4M[KL]*U[IJ];
                      }

           }
        //c.......... (KLM)NN Augers ....................................
       if((K != 0) && (NEN <= 30))
                {
                int IJ = f_numPP(L+199);
                int KL = IJ-326;
                UP[N] += AKLM4[KL]*U[IJ];
                }
   } //continue; L500:
//c******************************************************************
CalcOlegSum(U);
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void INI_F4() {
//c...............................................................
//     Construction of code table ICO(n)
//c...............................................................
      for(/*L21*/int m=1; m<=733;  m++) NUM1[m]=0;
      for(/*L22*/int m=1; m<=1911; m++) NUM2[m]=0;

      int n1=0;
      for(/*L31*/int i=0; i<=6; i++) {
      for(/*L32*/int j=0; j<=2; j++) {
      for(/*L33*/int k=0; k<=2; k++) {
              n1++;
              ICO1[n1]=i*100+j*10+k;
              int ICOD1=ICO1[n1]+1;
              NUM1[ICOD1]=n1;
        }}}

      int n2=84;
      for(/*L41*/ int ii=0; ii<=18; ii++) {
      for(/*L42*/ int jj=0; jj<=10; jj++) {
              n2++;
              ICO2[n2]=ii*100+jj;
              int ICOD2=ICO2[n2]+1;
              NUM2[ICOD2]=n2;
        }}


      int n3=326;
      for(/*L51*/ int ii=0; ii<=32; ii++) {
      for(/*L52*/ int jj=0; jj<=28; jj++) {
              n3++;
              ICO3[n3]=ii*100+jj;
              int ICOD3=ICO3[n3]+1;
              NUM3[ICOD3]=n3;
        }}
return;
}

//c****************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int f_num(int NCO)
{
int num=NUM1[NCO+1];
return num;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int f_numP(int NCO)
{
int numP=NUM2[NCO+1];
return numP;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int f_II(int num)
{
int f_II=ICO1[num]/100;
return f_II;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int f_JJ(int I, int num)
{
int f_JJ=(ICO1[num]-I*100)/10;
return f_JJ;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int f_KK(int I,int J,int num)
{
int KK=(ICO1[num]-I*100-J*10);
return KK;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int f_IM(int num)
{
int IM=ICO2[num]/100;
return IM;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int f_IKL(int M,int num)
{
int IKL=ICO2[num]-M*100;
return IKL;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int f_numPP(int NCO)
{
int numPP=NUM3[NCO+1];
return numPP;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int f_IN(int num)
{
int f_IN=ICO3[num]/100;
return f_IN;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int f_IKM(int N, int num)
{
int f_IKM=ICO3[num]-N*100;
return f_IKM;
}
