#include "../win/e_myextern.h"
#include "../win/e_Constant.h"

extern int f_num(int NCO);
extern int f_numP(int NCO);
extern int f_numPP(int NCO);
extern int f_II(int num);
extern int f_JJ(int I, int num);
extern int f_KK(int I,int J,int num);
extern int f_IM(int num);
extern int f_IKL(int M,int num);
extern int f_IN(int num);
extern int f_IKM(int N, int num);



void SecMean(double *Y);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SecMean(double *Y)
{
      double PKLM[210],C12_loc[210],D12_loc[210],C3_loc[210],D3_loc[210],EKLM4[210],DEKLM4[210];
/*
C12[] -> C12_loc[]
C12e[] -> C12[]    in order to correspond to other common blocks



      common /etats/ICO1[63],ICO2[294],ICO3[1284],NUM1[733],NUM2[1911],
                NUM3[3329];
      common/seceff/sec[34],cor[48],secs[48];

      common/secKLM/C12[209],D12[209],RAD2[209],AKL2[209],PA2[209],
           AKM3[209],RAD3[209],ALM3[209],C3[209],D3[209],E3[209],DE3[209],
           PA3[209],PA23[209];

      common/sec1234/P14[957],C13[957],D13[957],C4M[957],D4M[957],
                E4M[957],DE4M[957],PA4[957],AKLM4[957],PA123[957];
      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
      common/aug/AKLL,AKLM,ALMM,AM4;
      */
//c    ***********************************************
//c    calculation of actual (mean) cross sections
//c    ***********************************************
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

//not used      double e1s2s =GSEC[15];
//not used      double e1s2p =GSEC[16];

      double e1s3s =GSEC[17];
      double e1s3p =GSEC[18];
      double e1s3d =GSEC[19];
      double e1s4  =GSEC[20];
      double e2s3s =GSEC[21];
      double e2s3p =GSEC[22];
      double e2s3d =GSEC[23];
      double e2s4  =GSEC[24];
      double e2p3s =GSEC[25];
      double e2p3p =GSEC[26];
      double e2p3d =GSEC[27];
      double e2p4  =GSEC[28];
      double e3s4  =GSEC[29];
      double e3p4  =GSEC[30];
      double e3d4  =GSEC[31];
//not used      double e2s2p =GSEC[32];
//not used      double e3s3p=GSEC[33];
//not used      double e3p3d=GSEC[34];
//     add exchanges with n=4 for (1+2)*(3) fractions
      double C1se = C1s + yN*e1s4;
      double C2se = C2s + yN*e2s4;
      double C2pe = C2p + yN*e2p4;
//     ynn is the mean value of n(n-1) in n=4 shell...
//     AM4 is about 3.*ALMM
      double C3se = C3s + yN*e3s4 + AM4*yNn + Rad4*yN;
      double C3pe = C3p + yN*e3p4 + AM4*yNn + Rad4*yN;
      double C3de = C3d + yN*e3d4 + AM4*yNn + Rad4*yN;

      double D1se = D1s + (32.-yN)*e1s4;
      double D2se = D2s + (32.-yN)*e2s4;
      double D2pe = D2p + (32.-yN)*e2p4;
      double D3se = D3s + (32.-yN)*e3s4;
      double D3pe = D3p + (32.-yN)*e3p4;
      double D3de = D3d + (32.-yN)*e3d4;

for(int N=0; N<=209; N++)
        {
        PKLM[N]=0.;
        C12_loc[N]=0.;            D12_loc[N]=0.;
        C12[N]=0.;                D12[N]=0.;
        RAD2[N]=0.;               AKL2[N]=0.;
        PA2[N]=0.;
        RAD3[N]=0.;
        AKM3[N]=0.;               ALM3[N]=0.;
        C3[N]=0.;                 D3[N]=0.;
        C3_loc[N]=0.;             D3_loc[N]=0.;
        EKLM4[N]=0.;              DEKLM4[N]=0.;
                        //     E3, DE3 = 1+2 <-> 3
        E3[N]=0.;                DE3[N]=0.;
                //     (KL)MM Autoionization :
        PA3[N]=0.;
                //     KLM Autoionization :
        PA23[N]=0.;
        }

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW cicle L4 start
//     capture and ionization cross section in n=3,
//     as well as intershell exchanges, radiatives and Auger transitions
//     ne12 = (electron nb in n=1+2)+1
//     ne3 = (electron nb in n=3)+1
//     NN = state code (ne12,ne3)-84, in between 1 and 209 (=11*19)
//     a loop is made on a given number of electrons in 1s,2s,2p subshells
//     (n12 code in between 1 and 63), and a given number of electrons
//     in 3s,3p and 3d (L,LL,LLL): this determines one of the possible configurations
//     of a (ne12,ne3) state whose probability is given by Y(N)

for(int N12=1; N12<=63; N12++)
 for(int L=64; L<=66; L++)           // 3s
  for(int LL=67; LL<=73; LL++)       // 3p
    for(int LLL=74; LLL<=84; LLL++)  // 3d
      {
      int I = f_II(N12);      //  2p
      int J = f_JJ(I,N12);    //  2s
      int K = f_KK(I,J,N12);  //  1s

      int ne12  = I+J+K        +1;
      int ne3   = L+LL+LLL-205 +1;
      int ICOD2 = 100*(ne3-1)+(ne12-1);

      int NN = f_numP(ICOD2)-84;

      double PT = Y[N12]*Y[L]*Y[LL]*Y[LLL];
      PKLM[NN] += PT;

        //c--------------------------------------------------------------
        //c    Cross sections for n=1+2
        //c--------------------------------------------------------------
      C12_loc[NN] += PT*((2-K)*C1s +(2-J)*C2s+ (6-I)*C2p);
      C12[NN]     += PT*((2-K)*C1se+(2-J)*C2se+(6-I)*C2pe);

      D12_loc[NN] += PT*(K*D1s  +J*D2s  +I*D2p);
      D12[NN]     += PT*(K*D1se +J*D2se +I*D2pe);

      if (K < 2) {
               RAD2[NN] += PT*(2-K)*I*Rad;
               if ((I+J) > 1) {
                        AKL2[NN] += PT* (2-K)*(I+J)*(I+J-1)*AKLL;
                        PA2[NN]  += PT*((2-K)*(I+J)*(I+J-1)*AKLL)/ ((2-K)*(I+J)*(I+J-1)*AKLL+(2-K)*I*Rad);
                        }
                }
        //c--------------------------------------------------------------
        //c    Cross sections for n=3
        //c--------------------------------------------------------------
      C3_loc[NN] += PT*((66-L)*C3s +(73-LL)*C3p   +(84-LLL)*C3d);
      C3[NN]     += PT*((66-L)*C3se+(73-LL)*C3pe  +(84-LLL)*C3de);
      D3_loc[NN] += PT*((L-64)*D3s+ (LL-67)*D3p   +(LLL-74)*D3d);
      D3[NN]     += PT*((L-64)*D3se+(LL-67)*D3pe  +(LLL-74)*D3de);
        //c--------------------------------------------------------------
        //c    Cross sections for n=1+2 <-> n=3
        //c--------------------------------------------------------------
      E3[NN] += PT*((66-L)*(K*e1s3s+J*e2s3s+I*e2p3s)+(73-LL)*(K*e1s3p+J*e2s3p+I*e2p3p)+(84-LLL)*(K*e1s3d+J*e2s3d+I*e2p3d));
      DE3[NN]+= PT*(  (L-64)*((2-K)*e1s3s+(2-J)*e2s3s+(6-I)*(e2p3s+Rad3s))+
                      (LL-67)*((2-K)*(e1s3p+Rad3p1)+(2-J)*(e2s3p+Rad3p2)+(6-I)*e2p3p)+
                      (LLL-74)*((2-K)*e1s3d+(2-J)*e2s3d+(6-I)*(e2p3d+Rad3d)));

      double RAD3NN=((L-64)*(6-I)*Rad3s +(LL-67)*((2-K)*Rad3p1+(2-J)*Rad3p2) +(LLL-74)*(6-I)*Rad3d);
      double AKM3NN=AKLM*(ne3-1)*(I+J)*(2-K);
                //     ne3 is the (electron number in n=3) +1
                //     I,J and K are the real number of electrons in n=2 and 1
       double ALM3NN =0. ;
       if (ne3 > 2)  ALM3NN =ALMM*(ne3-1)*(ne3-2)*(8-J-I);


      RAD3[NN] += PT*RAD3NN;
      AKM3[NN] += PT*AKM3NN;
      ALM3[NN] += PT*ALM3NN;

      double ATOT=ALM3NN+AKM3NN+RAD3NN;
       if(ATOT > 0.000000000001) {
               PA23[NN] += PT*(AKM3NN)/ATOT;
               PA3[NN]  += PT*(ALM3NN)/ATOT;
               }
        //c--------------------------------------------------------------
        //c    Cross sections for n=1+2+3 <-> n=4
        //c    unwheighted for n=4 pop (but dependant on 1+2+3 pop)
        //c    PT is the probability of (n12,n3s,n3p,n3d) configuration
        //c--------------------------------------------------------------
      EKLM4[NN]  += PT*(K*e1s4+J*e2s4+I*e2p4 +(L-64)*e3s4+(LL-67)*e3p4+(LLL-74)*e3d4);
      DEKLM4[NN] += PT*( (2-K)*e1s4+(2-J)*e2s4+(6-I)*e2p4
                         +(66-L)*e3s4+(73-LL)*e3p4+(84-LLL)*e3d4
                         +(66-L)*Rad4+(73-LL)*Rad4+(84-LLL)*Rad4);
       } //continue; L4:
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW cicle L4 stop


//c---------------------------------------------------------------
//c    Normalization of wheighted cross sections
//c---------------------------------------------------------------
 for(int  N=1; N<=209; N++)  // L5
      {
      int NUE = N+84;
      int NM  = f_IM(NUE);
      int NKL = f_IKL(NM,NUE);
      double PKLMN=PKLM[N];

      if(PKLMN >= 0.000000000001) {
              C12_loc[N] /= PKLMN;
              D12_loc[N] /= PKLMN;
              C12[N]     /= PKLMN;
              D12[N]     /= PKLMN;
              RAD2[N]    /= PKLMN;
              AKL2[N]    /= PKLMN;
              PA2[N]     /= PKLMN;
              AKM3[N]    /= PKLMN;
              RAD3[N]    /= PKLMN;
              ALM3[N]    /= PKLMN;
              C3_loc[N]  /= PKLMN;
              D3_loc[N]  /= PKLMN;
              C3[N]      /= PKLMN;
              D3[N]      /= PKLMN;
              E3[N]      /= PKLMN;
              DE3[N]     /= PKLMN;
              PA3[N]     /= PKLMN;
              PA23[N]    /= PKLMN;
              EKLM4[N]   /= PKLMN;
              DEKLM4[N]  /= PKLMN;
              }
         else {
                //c    small probability configurations
              C12_loc[N]=(C1s *2.+C2s *2.+C2p *6.)*(10.-NKL)/10.;
              C12[N]    =(C1se*2.+C2se*2.+C2pe*6.)*(10.-NKL)/10.;

              D12_loc[N]=(D1s*2.+D2s*2.+D2p*6.)*(NKL)/10.;
              C3_loc[N] =(C3s*2.+C3p*6.+C3d*10.)*(18.-NM)/18.;
              C3    [N] =(C3se*2.+C3pe*6.+C3de*10.)*(18.-NM)/18.;
              D3_loc[N] =(D3s*2.+D3p*6.+D3d*10.)*(NM)/18.;

              D12[N]    =(D1se*2.+D2se*2.+D2pe*6.)*(NKL)/10.;
              D3[N]=(D3se*2.+D3pe*6.+D3de*10.)*(NM)/18.;
              E3[N]=2.*(e1s3s*2.+e1s3p*6.+e1s3d*10.);
              E3[N]=E3[N]+2.*(e2s3s*2.+e2s3p*6.+e2s3d*10.);
              E3[N]=(E3[N]+6.*(e2p3s*2.+e2p3p*6.+e2p3d*10.))*NKL*(18.-NM)/180.;
              D3_loc[N]=2.*(e1s3s*2.+(e1s3p+Rad3p1)*6.+e1s3d*10.);
              D3_loc[N]=D3_loc[N]+2.*(e2s3s*2.+(e2s3p+Rad3p2)*6.+e2s3d*10.);
              D3_loc[N]=D3_loc[N]+6.*((e2p3s+Rad3s)*2.+e2p3p*6.+(e2p3d+Rad3d)*10.);
              D3_loc[N]=D3_loc[N]*(10.-NKL)*(NM)/180.;
              DE3[N]=2.*(2.*e1s3s+2.*e2s3s+6.*(e2p3s+Rad3s))+
                   6.*(2.*(e1s3p+Rad3p1)+2.*(e2s3p+Rad3p2)+6.*e2p3p)
                   +10.*(2.*e1s3d+2.*e2s3d+6.*(e2p3d+Rad3d))*(10.-NKL)*(NM)/180.;
              RAD3[N]=(2.*(6.*Rad3s)+6.*(2.*Rad3p1+2.*Rad3p2)
                      +10.*6.*Rad3d)*(10.-NKL)*(NM)/180.;
              AKM3[N]=AKLM*NM*(y2s+y2p)*(2.-y1s);

                if (NM > 2)  ALM3[N]=ALMM*(NM)*(NM-1)*(8.-y2s-y2p) ;
                else         ALM3[N]=0.;

              RAD2[N]=0.;
              AKL2[N]=0.;
              PA2[N]=0.;
              PA23[N]=0.;
              PA3[N]=0.;
              EKLM4[N]=(y1s*e1s4+y2s*e2s4+y2p*e2p4
                        +y3s*e3s4+y3p*e3p4+y3d*e3d4);
              DEKLM4[N]=((2.-y1s)*e1s4+(2.-y2s)*e2s4+(6.-y2p)*e2p4
                        +(2.-y3s)*e3s4+(6.-y3p)*e3p4+(10.-y3d)*e3d4
                        +(2.-y3s)*Rad4+(6.-y3p)*Rad4+(10.-y3d)*Rad4);
              }
     } //continue;   L5
//==================================================================================
//     n=4 shell

 for(int N=0; N<=957; N++)
      {
      P14[N]=0.;
      C13[N]=0.;
      D13[N]=0.;
      C4M[N]=0.;
      D4M[N]=0.;
      E4M[N]=0.;
      DE4M[N]=0.;
      PA4[N]=0.;
      AKLM4[N]=0.;
      PA123[N]=0.;
      }

//     capture and ionization cross section in n=4 for correlated states (n123,n4),
//     as well as intershell exchanges, radiative and Auger transitions
//     ne123 = (electron nb in n=1+2+3)+1
//     ne4 = (electron nb in n=4)+1
//     N123 = state code {ne12,ne3} (in between 85 and 293)
//     NKLM = {ne12,ne3} state number (in between 1 and 209)
//     N4 = state code {ne4} (in between 294 and 326)
//     K = K+L electron number
//     M = M electron number
//     KLM = K+L+M electron number
//     NN = state number {ne123,ne4} (in between 1 and 957)
//     209 n12,n3 states, reduced to 29 constant n12+n3 states, combined to the
//     33 possible values of n4

if(EtachaVersion >= etacha_v4)
  {
  for(int  N123=85; N123<=293; N123++)   // L7
      {                  //      209 n1+2,n3 states, 29 possible values of ne123
      int M=f_IM(N123);
      int K=f_IKL(M,N123);
      int KLM=K+M;
      int NE123=K+M+1;
      // not used int L=100*M+K;
      int NKLM=N123-84;
      for(int N4=294; N4<=326; N4++)  // L7
           {                             //      33 states and 33 possible values of n4
           int NE4=N4-294+1;
           int ICOD3=100*(NE4-1)+(NE123-1);
           int NN = f_numPP(ICOD3)-326;

           double PT = Y[N123]*Y[N4];
           P14[NN] += PT;
                //c--------------------------------------------------------------
                //c    Cross sections for n=1+2+3,
                //c    and KLL,KLM,LMM autoionization probabilities (neKLM decreases by 1)
                //     NKLM is the {ne12,ne3} state number (in between 1 and 209)
                //c--------------------------------------------------------------
           C13[NN] += PT*(C12_loc[NKLM]+C3_loc[NKLM]);
           D13[NN] += PT*(D12_loc[NKLM]+D3_loc[NKLM]+AKL2[NKLM]+AKM3[NKLM]+ALM3[NKLM]);
           double PAKLMNN = PA2[NKLM]+PA23[NKLM]+PA3[NKLM];

           if (PAKLMNN > 1.)  PAKLMNN=1.;
           PA123[NN] += PT*PAKLMNN;
                //c--------------------------------------------------------------
                //c    Cross sections for n=4
                //     only depend on n4 value
                //     -> no need for normalization
                //c--------------------------------------------------------------
              C4M[NN] = (33-NE4)*C4;
              D4M[NN] = (NE4-1)*D4;
                //c--------------------------------------------------------------
                //c    Cross sections for n=1+2+3 <-> n=4
                //c--------------------------------------------------------------
              E4M[NN]  += PT*( EKLM4[NKLM]*(33-NE4));
              DE4M[NN] += PT*(DEKLM4[NKLM]*(NE4-1));
                //     (K or L or M) NN Augers : only depend on n4 and KLM values
                //     -> no need for normalization
               if ((NE4 > 2) && (KLM < 28))
                      {
                      AKLM4[NN]=AM4*(NE4-1)*(NE4-2)*(28-KLM);
                        //     autoionization probabilities of states (at the exit of target)
                      PA4[NN]  =AM4*(NE4-2)/(AM4*(NE4-2)+Rad4);
                      }
              else    {
                      PA4[NN]=0.;
                      }
          }} //continue; l7

        //c---------------------------------------------------------------
        //c    Normalization of weighted cross sections
        //c---------------------------------------------------------------
      for(int  N=1; N<=957; N++) //L8
        {
        int NUE=N+326;
        int NN = f_IN(NUE);
        int NKL= f_IKM(NN,NUE);
        double PKLMN = P14[N];
        if(PKLMN >= 0.000000000001)  {
              C13[N]   /= PKLMN;
              D13[N]   /= PKLMN;
              E4M[N]   /= PKLMN;
              DE4M[N]  /= PKLMN;
              PA123[N] /= PKLMN;
              }
        else {
              C13[N] =(C1s*2.+C2s*2.+C2p*6.+C3s*2.+C3p*6.+C3d*10.)*(28-NKL)/28.;
              D13[N] =(D1s*2.+D2s*2.+D2p*6.+D3s*2.+D3p*6.+D3d*10.)*(NKL)/28.;
              E4M[N] =(e1s4*2.+e2s4*2.+e2p4*6.+e3s4*2.+e3p4*6.+e3d4*10.)*NKL*(32-NN)/28.;
              DE4M[N]=(e1s4*2.+e2s4*2.+e2p4*6.+e3s4*2.+e3p4*6.+e3d4*10.)*(28-NKL)*NN/28.;
              PA123[N]=0.;
              }
        }
    }

}
//c-------------------------------------------------------------------
