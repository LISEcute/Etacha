//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void Auger(Y,Zp,PR,YTOT) {
//c..........................................
//c    Calculates final charge state distribution for each target thickness
//c    of the output files;
//c    takes into account autoionisation at the exit of target
//c..........................................
      double PR[62],P1234[29][33],Y[1283],YA[1283];
      common /etats/ICO1[63],ICO2[294],ICO3[1284],NUM1[733],NUM2[1911],
                NUM3[3329];
      common/secKLM/C12[209],D12[209],RAD2[209],AKL2[209],PA2[209],
           AKM3[209],RAD3[209],ALM3[209],C3[209],D3[209],E3[209],DE3[209],
           PA3[209],PA23[209];
      common/sec1234/P14[957],C13[957],D13[957],C4M[957],D4M[957],
                E4M[957],DE4M[957],PA4[957],AKLM4[957],PA123[957];
      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
      common/aug/AKLL,AKLM,ALMM,AM4;
      YTOT=0.;
      for(/*L1*/ N=327; N<=1283; N++) {
      YA[N]=Y[N];
      YTOT=YTOT+Y[N];
L1:   } //continue;
//     Calculation of correlated states (n123,n4) probabilities after autoionisation
      for(/*L3*/ N1=1; N1<=29; N1++) {
      for(/*L3*/ N2=1; N2<=33; N2++) {
      P1234[N1][N2]=0.;
L3:   } //continue;
//.....................
      for(/*L4*/ N1=1; N1<=29; N1++) {
      for(/*L4*/ N4=33; N4<=1; N4+-1) {
      I=N1-1;
      M=N4-1;
      INM=100*M+I;
      N=NUMPP[INM];
//     N stands for indices of correlated states (between 327 and 1283)
//     effect of KLL,KLM and LMM (PA123) Auger effects on these states
//     (gives rise to very little change in most cases)
      L=N-326;
       if(I > 2) {
            P1234[I+1][M+1]=P1234[I+1][M+1]+(1.-PA123[L])*YA[N];
            P1234[I][M+1]=P1234[I][M+1]+PA123[L]*YA[N];
      }
 else {
            P1234[I+1][M+1]=P1234[I+1][M+1]+YA[N];
      }
L4:   } //continue;
//....................
//     effet of (KLM)NN (PA4)) Augers effect on P1234 fractions (as calculated above)
//     fractions with N-1 e- KLM and M-1 e- N initially
//     fractions with 28 e- KLM (N=29) or no N shell e- (M=1) are unchanged
      for(/*L5*/ N=1; N<=28; N++) {
//     (states with at least one hole in K or L or M shells)
      for(/*L5*/ M=33; M<=2; M+-1) {
        INM=100*(M-1)+N-1;
        NM=NUMPP[INM]-326;
         if (M > 2) {
//     radiative decay 
            P1234[N+1][M-1]=P1234[N+1][M-1]+(1.-PA4[NM])*P1234[N][M];
//     Auger decay
            P1234[N+1][M-2]=P1234[N+1][M-2]+PA4[NM]*P1234[N][M];
//     NO Auger decay (as a test)
//           P1234(N+1,M-1)=P1234(N+1,M-1)+PA4(NM)*P1234(N,M)

        P1234[N][M]=0.;

        }
 else if (M == 2) {
//     only one N shell electron :
            P1234[N+1][M-1]=P1234[N+1][M-1]+P1234[N][M];
        P1234[N][M]=0.;
        }
L5:   } //continue;
//..................
//     conbine configurations and calculates final probabilities
      for(/*L6*/ N=1; N<=62; N++) {
       PR[N]=0.;
L6:   } //continue;
      for(/*L7*/ N=1; N<=29; N++) {
      for(/*L7*/ M=1; M<=33; M++) {
      Ne1=N+M-1;
       if (Ne1 <= (Zp+1)) {
      PR[Ne1]=PR[Ne1]+P1234[N][M];
      }
 else {
      PR[Zp+1]=PR[Zp+1]+P1234[N][M];
      }
L7:   } //continue;
      for(/*L100*/ N=1; N<=61; N++) {
      PP=Pr[N];
       if(PP < 0.0000001) {
      PP=0.;
      PR[N]=PP;
      }
L100: } //continue;
      return;
      }

