//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SecMean(Y,T) {
      double PKLM[209],Y[1283],C12[209],D12[209],C3[209],D3[209],
           EKLM4[209],DEKLM4[209];
      common /etats/ICO1[63],ICO2[294],ICO3[1284],NUM1[733],NUM2[1911],
                NUM3[3329];
      common/seceff/sec[34],cor[48],secs[48];
      common/secKLM/C12e[209],D12e[209],RAD2[209],AKL2[209],PA2[209],
           AKM3[209],RAD3[209],ALM3[209],C3e[209],D3e[209],E3[209],DE3[209],
           PA3[209],PA23[209];
      common/sec1234/P14[957],C13[957],D13[957],C4M[957],D4M[957],
                E4M[957],DE4M[957],PA4[957],AKLM4[957],PA123[957];
      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
      common/aug/AKLL,AKLM,ALMM,AM4;
//c    ***********************************************
//c    calculation of actual (mean) cross sections
//c    ***********************************************
      C1s=sec[1];
      C2s=sec[2];
      C2p=sec[3];
      C3s=sec[4];
      C3p=sec[5];
      C3d=sec[6];
      C4=sec[7];
      D1s=sec[8];
      D2s=sec[9];
      D2p=sec[10];
      D3s=sec[11];
      D3p=sec[12];
      D3d=sec[13];
      D4=sec[14];
      e2s=sec[15];
      e2p=sec[16];
      e3s=sec[17];
      e3p=sec[18];
      e3d=sec[19];
      e4=sec[20];
      s3s=sec[21];
      s3p=sec[22];
      s3d=sec[23];
      s4=sec[24];
      p3s=sec[25];
      p3p=sec[26];
      p3d=sec[27];
      p4=sec[28];
      e3s4=sec[29];
      e3p4=sec[30];
      e3d4=sec[31];
      esp=sec[32];
      esp3=sec[33];
      epd3=sec[34];
//     add exchanges with n=4 for (1+2)*(3) fractions
      C1se=C1s+yn*e4;
      C2se=C2s+yn*s4;
      C2pe=C2p+yn*p4;
//     ynn is the mean value of n(n-1) in n=4 shell...
//     AM4 is about 3.*ALMM 
      C3se=C3s+yn*e3s4+AM4*ynn+Rad4*yn;
      C3pe=C3p+yn*e3p4+AM4*ynn+Rad4*yn;
      C3de=C3d+yn*e3d4+AM4*ynn+Rad4*yn;
      D1se=D1s+(32.-yn)*e4;
      D2se=D2s+(32.-yn)*s4;
      D2pe=D2p+(32.-yn)*p4;
      D3se=D3s+(32.-yn)*e3s4;
      D3pe=D3p+(32.-yn)*e3p4;
      D3de=D3d+(32.-yn)*e3d4;
      for(/*L3*/ N=1; N<=209; N++) {
      PKLM[N]=0.;
      C12[N]=0.;
      D12[N]=0.;
      C12e[N]=0.;
      D12e[N]=0.;
      RAD2[N]=0.;
      AKL2[N]=0.;
      PA2[N]=0.;
      RAD3[N]=0.;
      AKM3[N]=0.;
      ALM3[N]=0.;
      C3e[N]=0.;
      D3e[N]=0.;
      C3[N]=0.;
      D3[N]=0.;
      EKLM4[N]=0.;
      DEKLM4[N]=0.;
//     E3, DE3 = 1+2 <-> 3
      E3[N]=0.;
      DE3[N]=0.;
//     (KL)MM Autoionization :
      PA3[N]=0.;
//     KLM Autoionization :
      PA23[N]=0.;
L3:   } //continue;

//     capture and ionization cross section in n=3,
//     as well as intershell exchanges, radiatives and Auger transitions
//     ne12 = (electron nb in n=1+2)+1
//     ne3 = (electron nb in n=3)+1
//     NN = state code (ne12,ne3)-84, in between 1 and 209 (=11*19)
//     a loop is made on a given number of electrons in 1s,2s,2p subshells
//     (n12 code in between 1 and 63), and a given number of electrons 
//     in 3s,3p and 3d (L,LL,LLL): this determines one of the possible configurations
//     of a (ne12,ne3) state whose probability is given by Y(N)

      for(/*L4*/ N12=1; N12<=63; N12++) {
      for(/*L4*/ L=64; L<=66; L++) {
      for(/*L4*/ LL=67; LL<=73; LL++) {
      for(/*L4*/ LLL=74; LLL<=84; LLL++) {
      I=II[N12];
      J=JJ[I][N12];
      K=KK[I][J][N12];
      NE12=I+J+K+1;
      NE3=L+LL+LLL-205+1;
      ICOD2=100*(NE3-1)+(NE12-1);
      NN=NUMP[ICOD2]-84;

      PT=Y[N12]*Y[L]*Y[LL]*Y[LLL];
      PKLM[NN]=PKLM[NN]+PT;
//c--------------------------------------------------------------
//c    Cross sections for n=1+2
//c--------------------------------------------------------------
      C12[NN]=C12[NN]+PT*((2-K)*C1S+(2-J)*C2s+
           [6-I]*C2p);
      C12e[NN]=C12e[NN]+PT*((2-K)*C1Se+(2-J)*C2se+
           [6-I]*C2pe);
      D12[NN]=D12[NN]+PT*(K*D1s+J*D2s+I*D2p);
      D12e[NN]=D12e[NN]+PT*(K*D1se+J*D2se+I*D2pe);
        if (K < 2) {
       RAD2[NN]=RAD2[NN]+PT*(2-K)*I*Rad;
         if ((I+J) > 1) {
        AKL2[NN]=AKL2[NN]
                +PT*(2-K)*(I+J)*(I+J-1)*AKLL;
        PA2[NN]=PA2[NN]+PT*((2-K)*(I+J)*(I+J-1)*AKLL)/
                [(2-K]*(I+J)*(I+J-1)*AKLL+(2-K)*I*Rad);
        }
       }
//c--------------------------------------------------------------
//c    Cross sections for n=3
//c--------------------------------------------------------------
      C3[NN]=C3[NN]+PT*((66-L)*C3s+(73-LL)*C3p
                +(84-LLL)*C3d);
      C3e[NN]=C3e[NN]+PT*((66-L)*C3se+(73-LL)*C3pe
                +(84-LLL)*C3de);
      D3[NN]=D3[NN]+PT*((L-64)*D3s+(LL-67)*D3p
                +(LLL-74)*D3d);
      D3e[NN]=D3e[NN]+PT*((L-64)*D3se+(LL-67)*D3pe
                +(LLL-74)*D3de);
//c--------------------------------------------------------------
//c    Cross sections for n=1+2 <-> n=3
//c--------------------------------------------------------------
      E3[NN]=E3[NN]+PT*((66-L)*(K*e3s+J*s3s+I*p3s)
                +(73-LL)*(K*e3p+J*s3p+I*p3p)+(84-LLL)*(K*e3d+J*s3d+I*p3d));
      DE3[NN]=DE3[NN]+PT*(
           [L-64]*((2-K)*e3s+(2-J)*s3s+(6-I)*(p3s+rad3s))+
           [LL-67]*((2-K)*(e3p+rad3p1)+(2-J)*(s3p+rad3p2)+(6-I)*p3p)
           +(LLL-74)*((2-K)*e3d+(2-J)*s3d+(6-I)*(p3d+rad3d)));
      RAD3NN=((L-64)*(6-I)*rad3s
                +(LL-67)*((2-K)*rad3p1+(2-J)*rad3p2)
              +(LLL-74)*(6-I)*rad3d);
      AKM3NN=AKLM*(NE3-1)*(I+J)*(2-K);
//     ne3 is the (electron number in n=3) +1
//     I,J and K are the real number of electrons in n=2 and 1
       if (NE3 > 2) {
       ALM3NN=ALMM*(NE3-1)*(NE3-2)*(8-J-I);
// 
      }
 else {
       ALM3NN=0. ;
      }
      RAD3[NN]=RAD3[NN]+PT*RAD3NN;
      AKM3[NN]=AKM3[NN]+PT*AKM3NN;
      ALM3[NN]=ALM3[NN]+PT*ALM3NN;
      ATOT=ALM3NN+AKM3NN+RAD3NN;
       if(ATOT > 0.000000000001) {
       PA23[NN]=PA23[NN]+PT*(AKM3NN)/ATOT;
       PA3[NN]=PA3[NN]+PT*(ALM3NN)/ATOT;
      }
//c--------------------------------------------------------------
//c    Cross sections for n=1+2+3 <-> n=4
//c    unwheighted for n=4 pop (but dependant on 1+2+3 pop)
//c    PT is the probability of (n12,n3s,n3p,n3d) configuration
//c--------------------------------------------------------------
      EKLM4[NN]=EKLM4[NN]+PT*(K*e4+J*s4+I*p4
                +(L-64)*e3s4+(LL-67)*e3p4+(LLL-74)*e3d4);
      DEKLM4[NN]=DEKLM4[NN]+PT*((2-K)*e4+(2-J)*s4+(6-I)*p4
                +(66-L)*e3s4+(73-LL)*e3p4+(84-LLL)*e3d4
                +(66-L)*Rad4+(73-LL)*Rad4+(84-LLL)*Rad4);
L4:   } //continue;
//c---------------------------------------------------------------
//c    Normalization of wheighted cross sections
//c---------------------------------------------------------------
      for(/*L5*/ N=1; N<=209; N++) {
      NUE=N+84;
      NM=IM[NUE];
      NKL=IKL[NM][NUE];
      PKLMN=PKLM[N];
       if(PKLMN >= 0.000000000001) {
      C12[N]=C12[N]/PKLMN;
      D12[N]=D12[N]/PKLMN;
      C12e[N]=C12e[N]/PKLMN;
      D12e[N]=D12e[N]/PKLMN;
      RAD2[N]=RAD2[N]/PKLMN;
      AKL2[N]=AKL2[N]/PKLMN;
      PA2[N]=PA2[N]/PKLMN;
      AKM3[N]=AKM3[N]/PKLMN;
      RAD3[N]=RAD3[N]/PKLMN;
      ALM3[N]=ALM3[N]/PKLMN;
      C3[N]=C3[N]/PKLMN;
      D3[N]=D3[N]/PKLMN;
      C3e[N]=C3e[N]/PKLMN;
      D3e[N]=D3e[N]/PKLMN;
      E3[N]=E3[N]/PKLMN;
      DE3[N]=DE3[N]/PKLMN;
      PA3[N]=PA3[N]/PKLMN;
      PA23[N]=PA23[N]/PKLMN;
      EKLM4[N]=EKLM4[N]/PKLMN;
      DEKLM4[N]=DEKLM4[N]/PKLMN;
      }
 else {
//c    small probability configurations
      C12[N]=(C1s*2.+C2s*2.+C2p*6.)*(10.-NKL)/10.;
      D12[N]=(D1s*2.+D2s*2.+D2p*6.)*(NKL)/10.;
      C3[N]=(C3s*2.+C3p*6.+C3d*10.)*(18.-NM)/18.;
      D3[N]=(D3s*2.+D3p*6.+D3d*10.)*(NM)/18.;
      C12e[N]=(C1se*2.+C2se*2.+C2pe*6.)*(10.-NKL)/10.;
      D12e[N]=(D1se*2.+D2se*2.+D2pe*6.)*(NKL)/10.;
      C3e[N]=(C3se*2.+C3pe*6.+C3de*10.)*(18.-NM)/18.;
      D3e[N]=(D3se*2.+D3pe*6.+D3de*10.)*(NM)/18.;
      E3[N]=2.*(e3s*2.+e3p*6.+e3d*10.);
      E3[N]=E3[N]+2.*(s3s*2.+s3p*6.+s3d*10.);
      E3[N]=(E3[N]+6.*(p3s*2.+p3p*6.+p3d*10.))*NKL*(18.-NM)/180.;
      D3[N]=2.*(e3s*2.+(e3p+rad3p1)*6.+e3d*10.);
      D3[N]=D3[N]+2.*(s3s*2.+(s3p+rad3p2)*6.+s3d*10.);
      D3[N]=D3[N]+6.*((p3s+rad3s)*2.+p3p*6.+(p3d+rad3d)*10.);
      D3[N]=D3[N]*(10.-NKL)*(NM)/180.;
      DE3[N]=2.*(2.*e3s+2.*s3s+6.*(p3s+rad3s))+
           6.*(2.*(e3p+rad3p1)+2.*(s3p+rad3p2)+6.*p3p)
           +10.*(2.*e3d+2.*s3d+6.*(p3d+rad3d))*(10.-NKL)*(NM)/180.;
      Rad3[N]=(2.*(6.*rad3s)+6.*(2.*rad3p1+2.*rad3p2)
              +10.*6.*rad3d)*(10.-NKL)*(NM)/180.;
      AKM3[N]=AKLM*NM*(y2s+y2p)*(2.-y1s);
        if (NM > 2) { 
          ALM3[N]=ALMM*(NM)*(NM-1)*(8.-y2s-y2p) ;
       }
 else {
          ALM3[N]=0.;
       }
      RAD2[N]=0.;
      AKL2[N]=0.;
      PA2[N]=0.;
      PA23[N]=0.;
      PA3[N]=0.;
      EKLM4[N]=(y1s*e4+y2s*s4+y2p*p4
                +y3s*e3s4+y3p*e3p4+y3d*e3d4);
      DEKLM4[N]=((2.-y1s)*e4+(2.-y2s)*s4+(6.-y2p)*p4
                +(2.-y3s)*e3s4+(6.-y3p)*e3p4+(10.-y3d)*e3d4
                +(2.-y3s)*Rad4+(6.-y3p)*Rad4+(10.-y3d)*Rad4);
      }
L5:   } //continue;

//     n=4 shell

      for(/*L6*/ N=1; N<=957; N++) {
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
L6:   } //continue;

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

      for(/*L7*/ N123=85; N123<=293; N123++) {
//      209 n1+2,n3 states, 29 possible values of ne123
      M=IM[N123];
      K=IKL[M][N123];
      KLM=K+M;
      NE123=K+M+1;
      L=100*M+K;
      NKLM=N123-84;
      for(/*L7*/ N4=294; N4<=326; N4++) {
//      33 states and 33 possible values of n4
      NE4=N4-294+1;
      ICOD3=100*(NE4-1)+(NE123-1);
      NN=NUMPP[ICOD3]-326;

      PT=Y[N123]*Y[N4];
      P14[NN]=P14[NN]+PT;
//c--------------------------------------------------------------
//c    Cross sections for n=1+2+3,
//c    and KLL,KLM,LMM autoionization probabilities (neKLM decreases by 1)
//     NKLM is the {ne12,ne3} state number (in between 1 and 209)
//c--------------------------------------------------------------
      C13[NN]=C13[NN]+PT*(C12[NKLM]+C3[NKLM]);
      D13[NN]=D13[NN]+PT*(D12[NKLM]+D3[NKLM]+AKL2[NKLM]+AKM3[NKLM]
                +ALM3[NKLM]);
      PAKLMNN=PA2[NKLM]+PA23[NKLM]+PA3[NKLM];
       if (PAKLMNN > 1.) {
      PAKLMNN=1.;
      }
      PA123[NN]=PA123[NN]+PT*PAKLMNN;
//c--------------------------------------------------------------
//c    Cross sections for n=4
//     only depend on n4 value
//     -> no need for normalization
//c--------------------------------------------------------------
      C4M[NN]=(33-ne4)*C4;
      D4M[NN]=(ne4-1)*D4;
//c--------------------------------------------------------------
//c    Cross sections for n=1+2+3 <-> n=4
//c--------------------------------------------------------------
      E4M[NN]=E4M[NN]+PT*(EKLM4[NKLM]*(33-ne4));
      DE4M[NN]=DE4M[NN]+PT*(DEKLM4[NKLM]*(ne4-1));
//     (K or L or M) NN Augers : only depend on n4 and KLM values
//     -> no need for normalization
       if ((NE4 > 2) && [KLM < 28]) { 
      AKLM4[NN]=AM4*(NE4-1)*(NE4-2)*(28-KLM);
//     autoionization probabilities of states (at the exit of target)
      PA4[NN]=AM4*(NE4-2)/(AM4*(NE4-2)+Rad4);
      }
 else {
      PA4[NN]=0.;
      }
L7:   } //continue;
//c---------------------------------------------------------------
//c    Normalization of wheighted cross sections
//c---------------------------------------------------------------
      for(/*L8*/ N=1; N<=957; N++) {
      NUE=N+326;
      NN=IN[NUE];
      NKL=IKM[NN][NUE];
      PKLMN=P14[N];
       if(PKLMN >= 0.000000000001) {
      C13[N]=C13[N]/PKLMN;
      D13[N]=D13[N]/PKLMN;
      E4M[N]=E4M[N]/PKLMN;
      DE4M[N]=DE4M[N]/PKLMN;
      PA123[N]=PA123[N]/PKLMN;
      }
 else {
      C13[N]=(C1s*2.+C2s*2.+C2p*6.+C3s*2.+C3p*6.+C3d*10.)*(28-NKL)/28.;
      D13[N]=(D1s*2.+D2s*2.+D2p*6.+D3s*2.+D3p*6.+D3d*10.)*(NKL)/28.;
      E4M[N]=(e4*2.+s4*2.+p4*6.+e3s4*2.+e3p4*6.+e3d4*10.)
                *NKL*(32-NN)/28.;
      DE4M[N]=(e4*2.+s4*2.+p4*6.+e3s4*2.+e3p4*6.+e3d4*10.)
                *(28-NKL)*NN/28.;
      PA123[N]=0.;
      }
L8:   } //continue;
      return;
      }
//c-------------------------------------------------------------------