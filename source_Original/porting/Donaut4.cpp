//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void donaut(QP,Z1,AP,EP,Z2,AC,RHO,EPM,iflag,E1,S1,E2,S2,
                istp,iprt,ilgn); {
      double*8 zp8,zt8,E8;
      int ninis;
      double sec(34),cor(48),secs(48),seci(48),SeSE(66),StSE(12);
      common/seceff/sec,cor,secs;
      common/SecSE/SeSE,StSE;
      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
      common/tol/ep0,ep1,erel,erabs;
      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
      common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/don/zp,E,zt;
      common/corr/ibin;
      char kk*2,ll*2,kc*2, nomf*24;
//_______________________________________________________________________

//     table of cross sections (units 10-20cm2):
//     indices                 secs              indices                 sec

//           1                 MeC1s SEIK        1                 C1s
//           2                 MeC2s SEIK        2                 C2s
//           3                 MeC2p SEIK        3                 C2p
//           4                 MeC3s SEIK        4                 C3s
//           5                 MeC3p SEIK        5                 C3p
//           6                 MeC3d SEIK        6                 C3d
//           7                 MeC4  SEIK        7                 C4

//           8                 ReC1s                   8                 Ion1s
//           9                 ReC2s                   9                 Ion2s
//           10                ReC2p                   10                Ion2p
//           11                ReC3s                   11                Ion3s
//           12                ReC3p                   12                Ion3p
//           13                ReC3d                   13                Ion3d
//           14                ReC4                    14                Ion4

//           15                Ion1s CDW   EIS         15                Se1s2s
//           16                Ion2s CDW   EIS         16                Se1s2p
//           17                Ion2p CDW   EIS         17                Se1s3s
//           18                Ion3s Scaling           18                Se1s3p
//           19                Ion3p Scaling           19                Se1s3d
//           20                Ion3d Scaling           20                Se1s4
//           21                Ion4  Scaling 

//     Scaling means Ion(3l,Z) = Ion(1s,Z/3) and ion(4,Z) = Ion(2p,Z/2)
// 
//           22                Se1s5 SE          21                Se2s3s
//           23                Se2s5 SE          22                Se2s3p
//           24                Se2p5 SE          23                Se2s3d
//           25                Se3s5 SE          24                Se2s4
//           26                Se3p5 SE 
//           27                Se3d5 Scaling SE  25                Se2p3s
//           28                Se456 Scaling SE  26                Se2p3p
//                                                           27                Se2p3d
//           29                Se1s2s      SE          28                Se2p4
//           30                Se1s2p      SE 
//           31                Se1s3s      SE                29                Se3s4
//           32                Se1s3p      SE                30                Se3p4
//           33                Se1s3d      SE            31              Se3d4
//           34                Se1s4 SE 
// 
//           35                Se2s3s      SE                32                Se2s2p
//           36                Se2s3p      SE            33              Se3s3p
//           37                Se2s3d      SE          34                Se3p3d
//           38                Se2s4 SE 

//           39                Se2p3s      SE
//           40                Se2p3p      SE
//           41                Se2p3d      SE
//           42                Se2p4 SE

//           43                Se3s4 SE
//           44                Se3p4 SE
//           45                Se3d4 Scaling SE

//           46                Se2s2p
//           47                Se3s3p
//           48                Se3p3d
//_________________________________________________________________________________
//...... ibin=0 empirical saturation correction for born .....
//...... ibin=1 binding correction included in born (not recommended) .....
//...... ibin=2 no empirical correction and no binding correction
//...... correction are only for PWBA (BORN1) cross sections
//...... as calculated in PION, SEXI, SEX2 and SEX3,
//...... are only used for the evolution of cross sections with effective charge,
//...... and should not make much a change 

//...... Excitation to n=5 (secs(i), i=22,28) is not added to ionisation and set 
//...... to 0

//...... CDW is not used for capture cross sections 
//________________________________________________________________________________
      ibin=0;
      iflag=1;
      iter=0;
      y1s=0.;
      y2s=0.;
      y2p=0.;
      y3s=0.;
      y3p=0.;
      y3d=0.;
      y4=0.;
      yl=0.;
      ym=0.;
      ymp=0.;
      ymm=0.;
      ykm=0.;
      yl1m=0.;
      yl2m=0.;
      ym1m=0.;
      ym2m=0.;
      tetak=1.;
      tetal1=1.;
      tetal2=1.;
      tetam1=1.;
      tetam2=1.;
      tetan=1.;
      nomf = 'files\etadon.etacha';
      fopen(UNIT=13,ERR=1003,STATUS='old',FILE=nomf);
      fscanf(13,11,err=1000) QP,ZP,AP,zt,AC;
//L11:  format(4(1x,f4.0),1x,f6.2/);
      fscanf(13,12,err=1000) E,RHO,EPM;
//L12:  format(2(1x,F8.3),1x,F10.3,/);
      fscanf(13,15,err=1000) E1,S1,E2,S2,ISTP;
//L15:  format(4(1x,F8.3),1xI2,/);
      fscanf(13,16,err=1000) ep0,ep1,erel,erabs;
//L16:  format(2(1x,g10.3),2(1x,g11.4),/);
      fscanf(13,17,err=1000) iprt,ilgn;
//L17:  format(2(1x,I2),/);
      CLOSE(UNIT=13);
      printf(6,*)'      double /*data*/ used in previous calculation :';
      printf(6,*)' ';
L500: printf(6,*)'      General double /*data*/ :';
      printf(6,*)' ';
      printf (6,20)zp,qp,ap,e,zt,ac,rho;
//L20:  format(2x26hPROJECTILE: atomic number=,f4.0,2x11hincident ch,
                5harge=,f4.0,2x12hatomic mass=f4.0,/,12x16hincident energy=,
                f8.3,6h MeV/u,/,6x22hTARGET: atomic number=,f4.0,2x6hatomic,
                6h mass=,f4.0,2x8hdensity=,f6.3,6h g/cm3);
      printf (6,21)epm,ep0,ep1,erabs,erel;
//L21:  format(2x35hmaximum target thickness (mg/cm2) =,f10.3,/,
                2x23hminimum step (µg/cm2) =,g10.3,2x13hmaximum step=,g10.3,
                /,2x35hnumerical uncertainties on output :,/,2x9habsolute=,
                g11.4,2x9hrelative=,g11.4,/);
      printf (6,*) ' Want to change any of these value? (No/y)';
      printf (6,*) ' (Type ''y'' to change';
      printf (6,*) '  Type ''return'' not to change) ';
      fscanf (5,'(A)') KK;
//...........................................
       if ((KK == 'Y') || (KK == 'y')) {
L501:   printf (6,22) zp;
//L22:  format(2x26hProjectile atomic number =,f4.0," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          iter=0;
          printf (6,*) ' New value?';
          fscanf (5,*) Zp;
        }
        printf (6,23) qp;
//L23:  format(2x19hProjectile charge =,f4.0," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) QP;
        }
        printf (6,24) Ap;
//L24:  format(2x17hProjectile mass =,f4.0," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          iter=0;
          printf (6,*) ' New value? ';
          fscanf (5,*) AP;
        }
        printf (6,25) E;
//L25:  format(2x25hIncident energy (MeV/u) =,f8.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          iter=0;
          printf (6,*) ' New value? ';
          fscanf (5,*) E;
        }
        printf (6,26) Zt;
//L26:  format(2x22hTarget atomic number =,f4.0," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          iter=0;
          printf (6,*) ' New value? ';
          fscanf (5,*) zt;
        }
        printf (6,27) Ac;
//L27:  format(2x20hTarget atomic mass =,f4.0," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          iter=0;
          printf (6,*) ' New value? ';
          fscanf (5,*) Ac;
        }
        printf (6,28) Rho;
//L28:  format(2x24hTarget density (g/cm3) =,f6.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) RHO;
        }
        printf (6,29) Epm;
//L29:  format(2x35hmaximum target thickness (mg/cm2) =,f10.3,
                " Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) EPM;
        }
        printf (6,30) Ep0;
//L30:  format(2x23hminimum step (µg/cm2) =,g10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Ep0;
        }
        printf (6,31) Ep1;
//L31:  format(2x23hmaximum step (µg/cm2) =,g10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Ep1;
        }
      printf(6,32) Erabs,Erel;
//L32:  format(2x46hMaximun numerical uncertainties on populations,
                28h P(i) (=erel*(Max(P,erabs))),/,2x9habsolute=,
                g11.4,2x9hrelative=,g11.4," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value for erabs? ';
          fscanf (5,*) ERABS;
          printf (6,*) ' New value for erel? ';
          fscanf (5,*) EREL;
        }
        printf(6,*)' New values :';
        printf(6,*)' ';
        printf (6,20)zp,qp,ap,e,zt,ac,rho;
        printf (6,21)epm,ep0,ep1,erabs,erel;
        printf (6,*) ' Want to change again any of these value?',
                  ' (No/y) ';
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          goto L501;
        }
      }
//.......................................................................
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
printf(6,*)' for( you want this PROGRAM to modify cross ') { {
      printf(6,*)'sections when the projectile ion loss';
      printf(6,*)' energy in thick targets? (No/y) ';
        fscanf (5,'(A)') KK;
//.................................................
       if ((KK == 'Y') || (KK == 'y')) {
      istp=1;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
printf(6,*)' This PROGRAM needs two values of energy',
                ' in order to perform this task'; {
      printf(6,*)' (You can choose one value close to the initial',
          ' energy and one close to ';
      printf(6,*)' the final one)';
      printf (6,*) ' Previous values :';
//c.........
      printf (6,33)E1,S1;
      printf (6,34)E2,S2;
//L33:  format(2x12hE1 (MeV/u) =,f8.3,2x17hS1 (MeV/mg/cm2) =,f8.3);
//L34:  format(2x12hE2 (MeV/u) =,f8.3,2x17hS2 (MeV/mg/cm2) =,f8.3);
      printf (6,*) ' Change these values ? (No/y) ';
      fscanf (5,'(A)') KK;
//c.........
       if ((KK == 'Y') || (KK == 'y')) {
L502:   printf (6,35)E1,S1;
//L35:    format(2x12hE1 (MeV/u) =,f8.3,2x17hS1 (MeV/mg/cm2) =,f8.3,
                  " Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value for E1 (MeV/u)? ';
          printf (6,*)' ( if your input for E1 is 0 there will be NO',
                  ' stopping correction)';
          fscanf (5,*) E10;
           if (E10 == 0.) {
          istp=0;
          goto L503;
          }
 else {
          E1=E10;
            call zstop(zp,E1,zt,S1);
          }
            printf (6,*) ' Calculated value : ';
            printf (6,351)E1,S1;
//L351:   format(2x12hE1 (MeV/u) =,f8.3,2x17hS1 (MeV/mg/cm2) =,f8.3);
          printf (6,*) ' Change value for S1? (No/y)';
            fscanf (5,'(A)') KK;
             if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value for S1 (MeV/mg/cm2)? ';
          fscanf (5,*) S1;
            }
        }
        printf (6,36)E2,S2;
//L36:    format(2x12hE2 (MeV/u) =,f8.3,2x17hS2 (MeV/mg/cm2) =,f8.3,
                  " Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value for E2 (MeV/u)? ';
          fscanf (5,*) E2;
            call zstop(zp,E2,zt,S2);
            printf (6,*) ' Calculated value : ';
            printf (6,361)E2,S2;
//L361:   format(2x12hE2 (MeV/u) =,f8.3,2x17hS2 (MeV/mg/cm2) =,f8.3);
          printf (6,*) ' Change value for S2? (No/y)';
            fscanf (5,'(A)') KK;
             if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value for S2 (MeV/mg/cm2)? ';
          fscanf (5,*) S2;
            }
        }
L503:   printf(6,*)' new values for stopping power correction:';
         if (istp == 0) {
        printf(6,*)' no stopping power correction ...';
        }
 else {
        printf (6,33)E1,S1;
        printf (6,34)E2,S2;
        }
        printf (6,*) ' Change these values ? (No/y) ';
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
        goto L502;
        }
      }
      }
 else {
        istp=0;
      }
//---------------------------------------------------
//     CROSS SECTIONS 
//---------------------------------------------------
      printf (6,*) ' ';
      printf (6,*) ' please, wait......';
      printf (6,*) ' Hydrogenic (reference) cross sections',
                ' calculation ....';
       if (iter == 0) {
        for(/*L50*/ I=1; I<=48; I++) {
          cor(I)=1.;
L50:    } //continue;
        call CSEC(E);
        for(/*L51*/ I=1; I<=48; I++) {
          seci(I)=secs(I);
L51:    } //continue;
//c
//     Preliminary calculation of scaling factors for excitation cross sections
//     of d and f states in SE approximation @ v=zp

      Ev=E;
      bet=zp/137.036;
      E=931.5*(1/pow((1-bet*bet),0.5)-1);
      printf (6,*)' E=',E;
      call sex3(e3s4,e3p4,e3d4,e3s5,e3p5,e3d5,e4s5,e4p5,e45,e456);
      R3d4=(e3d4)/(e3s4+e3p4);
      R3d5=(e3d5)/(e3s5+e3p5);
//     e456 is not corrected for 1/n**3 law anymore...
      R45=e456/(e4s5+e4p5);
      E=Ev;
//c
      printf (6,*) ' ';
      printf(6,*)' CDWEIS Ionisation Cross sections';
      ninis=10;
      call tceis(tot,ninis);
      secs(15)=tot*1.e20;
      ninis=20;
      call tceis(tot,ninis);
      secs(16)=tot*1.e20;
      ninis=21;
      call tceis(tot,ninis);
      secs(17)=tot*1.e20;

      zp3=zp;
      zp=zp3/3.;
      ninis=10;
      call tceis(tot,ninis);
      secs(18)=tot*1.e20;
      secs(19)=tot*1.e20;
      secs(20)=tot*1.e20;
      zp=zp3;
// 
      zp2=zp;
      zp=zp2/2.;
      ninis=21;
      call tceis(tot,ninis);
      secs(21)=tot*1.e20;
      zp=zp2;
// 
      cor(15)=secs(15)/seci(15);
      cor(16)=secs(16)/seci(16);
      cor(17)=secs(17)/seci(17);
      cor(18)=secs(18)/seci(18);
      cor(19)=secs(19)/seci(19);
      cor(20)=secs(20)/seci(20);
      cor(21)=secs(21)/seci(21) ;
// 
      printf(6,*)' SE Excitation Cross sections';
      zp8=dble(zp);
      zt8=dble(zt);
      E8=dble(E);
      call SEnlm(zp8,zt8,E8);

//     Excitation to n>=5 set to 0 in the present version
// 
//     1s - 5
      secs(22)=StSE(6)*1.e20*0.;
//     2s,2p - 5 
      secs(23)=StSE(7)*1.e20*0.;
      secs(24)=StSE(8)*1.e20*0.;
//     3s,3p,3d - 5 
      secs(25)=StSE(9)*1.e20*0.;
      secs(26)=StSE(10)*1.e20*0.;
      secs(27)=R3d5*(secs(25)+secs(26));
//     R45 does not include correction for the sum on 1/n**3 anymore
      secs(28)=R45*(StSE(11)+StSE(12))*1.e20*0.;
// 
      secs(29)=SeSE(1)*1.e20;
      secs(30)=SeSE(2)*1.e20;
      secs(31)=SeSE(3)*1.e20;
      secs(32)=SeSE(4)*1.e20;
      secs(33)=SeSE(5)*1.e20;
      secs(34)=StSE(1)*1.e20;
      secs(35)=SeSE(15)*1.e20;
      secs(36)=SeSE(16)*1.e20;
      secs(37)=SeSE(17)*1.e20;
      secs(38)=StSE(2)*1.e20;
      secs(39)=SeSE(27)*1.e20;
      secs(40)=SeSE(28)*1.e20;
      secs(41)=SeSE(29)*1.e20;
      secs(42)=StSE(3)*1.e20;
      secs(43)=StSE(4)*1.e20;
      secs(44)=StSE(5)*1.e20;
      secs(45)=R3d4*(secs(43)+secs(44));

//     Excitation to n=5 and 6 set to 0 

      cor(22)=0.;
      cor(23)=0.;
      cor(24)=0.;
      cor(25)=0.;
      cor(26)=0.;
      cor(27)=0.;
      cor(28)=0.;
// 
      cor(29)=secs(29)/seci(29);
      cor(30)=secs(30)/seci(30);
      cor(31)=secs(31)/seci(31);
      cor(32)=secs(32)/seci(32);
      cor(33)=secs(33)/seci(33);
      cor(34)=secs(34)/seci(34);
      cor(35)=secs(35)/seci(35);
      cor(36)=secs(36)/seci(36);
      cor(37)=secs(37)/seci(37);
      cor(38)=secs(38)/seci(38);
      cor(39)=secs(39)/seci(39);
      cor(40)=secs(40)/seci(40);
      cor(41)=secs(41)/seci(41);
      cor(42)=secs(42)/seci(42);
      cor(43)=secs(43)/seci(43);
      cor(44)=secs(44)/seci(44);
      cor(45)=secs(45)/seci(45);

//     note that secs(22) -> secs(28) are set to 0...
// 
      sec(8)=secs(15)+secs(22);
      sec(9)=secs(16)+secs(23);
      sec(10)=secs(17)+secs(24);
      sec(11)=secs(18)+secs(25);
      sec(12)=secs(19)+secs(26);
      sec(13)=secs(20)+secs(27);
      sec(14)=secs(21)+secs(28);
// 
      sec(15)=secs(29);
      sec(16)=secs(30);
      sec(17)=secs(31);
      sec(18)=secs(32);
      sec(19)=secs(33);
      sec(20)=secs(34);
      sec(21)=secs(35);
      sec(22)=secs(36);
      sec(23)=secs(37);
      sec(24)=secs(38);
      sec(25)=secs(39);
      sec(26)=secs(40);
      sec(27)=secs(41);
      sec(28)=secs(42);
      sec(29)=secs(43);
      sec(30)=secs(44);
      sec(31)=secs(45);
        iter=1;
      }
      printf (6,*) ' ';
      printf (6,40) (sec(I),I=1,33);
//L40:  format(2x49hCapture cross sections into fully stripped projec,
                8htile ...,/,2x43hIonization and excitation cross sections of,
                28h one electron projectile ...,/,2x22h(corrected for "satura,
                36htion" and screening + antiscreening),/,2x14h(all cross sec,
                21htions in 10E-20 cmý ),/,2x11hCAPTURE TO:,5x12h(capture inc,
                17hludes MEC + REC ),/,21x4h1s =,G10.3,2x4h2s =,G10.3,2x4h2p =,
                G10.3,/,21x4h3s =,G10.3,2x4h3p =,G10.3,2x4h3d =,G10.3,/,
                21x7h(n=4) =,G10.3,/,
             2x57hIONIZATION OF:  (ionization includes excitation to n > 4),
             /,21x4h1s =,G10.3,2x4h2s =,G10.3,2x4h2p =,G10.3,/,21x4h3s =,
                G10.3,2x4h3p =,G10.3,2x4h3d =,G10.3,/,21x7h(n=4) =,G10.3,/,
                2x17h1s EXCITATION TO:,2x4h2s =,
                G10.3,2x4h2p =,G10.3,/,21x4h3s =,G10.3,2x4h3p =,G10.3,2x4h3d =,
                G10.3,/,21x7h(n=4) =,G10.3,/,
                2x17h2s EXCITATION TO:,2x4h3s =,G10.3,2x4h3p =,
                G10.3,2x4h3d =,G10.3,/,21x7h(n=4) =,G10.3,/,
                2x17h2p EXCITATION TO:,2x4h3s =,G10.3,
                2x4h3p =,G10.3,2x4h3d =,G10.3,/,21x7h(n=4) =,G10.3,/,
                2x23hEXCITATION TO n=4 from:,2x4h3s =,G10.3,2x4h3p =,
                G10.3,2x4h3d =,G10.3,/,
                2x22hINTRASHELL EXCITATION:,
                2x10h2s to 2p =,G10.3,2x10h3s to 3p =,G10.3);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
      printf(6,*)'  During charge state calculation, cross sections',
                ' will be periodicaly adjusted as a function of mean charge',
                ' state ( if necessary) and as a function of energy  (if',
                ' desired).  if however you are not pleased with some of the',
                ' above values,you may now enter your own ones',
                '(fully stripped or one electron ion cross sections)',
                '   The PROGRAM will in this case scale his calculat',
                'ions to yours.'; {
      printf(6,*)'     for(/*L*/ you want to enter such corrected') {
      printf(6,*)' values ? (No/y) ';
      fscanf (5,'(A)') KK;
// ...................... debut des modifs ..........................
       if ((KK == 'Y') || (KK == 'y')) {
L504:   printf (6,*) ' Enter your cross sections in 10E-20 cm2';
        printf(6,*)'  CAPTURE to substates of fully stripped ions :';
      printf(6,*)' Want to change ? (No/y) ';
      fscanf (5,'(A)') KC;
       if ((KC == 'Y') || (KC == 'y')) {
        printf (6,101) secs(1);
//L101:   format(1x18h1s capture (MEC) =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
         printf (6,*) ' New value? ';
         fscanf (5,*) C1s;
         cor(1)=C1s/seci(1);
        }
        printf (6,102) secs(2);
//L102:   format(1x18h2s capture (MEC) =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) C2s;
          cor(2)=C2s/seci(2);
        }
        printf (6,103) secs(3);
//L103:   format(1x18h2p capture (MEC) =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) C2p;
          cor(3)=C2p/seci(3);
        }
        printf (6,104) secs(4);
//L104:   format(1x18h3s capture (MEC) =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) C3s;
          cor(4)=C3s/seci(4);
        }
        printf (6,105) secs(5);
//L105:   format(1x18h3p capture (MEC) =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) C3p;
          cor(5)=C3p/seci(5);
        }
        printf (6,106) secs(6);
//L106:   format(1x18h3d capture (MEC) =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) C3d;
          cor(6)=C3d/seci(6);
        }
        printf (6,107) secs(7);
//L107:   format(1x19hn=4 capture (MEC) =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
         printf (6,*) ' New value? ';
         fscanf (5,*) C4;
         cor(7)=C4/seci(7);
        }
        printf (6,108) secs(8);
//L108:   format(1x8h1s REC =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
         printf (6,*) ' New value? ';
         fscanf (5,*) SRK;
         cor(8)=SRK/seci(8);
        }
        printf (6,109) secs(9);
//L109:   format(1x8h2s REC =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) SRls;
          cor(9)=SRls/seci(9);
        }
        printf (6,110) secs(10);
//L110:   format(1x8h2p REC =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) SRlp;
          cor(10)=SRlp/seci(10);
        }
        printf (6,111) secs(11);
//L111:   format(1x8h3s REC =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) SR3s;
          cor(11)=SR3s/seci(11);
        }
        printf (6,112) secs(12);
//L112:   format(1x8h3p REC =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) SR3p;
          cor(12)=SR3p/seci(12);
        }
        printf (6,113) secs(13);
//L113:   format(1x8h3d REC =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) SR3d;
          cor(13)=SR3d/seci(13);
        }
        printf (6,114) secs(14);
//L114:   format(1x9hn=4 REC =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) SR4;
          cor(14)=SR4/seci(14);
        }
      }
        printf(6,*)'  IONIZATION of hydrogenlike ions :';
      printf(6,*)' Want to change ? (No/y) ';
      fscanf (5,'(A)') KC;
       if ((KC == 'Y') || (KC == 'y')) {
        printf (6,115) secs(15);
//L115:   format(1x15h1s Ionization =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Dk;
          cor(15)=Dk/seci(15);
        }
        printf (6,116) secs(16);
//L116:   format(1x15h2s Ionization =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Dl1;
          cor(16)=Dl1/seci(16);
        }
        printf (6,117) secs(17);
//L117:   format(1x15h2p Ionization =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Dl2;
          cor(17)=Dl2/seci(17);
        }
        printf (6,118) secs(18);
//L118:   format(1x15h3s Ionization =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Dm1;
          cor(18)=Dm1/seci(18);
        }
        printf (6,119) secs(19);
//L119:   format(1x15h3p Ionization =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Dm2;
          cor(19)=Dm2/seci(19);
        }
        printf (6,120) secs(20);
//L120:   format(1x15h3d Ionization =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Dm3;
          cor(20)=Dm3/seci(20);
        }
        printf (6,121) secs(21);
//L121:   format(1x16hn=4 Ionization =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Dn;
          cor(21)=Dn/seci(21);
        }
      }
        printf (6,*)' excitation to n>4 to be added to ionization :';
      printf(6,*)' Want to change ? (No/y) ';
      fscanf (5,'(A)') KC;
       if ((KC == 'Y') || (KC == 'y')) {
        printf (6,122) secs(22);
//L122:   format(1x22h1s -> n>4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) es5;
          cor(22)=es5/seci(22);
        }
        printf (6,123) secs(23);
//L123:   format(1x22h2s -> n>4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e2s5;
          cor(23)=e2s5/seci(23);
        }
        printf (6,124) secs(24);
//L124:   format(1x22h2p -> n>4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e2p5;
          cor(24)=e2p5/seci(24);
        }
        printf (6,125) secs(25);
//L125:   format(1x22h3s -> n>4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e3s5;
          cor(25)=e3s5/seci(25);
        }
        printf (6,126) secs(26);
//L126:   format(1x22h3p -> n>4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e3p5;
          cor(26)=e3p5/seci(26);
        }
        printf (6,127) secs(27);
//L127:   format(1x22h3d -> n>4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e3d5;
          cor(27)=e3d5/seci(27);
        }
        printf (6,128) secs(28);
//L128:   format(1x23hn=4 -> n>4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e45;
          cor(28)=e45/seci(28);
        }
      }
        printf(6,*)'  EXCITATION of hydrogenlike ions :';
        printf(6,*)' Excitation from 1s :';
      printf(6,*)' Want to change ? (No/y) ';
      fscanf (5,'(A)') KC;
       if ((KC == 'Y') || (KC == 'y')) {
        printf (6,129) secs(29);
//L129:   format(1x21h1s -> 2s excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) E2s;
          cor(29)=E2s/seci(29);
        }
        printf (6,130) secs(30);
//L130:   format(1x21h1s -> 2p excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) E2p;
          cor(30)=E2p/seci(30);
        }
        printf (6,131) secs(31);
//L131:   format(1x21h1s -> 3s excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) E3s;
          cor(31)=E3s/seci(31);
        }
        printf (6,132) secs(32);
//L132:   format(1x21h1s -> 3p excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) E3p;
          cor(32)=E3p/seci(32);
        }
        printf (6,133) secs(33);
//L133:   format(1x21h1s -> 3d excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) E3d;
          cor(33)=E3d/seci(33);
        }
        printf (6,134) secs(34);
//L134:   format(1x22h1s -> n=4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) es4;
          cor(34)=es4/seci(34);
        }
      }
        printf(6,*)' Excitation from n=2 :';
      printf(6,*)' Want to change ? (No/y) ';
      fscanf (5,'(A)') KC;
       if ((KC == 'Y') || (KC == 'y')) {
        printf (6,135) secs(35);
//L135:   format(1x21h2s -> 3s excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) s3s;
          cor(35)=s3s/seci(35);
        }
        printf (6,136) secs(36);
//L136:   format(1x21h2s -> 3p excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) s3p;
          cor(36)=s3p/seci(36);
        }
        printf (6,137) secs(37);
//L137:   format(1x21h2s -> 3d excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) s3d;
          cor(37)=s3d/seci(37);
        }
        printf (6,138) secs(38);
//L138:   format(1x22h2s -> n=4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e2s4;
          cor(38)=e2s4/seci(38);
        }
        printf (6,139) secs(39);
//L139:   format(1x21h2p -> 3s excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) p3s;
          cor(39)=p3s/seci(39);
        }
        printf (6,140) secs(40);
//L140:   format(1x21h2p -> 3p excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) p3p;
          cor(40)=p3p/seci(40);
        }
        printf (6,141) secs(41);
//L141:   format(1x21h2p -> 3d excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) p3d;
          cor(41)=p3d/seci(41);
        }
        printf (6,142) secs(42);
//L142:   format(1x22h2p -> n=4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e2p4;
          cor(42)=e2p4/seci(42);
        }
      }
        printf(6,*)' Excitation from n=3 :';
      printf(6,*)' Want to change ? (No/y) ';
      fscanf (5,'(A)') KC;
       if ((KC == 'Y') || (KC == 'y')) {
        printf (6,143) secs(43);
//L143:   format(1x22h3s -> n=4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e3s4;
          cor(43)=e3s4/seci(43);
        }
        printf (6,144) secs(44);
//L144:   format(1x22h3p -> n=4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e3p4;
          cor(44)=e3p4/seci(44);
        }
        printf (6,145) secs(45);
//L145:   format(1x22h3d -> n=4 excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) e3d4;
          cor(45)=e3d4/seci(45);
        }
      }
        printf(6,*)'  INTRASHELL EXCITATION of hydrogenlike ions :';
      printf(6,*)' Want to change ? (No/y) ';
      fscanf (5,'(A)') KC;
       if ((KC == 'Y') || (KC == 'y')) {
        printf (6,146) secs(46);
//L146:   format(1x21h2s -> 2p excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Esp;
          cor(46)=Esp/seci(46);
        }
        printf (6,147) secs(47);
//L147:   format(1x21h3s -> 3p excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Esp3;
          cor(47)=Esp3/seci(47);
        }
        printf (6,148) secs(48);
//L148:   format(1x21h3p -> 3d excitation =,G10.3," Change? (No/y) ");
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          printf (6,*) ' New value? ';
          fscanf (5,*) Epd3;
          cor(48)=Epd3/seci(48);
        }
      }
//..........................................
        for(/*L52*/ I=1; I<=48; I++) {
        secs(I)=cor(I)*seci(I);
L52:    } //continue;
        for(/*L53*/ I=1; I<=7; I++) {
//     add MEC and REC
        sec(i)=secs(i)+secs(i+7);
//     add ionisation and nl -> 5 excitation
        sec(i+7)=secs(i+14)+secs(i+21);
L53:    } //continue;
//     shift index for other processes (excitation)
      for(/*L54*/ i=1; i<=20; i++) {
        sec(i+14)=secs(i+28);
L54:  } //continue;
//..........................................
        printf (6,*) ' ';
        printf (6,40) (sec(I),I=1,33);
        printf (6,*) ' Want to change again any of these value? ';
        printf (6,*) '(No/y) ';
        fscanf (5,'(A)') KK;
         if ((KK == 'Y') || (KK == 'y')) {
          goto L504;
        }
      }
// ...................... fin des modifs ..........................
      printf (6,*) ' Want to recheck all input values? (No/y) ';
      fscanf (5,'(A)') KK;
       if ((KK == 'Y') || (KK == 'y')) {
        goto L500;
      }
      fopen(UNIT=13,ERR=1001,STATUS='old',FILE='files\etadon.etacha');
      printf(13,11,err=1001) QP,ZP,AP,ZT,AC;
      printf(13,12,err=1001) E,RHO,EPM;
      printf(13,15,err=1001) E1,S1,E2,S2,ISTP;
      printf(13,16,err=1001) ep0,ep1,erel,erabs;
      printf(13,17,err=1001) iprt,ilgn;
      close(UNIT=13);
      z1=zp;
      ep=e;
      z2=zt;
      Epm=1000.*EPM;
      goto L1002;
L1003: printf(6,*) 'Problem in opening double /*data*/ file... check', nomf;
      iflag=2;
      goto L1002;
L1000: printf(6,*) 'Problem in reading data... check', nomf;
      iflag=3;
      goto L1002;
L1001: printf(6,*) 'Problem in writing data... check', nomf;
      iflag=2;
L1002: } //continue;
      return;
      }
//c      ************** fin de donnees ***************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void CSEC(EF) {
      double sec(34),cor(48),secs(48);
      common/seceff/sec,cor,secs;
      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
      common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn;
      common/don/zp,E,zt;
      common/corr/ibin;
      QM=Zp-(y1s+yl+ym+yn);
      E=EF;
      VC=pow((1-pow((1+E/931.5),(-2))),(0.5));
      vu=vc*137.036;
      printf(6,2000)E,QM;
//L2000: format(33h Cross sections calculation at E=,g12.5,7h(MeV/u),
                3x8h, Qmean=,f7.3);
      printf(6,*)'ionization';
      call PION(dk,dl1,dl2,dm1,dm2,dn);
      printf(6,*)'capture';
      call seik(ck,cl,cm,cn,ct,cs);
      sig14=ck+cl+cm+cn;
       if (sig14 <= cs) {
      sig530=cs-sig14;
      }
 else {

//     ct is the sum of seik from 5 to 30
// 
      sig530=ct;
      }
      call rec(srk,srls,srlp);
      srm=8./27.*(srls+srlp);
      stm=cm+srm;
      srn=1./8.*(srls+srlp);
//     MEC
      secs(1)=ck;
      secs(2)=.25*cl;
      secs(3)=.75*cl;
      secs(4)=cm*1./9.;
      secs(5)=cm*3./9.;
      secs(6)=cm*5./9.;
//************************ 
      secs(7)=cn;
//     if (zp.gt.60.) then
//     secs(7)=secs(7)+sig530
//     endif
//************************ 
//     REC
      secs(8)=srk;
      secs(9)=srls;
      secs(10)=srlp;
      secs(11)=srls*8./27.;
      secs(12)=srlp*8./27.;
      secs(13)=srlp*8./45.;
      secs(14)=srn;
//     EXCITATION
      printf(6,*)'excitation 1s-nl';
      call sexi(e2s,e2p,e3s,e3p,e3d,es4,es5);
      printf(6,*)'excitation 2l-nl';
      call sex2(s3s,s3p,s3d,p3s,p3p,p3d,e2s4,e2p4,e2s5,e2p5);
      printf(6,*)'excitation 3l-nl';
      call sex3(e3s4,e3p4,e3d4,e3s5,e3p5,e3d5,e4s5,e4p5,e45,e456);
      printf(6,*)'excitation nl-nl';
      call snl(esp,esp3,epd3);
       if ((ibin == 1) || (ibin == 2)) {
//c.......... binding correction included (1) or no correction (2) ..........
      secs(15)=Dk;
      secs(16)=Dl1;
      secs(17)=Dl2;
      secs(18)=Dm1;
      secs(19)=Dm1;
      secs(20)=Dm2;
      secs(21)=Dn;
//c excitation to n=5 will be added to ionization
//     set to 0
      secs(22)=es5*0.;
//     (es5 includes a sum up to infinity)
//     secs(23)=e2s5*3.05
//     secs(24)=e2p5*3.05
      secs(23)=e2s5*0.;
      secs(24)=e2p5*0.;
//     (theoretical factor in 1/n**3 = 3.049358 for n=5)
//     secs(25)=e3s5*2.5
//     secs(26)=e3p5*2.5
//     secs(27)=e3d5*2.5
      secs(25)=e3s5*0.;
      secs(26)=e3p5*0.;
      secs(27)=e3d5*0.;
//     (3 ->5 => 2.5 instead of 3.05)
      secs(28)=e456*0.;
//     (e456=e45+e46 instead of 
//     e456=e45+2.5*e46 -> 2.5 au lieu de 3.5)
//.....................
//c excitation from 1s
      secs(29)=e2s;
      secs(30)=e2p;
      secs(31)=e3s;
      secs(32)=e3p;
      secs(33)=e3d;
      secs(34)=es4;
//c excitation from n=2
      secs(35)=s3s;
      secs(36)=s3p;
      secs(37)=s3d;
      secs(38)=e2s4;
      secs(39)=p3s;
      secs(40)=p3p;
      secs(41)=p3d;
      secs(42)=e2p4;
//c excitation from n=3 to n=4
      secs(43)=e3s4;
      secs(44)=e3p4;
      secs(45)=e3d4;
//...........................
      }
 else {
//c....... empirical saturation correction (09/09/94) ...........
//c    not used anymore since a long time ...
      cse=1.06;
      csi=0.735;
      zs1=pow((10.96*pow(vu,2)*exp(0.111*zp)/pow((zp),1.946)),0.5);
      zs2=pow((10.96*pow(vu,2)*exp(0.111*zp/1.5)/pow((zp/2.),1.946)),0.5);
      zs3=pow((10.96*pow(vu,2)*exp(0.111*zp/1.5)/pow((zp/3.),1.946)),0.5);
      sate1=pow((zs1/zt),2)*exp(-cse*zt/pow(vu,2.1));
      sate1=sate1*pow((1.-exp(-pow((zt/zs1),2.5))),0.8);
      sate2=pow((zs2/zt),2)*exp(-cse*zt/pow(vu,2.1));
      sate2=sate2*pow((1.-exp(-pow((zt/zs2),2.5))),0.8);
      sate3=pow((zs3/zt),2)*exp(-cse*zt/pow(vu,2.1));
      sate3=sate3*pow((1.-exp(-pow((zt/zs3),2.5))),0.8);
      sati1=pow((zs1/zt),2)*exp(-csi*zt/pow(vu,2.1));
      sati1=sati1*pow((1.-exp(-pow((zt/zs1),2.5))),0.8);
      sati2=pow((zs2/zt),2)*exp(-csi*zt/pow(vu,2.1));
      sati2=sati2*pow((1.-exp(-pow((zt/zs2),2.5))),0.8);
      sati3=pow((zs3/zt),2)*exp(-csi*zt/pow(vu,2.1));
      sati3=sati3*pow((1.-exp(-pow((zt/zs3),2.5))),0.8);
//c......................................................
      secs(15)=Dk*sati1;
      secs(16)=Dl1*sati2;
      secs(17)=Dl2*sati2;
      secs(18)=Dm1*sati3;
      secs(19)=Dm1*sati3;
      secs(20)=Dm2*sati3;
      secs(21)=Dn;
//c excitation to n=5 to be added to ionization
//     set to 0
      secs(22)=es5*sate1*0.;
//     (es5 includes a sum up to infinity)
      secs(23)=e2s5*3.05*sate2*0.;
      secs(24)=e2p5*3.05*sate2*0.;
//     (theoretical factor in 1/n**3 = 3.049358 for n=5)
      secs(25)=e3s5*2.5*sate3*0.;
      secs(26)=e3p5*2.5*sate3*0.;
      secs(27)=e3d5*2.5*sate3*0.;
//     (3 ->5 => 2.5 in place of 3.05)
      secs(28)=e456*0.;
//     (e456=e45+2.5*e46 -> 2.5 in place of 3.5)
//c excitation from 1s
      secs(29)=e2s*sate1;
      secs(30)=e2p*sate1;
      secs(31)=e3s*sate1;
      secs(32)=e3p*sate1;
      secs(33)=e3d*sate1;
      secs(34)=es4*sate1;
//c excitation a partir de n=2
      secs(35)=s3s*sate2;
      secs(36)=s3p*sate2;
      secs(37)=s3d*sate2;
      secs(38)=e2s4*sate2;
      secs(39)=p3s*sate2;
      secs(40)=p3p*sate2;
      secs(41)=p3d*sate2;
      secs(42)=e2p4*sate2;
//c excitation a partir de n=3 vers n=4
      secs(43)=e3s4*sate3;
      secs(44)=e3p4*sate3;
      secs(45)=e3d4*sate3;
//c
      }
//c intracouche
      secs(46)=esp;
      secs(47)=esp3;
      secs(48)=epd3;
//......................
      for(/*L2001*/ I=1; I<=48; I++) {
        secs(I)=cor(I)*secs(I);
L2001:   } //continue;
      for(/*L2002*/ I=1; I<=7; I++) {
//  capture totale (1 a 7)
        sec(i)=secs(i)+secs(i+7);
//  perte totale (8 a 14)
        sec(i+7)=secs(i+14)+secs(i+21);
L2002: } //continue;
//  excitation (15 a 34)
      for(/*L2003*/ i=1; i<=20; i++) {
        sec(i+14)=secs(i+28);
L2003: } //continue;
      return;
      }
//c	************** fin de Csec ***************************