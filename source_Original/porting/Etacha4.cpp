//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
PROGRAM etacha {
//c*****************************************************************************
//c***     This version of ETACHA is for ions with up to 60 électrons        ***
//c***     First PC version of ETACHA was 18/12/91                           ***
//c***     First distributed PC version of ETACHA was 07/96                  ***
//c*****************************************************************************
//c
//c   05/2015 version going up to nmax=4 with partially correlated states
//c   Last corrections 09/2015 
//c 
//c   NOTES :    1) COA = 1-(1.75*qmin/vu)**2
//c                  2) CDW calculations are not used any more, SEIK instead
//c                  3) CDWEIS +SC & ASC calculations for ionisation 
//c                  4) SE + SC & ASC calculations for excitation
//c              5) Excitation to n >= 5 not added to ionization anymore    // Oleg may be make it optional?
//c                 (See Donaut4.for)
//c
//c   Contact: J.P.ROZET  eMail: rozet@insp.jussieu.fr
//c            D. Vernhet eMail: dominique.vernhet@insp.jussieu.fr
//c            E. Lamour  eMail: lamour@insp.jussieu.fr
//c
//c       Calculates charge states of swift ions and their evolution as a fonction
//c    of traversed target thickness for a number of initial electrons =< 60.
//c    In a first step, only 1s,2s et 2p states are considered (as in oldest ETACHA),
//c    and 3s,3p and 3d state populations independently estimated.
//c    This program then solves first a set of 84 coupled differential equations : 
//c    63 are for the evolution of the 63 "correlated states of the type Y(i,j,k),
//c    where i,j and k stand for the number of 2p,2s and 1s electrons
//c    (i between 0 et 6, j and k between 0 and 2), and 21 are for the evolution
//c    of 3s (3 equations), 3p (7 equations) and 3d (11 equations) substates.
//c    Differential equations are integrated using a Runge-Kutta type method and associated subroutines
//c    (see EQDIFF.for, an improvement of Adams code - a variable order predictor-corrector method)
//c 
//c    In the actual version, the evolution of n<4 states is calculated in an improved way, 
//c    by considering the  11*19=209 Y(n12,n3) type states, where n12 is the electron number in n=1 and 2,
//c    and n3 the electron number in n=3. For this calculation, actual (averadged) cross sections,
//c    functions of n=n12+n3 are used (see SecMean.for).
//c    This corresponds to a total of 84+209=293 equations (this is the ETACHA3 version).
//c
//c    A final evolution is calculated by considering first in an "indépendant" way the n=4 shell (33 equations),
//c    and then finally  the 29*33 Y(n123,n4)type states, where n123 is the total electron number in n<4 and n4 is
//c    the electron number in n=4.
//c
//c__________________________________________________________________
//c          We then have eventually 293+33+29*33=293+990=1283 equations.
//c__________________________________________________________________
//c 
//c    Files to be included for building the exe file (15 files):
//c          -Etacha4.for (present file)
//c          -Auger4.for
//c          -Donaut4.for
//c          -EQDIF.for
//c          -F4.for
//c          -INTG.for
//c          -Pion.for (PWBA ionisation cross sections)
//c          -SecMean4.for
//c          -SEIK.for (capture cross sections in the SE approximation + REC)
//c          -senlm.for (excitation cross sections in teh SE approximation)
//c          -Sex2.for (PWBA excitation cross sections from n=2 and n=3)
//c          -Sexi.for (PWBA excitation cross sections from n=1)
//c          -Snl.for (PWBA intrashell excitation cross sections)
//c          -Tceis.for (CDW-EIS ionisation cross sections for n=1 and 2)
//c          -Zstop.for (Ziegler's SRIM derived stopping power)
//c    Files to be included as data files (4 files)
//c          -Etadon.dat   ==> etadon.etacha
//c          -SCOEFgas.95, SCOEF.95A, SCOEF.95B
//c    Output files :
//c          -ETA09.prn, ETA1019.prn, (ETA2029.prn, ETA3039.prn, ETA4049.prn and ETA5059.prn if usefull)
//c          -ETAPied.prn, POPMean.prn, seceff.prn
//c
//c
//c    ***********************************************************
      double Y[1284],YX[1284],WORK[1680000],IWORK[1350],PR[62],
                P1s[2],P2s[2],P2p[2],PRF[62];
//     LENW=N*N+17*N+204=1670688
//     LENIW=N+21=1305 
      char nomf*24;
      common /etats/ICO1[63],ICO2[294],ICO3[1284],NUM1[733],NUM2[1911],
                NUM3[3329];
      common/seceff/sec[34],cor[48],secs[48];
      common/secKLM/C12[209],D12[209],RAD2[209],AKL2[209],PA2[209],
           AKM3[209],RAD3[209],ALM3[209],C3[209],D3[209],E3[209],DE3[209],
           PA3[209],PA23[209];
      common/sec1234/P14[957],C13[957],D13[957],C4M[957],D4M[957],
                E4M[957],DE4M[957],PA4[957],AKLM4[957],PA123[957];
      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
      common/tol/ep0,ep1,erel,erabs;
      common/aug/AKLL,AKLM,ALMM,AM4;
      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
      common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m;
      common/popi/n02s,n02p,n03s,n03p,n03d,nl0,nm0;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn;
      common/proj/zp;
      external F;
//c    ***********************************************************
//c    initializations and data input procedures
//c    ***********************************************************
      ilgn=64;
      iprt=1;
      NEQ=1283;
      iflag=1;
      icor=0;
      call ini[ICO1][ICO2][ICO3][NUM1][NUM2][NUM3];
L20:  call donaut[QP][Z1][AP][EP][Z2][AC][RHO][EPM][iflag][E1][S1][E2][S2][
                istp][iprt][ilgn];
       if (iflag == 2) {
      printf (6,*) 'iflag=2 !! Input double /*data*/ file problem';
      goto L5000;
      }
      ilp=ilgn-55;
      for(/*L21*/ i=1; i<=38; i++) {
         if (cor[i] != 1.) {
          icor=1;
        }
L21:  } //continue;
       if (icor == 1) {
      ilp=ilp-2;
      }
      pas0=20.*ep0;
      zp=z1;
      ip=zp;
      e=ep;
      zc=z2;
      VC=pow((1-pow((1+E/931.5),(-2))),(0.5));
      VU=137.036*VC;
      ok=omk[zp];
      ol=oml[zp];
//c....... slowing down coefficients ..........
       if(istp == 0) {
        tc=1.e9;
        as=0.;
        bs=e;
//c.......!!!
      }
 else {
        as=(s2-s1)/((e2-e1)*ap);
        bs=(s2*e1-s1*e2)/(s2-s1);
        stp=s1+(e-E1)*(s2-s1)/(e2-e1);
         if(stp < 0.) {
        printf(6,*)' ';
        printf(6,*)' ';
        printf(6,*)' stp=',stp;
        STOP' problem with stopping power values, } of calculation';
        }
        stp1=log10[ap*e/stp];
        astp=pow(10.,aint[stp1]);
         if (stp1 < 0.) {
          astp=astp/10.;
        }
        sstp=ap*e/(stp*astp);
         if (sstp <= 1.5) {
          tc=2.5*astp;
        }
 else if (sstp <= 3.5) {
          tc=5.*astp;
        }
 else if (sstp < 7.) {
          tc=10.*astp;
        }
 else {
          tc=25.*astp;
        }
      }
       if (zp <= 6) {
        ddq=8.;
         if (erel > 1.e-4) {
          ddq=ddq/2.;
        }
      }
 else if (zp <= 20) {
        ddq=4.;
         if (erel > 1.e-4) {
          ddq=ddq/2.;
        }
      }
 else if (zp <= 40) {
        ddq=2.;
         if (erel > 1.e-4) {
          ddq=ddq/2.;
        }
      }
 else {
        ddq=1.;
      }
      printf (6,*) 'initialization of populations';
//c    ***** number of electrons in shells of neutral atom ***********
      n02s=0;
      n02p=0;
      n03s=0;
      n03p=0;
      n03d=0;
      n04s=0;
      n04p=0;
      n04d=0;
      n04f=0;
       if (IP <= 2) {
      } //continue;
      }
 else if (IP <= 4) {
      n02s=IP-2;
      }
 else if (IP <= 10) {
      n02s=2;
      n02p=IP-4;
      }
 else if (IP <= 12) {
      n02s=2;
      n02p=6;
      n03s=IP-10;
      }
 else if (IP <= 18) {
      n02s=2;
      n02p=6;
      n03s=2;
      n03p=IP-12;
      }
 else if (IP <= 20) {
      n02s=2;
      n02p=6;
      n03s=2;
      n03p=6;
      n04s=IP-18;
      }
 else if (IP <= 28) {
      n02s=2;
      n02p=6;
      n03s=2;
      n03p=6;
      n04s=2;
      n03d=IP-20;
      }
 else if (IP <= 36) {
      n02s=2;
      n02p=6;
      n03s=2;
      n03p=6;
      n04s=2;
      n03d=10;
      n04p=IP-30;
      }
 else if (IP <= 38) {
      n02s=2;
      n02p=6;
      n03s=2;
      n03p=6;
      n04s=2;
      n03d=10;
      n04p=6;
      n05s=IP-36;
      }
 else if (IP <= 48) {
      n02s=2;
      n02p=6;
      n03s=2;
      n03p=6;
      n04s=2;
      n03d=10;
      n04p=6;
      n05s=2;
      n04d=IP-38;
      }
 else if (IP <= 54) {
      n02s=2;
      n02p=6;
      n03s=2;
      n03p=6;
      n04s=2;
      n03d=10;
      n04p=6;
      n05s=2;
      n04d=10;
      n05p=IP-48;
      }
 else if (IP <= 56) {
      n02s=2;
      n02p=6;
      n03s=2;
      n03p=6;
      n04s=2;
      n03d=10;
      n04p=6;
      n05s=2;
      n04d=10;
      n05p=6;
      n06s=IP-54;
      }
 else {
      n02s=2;
      n02p=6;
      n03s=2;
      n03p=6;
      n04s=2;
      n03d=10;
      n04p=6;
      n05s=2;
      n04d=10;
      n05p=6;
      n06s=2;
      n04f=IP-56;
      }
      nl0=n02s+n02p;
      nm0=n03s+n03p+n03d;
      nn0=n04s+n04p+n04d+n04f;
//c    ***** Y(n) initialization *****************************
      for(/*L1*/ n=1; n<=1283; n++) {
      Y[n]=0.;
L1:   } //continue;
//c    ***** number of electrons in shells of incident ion ***********
      NEL=Zp-Qp;
      n1s=0;
      n2s=0;
      n2p=0;
      n3s=0;
      n3p=0;
      n3d=0;
      n4s=0;
      n4p=0;
      n4d=0;
      n4f=0;
       if (NEL <= 2) {
      n1s=NEL;
      }
 else if (NEL <= 4) {
      n1s=2;
      n2s=NEL-2;
      }
 else if (NEL <= 10) {
      n1s=2;
      n2s=2;
      n2p=NEL-4;
      }
 else if (NEL <= 12) {
      n1s=2;
      n2s=2;
      n2p=6;
      n3s=NEL-10;
      }
 else if (NEL <= 18) {
      n1s=2;
      n2s=2;
      n2p=6;
      n3s=2;
      n3p=NEL-12;
      }
 else if (NEL <= 28) {
      n1s=2;
      n2s=2;
      n2p=6;
      n3s=2;
      n3p=6;
      n3d=NEL-18;
      }
 else if (NEL <= 30) {
      n1s=2;
      n2s=2;
      n2p=6;
      n3s=2;
      n3p=6;
      n3d=10;
      n4s=NEL-28;
      }
 else if (NEL <= 36) {
      n1s=2;
      n2s=2;
      n2p=6;
      n3s=2;
      n3p=6;
      n3d=10;
      n4s=2;
      n4p=NEL-30;
      }
 else if (NEL <= 46) {
      n1s=2;
      n2s=2;
      n2p=6;
      n3s=2;
      n3p=6;
      n3d=10;
      n4s=2;
      n4p=6;
      n4d=NEL-36;
      }
 else if (NEL <= 60) {
      n1s=2;
      n2s=2;
      n2p=6;
      n3s=2;
      n3p=6;
      n3d=10;
      n4s=2;
      n4p=6;
      n4d=10;
      n4f=NEL-46;
      }
 else {
      printf (6,*) ' ';
      printf (6,*) '**************************************************';
      printf (6,*) ' too much electrons ...';
      printf (6,*) ' try again with ZP-Q < 61';
      printf (6,*) '**************************************************';
      printf (6,*) ' ';
      goto L20;
      }
      nkl=n1s+n2s+n2p;
      nm=n3s+n3p+n3d;
      nklm=nkl+nm;
      n4l=n4s+n4p+n4d+n4f;
      Ic1=n2p*100+n2s*10+n1s;
      Ic2=nm*100+nkl;
      n1=num[Ic1];
      n2=numP[Ic2];
      printf (6,*) 'initialization of populations for actual charge';
      Y[n1]=1.;
      Y[64+n3s]=1.;
      Y[67+n3p]=1.;
      Y[74+n3d]=1.;
      Y[n2]=1.;
      Y[294+n4l]=1.;
      n1234=327+nklm+29*n4l;
      Y[n1234]=1.;
//c    ******************************************************************
//c    *** radiative and Auger "cross sections" (units 10e-20 cm²) ***
      conv=6.023E-3/AC;
      RAD=3.461E-6*pow((Zk*tetak),4)*AC/(RHO*VC)/2.;
      RAD3s=.034891E-6*pow((Zk*tetak),4)*AC/(RHO*VC)/6.;
      RAD3p=1.0301E-6*pow((Zk*tetak),4)*AC/(RHO*VC)/2.;
      RAD3d=.3544E-6*pow((Zk*tetak),4)*AC/(RHO*VC)/6.;
      RAD4=0.0454E-6*pow((Zk*tetak),4)*AC/(RHO*VC)/18.;
       if (n02p >= 1) {
      AKLL=rad*(1./ok-1.)*n02p/(nl0*(nl0-1));
      }
 else {
      AKLL=7.6E-2*AC/(RHO*VC);
      }
       if (nm0 >= 1) {
      AKLM=rad3p*0.882*(1./ok-1.)*n02p/(nl0*nm0);
      }
 else {
      AKLM=0.3*AC/(RHO*VC);
      }
       if (nm0 >= 2) {
      ALMM=(rad3s*n03s+0.118*rad3p*n03p+rad3d*n03d)*(1./ol-1.)
           /(nm0*(nm0-1));
      }
 else {
      ALMM=0.05*AC/(RHO*VC);
      }
      AM4=3.*ALMM;
//c      ***************************************************************
//c      Open and write headings of output files
//c      ***************************************************************
      call datetime[iyr][imon][iday][ihour][imin][isec][imil];
//c ...................................................................
//L200:    format(2x,i2.2,1h/,i2.2,1h/,i4,21x,7hEtacha4,13x,
                3h(C),2x,9hINSP-ASUR,2x,11hJPR 05/2015);
//L202:    format(89('*'));
//L203:    format(118('*'));
//L204:    format(2x,"PROJECTILE: atomic number=",f4.0,2x,"incident ch",
             "harge=",f4.0,2x,"atomic mass=",f4.0,/,12x,"incident energy=",
             f8.3," MeV/u",2x,"velocity=",f8.3," (au)",/,
             6x,"TARGET: atomic number=",f4.0,2x,"atomic",
             " mass=",f4.0,2x,"density=",f6.3," g/cm3");
//L205:    format(2x,"relative error=",g11.4,3x,"absolute=",g11.4);
//L206:    format(2x,"No stopping power correction",/);
//L207:    format(2x,"Energy loss: S1=",f8.3," MeV/mg/cm² at E1=",f8.3,
           "MeV/u",/,10x,"and: S2=",f8.3," MeV/mg/cm² at E2=",f8.3," MeV/u") ;
//L208: format(2x,"Reference ('hydrogenic') cross sections in 10E-20 cm²:",
              /,"  CAPTURE (MEC+REC) TO:  ",
              " 1s =",G10.3,"  2s =",G10.3,"  2p =",G10.3,/,25x,
              " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
              "  (n=4) =",G10.3,/,
                "  IONIZATION OF:",9x,
              " 1s =",G10.3,"  2s =",G10.3,"  2p =",G10.3,/,25x,
              " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
              "  (n=4) =",G10.3,/,
                2x,"1s EXCITATION TO:",
              23x,            "2s =",G10.3,"  2p =",G10.3,/,25x,
              " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
                "  (n=4) =",G10.3,/,
                2x,"2s EXCITATION TO:",6x,
                " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
                "  (n=4) =",G10.3,/,
              2x,"2p EXCITATION TO:",6x,
              " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
              "  (n=4) =",G10.3,/,
                2x,"EXCITATION To n=4 from:",
                " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,/,
                2x,"INTRASHELL EXCITATION:",
                "  2s to 2p =",G10.3,2x,"3s to 3p =",G10.3);
//L309: format (\,"  Binding Energy correction factors =",6(2x,F8.4));
//L310: format (\,"  Effective charges =",6(2x,F8.4)) ;
//L311: format ("  Binding energies  =",6(2x,F8.1)) ;
//c ....................................................................
//L210:    format('Tar.thick.    0e-       1e-       2e-       3e-',
            '       4e-       5e-       6e-       7e-       8e-       9e-',
            '     Eout');
//L211:    format(1x,"(ug/cm²) ",10(3x,"(",i2,"+)  "));
//L215:    format(1x,"(ug/cm²) ",10(3x,"(",i2,"+)  "));
//L212:    format('Tar.thick.    10e-      11e-      12e-      13e-',
            '      14e-      15e-      16e-      17e-      18e-      19e-');
//L213:    format('Tar.thick.    20e-      21e-      22e-      23e-',
            '      24e-      25e-      26e-      27e-      28e-      29e-');
//L214:    format('Tar.thick.    30e-      31e-      32e-      33e-',
            '      34e-      35e-      36e-      37e-      38e-      39e-');
//L218:    format('Tar.thick.    40e-      41e-      42e-      43e-',
            '      44e-      45e-      46e-      47e-      48e-      49e-');
//L219:    format('Tar.thick.    50e-      51e-      52e-      53e-',
            '      54e-      55e-      56e-      57e-      58e-      59e-');
//L216:    format('T (ug/cm²)  bare    1s     2s     2p    1s²',
            '    1s2s   1s2p  1s²2s  1s²2p   tot');
//L217:    format('T (ug/cm²)  y1s       y2s      y2p      ym',
            '      yn      Qm      Qm in   Qm out    PTOT');
//c ...................................................................
      nomf='results\ETA0009.PRN';
      fopen (10,file=nomf);
      printf (10,200) imon,iday,iyr;
      printf (10,202);
      printf (10,204) zp,qp,ap,e,vu,zc,ac,rho;
      printf (10,205) erel,erabs;
       if (istp == 0) {
          printf (10,206);
      }
 else {
          printf (10,207) S1,E1,S2,E2;
      }
      printf (10,208) (sec[I],I=1,33);
      printf (10,203);
      printf (10,210);
      printf (10,211)ip,ip-1,ip-2,ip-3,ip-4,ip-5,ip-6,ip-7,ip-8,ip-9;
      printf (10,203);

      nomf='results\ETA1019.PRN';
      fopen (11,file=nomf);
      printf (11,200) imon,iday,iyr;
      printf (11,202);
      printf (11,204) zp,qp,ap,e,vu,zc,ac,rho;
      printf (11,205) erel,erabs;
       if (istp == 0) {
          printf (11,206);
      }
 else {
        printf (11,207) S1,E1,S2,E2;
      }
      printf (11,208) (sec[I],I=1,33);
      printf (11,203);
      printf (11,212);
       if (ip >= 19) {
        printf (11,211)ip-10,ip-11,ip-12,ip-13,ip-14,ip-15,ip-16,ip-17,ip-18,
           ip-19;
      }
 else {
        printf (11,*) ' ';
      }
      printf (11,203);

       if (ip >= 20) {
      nomf='results\ETA2029.PRN';
      fopen (12,file=nomf);
      printf (12,200) imon,iday,iyr;
      printf (12,202);
      printf (12,204) zp,qp,ap,e,vu,zc,ac,rho;
      printf (12,205) erel,erabs;
       if (istp == 0) {
          printf (12,206);
      }
 else {
          printf (12,207) S1,E1,S2,E2;
      }
      printf (12,208) (sec[I],I=1,33);
      printf (12,203);
      printf (12,213);
       if (ip >= 29) {
        printf (12,211)ip-20,ip-21,ip-22,ip-23,ip-24,ip-25,ip-26,ip-27,ip-28,
           ip-29;
      }
 else {
        printf (12,*) ' ';
      }
      printf (12,203);
      }

       if (ip >= 30) {
      nomf='results\ETA3039.PRN';
      fopen (13,file=nomf);
      printf (13,200) imon,iday,iyr;
      printf (13,202);
      printf (13,204) zp,qp,ap,e,vu,zc,ac,rho;
      printf (13,205) erel,erabs;
       if (istp == 0) {
          printf (13,206);
      }
 else {
          printf (13,207) S1,E1,S2,E2;
      }
      printf (13,208) (sec[I],I=1,33);
      printf (13,203);
      printf (13,214);
       if (ip >= 39) {
        printf (13,215)ip-30,ip-31,ip-32,ip-33,ip-34,ip-35,ip-36,ip-37,ip-38,
           ip-39;
      }
 else {
        printf (13,*) '   ';
      }
      printf (13,203);
      }
// 
       if (ip >= 40) {
      nomf='results\ETA4049.PRN';
      fopen (20,file=nomf);
      printf (20,200) imon,iday,iyr;
      printf (20,202);
      printf (20,204) zp,qp,ap,e,vu,zc,ac,rho;
      printf (20,205) erel,erabs;
       if (istp == 0) {
          printf (20,206);
      }
 else {
          printf (20,207) S1,E1,S2,E2;
      }
      printf (20,208) (sec[I],I=1,33);
      printf (20,203);
      printf (20,218);
       if (ip >= 49) {
        printf (20,215)ip-40,ip-41,ip-42,ip-43,ip-44,ip-45,ip-46,ip-47,ip-48,
           ip-49;
      }
 else {
        printf (20,*) '   ';
      }
      printf (20,203);
      }
// 
       if (ip >= 50) {
      nomf='results\ETA5059.PRN';
      fopen (21,file=nomf);
      printf (21,200) imon,iday,iyr;
      printf (21,202);
      printf (21,204) zp,qp,ap,e,vu,zc,ac,rho;
      printf (21,205) erel,erabs;
       if (istp == 0) {
          printf (21,206);
      }
 else {
          printf (21,207) S1,E1,S2,E2;
      }
      printf (21,208) (sec[I],I=1,33);
      printf (21,203);
      printf (21,219);
       if (ip >= 59) {
        printf (21,215)ip-50,ip-51,ip-52,ip-53,ip-54,ip-55,ip-56,ip-57,ip-58,
           ip-59;
      }
 else {
        printf (21,*) '   ';
      }
      printf (21,203);
      }

      nomf='results\ETAPIED.PRN';
      fopen (14,file=nomf);
      printf (14,200) imon,iday,iyr;
      printf (14,202);
      printf (14,204) zp,qp,ap,e,vu,zc,ac,rho;
      printf (14,205) erel,erabs;
       if (istp == 0) {
        printf (14,206);
      }
 else {
        printf (14,207) S1,E1,S2,E2;
      }
      printf (14,208) (sec[I],I=1,33);
      printf (14,203);
      printf (14,*) ' ';
      printf (14,216);
      printf (14,203);

      nomf='results\POPMEAN.PRN';
      fopen (15,file=nomf);
      printf (15,200) imon,iday,iyr;
      printf (15,202);
      printf (15,204) zp,qp,ap,e,vu,zc,ac,rho;
      printf (15,205) erel,erabs;
       if (istp == 0) {
        printf (15,206);
      }
 else {
        printf (15,207) S1,E1,S2,E2;
      }
      printf (15,208) (sec[I],I=1,33);
      printf (15,202);
      printf (15,*) '  ';
      printf (15,217);
      printf (15,203);
//c
      nomf='results\seceff.PRN';
      fopen (16,file=nomf);
      printf (16,200) imon,iday,iyr;
      printf (16,202);
      printf (16,204) zp,qp,ap,e,vu,zc,ac,rho;
      printf (16,205) erel,erabs;
       if (istp == 0) {
        printf (16,206);
      }
 else {
        printf (16,207) S1,E1,S2,E2;
      }
      printf (16,208) (sec[I],I=1,33);
      printf (16,309) tetak,tetal1,tetal2,tetam1,tetam2,tetan;
      printf (16,310) zk,zl1,zl2,zm1,zm2,zn;
      printf (16,311) Bk,Bl1,Bl2,Bm1,Bm2,Bn;
      printf (16,202);
//c
//c    *******************************************************
//c    calculates cross sections (incident ion) in reduced units (µg/cm²)-1
//c    *******************************************************
      call POPMEAN[Y];
      QM=Zp-(y1s+yl+ym+yn);
      call CSEC[E];
      printf (16,*) 'QM=',QM,'E=',E;
      printf (16,208) (sec[I],I=1,33);
      printf (16,309) tetak,tetal1,tetal2,tetam1,tetam2,tetan;
      printf (16,310) zk,zl1,zl2,zm1,zm2,zn;
      printf (16,311) Bk,Bl1,Bl2,Bm1,Bm2,Bn;
      printf (16,203);
//L220: format(1x,f8.3,1x,f7.3,6(1x,G10.3));
      for(/*L2*/ I=1; I<=34; I++) {
      Sec[I]=Sec[I]*conv;
L2:   } //continue;
//c
//c    cross sections per electron and per hole
//c
      Sec[1]=sec[1]/2.;
      sec[2]=sec[2]/2.;
      sec[3]=sec[3]/6.;
      sec[4]=sec[4]/2.;
      sec[5]=sec[5]/6.;
      sec[6]=sec[6]/10.;
      sec[7]=sec[7]/32.;
      sec[15]=sec[15]/2.;
      sec[16]=sec[16]/6.;
      sec[17]=sec[17]/2.;
      sec[18]=sec[18]/6.;
      sec[19]=sec[19]/10.;
      sec[20]=sec[20]/32.;
      sec[21]=sec[21]/2.;
      sec[22]=sec[22]/6.;
      sec[23]=sec[23]/10.;
      sec[24]=sec[24]/32.;
      sec[25]=sec[25]/2.;
      sec[26]=sec[26]/6.;
      sec[27]=sec[27]/10.;
      sec[28]=sec[28]/32.;
      sec[29]=sec[29]/32.;
      sec[30]=sec[30]/32.;
      sec[31]=sec[31]/32.;
      sec[32]=sec[32]/6.;
      sec[33]=sec[33]/6.;
      sec[34]=sec[34]/10.;
//c 
      Rad=rad*conv;
      Rad3s=Rad3s*conv;
      Rad3p1=Rad3p*0.882*conv;
      Rad3p2=Rad3p*0.118*conv;
      Rad3d=Rad3d*conv;
      Rad4=Rad4*conv;
      rad3s=rad3s+sec[25];
      rad3p1=rad3p1+sec[18];
      rad3p2=rad3p2+sec[22];
      rad3d=rad3d+sec[27];
      AKLL=AKLL*conv;
      AKLM=AKLM*conv;
      ALMM=ALMM*conv;
      AM4=AM4*conv;
//c    ***********************************************
//c    actual cross sections
//c    ***********************************************
      T=0.;
      Call SecMean[Y][T];
//c    *******************************************************
//c    start integration
//c    *******************************************************
      call datetime[iyr][imon][iday][ih0][imi0][isec0][it0];
      printf (6,160) ih0,imi0,isec0;
//L160: format (2x,12hbegin time =,2x,i2,3h h ,i2,4h mn ,i2,2h s);
      t0=(ih0*3600+imi0*60+isec0+it0/1000.);
      iboucle=0;
      MINT=3;
      NROOT=0;
      T=0.;
      TOUT=ep0;
      toutc=0.;
      QM0=QM;
L100: dtout=(tout-toutc)/tc;
      DQM=fabs((QM-QM0)*ddq);
       if((dtout >= 1.) || [DQM >= 1.]) {
        dt=tout-toutc;
        call chgt[e][zp][ac][rho][as][bs][dt][ef][qm];
        e=ef;
        call popmean[Y];
        QM=Zp-(y1s+yl+ym+yn);
        Call SecMean[Y][T];
        toutc=tout;
        QM0=QM;
        call datetime[iyr][imon][iday][ih1][imi1][isec1][it1];
        t1=(ih1*3600+imi1*60+isec1+it1/1000.);
        ttot=t1-t0;
        iht=ttot/3600;
        imt=(ttot-3600*iht)/60;
        tst=(ttot-3600*iht-60*imt);
        printf (6,105) iht,imt,tst;
      }
//L105: format (2x23hrestart integration ...,14h(elapsed time=,
                2x,i2,3h h ,i2,4h mn ,f5.2,3h s));
      MSTATE=1;
      EPS=erel;
      EWT=erabs;
      LENW=1670000;
      LENIW=1310;
      X=T;
      nequ=0;
      YPMax=0.;
//c    ****************************************************
      call EQDIF[NEQ][T][Y][F][TOUT][MSTATE][NROOT][EPS][EWT][MINT][WORK][
                LENW][IWORK][LENIW][F];
//c    *************** IDID=1 as long as T <= TOUT **********
       if (MSTATE > 2){
      call message[mstate];
      goto L5000;
      }
      NSTEP=IWORK[3];
      iboucle=NSTEP;
      AVGORD=WORK[3];
      for(/*L101*/ IN=1; IN<=1283; IN++) {
       if (Y[IN] < 0.) {
      Y[IN]=0.;
      }
      YX[IN]=Y[IN]*100.;
L101: } //continue;
//c    **************** end of one step calculation ************
      X=T;
//c    ***********************************************
//c    charge state probabilities calculation
//c    ***********************************************
      call popmean[Y];
      QM=Zp-(y1s+yl+ym+yn);
//     charge states before autoionization
      for(/*L500*/ N=1; N<=61; N++) {
      PRF[N]=0.;
L500: } //continue;
      for(/*L600*/ N1=1; N1<=29; N1++) {
      for(/*L600*/ N4=1; N4<=33; N4++) {
      I=N1-1;
      N=N4-1;
      INM=100*N+I;
      NN=NUMPP[INM];
      Nel=I+N4;
      PRF[Nel]=PRF[Nel]+Y[NN];
L600: } //continue;
      Qin=0.;
      PTF=0.;
      for(/*L700*/ M=1; M<=61; M++) {
      RM=double[M]-1.;
      PTF=PTF+PRF[M];
      Qin=Qin+(Zp-RM)*PRF[M];
L700: } //continue;
//c    ***********************************************
//c    re-calculate cross sections for actual charge state
//c    ***********************************************
      Call SecMean[Y][T];
//c-------- K shell autoionization and K+L  populations --------
      Call Auger[Y][Zp][PR][YTOT];
      PT=0.;
      QF=0.;
      for(/*L9*/ M=1; M<=62; M++) {
      RM=double[M]-1.;
      PR[M]=100.*PR[M];
      PT=PT+PR[M];
      QF=QF+(Zp-RM)*PR[M]/100.;
L9:   } //continue;
      for(/*L50*/ K=1; K<=2; K++) {
      P1s[K]=0;
L50:  } //continue;
      for(/*L51*/ K=1; K<=2; K++) {
      for(/*L51*/ J=0; J<=2; J++) {
      for(/*L51*/ I=0; I<=6; I++) {
      L=100*I+10*J+K;
      N=NUM[L];
      P1s[K]=P1s[K]+Y[N];
L51:  } //continue;
      for(/*L52*/ J=1; J<=2; J++) {
      P2s[J]=0;
L52:  } //continue;
      for(/*L53*/ J=1; J<=2; J++) {
      for(/*L53*/ K=0; K<=2; K++) {
      for(/*L53*/ I=0; I<=6; I++) {
      L=100*I+10*J+K;
      N=NUM[L];
      P2s[J]=P2s[J]+Y[N];
L53:  } //continue;
      for(/*L54*/ I=1; I<=2; I++) {
      P2p[I]=0;
L54:  } //continue;
      for(/*L55*/ I=1; I<=2; I++) {
      for(/*L55*/ K=0; K<=2; K++) {
      for(/*L55*/ J=0; J<=2; J++) {
      L=100*I+10*J+K;
      N=NUM[L];
      P2p[I]=P2p[I]+Y[N];
L55:  } //continue;
//c    ********* prints probabilities at each step *****
//     write (10,125) T,(PR(I),I=1,10),e,PT
// 125 format(F10.3,F9.5,9(1x,F9.5),2x,2(1x,g10.5))
      printf (10,125) T,(PR[I],I=1,10),e;
      printf (11,25) T,(PR[I],I=11,20);
       if (ip >= 20) {
      printf (12,25) T,(PR[I],I=21,30);
      }
       if (ip >= 30) {
      printf (13,26) T,(PR[I],I=31,40);
      }
       if (ip >= 40) {
      printf (20,26) T,(PR[I],I=41,50);
      }
       if (ip >= 50) {
      printf (21,26) T,(PR[I],I=51,60);
      }
      tots=yx[1]+yx[2]+yx[4]+yx[10]+yx[3]+yx[5]+yx[11]+yx[6]+yx[12];
      printf (14,27) T,yx[1],yx[2],yx[4],yx[10],yx[3],yx[5],yx[11],
                yx[6],yx[12],tots;
// 
      printf (15,28) T,y1s,y2s,y2p,ym,yn,Qm,Qin,QF,PTF ;
//L25:  format(F10.3,F9.5,9(1x,F9.5));
//L125: format(F10.3,F9.5,9(1x,F9.5),3x,g12.5);
//L26:  format(F9.2,F9.5,9(1x,F9.5));
//L27:  format(F9.2,10(1x,F6.2));
//L28:  format(F9.2,3(1x,F7.5),6(1x,F8.5));
      T=TOUT;
      printf(6,29) T,iboucle,nequ;
//L29:  format(2x,19hachieved thickness=,f9.2,4x,17h(iteration # for ,
                10hthis step=,i5,1x,i5,1h));
      iboucle=0;
//c********* calculation goes on as long as T <= EPM *****
       if (TOUT < EPM) {
             if (TOUT < pas0){
//c*************** (pas0=20*ep0) ******************
              TOUT=TOUT+ep0;
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT <= 0.2) {
              TOUT=TOUT+ep0;
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 1.4999) {
               if (ep1 > 0.05) {
               TOUT=TOUT+0.05;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 1.9999) {
               if (ep1 > 0.1) {
               TOUT=TOUT+0.1;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 4.9999) {
               if (ep1 > 0.2) {
               TOUT=TOUT+0.2;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 14.999) {
               if (ep1 > 0.5) {
               TOUT=TOUT+0.5;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 19.99) {
               if (ep1 > 1.) {
               TOUT=TOUT+1.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 49.99) {
               if (ep1 > 2.) {
               TOUT=TOUT+2.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 149.99) {
               if (ep1 > 5.) {
               TOUT=TOUT+5.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 199.99) {
               if (ep1 > 10.) {
               TOUT=TOUT+10.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 499.99) {
               if (ep1 > 20.) {
               TOUT=TOUT+20.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 1499.9) {
               if (ep1 > 50.) {
               TOUT=TOUT+50.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 1999.) {
               if (ep1 > 100.) {
               TOUT=TOUT+100.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 4999.) {
               if (ep1 > 200.) {
               TOUT=TOUT+200.;
              }
 else {
               TOUT=TOUT+ep1;
              }
                if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 14999.) {
               if (ep1 > 500.) {
               TOUT=TOUT+500.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 19999.) {
               if (ep1 > 1000.) {
               TOUT=TOUT+1000.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else if (TOUT < 49999.) {
               if (ep1 > 2000.) {
               TOUT=TOUT+2000.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
 else {
               if (ep1 > 5000.) {
               TOUT=TOUT+5000.;
              }
 else {
               TOUT=TOUT+ep1;
              }
               if (TOUT > EPM) TOUT=EPM;
              iboucle=0.;
              goto L100;
            }
      }
//c ************* end of output files **************************
      call datetime[iyr][imon][iday][ih1][imi1][isec1][it1];
      printf (6,161) ih1,imi1,isec1;
      t1=(ih1*3600+imi1*60+isec1+it1/1000.);
      ttot=t1-t0;
      iht=ttot/3600;
      imt=(ttot-3600*iht)/60;
      tst=(ttot-3600*iht-60*imt);
      printf (6,162) iht,imt,tst;
//L161: format ("     end time = ",i2," h ",i2," mn ",i2," s");
//L162: format (" elapsed time = ",i2," h ",i2," mn ",f5.2," s");
       if(tc < epm) {
        dt=tout-toutc;
         if (dt > 0.) {
          ef=e+(bs-e)*(1.-exp(-0.001*tc*as));
          e=ef;
        }
      }
      printf (6,163) e;
//L163: format (35x,"final energy : ",g11.4,"(MeV/u)");
      printf (10,203);
      printf (10,210);
       if (ip >= 9) {
      printf (10,211)ip,ip-1,ip-2,ip-3,ip-4,ip-5,ip-6,ip-7,ip-8,ip-9;
      }
      printf (10,203);
      printf (10,160) ih0,imi0,isec0;
      printf (10,161) ih1,imi1,isec1;
      printf (10,163) e ;
      close[unit=10];
// 
      printf (11,203);
      printf (11,212);
       if (ip >= 19) {
        printf (11,211)ip-10,ip-11,ip-12,ip-13,ip-14,ip-15,ip-16,ip-17,ip-18,
           ip-19;
      }
      printf (11,203);
      printf (11,160) ih0,imi0,isec0;
      printf (11,161) ih1,imi1,isec1;
      printf (11,163) e;
      close[unit=11];

       if (ip >= 20) {
       printf (12,203);
       printf (12,213);
        if (ip >= 29) {
       printf (12,211)ip-20,ip-21,ip-22,ip-23,ip-24,ip-25,ip-26,ip-27,ip-28,
           ip-29;
       }
        printf (12,203);
        printf (12,160) ih0,imi0,isec0;
        printf (12,161) ih1,imi1,isec1;
        printf (12,163) e;
        close[unit=12];
      }
// 
       if (ip >= 30) {
        printf (13,203);
        printf (13,214);
         if (ip >= 39) {
        printf (13,215)ip-30,ip-31,ip-32,ip-33,ip-34,ip-35,ip-36,ip-37,ip-38,
            ip-39;
        }
        printf (13,203);
        printf (13,160) ih0,imi0,isec0;
        printf (13,161) ih1,imi1,isec1;
        printf (13,163) e;
        close[unit=13];
      }
// 
       if (ip >= 40) {
        printf (20,203);
        printf (20,218);
         if (ip >= 49) {
        printf (20,215)ip-40,ip-41,ip-42,ip-43,ip-44,ip-45,ip-46,ip-47,ip-48,
            ip-49;
        }
        printf (20,203);
        printf (20,160) ih0,imi0,isec0;
        printf (20,161) ih1,imi1,isec1;
        printf (20,163) e;
        close[unit=20];
      }

       if (ip >= 50) {
        printf (21,203);
        printf (21,219);
         if (ip >= 59) {
        printf (21,215)ip-50,ip-51,ip-52,ip-53,ip-54,ip-55,ip-56,ip-57,ip-58,
           ip-59;
        }
        printf (21,203);
        printf (21,160) ih0,imi0,isec0;
        printf (21,161) ih1,imi1,isec1;
        printf (21,163) e;
        close[unit=21];
      }
// 
      printf (14,203);
      printf (14,216);
      printf (14,203);
      printf (14,160) ih0,imi0,isec0;
      printf (14,161) ih1,imi1,isec1;
      printf (14,163) e;
// 
      for(/*L*/ n3=327; n3<=1284; n3++) {
       if (Y[n3] >= 0.01) {
      printf(14,*) 'Code=',Ico3[n3], 'pop=',Y[n3];
      }
      } for(/*L*/ ) {
      close[unit=14];
// 
      printf (15,203);
      printf (15,217);
      printf (15,203);
      printf (15,160) ih0,imi0,isec0;
      printf (15,161) ih1,imi1,isec1;
      printf (15,163) e;
      close[unit=15];

      close[unit=16];
      printf(6,*) 'output double /*data*/ in files:';
      printf(6,*) '   0 to  9 e- charge states    in ETA09.PRN';
      printf(6,*) '  10 to 19 e- charge states    in ETA1019.PRN';
       if (ip >= 20) {
      printf(6,*) '  20 to 29 e- charge states    in ETA2029.PRN';
      }
       if (ip >= 30) {
      printf(6,*) '  30 to 39 e- charge states    in ETA3039.PRN';
      }
       if (ip >= 40) {
      printf(6,*) '  40 to 49 e- charge states    in ETA4049.PRN';
      }
       if (ip >= 50) {
      printf(6,*) '  50 to 59 e- charge states    in ETA5059.PRN';
      }
      printf(6,*) ' bare,1s,2s,2p,1sý,1s2s,1s2p,1sý2s,1sý2p ions ';
      printf(6,*) ' and sum of these             in ETAPIED.PRN';
      printf(6,*) ' mean 1s,2s,2p,3s,3p and 3d populations';
      printf(6,*) '                              in POPMEAN.PRN';
      printf(6,*) '  ';
      printf(6,*)'WARNING! Next calculation will overwrite these files';
      printf(6,*)'Consider saving or renaming these results ! ';
L5000: } //continue;
      }
//c    ***************************************************
//c    *********** END of etacha *************************
//c    ***************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void chgt(e,zp,ac,rho,as,bs,tc,ef,qm) {
      common/seceff/sec[34],cor[48],secs[48];
      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn;
      common/aug/AKLL,AKLM,ALMM,AM4;
      common/popi/n02s,n02p,n03s,n03p,n03d,nl0,nm0;
      ef=e+(bs-e)*(1.-exp(-0.001*tc*as));
      e=ef;
      call CSEC[ef];
//L220: format(1x,f8.3,1x,f7.3,6(1x,G10.3));
//L203: format(89('*'));
//L308: format(2x,"Reference ('hydrogenic') cross sections in 10E-20 cm²:",
              /,"  CAPTURE (MEC+REC) TO:  ",
              " 1s =",G10.3,"  2s =",G10.3,"  2p =",G10.3,/,25x,
              " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
              "  (n=4) =",G10.3,/,
                "  IONIZATION OF:",9x,
              " 1s =",G10.3,"  2s =",G10.3,"  2p =",G10.3,/,25x,
              " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
              "  (n=4) =",G10.3,/,
                2x,"1s EXCITATION TO:",
              23x,            "2s =",G10.3,"  2p =",G10.3,/,25x,
              " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
                "  (n=4) =",G10.3,/,
                2x,"2s EXCITATION TO:",6x,
                " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
                "  (n=4) =",G10.3,/,
              2x,"2p EXCITATION TO:",6x,
              " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
              "  (n=4) =",G10.3,/,
                2x,"EXCITATION To n=4 from:",
                " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,/,
                2x,"INTRASHELL EXCITATION:",
                "  2s to 2p =",G10.3,2x,"3s to 3p =",G10.3);
//L309: format (\,"  Binding Energy correction factors =",6(2x,F8.4));
//L310: format (\,"  Effective charges =",6(2x,F8.4)) ;
//L311: format ("  Binding energies  =",6(2x,F8.1)) ;
      printf (16,*) 'QM=',QM,'E=',E;
      printf (16,308) (sec[I],I=1,33);
      printf (16,309) tetak,tetal1,tetal2,tetam1,tetam2,tetan;
      printf (16,310) zk,zl1,zl2,zm1,zm2,zn;
      printf (16,311) Bk,Bl1,Bl2,Bm1,Bm2,Bn;
      printf (16,203);
      dum=qm;
      e=ef;
      conv=6.023E-3/AC;
      printf(6,*)'Radiative and Auger yields';
      VC=pow((1-pow((1+E/931.5),(-2))),(0.5));
      RAD=3.461E-6*pow((Zk*tetak),4)*AC/(RHO*VC)/2.;
      RAD3s=.034891E-6*pow((Zk*tetak),4)*AC/(RHO*VC)/6.;
      RAD3p=1.0301E-6*pow((Zk*tetak),4)*AC/(RHO*VC)/2.;
      RAD3d=.3544E-6*pow((Zk*tetak),4)*AC/(RHO*VC)/6.;
      RAD4=0.0454E-6*pow((Zk*tetak),4)*AC/(RHO*VC)/18.;
       if (n02p >= 1) {
      AKLL=rad*(1./ok-1.)*n02p/(nl0*(nl0-1));
      }
 else {
      AKLL=7.6E-2*AC/(RHO*VC);
      }
       if (nm0 >= 1) {
      AKLM=rad3p*0.882*(1./ok-1.)*n02p/(nl0*nm0);
      }
 else {
      AKLM=0.3*AC/(RHO*VC);
      }
       if (nm0 >= 2) {
      ALMM=(rad3s*n03s+0.118*rad3p*n03p+rad3d*n03d)*(1./ol-1.)
           /(nm0*(nm0-1));
      }
 else {
      ALMM=0.05*AC/(RHO*VC);
      }
      AM4=3.*ALMM;
      for(/*L1*/ I=1; I<=34; I++) {
      Sec[I]=Sec[I]*conv;
L1:   } //continue;
      Sec[1]=sec[1]/2.;
      sec[2]=sec[2]/2.;
      sec[3]=sec[3]/6.;
      Sec[4]=sec[4]/2.;
      sec[5]=sec[5]/6.;
      sec[6]=sec[6]/10.;
      sec[7]=sec[7]/32.;
      sec[15]=sec[15]/2.;
      sec[16]=sec[16]/6.;
      sec[17]=sec[17]/2.;
      sec[18]=sec[18]/6.;
      sec[19]=sec[19]/10.;
      sec[20]=sec[20]/32.;
      sec[21]=sec[21]/2.;
      sec[22]=sec[22]/6.;
      sec[23]=sec[23]/10.;
      sec[24]=sec[24]/32.;
      sec[25]=sec[25]/2.;
      sec[26]=sec[26]/6.;
      sec[27]=sec[27]/10.;
      sec[28]=sec[28]/32.;
      sec[29]=sec[29]/32.;
      sec[30]=sec[30]/32.;
      sec[31]=sec[31]/32.;
      sec[32]=sec[32]/6.;
      sec[33]=sec[33]/6.;
      sec[34]=sec[34]/10.;
      Rad=rad*conv;
      Rad3s=Rad3s*conv;
      Rad3p1=Rad3p*0.882*conv;
      Rad3p2=Rad3p*0.118*conv;
      Rad3d=Rad3d*conv;
      Rad4=Rad4*conv;
      rad3s=rad3s+sec[25];
      rad3p1=rad3p1+sec[18];
      rad3p2=rad3p2+sec[22];
      rad3d=rad3d+sec[27];
      AKLL=AKLL*conv;
      AKLM=AKLM*conv;
      ALMM=ALMM*conv;
      AM4=AM4*conv;
      return;
      }
//c    ***************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void popmean(u) {
      double u[1283];
      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
      common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m;
//c    mean number of electrons in each subshell
      y3s=U[65]+2.*U[66];
      y3p=U[68]+2.*U[69]+3.*U[70]+4.*U[71]+5.*U[72]+6.*U[73];
      y3d=U[75]+2.*U[76]+3.*U[77]+4.*U[78]+5.*U[79]+6.*U[80];
      y3d=y3d+7.*U[81]+8.*U[82]+9.*U[83]+10.*U[84];
      ym=y3s+y3p+y3d;
      yn=0.;
      ynn=0.;
      for(/*L10*/ n=1; n<=32; n++) {
      rn=double[n];
      yn=yn+rn*U[294+n];
       if (n >= 2) {
      ynn=ynn+rn*(rn-1.)*U[294+n];
      }
L10:  } //continue;
// -------- see F.for (ALMM) ------------------------------
      ymp=0.;
      ymm=0.;
      for(/*L1*/ I=0; I<=2; I++) {
      for(/*L1*/ J=0; J<=6; J++) {
      for(/*L1*/ K=0; K<=10; K++) {
      N=I+J+K;
      ymp=ymp+U[64+I]*U[67+J]*U[74+K]*N;
       if (N >= 2) {
      ymm=ymm+U[64+I]*U[67+J]*U[74+K]*N*(N-1);
      }
L1:   } //continue;
      ym1m=0.;
      for(/*L2*/ I=0; I<=2; I++) {
      for(/*L2*/ J=0; J<=6; J++) {
      N=I+J-1;
       if (N >= 2) {
      ym1m=ym1m+U[64+I]*U[67+J]*(N-1);
      }
L2:   } //continue;
      ym2m=U[76]+2.*U[77]+3.*U[78]+4.*U[79]+5.*U[80];
      ym2m=ym2m+6.*U[81]+7.*U[82]+8.*U[83]+9.*U[84];
//------------------------------------------------------
//.....mean values for 1s, 2s, 2p ......
      y1s=0.;
      y2s=0.;
      y2p=0.;
      for(/*L3*/ N=1; N<=63; N++) {
      I=II[N];
      J=JJ[I][N];
      K=KK[I][J][N];
      y1s=y1s+K*U[N];
      y2s=y2s+J*U[N];
      y2p=y2p+I*U[N];
L3:   } //continue;
      yl=y2s+y2p;
//.....mean values for 1s², 2s², 2p² (in case is needed) ......
      ykm=0.;
      yl1m=0.;
      yl2m=0.;
      for(/*L4*/ N=1; N<=63; N++) {
      I=II[N];
      J=JJ[I][N];
      K=KK[I][J][N];
       if(K >= 2) {
      ykm=ykm+U[N];
      }
       if(J >= 2) {
      yl1m=yl1m+U[N];
      }
       if(I >= 2) {
      yl2m=yl2m+(I-1)*U[N];
      }
L4:   } //continue;
      return;
      }
//c*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function omk(z) {
      double ok[92];
      double /*data*/ (ok[i],i=1,92) / 0.,0.,0.001,0.0012,
                .170E-02,.280E-02,.520E-02,.830E-02,.130E-01,.180E-01,.230E-01,
                .300E-01,.390E-01,.500E-01,.630E-01,.780E-01,.970E-01,.118,
                .140   ,.163        ,.188    ,.214    ,.243    ,.275    ,.308,
                .340   ,.373        ,.406    ,.440    ,.474    ,.507    ,.535,
                .562   ,.589        ,.618    ,.643    ,.667    ,.690    ,.710,
                .730   ,.747        ,.765    ,.780    ,.794    ,.808    ,.820,
                .831   ,.843        ,.853    ,.862    ,.870    ,.877    ,.884,
                .891   ,.897        ,.902    ,.907    ,.912    ,.917    ,.921,
                .925   ,.929        ,.932    ,.935    ,.938    ,.941    ,.944,
                .947   ,.949        ,.951    ,.953    ,.955    ,.957    ,.958,
                .959   ,.961        ,.962    ,.963    ,.964    ,.965    ,.966,
                .967   ,.968        ,.968    ,.969    ,.969    ,.970    ,.970,
                .971   ,.971        ,.972        ,.972 /;
      i=Z;
      omk=ok[i];
      return;
      }
//c*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function oml(z) {
      double ol[92];
      double /*data*/ (ol[i],i=1,92) / 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                .120E-02,.750E-03,.380E-03,.310E-03,.260E-03,.240E-03,.220E-03,
            .270E-03,.330E-03,.840E-03,.150E-02,.260E-02,.370E-02,.500E-02,
            .630E-02,.770E-02,.930E-02,.110E-01,.120E-01,.130E-01,.150E-01,
            .160E-01,.180E-01,.200E-01,.220E-01,.240E-01,.260E-01,.280E-01,
            .310E-01,.340E-01,.370E-01,.400E-01,.430E-01,.460E-01,.490E-01,
            .520E-01,.560E-01,.600E-01,.640E-01,.690E-01,.740E-01,.790E-01,
            .850E-01,.910E-01,.970E-01,.104    ,.111    ,.118    ,.125    ,
                .132  ,.139  ,.147        ,.155    ,.164    ,.174    ,.182    ,
            .192    ,.201    ,.210    ,.220    ,.231    ,.243    ,.255    ,
            .268    ,.281    ,.294    ,.306    ,.320    ,.333    ,.347    ,
            .360    ,.373    ,.386    ,.399    ,.411    ,.424    ,.437    ,
                .450  ,.463  ,.476        ,.489 /;
      i=Z;
      oml=ol[i];
      return;
      }
//c*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void message(idid) {
      printf (6,*)'problem integrating with QDIFF';
       if (idid == 3) {
      printf (6,*)' more than 1000 iterations have been atempted ...';
      printf (6,*)' try reducing maximum step size';
      }
       if (idid == 4) {
      printf (6,*)' you are probably asking too much accuracy ...';
      printf (6,*)' try again with larger uncertainties';
      }
       if (idid > 4) {
      printf (6,*)' for some reason, the problem is very stiff and',
                ' cannot be solved with the present integration routine';
      }
      return;
      }
//c*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void datetime(iyr,imon,iday,ihour,imin,isec,imil) {
      int Datime[8];
      char*12 Dum[3];
      Call date_and_time[Dum[1]][Dum[2]][Dum[3]][Datime];
      iyr=datime[1];
      imon=datime[2];
      iday=datime[3];
      ihour=datime[5];
      imin=datime[6];
      isec=datime[7];
      imil=datime[8];
      return;
      }