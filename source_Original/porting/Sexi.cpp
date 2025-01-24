//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void sexi(e2s,e2p,e3s,e3p,e3d,es4,es5) {
//c         Version du 02/03/94 simplifiee pour eta3
//c           dernieres corrections : 31/08/2000 (COA)
//c    version avec makefile utilisant la routine d'integration INTG
//c*********************************************************************
//c        Calcule les sections efficaces d'excitation 1s-2s, 1s-2p et
//c    1s-n dans l'approximation PWBA avec:
//c        -facteurs de forme analytiques d'Anholt
//c        -facteurs d'ecrantage et d'antiecrantage dans le modŠle
//c         de Thomas-Fermi
//c    Les facteurs de forme sont calcules en fonction de la variable
//c    reduite k =  q/Zp (q en u.a.)
//c
//c    Les resultats sont donnes en unites de 10-20 cmý
//c*********************************************************************
      int LIMIT,LENW;
      PARAMETER[LIMIT=100][LENW=LIMIT*4];
      double BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK[LENW];
      double*8 sig1,rn,qmin,scf,coa;
      int INF,NEVAL,IER,LAST,IWORK[LIMIT];
      EXTERNAL FEX;
//**********************************************************************
      char erfi*12;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common /don/zp,E,zt;
      common /nn/ rn,coa,scf,ll;
      common/corr/ibin;
      erfi='messerr.fil';
      bet=pow((1.-pow((1.+E/931.5),(-2.))),(0.5));
      vu=137.036*bet;
      fs1=1.;
       if (zt == 1.) fs1=1.23;
       if (zt == 2.) fs1=1.526;
       if (zt == 6.) fs1=0.78;
       if (zt == 7.) fs1=0.85;
       if (zt == 10.) fs1=1.04;
       if (zt == 13.) fs1=0.58;
       if (zt == 14.) fs1=0.59;
       if (zt == 18.) fs1=0.68;
       if (zt == 29.) fs1=0.672;
       if (zt == 36.) fs1=0.61;
       if (zt == 54.) fs1=0.535;
      scf=dble[fs1*1.13*pow(zt,(1./3.]));
      nmin=2;
      nmax=5;
      som1=0.;
      so51=0.;
//c****************************
//c******integration***********
      for(/*L20*/ n=nmin; n<=nmax; n++) {
       if (n == 2) {
          ll=1;
      }
       if (n == 3) {
          ll=4;
      }
       if (n > 3) {
          ll=7;
      }
L21:  rn=dble[n];
//c**********************************
//c correction pour effet energie liaison:
//c**********************************
       if(ibin == 1) {
      y=2.*vu/(zk*tetak);
      g=(1.+5.*y+7.14*pow(y,2.)+4.27*pow(y,3.)+0.947*pow(y,4.))/pow((1.+y),5.);
      epsi=pow((1+zt*g/(zk*tetak)),2.);
      }
 else {
      epsi=1.;
      }
//c****************************
//     pmin=deltaE/2v => qmin=deltaE/2v/Zp

      qmin=dble[epsi*(zk*tetak/(vu*2.])*(1.-1./(rn*rn)));
//     calcul du coefficient d'antiscreening (08/2000)
      vs=1.75*epsi*zk*tetak*qmin;
       if (vu > vs) {
      coas=1.-pow((vs/vu),2);
      }
 else {
      coas=0.;
      }
      coa=dble[coas];
      BOUND = qmin;
      INF = 1;
      EPSABS = 0.0;
      EPSREL = 1.E-4;
//*************************************************
      CALL INTG[FEX][BOUND][INF][EPSABS][EPSREL][RESULT][ABSERR][NEVAL][
                  IER][LIMIT][LENW][LAST][IWORK][WORK];
      sig1=RESULT;
       if (ier > 0) {
      fopen(unit=9,file=erfi,status='unknown');
      printf(9,*)' n= ',n,'  ier= ',ier;
      printf(9,*)' result= ',result,' abserr= ',abserr;
      close[unit=9];
      }
//*************************************************
      sec1=sig1*(1.874)*zt*zt/(zk*zk*bet*bet);
       if(ll == 3 || ll == 7) {
      som1=som1+sec1;
      }
       if (ll < 3) {
           if(ll == 1) e2s=sec1;
           if(ll == 2) e2p=sec1;
      ll=ll+1;
      goto L21;
      }
       if (ll > 3 && ll < 7) {
           if(ll == 4) e3s=sec1;
           if(ll == 5) e3p=sec1;
           if(ll == 6) e3d=sec1;
      ll=ll+1;
      goto L21;
      }
       if (n == 4) {
      es4=sec1;
      }
       if (n > 4) {
      so51=so51+sec1;
      }
L20:  } //continue;
//     so51=so51+4.52*sec1
      es5=so51;
      return;
      }
//c*************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function FEX[x] {
      double*8 q,scf,stf,coa,corr,fact,bp,rn;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common /don/zp,E,zt;
      common /nn/ rn,coa,scf,ll;
      q=dble[x*zk];
      stf=(scf)/(q);
      atf=(q)/(scf);
      corr=dble[1./(1.+pow(pow(stf,2.]),2.));
      corr=corr+coa*dble[(1.-(1./(1.+pow(pow(atf,2.]),2.)))/zt);
//c*********************************************************************
      bp=dble[x];
       if (ll == 1) {
         fact=64.*bp/pow((bp*bp+9./4.),6.);
      }
 else if (ll == 2) {
         fact=144./(bp*pow((bp*bp+9./4.),6.));
      }
 else if (ll == 4) {
         fact=1119744.*bp*pow((16.+27.*bp*bp),2.)/pow((9.*bp*bp+16.),8.);
      }
 else if (ll == 5) {
         fact=2985984.*pow((16.+27.*bp*bp),2.)/(bp*pow((9.*bp*bp+16.),8.));
      }
 else if (ll == 6) {
         fact=573308928.*bp/pow((9.*bp*bp+16.),8.);
      }
 else {
         fact=(512./(3.*rn*rn*rn))*(3.*bp*bp+1.-1./(rn*rn));
         fact=fact*pow((bp*bp+pow((1.-1./rn),2.)),(rn-3.));
         fact=fact/(bp*pow((bp*bp+pow((1.+1./rn),2.)),(rn+3.)));
      }
      FEX=fact*corr;
      return;
      }
//c*********************************************************************
