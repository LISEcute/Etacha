//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void snl(esp,esp3,epd3) {
//c*********************************************************************
//c        Calcule les sections efficaces d'excitation ns-np  dans
//c    l'approximation PWBA avec:
//c        -facteurs de forme analytiques utilisant la forme recurrente
//c        etablie a partir des calculs n=2 a n=6
//c        -des facteurs d'ecrantage et d'antiecrantage Thomas-Fermi
//c    Les facteurs de forme sont calcules en fonction de la variable
//c    reduite k =  q/Zp (q en u.a.)
//c
//c    Les resultats sont donnes ici en unites de 10-20 cm2
//c*********************************************************************
      double secnl[3];
      int LIMIT,LENW;
      PARAMETER[LIMIT=100][LENW=LIMIT*4];
      double BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK[LENW];
      int INF,NEVAL,IER,LAST,IWORK[LIMIT];
      EXTERNAL FNL;
      char erfi*12;
      common/don/zp,E,zt;
      common /nsp/ n;
      erfi='messerr.fil';
//c******integration***********
      bet=pow((1-pow((1+E/931.5),(-2))),(0.5));
      for(/*L30*/ n=1; n<=3; n++) {
      BOUND = 1.E-5/zp;
      INF = 1;
      EPSABS = 0.0;
      EPSREL = 1.E-4;
//*************************************************
      CALL INTG[FNL][BOUND][INF][EPSABS][EPSREL][RESULT][ABSERR][NEVAL][
                  IER][LIMIT][LENW][LAST][IWORK][WORK];
      sig1=RESULT;
       if (ier > 0) {
      fopen(unit=9,file=erfi,status='unknown');
      printf(9,*)' n= ',n,'  ier= ',ier;
      printf(9,*)' result= ',result,' abserr= ',abserr;
      close[unit=9];
      }
//*************************************************
      secnl[n]=sig1*(1.874)*zt*zt/(zp*zp*bet*bet);
L30:  } //continue;
      esp=secnl[1];
      esp3=secnl[2];
      epd3=secnl[3];
      return;
      }
//c*************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function fnl[xx] {
      double*8 q,corr,scf,stf,atf,fac,bp;
      common/don/zp,E,zt;
      common /nsp/ n;
      scf=dble[1.13*pow(zt,(1./3.]));
      q=dble[xx*zp];
      stf=(scf)/(q);
      atf=(q)/(scf);
      corr=dble[1./(1.+pow(pow(stf,2.]),2.));
      corr=corr+dble[(1.-(1./(1.+pow(pow(atf,2.]),2.)))/zt);
//c*********************************************************************
      bp=dble[xx];
       if (n == 1) {
// 2s-2p
      fac=1.-pow(bp,2.);
      fac=18.*pow(fac,2.)/(bp*pow((1.+pow(bp,2.)),8.));
      }
 else if (n == 2) {
// 3s-3p
      fac=128.-1392.*pow(bp,2.)+3888.*pow(bp,4.)-2187.*pow(bp,6.);
      fac=110592.*pow(fac,2.)/(bp*pow((4.+9.*pow(bp,2.)),12.));
      }
 else if (n == 3) {
// 3p-3d
      fac=256.-3840.*pow(bp,2.)+30816.*pow(bp,4.)-81648.*pow(bp,6.);
      fac=fac+76545.*pow(bp,8.);
      fac=2949120.*fac/(bp*pow((4.+9.*pow(bp,2.)),12.));
      }
      fnl=fac*corr;
      return;
      }
//c*********************************************************************
