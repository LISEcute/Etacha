//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void sex2(s3s,s3p,s3d,p3s,p3p,p3d,e2s4,e2p4,e2s5,e2p5) {
//c         Version du 21/06/2000 simplifiée pour eta4
//c           dernieres corrections : 31/08/2000 (COA et certains facteurs)
//c    version utilisant la routine d'intégration INTG
//c*********************************************************************
//c        Calcule les sections efficaces d'excitation 2l-3l' et 2l-4l'
//c        dans l'approximation PWBA avec facteurs d'écrantage et
//c        d'antiécrantage dans le modèle de Thomas-Fermi
//c    Les facteurs de forme sont calculés en fonction de la variable
//c    réduite k = q/Zp (q en u.a.)
//c
//c    Les résultats sont donnés en unités de 10-20 cm2
//c*********************************************************************
      int LIMIT,LENW;
      PARAMETER[LIMIT=100][LENW=LIMIT*4];
      double BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK[LENW];
      double*8 sig1,rn,qmin,scf,coa;
      int INF,NEVAL,IER,LAST,IWORK[LIMIT];
      EXTERNAL FEX2;
//**********************************************************************
      double sec[10];
      char erfi*12;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common /don/zp,E,zt;
      common/nn/rn,coa,scf,ll;
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
//c****************************
//c******integration***********
      for(/*L20*/ n=1; n<=10; n++) {
      rn=dble[n];
//c**********************************
//c correction pour effet energie liaison:
//c**********************************
       if (ibin == 1) {
      y=4.*vu/(zl2*tetal2);
      g=1.+9.*y+34.7*pow(y,2.)+81.2*pow(y,3.)+112.*pow(y,4.);
      g=g+93.5*pow(y,5.)+46.6*pow(y,6.)+12.9*pow(y,7.)+1.549*pow(y,8.);
      g=g/pow((1.+y),9.);
      epsi=pow((1+zt*g/(zl2*tetal2)),2.);
      }
 else {
      epsi=1.;
      }
//c****************************
       if(n <= 6) {
//      2 -> 3
      qmin=dble[epsi*(zl2*tetal2/(vu*2.])*(5./36.));
      }
 else if(n <= 8) {
//      2 -> 4
      qmin=dble[epsi*(zl2*tetal2/(vu*2.])*(3./16.));
      }
 else {
//      2 -> 5
      qmin=dble[epsi*(zl2*tetal2/(vu*2.])*(21./100.));
      }
//     calcul du coefficient d'antiscreening (08/2000)
      vs=1.75*epsi*zl2*tetal2*qmin;
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
      CALL INTG[FEX2][BOUND][INF][EPSABS][EPSREL][RESULT][ABSERR][NEVAL][
                  IER][LIMIT][LENW][LAST][IWORK][WORK];
      sig1=RESULT;
       if (ier > 0) {
      fopen(unit=9,file=erfi,status='unknown');
      printf(9,*)' n= ',n,'  ier= ',ier;
      printf(9,*)' result= ',result,' abserr= ',abserr;
      close[unit=9];
      }
//*************************************************
      sec[n]=sig1*(1.874)*zt*zt/(zl2*zl2*bet*bet);
L20:  } //continue;
      s3s=sec[1];
      s3p=sec[2];
      s3d=sec[3];
      p3s=sec[4];
      p3p=sec[5];
      p3d=sec[6];
      e2s4=sec[7];
      e2p4=sec[8];
      e2s5=sec[9];
      e2p5=sec[10];
      return;
      }
//c*************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function FEX2[x] {
      double*8 q,scf,stf,atf,corr,fact,bp,rn,coa;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common /don/zp,E,zt;
      common/nn/rn,coa,scf,ll;
      q=dble[x*zl2];
      stf=(scf)/(q);
      atf=(q)/(scf);
      corr=dble[1./(1.+pow(pow(stf,2.]),2.));
      corr=corr+coa*((1.-(1./pow((1.+pow(atf,2.)),2.)))/zt);
//c*********************************************************************
      bp=dble[x];
       if (rn == 1) {
//c 2s->3s
      fact=pow(2.,18.)*pow(3.,7.)*bp/pow((25.+36.*pow(bp,2.)),10.);
      fact=fact*pow((2875.-6984.*pow(bp,2.)+3888.*pow(bp,4.)),2.);
      }
 else if (rn == 2) {
//c 2s->3p
      fact=2.*pow(27648.,2.)/(bp*pow((25.+36.*pow(bp,2.)),10.));
      fact=fact*pow((-625.+5040.*pow(bp,2.)-3888.*pow(bp,4.)),2.);
      }
 else if (rn == 3) {
//c 2s->3d
      fact=2.*pow(65536.,2.)*pow(3.,7.)*bp/pow((25.+36.*pow(bp,2.)),10.);
      fact=fact*pow((-25.+18.*pow(bp,2.)),2.);
      }
 else if (rn == 4) {
//c 2p->3s
//c corrigé le 19/08/2009 : 3.**7 remplacé par 3.**6
      fact=pow(2.,16.)*pow(3.,6.)/(bp*pow((25.+36.*pow(bp,2.)),10.));
      fact=fact*pow((625.-14040.*pow(bp,2.)+11664.*pow(bp,4.)),2.);
      }
 else if (rn == 5) {
//c 2p->3p
      fact=pow(2.,25.)*pow(3.,9.)*bp/pow((25.+36.*pow(bp,2.)),10.);
      fact=fact*(6875.-23400.*pow(bp,2.)+34992.*pow(bp,4.));
      }
 else if (rn == 6) {
//c 2p->3d
      fact=pow(2.,24.)*pow(3.,6.)*pow(5.,2.)/(bp*pow((25.+36.*pow(bp,2.)),10.));
      fact=fact*(3125.-5400.*pow(bp,2.)+27216.*pow(bp,4.));
      }
 else if (rn == 7) {
//c 2s->4
      fact=405.+5120.*pow(bp,2.)+98816.*pow(bp,4.);
      fact=fact-163840.*pow(bp,6.)+65536.*pow(bp,8.);
      fact=fact*pow(2.,20.)/(bp*pow((9.+16.*pow(bp,2.)),9.));
      }
 else if (rn == 8) {
//c 2p->4
      fact=4194304.*(123. + 3728.*pow(bp,2 )- 3840.*pow(bp,4 )+
                12288.*pow(bp,6))/(bp*pow((9.+16.*pow(bp,2)),9));
      }
 else if (rn == 9) {
//c 2s->5
      fact=(20480000000.*(1555848.+29336825.*pow(bp,2)+376570000.*pow(bp,4)+ 
                3689500000.*pow(bp,6)-6300000000.*pow(bp,8)+2500000000.*pow(bp,10)))/
            [bp*(49.+100.*pow(pow(bp,2]),10));
      }
 else {
//c 2p->5
      fact=(128000000000.*(1643158523433.+74683227547040.*pow(bp,2)+ 
                993723878526400.*pow(bp,4)+4820995175360000.*pow(bp,6 )+ 
                14083087256000000.*pow(bp,8)+45107440000000000.*pow(bp,10)+ 
                81423200000000000.*pow(bp,12)-79200000000000000.*pow(bp,14)+ 
                90000000000000000.*pow(bp,16)))/(bp*pow((49.+100.*pow(bp,2)),14));
      }
      FEX2=fact*corr;
      return;
      }
//c*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void sex3(e3s4,e3p4,e3d4,e3s5,e3p5,e3d5,e4s5,e4p5,e45,e456) {
//c         Version du 25/01/95 simplifiée pour eta4
//c           dernieres corrections : 31/08/2000 (COA)
//c    version utilisant la routine d'intégration INTG
//c*********************************************************************
//c      Calcule les sections efficaces d'excitation 3l-4l'
//c        dans l'approximation PWBA avec facteurs d'ecrantage et
//c        d'antiecrantage dans le modŠle de Thomas-Fermi
//c    Les facteurs de forme sont calcules en fonction de la variable
//c    reduite k = q/Zp (q en u.a.)
//c
//c    Les resultats sont donnes en unites de 10-20 cm2
//c*********************************************************************
      int LIMIT,LENW;
      PARAMETER[LIMIT=100][LENW=LIMIT*4];
      double BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK[LENW];
      double*8 sig1,rn,qmin,scf,coa;
      int INF,NEVAL,IER,LAST,IWORK[LIMIT];
      EXTERNAL FEX3;
//**********************************************************************
      double sec[10];
      char erfi*12;
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common /don/zp,E,zt;
      common/nn/rn,coa,scf,ll;
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
//c******integration***********
      for(/*L20*/ n=1; n<=10; n++) {
      rn=dble[n];
//c**********************************
//c correction pour effet energie liaison:
//c**********************************
       if (ibin == 1) {
      y=18.*vu/(zm2*tetam2);
      g=(1.+5.*y+7.14*pow(y,2.)+4.27*pow(y,3.)+0.947*pow(y,4.))/pow((1.+y),5.);
      epsi=pow((1+zt*g/(zm2*tetam2)),2.);
      }
 else {
      epsi=1.;
      }
//c****************************
       if (n <= 3) {
//     3 vers 4 
      qmin=dble[epsi*(zm2*tetam2/(vu*2.])*(7./144.));
      vs=1.75*epsi*zm2*tetam2*qmin ;
      }
 else if (n <= 6) {
//     3 vers 5 
      qmin=dble[epsi*(zm2*tetam2/(vu*2.])*(16./225.));
      vs=1.75*epsi*zm2*tetam2*qmin ;
      }
 else if (n <= 9) {
//     4 vers 5 
      qmin=dble[(zn*tetan/(vu*2.])*(9./400.));
      vs=1.75*epsi*zn*tetan*qmin ;
      }
 else {
//     4 vers 6 
      qmin=dble[(zn*tetan/(vu*2.])*(5./144.));
      vs=1.75*epsi*zn*tetan*qmin;
      }
//     calcul du coefficient d'antiscreening (08/2000)
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
      CALL INTG[FEX3][BOUND][INF][EPSABS][EPSREL][RESULT][ABSERR][NEVAL][
                  IER][LIMIT][LENW][LAST][IWORK][WORK];
      sig1=RESULT;
       if (ier > 0) {
      fopen(unit=9,file=erfi,status='unknown');
      printf(9,*)' n= ',n,'  ier= ',ier;
      printf(9,*)' result= ',result,' abserr= ',abserr;
      close [unit=9];
      }
//*************************************************
       if (n >= 7) {
      sec[n]=sig1*(1.874)*zt*zt/(zn*zn*bet*bet);
      }
 else {
      sec[n]=sig1*(1.874)*zt*zt/(zm2*zm2*bet*bet);
      }
L20:  } //continue;
      e3s4=sec[1];
      e3p4=sec[2];
      e3d4=sec[3];
//     (facteur theorique en 1/n**3 = 3.049358 pour n=5)
//     les facteurs sont dans CSEC(EF) et valent 2.5 (13/09/01)
      e3s5=sec[4];
      e3p5=sec[5];
      e3d5=sec[6];
      e4s5=sec[7];
      e4p5=sec[8];
      e45=sec[9];
//     on suppose que la loi en 1/n**3 n'est pas etablie des n=6
//     (facteur theorique en 1/n**3 = 3.5412910805 pour n=6)
//     e456=sec(9)+2.5*sec(10)
      e456=sec[9]+sec[10];
      return;
      }
//c*************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function FEX3[x] {
      double*8 q,scf,stf,atf,corr,fact,bp,rn,coa;
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/don/zp,E,zt;
      common/nn/rn,coa,scf,ll;
       if (rn <= 6.) {
      q=dble[x*zm2];
      }
 else {
      q=dble[x*zn];
      }
      stf=(scf)/(q);
      atf=(q)/(scf);
      corr=dble[1./(1.+pow(pow(stf,2.]),2.));
      corr=corr+coa*((1.-(1./pow((1.+pow(atf,2.)),2.)))/zt);
//c*********************************************************************
      bp=dble[x];
       if (rn == 1) {
//c 3s->4
      fact=4250070125.+110475272992.*pow(bp,2.);
      fact=fact+251535456000.*pow(bp,4.);
      fact=fact-8858287816704.*pow(bp,6.);
      fact=fact+30668348915712.*pow(bp,8.);
      fact=fact-29926726041600.*pow(bp,10.);
      fact=fact+8916100448256.*pow(bp,12.);
      fact=fact*18345885696./(bp*pow((49.+144.*pow(bp,2.)),11.));
      }
 else if (rn == 2) {
//c 3p->4
      fact=401065441.+5869103184.*pow(bp,2)-57396708864.*pow(bp,4);
      fact=fact+313343188992.*pow(bp,6)-651422269440.*pow(bp,8);
      fact=fact+557256278016.*pow(bp,10);
      fact=fact*pow(2.,30)*pow(3.,5)/(bp*pow((49.+144.*pow(bp,2)),11));
      }
 else if (rn == 3) {
//c 3d->4
      fact=11008585.-92829520.*pow(bp,2)+552524544.*pow(bp,4);
      fact=fact+417042432.*pow(bp,6)+1528823808.*pow(bp,8);
      fact=fact*pow(2.,35)*pow(3.,7)/(5.*bp*pow((49.+144.*pow(bp,2)),11));
      }
 else if (rn == 4) {
//c 3s->5
      fact=87480000000.*(183744069632.+7768401510400.*pow(bp,2)+
                227789475840000.*pow(bp,4)+3915483300000000.*pow(bp,6)-
                41894536500000000.*pow(bp,8)+115844706562500000.*pow(bp,10)-
                103451080078125000.*pow(bp,12)+29192926025390625.*pow(bp,14))/
                [bp*(64.+225*pow(pow(bp,2]),12));
      }
 else if (rn == 5) {
//c 3p->5
      fact=1944.d9*(10020192256.d0+4232282112.d2*pow(bp,2)+
                1934093376.d4*pow(bp,4)-1692872865.d5*pow(bp,6)+
                79597231875.d4*pow(bp,8)-152235703125.d4*pow(bp,10)+
                1167717041015625.d0*pow(bp,12))/(bp*pow((64.+225.*pow(bp,2)),12));
      }
 else if (rn == 6) {
//c 3d->5
      fact=69984.d11*(3014656.d0+382644224.d0*pow(bp,2)-33788808.d2*pow(bp,4
             )+1613793375.d1*pow(bp,6)+10308515625.d0*pow(bp,8)+512578125.d2*pow(bp,10))/
             [bp*(64.+225.*pow(pow(bp,2]),12));
      }
 else if (rn == 7) {
//c 4s->5
      fact=(2414107905106248.d0+122375509655954625.d0*pow(bp,2)+
                141571857480192.d4*pow(bp,4)-42952190122096.d6*pow(bp,6)-
                3750709452544.d8*pow(bp,8)+898924011392.d10*pow(bp,10)-
                446969856.d14*pow(bp,12)+865128448.d14*pow(bp,14)-
                6488064.d16*pow(bp,16)+16384.d18*pow(bp,18));
      fact=fact*pow(2.,27)*pow(5.,7)/(bp*pow((81.+400.*pow(bp,2)),14));
      }
 else if (rn == 8) {
//c 4p->5
      fact=376041246141627.d0+143031288139776.d2*pow(bp,2);
      fact=fact-8084779967808.d4*pow(bp,4)-4296764470272.d6*pow(bp,6);
      fact=fact+870387186176.d8*pow(bp,8)-63160877056.d10*pow(bp,10);
      fact=fact+22278144.d14*pow(bp,12)-33816576.d14*pow(bp,14);
      fact=fact+196608.d16*pow(bp,16);
      fact=fact*pow(2.,23)*pow(5.,10)/(bp*pow((81 + 400*pow(bp,2)),14)) ;
      }
 else if (rn == 9) { 
//c 4->5
      fact=2217735398973.d0-546176812856.d2*pow(bp,2);
      fact=fact+7931167704.d5*pow(bp,4)-37782745600.d5*pow(bp,6);
      fact=fact+97729280000.d5*pow(bp,8)-71884800000.d5*pow(bp,10);
      fact=fact+40960000000.d5*pow(bp,12);
      fact=fact*pow(2.,19)*pow(5.,7)/(bp*pow((81.+400.*pow(bp,2)),11));
      }
 else if (rn == 10) {
//c 4->6
      fact=12881310395.d0+2722585301712.d0*pow(bp,2);
      fact=fact-66947504285952.d0*pow(bp,4)+87094878234624.d1*pow(bp,6);
      fact=fact-396354332491776.d1*pow(bp,8)+960059690975232.d1*pow(bp,10);
      fact=fact-7016971052777472.d0*pow(bp,12)+3851755393646592.d0*pow(bp,14);
      fact=fact*47775744.d0/(bp*pow((25.+144.*pow(bp,2)),12));
      }
      FEX3=fact*corr;
      return;
      }
//c*********************************************************************
