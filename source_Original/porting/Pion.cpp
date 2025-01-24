//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void pion(tck,tcl1,tcl2,tcm1,tcm2,tcn) {
//------------------------------------------------PION.FOR-----------------
// SECTIONS EFFICACES D'IONISATION PWBA                Version du 21/04/94
// d'un projectile à rk,rl1,rl2,rm1 et rm2          -----------------------
// électrons K,2s,2p,3s,p et 3d                  DEBUT DU PROGRAMME PRINCIPAL.
//-------------------------------------------------------------------------
      double DCKE[80],DCL1E[80],DCL2E[80],DCM1E[80],
                DCM2E[80],T[80],WK[80],WL1[80],WL2[80],WM1[80],WM2[80],
                sk[80],sl1[80],sl2[80],sm1[80],sm2[80],
                DCNE[80],Wn[80],sn[80];
      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
      common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m;
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn;
      common/don/zp,E,zt;
      z1=zp;
      z2=zt;
      epo=E;
      rk=(y1s);
      rl1=(y2s);
      rl2=(y2p);
      rm1=(y3s+y3p);
      rm2=(y3d);
      rn=(yn);
      rk0=max[rk][1.];
      rl10=max[rl1][1.];
      rl20=max[rl2][1.];
      rm10=max[rm1][1.];
      rm20=max[rm2][1.];
      rn0=max[rn][1.];
      rkm=rk0-1.;
      rl1m=rl10-1.;
      rl2m=rl20-1.;
      rm1m=rm10-1.;
      rm2m=rm20-1.;
      rnm=rn0-1.;
      Bk=BTot[z1][rk0][rl1][rl2][rm1][rm2][rn]
                -BTot[z1][rkm][rl1][rl2][rm1][rm2][rn];
      Bl1=BTot[z1][rk][rl10][rl2][rm1][rm2][rn]
                -BTot[z1][rk][rl1m][rl2][rm1][rm2][rn];
      Bl2=BTot[z1][rk][rl1][rl20][rm1][rm2][rn]
                -BTot[z1][rk][rl1][rl2m][rm1][rm2][rn];
      Bm1=BTot[z1][rk][rl1][rl2][rm10][rm2][rn]
                -BTot[z1][rk][rl1][rl2][rm1m][rm2][rn];
      Bm2=BTot[z1][rk][rl1][rl2][rm1][rm20][rn]
                -BTot[z1][rk][rl1][rl2][rm1][rm2m][rn];
      Bn=BTot[z1][rk][rl1][rl2][rm1][rm2][rn0]
                -BTot[z1][rk][rl1][rl2][rm1][rm2][rnm];
       if (bl1 < 3.4) bl1=3.4;
       if (bl2 < 3.4) bl2=3.4;
       if(Bm1 <= 1.51) Bm1=1.51;
       if(Bm2 <= 1.51) Bm2=1.51;
       if(Bn <= 0.85) Bn=0.85;
// ************** calcul des constantes *******************************
      betal=(1.-(1./pow((1.+epo/931.5),2.)))*pow((137.036),2);
      CSTE=3.519D4*pow(z2,2.);
//-------- charges effectives ecrantage Slater --------
       if(rk > 1.) {
        zk=z1-0.3;
      }
 else {
        zk=z1;
      }
      zl1=z1-(0.8*rk+0.3*max[(rl1-1.],0.));
      zl2=z1-(rk+0.75*rl1+0.35*max[(rl2-1.],0.));
      zm1=z1-(rk+0.85*(rl1+rl2)+0.35*max[(rm1-1.],0.));
      zm2=z1-(rk+rl1+rl2+rm1+0.35*max[(rm2-1.],0.));
      zn=z1-(rk+rl1+rl2+rm1+rm2+0.35*max[(rn-1.],0.));
       if(zl1 <= 1.) zl1=1.;
       if(zl2 <= 1.) zl2=1.;
       if(zm1 <= 1.) zm1=1.;
       if(zm2 <= 1.) zm2=1.;
       if(zn <= 1.) zn=1.;
// -------------------------------------------------------
//     factors for screening and antiscreening corrections
//........................................................
      vu=137.036*pow((1.-pow((1.+epo/931.5),(-2.))),(0.5));
      zp1=0.9354*zk;
      zp21=0.9354*zl1/2.;
      zp22=0.9354*zl2/2.;
      zp31=0.9354*zm1/3.;
      zp32=0.9354*zm2/3.;
      zp4=0.9534*zn/4.;
       if (vu > zp1) {
      coak=1.-pow((zp1/vu),4);
      }
 else {
      coak=0.;
      }
       if (vu > zp21) {
      coal1=1.-pow((zp21/vu),4);
      }
 else {
      coal1=0.;
      }
       if (vu > zp22) {
      coal2=1.-pow((zp22/vu),4);
      }
 else {
      coal2=0.;
      }
       if (vu > zp31) {
      coam1=1.-pow((zp31/vu),4);
      }
 else {
      coam1=0.;
      }
       if (vu > zp32) {
      coam2=1.-pow((zp32/vu),4);
      }
 else {
      coam2=0.;
      }
       if (vu > zp4) {
      coan=1.-pow((zp4/vu),4);
      }
 else {
      coan=0.;
      }
// ----------------------------------------------------
//     l'ionisation en n=4 est calculée comme étant
//     l'ionisation 2p de z/2
// ----------------------------------------------------
      zn=zn/2.;

      etak=betal/pow((zk),2.);
      etal1=betal/pow((zl1),2.);
      etal2=betal/pow((zl2),2.);
      etam1=betal/pow((zm1),2.);
      etam2=betal/pow((zm2),2.);
      etan=betal/pow((zn),2.);
      tetak=bk/(13.6058*pow(zk,2.));
      tetal1=4*bl1/(13.6058*pow(zl1,2.));
      tetal2=4*bl2/(13.6058*pow(zl2,2.));
      tetam1=9*bm1/(13.6058*pow(zm1,2.));
      tetam2=9*bm2/(13.6058*pow(zm2,2.));
      tetan=16*bn/(13.6058*pow((2.*zn),2.));
//     tetak=1.
//     tetal1=1.
//     tetal2=1.
//     tetam1=1.
//     tetam2=1.
//     tetan=1.
      cstk=cste/(pow(zk,4.));
      cstl1=cste/(pow(zl1,4.));
      cstl2=cste/(pow(zl2,4.));
      cstm1=cste/(pow(zm1,4.));
      cstm2=cste/(pow(zm2,4.));
      cstn=cste/(pow(zn,4.));

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
      scfk=pow(fs1,2)*1.277*pow(z2,(2./3.))/pow(zk,2.);
      scfl1=pow(fs1,2)*1.277*pow(z2,(2./3.))/pow(zl1,2.);
      scfl2=pow(fs1,2)*1.277*pow(z2,(2./3.))/pow(zl2,2.);
      scfm1=pow(fs1,2)*1.277*pow(z2,(2./3.))/pow(zm1,2.);
      scfm2=pow(fs1,2)*1.277*pow(z2,(2./3.))/pow(zm2,2.);
      scfn=pow(fs1,2)*1.277*pow(z2,(2./3.))/pow(zn,2.);
// ...................................................................
//      correction energie de liaison
// ....................................................................
       if (ibin == 1) {
      xk=2.*vu/(zk*tetak);
      gk=(1.+5.*xk+7.14*pow(xk,2.)+4.27*pow(xk,3.)+0.947*pow(xk,4.))/pow((1.+xk),5.);
      ek=pow((1.+zt*gk/(zk*tetak)),2.);
      xl1=4.*vu/(zl1*tetal1);
      gl1=1.+9.*xl1+30.2*pow(xl1,2.)+66.8*pow(xl1,3.)+100.*pow(xl1,4.);
      gl1=gl1+94.1*pow(xl1,5.)+51.3*pow(xl1,6.)+15.2*pow(xl1,7.)+1.891*pow(xl1,8.);
      gl1=gl1/pow((1.+xl1),9.);
      el1=pow((1.+zt*gl1/(zl1*tetal1)),2.);
      xl2=4.*vu/(zl2*tetal2);
      gl2=1.+9.*xl2+34.7*pow(xl2,2.)+81.2*pow(xl2,3.)+112.*pow(xl2,4.);
      gl2=gl2+93.5*pow(xl2,5.)+46.6*pow(xl2,6.)+12.9*pow(xl2,7.)+1.549*pow(xl2,8.);
      gl2=gl2/pow((1.+xl2),9.);
      el2=pow((1.+zt*gl2/(zl2*tetal2)),2.);
      xm1=18.*vu/(zm1*tetam1);
      gm1=(1.+5.*xm1+7.14*pow(xm1,2.)+4.27*pow(xm1,3.)+0.947*pow(xm1,4.));
      gm1=gm1/pow((1.+xm1),5.);
      em1=pow((1.+zt*gm1/(zm1*tetam1)),2.);
      xm2=18.*vu/(zm2*tetam2);
      gm2=(1.+5.*xm2+7.14*pow(xm2,2.)+4.27*pow(xm2,3.)+0.947*pow(xm2,4.));
      gm2=gm2/pow((1.+xm2),5.);
      em2=pow((1.+zt*gm2/(zm2*tetam2)),2.);
//     xn=18.*vu/(zn*tetan)
//     gn=(1.+5.*xn+7.14*xn**2.+4.27*xn**3.+0.947*xn**4.)
//     gn=gn/(1.+xn)**5.
//     en=(1.+zt*gn/(zn*tetan))**2.
      en=1.;
      }
 else {
// ...................................................................
      ek=1.;
      el1=1.;
      el2=1.;
      em1=1.;
      em2=1.;
      en=1.;
      }
// --------------------------------------------------------------------
//   Table des energies cinetiques de l'electron
// ---------------------------------------------
      TMAX1=2194.32*epo;
      TMAX=max[tmax1][4.*bk];
// ---------------------------------------------
       if (TMAX <= 1000.) {
      IMAX=39;
      }
 else if (TMAX <= 2000.) {
      IMAX=41;
      }
 else if (TMAX <= 5000.) {
      IMAX=46;
      }
 else if (TMAX <= 10000.) {
      IMAX=51;
      }
 else if (TMAX <= 20000.) {
      IMAX=53;
      }
 else if (TMAX <= 50000.) {
      IMAX=56;
      }
 else if (TMAX <= 100000.) {
      IMAX=58;
      }
 else if (TMAX <= 200000.) {
      IMAX=60;
      }
 else if (TMAX <= 500000.) {
      IMAX=63;
      }
 else if (TMAX <= 1000000.) {
      IMAX=65;
      }
 else if (TMAX <= 2000000.) {
      IMAX=67;
      }
 else if (TMAX <= 5000000.) {
      IMAX=70;
      }
 else {
      IMAX=73;
      }
      for(/*L101*/ i=1; i<=imax; i++) {
      T[i]=enel[i];
L101: } //continue;
//*********** Energies transferees en unites reduites ******
      for(/*L106*/ I=1; I<=IMAX; I++) {
      WK[I]=ek*tetak+T[I]/(13.6058*pow(zk,2.));
      WL1[I]=el1*tetal1/4.+T[I]/(13.6058*pow(zl1,2.));
      WL2[I]=el2*tetal2/4.+T[I]/(13.6058*pow(zl2,2.));
      WM1[I]=em1*tetam1/9.+T[I]/(13.6058*pow(zm1,2.));
      WM2[I]=em2*tetam2/9.+T[I]/(13.6058*pow(zm2,2.));
      WN[I]=en*tetan/4.+T[I]/(13.6058*pow(zn,2.));
L106: } //continue;
//      **************************************************** 
//      *****   Calcul des DCS et des pertes d'energie *****
//      ****************************************************
      CALL SEKE[WK][IMAX][DCKE];
      CALL SEL1E[WL1][IMAX][DCL1E];
      CALL SEL2E[WL2][IMAX][DCL2E];
      CALL SEM1E[WM1][IMAX][DCM1E];
      CALL SEM2E[WM2][IMAX][DCM2E];
      CALL SENE[WN][IMAX][DCNE];
// --------traitement des donnees du calcul : ---------------------------
// ----------------------------------------------------debut de la boucle
      for(/*L107*/ I=1; I<=IMAX; I++) {
// -----------------------------------
// Sections efficaces differentielles:
// -----------------------------------
      sk[i]=cstk*dcke[i]/(etak*13.6058*pow(zk,2.));
      sl1[i]=cstl1*dcl1e[i]/(etal1*13.6058*pow(zl1,2.));
      sl2[i]=cstl2*dcl2e[i]/(etal2*13.6058*pow(zl2,2.));
      sm1[i]=cstm1*DCM1E[i]/(etam1*13.6058*pow(zm1,2.));
      sm2[i]=cstm2*dcm2e[i]/(etam2*13.6058*pow(zm2,2.));
      sn[i]=cstn*dcne[i]/(etan*13.6058*pow(zn,2.));
// --------------------------------------------------
L107: } //continue;
//----------------------------------------------------------------
// Calcul de la section efficace totale:
// --------------------------------------
      TCk=quad1[Sk][T][IMAX];
      TCl1=quad1[Sl1][T][IMAX];
      TCl2=quad1[Sl2][T][IMAX];
      TCm1=quad1[Sm1][T][IMAX];
      TCm2=quad1[Sm2][T][IMAX];
      TCn=quad1[Sn][T][IMAX];
//     write(16,*)'Bk Bl1 BM1 Bn =', Bk,Bl1,BM1,Bn
//     write(16,*) Tck,TCL1,TCm1,TCn
// ------------------------------------
// on retabli la valeur de zn
      zn=2.*zn;
//--------------------------------------
      return;
        }
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//                                   FIN de BORN.FOR
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
// -----------------------------------
// fonction d'integration sur les T(I):
// -----------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function quad1(SOM,T,IMAX) {
      double SOM[1],T[1];
      quad1=0.;
      for(/*L200*/ I=1; I<=IMAX-1; I++) {
      quad1=quad1+(SOM[I+1]+SOM[I])*(T[I+1]-T[I])/2;
L200: } //continue;
      return;
      }
//-------------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function enel(i) {
      double el[73];
      double /*data*/ (el[i],i=1,73) / 0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,
                6.0,7.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,
                70.0,80.0,90.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,
                600.0,700.0,800.0,900.0,1000.0,1500.0,2000.0,2500.0,3000.0,
                3500.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0,
                15000.0,20000.0,30000.0,40000.0,50000.0,70000.0,100000.0,
                150000.0,200000.0,300000.0,400000.0,500000.0,700000.0,1000000.0,
                1500000.0,2000000.0,3000000.0,4000000.0,5000000.0,7000000.0,
                10000000.0,20000000.0 /;
      enel=el[i];
      return;
      }
// -----------------------------------
// fonction d'integration sur les Q(I):
// -----------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function quad2(SOM2,Q) {
      double SOM2[1],Q[1];
      quad2=0.;
      for(/*L200*/ I=1; I<=57; I++) {
      quad2=quad2+(SOM2[I+1]+SOM2[I])*(Q[I+1]-Q[I])/2.;
L200: } //continue;
      return;
      }
//-------------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function qred(i) {
      double red[58];
      double /*data*/ (red[i],i=1,58) / 0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,
                6.0,7.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,
                70.0,80.0,90.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,
                600.0,700.0,800.0,900.0,1000.0,1500.0,2000.0,2500.0,3000.0,
                3500.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0,
                15000.0,20000.0,30000.0,40000.0,50000.0,70000.0,100000.0 /;
      qred=red[i];
      return;
      }
// -----------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function BTOT(z1,rk,rl1,rl2,rm1,rm2,rn) {
      rks=max[(rk-1.],0.);
      rl1s=max[(rl1-1.],0.);
      rl2s=max[(rl2-1.],0.);
      rm1s=max[(rm1-1.],0.);
      rm2s=max[(rm2-1.],0.);
      rm2s=max[(rm2-1.],0.);
      rns=max[(rn-1.],0.);
      B0=rk*pow((z1-0.3125*rks),2.);
      B1=0.25*rl1*pow((z1-0.8*rk-0.3*rl1s),2.);
      B2=0.25*rl2*pow((z1-rk-0.75*rl1-0.35*rl2s),2.);
      B3=1./9.*rm1*pow((z1-rk-0.85*(rl1+rl2)-0.35*rm1s),2.);
      B4=1./9.*rm2*pow((z1-(rk+rl1+rl2+rm1)-0.35*rm2s),2.);
      B5=1./16.*rn*pow((z1-(rk+rl1+rl2+rm1+rm2)-0.35*rns),2.);
      Btot=13.6058*(b0+b1+b2+b3+b4+b5);
      return;
      }
// -----------------------------------------

//                       Debut des routines de sections efficaces.
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSk(T)/dT                    pour les electrons 1S1/2
// -----------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEKE(WK,IMAX,DCKE) {
      double WK[1],DCKE[1];
      double QP[58],FF[58];
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parf1s/W;
      printf(6,*) '1s ionization';
      QP=0.;
      FF=0.;
      for(/*L300*/ i=1; i<=IMAX; i++) {
        W=WK[i];
        Qm=pow(w,2)/(4.*etak);
//------ recherche du qmax de convergence -> pas d'integration------
        F0=F1s[Qm];
        Q=Qm;
L299:     Q=2.*Q;
          Fn=ANINT[100.*F1s[Q]/F0];
           if (Fn >= 1.) goto L299;
        PQ=Q/1000.;
//------------------------------------------------------------------
        for(/*L301*/ j=1; j<=39; j++) {
          QP[j]=qred[j]*PQ;
          Q=Qm+QP[j];
          FF[j]=F1s[Q];
L301:   } //continue;
        DCKE[i]=QUAD2[FF][QP];
L300: } //continue;
        return;
      }
// ~~~~~~~~~~~~~~~~~~~~~~~
// FONCTION F1S FACTEUR DE FORME POUR ELECTRON 1S
// ~~~~~~~~~~~~~~~~~~~~~~~
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function F1S(Q) {
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parf1s/W;
      Tol0=-1D-3;
      Tol1=1D-3;
      Pi=4.*atan(1.);
      sq=pow((1.+scfk/Q),(-2.))+coak*(1.-pow((1+Q/scfk),(-2.)))/z2;
      AK2=W-1;
      AS=Q+(AK2)/3.+(1./3.);
            if(AK2 < Tol0) {
      ACC=sqrt(fabs(AK2));
      BS=pow(((Q+pow((1.-ACC),2.))/(Q+pow((1.+ACC),2.))),(1./ACC));
      CS=1.;
          }
 else if(AK2 < 0) {
      ACC=sqrt(fabs(AK2));
      BS=exp(-4./(Q+pow((1.+ACC),2)));
      CS=1.;
          }
 else if(AK2 < Tol1) {
      BS=exp(-4./(Q+1.-AK2));
      CS=1.;
          }
 else {
      ACC=sqrt(AK2);
      BS=exp((-2./ACC)*(ATAN2[(2.*ACC],(Q-AK2+1.))));
      CS=(1.-exp((-2.*Pi)/ACC));
             }
      DS=pow((pow((Q-AK2+1.),2.)+4.*AK2),3.);
      ES=pow(2.,7.);
      F1S=sq*(ES*AS*BS)/(CS*DS*Q);
        return;
      }
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSl1(T)/dT                   pour les electrons 2S
// -----------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEL1E(WL1,IMAX,DCL1E) {
      double WL1[1],DCL1E[1];
      double QP[58],FF[58];
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parf2s/W1;
      printf(6,*) '2s ionization';
      for(/*L400*/ i=1; i<=IMAX; i++) {
        W1=WL1[i];
        Qm=pow(W1,2)/(4.*etal1);
//------ recherche du qmax de convergence -> pas d'integration------
        F0=F2s[Qm];
        Q=Qm;
L399:     Q=2.*Q;
          Fn=ANINT[100.*F2s[Q]/F0];
           if (Fn >= 1.) goto L399;
        PQ=Q/1000.;
//------------------------------------------------------------------
        for(/*L401*/ j=1; j<=39; j++) {
          QP[j]=qred[j]*PQ;
          Q=Qm+QP[j];
          FF[j]=F2s[Q];
L401:   } //continue;
        DCL1E[i]=QUAD2[FF][QP];
L400: } //continue;
        return;
      }
// ~~~~~~~~~~~~~~~~~~~~~~~
// FONCTION F2S FACTEUR DE FORME POUR ELECTRON 2S
// ~~~~~~~~~~~~~~~~~~~~~~~
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function F2S(Q) {
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parf2s/W1;
      Tol0=-1D-3;
      Tol1=1D-3;
      Pi=4.*atan(1.);
      sq=pow((1.+scfl1/Q),(-2.))+coal1*(1.-pow((1+Q/scfl1),(-2.)))/z2;
      AK2=W1-0.25;
            if(AK2 < Tol0) {
      ACC=sqrt(fabs(AK2));
      BS=pow(((Q+pow((0.5-ACC),2.))/(Q+pow((0.5+ACC),2.))),(1./ACC));
      CS=1.;
          }
 else if(AK2 < 0) {
      ACC=sqrt(fabs(AK2));
      BS=exp(-2./(Q+pow((0.5+ACC),2)));
      CS=1.;
          }
 else if(AK2 < Tol1) {
      BS=exp(-2./(Q+0.25-AK2));
      CS=1.;
          }
 else {
      ACC=sqrt(AK2);
      BS=exp((-2./ACC)*(ATAN2[(ACC],(Q-AK2+0.25))));
      CS=(1.-exp((-2.*Pi)/ACC));
           }
      AS=pow((Q-AK2+0.25),2.)+AK2;
      ES=pow(2.,4.);
      AL1=BS*ES/(CS*pow(AS,5.)*Q);
      A5=pow(Q,5.);
      A4=-(8./3.+11./3.*AK2)*pow(Q,4.);
      A3=(41./24.+6.*AK2+14./3.*pow(AK2,2.))*pow(Q,3.);
      A2=(5./48.-31./24.*AK2-10./3.*pow(AK2,2.)-2.*pow(AK2,3.))*pow(Q,2.);
      A1=(47./3840.-41./120.*pow(AK2,2.)-2./3.*pow(AK2,3.)-1./3.*pow(AK2,4))*Q;
      A0=1./768.+17./768.*AK2+7./48.*pow(AK2,2.)+11./24.*pow(AK2,3.
                )+2./3.*pow(AK2,4.)+1./3.*pow(AK2,5.);
      F2S=sq*AL1*(A5+A4+A3+A2+A1+A0);
        return;
      }
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSl2(T)/dT                   pour les electrons 2P
// -----------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEL2E(WL2,IMAX,DCL2E) {
      double WL2[1],DCL2E[1];
      double QP[58],FF[58];
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parf2p/W2;
      printf(6,*) '2p ionization';
      for(/*L500*/ i=1; i<=IMAX; i++) {
        W2=WL2[i];
        Qm=pow(w2,2)/(4.*etal2);
//------ recherche du qmax de convergence -> pas d'integration------
        F0=F2p[Qm];
        Q=Qm;
L499:     Q=2.*Q;
          Fn=ANINT[100.*F2p[Q]/F0];
           if (Fn >= 1.) goto L499;
        PQ=Q/1000.;
//------------------------------------------------------------------
        for(/*L501*/ j=1; j<=39; j++) {
          QP[j]=qred[j]*PQ;
          Q=Qm+QP[j];
          FF[j]=F2p[Q];
L501:   } //continue;
        DCL2E[i]=QUAD2[FF][QP];
L500: } //continue;
        return;
      }
// ~~~~~~~~~~~~~~~~~~~~~~~~~
//    FONCTION F2P FACTEUR DE FORME POUR ELECTRON 2P
// ~~~~~~~~~~~~~~~~~~~~~~~~~
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function F2P(Q) {
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parf2p/W2;
      Tol0=-1D-3;
      Tol1=1D-3;
      Pi=4.*atan(1.);
      sq=pow((1.+scfl2/Q),(-2.))+coal2*(1.-pow((1+Q/scfl2),(-2.)))/z2;
      AK2=W2-0.25;
            if(AK2 < Tol0) {
      ACC=sqrt(fabs(AK2));
      BS=pow(((Q+pow((0.5-ACC),2.))/(Q+pow((0.5+ACC),2.))),(1./ACC));
      CS=1.;
          }
 else if(AK2 < 0) {
      ACC=sqrt(fabs(AK2));
      BS=exp(-2./(Q+pow((0.5+ACC),2)));
      CS=1.;
          }
 else if(AK2 < Tol1) {
      BS=exp(-2./(Q+0.25-AK2));
      CS=1.;
          }
 else {
      ACC=sqrt(AK2);
      BS=exp((-2./ACC)*(ATAN2[(ACC],(Q-AK2+0.25))));
      CS=(1.-exp((-2.*Pi)/ACC));
           }
      AS=pow((Q-AK2+0.25),2.)+AK2;
      ES=pow(2.,4.);
      AL2=BS*ES/(CS*pow(AS,5.)*Q);
      A4=9./4.*pow(Q,4.);
      A3=-(0.75+3.*AK2)*pow(Q,3.);
      A2=(19./32.-0.75*AK2-0.5*pow(AK2,2.))*pow(Q,2.);
      A1=(107./960.+41./48.*AK2+113./60.*pow(AK2,2.)+pow(AK2,3.))*Q;
      A0=1./4.*pow(AK2,4.)+5./12.*pow(AK2,3.)+7./32.*pow(AK2,2.)+3./64.*AK2
                +11./3072.;
      F2P=sq*AL2*(A4+A3+A2+A1+A0)*(1./3.);
      return;
      }
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSm1(T)/dT                   pour les electrons 3s,p
// -----------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEM1E(WM1,IMAX,DCM1E) {
      double WM1[1],DCM1E[1];
      double QP[58],FF[58];
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parfm1/W1M;
      printf(6,*) '3s,p ionization';
      for(/*L600*/ i=1; i<=IMAX; i++) {
        W1M=WM1[i];
        Qm=pow(w1m,2)/(4.*etam1);
//------ recherche du qmax de convergence -> pas d'integration------
        F0=Fm1[Qm];
        Q=Qm;
L599:     Q=2.*Q;
          Fn=ANINT[100.*Fm1[Q]/F0];
           if (Fn >= 1.) goto L599;
        PQ=Q/1000.;
//------------------------------------------------------------------
        for(/*L601*/ j=1; j<=39; j++) {
          QP[j]=qred[j]*PQ;
          Q=Qm+QP[j];
          FF[j]=Fm1[Q];
L601:   } //continue;
        DCM1E[i]=QUAD2[FF][QP];
L600: } //continue;
        return;
      }
// ~~~~~~~~~~~~~~~~~~~~~~~~~
// FONCTION Fm1 FACTEUR DE FORME POUR ELECTRON 3s,p
// ~~~~~~~~~~~~~~~~~~~~~~~~~
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function Fm1(Q) {
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parfm1/W1M;
      Tol0=-1D-3;
      Tol1=1D-3;
      Pi=4.*atan(1.);
      sq=pow((1.+scfm1/Q),(-2.))+coam1*(1.-pow((1+Q/scfm1),(-2.)))/z2;
      AK2=W1M-1/9.;
            if(AK2 < Tol0) {
      ACC=sqrt(fabs(AK2));
      BS=pow(((Q+pow((1./3.-ACC),2.))/(Q+pow((1./3.+ACC),2.))),(1./ACC));
      CS=1.;
          }
 else if(AK2 < 0) {
      ACC=sqrt(fabs(AK2));
      BS=exp(-4./3./(Q+pow((1./3.+ACC),2)));
      CS=1.;
          }
 else if(AK2 < Tol1) {
      BS=exp(-4./3./(Q+1./9.-AK2));
      CS=1.;
          }
 else {
      ACC=sqrt(AK2);
      BS=exp((-2./ACC)*(ATAN2[(2.*ACC/3.],(Q-AK2+1./9.))));
      CS=(1.-exp((-2.*Pi)/ACC));
           }
      AS=pow((Q-AK2+1./9.),2.)+4.*AK2/9.;
      ES=pow(2.,7.)/27.;
      AL3=BS*ES/(Q*CS*pow(AS,5.));
      A5=pow(Q,5.);
      A4=-(43./27.+11./3.*AK2)*pow(Q,4.);
      A3=(518./243.+412./81.*AK2+14./3.*pow(AK2,2.))*pow(Q,3.);
      A2=-(442./729.+310./81.*AK2+122./27.*pow(AK2,2.)+2.*pow(AK2,3.))*pow(Q,2.);
      A1=(1943./pow(3.,9.)+1460./pow(3.,7.)*AK2+290./243.*pow(AK2,2.
                )+4./27.*pow(AK2,3.)-1./3.*pow(AK2,4.))*Q;
      A0=1./3.*pow(AK2,5.)+71./81.*pow(AK2,4.)+62./81.*pow(AK2,3.
                )+1790./pow(3.,8.)*pow(AK2,2.)+2431./pow(3.,10.)*AK2+377./pow(3.,11.);
      Fm1=sq*AL3*(A5+A4+A3+A2+A1+A0)*(1./9.);
      return;
      }
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSm2(T)/dT                   pour les electrons 3d
// -----------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEM2E(WM2,IMAX,DCM2E) {
      double WM2[1],DCM2E[1];
      double QP[58],FF[58];
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parfm2/W2M;
      printf(6,*) '3d ionization';
      for(/*L600*/ i=1; i<=IMAX; i++) {
        W2M=WM2[i];
        Qm=pow(w2m,2)/(4.*etam2);
//------ recherche du qmax de convergence -> pas d'integration------
        F0=Fm2[Qm];
        Q=Qm;
       if(f0 > 0.) {
L599:     Q=2.*Q;
          Fn=ANINT[100.*Fm2[Q]/F0];
           if (Fn >= 1.) goto L599;
      }
        PQ=Q/1000.;
//------------------------------------------------------------------
        for(/*L601*/ j=1; j<=39; j++) {
          QP[j]=qred[j]*PQ;
          Q=Qm+QP[j];
          FF[j]=Fm2[Q];
L601:   } //continue;
        DCM2E[i]=QUAD2[FF][QP];
L600: } //continue;
        return;
      }
// ~~~~~~~~~~~~~~~~~~~~~~~~~
// FONCTION Fm2 FACTEUR DE FORME POUR ELECTRON 3d
// ~~~~~~~~~~~~~~~~~~~~~~~~~
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function Fm2(Q) {
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parfm2/W2M;
      Tol0=-1D-3;
      Tol1=1D-3;
      Pi=4.*atan(1.);
      sq=pow((1.+scfm2/Q),(-2.))+coam2*(1.-pow((1+Q/scfm2),(-2.)))/z2;
      AK2=W2M-1/9.;
            if(AK2 < Tol0) {
      ACC=sqrt(fabs(AK2));
      BS=pow(((Q+pow((1./3.-ACC),2.))/(Q+pow((1./3.+ACC),2.))),(1./ACC));
      CS=1.;
          }
 else if(AK2 < 0) {
      ACC=sqrt(fabs(AK2));
      BS=exp(-4./3./(Q+pow((1./3.+ACC),2)));
      CS=1.;
          }
 else if(AK2 < Tol1) {
      BS=exp(-4./3./(Q+1./9.-AK2));
      CS=1.;
          }
 else {
      ACC=sqrt(AK2);
      BS=exp((-2./ACC)*(ATAN2[(2.*ACC/3.],(Q-AK2+1./9.))));
      CS=(1.-exp((-2.*Pi)/ACC));
           }
      AS=pow((Q-AK2+1./9.),2.)+4.*AK2/9.;
      ES=pow(2.,7.)/27.;
      AL3=BS*ES/(Q*CS*pow(AS,5.));
      A5=pow(Q,5.);
      A4=-(43./27.+11./3.*AK2)*pow(Q,4.);
      A3=(518./243.+412./81.*AK2+14./3.*pow(AK2,2.))*pow(Q,3.);
      A2=-(442./729.+310./81.*AK2+122./27.*pow(AK2,2.)+2.*pow(AK2,3.))*pow(Q,2.);
      A1=(1943./pow(3.,9.)+1460./pow(3.,7.)*AK2+290./243.*pow(AK2,2.
                )+4./27.*pow(AK2,3.)-1./3.*pow(AK2,4.))*Q;
      A0=1./3.*pow(AK2,5.)+71./81.*pow(AK2,4.)+62./81.*pow(AK2,3.
                )+1790./pow(3.,8.)*pow(AK2,2.)+2431./pow(3.,10.)*AK2+377./pow(3.,11.);
      Fm2=sq*AL3*(A5+A4+A3+A2+A1+A0)*(1./9.);
      return;
      }
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
//     dSN(T)/dT                    pour les electrons N -> "2P"
// -----------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SENE(WN,IMAX,DCNE) {
      double WN[1],DCNE[1];
      double QP[58],FF[58];
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parfn/W4;
      printf(6,*) 'n=4 ionization';
      for(/*L500*/ i=1; i<=IMAX; i++) {
        W4=Wn[i];
        Qm=pow(w4,2)/(4.*etan);
//------ recherche du qmax de convergence -> pas d'integration------
        F0=FN[Qm];
        Q=Qm;
L499:     Q=2.*Q;
          F4=ANINT[100.*FN[Q]/F0];
           if (F4 >= 1.) goto L499;
        PQ=Q/1000.;
//------------------------------------------------------------------
        for(/*L501*/ j=1; j<=39; j++) {
          QP[j]=qred[j]*PQ;
          Q=Qm+QP[j];
          FF[j]=FN[Q];
L501:   } //continue;
        DCNE[i]=QUAD2[FF][QP];
L500: } //continue;
        return;
      }
// ~~~~~~~~~~~~~~~~~~~~~~~~~
//    FONCTION FN FACTEUR DE FORME POUR ELECTRON 2P -> n=4
// ~~~~~~~~~~~~~~~~~~~~~~~~~
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function FN(Q) {
      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
      common/par3/coak,coal1,coal2,coam1,coam2,coan;
      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
           tetan;
      common/parfn/W4;
      Tol0=-1D-3;
      Tol1=1D-3;
      Pi=4.*atan(1.);
      sq=pow((1.+scfn/Q),(-2.))+coan*(1.-pow((1+Q/scfn),(-2.)))/z2;
      AK2=W4-0.25;
            if(AK2 < Tol0) {
      ACC=sqrt(fabs(AK2));
      BS=pow(((Q+pow((0.5-ACC),2.))/(Q+pow((0.5+ACC),2.))),(1./ACC));
      CS=1.;
          }
 else if(AK2 < 0) {
      ACC=sqrt(fabs(AK2));
      BS=exp(-2./(Q+pow((0.5+ACC),2)));
      CS=1.;
          }
 else if(AK2 < Tol1) {
      BS=exp(-2./(Q+0.25-AK2));
      CS=1.;
          }
 else {
      ACC=sqrt(AK2);
      BS=exp((-2./ACC)*(ATAN2[(ACC],(Q-AK2+0.25))));
      CS=(1.-exp((-2.*Pi)/ACC));
           }
      AS=pow((Q-AK2+0.25),2.)+AK2;
      ES=pow(2.,4.);
      AL2=BS*ES/(CS*pow(AS,5.)*Q);
      A4=9./4.*pow(Q,4.);
      A3=-(0.75+3.*AK2)*pow(Q,3.);
      A2=(19./32.-0.75*AK2-0.5*pow(AK2,2.))*pow(Q,2.);
      A1=(107./960.+41./48.*AK2+113./60.*pow(AK2,2.)+pow(AK2,3.))*Q;
      A0=1./4.*pow(AK2,4.)+5./12.*pow(AK2,3.)+7./32.*pow(AK2,2.)+3./64.*AK2
                +11./3072.;
      FN=sq*AL2*(A4+A3+A2+A1+A0)*(1./3.);
      return;
      }
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//                                  FIN DES ROUTINES DES SECTIONS EFFICACES.
//                           ----------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//           FIN DU PROGRAMME
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
