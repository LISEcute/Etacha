//********1*********1*********1*********1*********1*********1********1**
// Pour ETACHA4
//.....calcula la secc. eficaz total de excitacion desde los
//.....estados ns, np1 y np-1 al estado final nlm, por S.Eikonal
//.....Autor: Cesar A. Ramirez            Fecha: 10-12-2012
//         zt=carga del nucleo blanco
//         zp=carga del nucleo proyectil 
//         ei=energia del estado ligado inicial
//         ef=energia del estado ligado final
//         n,l,m = numeros cuanticos del estado final
//         eta=componente normal del momento transferido
//**********************************************************************
//----- Program to calculte total cross sections of 
//----- monoelectronic atoms excitation, for colision with nude projectiles 
//----- in the SE (Symmetric-Eikonal) aproximation.
//--------------------------------------------------
//......changes in the present version (JPR 02/2012)
//----- zp: Projectile charge (to be excited)
//----- zt: Target charge
//......SC and ASC corrections
//-------------------------------------------------
//----- (n0,l0,m0) quantum numbers of initial state 
//----- (n,l,m) quantum numbers of final state
//----- SUBROUTINES:
//----- * rns, rnp0, rnp1, rnp1m: excitation from ns, np0, np1 np-1 to n'lm
//-----                    respectively with n'lm arbitraries (n'>n)
//----- * intdef: for definite integrals from 0 to infinity
//----- * fac: factorial function
//----- * cg: Clebsch-Gordan coeficients
//----- * armonic: spherical harmonic function
//----- * plgndr: Legendre polinomials
//----- * gamac: gama function
//----- * cf21d, chyper, hypgfx: hypergeometric functions F21(ca,cb,cc,cx)
//**********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEnlm(zp8,zt8,E8) {
      // implicit double*8[a-b][d-h][o-z],complex*16[c];
      double*4 SeSE,StSE;
      char kk*2;
//-----------------------------------------------------------------
      common/m1/zp,zt;
      common/m2/epsi,v,pi;
      common/m3/ca1,ca2,ca3,ca4,calfa;
      common/m4/ci,cinu,cuno,cdos,con0;
      common/m5/n0,l0,m0,n,l,m;
      common/asc/sc,scoa,icor;
      common/SecSE/SeSE,StSE;
      double sexnlm[250];
      double signlm[4][2][3][5][5][6],signl[4][2][5][5],signt[4][2][5];
      double SeSE[66],StSE[12];
      external rns;
      external rnp0;
      external rnp1;
      external rnp1m;
//-----------------------------------------------------------------
      zp=zp8;
      zt=zt8;
      E=E8;
      bet=dble[(1.-pow((1.+E/pow(931.5],(-2.))),(0.5)));
      v=137.036d0*bet;
//----------------------------------------------------------------------
      scf=dble[1.13*pow(zt,(1./3.]));
      fs1=1.d0;
       if (zt == 1.) fs1=1.23d0;
       if (zt == 2.) fs1=1.526d0;
       if (zt == 6.) fs1=0.78d0;
       if (zt == 7.) fs1=0.85d0;
       if (zt == 10.) fs1=1.04d0;
       if (zt == 13.) fs1=0.58d0;
       if (zt == 14.) fs1=0.59d0;
       if (zt == 18.) fs1=0.68d0;
       if (zt == 29.) fs1=0.672d0;
       if (zt == 36.) fs1=0.61d0;
       if (zt == 54.) fs1=0.535d0;
      sc=scf*fs1;
      icor=1;
//----------------------------------------------------------------------
      pi=dacos[-1.d0] ;
      xnu=zt/v;

      ci=(0.d0,1.d0);
      cuno=(1.d0,0.d0);
      cdos=(2.d0,0.d0);
      cinu=ci*xnu;
      ca1=1.d0+cinu;
      ca2=1.d0-cinu;
      ca3=1.d0+2.d0*cinu;
      ca4=1.d0-2.d0*cinu;
      senhy=dsinh[pi*xnu];
        test1=1.d-6;
        test2=1.d-6;
      call gamac[cinu][1][test1][test2][nt][cloga];
        cgama=cdexp[cloga];
//----------------------------------------------------------------------- 
//     open(40,file='tnlm4.dat',status='unknown')
//      write(40,200)zp,zt,e
//     if (icor.eq.1) then
//     write(40,*)' with SC+ASC corrections'
//     else
//     write(40,*)' no   SC+ASC corrections'
//     endif
//200   format(' zp=',d15.5,' zt=',d15.5,' E=',e15.5)
//L201: format(' initial=',3i3,'        final=',3i3);
//-----------------------------------------------------------------------
      index=1;
      for(/*L100*/ n0=1; n0<=4; n0++) {
      l0m=min[n0-1][1];
      for(/*L100*/ l0=0; l0<=l0m; l0++) {
      for(/*L100*/ m0=0; m0<=l0; m0++) {
      nfm=n0+1;
      for(/*L110*/ n=nfm; n<=5; n++) {
      for(/*L120*/ l=0; l<=n-1; l++) {
      for(/*L130*/ m=0; m<=l; m++) {
      printf(*,201)n0,l0,m0,n,l,m;
//----------------------------------------------------------------------
      n10=n0-l0-1;
      n20=n0+l0;
      fac10=fac[n10];
      fac20=fac[n20];
      xn0=2.d0/pow(n0,2)*pow(sqrtl(zp),3)*dsqrt[fac20*fac10]*pow((2*zp/n0),l0);
      n1=n-l-1;
      n2=n+l;
      fac1=fac[n1];
      fac2=fac[n2];
      xnl=2.d0/pow(n,2)*pow(sqrtl(zp),3)*sqrtl(fac2*fac1)*pow((2*zp/n),l);
//----------------------------------------------------------------------
      ei=-pow((zp/n0),2)/2.d0;
      ef=-pow((zp/n),2)/2.d0;
          epsi=(ef-ei)/v;
           if(epsi == 0.d0)epsi=1.d-6;
//.............................
       if (v > 1.75d0*epsi) {
      scoa=1.d0-pow((1.75d0*epsi/v),2);
      }
 else {
      scoa=0.0d0;
      }
//.............................
          calfa=2.d0*ci*v*epsi;
          beta=pow((2.d0*v*epsi),2);
//----------------------------
      con0=4.d0*pi*v*pow((pi*xnu/cgama/senhy),2)/beta*xn0*xnl;
//     *        /cdexp(cinu*dlog(beta))         cte. de módulo=1
//----------------------------
      ymin=0.d0;
      ymax=3.d0;
        a02=2.8d-17;
        pre=1.d-3;
//---------------------------------------------
       if(l0 == 0){
      call intdef[pre][ymin][ymax][rns][xint];
      tnlm=2*pi*xint*a02;
      }
 else {
         if(m0 == 0){
          call intdef[pre][ymin][ymax][rnp0][xint];
          tnlm=2*pi*xint*a02;
        }
 else {
          call intdef[pre][ymin][ymax][rnp1][xint];
          tp1=2*pi*xint*a02;
           if(m != 0)tp1=2*tp1;
          call intdef[pre][ymin][ymax][rnp1m][xint];
          tp1m=2*pi*xint*a02;
           if(m != 0)tp1m=2*tp1m;
          tnlm=tp1+tp1m;

        }
      }
//--------------------------------------------------------------------
//---- Se multiplica por 2 la sección eficaz correspondiente a 
//---- la transición desde un estado inicial con m=0 
//---- a un estado final con m distinto a cero, para tener en cuenta los 
//---- dos casos simétricos +m y -m. 
//--------------------------------------------------------------------
       if(m0 == 0 && m != 0)tnlm=2*tnlm;
//      write(40,300)n0,l0,m0,n,l,m,tnlm,index,sc,scoa
      sexnlm[index]=tnlm;
      signlm[n0][l0+1][m0+1][n][l+1][m+1]=tnlm;
      index=index+1;
L130: } //continue;
L120: } //continue;
L110: } //continue;
L100: } //continue;
//     somme sur les m0 et m
//     write(40,*)' '
      index=1;
      for(/*L140*/ n0=1; n0<=4; n0++) {
      l0m=min[n0-1][1];
      nfm=n0+1;
      for(/*L140*/ l0=0; l0<=l0m; l0++) {
      for(/*L140*/ n=nfm; n<=5; n++) {
      for(/*L140*/ l=0; l<=n-1; l++) {
      sigt=0.0;
      for(/*L150*/ m0=0; m0<=l0m; m0++) {
      for(/*L150*/ m=0; m<=l; m++) {
      sigt=sigt+signlm[n0][l0+1][m0+1][n][l+1][m+1];
L150: } //continue;
       if(l0 == 1) sigt=sigt/3.;
      signl[n0][l0+1][n][l+1]=sigt;
      SeSE[index]=sigt;
      index=index+1;
//      write(40,310)n0,l0,n,l,sigt
L140: } //continue;
//     write(40,*)' '
//     somme sur les l pour n=4 et 5
      index=1;
      for(/*L160*/ n=4; n<=5; n++) {
      nim=n-1;
      for(/*L160*/ n0=1; n0<=nim; n0++) {
      l0m=min[n0-1][1];
      for(/*L160*/ l0=0; l0<=l0m; l0++) {
      sigtn=0.;
      for(/*L170*/ l=0; l<=n-1; l++) {
      sigtn=sigtn+signl[n0][l0+1][n][l+1];
L170: } //continue;
      signt[n0][l0+1][n]=sigtn;
      StSE[index]=sigtn;
      index=index+1;
//      write(40,311)n0,l0,n,sigtn
L160: } //continue ;
//      close(unit=40)
//-----
//L300: format(' inicial=',3i3,'   final=',3i3,'  tot = ',d15.5,
           '   index=',i3,' coef ASC=',2d15.5);
//L310: format(' inicial=',2i3,'      final=',2i3,'     tot = ',d15.5);
//L311: format(' inicial=',2i3,'      final=',i3,'        tot = ',d15.5);
      return;
      }
//*************************************** fin  programa principal ******

//     SE    ns-nlm.for

//***********************************************************************
//.....function=eta * /R(eta)/
          2      en la aprox. SE 
//.....desde ns a nlm
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void rns(eta,xf) {
      // implicit double*8[a-b][d-h][o-z], complex*16[c];
      common/m1/zp,zt;
      common/m2/epsi,v,pi;
      common/m3/ca1,ca2,ca3,ca4,calfa;
      common/m4/ci,cinu,cuno,cdos,con0;
      common/m5/n0,l0,m0,n,l,m;
      common/asc/sc,scoa,icor;
//-----------------------------------------------------------------------
      cero=(0.d0,0.d0);
      gama=pow(eta,2)+pow(epsi,2);
      xk=sqrtl(gama);
      znz=zp/n+zp/n0;
      vari=-pow((eta/epsi),2);
      call cf21d[cinu][cinu][cdos][vari][chyp1];
      call cf21d[ca1][ca1][cdos][vari][chyp2];
      chyp1=(ca2*chyp1-cinu*(1.d0-vari)*chyp2)/ca4;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      csum2=cero;
      csum3=cero;
       if(l == 0){
         minla=0;
      }
 else {
         minla=l-1;
      }
         maxla=l+1;
//-----
      for(/*L*/ la=minla; la<=maxla; la++) {
        ga=sqrtl(pi);
        for(/*L*/ j=0; j<=la; j++) {
        ga=ga*(0.5d0+j);
        } for(/*L*/ ) {
      xla=sqrtl(pi)*pow(xk,la)/pow(2.d0,(la+1))/ga;
        sumq1=cero;
        sumq=cero ;
        maxiq=n-l-1;
      for(/*L*/ iq=0; iq<=maxiq; iq++) {
         n4=n-l-1-iq;
         n5=iq;
         n6=2*l+1+iq;
      fac4=fac[n4];
      fac5=fac[n5];
      fac6=fac[n6];
      xq=pow( (-2.d0*zp/n),iq)/fac4/fac5/fac6;
        sump1=cero ;
        sump=cero ;
        maxip=n0-1;
      for(/*L*/ ip=0; ip<=maxip; ip++) {
        n40=n0-1-ip;
        n50=ip;
        n60=1+ip;
      fac40=fac[n40];
      fac50=fac[n50];
      fac60=fac[n60];
      xp =pow( (-2.d0*zp/n0),ip)/fac40/fac50/fac60;
//-------------------------------------------------
      xmu=2.5d0+l+iq+ip;
      xnu=la+0.5d0;
      xa=(xmu+xnu)/2.d0;
      xc=xnu+1.d0;
      xz=-gama/pow(znz,2);
      xc1=2.d0*xa;
      xc2=2.d0*xa-xc+1.d0;
      yz=(1.d0-sqrtl(1.d0-xz))/(1.d0+dsqrt[1.d0-xz])*cuno ;
      n3 = 3+l+la+iq+ip;
      h1=pow((2.d0/(1.d0+sqrtl(1.d0-xz))),n3 );
      call hygfx[xc1][xc2][xc][yz][hyp31];
       nn2 = 2+l+la+iq+ip;
       fac3 = fac[nn2];
       xintx1 = xla*fac3/pow(znz,n3)*h1*hyp31;
       if(la == l){
      sump1 = sump1 + xp*xintx1;
      }
 else {
       xc1=2.d0*xa-1.d0 ;
       xc2=2.d0*xa-xc ;
       h2=pow((2.d0/(1.d0+sqrtl(1.d0-xz))),(n3-1));
      call hygfx[xc1][xc2][xc][yz][hyp32] ;
       xintx2 = xla*fac3/pow(znz,n3)*h2*hyp32;
      sump = sump + xp*(-zp/n0*xintx1 + znz*ip/(2.d0*xa-1.d0)*xintx2);
      }
      } for(/*L*/ ) {
       if(la == l){
        sumq1 = sumq1+sump1*xq;
      }
 else {
        sumq = sumq+sump*xq;
      }
      } for(/*L*/ ) {
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//              s=-1      indica que es un armonico conjugado
      s=-1.d0;
      xkz=-epsi;
      tk=dacos[xkz/xk];
      calga=calfa/gama;
      xl=sqrtl((2.d0*la+1.d0)/(2.d0*l+1.d0)/4.d0/pi);
       if(la == l){ 
        call armonic[s][l][m][tk][0.d0][cylam];
        csum1=pow( ci,la)*xl*sumq1*cylam*calfa*(ca2/cinu*chyp1+chyp2);
      }
 else {
        call cg[la][1][l][0][0][0][xc3];
        call cg[la][1][l][m+1][-1][m][xc4];
        call cg[la][1][l][m-1][1][m][xc5];
        call cg[la][1][l][m][0][m][xc6];
        call armonic[s][la][m+1][tk][0.d0][cylam1] ;
        call armonic[s][la][m-1][tk][0.d0][cylam2] ;
        call armonic[s][la][m][tk][0.d0][cylam3] ;
      csum2=csum2+pow(ci,la)*xl*sumq*xc3*(xc4*cylam1-xc5*cylam2)/sqrtl(2.d0) 
               *ci*eta*calga*ca2/cinu*chyp1 ;
      csum3=csum3 +pow( ci,la)*xl*sumq*xc3*xc6*cylam3 
               *(2.d0*v*chyp2-ci*epsi*calga*ca2/cinu*chyp1) ;
      }
      } for(/*L*/ ) {
//-----------------------------------------------------------------------
      csum= 4.d0*pi* (csum1 - 2.d0 * (csum2+csum3));
//-----------------------------------------------------------------------
      cf= ci/(2.d0*pi*v) * con0 * csum ;
//     *         * cdexp(2*cinu*dlog(gama))    cte de modulo=1
      xf=eta*cdabs[cf]*cdabs[cf];
//............................... 
//     SC and ASC corrections
       if (icor == 1) {
      sq=pow((1.d0+pow((sc/xk),2)),(-2))+scoa*(1.d0-pow((1.d0+pow((xk/sc),2)),(-2)))/zt;
      xf=xf*sq;
      }
//...............................
      return;
      } ;
//***********************************************************************

//   SE  np0-nlm.for

//***********************************************************************
//.....function=eta * /R(eta)/
          2      en la aprox. SE 
//.....desde np0 a nlm
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void rnp0(eta,xf)  {
      // implicit double*8[a-b][d-h][o-z],complex*16[c];
      common/m1/zp,zt;
      common/m2/epsi,v,pi;
      common/m3/ca1,ca2,ca3,ca4,calfa;
      common/m4/ci,cinu,cuno,cdos,con0;
      common/m5/n0,l0,m0,n,l,m;
      common/asc/sc,scoa,icor;
//-----------------------------------------------------------------------
      cero=(0.d0,0.d0);
      gama=pow(eta,2)+pow(epsi,2);
      xk=sqrtl(gama);
      znz=zp/n+zp/n0;
//-----
      vari=-pow(eta,2)/pow(epsi,2);
      call cf21d[cinu][cinu][cdos][vari][chyp1];
      call cf21d[ca1][ca1][cdos][vari][chyp2];
      chyp1=(ca2*chyp1-cinu*(1.d0-vari)*chyp2)/ca4;
      csum1=cero;
      csum2=cero;
      csum3=cero;
        minla2=abs(l-2);
        minla1=abs(l-1);
      minla=min[minla1][minla2];
         if(l == 0)minla=0;
        maxla=l+2;
      for(/*L*/ la=minla; la<=maxla; la++) {
      sumq1=0.d0 ;
      sumq2=0.d0 ;
      sumq=0.d0 ;
         maxiq=n-l-1;
      for(/*L*/ iq=0; iq<=maxiq; iq++) {
         n4=n-l-1-iq;
         n5=iq;
         n6=2*l+1+iq;
      fac4=fac[n4];
      fac5=fac[n5];
      fac6=fac[n6];
      xq=pow( (-2.d0*zp/n),iq)/fac4/fac5/fac6;
      sump1=0.d0 ;
      sump2=0.d0 ;
      sump=0.d0 ;
         maxip=n0-l0-1;
      for(/*L*/ ip=0; ip<=maxip; ip++) {
         n40=n0-l0-1-ip;
         n50=ip;
         n60=2*l0+1+ip;
      fac40=fac[n40];
      fac50=fac[n50];
      fac60=fac[n60];
      xp =pow( (-2*zp/n0),ip)/fac40/fac50/fac60;
//--------------------------------------------------
      nn3 = 3+l+l0+la+iq+ip;
      nn2 = 2+l+l0+la+iq+ip;
      fn2 =fac[nn2];
      xmu=  2.5d0+l+l0+iq+ip;
      xnu=  la+0.5d0;
      xc=la+1.5d0;
      xxa=(xnu+xmu)/2;
      xxb=(xnu-xmu)/2+0.5d0;
      xxz=gama/(gama+pow(znz,2));
      call hygfx[xxa][xxb][xc][xxz][hyp31];
      hyp31=hyp31/pow(sqrtl(gama+pow(znz,2)),nn3);
       if(la == l+1 || la == abs(l-1)){
      sump1 = sump1 + xp*fn2*hyp31;
      }
 else {
       xxa=xxa-0.5d0;
       xxb=xxb+0.5d0 ;
      call hygfx[xxa][xxb][xc][xxz][hyp32];
      hyp32=hyp32/pow(sqrtl(gama+pow(znz,2)),nn2);
      sump  = sump  + xp*fn2*( ip*hyp32/nn2 - zp/n0*hyp31 );
        if(la == l){
      sump2 = sump2 + xp*fn2*( (ip+3)*hyp32/nn2 - zp/n0*hyp31 );
       }
      }
       } for(/*L*/ ) {
       if(la == l+1 || la == abs(l-1)){
      sumq1 = sumq1 + xq*sump1;
      }
 else {
      sumq  = sumq + xq*sump;
        if(la == l){
      sumq2=sumq2+xq*sump2;
       }
      }
      } for(/*L*/ ) {
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//                           s=-1 (indica el conjugado)
      s=-1.d0;
        ga=sqrtl(pi);
        for(/*L*/ j=0; j<=la; j++) {
        ga=ga*(0.5d0+j);
        } for(/*L*/ ) {
      xkz=-epsi;
      tk=dacos[xkz/xk];
      calga=calfa/gama;
      xla=sqrtl(pi)*pow(xk,la)/pow(2.d0,(la+1))/ga
                           *sqrtl((2*la+1.d0)/(2*l+1.d0)/4.d0/pi);
       if(la == l+1 || la == abs(l-1)){
         call cg[la][1][l][0][0][0][xc1];
         call cg[la][1][l][m][0][m][xc2];
         call armonic[s][la][m][tk][0.d0][cylam1] ;
      csum1=csum1 +pow( ci,la)*xla*sqrtl(3.d0)*sumq1*xc1*xc2*cylam1;
      }
 else {
         call cg[la][2][l][0][0][0][xc3];
         call cg[la][2][l][m][0][m][xc4];
         call cg[la][2][l][m+1][-1][m][xc5];
         call cg[la][2][l][m-1][1][m][xc6];
         call armonic[s][la][m][tk][0.d0][cylam];
         call armonic[s][la][m+1][tk][0.d0][cylam2];
         call armonic[s][la][m-1][tk][0.d0][cylam3];
      csum2=csum2 +pow( ci,la)*xla*sumq*xc3
                                *(xc5*cylam2-xc6*cylam3)/sqrtl(2.d0);
      csum3=csum3 +pow( ci,la)*xla*sumq* xc3*xc4*cylam*2.d0/sqrtl(3.d0);
         if(la == l){
      csum4=pow(ci,l)*xla*sumq2/sqrtl(3.d0)*cylam;
        }
      }
      } for(/*L*/ ) {
      csum1=csum1*calfa*(ca2/cinu*chyp1+chyp2);
      csum2=csum2*ci*eta*calga*ca2/cinu*chyp1;
      csum3=csum3*(2*v*chyp2-ci*epsi*calga*ca2/cinu*chyp1);
      csum4=csum4*(2*v*chyp2-ci*epsi*calga*ca2/cinu*chyp1);
//-----------------------------------------------------------------------
      csum= 4*pi*(  csum1 - 2 *(csum2+csum3+csum4) );
//-----------------------------------------------------------------------
      cf= ci/(2.d0*pi*v) * con0 * csum ;
//     *          * cdexp(2.d0*cinu*dlog(gama))      cte de modulo=1
      xf=eta * cdabs[cf]*cdabs[cf];
//............................... 
//     SC and ASC corrections
       if (icor == 1) {
      sq=pow((1.d0+pow((sc/xk),2)),(-2))+scoa*(1.d0-pow((1.d0+pow((xk/sc),2)),(-2)))/zt;
      xf=xf*sq;
      }
//...............................
      return;
      }
//***********************************************************************

//     SE    np1-nlm.for

//***********************************************************************
//.....function=eta * /R(eta)/
          2      en la aprox. SE 
//.....desde np1 a nlm
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void rnp1(eta,xf)  {
      // implicit double*8[a-b][d-h][o-z],complex*16[c];
      common/m1/zp,zt;
      common/m2/epsi,v,pi;
      common/m3/ca1,ca2,ca3,ca4,calfa;
      common/m4/ci,cinu,cuno,cdos,con0;
      common/m5/n0,l0,m0,n,l,m;
      common/asc/sc,scoa,icor;
//-----------------------------------------------------------------------
      cero=(0.d0,0.d0);
      gama=pow(eta,2)+pow(epsi,2);
      xk=sqrtl(gama);
      znz=zp/n+zp/n0;
//-----
      vari=-pow(eta,2)/pow(epsi,2);
      call cf21d[cinu][cinu][cdos][vari][chyp1];
      call cf21d[ca1][ca1][cdos][vari][chyp2];
      chyp1=(ca2*chyp1-cinu*(1.d0-vari)*chyp2)/ca4;
//-----------------------------------------------------------------------
      csum1=cero;
      csum2=cero;
      csum3=cero;
        minla2=abs(l-2);
        minla1=abs(l-1);
      minla=min[minla1][minla2];
         if(l == 0)minla=0;
        maxla=l+2;
//-----
      for(/*L*/ la=minla; la<=maxla; la++) {
      sumq1=0.d0 ;
      sumq2=0.d0 ;
      sumq=0.d0 ;
         maxiq=n-l-1;
      for(/*L*/ iq=0; iq<=maxiq; iq++) {
         n4=n-l-1-iq;
         n5=iq;
         n6=2*l+1+iq;
      fac4=fac[n4];
      fac5=fac[n5];
      fac6=fac[n6];
      xq=pow( (-2.d0*zp/n),iq)/fac4/fac5/fac6;

      sump1=0.d0 ;
      sump2=0.d0 ;
      sump=0.d0 ;

         maxip=n0-l0-1;
      for(/*L*/ ip=0; ip<=maxip; ip++) {
         n40=n0-l0-1-ip;
         n50=ip;
         n60=2*l0+1+ip;
      fac40=fac[n40];
      fac50=fac[n50];
      fac60=fac[n60];
      xp =pow( (-2*zp/n0),ip)/fac40/fac50/fac60;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          nn3 = 3+l+l0+la+iq+ip;
          nn2 = 2+l+l0+la+iq+ip;
          fn2 =fac[nn2];
         xmu= 2.5d0+l+l0+iq+ip;
         xnu= la+0.5d0;
         xa= (xnu+xmu)/2;
         xb= xa+0.5d0;
         xc= xnu+1.d0;
         xxb= xc-xb;
         xxz= gama/(gama+pow(znz,2));
      call hygfx[xa][xxb][xc][xxz][hyp31];
      hyp31=hyp31/pow(sqrtl(gama+pow(znz,2)),nn3);
       if(la == l+1 || la == abs(l-1)){
      sump1 = sump1 + xp*fn2*hyp31;
      }
 else {
      xxa=xa-0.5d0;
      xxb=xc-xa ;
      call hygfx[xxa][xxb][xc][xxz][hyp32] ;
      hyp32=hyp32/pow(sqrtl(gama+pow(znz,2)),nn2);
      sump = sump + xp*fn2*( ip*hyp32/nn2 - zp/n0*hyp31 );
         if(la == l){
      sump2 = sump2 + xp*fn2*((ip+3)*hyp32/nn2 - zp/n0*hyp31 );
        }
      }
      } for(/*L*/ ) {

       if(la == l+1 || la == abs(l-1)){
      sumq1 = sumq1 + xq*sump1;
      }
 else {
      sumq  = sumq + xq*sump;
           if(la == l){
      sumq2=sumq2+xq*sump2;
          }
      }
      } for(/*L*/ ) {
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//                           s=-1 (indica el conjugado)
      s=-1.d0;
        ga=sqrtl(pi);
        for(/*L*/ j=0; j<=la; j++) {
        ga=ga*(0.5d0+j);
        } for(/*L*/ ) {
      xkz=-epsi;
      tk=dacos[xkz/xk];
      calga=calfa/gama;
      xla=sqrtl(pi)*pow(xk,la)/pow(2.d0,(la+1))/ga
                     *sqrtl((2*la+1.d0)/(2*l+1.d0)/4.d0/pi);
       if(la == l+1 || la == abs(l-1)){
         call cg[la][1][l][0][0][0][xc1];
         call cg[la][1][l][m-1][1][m][xc2];
         call armonic[s][la][m-1][tk][0.d0][cylam1] ;
      csum1=csum1 +pow( ci,la)*xla*sqrtl(3.d0)*sumq1*xc1*xc2*cylam1;
      }
 else {
         call cg[la][2][l][0][0][0][xc3];
         call cg[la][2][l][m][0][m][xc4];
         call cg[la][2][l][m-2][2][m][xc5];
         call cg[la][2][l][m-1][1][m][xc6];
         call armonic[s][la][m][tk][0.d0][cylam];
         call armonic[s][la][m-2][tk][0.d0][cylam2];
         call armonic[s][la][m-1][tk][0.d0][cylam3];
      csum2=csum2 +pow( ci,la)*xla*sumq*xc3/sqrtl(6.d0)
                         *(xc4*cylam-sqrtl(6.d0)*xc5*cylam2);
      csum3=csum3 +pow( ci,la)*xla*sumq*xc3*xc6*cylam3;
             if(la == l){
      csum4=-pow(ci,l)*xla*sumq2/sqrtl(6.d0)*cylam;
            }
      }
         } for(/*L*/ ) {
      csum1=csum1*calfa*(ca2/cinu*chyp1+chyp2);
      csum2=csum2*ci*eta*calga*ca2/cinu*chyp1;
      csum3=csum3*(2*v*chyp2-ci*epsi*calga*ca2/cinu*chyp1);
      csum4=csum4*ci*eta*calga*ca2/cinu*chyp1;
//-----------------------------------------------------------------------
      csum= 4*pi* ( csum1 - 2 *(csum2+csum3+csum4)  );
//-----------------------------------------------------------------------
      cf= ci/(2.d0*pi*v) * con0 * csum;
//     *           * cdexp(2.d0*cinu*dlog(gama))   cte. de modulo=1
      xf=eta*cdabs[cf]*cdabs[cf];
//............................... 
//     SC and ASC corrections
       if (icor == 1) {
      sq=pow((1.d0+pow((sc/xk),2)),(-2))+scoa*(1.d0-pow((1.d0+pow((xk/sc),2)),(-2)))/zt;
      xf=xf*sq;
      }
//...............................
      return;
      }
//***********************************************************************

//     SE    np1m-nlm.for

//***********************************************************************
//.....function=eta * /R(eta)/
          2      en la aprox. SE 
//.....desde np-1 a nlm 
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void rnp1m(eta,xf)  {
      // implicit double*8[a-b][d-h][o-z],complex*16[c];
      common/m1/zp,zt;
      common/m2/epsi,v,pi;
      common/m3/ca1,ca2,ca3,ca4,calfa;
      common/m4/ci,cinu,cuno,cdos,con0;
      common/m5/n0,l0,m0,n,l,m;
      common/asc/sc,scoa,icor;
//-----------------------------------------------------------------------
      cero=(0.d0,0.d0);
      gama=pow(eta,2)+pow(epsi,2);
      xk=sqrtl(gama);
      znz=zp/n+zp/n0;
//-----
      vari=-pow(eta,2)/pow(epsi,2);
      call cf21d[cinu][cinu][cdos][vari][chyp1];
      call cf21d[ca1][ca1][cdos][vari][chyp2];
      chyp1=(ca2*chyp1-cinu*(1.d0-vari)*chyp2)/ca4;
//-----------------------------------------------------------------------
      csum1=cero;
      csum2=cero;
      csum3=cero;
        minla2=abs(l-2);
        minla1=abs(l-1);
      minla=min[minla1][minla2];
         if(l == 0)minla=0;
        maxla=l+2;
//-----
      for(/*L*/ la=minla; la<=maxla; la++) {

      sumq1=0.d0 ;
      sumq2=0.d0 ;
      sumq=0.d0 ;

         maxiq=n-l-1;
      for(/*L*/ iq=0; iq<=maxiq; iq++) {
         n4=n-l-1-iq;
         n5=iq;
         n6=2*l+1+iq;
      fac4=fac[n4];
      fac5=fac[n5];
      fac6=fac[n6];
      xq=pow( (-2.d0*zp/n),iq)/fac4/fac5/fac6;

      sump1=0.d0 ;
      sump2=0.d0 ;
      sump=0.d0 ;

         maxip=n0-l0-1;
      for(/*L*/ ip=0; ip<=maxip; ip++) {
         n40=n0-l0-1-ip;
         n50=ip;
         n60=2*l0+1+ip;
      fac40=fac[n40];
      fac50=fac[n50];
      fac60=fac[n60];
      xp =pow( (-2*zp/n0),ip)/fac40/fac50/fac60;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          nn3 = 3+l+l0+la+iq+ip;
          nn2 = 2+l+l0+la+iq+ip;
          fn2 =fac[nn2];
         xmu= 2.5d0+l+l0+iq+ip;
         xnu= la+0.5d0;
         xa= (xnu+xmu)/2;
         xb= xa+0.5d0;
         xc= xnu+1.d0;
         xxb= xc-xb;
         xxz= gama/(gama+pow(znz,2));
      call hygfx[xa][xxb][xc][xxz][hyp31];
      hyp31=hyp31/pow(sqrtl(gama+pow(znz,2)),nn3);
       if(la == l+1 || la == abs(l-1)){
      sump1 = sump1 + xp*fn2*hyp31;
      }
 else {
      xxa=xa-0.5d0;
      xxb=xc-xa ;
      call hygfx[xxa][xxb][xc][xxz][hyp32] ;
      hyp32=hyp32/pow(sqrtl(gama+pow(znz,2)),nn2);
      sump = sump + xp*fn2*( ip*hyp32/nn2 - zp/n0*hyp31 );
         if(la == l){
      sump2 = sump2 + xp*fn2*((ip+3)*hyp32/nn2 - zp/n0*hyp31 );
        }
      }
      } for(/*L*/ ) {

       if(la == l+1 || la == abs(l-1)){
      sumq1 = sumq1 + xq*sump1;
      }
 else {
      sumq  = sumq + xq*sump;
           if(la == l){
      sumq2=sumq2+xq*sump2;
          }
      }
      } for(/*L*/ ) {
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//                           s=-1 (indica el conjugado)
      s=-1.d0;
        ga=sqrtl(pi);
        for(/*L*/ j=0; j<=la; j++) {
        ga=ga*(0.5d0+j);
        } for(/*L*/ ) {
      xkz=-epsi;
      tk=dacos[xkz/xk];
      calga=calfa/gama;
      xla=sqrtl(pi)*pow(xk,la)/pow(2.d0,(la+1))/ga
                     *sqrtl((2*la+1.d0)/(2*l+1.d0)/4.d0/pi);
       if(la == l+1 || la == abs(l-1)){
         call cg[la][1][l][0][0][0][xc1];
         call cg[la][1][l][m+1][-1][m][xc2];
         call armonic[s][la][m+1][tk][0.d0][cylam1] ;
      csum1=csum1 +pow( ci,la)*xla*sqrtl(3.d0)*sumq1*xc1*xc2*cylam1;
      }
 else {
         call cg[la][2][l][0][0][0][xc3];
         call cg[la][2][l][m][0][m][xc4];
         call cg[la][2][l][m+2][-2][m][xc5];
         call cg[la][2][l][m+1][-1][m][xc6];
         call armonic[s][la][m][tk][0.d0][cylam];
         call armonic[s][la][m+2][tk][0.d0][cylam2];
         call armonic[s][la][m+1][tk][0.d0][cylam3];
      csum2=csum2 +pow( ci,la)*xla*sumq*xc3
                         *(-xc4*cylam/sqrtl(6.d0)+xc5*cylam2);
      csum3=csum3 +pow( ci,la)*xla*sumq*xc3*xc6*cylam3;
          if(la == l){
      csum4=pow(ci,l)*xla*sumq2/sqrtl(6.d0)*cylam;
         }
      }
         } for(/*L*/ ) {
      csum1=csum1*calfa*(ca2/cinu*chyp1+chyp2);
      csum2=csum2*ci*eta*calga*ca2/cinu*chyp1;
      csum3=csum3*(2*v*chyp2-ci*epsi*calga*ca2/cinu*chyp1);
      csum4=csum4*ci*eta*calga*ca2/cinu*chyp1;
//-----------------------------------------------------------------------
      csum= 4*pi*( csum1 - 2 *(csum2+csum3+csum4)  );
//-----------------------------------------------------------------------
      cf= ci/(2.d0*pi*v) * con0 * csum;
//     *   * cdexp(2.d0*cinu*dlog(gama))  factor de módulo=1
      xf=eta* cdabs[cf]*cdabs[cf];
//............................... 
//     SC and ASC corrections
       if (icor == 1) {
      sq=pow((1.d0+pow((sc/xk),2)),(-2))+scoa*(1.d0-pow((1.d0+pow((xk/sc),2)),(-2)))/zt;
      xf=xf*sq;
      }
//...............................
      return;
      }
//***********************************************************************

//     Subrutinas.for

//***********************************************************************
//     Armonicos esfericos 
//     Autor: Cesar Ramirez          fecha: 10-4-2000
//     l,/m/:   numeros enteros
//     tk     angulo teta
//     s  +1: real  -1: imaginario
//     ma  valor absoluto de m
//-----------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void armonic(s,l,m,tk,fik,cylm) {
//     ------------------------------------------------------------------
      // implicit double*8[a-b][d-h][o-z],complex*16[c];
        ci=(0.d0,1.d0);
        pi=dacos[-1.d0];
        ma=abs(m);
      cylm=(0.d0,0.d0);
       if(ma > l)return;
        n1=l-ma;
        n2=l+ma;
        fac11=fac[n1];
        fac12=fac[n2];
      x1= sqrtl((2.d0*l+1.d0)/(4.d0*pi)*fac11/fac12);
      x=dcos[tk];
        xlm=plgndr[l][ma][x];
       if(m >= 0){
       cylm=x1*cdexp[s*ci*ma*fik]*xlm;
      }
 else {
       cylm=pow((-1.d0),ma)*x1*cdexp[-s*ci*ma*fik]*xlm;
      }
        return;
      }
//***********************************************************************
//     Legendre polynomials (Numerical Recipes Software: www.nr.com)
//     x=cos(t)
//     ------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function PLGNDR(L,M,X) {
      double*8 x,plgndr;
       if(M < 0 || M > L || dABS[X] > 1.){
      printf(*,*)' m<0,   or   m>l,   or   dabs[x]>1';
      pause 'bad arguments';
      }
      PMM=1.;
       if(M > 0) {
        SOMX2=sqrt((1.-X)*(1.+X));
        FACT=1.;
        for(/*L11*/ I=1; I<=M; I++) {
          PMM=-PMM*FACT*SOMX2;
          FACT=FACT+2.;
L11:    } //continue;
      }
       if(L == M) {
        PLGNDR=PMM;
      }
 else {
        PMMP1=X*(2*M+1)*PMM;
         if(L == M+1) {
          PLGNDR=PMMP1;
        }
 else {
          for(/*L12*/ LL=M+2; LL<=L; LL++) {
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M);
            PMM=PMMP1;
            PMMP1=PLL;
L12:      } //continue;
          PLGNDR=PLL;
        }
      }
//      write(*,*)x,plgndr
      return;
      }
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void gamac(za,ny,test1,test2,nt,logam)  {
//     ******************************************************************
//          auteur: r. gayet              date inconnue 
//          version 2    (modif a. salin)      30/10/81 
//     log de gamma(z) ecrit si ny=2 
//     test1 est la difference minimum admise entre reel(z) et un pole 
//     pour distinguer reel(z) du pole 
//     test2 est la valeur minimum admise pour imag(z) quand la 
//     difference entre reel(z) et un pole est inferieure a test1 
//     nt=2 si z est un pole de gamma:dans ce cas logam=0 arbitrairement-
//     sinon nt=1 
      complex*16 z,z2,logam,c,za ;
      double*8 a[2],pi,d,test1,test2,test10 ;
      equivalence[z][a[1]] ;
      double /*data*/ pi/.918938533204673d0/ ;
      z=za ;
      logam=(0.d0,0.d0) ;
      irz=idint[a[1]] ;
      test10=test1*2.d0 ;
       if(a[1]-test10)2,2,1 ;
L2:   d=dabs[a[1]-dfloat[irz]] ;
       if(d > test1 ) goto L1 ;
       if(dabs[a[2]] > test2)  goto L1 ;
      printf(6,3) d,z ;
//L3:   format(//,4x,'reel(z)=',1pd9.2,'+',1pd22.15,/,4x,'imag(z)=',1pd22.
          15,//,4x,'z est considere comme un pole de la fonction gamma') ;
      nt=2 ;
      goto L100 ;
L1:   } //continue ;
      nt=1 ;
      ng=10-irz ;
       if(ng)4,4,5 ;
L5:   c=dcmplx[dfloat[ng]][0.d0] ;
      z=z+c ;
L4:   z2=z*z ;
      d=cdabs[z] ;
       if(d >= 1.d+02) goto L6 ;
      logam=1.d0/156.d0/z2-691.d0/360360.d0 ;
      logam=logam/z2+1.d0/1188.d0 ;
      logam=logam/z2-1.d0/1680.d0 ;
L6:    if(d >= 1.d+04)goto L7 ;
      logam=logam/z2+1.d0/1260.d0 ;
      logam=logam/z2-1.d0/360.d0 ;
L7:    if(d >= 1.d+07) goto L8 ;
      logam=(logam/z2+1.d0/12.d0)/z ;
L8:   logam=logam+pi-z+(z-0.5d0)*cdlog[z] ;
       if(ng)100,100,9 ;
L9:   c=(1.d0,0.d0) ;
      for(/*L10*/ i=1; i<=ng; i++) {
      z=z-c ;
      logam=logam-cdlog[z] ;
L10:  } //continue ;
L100: } //continue ;
       if(ny == 2) printf(ny,200) z,logam;
//L200: format(//,2x,'log de gamma(',2(1pd22.15),') = ',2(1pd22.15)) ;
      return ;
      } ;
//***********************************************************************
//***********************************************************************
//  Rutina para calcular el factorial de n
//  n! = n(n-1)(n-2)...2 ,
//  0! = 1
//----------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
function fac(n) {
      // implicit none;
      double*8 xn,xfac,fac,pi;
      int n,j;

      pi=dacos[-1.d0];
      xn=n*1.d0;
       if(n == 0 || n == 1){
          fac=1.d0;
      }
 else {
           if(n > 90){
          fac=sqrtl(pi*(2*xn+1/3.d0))*pow(xn,xn)/expl(xn);
          }
 else {
          xfac=0.d0;
          for(/*L*/ j=2; j<=n; j++) {
          xfac=xfac+logl(dfloat[j]);
          } for(/*L*/ ) {
          fac=expl(xfac);
          }
      }
      return;
      }
//*********************************************************************
//**************************** l o g f a c ****************************** 
//     calculates ln((n-1)!).
//     beware: fac(n) = ln((n-1)!) = ln(gamma(n))
//-----------------------------------------------------------------------
      block double /*data*/ faclog;
      parameter [idf=140] ;
      // implicit double*8 [a-h][o-z];
      logical init ;
      common/fact/fac[idf],init ;
      double /*data*/ init/.false./ ;
      } ;
//-----------------------------------------------------------------------
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void logfac  {
      parameter [idf=140] ;
      // implicit double*8 [a-h][o-z];
      logical init ;
      common/fact/fac[idf],init ;
      init=.true. ;
      fac[1]=0.d0 ;
      for(/*L*/ i=2; i<=idf; i++) {
        fac[i]=fac[i-1]+logl(dfloat[i-1]) ;
      } for(/*L*/ ) {
      } ;
//*********************************************************************
//.......... cg Calcula los coeficientes de Clebsch-Gordan
//.......... Debe llamarse antes al logfac.
//.......... logcle:  simbolos 3j
//.......... xcg:     coef de C-G
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void cg(l1,l2,l3,m1,m2,m3,z)  {
      // implicit double*8 [a-h][o-z];
      call logfac;
      call logcle[l1][l2][l3][m1][m2][-m3][x3j] ;
      z=pow((-1.d0),abs(l1-l2-m3))*sqrtl(2.d0*l3+1.d0)*x3j;
      }
//**************************** l o g c l e ****************************** 
//         author:  a.salin     version 2  23/2/77 - 5/2/96

//     calculation of WIGNER's 3j  (definition of MESSIAH)
//     restriction: INTEGER MOMENTS ONLY.

//     before the first run of logcle, a call to logfac should be
//     performed. logfac calculates ln((n-1)!) for n going from 1 to idf
//     and stores the results in the common fact. 
//     for large values of the l's, increase idf. 
//        input 
//                l1, l2, l3, m1, m2, m3: momenta and their components. 
//        output
//                z: 3j coefficient.
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void logcle(l1,l2,l3,m1,m2,m3,z) {
      // implicit double*8 [a-h][o-z];
      parameter[idf=140] ;
      logical init ;
      double ac[idf] ;
      common/fact/fac[idf],init ;
// 
       if(.not.init) { 
        printf(*,98) ;
//L98:    format(1x,' logfac not called before logcle') ;
        return ;
      }

      i4=l1+l2+l3+2 ;
       if(i4 > idf) { 
        printf(*,99) ;
//L99:    format(1x,'logcle: value of idf too small') ;
        return;
      }

       if(m1+m2+m3 != 0) { 
        z=0.d0 ;
        return ;
      } ;
// 
      izmax=min[l1+l2-l3][l1-m1][l2+m2]+1 ;
      izmin=max[0][l2-l3-m1][l1+m2-l3]+1 ;
       if(izmax < izmin) { 
        z=0.d0 ;
        return ;
      } ;

      i1=l1+l2-l3+1 ;
      i2=l1-m1+1 ;
      i3=l2+m2+1 ;
      abra=0.5d0*(fac[i1]+fac[l3+l1-l2+1]+fac[l3+l2-l1+1]-fac[i4]
               +fac[l1+m1+1]+fac[i2]+fac[i3]+fac[l2-m2+1]+fac[l3+m3+1]
               +fac[l3-m3+1]);
      k1=l3-l2+m1+1;
      k2=l3-l1-m2+1 ;
      gros=250.d0 ;
      for(/*L8*/ ii=izmin; ii<=izmax; ii++) {
        i=ii-1 ;
        ac[ii]=fac[i+1]+fac[i1-i]+fac[i2-i]
                   +fac[i3-i]+fac[k1+i]+fac[k2+i];
         if(ac[ii] < gros)  gros=ac[ii];
L8:   } //continue ;
      accu=0.d0;
      sig=pow((-1.d0),izmin );
      for(/*L9*/ ii=izmin; ii<=izmax; ii++) {
        sig=-sig ;
        ac[ii]=ac[ii]-gros ;
        accu=accu+sig*expl(-ac[ii]) ;
L9:   } //continue ;
      z=pow((-1.d0),abs(l1-l2-m3))*expl(abra-gros)*accu ;
      return ;
      } ;
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void cf21d[ca][cb][cc][x][cf]
          ********************************
          gram: hyper.f
          sion: 1.6 [see readme_hyper] 24/03/1997
          hor: Pablof [pablof@cab.cnea.edu.ar]
          ergeometric function for CDW-EIS calculations. 
          ble precision, see comment for V 1.4; {
*  We assume that x is double
          ~~~~~~~~~~
          ********************************;
      // implicit none;
      complex*16 ca,cb,cc,cf;
      double*8 x;
      complex*16 calgam,c0,c1;
      parameter[c0=(0.d0][0.d0], c1=(1.d0,0.d0));

       if(ca == c0 || cb == c0) {
        cf=c1;
        return;
      }
       if(x == 0.d0) {
        cf=c1;
        return;
      }
 else if(x == 1.d0) {
        cf=cdexp[calgam[cc]-calgam[cc-ca-cb]
              -calgam[cc-ca]-calgam[cc-cb]];
        return;
      }
 else {
        call hyper[ca][cb][cc][x][cf];
      }
      }
          ****;
//       ====================================================
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  void HYGFX[A][B][C][X][HF] {
//       ====================================================
//       Purpose: Compute hypergeometric function F(a,b,c,x)
//       Input :  a --- Parameter 
//                b --- Parameter
//                c --- Parameter, c <> 0,-1,-2,...
//                x --- Argument   ( x < 1 )
//       Output:  HF --- F(a,b,c,x)
//       Routines called:
//            (1) GAMMA for computing gamma function
//            (2) PSI for computing psi function
//       ====================================================
        // implicit DOUBLE PRECISION [A-H][O-Z];
        LOGICAL L0,L1,L2,L3,L4,L5;
        PI=3.141592653589793D0;
        EL=.5772156649015329D0;
        L0=C == INT[C] && C < 0.0;
        L1=1.0D0-X < 1.0D-15 && C-A-B <= 0.0;
        L2=A == INT[A] && A < 0.0;
        L3=B == INT[B] && B < 0.0;
        L4=C-A == INT[C-A] && C-A <= 0.0;
        L5=C-B == INT[C-B] && C-B <= 0.0;
         if (L0 || L1) {
           printf(*,*)'The hypergeometric series is divergent';
           return;
        }
        EPS=1.0D-15;
         if (X > 0.95) EPS=1.0D-8;
         if (X == 0.0 || A == 0.0 || B == 0.0) {
           HF=1.0D0;
           return;
        }
 else if (1.0D0-X == EPS && C-A-B > 0.0) {
           CALL GAMMA[C][GC];
           CALL GAMMA[C-A-B][GCAB];
           CALL GAMMA[C-A][GCA];
           CALL GAMMA[C-B][GCB];
           HF=GC*GCAB/(GCA*GCB);
           return;
        }
 else if (1.0D0+X <= EPS && DABS[C-A+B-1.0] <= EPS) {
           G0=sqrtl(PI)*pow(2.0D0,(-A));
           CALL GAMMA[C][G1];
           CALL GAMMA[1.0D0+A/2.0-B][G2];
           CALL GAMMA[0.5D0+0.5*A][G3];
           HF=G0*G1/(G2*G3);
           return;
        }
 else if (L2 || L3) {
            if (L2) NM=INT[fabs(A)];
            if (L3) NM=INT[fabs(B)];
           HF=1.0D0;
           R=1.0D0;
           for(/*L10*/ K=1; K<=NM; K++) {
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X;
L10:          HF=HF+R;
           return;
        }
 else if (L4 || L5) {
            if (L4) NM=INT[fabs(C-A)];
            if (L5) NM=INT[fabs(C-B)];
           HF=1.0D0;
           R=1.0D0;
           for(/*L15*/ K=1; K<=NM; K++) {
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X;
L15:          HF=HF+R;
           HF=pow((1.0D0-X),(C-A-B))*HF;
           return;
        }
        AA=A;
        BB=B;
        X1=X;
         if (X < 0.0D0) {
           X=X/(X-1.0D0);
            if (C > A && B < A && B > 0.0) {
              A=BB;
              B=AA;
           }
           B=C-B;
        }
         if (X >= 0.75D0) {
           GM=0.0D0;
            if (DABS[C-A-B-INT[C-A-B]] < 1.0D-15) {
              M=INT[C-A-B];
              CALL GAMMA[A][GA];
              CALL GAMMA[B][GB];
              CALL GAMMA[C][GC];
              CALL GAMMA[A+M][GAM];
              CALL GAMMA[B+M][GBM];
              CALL PSI[A][PA];
              CALL PSI[B][PB];
               if (M != 0) GM=1.0D0;
              for(/*L30*/ J=1; J<=fabs[M; J++)-1) {
L30:             GM=GM*J;
              RM=1.0D0;
              for(/*L35*/ J=1; J<=fabs[M; J++)) {
L35:             RM=RM*J;
              F0=1.0D0;
              R0=1.0D0;
              R1=1.0D0;
              SP0=0.D0;
              SP=0.0D0;
               if (M >= 0) {
                 C0=GM*GC/(GAM*GBM);
                 C1=-GC*pow((X-1.0D0),M)/(GA*GB*RM);
                 for(/*L40*/ K=1; K<=M-1; K++) {
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X);
L40:                F0=F0+R0;
                 for(/*L45*/ K=1; K<=M; K++) {
L45:                SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K;
                 F1=PA+PB+SP0+2.0D0*EL+logl(1.0D0-X);
                 for(/*L55*/ K=1; K<=250; K++) {
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0));
                    SM=0.0D0;
                    for(/*L50*/ J=1; J<=M; J++) {
L50:                   SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
                              [B+J+K-1.0];
                    RP=PA+PB+2.0D0*EL+SP+SM+logl(1.0D0-X);
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X);
                    F1=F1+R1*RP;
                     if (DABS[F1-HW] < DABS[F1]*EPS) goto L60;
L55:                HW=F1;
L60:             HF=F0*C0+F1*C1;
              }
 else if (M < 0) {
                 M=-M;
                 C0=GM*GC/(GA*GB*pow((1.0D0-X),M));
                 C1=-pow((-1),M)*GC/(GAM*GBM*RM);
                 for(/*L65*/ K=1; K<=M-1; K++) {
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X);
L65:                F0=F0+R0;
                 for(/*L70*/ K=1; K<=M; K++) {
L70:                SP0=SP0+1.0D0/K;
                 F1=PA+PB-SP0+2.0D0*EL+logl(1.0D0-X);
                 for(/*L80*/ K=1; K<=250; K++) {
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0));
                    SM=0.0D0;
                    for(/*L75*/ J=1; J<=M; J++) {
L75:                   SM=SM+1.0D0/(J+K);
                    RP=PA+PB+2.0D0*EL+SP-SM+logl(1.0D0-X);
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X);
                    F1=F1+R1*RP;
                     if (DABS[F1-HW] < DABS[F1]*EPS) goto L85;
L80:                HW=F1;
L85:             HF=F0*C0+F1*C1;
              }
           }
 else {
              CALL GAMMA[A][GA];
              CALL GAMMA[B][GB];
              CALL GAMMA[C][GC];
              CALL GAMMA[C-A][GCA];
              CALL GAMMA[C-B][GCB];
              CALL GAMMA[C-A-B][GCAB];
              CALL GAMMA[A+B-C][GABC];
              C0=GC*GCAB/(GCA*GCB);
              C1=GC*GABC/(GA*GB)*pow((1.0D0-X),(C-A-B));
              HF=0.0D0;
              R0=C0;
              R1=C1;
              for(/*L90*/ K=1; K<=250; K++) {
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X);
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
                        *(1.0-X);
                 HF=HF+R0+R1;
                  if (DABS[HF-HW] < DABS[HF]*EPS) goto L95;
L90:             HW=HF;
L95:          HF=HF+C0+C1;
           }
        }
 else {
           A0=1.0D0;
            if (C > A && C < 2.0D0*A && 
                   C > B && C < 2.0D0*B) {
              A0=pow((1.0D0-X),(C-A-B));
              A=C-A;
              B=C-B;
           }
           HF=1.0D0;
           R=1.0D0;
           for(/*L100*/ K=1; K<=250; K++) {
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X;
              HF=HF+R;
               if (DABS[HF-HW] <= DABS[HF]*EPS) goto L105;
L100:         HW=HF;
L105:      HF=A0*HF;
        }
         if (X1 < 0.0D0) {
           X=X1;
           C0=1.0D0/pow((1.0D0-X),AA);
           HF=C0*HF;
        }
        A=AA;
        B=BB;
         if (K > 120) printf(*,115);
//L115:   format(1X,'Warning! You should check the accuracy');
        return;
        }


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  void GAMMA[X][GA] {

//       ==================================================
//       Purpose: Compute gamma function â(x)
//       Input :  x  --- Argument of â(x)
//                       ( x is not equal to 0,-1,-2,úúú)
//       Output:  GA --- â(x)
//       ==================================================

        // implicit DOUBLE PRECISION [A-H][O-Z];
        double G[26];
        PI=3.141592653589793D0;
         if (X == INT[X]) {
            if (X > 0.0D0) {
              GA=1.0D0;
              M1=X-1;
              for(/*L10*/ K=2; K<=M1; K++) {
L10:             GA=GA*K;
           }
 else {
              GA=1.0D+300;
           }
        }
 else {
            if (DABS[X] > 1.0D0) {
              Z=DABS[X];
              M=INT[Z];
              R=1.0D0;
              for(/*L15*/ K=1; K<=M; K++) {
L15:             R=R*(Z-K);
              Z=Z-M;
           }
 else {
              Z=X;
           }
           double /*data*/ G/1.0D0,0.5772156649015329D0,
                    -0.6558780715202538D0, -0.420026350340952D-1,
                    0.1665386113822915D0,-.421977345555443D-1,
                    -.96219715278770D-2, .72189432466630D-2,
                    -.11651675918591D-2, -.2152416741149D-3,
                    .1280502823882D-3, -.201348547807D-4,
                    -.12504934821D-5, .11330272320D-5,
                    -.2056338417D-6, .61160950D-8,
                    .50020075D-8, -.11812746D-8,
                    .1043427D-9, .77823D-11,
                    -.36968D-11, .51D-12,
                    -.206D-13, -.54D-14, .14D-14, .1D-15/;
           GR=G[26];
           for(/*L20*/ K=25; K<=1; K+-1) {
L20:          GR=GR*Z+G[K];
           GA=1.0D0/(GR*Z);
            if (DABS[X] > 1.0D0) {
              GA=GA*R;
               if (X < 0.0D0) GA=-PI/(X*GA*DSIN[PI*X]);
           }
        }
        return;
        }


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  void PSI[X][PS] {

//       ======================================
//       Purpose: Compute Psi function
//       Input :  x  --- Argument of psi(x)
//       Output:  PS --- psi(x)
//       ======================================

        // implicit DOUBLE PRECISION [A-H][O-Z];
        XA=DABS[X];
        PI=3.141592653589793D0;
        EL=.5772156649015329D0;
        S=0.0D0;
         if (X == INT[X] && X <= 0.0) {
           PS=1.0D+300;
           return;
        }
 else if (XA == INT[XA]) {
           N=XA;
           for(/*L10*/ K=1; K<=N-1; K++) {
L10:          S=S+1.0D0/K;
           PS=-EL+S;
        }
 else if (XA+.5 == INT[XA+.5]) {
           N=XA-.5;
           for(/*L20*/ K=1; K<=N; K++) {
L20:          S=S+1.0/(2.0D0*K-1.0D0);
           PS=-EL+2.0D0*S-1.386294361119891D0;
        }
 else {
            if (XA < 10.0) {
              N=10-INT[XA];
              for(/*L30*/ K=0; K<=N-1; K++) {
L30:             S=S+1.0D0/(XA+K);
              XA=XA+N;
           }
           X2=1.0D0/(XA*XA);
           A1=-.8333333333333D-01;
           A2=.83333333333333333D-02;
           A3=-.39682539682539683D-02;
           A4=.41666666666666667D-02;
           A5=-.75757575757575758D-02;
           A6=.21092796092796093D-01;
           A7=-.83333333333333333D-01;
           A8=.4432598039215686D0;
           PS=logl(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
                  A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1);
           PS=PS-S;
        }
         if (X < 0.0) PS=PS-PI*DCOS[PI*X]/DSIN[PI*X]-1.0D0/X;
        return;
        }
//************************************************************

//     Integrador.for

//*************************************** fin  programa principal ******
//     INTEGRADOR POR SIMPSO en [ymin, infinito)
//     AJUSTANDO LA PRESICION EN CADA INTERVALO Y LUEGO
//     entre el ultimo intervalo Y LA SUMA TOTAL
//     Para funciones Reales 
//      (10-8-94)                        autor:  Cesar Ramirez

//*****
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void INTDEF(PRE,Ymin,Ymax,cfun,cint) {
      // implicit double*8[A-H][O-Z];
      preci=pre;
      x1=Ymin;
      x2=Ymax;
      csum1=0.d0;
      cint=0.d0;
L39:  n=1;
      csum=0.d0;
      h=x2-x1;
      call cfun[x1][cf1];
      call cfun[x2][cf2];
      cf12=cf1+cf2;
      cfmc=0.d0;
      cfm=0.d0;
L40:  csum1=csum;
      d=h/n;
      csum=(cf12+2.d0*cfmc)*d/6.d0;
      for(/*L50*/ i=1; i<=n; i++) {
      xm=h/n*(i-0.5d0)+x1;
      call cfun[xm][cfm] ;
      csum=csum+4.d0*cfm*d/6.d0 ;
      cfmc=cfmc+cfm;
L50:  } //continue;
       if(dabs[csum-csum1] > preci*dabs[csum1]){
      n=n*2.d0;
      goto L40;
      }
 else {
      x1=x2;
      x2=x1+h;
      cint=cint+csum;
      }
       if(dabs[csum] > preci*dabs[cint]){
      goto L39;
      }
 else {
      } //continue;
      } ;
      return;
      }
//***********************************************************************