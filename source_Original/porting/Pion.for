	subroutine pion(tck,tcl1,tcl2,tcm1,tcm2,tcn)
C------------------------------------------------PION.FOR-----------------
C SECTIONS EFFICACES D'IONISATION PWBA	            Version du 21/04/94
C d'un projectile à rk,rl1,rl2,rm1 et rm2          -----------------------
C électrons K,2s,2p,3s,p et 3d		      DEBUT DU PROGRAMME PRINCIPAL.
C-------------------------------------------------------------------------
	DIMENSION DCKE(80),DCL1E(80),DCL2E(80),DCM1E(80),
     s	DCM2E(80),T(80),WK(80),WL1(80),WL2(80),WM1(80),WM2(80),
     s	sk(80),sl1(80),sl2(80),sm1(80),sm2(80),
     s	DCNE(80),Wn(80),sn(80)
	common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn
	common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn
	common/don/zp,E,zt
	z1=zp
	z2=zt
	epo=E
	rk=(y1s)
	rl1=(y2s)
	rl2=(y2p)
	rm1=(y3s+y3p)
	rm2=(y3d)
	rn=(yn)
	rk0=max(rk,1.)
	rl10=max(rl1,1.)
	rl20=max(rl2,1.)
	rm10=max(rm1,1.)
	rm20=max(rm2,1.)
	rn0=max(rn,1.)
	rkm=rk0-1.
	rl1m=rl10-1.
	rl2m=rl20-1.
	rm1m=rm10-1.
	rm2m=rm20-1.
	rnm=rn0-1.
	Bk=BTot(z1,rk0,rl1,rl2,rm1,rm2,rn)
     s	-BTot(z1,rkm,rl1,rl2,rm1,rm2,rn)
	Bl1=BTot(z1,rk,rl10,rl2,rm1,rm2,rn)
     s 	-BTot(z1,rk,rl1m,rl2,rm1,rm2,rn)
	Bl2=BTot(z1,rk,rl1,rl20,rm1,rm2,rn)
     s	-BTot(z1,rk,rl1,rl2m,rm1,rm2,rn)
	Bm1=BTot(z1,rk,rl1,rl2,rm10,rm2,rn)
     s	-BTot(z1,rk,rl1,rl2,rm1m,rm2,rn)
	Bm2=BTot(z1,rk,rl1,rl2,rm1,rm20,rn)
     s	-BTot(z1,rk,rl1,rl2,rm1,rm2m,rn)
	Bn=BTot(z1,rk,rl1,rl2,rm1,rm2,rn0)
     s	-BTot(z1,rk,rl1,rl2,rm1,rm2,rnm)
	if (bl1.lt.3.4) bl1=3.4
	if (bl2.lt.3.4) bl2=3.4
	if(Bm1.le.1.51) Bm1=1.51
	if(Bm2.le.1.51) Bm2=1.51
	if(Bn.le.0.85) Bn=0.85
c ************** calcul des constantes *******************************
	betal=(1.-(1./(1.+epo/931.5)**2.))*(137.036)**2
	CSTE=3.519D4*z2**2.
c-------- charges effectives ecrantage Slater --------
	if(rk.gt.1.) then
	  zk=z1-0.3
	else
	  zk=z1
	end if
	zl1=z1-(0.8*rk+0.3*max((rl1-1.),0.))
	zl2=z1-(rk+0.75*rl1+0.35*max((rl2-1.),0.))
	zm1=z1-(rk+0.85*(rl1+rl2)+0.35*max((rm1-1.),0.))
	zm2=z1-(rk+rl1+rl2+rm1+0.35*max((rm2-1.),0.))
	zn=z1-(rk+rl1+rl2+rm1+rm2+0.35*max((rn-1.),0.))
	if(zl1.le.1.) zl1=1.
	if(zl2.le.1.) zl2=1.
	if(zm1.le.1.) zm1=1.
	if(zm2.le.1.) zm2=1.
	if(zn.le.1.) zn=1.
c -------------------------------------------------------
c     factors for screening and antiscreening corrections
c........................................................
	vu=137.036*(1.-(1.+epo/931.5)**(-2.))**(0.5)
	zp1=0.9354*zk
	zp21=0.9354*zl1/2.
	zp22=0.9354*zl2/2.
	zp31=0.9354*zm1/3.
	zp32=0.9354*zm2/3.
	zp4=0.9534*zn/4.
	if (vu.gt.zp1) then
	coak=1.-(zp1/vu)**4
	else
	coak=0.
	endif
	if (vu.gt.zp21) then
	coal1=1.-(zp21/vu)**4
	else
	coal1=0.
	endif
	if (vu.gt.zp22) then
	coal2=1.-(zp22/vu)**4
	else
	coal2=0.
	endif
	if (vu.gt.zp31) then
	coam1=1.-(zp31/vu)**4
	else
	coam1=0.
	endif
	if (vu.gt.zp32) then
	coam2=1.-(zp32/vu)**4
	else
	coam2=0.
	endif
	if (vu.gt.zp4) then
	coan=1.-(zp4/vu)**4
	else
	coan=0.
	endif
c ----------------------------------------------------
c	l'ionisation en n=4 est calculée comme étant
c	l'ionisation 2p de z/2
c ----------------------------------------------------
	zn=zn/2.
c
	etak=betal/(zk)**2.
	etal1=betal/(zl1)**2.
	etal2=betal/(zl2)**2.
	etam1=betal/(zm1)**2.
	etam2=betal/(zm2)**2.
	etan=betal/(zn)**2.
	tetak=bk/(13.6058*zk**2.)
	tetal1=4*bl1/(13.6058*zl1**2.)
	tetal2=4*bl2/(13.6058*zl2**2.)
	tetam1=9*bm1/(13.6058*zm1**2.)
	tetam2=9*bm2/(13.6058*zm2**2.)
	tetan=16*bn/(13.6058*(2.*zn)**2.)
c	tetak=1.
c	tetal1=1.
c	tetal2=1.
c	tetam1=1.
c	tetam2=1.
c	tetan=1.
	cstk=cste/(zk**4.)
	cstl1=cste/(zl1**4.)
	cstl2=cste/(zl2**4.)
	cstm1=cste/(zm1**4.)
	cstm2=cste/(zm2**4.)
	cstn=cste/(zn**4.)
c
	fs1=1.
	if (zt.eq.1.) fs1=1.23
	if (zt.eq.2.) fs1=1.526
	if (zt.eq.6.) fs1=0.78
	if (zt.eq.7.) fs1=0.85
	if (zt.eq.10.) fs1=1.04
	if (zt.eq.13.) fs1=0.58
	if (zt.eq.14.) fs1=0.59
	if (zt.eq.18.) fs1=0.68
	if (zt.eq.29.) fs1=0.672
	if (zt.eq.36.) fs1=0.61
	if (zt.eq.54.) fs1=0.535
	scfk=fs1**2*1.277*z2**(2./3.)/zk**2.
	scfl1=fs1**2*1.277*z2**(2./3.)/zl1**2.
	scfl2=fs1**2*1.277*z2**(2./3.)/zl2**2.
	scfm1=fs1**2*1.277*z2**(2./3.)/zm1**2.
	scfm2=fs1**2*1.277*z2**(2./3.)/zm2**2.
	scfn=fs1**2*1.277*z2**(2./3.)/zn**2.
c ...................................................................
c      correction energie de liaison
c ....................................................................
	if (ibin.eq.1) then
	xk=2.*vu/(zk*tetak)
	gk=(1.+5.*xk+7.14*xk**2.+4.27*xk**3.+0.947*xk**4.)/(1.+xk)**5.
	ek=(1.+zt*gk/(zk*tetak))**2.
	xl1=4.*vu/(zl1*tetal1)
	gl1=1.+9.*xl1+30.2*xl1**2.+66.8*xl1**3.+100.*xl1**4.
	gl1=gl1+94.1*xl1**5.+51.3*xl1**6.+15.2*xl1**7.+1.891*xl1**8.
	gl1=gl1/(1.+xl1)**9.
	el1=(1.+zt*gl1/(zl1*tetal1))**2.
	xl2=4.*vu/(zl2*tetal2)
	gl2=1.+9.*xl2+34.7*xl2**2.+81.2*xl2**3.+112.*xl2**4.
	gl2=gl2+93.5*xl2**5.+46.6*xl2**6.+12.9*xl2**7.+1.549*xl2**8.
	gl2=gl2/(1.+xl2)**9.
	el2=(1.+zt*gl2/(zl2*tetal2))**2.
	xm1=18.*vu/(zm1*tetam1)
	gm1=(1.+5.*xm1+7.14*xm1**2.+4.27*xm1**3.+0.947*xm1**4.)
	gm1=gm1/(1.+xm1)**5.
	em1=(1.+zt*gm1/(zm1*tetam1))**2.
	xm2=18.*vu/(zm2*tetam2)
	gm2=(1.+5.*xm2+7.14*xm2**2.+4.27*xm2**3.+0.947*xm2**4.)
	gm2=gm2/(1.+xm2)**5.
	em2=(1.+zt*gm2/(zm2*tetam2))**2.
c	xn=18.*vu/(zn*tetan)
c	gn=(1.+5.*xn+7.14*xn**2.+4.27*xn**3.+0.947*xn**4.)
c	gn=gn/(1.+xn)**5.
c	en=(1.+zt*gn/(zn*tetan))**2.
      en=1.
	else
c ...................................................................
	ek=1.
	el1=1.
	el2=1.
	em1=1.
	em2=1.
	en=1.
	end if
C --------------------------------------------------------------------
C   Table des energies cinetiques de l'electron
c ---------------------------------------------
	TMAX1=2194.32*epo
	TMAX=max(tmax1,4.*bk)
c ---------------------------------------------
	IF (TMAX.LE.1000.) THEN
	IMAX=39
	ELSE IF (TMAX.LE.2000.) THEN
	IMAX=41
	ELSE IF (TMAX.LE.5000.) THEN
	IMAX=46
	ELSE IF (TMAX.LE.10000.) THEN
	IMAX=51
	ELSE IF (TMAX.LE.20000.) THEN
	IMAX=53
	ELSE IF (TMAX.LE.50000.) THEN
	IMAX=56
	ELSE IF (TMAX.LE.100000.) THEN
	IMAX=58
	ELSE IF (TMAX.LE.200000.) THEN
	IMAX=60
	ELSE IF (TMAX.LE.500000.) THEN
	IMAX=63
	ELSE IF (TMAX.LE.1000000.) THEN
	IMAX=65
	ELSE IF (TMAX.LE.2000000.) THEN
	IMAX=67
	ELSE IF (TMAX.LE.5000000.) THEN
	IMAX=70
	ELSE
	IMAX=73
	END IF
	do 101 i=1,imax
	T(i)=enel(i)
101	continue
c*********** Energies transferees en unites reduites ******
	do 106 I=1,IMAX
	WK(I)=ek*tetak+T(I)/(13.6058*zk**2.)
	WL1(I)=el1*tetal1/4.+T(I)/(13.6058*zl1**2.)
	WL2(I)=el2*tetal2/4.+T(I)/(13.6058*zl2**2.)
	WM1(I)=em1*tetam1/9.+T(I)/(13.6058*zm1**2.)
	WM2(I)=em2*tetam2/9.+T(I)/(13.6058*zm2**2.)
	WN(I)=en*tetan/4.+T(I)/(13.6058*zn**2.)
106	CONTINUE
C      ****************************************************	
C      *****   Calcul des DCS et des pertes d'energie *****
C      ****************************************************
	CALL SEKE(WK,IMAX,DCKE)
	CALL SEL1E(WL1,IMAX,DCL1E)
	CALL SEL2E(WL2,IMAX,DCL2E)
	CALL SEM1E(WM1,IMAX,DCM1E)
	CALL SEM2E(WM2,IMAX,DCM2E)
	CALL SENE(WN,IMAX,DCNE)
c --------traitement des donnees du calcul : ---------------------------
c ----------------------------------------------------debut de la boucle
	DO 107 I=1,IMAX
c -----------------------------------
c Sections efficaces differentielles:
c -----------------------------------
	sk(i)=cstk*dcke(i)/(etak*13.6058*zk**2.)
	sl1(i)=cstl1*dcl1e(i)/(etal1*13.6058*zl1**2.)
	sl2(i)=cstl2*dcl2e(i)/(etal2*13.6058*zl2**2.)
	sm1(i)=cstm1*DCM1E(i)/(etam1*13.6058*zm1**2.)
	sm2(i)=cstm2*dcm2e(i)/(etam2*13.6058*zm2**2.)
	sn(i)=cstn*dcne(i)/(etan*13.6058*zn**2.)
c --------------------------------------------------
107	continue
c----------------------------------------------------------------
c Calcul de la section efficace totale:
c --------------------------------------
	TCk=quad1(Sk,T,IMAX)
	TCl1=quad1(Sl1,T,IMAX)
	TCl2=quad1(Sl2,T,IMAX)
	TCm1=quad1(Sm1,T,IMAX)
	TCm2=quad1(Sm2,T,IMAX)
	TCn=quad1(Sn,T,IMAX)
c	write(16,*)'Bk Bl1 BM1 Bn =', Bk,Bl1,BM1,Bn
c	write(16,*) Tck,TCL1,TCm1,TCn
c ------------------------------------
c on retabli la valeur de zn
	zn=2.*zn
c--------------------------------------
	return
        END
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c						FIN de BORN.FOR
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c -----------------------------------
c fonction d'integration sur les T(I):
c -----------------------------------
	FUNCTION quad1(SOM,T,IMAX)
	DIMENSION SOM(1),T(1)
	quad1=0.
	DO 200 I=1,IMAX-1
	quad1=quad1+(SOM(I+1)+SOM(I))*(T(I+1)-T(I))/2
200	CONTINUE
	RETURN
	END
C-------------------------------------------------------------------------
	function enel(i)
	dimension el(73)
	data (el(i),i=1,73) / 0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,
     s	6.0,7.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,
     s	70.0,80.0,90.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,
     s	600.0,700.0,800.0,900.0,1000.0,1500.0,2000.0,2500.0,3000.0,
     s	3500.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0,
     s	15000.0,20000.0,30000.0,40000.0,50000.0,70000.0,100000.0,
     s	150000.0,200000.0,300000.0,400000.0,500000.0,700000.0,1000000.0,
     s	1500000.0,2000000.0,3000000.0,4000000.0,5000000.0,7000000.0,
     s	10000000.0,20000000.0 /
	enel=el(i)
	return
	end
c -----------------------------------
c fonction d'integration sur les Q(I):
c -----------------------------------
	FUNCTION quad2(SOM2,Q)
	DIMENSION SOM2(1),Q(1)
	quad2=0.
	DO 200 I=1,57
	quad2=quad2+(SOM2(I+1)+SOM2(I))*(Q(I+1)-Q(I))/2.
200	CONTINUE
	RETURN
	END
C-------------------------------------------------------------------------
	function qred(i)
	dimension red(58)
	data (red(i),i=1,58) / 0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,
     s	6.0,7.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,
     s	70.0,80.0,90.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,
     s	600.0,700.0,800.0,900.0,1000.0,1500.0,2000.0,2500.0,3000.0,
     s	3500.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0,
     s	15000.0,20000.0,30000.0,40000.0,50000.0,70000.0,100000.0 /
	qred=red(i)
	return
	end
C -----------------------------------------
	Function BTOT(z1,rk,rl1,rl2,rm1,rm2,rn)
	rks=max((rk-1.),0.)
	rl1s=max((rl1-1.),0.)
	rl2s=max((rl2-1.),0.)
	rm1s=max((rm1-1.),0.)
	rm2s=max((rm2-1.),0.)
	rm2s=max((rm2-1.),0.)
	rns=max((rn-1.),0.)
	B0=rk*(z1-0.3125*rks)**2.
	B1=0.25*rl1*(z1-0.8*rk-0.3*rl1s)**2.
	B2=0.25*rl2*(z1-rk-0.75*rl1-0.35*rl2s)**2.
	B3=1./9.*rm1*(z1-rk-0.85*(rl1+rl2)-0.35*rm1s)**2.
	B4=1./9.*rm2*(z1-(rk+rl1+rl2+rm1)-0.35*rm2s)**2.
	B5=1./16.*rn*(z1-(rk+rl1+rl2+rm1+rm2)-0.35*rns)**2.
	Btot=13.6058*(b0+b1+b2+b3+b4+b5)
	return
	end
C -----------------------------------------
c
c				Debut des routines de sections efficaces.
coooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
c	dSk(T)/dT			     pour les electrons	1S1/2
c -----------------------------------------------------------------------
	SUBROUTINE SEKE(WK,IMAX,DCKE)
	DIMENSION WK(1),DCKE(1)
	Dimension QP(58),FF(58)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parf1s/W
	write(6,*) '1s ionization'
	QP=0.
	FF=0.
	do 300 i=1,IMAX
	  W=WK(i)
	  Qm=w**2/(4.*etak)
c------ recherche du qmax de convergence -> pas d'integration------
	  F0=F1s(Qm)
	  Q=Qm
299	    Q=2.*Q
	    Fn=ANINT(100.*F1s(Q)/F0)
	    If (Fn.GE.1.) go to 299
	  PQ=Q/1000.
c------------------------------------------------------------------
	  do 301 j=1,39
	    QP(j)=qred(j)*PQ
	    Q=Qm+QP(j)
	    FF(j)=F1s(Q)
301	  continue
	  DCKE(i)=QUAD2(FF,QP)
300	continue
        RETURN
	END
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C FONCTION F1S FACTEUR DE FORME POUR ELECTRON 1S
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	FUNCTION F1S(Q)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parf1s/W
	Tol0=-1D-3
	Tol1=1D-3
	Pi=4.*ATAN(1.)
	sq=(1.+scfk/Q)**(-2.)+coak*(1.-(1+Q/scfk)**(-2.))/z2
	AK2=W-1
	AS=Q+(AK2)/3.+(1./3.)
	     IF(AK2.LT.Tol0) THEN
	ACC=SQRT(ABS(AK2))
	BS=((Q+(1.-ACC)**2.)/(Q+(1.+ACC)**2.))**(1./ACC)
	CS=1.
	    ELSE IF(AK2.LT.0) THEN
	ACC=SQRT(ABS(AK2))
	BS=EXP(-4./(Q+(1.+ACC)**2))
	CS=1.
	    ELSE IF(AK2.LT.Tol1) THEN
	BS=EXP(-4./(Q+1.-AK2))
	CS=1.
	    ELSE
	ACC=SQRT(AK2)
	BS=EXP((-2./ACC)*(ATAN2((2.*ACC),(Q-AK2+1.))))
	CS=(1.-EXP((-2.*Pi)/ACC))
             ENDIF
	DS=((Q-AK2+1.)**2.+4.*AK2)**3.
	ES=2.**7.
	F1S=sq*(ES*AS*BS)/(CS*DS*Q)
        RETURN
	END
coooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
c	dSl1(T)/dT			     pour les electrons	2S
c -----------------------------------------------------------------------
	SUBROUTINE SEL1E(WL1,IMAX,DCL1E)
	DIMENSION WL1(1),DCL1E(1)
	Dimension QP(58),FF(58)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parf2s/W1
	write(6,*) '2s ionization'
	do 400 i=1,IMAX
	  W1=WL1(i)
	  Qm=W1**2/(4.*etal1)
c------ recherche du qmax de convergence -> pas d'integration------
	  F0=F2s(Qm)
	  Q=Qm
399	    Q=2.*Q
	    Fn=ANINT(100.*F2s(Q)/F0)
	    If (Fn.GE.1.) go to 399
	  PQ=Q/1000.
c------------------------------------------------------------------
	  do 401 j=1,39
	    QP(j)=qred(j)*PQ
	    Q=Qm+QP(j)
	    FF(j)=F2s(Q)
401	  continue
	  DCL1E(i)=QUAD2(FF,QP)
400	continue
        RETURN
	END
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c FONCTION F2S FACTEUR DE FORME POUR ELECTRON 2S
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	FUNCTION F2S(Q)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parf2s/W1
	Tol0=-1D-3
	Tol1=1D-3
	Pi=4.*ATAN(1.)
	sq=(1.+scfl1/Q)**(-2.)+coal1*(1.-(1+Q/scfl1)**(-2.))/z2
	AK2=W1-0.25
	     IF(AK2.LT.Tol0) THEN
	ACC=SQRT(ABS(AK2))
	BS=((Q+(0.5-ACC)**2.)/(Q+(0.5+ACC)**2.))**(1./ACC)
	CS=1.
	    ELSE IF(AK2.LT.0) THEN
	ACC=SQRT(ABS(AK2))
	BS=EXP(-2./(Q+(0.5+ACC)**2))
	CS=1.
	    ELSE IF(AK2.LT.Tol1) THEN
	BS=EXP(-2./(Q+0.25-AK2))
	CS=1.
	    ELSE
	ACC=SQRT(AK2)
	BS=EXP((-2./ACC)*(ATAN2((ACC),(Q-AK2+0.25))))
	CS=(1.-EXP((-2.*Pi)/ACC))
	     ENDIF
	AS=(Q-AK2+0.25)**2.+AK2
	ES=2.**4.
	AL1=BS*ES/(CS*AS**5.*Q)
	A5=Q**5.
	A4=-(8./3.+11./3.*AK2)*Q**4.
	A3=(41./24.+6.*AK2+14./3.*AK2**2.)*Q**3.
	A2=(5./48.-31./24.*AK2-10./3.*AK2**2.-2.*AK2**3.)*Q**2.
	A1=(47./3840.-41./120.*AK2**2.-2./3.*AK2**3.-1./3.*AK2**4)*Q
	A0=1./768.+17./768.*AK2+7./48.*AK2**2.+11./24.*AK2**3.
     S	+2./3.*AK2**4.+1./3.*AK2**5.
	F2S=sq*AL1*(A5+A4+A3+A2+A1+A0)
        RETURN
	END
coooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
c	dSl2(T)/dT			     pour les electrons	2P
c -----------------------------------------------------------------------
	SUBROUTINE SEL2E(WL2,IMAX,DCL2E)
	DIMENSION WL2(1),DCL2E(1)
	Dimension QP(58),FF(58)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parf2p/W2
	write(6,*) '2p ionization'
	do 500 i=1,IMAX
	  W2=WL2(i)
	  Qm=w2**2/(4.*etal2)
c------ recherche du qmax de convergence -> pas d'integration------
	  F0=F2p(Qm)
	  Q=Qm
499	    Q=2.*Q
	    Fn=ANINT(100.*F2p(Q)/F0)
	    If (Fn.GE.1.) go to 499
	  PQ=Q/1000.
c------------------------------------------------------------------
	  do 501 j=1,39
	    QP(j)=qred(j)*PQ
	    Q=Qm+QP(j)
	    FF(j)=F2p(Q)
501	  continue
	  DCL2E(i)=QUAD2(FF,QP)
500	continue
        RETURN
	END
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C    FONCTION F2P FACTEUR DE FORME POUR ELECTRON 2P
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	FUNCTION F2P(Q)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parf2p/W2
	Tol0=-1D-3
	Tol1=1D-3
	Pi=4.*ATAN(1.)
	sq=(1.+scfl2/Q)**(-2.)+coal2*(1.-(1+Q/scfl2)**(-2.))/z2
	AK2=W2-0.25
	     IF(AK2.LT.Tol0) THEN
	ACC=SQRT(ABS(AK2))
	BS=((Q+(0.5-ACC)**2.)/(Q+(0.5+ACC)**2.))**(1./ACC)
	CS=1.
	    ELSE IF(AK2.LT.0) THEN
	ACC=SQRT(ABS(AK2))
	BS=EXP(-2./(Q+(0.5+ACC)**2))
	CS=1.
	    ELSE IF(AK2.LT.Tol1) THEN
	BS=EXP(-2./(Q+0.25-AK2))
	CS=1.
	    ELSE
	ACC=SQRT(AK2)
	BS=EXP((-2./ACC)*(ATAN2((ACC),(Q-AK2+0.25))))
	CS=(1.-EXP((-2.*Pi)/ACC))
	     ENDIF
	AS=(Q-AK2+0.25)**2.+AK2
	ES=2.**4.
	AL2=BS*ES/(CS*AS**5.*Q)
	A4=9./4.*Q**4.
	A3=-(0.75+3.*AK2)*Q**3.
	A2=(19./32.-0.75*AK2-0.5*AK2**2.)*Q**2.
	A1=(107./960.+41./48.*AK2+113./60.*AK2**2.+AK2**3.)*Q
	A0=1./4.*AK2**4.+5./12.*AK2**3.+7./32.*AK2**2.+3./64.*AK2
     s	+11./3072.
	F2P=sq*AL2*(A4+A3+A2+A1+A0)*(1./3.)
	RETURN
	END
coooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
c	dSm1(T)/dT			     pour les electrons	3s,p
c -----------------------------------------------------------------------
	SUBROUTINE SEM1E(WM1,IMAX,DCM1E)
	DIMENSION WM1(1),DCM1E(1)
	Dimension QP(58),FF(58)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parfm1/W1M
	write(6,*) '3s,p ionization'
	do 600 i=1,IMAX
	  W1M=WM1(i)
	  Qm=w1m**2/(4.*etam1)
c------ recherche du qmax de convergence -> pas d'integration------
	  F0=Fm1(Qm)
	  Q=Qm
599	    Q=2.*Q
	    Fn=ANINT(100.*Fm1(Q)/F0)
	    If (Fn.GE.1.) go to 599
	  PQ=Q/1000.
c------------------------------------------------------------------
	  do 601 j=1,39
	    QP(j)=qred(j)*PQ
	    Q=Qm+QP(j)
	    FF(j)=Fm1(Q)
601	  continue
	  DCM1E(i)=QUAD2(FF,QP)
600	continue
        RETURN
	END
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c FONCTION Fm1 FACTEUR DE FORME POUR ELECTRON 3s,p
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	FUNCTION Fm1(Q)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parfm1/W1M
	Tol0=-1D-3
	Tol1=1D-3
	Pi=4.*ATAN(1.)
	sq=(1.+scfm1/Q)**(-2.)+coam1*(1.-(1+Q/scfm1)**(-2.))/z2
	AK2=W1M-1/9.
	     IF(AK2.LT.Tol0) THEN
	ACC=SQRT(ABS(AK2))
	BS=((Q+(1./3.-ACC)**2.)/(Q+(1./3.+ACC)**2.))**(1./ACC)
	CS=1.
	    ELSE IF(AK2.LT.0) THEN
	ACC=SQRT(ABS(AK2))
	BS=EXP(-4./3./(Q+(1./3.+ACC)**2))
	CS=1.
	    ELSE IF(AK2.LT.Tol1) THEN
	BS=EXP(-4./3./(Q+1./9.-AK2))
	CS=1.
	    ELSE
	ACC=SQRT(AK2)
	BS=EXP((-2./ACC)*(ATAN2((2.*ACC/3.),(Q-AK2+1./9.))))
	CS=(1.-EXP((-2.*Pi)/ACC))
	     ENDIF
	AS=(Q-AK2+1./9.)**2.+4.*AK2/9.
	ES=2.**7./27.
	AL3=BS*ES/(Q*CS*AS**5.)
	A5=Q**5.
	A4=-(43./27.+11./3.*AK2)*Q**4.
	A3=(518./243.+412./81.*AK2+14./3.*AK2**2.)*Q**3.
	A2=-(442./729.+310./81.*AK2+122./27.*AK2**2.+2.*AK2**3.)*Q**2.
	A1=(1943./3.**9.+1460./3.**7.*AK2+290./243.*AK2**2.
     s	+4./27.*AK2**3.-1./3.*AK2**4.)*Q
	A0=1./3.*AK2**5.+71./81.*AK2**4.+62./81.*AK2**3.
     s	+1790./3.**8.*AK2**2.+2431./3.**10.*AK2+377./3.**11.
	Fm1=sq*AL3*(A5+A4+A3+A2+A1+A0)*(1./9.)
	RETURN
	END
coooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
c	dSm2(T)/dT			     pour les electrons	3d
c -----------------------------------------------------------------------
	SUBROUTINE SEM2E(WM2,IMAX,DCM2E)
	DIMENSION WM2(1),DCM2E(1)
	Dimension QP(58),FF(58)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parfm2/W2M
	write(6,*) '3d ionization'
	do 600 i=1,IMAX
	  W2M=WM2(i)
	  Qm=w2m**2/(4.*etam2)
c------ recherche du qmax de convergence -> pas d'integration------
	  F0=Fm2(Qm)
	  Q=Qm
	if(f0.gt.0.) then
599	    Q=2.*Q
	    Fn=ANINT(100.*Fm2(Q)/F0)
	    If (Fn.GE.1.) go to 599
	end if
	  PQ=Q/1000.
c------------------------------------------------------------------
	  do 601 j=1,39
	    QP(j)=qred(j)*PQ
	    Q=Qm+QP(j)
	    FF(j)=Fm2(Q)
601	  continue
	  DCM2E(i)=QUAD2(FF,QP)
600	continue
        RETURN
	END
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c FONCTION Fm2 FACTEUR DE FORME POUR ELECTRON 3d
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	FUNCTION Fm2(Q)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parfm2/W2M
	Tol0=-1D-3
	Tol1=1D-3
	Pi=4.*ATAN(1.)
	sq=(1.+scfm2/Q)**(-2.)+coam2*(1.-(1+Q/scfm2)**(-2.))/z2
	AK2=W2M-1/9.
	     IF(AK2.LT.Tol0) THEN
	ACC=SQRT(ABS(AK2))
	BS=((Q+(1./3.-ACC)**2.)/(Q+(1./3.+ACC)**2.))**(1./ACC)
	CS=1.
	    ELSE IF(AK2.LT.0) THEN
	ACC=SQRT(ABS(AK2))
	BS=EXP(-4./3./(Q+(1./3.+ACC)**2))
	CS=1.
	    ELSE IF(AK2.LT.Tol1) THEN
	BS=EXP(-4./3./(Q+1./9.-AK2))
	CS=1.
	    ELSE
	ACC=SQRT(AK2)
	BS=EXP((-2./ACC)*(ATAN2((2.*ACC/3.),(Q-AK2+1./9.))))
	CS=(1.-EXP((-2.*Pi)/ACC))
	     ENDIF
	AS=(Q-AK2+1./9.)**2.+4.*AK2/9.
	ES=2.**7./27.
	AL3=BS*ES/(Q*CS*AS**5.)
	A5=Q**5.
	A4=-(43./27.+11./3.*AK2)*Q**4.
	A3=(518./243.+412./81.*AK2+14./3.*AK2**2.)*Q**3.
	A2=-(442./729.+310./81.*AK2+122./27.*AK2**2.+2.*AK2**3.)*Q**2.
	A1=(1943./3.**9.+1460./3.**7.*AK2+290./243.*AK2**2.
     s	+4./27.*AK2**3.-1./3.*AK2**4.)*Q
	A0=1./3.*AK2**5.+71./81.*AK2**4.+62./81.*AK2**3.
     s	+1790./3.**8.*AK2**2.+2431./3.**10.*AK2+377./3.**11.
	Fm2=sq*AL3*(A5+A4+A3+A2+A1+A0)*(1./9.)
	RETURN
	END
coooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c Creation d'une table de SECTIONS EFFICACES DIFFERENTIELLES
c	dSN(T)/dT			     pour les electrons	N -> "2P"
c -----------------------------------------------------------------------
	SUBROUTINE SENE(WN,IMAX,DCNE)
	DIMENSION WN(1),DCNE(1)
	Dimension QP(58),FF(58)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parfn/W4
	write(6,*) 'n=4 ionization'
	do 500 i=1,IMAX
	  W4=Wn(i)
	  Qm=w4**2/(4.*etan)
c------ recherche du qmax de convergence -> pas d'integration------
	  F0=FN(Qm)
	  Q=Qm
499	    Q=2.*Q
	    F4=ANINT(100.*FN(Q)/F0)
	    If (F4.GE.1.) go to 499
	  PQ=Q/1000.
c------------------------------------------------------------------
	  do 501 j=1,39
	    QP(j)=qred(j)*PQ
	    Q=Qm+QP(j)
	    FF(j)=FN(Q)
501	  continue
	  DCNE(i)=QUAD2(FF,QP)
500	continue
        RETURN
	END
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C    FONCTION FN FACTEUR DE FORME POUR ELECTRON 2P -> n=4
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	FUNCTION FN(Q)
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn
	common/par3/coak,coal1,coal2,coam1,coam2,coan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/parfn/W4
	Tol0=-1D-3
	Tol1=1D-3
	Pi=4.*ATAN(1.)
	sq=(1.+scfn/Q)**(-2.)+coan*(1.-(1+Q/scfn)**(-2.))/z2
	AK2=W4-0.25
	     IF(AK2.LT.Tol0) THEN
	ACC=SQRT(ABS(AK2))
	BS=((Q+(0.5-ACC)**2.)/(Q+(0.5+ACC)**2.))**(1./ACC)
	CS=1.
	    ELSE IF(AK2.LT.0) THEN
	ACC=SQRT(ABS(AK2))
	BS=EXP(-2./(Q+(0.5+ACC)**2))
	CS=1.
	    ELSE IF(AK2.LT.Tol1) THEN
	BS=EXP(-2./(Q+0.25-AK2))
	CS=1.
	    ELSE
	ACC=SQRT(AK2)
	BS=EXP((-2./ACC)*(ATAN2((ACC),(Q-AK2+0.25))))
	CS=(1.-EXP((-2.*Pi)/ACC))
	     ENDIF
	AS=(Q-AK2+0.25)**2.+AK2
	ES=2.**4.
	AL2=BS*ES/(CS*AS**5.*Q)
	A4=9./4.*Q**4.
	A3=-(0.75+3.*AK2)*Q**3.
	A2=(19./32.-0.75*AK2-0.5*AK2**2.)*Q**2.
	A1=(107./960.+41./48.*AK2+113./60.*AK2**2.+AK2**3.)*Q
	A0=1./4.*AK2**4.+5./12.*AK2**3.+7./32.*AK2**2.+3./64.*AK2
     s	+11./3072.
	FN=sq*AL2*(A4+A3+A2+A1+A0)*(1./3.)
	RETURN
	END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                  FIN DES ROUTINES DES SECTIONS EFFICACES.
C				    ----------------------------------------
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c		FIN DU PROGRAMME
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
