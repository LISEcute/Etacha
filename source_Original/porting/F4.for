	subroutine F(NEQ,X,U,UP)
cc	*********************************************************
cc	création des équations différentielles (fonction F)
cc	*********************************************************
	dimension U(1284),UP(1284)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
cc--------------------------------------------------------------
cc	ICO1=100*I+10*J+K
cc	ICO2=100*n3+n12
cc	ICO3=100*n4+n123
cc
cc	I= nb elec 2p , J= nb elec 2s , K= nb elec 1s
cc	n12=nb elec 1+2 , n3=nb elec 3, n123=nb elec 1+2+3 , n4=nb elec 4
cc
cc	y1s= nb moyen e- 1s , y2s= ...etc
cc
cc	This version takes into account in "partially" correlated way
cc	M shell lectrons, and similarly N shell electrons :
cc
cc	- for n=1 et 2 , 63 correlated states of the type {I,J,K}
cc	(fractions of ions with I e- 2p AND J e- 2s AND K e- 1s)
cc
cc	- for n=3 , 3+7+11 = 21 independent states
cc	(fractions of ions with 0,1,2 (,..) e- in 3s,3p and 3d)
cc	- for n=1+2+3 , 11*19=209 correlated states of the type {n12,n3}
cc
cc	- for n=4 , 33 independent states
cc	- for n=1+2+3+4 , 29*33=957 correlated states of the type {n123,n4}
cc
cc	corresponding to 63+21+209+33+957=1283 states and equations...
cc
cc--------------------------------------------------------------
	common/seceff/sec(34),cor(48),secs(48)
	common/secKLM/C12(209),D12(209),RAD2(209),AKL2(209),PA2(209),
     s AKM3(209),RAD3(209),ALM3(209),C3(209),D3(209),E3(209),DE3(209),
     s PA3(209),PA23(209)
	common/sec1234/P14(957),C13(957),D13(957),C4M(957),D4M(957),
     s	E4M(957),DE4M(957),PA4(957),AKLM4(957),PA123(957)
	common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4
	common/aug/AKLL,AKLM,ALMM,AM4
	common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn
	common/proj/zp
cc--------------------------------------------------------------
cc	Cross sections ...
cc--------------------------------------------------------------
	nequ=nequ+1
	dum=x
	Ndum=NEQ
	nzp=int(zp)
	C1s=sec(1)
	C2s=sec(2)
	C2p=sec(3)
	C3s=sec(4)
	C3p=sec(5)
	C3d=sec(6)
	C4=sec(7)
	D1s=sec(8)
	D2s=sec(9)
	D2p=sec(10)
	D3s=sec(11)
	D3p=sec(12)
	D3d=sec(13)
	D4=sec(14)
	e2s=sec(15)
	e2p=sec(16)
	e3s=sec(17)
	e3p=sec(18)
	e3d=sec(19)
	e4=sec(20)
	s3s=sec(21)
	s3p=sec(22)
	s3d=sec(23)
	s4=sec(24)
	p3s=sec(25)
	p3p=sec(26)
	p3d=sec(27)
	p4=sec(28)
	e3s4=sec(29)
	e3p4=sec(30)
	e3d4=sec(31)
	esp=sec(32)
	esp3=sec(33)
	epd3=sec(34)
cc--------------------------------------------------------------
cc	actual cross sections for 1s,2s,2p states
cc	excitation included in rad3p1, rad3p2, rad 3s & rad 3d
cc	exchanges with n=4 included (wheighted with yn or 32-yn)
cc--------------------------------------------------------------
	C1se=C1S+y3s*e3s+y3p*Rad3p1+y3d*e3d
	C1se=C1se+yn*e4
	C2se=C2S+y3s*s3s+y3p*Rad3p2+y3d*s3d+ALMM*ymm
	C2se=C2se+yn*s4
	C2pe=C2P+y3s*Rad3s+y3p*p3p+y3d*Rad3d+ALMM*ymm
	C2pe=C2pe+yn*p4
	D1se=D1s+(2.-y3s)*e3s+(6.-y3p)*e3p+(10.-y3d)*e3d
	D1se=D1se+(32.-yn)*e4
	D2se=D2s+(2.-y3s)*s3s+(6.-y3p)*s3p+(10.-y3d)*s3d
	D2se=D2se+(32.-yn)*s4
	D2pe=D2p+(2.-y3s)*p3s+(6.-y3p)*p3p+(10.-y3d)*p3d
	D2pe=D2pe+(32.-yn)*p4
cc--------------------------------------------------------------
cc--------------------------------------------------------------
cc	Loop for calculation of differential equations
cc	for the 63 correlated states {1s,2s,2p}
cc--------------------------------------------------------------
	do 200 N=1,63
	I=II(N)
	J=JJ(I,N)
	K=KK(I,J,N)
	L=100*I+10*J+K
	nelec=I+J+K
cc------------------------------------------------------------------
cc electron loss processes (exchanges with n=3 & 4 included)
cc------------------------------------------------------------------
	UP(N)=-U(N)*(I*D2Pe+J*D2Se+K*D1Se)
	UP(N)=UP(N)-U(N)*((2-K)*I*RAD+2*(K+J-K*J)*E2S)
	UP(N)=UP(N)-U(N)*(2*(3*K+I-K*I)*E2P)
	UP(N)=UP(N)-U(N)*(2*(3*J+I-I*J)*ESP)
	if (nelec.lt.nzp)then
	UP(N)=UP(N)-U(N)*((6-I)*C2Pe+(2-J)*C2Se+(2-K)*C1Se)
	end if
cc.......... KLL & KLM Auger ....................................
cc	(electron loss)
cc.................................................................
	if (K.lt.2) then
	  if (I.ge.2) then
	    UP(N)=UP(N)-U(N)*AKLL*I*(I-1)*(2-K)
	  end if
	  if (I.ge.1) then
	    UP(N)=UP(N)-U(N)*AKLM*I*ym*(2-K)
	      if (J.ge.1) then
		UP(N)=UP(N)-U(N)*AKLL*2.*I*J*(2-K)
	      end if
	  end if
	  if (J.ge.1) then
	    UP(N)=UP(N)-U(N)*AKLM*J*ym*(2-K)
	  end if
	  if (J.ge.2) then
	    UP(N)=UP(N)-U(N)*AKLL*J*(J-1)*(2-K)
	  end if
	end if
cc.................................................................
cc	electron gain processses
cc.................................................................
	if (K.ne.0) then
cc........ Augers KLL ........................
	  if (I.le.4) then
	    M=NUM(L+199)
	    UP(N)=UP(N)+(3-K)*(I+2)*(I+1)*AKLL*U(M)
	  end if
	  if ((I.le.5).and.(J.le.1)) then
	    M=NUM(L+109)
	    UP(N)=UP(N)+(3-K)*2*(I+1)*(J+1)*AKLL*U(M)
	  end if
	  if (J.eq.0) then
	    M=NUM(L+19)
	    UP(N)=UP(N)+(3-K)*2.*AKLL*U(M)
	  end if
cc........ Augers KLM ........................
	  if (I.le.5) then
	    M=NUM(L+99)
	    UP(N)=UP(N)+(3-K)*(I+1)*ym*AKLM*U(M)
	  end if
	  if (J.le.1) then
	    M=NUM(L+9)
	    UP(N)=UP(N)+(3-K)*(J+1)*ym*AKLM*U(M)
	  end if
cc............................................
	end if
cc--------------------------------------------------------------
cc States are populated by electron gain if the considered states possess
cc electrons, 
cc and intershell processes if states possess holes
cc--------------------------------------------------------------
cc    test on unpopulated states if nelec > zp
      if(nelec.le.nzp) then
cc..............................................................
cc		    2p state
cc..............................................................
	if (I.ne.0) then
	M=NUM(L-100)
	UP(N)=UP(N)+(7-I)*C2Pe*U(M)
	if (K.ne.2) then
	M=NUM(L-99)
	UP(N)=UP(N)+(7-I)*(K+1)*E2P*U(M)
	end if
	if (J.ne.2) then
	M=NUM(L-90)
	UP(N)=UP(N)+(7-I)*(J+1)*ESP*U(M)
	end if
	end if
cc..............................................................
cc		    2s state
cc..............................................................
	if (J.ne.0) then
	M=NUM(L-10)
	UP(N)=UP(N)+(3-J)*C2Se*U(M)
	if (K.ne.2) then
	M=NUM(L-9)
	UP(N)=UP(N)+(3-J)*(K+1)*E2S*U(M)
	end if
	if (I.ne.6) then
	M=NUM(L+90)
	UP(N)=UP(N)+(3-J)*(I+1)*ESP*U(M)
	end if
	end if
cc..............................................................
cc		    1s state
cc..............................................................
	if (K.ne.0) then
	M=NUM(L-1)
	UP(N)=UP(N)+(3-K)*C1Se*U(M)
	if (I.ne.6) then
	M=NUM(L+99)
	UP(N)=UP(N)+(3-K)*(I+1)*(RAD+E2P)*U(M)
	end if
	if (J.ne.2) then
	M=NUM(L+9)
	UP(N)=UP(N)+(3-K)*(J+1)*E2S*U(M)
	end if
	end if
cc--------------------------------------------------------------
cc    end of test on unpopulated states if nelec > zp
      end if
cc--------------------------------------------------------------
cc Sates are populated by ionization if if the considered states possess
cc holes 
cc--------------------------------------------------------------
cc
cc		    2p state
	if (I.ne.6) then
	M=NUM(L+100)
	UP(N)=UP(N)+(1+I)*D2Pe*U(M)
	end if
cc		    2s state
	if (J.ne.2) then
	M=NUM(L+10)
	UP(N)=UP(N)+(1+J)*D2Se*U(M)
	end if
cc		    1s state
	if (K.ne.2) then
	M=NUM(L+1)
	UP(N)=UP(N)+(K+1)*D1Se*U(M)
	end if
 200	continue
cc--------------------------------------------------------------
cc    actual cross sections for 3s,3p,3d states
cc	exchanges with n=4 included (wheighted with yn or 32-yn)
cc--------------------------------------------------------------
	SLS=d3s+(6.-y2p)*rad3s+(2.-y2s)*s3s+(2.-y1s)*e3s+(6.-y3p)*esp3
	SLS=SLS+(2.-y1s)*yl*AKLM
	SLS=SLS+(32.-yn)*e3s4
	SLP=d3p+(6.-y2p)*p3p+(2.-y2s)*rad3p2+(2.-y1s)*rad3p1
	SLP=SLP+(2.-y3s)*esp3+(10.-y3d)*epd3
	SLP=SLP+(2.-y1s)*yl*AKLM
	SLP=SLP+(32.-yn)*e3p4
	SLD=d3d+(6.-y2p)*rad3d+(2.-y2s)*s3d+(2.-y1s)*e3d+(6.-y3p)*epd3
	SLD=SLD+(2.-y1s)*yl*AKLM
	SLD=SLD+(32.-yn)*e3d4
	SGS=c3s+y2p*p3s+y2s*s3s+y1s*e3s+y3p*esp3
	SGS=SGS+yn*e3s4+AM4*ynn+Rad4*yn
	SGP=c3p+y2p*p3p+y2s*s3p+y1s*e3p+y3s*esp3+y3d*epd3
	SGP=SGP+yn*e3p4+AM4*ynn+Rad4*yn
	SGD=c3d+y2p*p3d+y2s*s3d+y1s*e3d+y3p*epd3
	SGD=SGD+yn*e3d4+AM4*ynn+Rad4*yn
cc--------------------------------------------------------------
cc	differential equations for 3s,3p et 3d populations :
cc	(no loop)
cc			64 to 66 -> 3s states with 0,1 or 2 e-
cc			67 to 73 -> 3p states with 0,1 ... 6 e-
cc			74 to 84 -> 3d states with 0,1 ... 10 e-
cc--------------------------------------------------------------
	UP(64)=-U(64)*2.*SGS+U(65)*SLS
	UP(65)=-U(65)*(SGS+SLS)+U(64)*2.*SGS +U(66)*2.*SLS
	UP(66)=-U(66)*2.*SLS+U(65)*SGS
cc
	UP(67)=-U(67)*(6.*SGP)+U(68)*1.*SLP
	UP(68)=-U(68)*(5.*SGP+1.*SLP)+U(67)*6.*SGP +U(69)*2.*SLP
	UP(69)=-U(69)*(4.*SGP+2.*SLP)+U(68)*5.*SGP +U(70)*3.*SLP
	UP(70)=-U(70)*(3.*SGP+3.*SLP)+U(69)*4.*SGP +U(71)*4.*SLP
	UP(71)=-U(71)*(2.*SGP+4.*SLP)+U(70)*3.*SGP +U(72)*5.*SLP
	UP(72)=-U(72)*(1.*SGP+5.*SLP)+U(71)*2.*SGP +U(73)*6.*SLP
	UP(73)=-U(73)*(6.*SLP)+U(72)*1.*SGP
cc
	UP(74)=-U(74)*(10.*SGD)+U(75)*1.*SLD
	UP(75)=-U(75)*(9.*SGD+1.*SLD)+U(74)*10.*SGD +U(76)*2.*SLD
	UP(76)=-U(76)*(8.*SGD+2.*SLD)+U(75)*9.*SGD +U(77)*3.*SLD
	UP(77)=-U(77)*(7.*SGD+3.*SLD)+U(76)*8.*SGD +U(78)*4.*SLD
	UP(78)=-U(78)*(6.*SGD+4.*SLD)+U(77)*7.*SGD +U(79)*5.*SLD
	UP(79)=-U(79)*(5.*SGD+5.*SLD)+U(78)*6.*SGD +U(80)*6.*SLD
	UP(80)=-U(80)*(4.*SGD+6.*SLD)+U(79)*5.*SGD +U(81)*7.*SLD
	UP(81)=-U(81)*(3.*SGD+7.*SLD)+U(80)*4.*SGD +U(82)*8.*SLD
	UP(82)=-U(82)*(2.*SGD+8.*SLD)+U(81)*3.*SGD +U(83)*9.*SLD
	UP(83)=-U(83)*(1.*SGD+9.*SLD)+U(82)*2.*SGD +U(84)*10.*SLD
	UP(84)=-U(84)*(10.*SLD)+U(83)*1.*SGD
cc
cc........... AUGERS LMM .......................................
	AL=ALMM*(8.-yl)
	ysp=2.*(y3s+y3p)
	ysd=2.*(y3s+y3d)
	ypd=2.*(y3p+y3d)
cc.................................................
c
	UP(64)=UP(64)+AL*(+U(65)*1.*ypd+U(66)*2.)
	UP(65)=UP(65)+AL*(-U(65)*1.*ypd+U(66)*2.*ypd)
	UP(66)=UP(66)+AL*(-U(66)*2.*(ypd+1.))
c
	UP(67)=UP(67)+AL*(+U(68)*1.*ysd+U(69)*2.)
	UP(68)=UP(68)+AL*(-U(68)*1.*(ysd+0.)+U(69)*2.*ysd+U(70)*6. )
	UP(69)=UP(69)+AL*(-U(69)*2.*(ysd+1.)+U(70)*3.*ysd+U(71)*12.)
	UP(70)=UP(70)+AL*(-U(70)*3.*(ysd+2.)+U(71)*4.*ysd+U(72)*20.)
	UP(71)=UP(71)+AL*(-U(71)*4.*(ysd+3.)+U(72)*5.*ysd+U(73)*30.)
	UP(72)=UP(72)+AL*(-U(72)*5.*(ysd+4.)+U(73)*6.*ysd)
	UP(73)=UP(73)+AL*(-U(73)*6.*(ysd+5.))
c
	UP(74)=UP(74)+AL*(+U(75)*1.*ysp+U(76)*2.)
	UP(75)=UP(75)+AL*(-U(75)*1.*(ysp+0.)+U(76)*2.*ysp+U(77)*6. )
	UP(76)=UP(76)+AL*(-U(76)*2.*(ysp+1.)+U(77)*3.*ysp+U(78)*12.)
	UP(77)=UP(77)+AL*(-U(77)*3.*(ysp+2.)+U(78)*4.*ysp+U(79)*20.)
	UP(78)=UP(78)+AL*(-U(78)*4.*(ysp+3.)+U(79)*5.*ysp+U(80)*30.)
	UP(79)=UP(79)+AL*(-U(79)*5.*(ysp+4.)+U(80)*6.*ysp+U(81)*42.)
	UP(80)=UP(80)+AL*(-U(80)*6.*(ysp+5.)+U(81)*7.*ysp+U(82)*56.)
	UP(81)=UP(81)+AL*(-U(81)*7.*(ysp+6.)+U(82)*8.*ysp+U(83)*72.)
	UP(82)=UP(82)+AL*(-U(82)*8.*(ysp+7.)+U(83)*9.*ysp+U(84)*90.)
	UP(83)=UP(83)+AL*(-U(83)*9.*(ysp+8.)+U(84)*10.*ysp)
	UP(84)=UP(84)+AL*(-U(84)*10.*(ysp+9.))
cc    **********************************************************
cc                  correlated states	
cc--------------------------------------------------------------
cc	Loop for calculation of differential equations
cc	for the 209 correlated states {(1+2)*N,3l*M}
cc    Actual cross sections as a function of state code are
cc    calculated in secmean.for
cc--------------------------------------------------------------
	Do 300 N=85,293
	M=IM(N)
	K=IKL(M,N)
	L=100*M+K
	nelec=M+K
	NN=N-84
cc------------------------------------------------------------------
cc    electron loss processes (exchanges with n=3 & 4 included)
cc	E3=2->3,DE3=3->2
cc	!! RAD3 is included in DE3, and RAD2 do not change populations !!
cc------------------------------------------------------------------
	UP(N)=-U(N)*(E3(NN)+DE3(NN)+AKM3(NN)+ALM3(NN))
	UP(N)=UP(N)-U(N)*(D12(NN)+AKL2(NN)+D3(NN))
	if(nelec.lt.nzp) then
	UP(N)=UP(N)-U(N)*(C12(NN)+C3(NN))
	end if
cc.................................................................
cc	electron gain processes
cc.................................................................
	if(nelec.le.nzp) then
c............................	
	if (K.ne.0) then
	  IJ=NUMP(L-1)
	  KL=IJ-84
	  UP(N)=UP(N)+C12(KL)*U(IJ)
	  if (M.ne.18) then
	    IJ=NUMP(L+99)
	    KL=IJ-84
	    UP(N)=UP(N)+DE3(KL)*U(IJ)
	  end if
	end if
c	
	if (K.ne.10) then
	  IJ=NUMP(L+1)
	  KL=IJ-84
	  UP(N)=UP(N)+D12(KL)*U(IJ)
	  if (M.ne.0) then
	    IJ=NUMP(L-99)
	    KL=IJ-84
	    UP(N)=UP(N)+E3(KL)*U(IJ)
	  end if
	end if
c	
	if (M.ne.0) then
	  IJ=NUMP(L-100)
	  KL=IJ-84
	  UP(N)=UP(N)+C3(KL)*U(IJ)
	end if
	if (M.ne.18) then
	  IJ=NUMP(L+100)
	  KL=IJ-84
	  UP(N)=UP(N)+D3(KL)*U(IJ)
	end if	
cc.......... Augers   ....................................
	if (K.le.8) then
	    IJ=NUMP(L+1)
		KL=IJ-84
	    UP(N)=UP(N)+AKL2(KL)*U(IJ)
	end if
	If((K.ne.0).and.(M.le.17)) then
	    IJ=NUMP(L+100)
		KL=IJ-84
		UP(N)=UP(N)+AKM3(KL)*U(IJ)
	End If
	If((K.ne.0).and.(M.le.16)) then
	    IJ=NUMP(L+199)
		KL=IJ-84
		UP(N)=UP(N)+ALM3(KL)*U(IJ)
	End If
c.....................
      end if
300	continue
cc--------------------------------------------------------------
cc	actual cross sections for n=4 states
cc--------------------------------------------------------------
	SG4=C4+y1s*e4+y2s*s4+y2p*p4+y3s*e3s4+y3p*e3p4+y3d*e3d4
	SL4=D4+(2.-y1s)*e4+(2.-y2s)*s4+(6.-y2p)*p4
	SL4=SL4+(2.-y3s)*e3s4+(6.-y3p)*e3p4+(10.-y3d)*e3d4
	SL4=SL4+((2.-y3s)+(6.-y3p)+(10.-y3d))*Rad4
	AM=AM4*(18.-ym)
cc--------------------------------------------------------------
cc	differential equations for n=4 population :
cc			294 to 326 -> 33 n=4 states with 0,1 .....or 32 e-
cc    (corrected 12/12)
cc--------------------------------------------------------------
	Do 400 i=1,33
	ri=real(i)
	UP1=0.
	UP2=0.
	UP3=0.
	UP4=0.
	UP5=0.
	UP6=0.
	if (i.ne.1)  UP1=(34.-ri)*SG4*U(i+292)
	if (i.ne.33) UP2=(ri)*SL4*U(i+294)
	if (i.ne.33)  UP3=-(33.-ri)*SG4*U(i+293)
	if (i.ne.1) UP4=-(ri-1.)*SL4*U(i+293)
	if (i.lt.32) UP5=AM*ri*(ri+1.)*U(i+295)
	if (i.gt.2)  UP6=-AM*(ri-1.)*(ri-2.)*U(i+293)
	UP(i+293)=UP1+UP2+UP3+UP4+UP5+UP6
400	continue
cc--------------------------------------------------------------
cc	Loop for calculation of differential equations
cc	for the 957 correlated {(1+2+3)*N1,4l*N4} states
cc--------------------------------------------------------------
	Do 500 NN=1,957
	N=NN+326
	NEN=IN(N)
	K=IKM(NEN,N)
	L=100*NEN+K
	nelec=NEN+K
cc------------------------------------------------------------------
cc electron loss processes (exchanges with 4 included)
cc	
cc	
cc------------------------------------------------------------------
	UP(N)=-U(N)*(D13(NN)+D4M(NN)+E4M(NN)+DE4M(NN)+AKLM4(NN))
	if(nelec.lt.nzp) then
	UP(N)=UP(N)-U(N)*(C13(NN)+C4M(NN))
	end if
cc.................................................................
cc	electron gain processes
cc.................................................................
	if(nelec.le.nzp) then
c............................	
	if (K.ne.0) then
c	123 capture and 4->123 desexcitation
	IJ=NUMPP(L-1)
	KL=IJ-326
	UP(N)=UP(N)+C13(KL)*U(IJ)
	if (NEN.ne.32) then
	IJ=NUMPP(L+99)
	KL=IJ-326
	UP(N)=UP(N)+DE4M(KL)*U(IJ)
	end if
	end if
c	
	if (K.ne.28) then
c	123 ionization and 123->4 excitation
	IJ=NUMPP(L+1)
	KL=IJ-326
	UP(N)=UP(N)+D13(KL)*U(IJ)
	if (NEN.ne.0) then
	IJ=NUMPP(L-99)
	KL=IJ-326
	UP(N)=UP(N)+E4M(KL)*U(IJ)
	end if
	end if
c	
	if (NEN.ne.0) then
c	n=4 capture
	IJ=NUMPP(L-100)
	KL=IJ-326
	UP(N)=UP(N)+C4M(KL)*U(IJ)
	end if
	if (NEN.ne.32) then
c	n=4 ionization
	IJ=NUMPP(L+100)
	KL=IJ-326
	UP(N)=UP(N)+D4M(KL)*U(IJ)
	end if
c.....................
      end if	
cc.......... (KLM)NN Augers ....................................
	If((K.ne.0).and.(NEN.le.30)) then
	    IJ=NUMPP(L+199)
		KL=IJ-326
		UP(N)=UP(N)+AKLM4(KL)*U(IJ)
	End If
500	continue
	end
cc******************************************************************
	subroutine ini(ICO1,ICO2,ICO3,NUM1,NUM2,NUM3)
cc...............................................................
c	Construction of code table ICO(n)
cc...............................................................
	dimension ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	DO 21 m=1,733
	NUM1(m)=0
 21	continue
	DO 22 m=1,1911
	NUM2(m)=0
 22	continue
	n1=0
	DO 31 i=0,6
	DO 32 j=0,2
	DO 33 k=0,2
	n1=n1+1
	ICO1(n1)=i*100+j*10+k
	ICOD1=ICO1(n1)+1
	NUM1(ICOD1)=n1
 33	continue
 32	continue
 31	continue
	n2=84
	DO 41 ii=0,18
	DO 42 jj=0,10
	n2=n2+1
	ICO2(n2)=ii*100+jj
	ICOD2=ICO2(n2)+1
	NUM2(ICOD2)=n2
 42	continue
 41	continue
	n3=326
	DO 51 ii=0,32
	DO 52 jj=0,28
	n3=n3+1
	ICO3(n3)=ii*100+jj
	ICOD3=ICO3(n3)+1
	NUM3(ICOD3)=n3
 52	continue
 51	continue
	return
	end
cc****************************************************************
	function num(NCO)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	num=NUM1(NCO+1)
	return
	end
CC      *********************************************************
	function II(num)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	II=ICO1(num)/100
	return
	end
CC      *********************************************************
	function JJ(I,num)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	JJ=(ICO1(num)-I*100)/10
	return
	end
CC      *********************************************************
	function KK(I,J,num)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	KK=(ICO1(num)-I*100-J*10)
	return
	end
cc****************************************************************
	function numP(NCO)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	numP=NUM2(NCO+1)
	return
	end
CC      *********************************************************
	function IM(num)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	IM=ICO2(num)/100
	return
	end
CC      *********************************************************
	function IKL(M,num)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	IKL=ICO2(num)-M*100
	return
	end
cc****************************************************************
	function numPP(NCO)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	numPP=NUM3(NCO+1)
	return
	end
CC      *********************************************************
	function IN(num)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	IN=ICO3(num)/100
	return
	end
CC      *********************************************************
	function IKM(N,num)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	IKM=ICO3(num)-N*100
	return
	end
