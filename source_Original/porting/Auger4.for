	Subroutine Auger(Y,Zp,PR,YTOT)
cc..........................................
cc    Calculates final charge state distribution for each target thickness
cc    of the output files;
cc    takes into account autoionisation at the exit of target
cc..........................................
	Dimension PR(62),P1234(29,33),Y(1283),YA(1283)
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	common/secKLM/C12(209),D12(209),RAD2(209),AKL2(209),PA2(209),
     s AKM3(209),RAD3(209),ALM3(209),C3(209),D3(209),E3(209),DE3(209),
     s PA3(209),PA23(209)
	common/sec1234/P14(957),C13(957),D13(957),C4M(957),D4M(957),
     s	E4M(957),DE4M(957),PA4(957),AKLM4(957),PA123(957)
	common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4
	common/aug/AKLL,AKLM,ALMM,AM4
	YTOT=0.
	Do 1 N=327,1283
	YA(N)=Y(N)
	YTOT=YTOT+Y(N)
1	Continue
c	Calculation of correlated states (n123,n4) probabilities after autoionisation
	Do 3 N1=1,29
	Do 3 N2=1,33
	P1234(N1,N2)=0.
3	Continue
c.....................
	Do 4 N1=1,29
	Do 4 N4=33,1,-1
	I=N1-1
	M=N4-1
	INM=100*M+I
	N=NUMPP(INM)
c	N stands for indices of correlated states (between 327 and 1283)
c	effect of KLL,KLM and LMM (PA123) Auger effects on these states
c     (gives rise to very little change in most cases)
	L=N-326
	if(I.gt.2) then
		P1234(I+1,M+1)=P1234(I+1,M+1)+(1.-PA123(L))*YA(N)
		P1234(I,M+1)=P1234(I,M+1)+PA123(L)*YA(N)
	else
		P1234(I+1,M+1)=P1234(I+1,M+1)+YA(N)
	end if
4	Continue
c....................
c	effet of (KLM)NN (PA4)) Augers effect on P1234 fractions (as calculated above)
c     fractions with N-1 e- KLM and M-1 e- N initially
c     fractions with 28 e- KLM (N=29) or no N shell e- (M=1) are unchanged
	Do 5 N=1,28
c	(states with at least one hole in K or L or M shells)
	Do 5 M=33,2,-1
	  INM=100*(M-1)+N-1
	  NM=NUMPP(INM)-326
	  If (M.gt.2) then
c     radiative decay 
		P1234(N+1,M-1)=P1234(N+1,M-1)+(1.-PA4(NM))*P1234(N,M)
c     Auger decay
		P1234(N+1,M-2)=P1234(N+1,M-2)+PA4(NM)*P1234(N,M)
c     NO Auger decay (as a test)
c		P1234(N+1,M-1)=P1234(N+1,M-1)+PA4(NM)*P1234(N,M)
c
	  P1234(N,M)=0.
c
	  else if (M.eq.2) then
c	only one N shell electron :
		P1234(N+1,M-1)=P1234(N+1,M-1)+P1234(N,M)
	  P1234(N,M)=0.
	  end if
5	Continue
c..................
c	conbine configurations and calculates final probabilities
	Do 6 N=1,62
	 PR(N)=0.
6	Continue
	Do 7 N=1,29
	Do 7 M=1,33
	Ne1=N+M-1
	If (Ne1.le.(Zp+1)) then
	PR(Ne1)=PR(Ne1)+P1234(N,M)
	Else
	PR(Zp+1)=PR(Zp+1)+P1234(N,M)
	End If
7	Continue
	Do 100 N=1,61
	PP=Pr(N)
	If(PP.lt.0.0000001) then
	PP=0.
	PR(N)=PP
	End If
100	Continue
	Return
	End
	
