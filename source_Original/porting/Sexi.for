	subroutine sexi(e2s,e2p,e3s,e3p,e3d,es4,es5)
cc	     Version du 02/03/94 simplifiee pour eta3
cc		 dernieres corrections : 31/08/2000 (COA)
cc	version avec makefile utilisant la routine d'integration INTG
cc*********************************************************************
cc	    Calcule les sections efficaces d'excitation 1s-2s, 1s-2p et
cc	1s-n dans l'approximation PWBA avec:
cc	    -facteurs de forme analytiques d'Anholt
cc	    -facteurs d'ecrantage et d'antiecrantage dans le modŠle
cc	     de Thomas-Fermi
cc	Les facteurs de forme sont calcules en fonction de la variable
cc	reduite k =  q/Zp (q en u.a.)
cc
cc	Les resultats sont donnes en unites de 10-20 cmý
cc*********************************************************************
	INTEGER LIMIT,LENW
	PARAMETER(LIMIT=100,LENW=LIMIT*4)
	REAL BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK(LENW)
	real*8 sig1,rn,qmin,scf,coa
	INTEGER INF,NEVAL,IER,LAST,IWORK(LIMIT)
	EXTERNAL FEX
c**********************************************************************
	character erfi*12
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common /don/zp,E,zt
	common /nn/ rn,coa,scf,ll
	common/corr/ibin
	erfi='messerr.fil'
	bet=(1.-(1.+E/931.5)**(-2.))**(0.5)
	vu=137.036*bet
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
	scf=dble(fs1*1.13*zt**(1./3.))
	nmin=2
	nmax=5
	som1=0.
	so51=0.
cc****************************
cc******integration***********
	do 20 n=nmin,nmax
	if (n.eq.2) then
	    ll=1
	end if
	if (n.eq.3) then
	    ll=4
	end if
	if (n.gt.3) then
	    ll=7
	end if
21	rn=dble(n)
cc**********************************
cc correction pour effet energie liaison:
cc**********************************
	if(ibin.eq.1) then
	y=2.*vu/(zk*tetak)
	g=(1.+5.*y+7.14*y**2.+4.27*y**3.+0.947*y**4.)/(1.+y)**5.
	epsi=(1+zt*g/(zk*tetak))**2.
	else
	epsi=1.
	end if
cc****************************
c     pmin=deltaE/2v => qmin=deltaE/2v/Zp
c
	qmin=dble(epsi*(zk*tetak/(vu*2.))*(1.-1./(rn*rn)))
c	calcul du coefficient d'antiscreening (08/2000)
	vs=1.75*epsi*zk*tetak*qmin
	if (vu.gt.vs) then
	coas=1.-(vs/vu)**2
	else
	coas=0.
	endif
	coa=dble(coas)
	BOUND = qmin
	INF = 1
	EPSABS = 0.0
	EPSREL = 1.E-4
c*************************************************
	CALL INTG(FEX,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
     *	  IER,LIMIT,LENW,LAST,IWORK,WORK)
	sig1=RESULT
	if (ier.gt.0) then
	open(unit=9,file=erfi,status='unknown')
	write(9,*)' n= ',n,'  ier= ',ier
	write(9,*)' result= ',result,' abserr= ',abserr
	close(unit=9)
	end if
c*************************************************
	sec1=sig1*(1.874)*zt*zt/(zk*zk*bet*bet)
	if(ll.eq.3.or.ll.eq.7) then
	som1=som1+sec1
	end if
	if (ll.lt.3) then
	    if(ll.eq.1) e2s=sec1
	    if(ll.eq.2) e2p=sec1
	ll=ll+1
	go to 21
	end if
	if (ll.gt.3.and.ll.lt.7) then
	    if(ll.eq.4) e3s=sec1
	    if(ll.eq.5) e3p=sec1
	    if(ll.eq.6) e3d=sec1
	ll=ll+1
	go to 21
	end if
	if (n.eq.4) then
	es4=sec1
	end if
	if (n.gt.4) then
	so51=so51+sec1
	end if
20	continue
c	so51=so51+4.52*sec1
	es5=so51
	return
	end
cc*************************************************************
	real function FEX(x)
	real*8 q,scf,stf,coa,corr,fact,bp,rn
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common /don/zp,E,zt
	common /nn/ rn,coa,scf,ll
	q=dble(x*zk)
	stf=(scf)/(q)
	atf=(q)/(scf)
	corr=dble(1./(1.+stf**2.)**2.)
	corr=corr+coa*dble((1.-(1./(1.+atf**2.)**2.))/zt)
cc*********************************************************************
	bp=dble(x)
	if (ll.eq.1) then
	   fact=64.*bp/(bp*bp+9./4.)**6.
	elseif (ll.eq.2) then
	   fact=144./(bp*(bp*bp+9./4.)**6.)
	elseif (ll.eq.4) then
	   fact=1119744.*bp*(16.+27.*bp*bp)**2./(9.*bp*bp+16.)**8.
	elseif (ll.eq.5) then
	   fact=2985984.*(16.+27.*bp*bp)**2./(bp*(9.*bp*bp+16.)**8.)
	elseif (ll.eq.6) then
	   fact=573308928.*bp/(9.*bp*bp+16.)**8.
	else
	   fact=(512./(3.*rn*rn*rn))*(3.*bp*bp+1.-1./(rn*rn))
	   fact=fact*(bp*bp+(1.-1./rn)**2.)**(rn-3.)
	   fact=fact/(bp*(bp*bp+(1.+1./rn)**2.)**(rn+3.))
	end if
	FEX=fact*corr
	return
	end
cc*********************************************************************
