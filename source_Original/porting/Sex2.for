	subroutine sex2(s3s,s3p,s3d,p3s,p3p,p3d,e2s4,e2p4,e2s5,e2p5)
cc	     Version du 21/06/2000 simplifiée pour eta4
cc		 dernieres corrections : 31/08/2000 (COA et certains facteurs)
cc	version utilisant la routine d'intégration INTG
cc*********************************************************************
cc	    Calcule les sections efficaces d'excitation 2l-3l' et 2l-4l'
cc	    dans l'approximation PWBA avec facteurs d'écrantage et
cc	    d'antiécrantage dans le modèle de Thomas-Fermi
cc	Les facteurs de forme sont calculés en fonction de la variable
cc	réduite k = q/Zp (q en u.a.)
cc
cc	Les résultats sont donnés en unités de 10-20 cm2
cc*********************************************************************
	INTEGER LIMIT,LENW
	PARAMETER(LIMIT=100,LENW=LIMIT*4)
	REAL BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK(LENW)
	real*8 sig1,rn,qmin,scf,coa
	INTEGER INF,NEVAL,IER,LAST,IWORK(LIMIT)
	EXTERNAL FEX2
c**********************************************************************
	dimension sec(10)
	character erfi*12
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common /don/zp,E,zt
	common/nn/rn,coa,scf,ll
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
cc****************************
cc******integration***********
	do 20 n=1,10
	rn=dble(n)
cc**********************************
cc correction pour effet energie liaison:
cc**********************************
	if (ibin.eq.1) then
	y=4.*vu/(zl2*tetal2)
	g=1.+9.*y+34.7*y**2.+81.2*y**3.+112.*y**4.
	g=g+93.5*y**5.+46.6*y**6.+12.9*y**7.+1.549*y**8.
	g=g/(1.+y)**9.
	epsi=(1+zt*g/(zl2*tetal2))**2.
	else
	epsi=1.
	end if
cc****************************
	if(n.le.6) then
c	 2 -> 3
	qmin=dble(epsi*(zl2*tetal2/(vu*2.))*(5./36.))
	elseif(n.le.8) then
c	 2 -> 4
	qmin=dble(epsi*(zl2*tetal2/(vu*2.))*(3./16.))
	else
c	 2 -> 5
	qmin=dble(epsi*(zl2*tetal2/(vu*2.))*(21./100.))
	endif
c	calcul du coefficient d'antiscreening (08/2000)
	vs=1.75*epsi*zl2*tetal2*qmin
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
	CALL INTG(FEX2,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
     *	  IER,LIMIT,LENW,LAST,IWORK,WORK)
	sig1=RESULT
	if (ier.gt.0) then
	open(unit=9,file=erfi,status='unknown')
	write(9,*)' n= ',n,'  ier= ',ier
	write(9,*)' result= ',result,' abserr= ',abserr
	close(unit=9)
	end if
c*************************************************
	sec(n)=sig1*(1.874)*zt*zt/(zl2*zl2*bet*bet)
20	continue
	s3s=sec(1)
	s3p=sec(2)
	s3d=sec(3)
	p3s=sec(4)
	p3p=sec(5)
	p3d=sec(6)
	e2s4=sec(7)
	e2p4=sec(8)
	e2s5=sec(9)
	e2p5=sec(10)
	return
	end
cc*************************************************************
	real function FEX2(x)
	real*8 q,scf,stf,atf,corr,fact,bp,rn,coa
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common /don/zp,E,zt
	common/nn/rn,coa,scf,ll
	q=dble(x*zl2)
	stf=(scf)/(q)
	atf=(q)/(scf)
	corr=dble(1./(1.+stf**2.)**2.)
	corr=corr+coa*((1.-(1./(1.+atf**2.)**2.))/zt)
cc*********************************************************************
	bp=dble(x)
	if (rn.eq.1) then
cc 2s->3s
	fact=2.**18.*3.**7.*bp/(25.+36.*bp**2.)**10.
	fact=fact*(2875.-6984.*bp**2.+3888.*bp**4.)**2.
	elseif (rn.eq.2) then
cc 2s->3p
	fact=2.*27648.**2./(bp*(25.+36.*bp**2.)**10.)
	fact=fact*(-625.+5040.*bp**2.-3888.*bp**4.)**2.
	elseif (rn.eq.3) then
cc 2s->3d
	fact=2.*65536.**2.*3.**7.*bp/(25.+36.*bp**2.)**10.
	fact=fact*(-25.+18.*bp**2.)**2.
	elseif (rn.eq.4) then
cc 2p->3s
cc corrigé le 19/08/2009 : 3.**7 remplacé par 3.**6
	fact=2.**16.*3.**6./(bp*(25.+36.*bp**2.)**10.)
	fact=fact*(625.-14040.*bp**2.+11664.*bp**4.)**2.
	elseif (rn.eq.5) then
cc 2p->3p
	fact=2.**25.*3.**9.*bp/(25.+36.*bp**2.)**10.
	fact=fact*(6875.-23400.*bp**2.+34992.*bp**4.)
	elseif (rn.eq.6) then
cc 2p->3d
	fact=2.**24.*3.**6.*5.**2./(bp*(25.+36.*bp**2.)**10.)
	fact=fact*(3125.-5400.*bp**2.+27216.*bp**4.)
	elseif (rn.eq.7) then
cc 2s->4
	fact=405.+5120.*bp**2.+98816.*bp**4.
	fact=fact-163840.*bp**6.+65536.*bp**8.
	fact=fact*2.**20./(bp*(9.+16.*bp**2.)**9.)
	elseif (rn.eq.8) then
cc 2p->4
	fact=4194304.*(123. + 3728.*bp**2 - 3840.*bp**4 +
     s	12288.*bp**6)/(bp*(9.+16.*bp**2)**9)
	elseif (rn.eq.9) then
cc 2s->5
	fact=(20480000000.*(1555848.+29336825.*bp**2+376570000.*bp**4+ 
     -      3689500000.*bp**6-6300000000.*bp**8+2500000000.*bp**10))/
     -  (bp*(49.+100.*bp**2)**10)
	else
cc 2p->5
	fact=(128000000000.*(1643158523433.+74683227547040.*bp**2+ 
     -      993723878526400.*bp**4+4820995175360000.*bp**6 + 
     -      14083087256000000.*bp**8+45107440000000000.*bp**10+ 
     -      81423200000000000.*bp**12-79200000000000000.*bp**14+ 
     -      90000000000000000.*bp**16))/(bp*(49.+100.*bp**2)**14)
	end if
	FEX2=fact*corr
	return
	end
cc*********************************************************************
	subroutine sex3(e3s4,e3p4,e3d4,e3s5,e3p5,e3d5,e4s5,e4p5,e45,e456)
cc	     Version du 25/01/95 simplifiée pour eta4
cc		 dernieres corrections : 31/08/2000 (COA)
cc	version utilisant la routine d'intégration INTG
cc*********************************************************************
cc	  Calcule les sections efficaces d'excitation 3l-4l'
cc	    dans l'approximation PWBA avec facteurs d'ecrantage et
cc	    d'antiecrantage dans le modŠle de Thomas-Fermi
cc	Les facteurs de forme sont calcules en fonction de la variable
cc	reduite k = q/Zp (q en u.a.)
cc
cc	Les resultats sont donnes en unites de 10-20 cm2
cc*********************************************************************
	INTEGER LIMIT,LENW
	PARAMETER(LIMIT=100,LENW=LIMIT*4)
	REAL BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK(LENW)
	real*8 sig1,rn,qmin,scf,coa
	INTEGER INF,NEVAL,IER,LAST,IWORK(LIMIT)
	EXTERNAL FEX3
c**********************************************************************
	dimension sec(10)
	character erfi*12
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common /don/zp,E,zt
	common/nn/rn,coa,scf,ll
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
cc******integration***********
	do 20 n=1,10
	rn=dble(n)
cc**********************************
cc correction pour effet energie liaison:
cc**********************************
	if (ibin.eq.1) then
	y=18.*vu/(zm2*tetam2)
	g=(1.+5.*y+7.14*y**2.+4.27*y**3.+0.947*y**4.)/(1.+y)**5.
	epsi=(1+zt*g/(zm2*tetam2))**2.
	else
	epsi=1.
	end if
cc****************************
	if (n.le.3) then
c     3 vers 4	
	qmin=dble(epsi*(zm2*tetam2/(vu*2.))*(7./144.))
	vs=1.75*epsi*zm2*tetam2*qmin 
	elseif (n.le.6) then
c     3 vers 5	
	qmin=dble(epsi*(zm2*tetam2/(vu*2.))*(16./225.))
	vs=1.75*epsi*zm2*tetam2*qmin 
	elseif (n.le.9) then
c     4 vers 5	
	qmin=dble((zn*tetan/(vu*2.))*(9./400.))
	vs=1.75*epsi*zn*tetan*qmin 
	else
c     4 vers 6	
	qmin=dble((zn*tetan/(vu*2.))*(5./144.))
	vs=1.75*epsi*zn*tetan*qmin
	end if
c	calcul du coefficient d'antiscreening (08/2000)
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
	CALL INTG(FEX3,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
     *	  IER,LIMIT,LENW,LAST,IWORK,WORK)
	sig1=RESULT
	if (ier.gt.0) then
	open(unit=9,file=erfi,status='unknown')
	write(9,*)' n= ',n,'  ier= ',ier
	write(9,*)' result= ',result,' abserr= ',abserr
	close (unit=9)
	end if
c*************************************************
	if (n.ge.7) then
	sec(n)=sig1*(1.874)*zt*zt/(zn*zn*bet*bet)
	else
	sec(n)=sig1*(1.874)*zt*zt/(zm2*zm2*bet*bet)
	end if
20	continue
	e3s4=sec(1)
	e3p4=sec(2)
	e3d4=sec(3)
c	(facteur theorique en 1/n**3 = 3.049358 pour n=5)
c	les facteurs sont dans CSEC(EF) et valent 2.5 (13/09/01)
	e3s5=sec(4)
	e3p5=sec(5)
	e3d5=sec(6)
	e4s5=sec(7)
	e4p5=sec(8)
	e45=sec(9)
c     on suppose que la loi en 1/n**3 n'est pas etablie des n=6
c	(facteur theorique en 1/n**3 = 3.5412910805 pour n=6)
c	e456=sec(9)+2.5*sec(10)
	e456=sec(9)+sec(10)
	return
	end
cc*************************************************************
	real function FEX3(x)
	real*8 q,scf,stf,atf,corr,fact,bp,rn,coa
	common/par1/epo,z1,z2,etak,etal1,etal2,etan
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/don/zp,E,zt
	common/nn/rn,coa,scf,ll
	if (rn.le.6.) then
	q=dble(x*zm2)
	else
	q=dble(x*zn)
	end if
	stf=(scf)/(q)
	atf=(q)/(scf)
	corr=dble(1./(1.+stf**2.)**2.)
	corr=corr+coa*((1.-(1./(1.+atf**2.)**2.))/zt)
cc*********************************************************************
	bp=dble(x)
	if (rn.eq.1) then
cc 3s->4
	fact=4250070125.+110475272992.*bp**2.
	fact=fact+251535456000.*bp**4.
	fact=fact-8858287816704.*bp**6.
	fact=fact+30668348915712.*bp**8.
	fact=fact-29926726041600.*bp**10.
	fact=fact+8916100448256.*bp**12.
	fact=fact*18345885696./(bp*(49.+144.*bp**2.)**11.)
	elseif (rn.eq.2) then
cc 3p->4
	fact=401065441.+5869103184.*bp**2-57396708864.*bp**4
	fact=fact+313343188992.*bp**6-651422269440.*bp**8
	fact=fact+557256278016.*bp**10
	fact=fact*2.**30*3.**5/(bp*(49.+144.*bp**2)**11)
	elseif (rn.eq.3) then
cc 3d->4
	fact=11008585.-92829520.*bp**2+552524544.*bp**4
	fact=fact+417042432.*bp**6+1528823808.*bp**8
	fact=fact*2.**35*3.**7/(5.*bp*(49.+144.*bp**2)**11)
	elseif (rn.eq.4) then
cc 3s->5
	fact=87480000000.*(183744069632.+7768401510400.*bp**2+
     s	227789475840000.*bp**4+3915483300000000.*bp**6-
     s	41894536500000000.*bp**8+115844706562500000.*bp**10-
     s	103451080078125000.*bp**12+29192926025390625.*bp**14)/
     s	(bp*(64.+225*bp**2)**12)
	elseif (rn.eq.5) then
cc 3p->5
	fact=1944.d9*(10020192256.d0+4232282112.d2*bp**2+
     s	1934093376.d4*bp**4-1692872865.d5*bp**6+
     s	79597231875.d4*bp**8-152235703125.d4*bp**10+
     s	1167717041015625.d0*bp**12)/(bp*(64.+225.*bp**2)**12)
	elseif (rn.eq.6) then
cc 3d->5
	fact=69984.d11*(3014656.d0+382644224.d0*bp**2-33788808.d2*bp**4
     s   +1613793375.d1*bp**6+10308515625.d0*bp**8+512578125.d2*bp**10)/
     s   (bp*(64.+225.*bp**2)**12)
	elseif (rn.eq.7) then
cc 4s->5
	fact=(2414107905106248.d0+122375509655954625.d0*bp**2+
     s	141571857480192.d4*bp**4-42952190122096.d6*bp**6-
     s	3750709452544.d8*bp**8+898924011392.d10*bp**10-
     s	446969856.d14*bp**12+865128448.d14*bp**14-
     s	6488064.d16*bp**16+16384.d18*bp**18)
	fact=fact*2.**27*5.**7/(bp*(81.+400.*bp**2)**14)
      elseif (rn.eq.8) then
cc 4p->5
	fact=376041246141627.d0+143031288139776.d2*bp**2
	fact=fact-8084779967808.d4*bp**4-4296764470272.d6*bp**6
	fact=fact+870387186176.d8*bp**8-63160877056.d10*bp**10
	fact=fact+22278144.d14*bp**12-33816576.d14*bp**14
	fact=fact+196608.d16*bp**16
	fact=fact*2.**23*5.**10/(bp*(81 + 400*bp**2)**14)	 
      elseif (rn.eq.9) then	
cc 4->5
	fact=2217735398973.d0-546176812856.d2*bp**2
	fact=fact+7931167704.d5*bp**4-37782745600.d5*bp**6
	fact=fact+97729280000.d5*bp**8-71884800000.d5*bp**10
	fact=fact+40960000000.d5*bp**12
	fact=fact*2.**19*5.**7/(bp*(81.+400.*bp**2)**11)
	elseif (rn.eq.10) then
cc 4->6
	fact=12881310395.d0+2722585301712.d0*bp**2
	fact=fact-66947504285952.d0*bp**4+87094878234624.d1*bp**6
	fact=fact-396354332491776.d1*bp**8+960059690975232.d1*bp**10
	fact=fact-7016971052777472.d0*bp**12+3851755393646592.d0*bp**14
	fact=fact*47775744.d0/(bp*(25.+144.*bp**2)**12)
	end if
	FEX3=fact*corr
	return
	end
cc*********************************************************************
