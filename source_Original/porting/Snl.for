	subroutine snl(esp,esp3,epd3)
cc*********************************************************************
cc	    Calcule les sections efficaces d'excitation ns-np  dans
cc	l'approximation PWBA avec:
cc	    -facteurs de forme analytiques utilisant la forme recurrente
cc	    etablie a partir des calculs n=2 a n=6
cc	    -des facteurs d'ecrantage et d'antiecrantage Thomas-Fermi
cc	Les facteurs de forme sont calcules en fonction de la variable
cc	reduite k =  q/Zp (q en u.a.)
cc
cc	Les resultats sont donnes ici en unites de 10-20 cm2
cc*********************************************************************
	dimension secnl(3)
	INTEGER LIMIT,LENW
	PARAMETER(LIMIT=100,LENW=LIMIT*4)
	REAL BOUND,EPSABS,EPSREL,ABSERR,RESULT,WORK(LENW)
	INTEGER INF,NEVAL,IER,LAST,IWORK(LIMIT)
	EXTERNAL FNL
	character erfi*12
	common/don/zp,E,zt
	common /nsp/ n
	erfi='messerr.fil'
cc******integration***********
	bet=(1-(1+E/931.5)**(-2))**(0.5)
	do 30 n=1,3
	BOUND = 1.E-5/zp
	INF = 1
	EPSABS = 0.0
	EPSREL = 1.E-4
c*************************************************
	CALL INTG(FNL,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
     *	  IER,LIMIT,LENW,LAST,IWORK,WORK)
	sig1=RESULT
	if (ier.gt.0) then
	open(unit=9,file=erfi,status='unknown')
	write(9,*)' n= ',n,'  ier= ',ier
	write(9,*)' result= ',result,' abserr= ',abserr
	close(unit=9)
	end if
c*************************************************
	secnl(n)=sig1*(1.874)*zt*zt/(zp*zp*bet*bet)
 30	continue
	esp=secnl(1)
	esp3=secnl(2)
	epd3=secnl(3)
	return
	end
cc*************************************************************
	real function fnl(xx)
	real*8 q,corr,scf,stf,atf,fac,bp
	common/don/zp,E,zt
	common /nsp/ n
	scf=dble(1.13*zt**(1./3.))
	q=dble(xx*zp)
	stf=(scf)/(q)
	atf=(q)/(scf)
	corr=dble(1./(1.+stf**2.)**2.)
	corr=corr+dble((1.-(1./(1.+atf**2.)**2.))/zt)
cc*********************************************************************
	bp=dble(xx)
	if (n.eq.1) then
c 2s-2p
	fac=1.-bp**2.
	fac=18.*fac**2./(bp*(1.+bp**2.)**8.)
	elseif (n.eq.2) then
c 3s-3p
	fac=128.-1392.*bp**2.+3888.*bp**4.-2187.*bp**6.
	fac=110592.*fac**2./(bp*(4.+9.*bp**2.)**12.)
	elseif (n.eq.3) then
c 3p-3d
	fac=256.-3840.*bp**2.+30816.*bp**4.-81648.*bp**6.
	fac=fac+76545.*bp**8.
	fac=2949120.*fac/(bp*(4.+9.*bp**2.)**12.)
	end if
	fnl=fac*corr
	return
	end
cc*********************************************************************
