	subroutine SEIK(ck,cl,cm,cn,ct,cs)
cc*********************************************************************	
cc			      Version du 20/01/98, revisitée en 2012, avec :
cc            - calcul du nombre d'electrons modifie
cc            - possibilite d'avoir Qp et Zp differents
cc*********************************************************************
cc	    Calcule les sections efficaces de capture non radiative
cc	     dans l'approximation eikonale
cc	(formules 8 p 3295  Meyerhof et al. Phys. Rev A 32 (1985) 3291 )
cc*********************************************************************
	dimension oc(3)
	dimension SIG(30)
	dimension ceik(30,3)
	common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn
	common/don/zp,E,zt
	do 10 i=1,30
	do 10 j=1,3
	ceik(i,j)=0.
10	continue
c     yo is the mean number of electrons in n=5
c     not yet used, but may be used one day...
      yo=0.
	QM=Zp-(y1s+yl+ym+yn+yo)
	qp=QM
c     .........	
	pi=4.*atan(1.)
	alph=1./137.036
	s0=2.8e3
	bet=(1-(1+E/931.5)**(-2))**(0.5)
	gam=1.+E/931.5
	d2=(gam-1.)/(gam+1.)
	d1=sqrt(d2)
	gag=gam/(gam+1.)
	eta=alph/bet
	y1=max(0.,y1s-1.)
	y2=max(0.,yl-1.)
	y3=max(0.,ym-1.)
	y4=max(0.,yn-1.)
	y5=max(0.,yo-1.)
cc****************************
cc calcul du nombre d'electrons par couche de la cible (version modifiee)
	jm=3
	tk=1.
	if(zt.lt.29.) then
	tm=zt-10.
	else
	tm=18.
	endif
	if(zt.lt.11.) then
	tm=0.
	tl=zt-2.
	jm=2
	else
	tl=8.
	endif
	if(zt.lt.3) then
	tk=zt-1.
	tl=0.
	jm=1
	else
	tk=1.
	endif	
	oc(1)=2.
	if(zt.eq.1) then
	oc(1)=1.
c	tk=0.8
	endif
c..............................
	oc(2)=tl
	oc(3)=tm
	do 1 jp=1,30
	if(jp.eq.1) zps=zp-0.15*y1
	if(jp.eq.2) zps=(zp-0.85*y1s-0.35*y2)/2.
	if(jp.eq.3) zps=(zp-y1s-0.85*yl-0.35*y3)/3.
	if(jp.eq.4) zps=(zp-y1s-yl-0.85*ym-0.35*y4)/4.
	if(jp.eq.5) zps=(zp-y1s-yl-ym-0.85*y4)/5.
	if(jp.eq.6) zps=(zp-y1s-yl-ym-y4-0.85*y5)/6.
	if(jp.ge.7) zps=qm/jp
	do 2 jt=1,jm
	if(jt.eq.1) zts=zt-0.3*tk
	if(jt.eq.2) zts=(zt-1.7-0.35*tl)/2.
	if(jt.eq.3) zts=(zt-8.8-0.35*tm)/3.
	if(jt.eq.1.and.jp.eq.1) then
	Z2PFAC=1.
	else
	Z2PFAC=1.16
	endif
	if(zps.le.zts) then
	z1=zps
	z2=zts
	z2q=zts*z2pfac
	else
	z1=zts
	z2=zps
	z2q=zps*z2pfac
	endif
	EI=sqrt(1.-alph*alph*z2*z2)
	EF=sqrt(1.-alph*alph*z1*z1)
	PM=eta*(ef/gam-ei)/(alph*alph)
	ezp=eta*z2q
	ezp2=ezp*ezp
	pez=pi*ezp
	expez=exp(pez)
	exzt=-2.*ezp*atan(-pm/z2)
	zfac=z1*z2/(z2*z2+pm*pm)
	obkfac=s0*128*pi*eta*eta*zfac**5./(5.*gag*gam)
	eikfac=2.*pez*exp(exzt)/(expez-1./expez)
	eik=1.+5.*ezp*pm/(4.*z2)+5.*ezp2*pm*pm/(12.*z2*z2)+ezp2/6.
	mag=-d2+5.*d2*d2/16.+5.*d2*gag*z2q/(8.*z2)+d2*ezp2/4.
	mag=mag+5.*d2*d2*ezp2/48.
	orb=5.*pi*d1*alph*(z1+z2)*(1.-d2/2.)/18.
	orb=orb-5.*d1*alph*z2*ezp*(1-d2/2.)/8.
	orb=orb-5.*pi*d1*gag*alph*z1*z2q/(18.*z2)
	orb=orb+5.*pi*d1*gag*gag*alph*z1*z2q*z2q/(28.*z2*z2)
	orb=orb-5.*pi*d1*gag*alph*(z1+z2-d2*z1)*z2q/(28*z2)
	ceik(jp,jt)=jp*jp*oc(jt)*obkfac*eikfac*(eik+mag+orb)
2	continue
	SIG(jp)=ceik(jp,1)+ceik(jp,2)+ceik(jp,3)
1	continue
cc*****************************************************************
	ck=SIG(1)
	cl=SIG(2)
	cm=SIG(3)
	cn=SIG(4)
c not used	
	co=SIG(5)
	cp=SIG(6)
c ........	
	ct=0.
	do i=5,30
	ct=ct+SIG(i)
	end do
c*******formule de schlachter****************************
	E1=1000.*E/((Zt**1.25)*(qp**0.7))
	sig1=(1.1e-8/(E1**4.8))*(1-exp(-0.037*(E1**2.2)))
	sig1=sig1*(1-exp(-2.44e-5*(E1**2.6)))
	sigs=(sig1*(qp**0.5))/(Zt**1.8)  
cc*****************************************************************
      cs=sigs
	return
	end
cc*************************************************************
	subroutine REC(sk,sl1,slp)
cc			      Version du 3/11/93
cc*********************************************************************
cc	    Calcule les sections efficaces de REC 1s,2s et 2p
cc	     dans l'approximation type Bethe et Salpeter
cc	 (formules 71.7 p 304 et 71.14 71.15 p 306 ; voir aussi
cc	  formule 75.6 p 322)
cc*********************************************************************
	common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn
	common/don/zp,E,zt
	QM=Zp-(y1s+yl+ym+yn)
	pi=4.*atan(1.)
	Tel=E*511003.4/931.5
	bet=(1-(1+E/931.5)**(-2))**(0.5)
	gam=1.+E/931.5
	azk=(zp-0.15*y1s)/137.036
	azl=(zp-0.85*y1s-0.35*yl)/137.036
	azm=(zp-y1s-0.85*yl-0.35*ym)/137.036
	az=azk
	BK=-511003.4*((1.+(az/(1.-az**2.)**0.5)**2.)**(-0.5)-1.)
	az=azl
	BL1=-511003.4*((1.+(az/(1.+(1.-az**2.)**0.5))**2.)**(-0.5)-1.)
	az=azm
	BL3=-511003.4*((1.+(az/(4.-az**2.)**0.5)**2.)**(-0.5)-1.)
	IBK=BK+0.5
	IBL1=BL1+0.5
	IBL3=BL3+0.5
cc****************************
	EK=BK+Tel
	IEK=EK+0.5
	ZK=(BK/Tel)**0.5
	sphk=256.*pi*3.5e3/3.*BK**3./EK**4.*exp(-4.*zk*atan(1./zk))
	sphk=sphk/(1.-exp(-2.*pi*zk))
	EL1=BL1+Tel
	IEL1=EL1+0.5
	ZL1=(BL1/Tel)**0.5
	sphl1=2048.*pi*3.5e3/3.*BL1**3./EL1**4.*(1.+3.*bl1/el1)
	sphl1=sphl1*exp(-8.*zl1*atan(1./zl1))/(1.-exp(-4.*pi*zl1))
	sphl2=2048.*pi*3.5e3/9.*BL1**4./EL1**5.*(3.+8.*bl1/el1)
	sphl2=sphl2*exp(-8.*zl1*atan(1./zl1))/(1.-exp(-4.*pi*zl1))
	EL3=BL3+Tel
	IEL3=EL3+0.5
	ZL3=(BL3/Tel)**0.5
	sphl3=4096.*pi*3.5e3/9.*BL3**4./EL3**5.*(3.+8.*bl3/el3)
	sphl3=sphl3*exp(-8.*zl3*atan(1./zl3))/(1.-exp(-4.*pi*zl3))
	sk=zt*sphk*(ek/(gam*bet*511003.4))**2.
	sl1=zt*sphl1*(el1/(gam*bet*511003.4))**2.
	sl2=zt*sphl2*(el1/(gam*bet*511003.4))**2.
	sl3=zt*sphl3*(el3/(gam*bet*511003.4))**2.
	slp=sl2+sl3
	sl=sl1+slp
cc*****************************************************************
	return
	end
cc*************************************************************

