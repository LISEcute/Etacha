c********1*********1*********1*********1*********1*********1********1**
c Pour ETACHA4
c.....calcula la secc. eficaz total de excitacion desde los
c.....estados ns, np1 y np-1 al estado final nlm, por S.Eikonal
c.....Autor: Cesar A. Ramirez            Fecha: 10-12-2012
c         zt=carga del nucleo blanco
c         zp=carga del nucleo proyectil 
c         ei=energia del estado ligado inicial
c         ef=energia del estado ligado final
c         n,l,m = numeros cuanticos del estado final
c         eta=componente normal del momento transferido
c**********************************************************************
c----- Program to calculte total cross sections of  
c----- monoelectronic atoms excitation, for colision with nude projectiles 
c----- in the SE (Symmetric-Eikonal) aproximation.
c--------------------------------------------------
c......changes in the present version (JPR 02/2012)
c----- zp: Projectile charge (to be excited)
c----- zt: Target charge
c......SC and ASC corrections
c-------------------------------------------------
c----- (n0,l0,m0) quantum numbers of initial state 
c----- (n,l,m) quantum numbers of final state
c----- SUBROUTINES:
c----- * rns, rnp0, rnp1, rnp1m: excitation from ns, np0, np1 np-1 to n'lm
c-----                    respectively with n'lm arbitraries (n'>n)
c----- * intdef: for definite integrals from 0 to infinity
c----- * fac: factorial function
c----- * cg: Clebsch-Gordan coeficients
c----- * armonic: spherical harmonic function
c----- * plgndr: Legendre polinomials
c----- * gamac: gama function
c----- * cf21d, chyper, hypgfx: hypergeometric functions F21(ca,cb,cc,cx)
c**********************************************************************
      subroutine SEnlm(zp8,zt8,E8)
      implicit real*8(a-b,d-h,o-z),complex*16(c)
      real*4 SeSE,StSE
	character kk*2
c-----------------------------------------------------------------
      common/m1/zp,zt
      common/m2/epsi,v,pi
      common/m3/ca1,ca2,ca3,ca4,calfa
      common/m4/ci,cinu,cuno,cdos,con0
      common/m5/n0,l0,m0,n,l,m
	common/asc/sc,scoa,icor
	common/SecSE/SeSE,StSE
	dimension sexnlm(250)
	dimension signlm(4,2,3,5,5,6),signl(4,2,5,5),signt(4,2,5)
	dimension SeSE(66),StSE(12)
      external rns
      external rnp0
      external rnp1
      external rnp1m
c-----------------------------------------------------------------
	zp=zp8
	zt=zt8
	E=E8
	bet=dble((1.-(1.+E/931.5)**(-2.))**(0.5))
	v=137.036d0*bet
c----------------------------------------------------------------------
	scf=dble(1.13*zt**(1./3.))
	fs1=1.d0
	if (zt.eq.1.) fs1=1.23d0
	if (zt.eq.2.) fs1=1.526d0
	if (zt.eq.6.) fs1=0.78d0
	if (zt.eq.7.) fs1=0.85d0
	if (zt.eq.10.) fs1=1.04d0
	if (zt.eq.13.) fs1=0.58d0
	if (zt.eq.14.) fs1=0.59d0
	if (zt.eq.18.) fs1=0.68d0
	if (zt.eq.29.) fs1=0.672d0
	if (zt.eq.36.) fs1=0.61d0
	if (zt.eq.54.) fs1=0.535d0
	sc=scf*fs1
	icor=1
c----------------------------------------------------------------------
      pi=dacos(-1.d0) 
      xnu=zt/v
c
      ci=(0.d0,1.d0)
      cuno=(1.d0,0.d0)
      cdos=(2.d0,0.d0)
      cinu=ci*xnu
      ca1=1.d0+cinu
      ca2=1.d0-cinu
      ca3=1.d0+2.d0*cinu
      ca4=1.d0-2.d0*cinu
      senhy=dsinh(pi*xnu)
        test1=1.d-6
        test2=1.d-6
      call gamac(cinu,1,test1,test2,nt,cloga)
        cgama=cdexp(cloga)
c-----------------------------------------------------------------------       
c     open(40,file='tnlm4.dat',status='unknown')
c      write(40,200)zp,zt,e
c	if (icor.eq.1) then
c	write(40,*)' with SC+ASC corrections'
c	else
c	write(40,*)' no   SC+ASC corrections'
c	end if
c200   format(' zp=',d15.5,' zt=',d15.5,' E=',e15.5)
201   format(' initial=',3i3,'        final=',3i3)
c-----------------------------------------------------------------------
	index=1
	do 100 n0=1,4
	l0m=min(n0-1,1)
	do 100 l0=0,l0m
	do 100 m0=0,l0
	nfm=n0+1
	do 110 n=nfm,5
	do 120 l=0,n-1
	do 130 m=0,l
	write(*,201)n0,l0,m0,n,l,m
c----------------------------------------------------------------------
      n10=n0-l0-1
      n20=n0+l0
      fac10=fac(n10)
      fac20=fac(n20)
      xn0=2.d0/n0**2*dsqrt(zp)**3*dsqrt(fac20*fac10)*(2*zp/n0)**l0
      n1=n-l-1
      n2=n+l
      fac1=fac(n1)
      fac2=fac(n2)
      xnl=2.d0/n**2*dsqrt(zp)**3*dsqrt(fac2*fac1)*(2*zp/n)**l
c----------------------------------------------------------------------
      ei=-(zp/n0)**2/2.d0
      ef=-(zp/n)**2/2.d0
          epsi=(ef-ei)/v
          if(epsi.eq.0.d0)epsi=1.d-6
c.............................
	if (v.gt.1.75d0*epsi) then
	scoa=1.d0-(1.75d0*epsi/v)**2
	else
	scoa=0.0d0
	end if
c.............................
          calfa=2.d0*ci*v*epsi
          beta=(2.d0*v*epsi)**2
c----------------------------
      con0=4.d0*pi*v*(pi*xnu/cgama/senhy)**2/beta*xn0*xnl
c     *        /cdexp(cinu*dlog(beta))         cte. de módulo=1
c----------------------------
      ymin=0.d0
      ymax=3.d0
        a02=2.8d-17
        pre=1.d-3
C---------------------------------------------
      if(l0.eq.0)then
      call intdef(pre,ymin,ymax,rns,xint)
      tnlm=2*pi*xint*a02
      else
        if(m0.eq.0)then
          call intdef(pre,ymin,ymax,rnp0,xint)
          tnlm=2*pi*xint*a02
        else
          call intdef(pre,ymin,ymax,rnp1,xint)
          tp1=2*pi*xint*a02
          if(m.ne.0)tp1=2*tp1
          call intdef(pre,ymin,ymax,rnp1m,xint)
          tp1m=2*pi*xint*a02
          if(m.ne.0)tp1m=2*tp1m
          tnlm=tp1+tp1m
c
        end if
      end if
C--------------------------------------------------------------------
C---- Se multiplica por 2 la sección eficaz correspondiente a 
C---- la transición desde un estado inicial con m=0 
C---- a un estado final con m distinto a cero, para tener en cuenta los 
C---- dos casos simétricos +m y -m. 
C--------------------------------------------------------------------
      if(m0.eq.0.and.m.ne.0)tnlm=2*tnlm
c      write(40,300)n0,l0,m0,n,l,m,tnlm,index,sc,scoa
	sexnlm(index)=tnlm
	signlm(n0,l0+1,m0+1,n,l+1,m+1)=tnlm
	index=index+1
130	Continue
120	Continue
110	Continue
100   Continue
c     somme sur les m0 et m
c	write(40,*)' '
      index=1
      do 140 n0=1,4
	l0m=min(n0-1,1)
	nfm=n0+1
	do 140 l0=0,l0m
	do 140 n=nfm,5
	do 140 l=0,n-1
	sigt=0.0
	do 150 m0=0,l0m
	do 150 m=0,l
	sigt=sigt+signlm(n0,l0+1,m0+1,n,l+1,m+1)
150   continue
	if(l0.eq.1) sigt=sigt/3.
      signl(n0,l0+1,n,l+1)=sigt
      SeSE(index)=sigt
      index=index+1
c      write(40,310)n0,l0,n,l,sigt
140   continue
c	write(40,*)' '
c     somme sur les l pour n=4 et 5
      index=1
      do 160 n=4,5
	nim=n-1
      do 160 n0=1,nim
	l0m=min(n0-1,1)
	do 160 l0=0,l0m
      sigtn=0.
      do 170 l=0,n-1
      sigtn=sigtn+signl(n0,l0+1,n,l+1)
170   continue
      signt(n0,l0+1,n)=sigtn
      StSE(index)=sigtn
      index=index+1
c      write(40,311)n0,l0,n,sigtn
160   continue     	
c      close(unit=40)
c-----
300   format(' inicial=',3i3,'   final=',3i3,'  tot = ',d15.5,
     s '   index=',i3,' coef ASC=',2d15.5)
310   format(' inicial=',2i3,'      final=',2i3,'     tot = ',d15.5)
311   format(' inicial=',2i3,'      final=',i3,'        tot = ',d15.5)
      return
      end
c*************************************** fin  programa principal ******
c
c     SE    ns-nlm.for
c
c***********************************************************************
c.....function=eta * /R(eta)/^2      en la aprox. SE  
c.....desde ns a nlm
c***********************************************************************
      subroutine rns(eta,xf)
      implicit real*8(a-b,d-h,o-z), complex*16(c)
      common/m1/zp,zt
      common/m2/epsi,v,pi
      common/m3/ca1,ca2,ca3,ca4,calfa
      common/m4/ci,cinu,cuno,cdos,con0
      common/m5/n0,l0,m0,n,l,m
	common/asc/sc,scoa,icor
c-----------------------------------------------------------------------
      cero=(0.d0,0.d0)
      gama=eta**2+epsi**2
      xk=dsqrt(gama)
      znz=zp/n+zp/n0
      vari=-(eta/epsi)**2
      call cf21d(cinu,cinu,cdos,vari,chyp1)
      call cf21d(ca1,ca1,cdos,vari,chyp2)
      chyp1=(ca2*chyp1-cinu*(1.d0-vari)*chyp2)/ca4
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      csum2=cero
      csum3=cero
      if(l.eq.0)then
         minla=0
      else
         minla=l-1
      end if
         maxla=l+1
c-----
      do la=minla,maxla
        ga=dsqrt(pi)
        do j=0,la
        ga=ga*(0.5d0+j)
        end do
      xla=dsqrt(pi)*xk**la/2.d0**(la+1)/ga
        sumq1=cero
        sumq=cero 
        maxiq=n-l-1
      do iq=0,maxiq
         n4=n-l-1-iq
         n5=iq
         n6=2*l+1+iq
      fac4=fac(n4)
      fac5=fac(n5)
      fac6=fac(n6)
      xq= (-2.d0*zp/n)**iq/fac4/fac5/fac6
        sump1=cero 
        sump=cero 
        maxip=n0-1
      do ip=0,maxip
        n40=n0-1-ip
        n50=ip
        n60=1+ip
      fac40=fac(n40)
      fac50=fac(n50)
      fac60=fac(n60)
      xp = (-2.d0*zp/n0)**ip/fac40/fac50/fac60
c-------------------------------------------------
      xmu=2.5d0+l+iq+ip
      xnu=la+0.5d0
      xa=(xmu+xnu)/2.d0
      xc=xnu+1.d0
      xz=-gama/znz**2
      xc1=2.d0*xa
      xc2=2.d0*xa-xc+1.d0
      yz=(1.d0-dsqrt(1.d0-xz))/(1.d0+dsqrt(1.d0-xz))*cuno 
      n3 = 3+l+la+iq+ip
      h1=(2.d0/(1.d0+dsqrt(1.d0-xz)))**n3 
      call hygfx(xc1,xc2,xc,yz,hyp31)
       nn2 = 2+l+la+iq+ip
       fac3 = fac(nn2)
       xintx1 = xla*fac3/znz**n3*h1*hyp31
      if(la.eq.l)then
      sump1 = sump1 + xp*xintx1
      else
       xc1=2.d0*xa-1.d0 
       xc2=2.d0*xa-xc 
       h2=(2.d0/(1.d0+dsqrt(1.d0-xz)))**(n3-1)
      call hygfx(xc1,xc2,xc,yz,hyp32) 
       xintx2 = xla*fac3/znz**n3*h2*hyp32
      sump = sump + xp*(-zp/n0*xintx1 + znz*ip/(2.d0*xa-1.d0)*xintx2)
      endif
      end do
      if(la.eq.l)then
        sumq1 = sumq1+sump1*xq
      else
        sumq = sumq+sump*xq
      endif
      end do
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c              s=-1      indica que es un armonico conjugado
      s=-1.d0
      xkz=-epsi
      tk=dacos(xkz/xk)
      calga=calfa/gama
      xl=dsqrt((2.d0*la+1.d0)/(2.d0*l+1.d0)/4.d0/pi)
      if(la.eq.l)then 
        call armonic(s,l,m,tk,0.d0,cylam)
        csum1= ci**la*xl*sumq1*cylam*calfa*(ca2/cinu*chyp1+chyp2)
      else
        call cg(la,1,l,0,0,0,xc3)
        call cg(la,1,l,m+1,-1,m,xc4)
        call cg(la,1,l,m-1,1,m,xc5)
        call cg(la,1,l,m,0,m,xc6)
        call armonic(s,la,m+1,tk,0.d0,cylam1) 
        call armonic(s,la,m-1,tk,0.d0,cylam2) 
        call armonic(s,la,m,tk,0.d0,cylam3) 
      csum2=csum2+ci**la*xl*sumq*xc3*(xc4*cylam1-xc5*cylam2)/dsqrt(2.d0) 
     *     *ci*eta*calga*ca2/cinu*chyp1 
      csum3=csum3 + ci**la*xl*sumq*xc3*xc6*cylam3 
     *     *(2.d0*v*chyp2-ci*epsi*calga*ca2/cinu*chyp1) 
      end if
      end do
c-----------------------------------------------------------------------
      csum= 4.d0*pi* (csum1 - 2.d0 * (csum2+csum3))
c-----------------------------------------------------------------------
      cf= ci/(2.d0*pi*v) * con0 * csum  
c     *         * cdexp(2*cinu*dlog(gama))    cte de modulo=1
      xf=eta*cdabs(cf)*cdabs(cf)
c...............................	
c	SC and ASC corrections
	if (icor.eq.1) then
	sq=(1.d0+(sc/xk)**2)**(-2)+scoa*(1.d0-(1.d0+(xk/sc)**2)**(-2))/zt
	xf=xf*sq
	end if
c...............................
      return
      end            
c***********************************************************************
c
c   SE  np0-nlm.for
c
c***********************************************************************
c.....function=eta * /R(eta)/^2      en la aprox. SE  
c.....desde np0 a nlm
c***********************************************************************
      subroutine rnp0(eta,xf) 
      implicit real*8(a-b,d-h,o-z),complex*16(c)
      common/m1/zp,zt
      common/m2/epsi,v,pi
      common/m3/ca1,ca2,ca3,ca4,calfa
      common/m4/ci,cinu,cuno,cdos,con0
      common/m5/n0,l0,m0,n,l,m
	common/asc/sc,scoa,icor
c-----------------------------------------------------------------------
      cero=(0.d0,0.d0)
      gama=eta**2+epsi**2
      xk=dsqrt(gama)
      znz=zp/n+zp/n0
c-----
      vari=-eta**2/epsi**2
      call cf21d(cinu,cinu,cdos,vari,chyp1)
      call cf21d(ca1,ca1,cdos,vari,chyp2)
      chyp1=(ca2*chyp1-cinu*(1.d0-vari)*chyp2)/ca4
      csum1=cero
      csum2=cero
      csum3=cero
        minla2=iabs(l-2)
        minla1=iabs(l-1)
      minla=min(minla1,minla2)
        if(l.eq.0)minla=0
        maxla=l+2
      do la=minla,maxla
      sumq1=0.d0 
      sumq2=0.d0 
      sumq=0.d0 
         maxiq=n-l-1
      do iq=0,maxiq
         n4=n-l-1-iq
         n5=iq
         n6=2*l+1+iq
      fac4=fac(n4)
      fac5=fac(n5)
      fac6=fac(n6)
      xq= (-2.d0*zp/n)**iq/fac4/fac5/fac6
      sump1=0.d0 
      sump2=0.d0 
      sump=0.d0 
         maxip=n0-l0-1
      do ip=0,maxip
         n40=n0-l0-1-ip
         n50=ip
         n60=2*l0+1+ip
      fac40=fac(n40)
      fac50=fac(n50)
      fac60=fac(n60)
      xp = (-2*zp/n0)**ip/fac40/fac50/fac60
c--------------------------------------------------
      nn3 = 3+l+l0+la+iq+ip
      nn2 = 2+l+l0+la+iq+ip
      fn2 =fac(nn2)
      xmu=  2.5d0+l+l0+iq+ip
      xnu=  la+0.5d0
      xc=la+1.5d0
      xxa=(xnu+xmu)/2
      xxb=(xnu-xmu)/2+0.5d0
      xxz=gama/(gama+znz**2)
      call hygfx(xxa,xxb,xc,xxz,hyp31)
      hyp31=hyp31/dsqrt(gama+znz**2)**nn3
      if(la.eq.l+1.or.la.eq.iabs(l-1))then
      sump1 = sump1 + xp*fn2*hyp31
      else
       xxa=xxa-0.5d0
       xxb=xxb+0.5d0 
      call hygfx(xxa,xxb,xc,xxz,hyp32)
      hyp32=hyp32/dsqrt(gama+znz**2)**nn2
      sump  = sump  + xp*fn2*( ip*hyp32/nn2 - zp/n0*hyp31 )
       if(la.eq.l)then
      sump2 = sump2 + xp*fn2*( (ip+3)*hyp32/nn2 - zp/n0*hyp31 )
       endif
      endif
       end do
      if(la.eq.l+1.or.la.eq.iabs(l-1))then
      sumq1 = sumq1 + xq*sump1
      else
      sumq  = sumq + xq*sump
       if(la.eq.l)then
      sumq2=sumq2+xq*sump2
       endif
      end if
      end do        
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c                           s=-1 (indica el conjugado)
      s=-1.d0
        ga=dsqrt(pi)
        do j=0,la
        ga=ga*(0.5d0+j)
        end do
      xkz=-epsi
      tk=dacos(xkz/xk)
      calga=calfa/gama
      xla=dsqrt(pi)*xk**la/2.d0**(la+1)/ga
     *                 *dsqrt((2*la+1.d0)/(2*l+1.d0)/4.d0/pi)
      if(la.eq.l+1.or.la.eq.iabs(l-1))then
         call cg(la,1,l,0,0,0,xc1)
         call cg(la,1,l,m,0,m,xc2)
         call armonic(s,la,m,tk,0.d0,cylam1) 
      csum1=csum1 + ci**la*xla*dsqrt(3.d0)*sumq1*xc1*xc2*cylam1
      else
         call cg(la,2,l,0,0,0,xc3)
         call cg(la,2,l,m,0,m,xc4)
         call cg(la,2,l,m+1,-1,m,xc5)
         call cg(la,2,l,m-1,1,m,xc6)
         call armonic(s,la,m,tk,0.d0,cylam)
         call armonic(s,la,m+1,tk,0.d0,cylam2)
         call armonic(s,la,m-1,tk,0.d0,cylam3)
      csum2=csum2 + ci**la*xla*sumq*xc3
     *                      *(xc5*cylam2-xc6*cylam3)/dsqrt(2.d0)
      csum3=csum3 + ci**la*xla*sumq* xc3*xc4*cylam*2.d0/dsqrt(3.d0)
        if(la.eq.l)then
      csum4=ci**l*xla*sumq2/dsqrt(3.d0)*cylam
        endif
      end if
      end do
      csum1=csum1*calfa*(ca2/cinu*chyp1+chyp2)
      csum2=csum2*ci*eta*calga*ca2/cinu*chyp1
      csum3=csum3*(2*v*chyp2-ci*epsi*calga*ca2/cinu*chyp1)
      csum4=csum4*(2*v*chyp2-ci*epsi*calga*ca2/cinu*chyp1)
c-----------------------------------------------------------------------
      csum= 4*pi*(  csum1 - 2 *(csum2+csum3+csum4) )
c-----------------------------------------------------------------------
      cf= ci/(2.d0*pi*v) * con0 * csum 
c     *          * cdexp(2.d0*cinu*dlog(gama))      cte de modulo=1
      xf=eta * cdabs(cf)*cdabs(cf)
c...............................	
c	SC and ASC corrections
	if (icor.eq.1) then
	sq=(1.d0+(sc/xk)**2)**(-2)+scoa*(1.d0-(1.d0+(xk/sc)**2)**(-2))/zt
	xf=xf*sq
	end if
c...............................
      return
      end
c***********************************************************************
c
c     SE    np1-nlm.for
c
c***********************************************************************
c.....function=eta * /R(eta)/^2      en la aprox. SE  
c.....desde np1 a nlm
c***********************************************************************
      subroutine rnp1(eta,xf) 
      implicit real*8(a-b,d-h,o-z),complex*16(c)
      common/m1/zp,zt
      common/m2/epsi,v,pi
      common/m3/ca1,ca2,ca3,ca4,calfa
      common/m4/ci,cinu,cuno,cdos,con0
      common/m5/n0,l0,m0,n,l,m
	common/asc/sc,scoa,icor
c-----------------------------------------------------------------------
      cero=(0.d0,0.d0)
      gama=eta**2+epsi**2
      xk=dsqrt(gama)
      znz=zp/n+zp/n0
c-----
      vari=-eta**2/epsi**2
      call cf21d(cinu,cinu,cdos,vari,chyp1)
      call cf21d(ca1,ca1,cdos,vari,chyp2)
      chyp1=(ca2*chyp1-cinu*(1.d0-vari)*chyp2)/ca4
c-----------------------------------------------------------------------
      csum1=cero
      csum2=cero
      csum3=cero
        minla2=iabs(l-2)
        minla1=iabs(l-1)
      minla=min(minla1,minla2)
        if(l.eq.0)minla=0
        maxla=l+2
c-----
      do la=minla,maxla
      sumq1=0.d0 
      sumq2=0.d0 
      sumq=0.d0 
         maxiq=n-l-1
      do iq=0,maxiq
         n4=n-l-1-iq
         n5=iq
         n6=2*l+1+iq
      fac4=fac(n4)
      fac5=fac(n5)
      fac6=fac(n6)
      xq= (-2.d0*zp/n)**iq/fac4/fac5/fac6

      sump1=0.d0 
      sump2=0.d0 
      sump=0.d0 

         maxip=n0-l0-1
      do ip=0,maxip
         n40=n0-l0-1-ip
         n50=ip
         n60=2*l0+1+ip
      fac40=fac(n40)
      fac50=fac(n50)
      fac60=fac(n60)
      xp = (-2*zp/n0)**ip/fac40/fac50/fac60
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          nn3 = 3+l+l0+la+iq+ip
          nn2 = 2+l+l0+la+iq+ip
          fn2 =fac(nn2)
         xmu= 2.5d0+l+l0+iq+ip
         xnu= la+0.5d0
         xa= (xnu+xmu)/2
         xb= xa+0.5d0
         xc= xnu+1.d0
         xxb= xc-xb
         xxz= gama/(gama+znz**2)
      call hygfx(xa,xxb,xc,xxz,hyp31)
      hyp31=hyp31/dsqrt(gama+znz**2)**nn3
      if(la.eq.l+1.or.la.eq.iabs(l-1))then
      sump1 = sump1 + xp*fn2*hyp31
      else
      xxa=xa-0.5d0
      xxb=xc-xa 
      call hygfx(xxa,xxb,xc,xxz,hyp32) 
      hyp32=hyp32/dsqrt(gama+znz**2)**nn2
      sump = sump + xp*fn2*( ip*hyp32/nn2 - zp/n0*hyp31 )
        if(la.eq.l)then
      sump2 = sump2 + xp*fn2*((ip+3)*hyp32/nn2 - zp/n0*hyp31 )
        endif
      endif
      end do

      if(la.eq.l+1.or.la.eq.iabs(l-1))then
      sumq1 = sumq1 + xq*sump1
      else
      sumq  = sumq + xq*sump
          if(la.eq.l)then
      sumq2=sumq2+xq*sump2
          endif
      end if
      end do        
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c                           s=-1 (indica el conjugado)
      s=-1.d0
        ga=dsqrt(pi)
        do j=0,la
        ga=ga*(0.5d0+j)
        end do
      xkz=-epsi
      tk=dacos(xkz/xk)
      calga=calfa/gama
      xla=dsqrt(pi)*xk**la/2.d0**(la+1)/ga
     *           *dsqrt((2*la+1.d0)/(2*l+1.d0)/4.d0/pi)
      if(la.eq.l+1.or.la.eq.iabs(l-1))then
         call cg(la,1,l,0,0,0,xc1)
         call cg(la,1,l,m-1,1,m,xc2)
         call armonic(s,la,m-1,tk,0.d0,cylam1) 
      csum1=csum1 + ci**la*xla*dsqrt(3.d0)*sumq1*xc1*xc2*cylam1
      else
         call cg(la,2,l,0,0,0,xc3)
         call cg(la,2,l,m,0,m,xc4)
         call cg(la,2,l,m-2,2,m,xc5)
         call cg(la,2,l,m-1,1,m,xc6)
         call armonic(s,la,m,tk,0.d0,cylam)
         call armonic(s,la,m-2,tk,0.d0,cylam2)
         call armonic(s,la,m-1,tk,0.d0,cylam3)
      csum2=csum2 + ci**la*xla*sumq*xc3/dsqrt(6.d0)
     *               *(xc4*cylam-dsqrt(6.d0)*xc5*cylam2)
      csum3=csum3 + ci**la*xla*sumq*xc3*xc6*cylam3
            if(la.eq.l)then
      csum4=-ci**l*xla*sumq2/dsqrt(6.d0)*cylam
            endif
      end if
         end do
      csum1=csum1*calfa*(ca2/cinu*chyp1+chyp2)
      csum2=csum2*ci*eta*calga*ca2/cinu*chyp1
      csum3=csum3*(2*v*chyp2-ci*epsi*calga*ca2/cinu*chyp1)
      csum4=csum4*ci*eta*calga*ca2/cinu*chyp1
c-----------------------------------------------------------------------
      csum= 4*pi* ( csum1 - 2 *(csum2+csum3+csum4)  )
c-----------------------------------------------------------------------
      cf= ci/(2.d0*pi*v) * con0 * csum
c     *           * cdexp(2.d0*cinu*dlog(gama))   cte. de modulo=1
      xf=eta*cdabs(cf)*cdabs(cf)
c...............................	
c	SC and ASC corrections
	if (icor.eq.1) then
	sq=(1.d0+(sc/xk)**2)**(-2)+scoa*(1.d0-(1.d0+(xk/sc)**2)**(-2))/zt
	xf=xf*sq
	end if
c...............................
      return
      end
c***********************************************************************
c
c     SE    np1m-nlm.for
c
c***********************************************************************
c.....function=eta * /R(eta)/^2      en la aprox. SE  
c.....desde np-1 a nlm 
c***********************************************************************
      subroutine rnp1m(eta,xf) 
      implicit real*8(a-b,d-h,o-z),complex*16(c)
      common/m1/zp,zt
      common/m2/epsi,v,pi
      common/m3/ca1,ca2,ca3,ca4,calfa
      common/m4/ci,cinu,cuno,cdos,con0
      common/m5/n0,l0,m0,n,l,m
	common/asc/sc,scoa,icor
c-----------------------------------------------------------------------
      cero=(0.d0,0.d0)
      gama=eta**2+epsi**2
      xk=dsqrt(gama)
      znz=zp/n+zp/n0
c-----
      vari=-eta**2/epsi**2
      call cf21d(cinu,cinu,cdos,vari,chyp1)
      call cf21d(ca1,ca1,cdos,vari,chyp2)
      chyp1=(ca2*chyp1-cinu*(1.d0-vari)*chyp2)/ca4
c-----------------------------------------------------------------------
      csum1=cero
      csum2=cero
      csum3=cero
        minla2=iabs(l-2)
        minla1=iabs(l-1)
      minla=min(minla1,minla2)
        if(l.eq.0)minla=0
        maxla=l+2
c-----
      do la=minla,maxla

      sumq1=0.d0 
      sumq2=0.d0 
      sumq=0.d0 

         maxiq=n-l-1
      do iq=0,maxiq
         n4=n-l-1-iq
         n5=iq
         n6=2*l+1+iq
      fac4=fac(n4)
      fac5=fac(n5)
      fac6=fac(n6)
      xq= (-2.d0*zp/n)**iq/fac4/fac5/fac6

      sump1=0.d0 
      sump2=0.d0 
      sump=0.d0 

         maxip=n0-l0-1
      do ip=0,maxip
         n40=n0-l0-1-ip
         n50=ip
         n60=2*l0+1+ip
      fac40=fac(n40)
      fac50=fac(n50)
      fac60=fac(n60)
      xp = (-2*zp/n0)**ip/fac40/fac50/fac60
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          nn3 = 3+l+l0+la+iq+ip
          nn2 = 2+l+l0+la+iq+ip
          fn2 =fac(nn2)
         xmu= 2.5d0+l+l0+iq+ip
         xnu= la+0.5d0
         xa= (xnu+xmu)/2
         xb= xa+0.5d0
         xc= xnu+1.d0
         xxb= xc-xb
         xxz= gama/(gama+znz**2)
      call hygfx(xa,xxb,xc,xxz,hyp31)
      hyp31=hyp31/dsqrt(gama+znz**2)**nn3
      if(la.eq.l+1.or.la.eq.iabs(l-1))then
      sump1 = sump1 + xp*fn2*hyp31
      else
      xxa=xa-0.5d0
      xxb=xc-xa 
      call hygfx(xxa,xxb,xc,xxz,hyp32) 
      hyp32=hyp32/dsqrt(gama+znz**2)**nn2
      sump = sump + xp*fn2*( ip*hyp32/nn2 - zp/n0*hyp31 )
        if(la.eq.l)then
      sump2 = sump2 + xp*fn2*((ip+3)*hyp32/nn2 - zp/n0*hyp31 )
        endif
      endif
      end do

      if(la.eq.l+1.or.la.eq.iabs(l-1))then
      sumq1 = sumq1 + xq*sump1
      else
      sumq  = sumq + xq*sump
          if(la.eq.l)then
      sumq2=sumq2+xq*sump2
          endif
      end if
      end do        
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c                           s=-1 (indica el conjugado)
      s=-1.d0
        ga=dsqrt(pi)
        do j=0,la
        ga=ga*(0.5d0+j)
        end do
      xkz=-epsi
      tk=dacos(xkz/xk)
      calga=calfa/gama
      xla=dsqrt(pi)*xk**la/2.d0**(la+1)/ga
     *           *dsqrt((2*la+1.d0)/(2*l+1.d0)/4.d0/pi)
      if(la.eq.l+1.or.la.eq.iabs(l-1))then
         call cg(la,1,l,0,0,0,xc1)
         call cg(la,1,l,m+1,-1,m,xc2)
         call armonic(s,la,m+1,tk,0.d0,cylam1) 
      csum1=csum1 + ci**la*xla*dsqrt(3.d0)*sumq1*xc1*xc2*cylam1
      else
         call cg(la,2,l,0,0,0,xc3)
         call cg(la,2,l,m,0,m,xc4)
         call cg(la,2,l,m+2,-2,m,xc5)
         call cg(la,2,l,m+1,-1,m,xc6)
         call armonic(s,la,m,tk,0.d0,cylam)
         call armonic(s,la,m+2,tk,0.d0,cylam2)
         call armonic(s,la,m+1,tk,0.d0,cylam3)
      csum2=csum2 + ci**la*xla*sumq*xc3
     *               *(-xc4*cylam/dsqrt(6.d0)+xc5*cylam2)
      csum3=csum3 + ci**la*xla*sumq*xc3*xc6*cylam3
         if(la.eq.l)then
      csum4=ci**l*xla*sumq2/dsqrt(6.d0)*cylam
         endif
      end if
         end do
      csum1=csum1*calfa*(ca2/cinu*chyp1+chyp2)
      csum2=csum2*ci*eta*calga*ca2/cinu*chyp1
      csum3=csum3*(2*v*chyp2-ci*epsi*calga*ca2/cinu*chyp1)
      csum4=csum4*ci*eta*calga*ca2/cinu*chyp1
c-----------------------------------------------------------------------
      csum= 4*pi*( csum1 - 2 *(csum2+csum3+csum4)  )
c-----------------------------------------------------------------------
      cf= ci/(2.d0*pi*v) * con0 * csum
c     *   * cdexp(2.d0*cinu*dlog(gama))  factor de módulo=1
      xf=eta* cdabs(cf)*cdabs(cf)
c...............................	
c	SC and ASC corrections
	if (icor.eq.1) then
	sq=(1.d0+(sc/xk)**2)**(-2)+scoa*(1.d0-(1.d0+(xk/sc)**2)**(-2))/zt
	xf=xf*sq
	end if
c...............................
      return
      end
c***********************************************************************
c
c     Subrutinas.for
c
c***********************************************************************
c     Armonicos esfericos 
c     Autor: Cesar Ramirez          fecha: 10-4-2000
c     l,/m/:   numeros enteros
c     tk     angulo teta
c     s  +1: real  -1: imaginario
c     ma  valor absoluto de m
c-----------------------------------------------------------------------
      subroutine armonic(s,l,m,tk,fik,cylm)
c     ------------------------------------------------------------------
      implicit real*8(a-b,d-h,o-z),complex*16(c)
        ci=(0.d0,1.d0)
        pi=dacos(-1.d0)
        ma=iabs(m)
      cylm=(0.d0,0.d0)
      if(ma.gt.l)return
        n1=l-ma
        n2=l+ma
        fac11=fac(n1)
        fac12=fac(n2)
      x1= dsqrt((2.d0*l+1.d0)/(4.d0*pi)*fac11/fac12)
      x=dcos(tk)
        xlm=plgndr(l,ma,x)
      if(m.ge.0)then
       cylm=x1*cdexp(s*ci*ma*fik)*xlm
      else
       cylm=(-1.d0)**ma*x1*cdexp(-s*ci*ma*fik)*xlm
      end if
        return
      end
c***********************************************************************
c     Legendre polynomials (Numerical Recipes Software: www.nr.com)
c     x=cos(t)
c     ------------------------------------------------------------------
      FUNCTION PLGNDR(L,M,X)
      real*8 x,plgndr
      IF(M.LT.0.OR.M.GT.L.OR.dABS(X).GT.1.)then
      write(*,*)' m<0,   or   m>l,   or   dabs(x)>1'
      pause 'bad arguments'
      endif
      PMM=1.
      IF(M.GT.0) THEN
        SOMX2=SQRT((1.-X)*(1.+X))
        FACT=1.
        DO 11 I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.
11      CONTINUE
      ENDIF
      IF(L.EQ.M) THEN
        PLGNDR=PMM
      ELSE
        PMMP1=X*(2*M+1)*PMM
        IF(L.EQ.M+1) THEN
          PLGNDR=PMMP1
        ELSE
          DO 12 LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          PLGNDR=PLL
        ENDIF
      ENDIF
c      write(*,*)x,plgndr
      RETURN
      END
c***********************************************************************
      subroutine gamac(za,ny,test1,test2,nt,logam)                      
c     ******************************************************************
c          auteur: r. gayet              date inconnue                  
c          version 2    (modif a. salin)      30/10/81                  
c     log de gamma(z) ecrit si ny=2                                     
c     test1 est la difference minimum admise entre reel(z) et un pole   
c     pour distinguer reel(z) du pole                                   
c     test2 est la valeur minimum admise pour imag(z) quand la          
c     difference entre reel(z) et un pole est inferieure a test1        
c     nt=2 si z est un pole de gamma:dans ce cas logam=0 arbitrairement-
c     sinon nt=1                                                        
      complex*16 z,z2,logam,c,za                                        
      real*8 a(2),pi,d,test1,test2,test10                               
      equivalence(z,a(1))                                               
      data pi/.918938533204673d0/                                       
      z=za                                                              
      logam=(0.d0,0.d0)                                                 
      irz=idint(a(1))                                                   
      test10=test1*2.d0                                                 
      if(a(1)-test10)2,2,1                                              
    2 d=dabs(a(1)-dfloat(irz))                                          
      if(d.gt.test1 ) go to 1                                           
      if(dabs(a(2)).gt.test2)  go to 1                                  
      write(6,3) d,z                                                    
    3 format(//,4x,'reel(z)=',1pd9.2,'+',1pd22.15,/,4x,'imag(z)=',1pd22.
     115,//,4x,'z est considere comme un pole de la fonction gamma')    
      nt=2                                                              
      go to 100                                                         
    1 continue                                                          
      nt=1                                                              
      ng=10-irz                                                         
      if(ng)4,4,5                                                       
    5 c=dcmplx(dfloat(ng),0.d0)                                         
      z=z+c                                                             
    4 z2=z*z                                                            
      d=cdabs(z)                                                        
      if(d.ge.1.d+02) go to 6                                           
      logam=1.d0/156.d0/z2-691.d0/360360.d0                             
      logam=logam/z2+1.d0/1188.d0                                       
      logam=logam/z2-1.d0/1680.d0                                       
    6 if(d.ge.1.d+04)go to 7                                            
      logam=logam/z2+1.d0/1260.d0                                       
      logam=logam/z2-1.d0/360.d0                                        
    7 if(d.ge.1.d+07) go to 8                                           
      logam=(logam/z2+1.d0/12.d0)/z                                     
    8 logam=logam+pi-z+(z-0.5d0)*cdlog(z)                               
      if(ng)100,100,9                                                   
    9 c=(1.d0,0.d0)                                                     
      do 10 i=1,ng                                                      
      z=z-c                                                             
      logam=logam-cdlog(z)                                              
   10 continue                                                          
  100 continue                                                          
      if(ny.eq.2) write(ny,200) z,logam
  200 format(//,2x,'log de gamma(',2(1pd22.15),') = ',2(1pd22.15))      
      return                                                            
      end                                                               
c***********************************************************************
c***********************************************************************
c  Rutina para calcular el factorial de n
c  n! = n(n-1)(n-2)...2 ,
c  0! = 1
c----------------------------------------
      function fac(n)
      implicit none
      real*8 xn,xfac,fac,pi
      integer n,j
c
      pi=dacos(-1.d0)
      xn=n*1.d0
      if(n.eq.0.or.n.eq.1)then
          fac=1.d0
      else
          if(n.gt.90)then
          fac=dsqrt(pi*(2*xn+1/3.d0))*xn**xn/dexp(xn)
          else
          xfac=0.d0
          do j=2,n
          xfac=xfac+dlog(dfloat(j))
          end do
          fac=dexp(xfac)
          end if
      end if
      return
      end
c*********************************************************************
c**************************** l o g f a c ******************************   
c     calculates ln((n-1)!).
c     beware: fac(n) = ln((n-1)!) = ln(gamma(n))
c-----------------------------------------------------------------------
      block data faclog
      parameter (idf=140)                                               
      implicit real*8 (a-h,o-z)
      logical init                                                      
      common/fact/fac(idf),init                                         
      data init/.false./                                                
      end                                                               
c-----------------------------------------------------------------------
      subroutine logfac                                                 
      parameter (idf=140)                                               
      implicit real*8 (a-h,o-z)
      logical init                                                      
      common/fact/fac(idf),init                                         
      init=.true.                                                       
      fac(1)=0.d0                                                       
      do i=2,idf                                                     
        fac(i)=fac(i-1)+dlog(dfloat(i-1))                               
      end do
      end                                                               
c*********************************************************************
c.......... cg Calcula los coeficientes de Clebsch-Gordan
c.......... Debe llamarse antes al logfac.
c.......... logcle:  simbolos 3j
c.......... xcg:     coef de C-G
c***********************************************************************
      Subroutine cg(l1,l2,l3,m1,m2,m3,z)   
      implicit real*8 (a-h,o-z)
      call logfac
      call logcle(l1,l2,l3,m1,m2,-m3,x3j) 
      z=(-1.d0)**iabs(l1-l2-m3)*dsqrt(2.d0*l3+1.d0)*x3j
      end
c**************************** l o g c l e ******************************      
c         author:  a.salin     version 2  23/2/77 - 5/2/96
c
c     calculation of WIGNER's 3j  (definition of MESSIAH)
c     restriction: INTEGER MOMENTS ONLY.
c
c     before the first run of logcle, a call to logfac should be
c     performed. logfac calculates ln((n-1)!) for n going from 1 to idf
c     and stores the results in the common fact. 
c     for large values of the l's, increase idf. 
c        input        
c                l1, l2, l3, m1, m2, m3: momenta and their components. 
c        output
c                z: 3j coefficient.
c***********************************************************************
      subroutine logcle(l1,l2,l3,m1,m2,m3,z)
      implicit real*8 (a-h,o-z)
      parameter(idf=140) 
      logical init 
      dimension ac(idf) 
      common/fact/fac(idf),init  
c     
      if(.not.init) then 
        write(*,98) 
 98     format(1x,' logfac not called before logcle')  
        stop  
      end if
c
      i4=l1+l2+l3+2 
      if(i4.gt.idf) then 
        write(*,99) 
 99     format(1x,'logcle: value of idf too small')  
        stop
      end if
c
      if(m1+m2+m3.ne.0) then 
        z=0.d0 
        return  
      end if  
c    
      izmax=min(l1+l2-l3,l1-m1,l2+m2)+1 
      izmin=max(0,l2-l3-m1,l1+m2-l3)+1 
      if(izmax.lt.izmin) then 
        z=0.d0   
        return 
      end if  
c
      i1=l1+l2-l3+1 
      i2=l1-m1+1 
      i3=l2+m2+1 
      abra=0.5d0*(fac(i1)+fac(l3+l1-l2+1)+fac(l3+l2-l1+1)-fac(i4)
     1     +fac(l1+m1+1)+fac(i2)+fac(i3)+fac(l2-m2+1)+fac(l3+m3+1)
     2     +fac(l3-m3+1))
      k1=l3-l2+m1+1
      k2=l3-l1-m2+1 
      gros=250.d0 
      do 8 ii=izmin,izmax
        i=ii-1 
        ac(ii)=fac(i+1)+fac(i1-i)+fac(i2-i)
     1         +fac(i3-i)+fac(k1+i)+fac(k2+i)
        if(ac(ii).lt.gros)  gros=ac(ii)
 8    continue 
      accu=0.d0
      sig=(-1.d0)**izmin 
      do 9 ii=izmin,izmax
        sig=-sig 
        ac(ii)=ac(ii)-gros 
        accu=accu+sig*dexp(-ac(ii))  
 9    continue 
      z=(-1.d0)**iabs(l1-l2-m3)*dexp(abra-gros)*accu  
      return  
      end 
c***********************************************************************
      subroutine cf21d(ca,cb,cc,x,cf)
**************************************
*  Program: hyper.f
*  Version: 1.6 (see readme_hyper) 24/03/1997
*  Author: Pablof (pablof@cab.cnea.edu.ar)
*  Hypergeometric function for CDW-EIS calculations. 
*  Double precision, see comment for V 1.4
*  We assume that x is REAL
*  ~~~~~~~~~~~~~~~~~~~~~~~~
**************************************
      implicit none
      complex*16 ca,cb,cc,cf
      real*8 x
      complex*16 calgam,c0,c1
      parameter(c0=(0.d0,0.d0), c1=(1.d0,0.d0))

      if(ca.eq.c0.or.cb.eq.c0) then
        cf=c1
        return
      endif
      if(x.eq.0.d0) then
        cf=c1
        return
      else if(x.eq.1.d0) then
        cf=cdexp(calgam(cc)-calgam(cc-ca-cb)
     *    -calgam(cc-ca)-calgam(cc-cb))
        return
      else
        call hyper(ca,cb,cc,x,cf)
      endif
      end
**********
c       ====================================================
        SUBROUTINE HYGFX(A,B,C,X,HF)
c       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter 
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI for computing psi function
C       ====================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA(C,G1)
           CALL GAMMA(1.0D0+A/2.0-B,G2)
           CALL GAMMA(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A+M,GAM)
              CALL GAMMA(B+M,GBM)
              CALL PSI(A,PA)
              CALL PSI(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              CALL GAMMA(C-A-B,GCAB)
              CALL GAMMA(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function â(x)
C       Input :  x  --- Argument of â(x)
C                       ( x is not equal to 0,-1,-2,úúú)
C       Output:  GA --- â(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        SUBROUTINE PSI(X,PS)
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
C************************************************************
c
c     Integrador.for
c
c*************************************** fin  programa principal ******
c     INTEGRADOR POR SIMPSO en [ymin, infinito)
c     AJUSTANDO LA PRESICION EN CADA INTERVALO Y LUEGO
c     entre el ultimo intervalo Y LA SUMA TOTAL
c     Para funciones Reales  
c      (10-8-94)                        autor:  Cesar Ramirez
 
C*****
      subroutine INTDEF(PRE,Ymin,Ymax,cfun,cint)
      IMPLICIT REAL*8(A-H,O-Z)
      preci=pre
      x1=Ymin
      x2=Ymax
      csum1=0.d0
      cint=0.d0
39    n=1
      csum=0.d0
      h=x2-x1
      call cfun(x1,cf1)
      call cfun(x2,cf2)
      cf12=cf1+cf2
      cfmc=0.d0
      cfm=0.d0
40    csum1=csum
      d=h/n
      csum=(cf12+2.d0*cfmc)*d/6.d0
      do 50 i=1,n
      xm=h/n*(i-0.5d0)+x1
      call cfun(xm,cfm)        
      csum=csum+4.d0*cfm*d/6.d0             
      cfmc=cfmc+cfm
50    continue
      if(dabs(csum-csum1).gt.preci*dabs(csum1))then
      n=n*2.d0
      go to 40
      else
      x1=x2
      x2=x1+h
      cint=cint+csum
      end if
      if(dabs(csum).gt.preci*dabs(cint))then
      go to 39
      else
      continue
      end if            
      return
      end
c***********************************************************************