********************
*Pour ETACHA4
*  Total cross sections for single ionization of Hydrogenic or Helium 
*  targets by bare ion impact using the CDW-EIS (RHF) approximation. 
*  (see Fainstein \etal, JPB 24 3091 (1991) and references therein)
*
*  Uses: AVINT : Integration Program - 2nd order polinomial fit
*        HYPER : Gauss Hypergeometric function
*        GALAG : Integration Program - Gauss-Laguerre
*        GQUAD : Integration Program - Gauss
*
*  Input Files  : gcdw.dat -->  NG   : no. of points for Gauss-Laguerre 
*                                      integration
*                               ALPHA: scaling factor for Gauss-Laguerre 
*                                      integration
*                               ND   : no. of points for Gauss integration
*                               NA   : no. of points for angular integration
*                               DTE  : electron angle (deg)
*
*
*  Author: PabloF  (pablof@cab.cnea.gov.ar)     Last Version: 14/06/2000
*          Hydrogenic 2s & 2p initial states by Mariel
********************
      subroutine tceis(tot,netat)
      implicit none
      integer i,j,in,na,nele,nang,ninis,ng,nd,nf,npedo,jconv,icor,netat
	real*4 zp1,zt1,E1,tot
      real*8 pi,a0,au,two,ei,zt,zp,epamu,ga2,alpha,
     *       zef,pv,rnu,eene,ek,de,fi,qmin,qmm,
     *       dte,eand,table,te,eanr,ep,ze,fip,
     *       galag,gquad,avint,d100,dd0,dd1,dcoufa,
     *       dd100,ddsia,ddsie,sdene,sdang,tcene,tcang,tcs,
     *		sc,coa,scf,fs1
      real*8 enstar,enstep,enelog
      real*8 frzt,srzt
      complex*16 c1,ci,ca0,ca1,cb0,cb1
      parameter(nele=100,nang=50)
      parameter(pi=3.1415926d0,a0=5.2917706e-9,au=27.2d0,two=2.d0)
      parameter(enstar=-2.0d0,enstep=0.125d0)
      parameter(c1=(1.d0,0.d0),ci=(0.d0,1.d0))
      external frzt,srzt,gauss
      common/pa1/zt,zef,te,ek,fip,sc,coa,zp
      common/pa2/pv,ep,de,fi,rnu,ze,alpha,qmin
      common/pa3/ninis,ng
      common/pa4/ca0,ca1,cb0,cb1
	common/pa5/icor
      dimension eand(nang),eanr(nang),eene(nele),table(nang)
      dimension ddsia(nang),ddsie(nele)
      dimension sdene(nele),sdang(nang),dd100(nele+1,nang+1)
	
	common/don/zp1,E1,zt1

      zt=dble(zp1)
	zp=dble(zt1)
	epamu=1.d03*dble(E1)
	ninis=netat

      jconv=0
	ng=20
	alpha=5.
	nd=15
	na=18
	DATA (table(i),i=1,18) / 0.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,
     s 60.0,70.0,80.0,90.0,100.0,110.0,130.0,150.0,170.0,180.0 /
	do 10 i=1,18
	eand(i)=table(i)
   10	continue

	ga2=(1.d0+epamu/931.5d3)*(1.d0+epamu/931.5d3)
	pv=1.3703599976D2*dsqrt(1.d0-1.d0/ga2)
c     pv=dsqrt(epamu/24.98d0)
      rnu=zp/pv

        if(ninis.eq.10) then
          ei=-zt*zt/2.d0
          zef=zt
          dd0=dexp(-two*pi*rnu)*dcoufa(rnu)*(zp/pv/pi)**2/2.d0
          dd0=zt**5*dd0
          write(*,311)
  311     format(1x,'Hydrogenic Target - 1s Initial State'/)
        else if(ninis.eq.20) then
          ei=-zt*zt/8.d0
          zef=zt
          dd0=dexp(-two*pi*rnu)*dcoufa(rnu)*(zp/pv/pi/4.d0)**2
          dd0=zt**5*dd0
          write(*,312)
  312     format(1x,'Hydrogenic Target - 2s Initial State'/)
        else if(ninis.eq.21) then
          ei=-zt*zt/8.d0
          zef=zt
          dd0=dexp(-two*pi*rnu)*dcoufa(rnu)*(zt*zp/pv/pi)**2/6.d0
          dd0=(zt/2.d0)**5*dd0
          write(*,313)
  313     format(1x,'Hydrogenic Target - 2p Initial State'/)
        endif

      ca0=ci*rnu
      ca1=ci*rnu+c1
	
	scf=dble(1.13*zp**(1./3.))
	fs1=1.d0
	if (zp.eq.1.) fs1=1.23d0
	if (zp.eq.2.) fs1=1.526d0
	if (zp.eq.6.) fs1=0.78d0
	if (zp.eq.7.) fs1=0.85d0
	if (zp.eq.10.) fs1=1.04d0
	if (zp.eq.13.) fs1=0.58d0
	if (zp.eq.14.) fs1=0.59d0
	if (zp.eq.18.) fs1=0.68d0
	if (zp.eq.29.) fs1=0.672d0
	if (zp.eq.36.) fs1=0.61d0
	if (zp.eq.54.) fs1=0.535d0
	sc=scf*fs1
	icor=1
     
*****Comienza el calculo de la DDCS.*******************************************
      i=1
      enelog=enstar
  310 eene(i)=10.d0**enelog

      ek=dsqrt(2.d0*eene(i)/au)
      de=0.5d0*ek*ek-ei
      fi=zef/ek
      qmin=de/pv
	qmm=-ei/pv
	if (pv.gt.1.75d0*qmin) then
	coa=1.d0-(1.75d0*qmin/pv)**2
	else
	coa=0.0d0
	end if
      do j=1,na
        dte=eand(j)
        if(dte.eq.0.d0) then
          te=0.d0
        else if(dte.eq.180.d0) then
          te=pi
        else
          te=dte*pi/180.d0 
        endif
        eanr(j)=te
        ep=dsqrt(ek*ek-2.d0*ek*pv*dcos(te)+pv*pv)
        ze=zp/ep
        cb0=ci*ze
        cb1=ci*ze+c1

        if(dte.eq.0.d0.or.dte.eq.180.d0) then
          fip=0.d0
          d100=two*pi*galag(ng,frzt,alpha,qmin)
        else
          d100=two*gquad(srzt,0.d0,pi,nd,npedo)
        end if

        dd1=dcoufa(fi)*dcoufa(ze)
        dd100(i,j)=ek*d100*dd0*dd1*a0**2/au 
        ddsia(j)=dsin(eanr(j))*dd100(i,j)
      enddo

      nf=i
      dd100(i,na+1)=2.d0*pi*avint(eanr,ddsia,na,eanr(1),eanr(na))
      sdene(i)=dd100(i,na+1)
      write(*,700) eene(i),sdene(i)
      if(sdene(i).lt.1.D-08*sdene(1)) go to 350
      enelog=enelog+enstep
      i=i+1
      if(i.gt.nele) then
        call error(1)
        call error(2)
        stop
      endif
      go to 310
      jconv=1

  350 continue
	do j=1,na
       do i=1,nf
          ddsie(i)=dd100(i,j)
        enddo
        dd100(nf+1,j)=avint(eene,ddsie,nf,eene(1),eene(nf))
        sdang(j)=dsin(eanr(j))*dd100(nf+1,j)
      enddo

*****Comienza el calculo de la TCS.********************************************
      tcene=avint(eene,sdene,nf,eene(1),eene(nf))
      tcang=2.d0*pi*avint(eanr,sdang,na,eanr(1),eanr(na))
      tcs=0.5d0*(tcene+tcang)
c	tcs=tcene
      if(zp.lt.1.d0) then
        tcs=tcs/zp**2
      endif
	write(*,*)' '
	if(icor.eq.1) then
	write(*,*)' SC and ASC corrected CDW-EIS cross sections'
	write(*,*)' '
	else
	write(*,*)' uncorrected CDW-EIS cross sections'
	write(*,*)' '
	end if

        if(ninis.eq.10) then
          write(*,301)
  301     format(1x,'Hydrogenic Target - 1s Initial State')
        else if(ninis.eq.20) then
          write(*,302)
  302     format(1x,'Hydrogenic Target - 2s Initial State')
        else if(ninis.eq.21) then
          write(*,303)
  303     format(1x,'Hydrogenic Target - 2p Initial State')
        endif


      write(*,1000) tcene,tcang,Zt,Zp,epamu,tcs 
      
      if(jconv.ne.0) then
        call error(1)
        call error(3)
      endif

	tot=tcs

 5000 continue 

  160	format(/,1x,' ion 1s=',1PE12.5,1x,'cm**2',/,
     *		 1x,' ion 2s=',1PE12.5,1x,'cm**2',/,
     *		 1x,' ion 2p=',1PE12.5,1x,'cm**2')

	

c  600 stop
  700 format(2x,'EEV=',1PE12.5,1x,'eV',
     *       5x,'SDCS(EEV)=',1PE12.5,1x,'cm**2/eV')
 1000 format(1x/,1x,'tcene(DTE->EEV)=',1PE12.5/,
     *           1x,'tcang(EEV->DTE)=',1PE12.5/,
     *       1x/,1x,'Zp=',0Pf4.1,1x,'Zt=',f4.1,
     *	       1x,'EPAMU=',1PE12.5,1x,'keV/amu',/,5x,
     *       'tcs = 0.5*(tcene+tcang) =',1PE12.5,1x,'cm**2')
	return
      end
******
      function srzt(x)
      implicit none
      real*8 x,srzt,sc,coa,zp
      real*8 frzt
      integer ninis,ng
      real*8 alpha,zt,zef,pv,rnu,ek,de,fi,qmin,te,ep,ze,fip,galag
      external frzt
      common/pa1/zt,zef,te,ek,fip,sc,coa,zp
      common/pa2/pv,ep,de,fi,rnu,ze,alpha,qmin
      common/pa3/ninis,ng

      fip=x
      srzt=galag(ng,frzt,alpha,qmin)
      return
      end
******
      function frzt(x)
      implicit none
      real*8 x,frzt,sc,coa,sq,zp
      integer ninis,nt,ng,icor
      real*8 zt,zef,te,ek,fip,pv,ep,de,fi,rnu,ze,alpha,qmin
      real*8 eta,al,ar,ag,ad
      complex*16 crcei,crx,cry,crz
      common/pa1/zt,zef,te,ek,fip,sc,coa,zp
      common/pa2/pv,ep,de,fi,rnu,ze,alpha,qmin
      common/pa3/ninis,ng
	common/pa5/icor

      eta=dsqrt(x*x-qmin*qmin)
      al=0.5d0*x*x
      ar=-ek*dsin(te)*eta*dcos(fip)-ek*dcos(te)*qmin
      ag=al+ar+de
      ad=-de+pv*(ek*dcos(te)-pv-ep)

c      nt=dint(zt)
c      if(zt/nt.eq.1.d0) then
        if(ninis.eq.10) then
          call hy1star(x,al,ar,ag,ad,crcei)
          frzt=x*(cdabs(crcei))**2
        else if(ninis.eq.20) then
          call hy2star(x,al,ar,ag,ad,crcei)
          frzt=x*(cdabs(crcei))**2
        else if(ninis.eq.21) then
          call hy2ptar(x,eta,al,ar,ag,ad,crx,cry,crz)
          frzt=x*((cdabs(crx))**2+(cdabs(cry))**2+(cdabs(crz))**2)
        endif
c      else
c          call hetar(x,al,ar,ag,ad,crcei)
c          frzt=x*(cdabs(crcei))**2
c      endif
c	SC and ASC corrections
	if (icor.eq.1) then
	sq=(1.d0+(sc/x)**2)**(-2)+coa*(1.d0-(1.d0+(x/sc)**2)**(-2))/zp
c	write (6,*) 'sq=', sq
c	pause
	frzt=frzt*sq
	end if
      return
      end
******
      subroutine hy1star(x,al,ar,ag,ad,crcei)
*
* R(\eta) for Hydrogenic Target 1s Initial State
*
      implicit none
      real*8 x,al,ar,ag,ad,sc,coa,zp
      complex*16 crcei
      integer ninis,ng
      real*8 zt,zef,te,ek,fip,pv,ep,de,fi,rnu,ze,alpha,qmin,zc
      complex*16 c1,c2,ci,ca0,ca1,cb0,cb1,cf1,cf2,cgg,cjj,cbc,ccc,cgc,
     *           crbn1,cauxs
      parameter(c1=(1.d0,0.d0),c2=(2.d0,0.d0),ci=(0.d0,1.d0))
      common/pa1/zt,zef,te,ek,fip,sc,coa,zp
      common/pa2/pv,ep,de,fi,rnu,ze,alpha,qmin
      common/pa3/ninis,ng
      common/pa4/ca0,ca1,cb0,cb1

      zc=1.d0+al*ad/de/ag
      call hypcei(ca0,cb0,c1,zc,cf1)
      call hypcei(ca1,cb1,c2,zc,cf2)

      cgg=(ar+ek*ek*(c1+ci*fi))/ag
      cjj=c1-cgg
      cbc=2.d0*al+ar*(c1+ci*fi)
      ccc=pv/ep*ag*cgg
     *   -(1.d0+pv/ep)*(-de+pv*ek*dcos(te)*(c1+ci*fi))
      cgc=-al*(ad*cbc+ag*ccc)/de/ag/cbc

      crbn1=cdexp(-ci*fi*cdlog(cjj))*cbc/cjj/al/ag**3
      cauxs=cf1-ci*rnu*cgc*cf2
      crcei=crbn1*cauxs
      return
      end
******
      subroutine hy2star(x,al,ar,ag,ad,crcei)
*
* R(\eta) for Hydrogenic Target 2s Initial State
*
*   Author: Mariel                     Last Version: 12/06/2000
*
      implicit none
      real*8 x,al,ar,ag,ad,sc,coa,zp
      complex*16 crcei
      integer ninis,ng
      real*8 zt,zef,te,ek,fip,pv,ep,de,fi,rnu,ze,alpha,qmin,zc
      real*8 bet,fibet,eta,gam
      complex*16 c0,c1,c2,ci,ca0,ca1,cb0,cb1
      complex*16 cf1,cf2,cgg,cjc,cbc,ccc,cgc,cr1,crh1s,cor,caux1,
     *           cfac,caux2,cr2,crh2s
      parameter(c0=(0.d0,0.d0), c1=(1.d0,0.d0), c2=(2.d0,0.d0),
     *          ci=(0.d0,1.d0))
      common/pa1/zt,zef,te,ek,fip,sc,coa,zp
      common/pa2/pv,ep,de,fi,rnu,ze,alpha,qmin
      common/pa3/ninis,ng
      common/pa4/ca0,ca1,cb0,cb1

      zc=1.d0+al*ad/de/ag
      call hypcei(ca0,cb0,c1,zc,cf1)
      call hypcei(ca1,cb1,c2,zc,cf2)

      bet=zt/2.d0
      fibet=bet/ek
      gam=0.5d0*(bet*bet+x*x+2.d0*ar+ek*ek)

      cgg=(ar+ek*ek*(c1+ci*fibet))/gam
      cjc=c1-cgg
      cbc=x*x+ar*(c1+ci*fibet)
      ccc=pv/ep*gam*cgg
     *    -(1.d0+pv/ep)*(-de+pv*ek*dcos(te)*(c1+ci*fibet))
      cgc=-al*(ad*cbc+ag*ccc)/de/ag/cbc
      cr1=cbc*cdexp(-ci*fi*cdlog(cjc))/cjc/gam**2
      crh1s=(cr1*cf1-ci*rnu*cr1*cgc*cf2)/(al*ag)
c
c***********************Estado de Slater 2s************************
c
      cor=-2.d0*bet*cjc+(c1+ci*fi)*(ci*ek-bet*cgg)
      caux1=cbc*cor+ci*cjc*gam*ar/ek
      cfac=-(al/de/ag)*(ad*ar/ek+ag*pv/ep*(ek-pv*dcos(te)-ep*dcos(te)))
      caux2=cbc*cor*cgc+ci*cjc*gam*cfac
      cr2=cdexp(-ci*fi*cdlog(cjc))/cjc**2/gam**3
      crh2s=cr2*(caux1*cf1-ci*rnu*caux2*cf2)/al/ag
c
c******************************************************************
c
      crcei=crh1s+bet*crh2s
      return
      end
******
      subroutine hy2ptar(x,eta,al,ar,ag,ad,crx,cry,crz)
*
* R(\eta) for Hydrogenic Target 2p Initial State
*
*   Author: Mariel                     Last Version: 12/06/2000
*
      implicit none
      real*8 x,al,ar,ag,ad,sc,coa,zp
      complex*16 crx,cry,crz
      integer ninis,ng,j
      real*8 zt,zef,te,ek,fip,pv,ep,de,fi,rnu,ze,alpha,qmin,zc
      real*8 bet,fibet,eta,gam
      real*8 mek(3),mpv(3),meta(3),mx(3),mep(3)
      complex*16 c0,c1,c2,ci,ca0,ca1,cb0,cb1
      complex*16 cf1,cf2,cgg,cjc,cbc,cor2p,cu2p,cm2p,cl2p,
     *           caux1,caux2,crsl1,cr2p(3)
      parameter(c0=(0.d0,0.d0), c1=(1.d0,0.d0), c2=(2.d0,0.d0),
     *          ci=(0.d0,1.d0))
      common/pa1/zt,zef,te,ek,fip,sc,coa,zp
      common/pa2/pv,ep,de,fi,rnu,ze,alpha,qmin
      common/pa3/ninis,ng
      common/pa4/ca0,ca1,cb0,cb1

      zc=1.d0+al*ad/de/ag
      call hypcei(ca0,cb0,c1,zc,cf1)
      call hypcei(ca1,cb1,c2,zc,cf2)

      mek(1)=ek*dsin(te)
      mek(2)=0.d0
      mek(3)=ek*dcos(te)
      mpv(1)=0.d0
      mpv(2)=0.d0
      mpv(3)=pv
      meta(1)=eta*dcos(fip)
      meta(2)=eta*dsin(fip)
      meta(3)=0.d0
      mep(1)=mek(1)-mpv(1)
      mep(2)=mek(2)-mpv(2)
      mep(3)=mek(3)-mpv(3)

      bet=zt/2.d0
      fibet=bet/ek
      gam=0.5d0*(bet*bet+x*x+2.d0*ar+ek*ek)

      cgg=(ar+ek*ek*(c1+ci*fibet))/gam
      cjc=c1-cgg
      cbc=x*x+ar*(c1+ci*fibet)
      cm2p=ar+ek*ek*(c1+ci*fibet)
     *    +(1.d0+ep/pv)*(de-pv*ek*dcos(te)*(c1+ci*fibet))
      crsl1=cdexp(-ci*fi*cdlog(cjc))/(cjc*cjc)/gam**3

      do j=1,3
        mx(j)=-meta(j)-de*mpv(j)/(pv*pv)
        cor2p=(c1-ci*fi)*(mx(j)+mek(j))*cjc+(c1+ci*fi)*mx(j)
        cu2p=-cbc*cor2p+gam*cjc*mx(j)
        cl2p=-cm2p*cor2p+gam*cjc*(mep(j)-ep/pv*mpv(j))
        caux1=cu2p*(cf1+ci*rnu*(al*ad/de/ag)*cf2)
        caux2=ci*rnu*(al*pv/de/ep)*cl2p*cf2
        cr2p(j)=-crsl1*(caux1+caux2)/(al*ag)
      enddo

      crx=cr2p(1)
      cry=cr2p(2)
      crz=cr2p(3)
      return
      end
******
      subroutine hetar(x,al,ar,ag,ad,crcei)
*
* R(\eta) for Hellium Target
*
      implicit none
      real*8 x,al,ar,ag,ad,sc,coa,zp
      complex*16 crcei
      integer ninis,ng,j
      real*8 zt,zef,te,ek,fip,pv,ep,de,fi,rnu,ze,alpha,qmin,zc,
     *       b100(5),d100(5),gam,fij
      complex*16 c0,c1,c2,ci,ca0,ca1,cb0,cb1,cf1,cf2,
     *           cgg,cjc,cbc,ccc,cgc,cs1,cs2,cr1,cr2
      parameter(c0=(0.d0,0.d0),c1=(1.d0,0.d0),
     *          c2=(2.d0,0.d0),ci=(0.d0,1.d0))
      common/pa1/zt,zef,te,ek,fip,sc,coa,zp
      common/pa2/pv,ep,de,fi,rnu,ze,alpha,qmin
      common/pa3/ninis,ng
      common/pa4/ca0,ca1,cb0,cb1
      data b100(1),b100(2),b100(3),b100(4),b100(5)/
     *     1.41714d0,2.37682d0,4.39628d0,6.52699d0,7.94252d0/,
     *     d100(1),d100(2),d100(3),d100(4),d100(5)/
     *     2.5925d0,1.6377d0,0.75254d0,-0.3315d0,0.103d0/
      
      zc=1.d0+al*ad/de/ag
      call hypcei(ca0,cb0,c1,zc,cf1)
      call hypcei(ca1,cb1,c2,zc,cf2)

      cs1=c0
      cs2=c0
      do j=1,5
        gam=0.5d0*(b100(j)**2+2.d0*al+ek*ek+2.d0*ar)
        fij=b100(j)/ek

        cgg=(ar+ek*ek*(c1+ci*fij))/gam
        cjc=c1-cgg
        cbc=2.d0*al+ar*(c1+ci*fij)
        ccc=pv/ep*gam*cgg
     *     -(1.d0+pv/ep)*(-de+pv*ek*dcos(te)*(c1+ci*fij))
        cgc=-al*(ad*cbc+ag*ccc)/de/ag/cbc

        cr1=d100(j)*cbc*cdexp(-ci*fi*cdlog(cjc))/cjc/gam**2
        cr2=cr1*cgc

        cs1=cs1+cr1
        cs2=cs2+cr2
      enddo

      crcei=(cs1*cf1-ci*rnu*cs2*cf2)/(al*ag)
      return
      end
******
      function dcoufa(x)
      implicit none
      real*8 x,dcoufa
      real*8 pi
      parameter(pi=3.1415926d0)
      dcoufa=2.d0*pi*x/(1.d0-dexp(-2.d0*pi*x))
      return
      end
******
      subroutine error(kod)
      implicit none
      integer kod

      if(kod.eq.1) then
        write(*,480)
      else if(kod.eq.2) then
        write(*,481)
      else if(kod.eq.3) then
        write(*,482)
      else if(kod.eq.4) then
        write(*,490)
      else if(kod.eq.5) then
        write(*,495)
      endif

  480 format(1x,/,
     *  1x,'I am sorry but this calculation did not converge',/,
     *  1x,'Please address your complain to pablof@cab.cnea.gov.ar')
  481 format(1x,'Error Code Number 1: vector dimension')
  482 format(1x,'Error Code Number 2: convergence')
  490 format(1x,'This initial state is not yet available')
  495 format(1x,'l > n-1 , please read a book on Quantum Mechanics !')
      return
      end
c***********************************************************************
	DOUBLE PRECISION FUNCTION GQUAD(F,AX,BX,NX,N)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL GAUSS,F
C     N-POINT GAUSSIAN QUADRATURE OF FUNCTION F OVER INTERVAL (AX,BX).
C
      COMMON /GQCOM/A(273),X(273),KTAB(96)
C
C-----TEST N
      N=NX
      ALPHA=0.5D0*(BX+AX)
      BETA=0.5D0*(BX-AX)
      N=MAX0(1,N)
      N=MIN0(96,N)
      IF(N.NE.1) GO TO 1
      GQUAD=(BX-AX)*F(ALPHA)
      RETURN
C
    1 IF(N.GT.16) N=4*(N/4)
      IF(N.GT.24) N=8*(N/8)
      IF(N.GT.48) N=16*(N/16)
C
C----- SET K EQUAL TO INITIAL SUBSCRIPT AND INTEGRATE
      K=KTAB(N)
      M=N/2
      SUM=0.0D0
      JMAX=K-1+M
      DO 2 J=K,JMAX
      DELTA=BETA*X(J)
      SUM=SUM+A(J)*(F(ALPHA+DELTA)+F(ALPHA-DELTA))
    2 CONTINUE
      IF((N-M-M).EQ.0) GO TO 3
      JMID=K+M
      SUM=SUM+A(JMID)*F(ALPHA)
    3 GQUAD=BETA*SUM
      RETURN
      END

      BLOCK DATA GAUSS
C     (DATA BLOCK FOR GQUAD AND GSET)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /GQCOM/A(273),X(273),KTAB(96)
C
C-----TABLE OF INITIAL SUBSCRIPTS FOR N=2(1)16(4)96
      DATA KTAB(2)/1/
      DATA KTAB(3)/2/
      DATA KTAB(4)/4/
      DATA KTAB(5)/6/
      DATA KTAB(6)/9/
      DATA KTAB(7)/12/
      DATA KTAB(8)/16/
      DATA KTAB(9)/20/
      DATA KTAB(10)/25/
      DATA KTAB(11)/30/
      DATA KTAB(12)/36/
      DATA KTAB(13)/42/
      DATA KTAB(14)/49/
      DATA KTAB(15)/56/
      DATA KTAB(16)/64/
      DATA KTAB(20)/72/
      DATA KTAB(24)/82/
      DATA KTAB(28)/82/
      DATA KTAB(32)/94/
      DATA KTAB(36)/94/
      DATA KTAB(40)/110/
      DATA KTAB(44)/110/
      DATA KTAB(48)/130/
      DATA KTAB(52)/130/
      DATA KTAB(56)/130/
      DATA KTAB(60)/130/
      DATA KTAB(64)/154/
      DATA KTAB(68)/154/
      DATA KTAB(72)/154/
      DATA KTAB(76)/154/
      DATA KTAB(80)/186/
      DATA KTAB(84)/186/
      DATA KTAB(88)/186/
      DATA KTAB(92)/186/
      DATA KTAB(96)/226/
C
C-----TABLE OF ABSCISSAE (X) AND WEIGHTS (A) FOR INTERVAL (-1,+1).
C
C-----N=2
      DATA X(1)/0.577350269189626  /, A(1)/1.000000000000000  /
C-----N=3
      DATA X(2)/0.774596669241483  /, A(2)/0.555555555555556  /
      DATA X(3)/0.000000000000000  /, A(3)/0.888888888888889  /
C-----N=4
      DATA X(4)/0.861136311594053  /, A(4)/0.347854845137454  /
      DATA X(5)/0.339981043584856  /, A(5)/0.652145154862546  /
C-----N=5
      DATA X(6)/0.906179845938664  /, A(6)/0.236926885056189  /
      DATA X(7)/0.538469310105683  /, A(7)/0.478628670499366  /
      DATA X(8)/0.000000000000000  /, A(8)/0.568888888888889  /
C-----N=6
      DATA X(9)/0.932469514203152  /, A(9)/0.171324492379170  /
      DATA X(10)/0.661209386466265 /, A(10)/0.360761573048139 /
      DATA X(11)/0.238619186083197 /, A(11)/0.467913934572691 /
C-----N=7
      DATA X(12)/0.949107912342759 /, A(12)/0.129484966168870 /
      DATA X(13)/0.741531185599394 /, A(13)/0.279705391489277 /
      DATA X(14)/0.405845151377397 /, A(14)/0.381830050505119 /
      DATA X(15)/0.000000000000000 /, A(15)/0.417959183673469 /
C-----N=8
      DATA X(16)/0.960289856497536 /, A(16)/0.101228536290376 /
      DATA X(17)/0.796666477413627 /, A(17)/0.222381034453374 /
      DATA X(18)/0.525532409916329 /, A(18)/0.313706645877887 /
      DATA X(19)/0.183434642495650 /, A(19)/0.362683783378362 /
C-----N=9
      DATA X(20)/0.968160239507626 /, A(20)/0.081274388361574 /
      DATA X(21)/0.836031107326636 /, A(21)/0.180648160694857 /
      DATA X(22)/0.613371432700590 /, A(22)/0.260610696402935 /
      DATA X(23)/0.324253423403809 /, A(23)/0.312347077040003 /
      DATA X(24)/0.000000000000000 /, A(24)/0.330239355001260 /
C-----N=10
      DATA X(25)/0.973906528517172 /, A(25)/0.066671344308688 /
      DATA X(26)/0.865063366688985 /, A(26)/0.149451349150581 /
      DATA X(27)/0.679409568299024 /, A(27)/0.219086362515982 /
      DATA X(28)/0.433395394129247 /, A(28)/0.269266719309996 /
      DATA X(29)/0.148874338981631 /, A(29)/0.295524224714753 /
C-----N=11
      DATA X(30)/0.978228658146057 /, A(30)/0.055668567116174 /
      DATA X(31)/0.887062599768095 /, A(31)/0.125580369464905 /
      DATA X(32)/0.730152005574049 /, A(32)/0.186290210927734 /
      DATA X(33)/0.519096129206812 /, A(33)/0.233193764591990 /
      DATA X(34)/0.269543155952345 /, A(34)/0.262804544510247 /
      DATA X(35)/0.000000000000000 /, A(35)/0.272925086777901 /
C-----N=12
      DATA X(36)/0.981560634246719 /, A(36)/0.047175336386512 /
      DATA X(37)/0.904117256370475 /, A(37)/0.106939325995318 /
      DATA X(38)/0.769902674194305 /, A(38)/0.160078328543346 /
      DATA X(39)/0.587317954286617 /, A(39)/0.203167426723066 /
      DATA X(40)/0.367831498998180 /, A(40)/0.233492536538355 /
      DATA X(41)/0.125233408511469 /, A(41)/0.249147045813403 /
C-----N=13
      DATA X(42)/0.984183054718588 /, A(42)/0.040484004765316 /
      DATA X(43)/0.917598399222978 /, A(43)/0.092121499837728 /
      DATA X(44)/0.801578090733310 /, A(44)/0.138873510219787 /
      DATA X(45)/0.642349339440340 /, A(45)/0.178145980761946 /
      DATA X(46)/0.448492751036447 /, A(46)/0.207816047536889 /
      DATA X(47)/0.230458315955135 /, A(47)/0.226283180262897 /
      DATA X(48)/0.000000000000000 /, A(48)/0.232551553230874 /
C-----N=14
      DATA X(49)/0.986283808696812 /, A(49)/0.035119460331752 /
      DATA X(50)/0.928434883663574 /, A(50)/0.080158087159760 /
      DATA X(51)/0.827201315069765 /, A(51)/0.121518570687903 /
      DATA X(52)/0.687292904811685 /, A(52)/0.157203167158194 /
      DATA X(53)/0.515248636358154 /, A(53)/0.185538397477938 /
      DATA X(54)/0.319112368927890 /, A(54)/0.205198463721296 /
      DATA X(55)/0.108054948707344 /, A(55)/0.215263853463158 /
C-----N=15
      DATA X(56)/0.987992518020485 /, A(56)/0.030753241996117 /
      DATA X(57)/0.937273392400706 /, A(57)/0.070366047488108 /
      DATA X(58)/0.848206583410427 /, A(58)/0.107159220467172 /
      DATA X(59)/0.724417731360170 /, A(59)/0.139570677926154 /
      DATA X(60)/0.570972172608539 /, A(60)/0.166269205816994 /
      DATA X(61)/0.394151347077563 /, A(61)/0.186161000015562 /
      DATA X(62)/0.201194093997435 /, A(62)/0.198431485327111 /
      DATA X(63)/0.000000000000000 /, A(63)/0.202578241925561 /
C-----N=16
      DATA X(64)/0.989400934991650 /, A(64)/0.027152459411754 /
      DATA X(65)/0.944575023073233 /, A(65)/0.062253523938648 /
      DATA X(66)/0.865631202387832 /, A(66)/0.095158511682493 /
      DATA X(67)/0.755404408355003 /, A(67)/0.124628971255534 /
      DATA X(68)/0.617876244402644 /, A(68)/0.149595988816577 /
      DATA X(69)/0.458016777657227 /, A(69)/0.169156519395003 /
      DATA X(70)/0.281603550779259 /, A(70)/0.182603415044924 /
      DATA X(71)/0.095012509837637 /, A(71)/0.189450610455069 /
C-----N=20
      DATA X(72)/0.993128599185094 /, A(72)/0.017614007139152 /
      DATA X(73)/0.963971927277913 /, A(73)/0.040601429800386 /
      DATA X(74)/0.912234428251325 /, A(74)/0.062672048334109 /
      DATA X(75)/0.839116971822218 /, A(75)/0.083276741576704 /
      DATA X(76)/0.746331906460150 /, A(76)/0.101930119817240 /
      DATA X(77)/0.636053680726515 /, A(77)/0.118194531961518 /
      DATA X(78)/0.510867001950827 /, A(78)/0.131688638449176 /
      DATA X(79)/0.373706088715419 /, A(79)/0.142096109318382 /
      DATA X(80)/0.227785851141645 /, A(80)/0.149172986472603 /
      DATA X(81)/0.076526521133497 /, A(81)/0.152753387130725 /
C-----N=24
      DATA X(82)/0.995187219997021 /, A(82)/0.012341229799987 /
      DATA X(83)/0.974728555971309 /, A(83)/0.028531388628933 /
      DATA X(84)/0.938274552002732 /, A(84)/0.044277438817419 /
      DATA X(85)/0.886415527004401 /, A(85)/0.059298584915436 /
      DATA X(86)/0.820001985973902 /, A(86)/0.073346481411080 /
      DATA X(87)/0.740124191578554 /, A(87)/0.086190161531953 /
      DATA X(88)/0.648093651936975 /, A(88)/0.097618652104113 /
      DATA X(89)/0.545421471388839 /, A(89)/0.107444270115965 /
      DATA X(90)/0.433793507626045 /, A(90)/0.115505668053725 /
      DATA X(91)/0.315042679696163 /, A(91)/0.121670472927803 /
      DATA X(92)/0.191118867473616 /, A(92)/0.125837456346828 /
      DATA X(93)/0.064056892862605 /, A(93)/0.127938195346752 /
C-----N=32
      DATA X(94)/0.997263861849481 /, A(94)/0.007018610009470 /
      DATA X(95)/0.985611511545268 /, A(95)/0.016274394730905 /
      DATA X(96)/0.964762255587506 /, A(96)/0.025392065309262 /
      DATA X(97)/0.934906075937739 /, A(97)/0.034273862913021 /
      DATA X(98)/0.896321155766052 /, A(98)/0.042835898022226 /
      DATA X(99)/0.849367613732569 /, A(99)/0.050998059262376 /
      DATA X(100)/0.794483795967942/, A(100)/0.058684093478535/
      DATA X(101)/0.732182118740289/, A(101)/0.065822222776361/
      DATA X(102)/0.663044266930215/, A(102)/0.072345794108848/
      DATA X(103)/0.587715757240762/, A(103)/0.078193895787070/
      DATA X(104)/0.506899908932229/, A(104)/0.083311924226946/
      DATA X(105)/0.421351276130635/, A(105)/0.087652093004403/
      DATA X(106)/0.331868602282127/, A(106)/0.091173878695763/
      DATA X(107)/0.239287362252137/, A(107)/0.093844399080804/
      DATA X(108)/0.144471961582796/, A(108)/0.095638720079274/
      DATA X(109)/0.048307665687738/, A(109)/0.096540088514727/
C-----N=40
      DATA X(110)/0.998237709710559/, A(110)/0.004521277098533/
      DATA X(111)/0.990726238699457/, A(111)/0.010498284531152/
      DATA X(112)/0.977259949983774/, A(112)/0.016421058381907/
      DATA X(113)/0.957916819213791/, A(113)/0.022245849194166/
      DATA X(114)/0.932812808278676/, A(114)/0.027937006980023/
      DATA X(115)/0.902098806968874/, A(115)/0.033460195282547/
      DATA X(116)/0.865959503212259/, A(116)/0.038782167974472/
      DATA X(117)/0.824612230833311/, A(117)/0.043870908185673/
      DATA X(118)/0.778305651426519/, A(118)/0.048695807635072/
      DATA X(119)/0.727318255189927/, A(119)/0.053227846983936/
      DATA X(120)/0.671956684614179/, A(120)/0.057439769099391/
      DATA X(121)/0.612553889667980/, A(121)/0.061306242492928/
      DATA X(122)/0.549467125095128/, A(122)/0.064804013456601/
      DATA X(123)/0.483075801686178/, A(123)/0.067912045815233/
      DATA X(124)/0.413779204371605/, A(124)/0.070611647391286/
      DATA X(125)/0.341994090825758/, A(125)/0.072886582395804/
      DATA X(126)/0.268152185007253/, A(126)/0.074723169057968/
      DATA X(127)/0.192697580701371/, A(127)/0.076110361900626/
      DATA X(128)/0.116084070675255/, A(128)/0.077039818164247/
      DATA X(129)/0.038772417506050/, A(129)/0.077505947978424/
C-----N=48
      DATA X(130)/0.998771007252426/, A(130)/0.003153346052305/
      DATA X(131)/0.993530172266350/, A(131)/0.007327553901276/
      DATA X(132)/0.984124583722826/, A(132)/0.011477234579234/
      DATA X(133)/0.970591592546247/, A(133)/0.015579315722943/
      DATA X(134)/0.952987703160430/, A(134)/0.019616160457355/
      DATA X(135)/0.931386690706554/, A(135)/0.023570760839324/
      DATA X(136)/0.905879136715569/, A(136)/0.027426509708356/
      DATA X(137)/0.876572020274247/, A(137)/0.031167227832798/
      DATA X(138)/0.843588261624393/, A(138)/0.034777222564770/
      DATA X(139)/0.807066204029442/, A(139)/0.038241351065830/
      DATA X(140)/0.767159032515740/, A(140)/0.041545082943464/
      DATA X(141)/0.724034130923814/, A(141)/0.044674560856694/
      DATA X(142)/0.677872379632663/, A(142)/0.047616658492490/
      DATA X(143)/0.628867396776513/, A(143)/0.050359035553854/
      DATA X(144)/0.577224726083972/, A(144)/0.052890189485193/
      DATA X(145)/0.523160974722233/, A(145)/0.055199503699984/
      DATA X(146)/0.466902904750958/, A(146)/0.057277292100403/
      DATA X(147)/0.408686481990716/, A(147)/0.059114839698395/
      DATA X(148)/0.348755886292160/, A(148)/0.060704439165893/
      DATA X(149)/0.287362487355455/, A(149)/0.062039423159892/
      DATA X(150)/0.224763790394689/, A(150)/0.063114192286254/
      DATA X(151)/0.161222356068891/, A(151)/0.063924238584648/
      DATA X(152)/0.097004699209462/, A(152)/0.064466164435950/
      DATA X(153)/0.032380170962869/, A(153)/0.064737696812683/
C-----N=64
      DATA X(154)/0.999305041735772/, A(154)/0.001783280721696/
      DATA X(155)/0.996340116771955/, A(155)/0.004147033260562/
      DATA X(156)/0.991013371476744/, A(156)/0.006504457968978/
      DATA X(157)/0.983336253884625/, A(157)/0.008846759826363/
      DATA X(158)/0.973326827789910/, A(158)/0.011168139460131/
      DATA X(159)/0.961008799652053/, A(159)/0.013463047896718/
      DATA X(160)/0.946411374858402/, A(160)/0.015726030476024/
      DATA X(161)/0.929569172131939/, A(161)/0.017951715775697/
      DATA X(162)/0.910522137078502/, A(162)/0.020134823153530/
      DATA X(163)/0.889315445995114/, A(163)/0.022270173808383/
      DATA X(164)/0.865999398154092/, A(164)/0.024352702568710/
      DATA X(165)/0.840629296252580/, A(165)/0.026377469715054/
      DATA X(166)/0.813265315122797/, A(166)/0.028339672614259/
      DATA X(167)/0.783972358943341/, A(167)/0.030234657072402/
      DATA X(168)/0.752819907260531/, A(168)/0.032057928354851/
      DATA X(169)/0.719881850171610/, A(169)/0.033805161837141/
      DATA X(170)/0.685236313054233/, A(170)/0.035472213256882/
      DATA X(171)/0.648965471254657/, A(171)/0.037055128540240/
      DATA X(172)/0.611155355172393/, A(172)/0.038550153178615/
      DATA X(173)/0.571895646202634/, A(173)/0.039953741132720/
      DATA X(174)/0.531279464019894/, A(174)/0.041262563242623/
      DATA X(175)/0.489403145707052/, A(175)/0.042473515123653/
      DATA X(176)/0.446366017253464/, A(176)/0.043583724529323/
      DATA X(177)/0.402270157963991/, A(177)/0.044590558163756/
      DATA X(178)/0.357220158337668/, A(178)/0.045491627927418/
      DATA X(179)/0.311322871990210/, A(179)/0.046284796581314/
      DATA X(180)/0.264687162208767/, A(180)/0.046968182816210/
      DATA X(181)/0.217423643740007/, A(181)/0.047540165714830/
      DATA X(182)/0.169644420423992/, A(182)/0.047999388596458/
      DATA X(183)/0.121462819296120/, A(183)/0.048344762234802/
      DATA X(184)/0.072993121787799/, A(184)/0.048575467441503/
      DATA X(185)/0.024350292663424/, A(185)/0.048690957009139/
C-----N=80
      DATA X(186)/0.999553822651630/, A(186)/0.001144950003186/
      DATA X(187)/0.997649864398237/, A(187)/0.002663533589512/
      DATA X(188)/0.994227540965688/, A(188)/0.004180313124694/
      DATA X(189)/0.989291302499755/, A(189)/0.005690922451403/
      DATA X(190)/0.982848572738629/, A(190)/0.007192904768117/
      DATA X(191)/0.974909140585727/, A(191)/0.008683945269260/
      DATA X(192)/0.965485089043799/, A(192)/0.010161766041103/
      DATA X(193)/0.954590766343634/, A(193)/0.011624114120797/
      DATA X(194)/0.942242761309872/, A(194)/0.013068761592401/
      DATA X(195)/0.928459877172445/, A(195)/0.014493508040509/
      DATA X(196)/0.913263102571757/, A(196)/0.015896183583725/
      DATA X(197)/0.896675579438770/, A(197)/0.017274652056269/
      DATA X(198)/0.878722567678213/, A(198)/0.018626814208299/
      DATA X(199)/0.859431406663111/, A(199)/0.019950610878141/
      DATA X(200)/0.838831473580255/, A(200)/0.021244026115782/
      DATA X(201)/0.816954138681463/, A(201)/0.022505090246332/
      DATA X(202)/0.793832717504605/, A(202)/0.023731882865930/
      DATA X(203)/0.769502420135041/, A(203)/0.024922535764115/
      DATA X(204)/0.744000297583597/, A(204)/0.026075235767565/
      DATA X(205)/0.717365185362099/, A(205)/0.027188227500486/
      DATA X(206)/0.689637644342027/, A(206)/0.028259816057276/
      DATA X(207)/0.660859898986119/, A(207)/0.029288369583267/
      DATA X(208)/0.631075773046871/, A(208)/0.030272321759557/
      DATA X(209)/0.600330622829751/, A(209)/0.031210174188114/
      DATA X(210)/0.568671268122709/, A(210)/0.032100498673487/
      DATA X(211)/0.536145920897131/, A(211)/0.032941939397645/
      DATA X(212)/0.502804111888784/, A(212)/0.033733214984611/
      DATA X(213)/0.468696615170544/, A(213)/0.034473120451753/
      DATA X(214)/0.433875370831756/, A(214)/0.035160529044747/
      DATA X(215)/0.398393405881969/, A(215)/0.035794393953416/
      DATA X(216)/0.362304753499487/, A(216)/0.036373749905835/
      DATA X(217)/0.325664370747701/, A(217)/0.036897714638276/
      DATA X(218)/0.288528054884511/, A(218)/0.037365490238730/
      DATA X(219)/0.250952358392272/, A(219)/0.037776364362001/
      DATA X(220)/0.212994502857666/, A(220)/0.038129711314477/
      DATA X(221)/0.174712291832646/, A(221)/0.038424993006959/
      DATA X(222)/0.136164022809143/, A(222)/0.038661759774076/
      DATA X(223)/0.097408398441584/, A(223)/0.038839651059051/
      DATA X(224)/0.058504437152420/, A(224)/0.038958395962769/
      DATA X(225)/0.019511383256793/, A(225)/0.039017813656306/
C-----N=96
      DATA X(226)/0.999689503883230/, A(226)/0.000796792065552/
      DATA X(227)/0.998364375863181/, A(227)/0.001853960788946/
      DATA X(228)/0.995981842987209/, A(228)/0.002910731817934/
      DATA X(229)/0.992543900323762/, A(229)/0.003964554338444/
      DATA X(230)/0.988054126329623/, A(230)/0.005014202742927/
      DATA X(231)/0.982517263563014/, A(231)/0.006058545504235/
      DATA X(232)/0.975939174585136/, A(232)/0.007096470791153/
      DATA X(233)/0.968326828463264/, A(233)/0.008126876925698/
      DATA X(234)/0.959688291448742/, A(234)/0.009148671230783/
      DATA X(235)/0.950032717784437/, A(235)/0.010160770535008/
      DATA X(236)/0.939370339752755/, A(236)/0.011162102099838/
      DATA X(237)/0.927712456722308/, A(237)/0.012151604671088/
      DATA X(238)/0.915071423120898/, A(238)/0.013128229566961/
      DATA X(239)/0.901460635315852/, A(239)/0.014090941772314/
      DATA X(240)/0.886894517402420/, A(240)/0.015038721026994/
      DATA X(241)/0.871388505909296/, A(241)/0.015970562902562/
      DATA X(242)/0.854959033434601/, A(242)/0.016885479864245/
      DATA X(243)/0.837623511228187/, A(243)/0.017782502316045/
      DATA X(244)/0.819400310737931/, A(244)/0.018660679627411/
      DATA X(245)/0.800308744139140/, A(245)/0.019519081140145/
      DATA X(246)/0.780369043867433/, A(246)/0.020356797154333/
      DATA X(247)/0.759602341176647/, A(247)/0.021172939892191/
      DATA X(248)/0.738030643744400/, A(248)/0.021966644438744/
      DATA X(249)/0.715676812348967/, A(249)/0.022737069658329/
      DATA X(250)/0.692564536642171/, A(250)/0.023483399085926/
      DATA X(251)/0.668718310043916/, A(251)/0.024204841792364/
      DATA X(252)/0.644163403784967/, A(252)/0.024900633222483/
      DATA X(253)/0.618925840125468/, A(253)/0.025570036005349/
      DATA X(254)/0.593032364777572/, A(254)/0.026212340735672/
      DATA X(255)/0.566510418561397/, A(255)/0.026826866725591/
      DATA X(256)/0.539388108324357/, A(256)/0.027412962726029/
      DATA X(257)/0.511694177154667/, A(257)/0.027970007616848/
      DATA X(258)/0.483457973920596/, A(258)/0.028497411065085/
      DATA X(259)/0.454709422167743/, A(259)/0.028994614150555/
      DATA X(260)/0.425478988407300/, A(260)/0.029461089958167/
      DATA X(261)/0.395797649828908/, A(261)/0.029896344136328/
      DATA X(262)/0.365696861472313/, A(262)/0.030299915420827/
      DATA X(263)/0.335208522892625/, A(263)/0.030671376123669/
      DATA X(264)/0.304364944354496/, A(264)/0.031010332586313/
      DATA X(265)/0.273198812591049/, A(265)/0.031316425596861/
      DATA X(266)/0.241743156163840/, A(266)/0.031589330770727/
      DATA X(267)/0.210031310460567/, A(267)/0.031828758894411/
      DATA X(268)/0.178096882367618/, A(268)/0.032034456231992/
      DATA X(269)/0.145973714654896/, A(269)/0.032206204794030/
      DATA X(270)/0.113695850110665/, A(270)/0.032343822568575/
      DATA X(271)/0.081297495464425/, A(271)/0.032447163714064/
      DATA X(272)/0.048812985136049/, A(272)/0.032516118713868/
      DATA X(273)/0.016276744849602/, A(273)/0.032550614492363/
C
      END
c************************************************************************
      subroutine hypcei(ca,cb,cc,x,cf)
**************************************
*  Program: hyper.f
*  Version: 1.6 (see readme_hyper) 24/03/1997
*  Author: Pablof (pablof@cab.cnea.gov.ar)
*  Hypergeometric function for CDW-EIS calculations. 
*  Double precision, see comment for V 1.4
*  We assume that x is REAL
*  ^^^^^^^^^^^^^^^^^^^^^^^^
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
      subroutine hyper(ca,cb,cc,x,cf)
      implicit none
      complex*16 ca,cb,cc,cf
      real*8 x
      complex*16 c1
      real*8 z,step1,step2,step3,step4
      parameter(step1=-200d0,step2=-0.618d0,step3=0.5d0,
     *          step4=1.1d0)
      parameter(c1=(1.d0,0.d0))

      if(x.lt.step1) then
	if(ca.eq.cb) then
          z=x/(x-1.d0)
	  call ca15310(cc-ca,cb,z,cf)
          cf=cdexp(-cb*dlog(1.d0-x))*cf
	  return
	else
          call ca1538(ca,cb,cc,x,cf)
	  return
	endif
      else if(x.lt.step2) then
	if(ca.eq.cb) then
          z=x/(x-1.d0)
          call chypser(cc-ca,cb,cc,z,cf)
          cf=cdexp(-cb*dlog(1.d0-x))*cf
	  return
	else
          call ca1538(ca,cb,cc,x,cf)
	  return
	endif
      else if(x.lt.step3) then
        call chypser(ca,cb,cc,x,cf)
        return
      else if(x.lt.step4) then
        if(cc.eq.ca+cb) then
	  call ca15310(ca,cb,x,cf)
	  return
	else
	  call ca1536(ca,cb,cc,x,cf)
          return
        endif
      else
	if(ca.eq.cb) then
          call ca15313(ca,cc,x,cf)
	  return
	else
          call ca1537(ca,cb,cc,x,cf)
	  return
	endif
      endif
      end
**********
      subroutine ca1536(ca,cb,cc,x,cf)
**************************************
*  Uses formula 15.3.6 from M. Abramowitz and I. Stegun
**************************************
      implicit none
      complex*16 ca,cb,cc,cf
      real*8 x,z,pi
      complex*16 c0,c1,cz,cca,ccb,ccab,calgam,cgam1,cf1,cgam2,cf2
      parameter(pi=3.141592653589793238462643d0)
      parameter(c0=(0.d0,0.d0), c1=(1.d0,0.d0))

      cca=cc-ca
      ccb=cc-cb
      ccab=cc-ca-cb
      if(ccab.eq.c0) pause 'There is a pole cc=ca+cb'
      
      z=1.d0-x
      cz=dcmplx(z,0.d0)
      cgam1=cdexp(calgam(cc)+calgam(ccab)-calgam(cca)-calgam(ccb))
      cgam2= cdexp(calgam(cc)+calgam(-ccab)-calgam(ca)-calgam(cb))*
     *       cdexp(ccab*cdlog(cz))
      call chypser(ca,cb,c1-ccab,z,cf1)
      call chypser(cca,ccb,c1+ccab,z,cf2)
      cf=cgam1*cf1+cgam2*cf2
      return
      end
**********
      subroutine ca1537(ca,cb,cc,x,cf)
**************************************
*  Uses formula 15.3.7 from M. Abramowitz and I. Stegun
**************************************
      implicit none
      complex*16 ca,cb,cc,cf
      real*8 x,z,pi
      complex*16 c0,c1,cca,ccb,cba,cx,calgam,cgam1,cf1,cgam2,cf2
      parameter(pi=3.141592653589793238462643d0)
      parameter(c0=(0.d0,0.d0), c1=(1.d0,0.d0))
      
      cca=cc-ca
      ccb=cc-cb
      cba=cb-ca
      if(cba.eq.c0) pause 'There is a pole ca=cb'

      cx=dcmplx(x,0.d0)
      cgam1=cdexp(calgam(cc)+calgam(cba)-calgam(cb)-calgam(cca))*
     *      cdexp(-ca*cdlog(-cx))
      cgam2=cdexp(calgam(cc)+calgam(-cba)-calgam(ca)-calgam(ccb))*
     *      cdexp(-cb*cdlog(-cx))
      z=1.d0/x
      call chypser(ca,c1-cca,c1-cba,z,cf1)
      call chypser(cb,c1-ccb,c1+cba,z,cf2)
      cf=cgam1*cf1+cgam2*cf2
      return
      end
**********
      subroutine ca1538(ca,cb,cc,x,cf)
**************************************
*  Uses formula 15.3.8 from M. Abramowitz and I. Stegun
**************************************
      implicit none
      complex*16 ca,cb,cc,cf
      real*8 x,z
      complex*16 c0,c1,cca,ccb,cba,calgam,cgam1,cf1,cgam2,cf2
      parameter(c0=(0.d0,0.d0), c1=(1.d0,0.d0))

      cca=cc-ca
      ccb=cc-cb
      cba=cb-ca
      if(cba.eq.c0) pause 'There is a pole ca=cb'

      cgam1=cdexp(calgam(cc)+calgam(cba)-calgam(cb)-calgam(cca))*
     *      cdexp(-ca*dlog(1.d0-x))
      cgam2=cdexp(calgam(cc)+calgam(-cba)-calgam(ca)-calgam(ccb))*
     *      cdexp(-cb*dlog(1.d0-x))
      z=1.d0/(1.d0-x)
      call chypser(ca,ccb,c1-cba,z,cf1)
      call chypser(cb,cca,c1+cba,z,cf2)
      cf=cgam1*cf1+cgam2*cf2
      return
      end
**********
      subroutine ca15310(ca,cb,x,cf)
**************************************
*  Uses formula 15.3.10 from M. Abramowitz and I. Stegun
**************************************
      implicit none
      complex*16 ca,cb,cf
      real*8 x
      complex*16 c1,caa,cbb,cc,psi1,psi2,psi3,
     *           cmx,coef,csum,cfac,ctemp,calgam
      real*8 pi
      integer n
      parameter(pi=3.141592653589793238462643d0)
      parameter(c1=(1.d0,0.d0))

      n=1
      caa=ca
      cbb=cb
      cc=caa+cbb
      call digam(c1,psi1)
      call digam(caa,psi2)
      call digam(cbb,psi3)

      cmx=dcmplx(1.d0-x,0.d0)
      coef=-cdlog(cmx)
      csum=coef+2.d0*psi1-psi2-psi3
      cfac=c1
      ctemp=cfac*csum

      psi1=psi1+1.d0
      psi2=psi2+1.d0/caa
      psi3=psi3+1.d0/cbb

100   continue
        cfac=cfac*caa*cbb*cmx/(n*n)
        csum=coef+2.d0*psi1-psi2-psi3
        cf=ctemp+cfac*csum
        if(cf.eq.ctemp) goto 200
        ctemp=cf
        n=n+1
        caa=caa+c1
        cbb=cbb+c1
        psi1=psi1+1.d0/n
        psi2=psi2+1.d0/caa
        psi3=psi3+1.d0/cbb
      go to 100
200   continue
      cf=cf*cdexp(calgam(cc)-calgam(ca)-calgam(cb))
      return
      end
**********
      subroutine ca15313(ca,cc,x,cf)
**************************************
*  Uses formula 15.3.13 from M. Abramowitz and I. Stegun
**************************************
      implicit none
      complex*16 ca,cc,cf
      real*8 x
      complex*16 c1,caa,cbb,psi1,psi2,psi3,
     *           cmx,coef,csum,cfac,ctemp,calgam
      real*8 pi
      integer n
      parameter(pi=3.141592653589793238462643d0)
      parameter(c1=(1.d0,0.d0))

      n=1
      caa=ca
      cbb=c1+ca-cc
      call digam(c1,psi1)
      call digam(caa,psi2)
      call digam(cbb,psi3)

      cmx=dcmplx(-x,0.d0)
      coef=cdlog(cmx)-pi*cdcos(pi*cbb)/cdsin(pi*cbb)
      csum=coef+2.d0*psi1-psi2-psi3
      cfac=c1
      ctemp=cfac*csum

      psi1=psi1+1.d0
      psi2=psi2+1.d0/caa
      psi3=psi3+1.d0/cbb

100   continue
        cfac=cfac*caa*cbb/(n*n*x)
        csum=coef+2.d0*psi1-psi2-psi3
        cf=ctemp+cfac*csum
        if(cf.eq.ctemp) goto 200
        ctemp=cf
        n=n+1
        caa=caa+c1
        cbb=cbb+c1
        psi1=psi1+1.d0/n
        psi2=psi2+1.d0/caa
        psi3=psi3+1.d0/cbb
      go to 100
200   continue
      cf=cf*cdexp(calgam(cc)-calgam(ca)-calgam(cc-ca))*
     *      cdexp(-ca*cdlog(cmx))
      return
      end
**********
      subroutine chypser(ca,cb,cc,x,cf)
**************************************
*  Hypergeometric series calculated to machine accuracy
**************************************
      implicit none
      complex*16 ca,cb,cc,cf
      real*8 x
      complex*16 c0,c1,cfac,ctemp,caa,cbb,ccc
      integer n
      parameter(c0=(0.d0,0.d0), c1=(1.d0,0.d0))
      
      if(ca.eq.c0.or.cb.eq.c0) then
        cf=c1
        return
      endif
      n=1
      caa=ca
      cbb=cb
      ccc=cc
      cfac=c1
      ctemp=cfac
100   continue
        cfac=((caa*cbb)/ccc)*cfac
        cfac=cfac*x/n
        cf=ctemp+cfac
        if(cf.eq.ctemp) go to 200
        ctemp=cf
        n=n+1
        caa=caa+c1
        cbb=cbb+c1
        ccc=ccc+c1
      go to 100
200   continue
      return
      end
**********
      subroutine digam(z,psi)
**************************************
*  Subroutine to calculate the DIGAMMA function of a complex
*  argument using formulas from M. Abramowitz and I. Stegun.
*
*  If z = integer --> uses A6.3.2  --> gives very good results
*  If z = complex --> uses A6.3.16 --> not so good ...
*
*  If Real(z).lt.0.d0 it uses the reflection formula A6.3.7
*  to avoid the poles in the series expansion.
*
*  Author: Pablof    V 0.8  15/12/1994
**************************************
      implicit none
      complex*16 z,psi
      complex*16 c1,diser
      real*8 xz,gama,pi
      parameter(gama=0.577215664901532860606512d0,
     *          pi=3.141592653589793238462643d0)
      parameter(c1=(1.d0,0.d0))
      if(z.eq.c1) then
        psi=-c1*gama
        return
      endif
      xz=dble(z)
      if(xz.lt.0.d0) then
        psi=diser(c1-z)-pi*cdcos(pi*z)/cdsin(pi*z)
        return
      else
        psi=diser(z)
        return
      endif
      end
**********
      complex*16 function diser(z)
      implicit none
      complex*16 z
      real*8 xz,yz,divx,gama
      integer nz,k
      parameter(gama=0.577215664901532860606512d0)
      nz=int(z)
      xz=dble(z)
      yz=dimag(z)
      if(xz.ne.0.d0) divx=nz/xz
      if(xz.eq.0.d0.or.divx.eq.1.d0) then
	if(yz.eq.0.d0) then
	  diser=-gama
	  do k=1,nz-1
            diser=diser+1.d0/k
          enddo
          return
        else
          call digi1(nz,xz,yz,diser)
	  return
        endif
      else
        call digi2(z,diser)
        return
      endif
      end
**********
      subroutine digi1(nz,xz,yz,diser)
      implicit none
      complex*16 diser
      real*8 xz,yz
      integer nz
      complex*16 z
      real*8 gama,fac,fac1,fac2,temp1,temp2,rpsi,cpsi
      integer n
      parameter(gama=0.577215664901532860606512d0)
      fac=1.d0/(1.d0+yz*yz)
      temp1=fac
      temp2=fac
      do n=2,10000
	fac2=1.d0/(n*n+yz*yz)
	fac1=fac2/n
	rpsi=temp1+fac1
	cpsi=temp2+fac2
	if(rpsi.eq.temp1.and.cpsi.eq.temp2) goto 200
	temp1=rpsi
	temp2=cpsi
      enddo
c     pause 'convergence failure in digi1'
200   continue
      rpsi=-gama+yz*yz*rpsi
      cpsi=yz*cpsi
      if(xz.eq.0.d0) then
        diser=dcmplx(rpsi,cpsi+1.d0/yz)
        return
      else if(xz.eq.1.d0) then
        diser=dcmplx(rpsi,cpsi)
        return
      else
        diser=dcmplx(rpsi,cpsi)
        do n=1,nz-1
	  z=dcmplx(dfloat(n),yz)
	  diser=diser+1.d0/z
	enddo
        return
      endif
      end
**********
      subroutine digi2(z,diser)
      implicit none
      complex*16 z,diser
      complex*16 c1,cfac,ctemp
      real*8 rn,gama,pi
      integer n
      parameter(gama=0.577215664901532860606512d0,
     *          pi=3.141592653589793238462643d0)
      parameter(c1=(1.d0,0.d0))
      cfac=-gama-1.d0/z+pi*pi/6.d0
      ctemp=cfac
      do n=1,10000
        rn=dfloat(n)
        cfac=z/(rn*rn+rn*z)-1.d0/(rn*rn)
        diser=ctemp+cfac
        if(diser.eq.ctemp) go to 200
        ctemp=diser
      enddo
c     pause 'convergence failure in digi2'
200   continue
      return
      end
**********
      complex*16 function calgam(cz)
**************************************
*  Subroutine tu calculate the logarithm of the Gamma function
*  It uses the program gammln from Numerical Recipes (2nd Ed.)
*  which is valid only when Real(cz).gt.0. Thus, the present
*  subroutine introduces appropiate formulas to calculate for
*  any value of cz.
*
*  Author: Pablof    V 1.0   9/12/1994
**************************************
      implicit none
      complex*16 cz
      complex*16 c0,c1,cz1,cgamln
      real*8 xz,yz,pi
      integer iw
      parameter(pi=3.141592654d0, c0=(0.d0,0.d0), c1=(1.d0,0.d0))
      data iw/6/
      if(cz.eq.c0) then
        write (iw,3) cz
	return
      endif
      xz=dble(cz)
      yz=dimag(cz)
      if(xz.eq.0.d0) then
        cz1=c1+cz
        calgam=cgamln(cz1)-cdlog(cz)
        return
      else if(xz.lt.0.d0) then
        if(int(xz)/dabs(xz).eq.-1.d0.and.yz.eq.0.d0) then
          write (iw,3) cz
	  calgam=c0
	  return
        endif
        cz1=c1-cz
        calgam=cdlog(pi/cdsin(pi*cz))-cgamln(cz1)
        return
      else
        calgam=cgamln(cz)
        return
      endif
3     format(2x,'z = ',2(1pe14.7,2x),'is considered as a pole')
      end
**********
      complex*16 function cgamln(cz)
      implicit none
      complex*16 cz
      complex*16 cx,cy,ctmp,cser
      real*8 rcof(6),stp
      integer j
      save rcof,stp
      data rcof,stp/76.18009172947146d0,-86.50532032941677d0,
     *              24.01409824083091d0,-1.231739572450155d0,
     *              0.1208650973866179d-2,-0.5395239384953d-5,
     *              2.5066282746310005d0/
      if(dble(cz).le.0.d0) stop 'real(z).le.0'
      cx=cz
      cy=cx
      ctmp=cx+5.5d0
      ctmp=(cx+0.5d0)*cdlog(ctmp)-ctmp
      cser=dcmplx(1.000000000190015d0,0.d0)
      do j=1,6
        cy=cy+1.d0
        cser=cser+rcof(j)/cy
      end do
      cgamln=ctmp+cdlog(stp*cser/cx)
      return
      end
c**********************************************************************
*
*    PabloF's library
*    Contents:           avint.f -----> Simpson
*                        galag.f -----> Gauss-Laguerre
*                        poch.f ------> Pochammer symbol
*
*-----------------------------------------------------------------------------*
*                                                                             *
*    This program computes the integral of y(x) in XLO < x < XUP .            *
*                                                                             *
*    X and Y are N dimensional arrays containing the abscissas x and the      *
*    ordinates y respectively.The x are usually unequally spaced and must     *
*    be given in ascending order.The y are usually data.XLO must be less      *
*    than XUP,but otherwise there are no restrictions on XLO and XUP except,  *
*    of course that of accuracy which requires that if XLO is less than x1,   *
*    it be close to x1,and if XUP is greater than xN,it be close to xN.       *
*                                                                             *
*    AVINT is adapted from P.E.Hennion,Algorithm 77,CACM 5,1962,p. 96.        *
*    (CACM = Communications of the Association for Computing Machinery)       *
*    AVINT was copied from P.J.Davis and P.Rabinowitz,NUMERICAL INTEGRATION,  *
*    Blaisdell Publishing Company,1967.                                       *
*                                                                             *
*    Some parts coded by PabloF (20/07/1999) (pablof@cab.cnea.gov.ar)         *
*-----------------------------------------------------------------------------*
      function avint(x,y,n,xlo,xup)
      implicit none
      integer n
      real*8 avint,x(n),y(n),xlo,xup
      integer i,j,ib,jm
      real*8 sum,syl,syu,x1,x2,x3,term1,term2,term3,a,b,c,ca,cb,cc

      sum=0.d0
      syl=xlo

      ib=2
      do i=1,n
        if(x(i).ge.xlo) then
          go to 170
        else
          ib=ib+1
        endif
      enddo

170   continue

      j=n
      do i=1,n
        if(xup.ge.x(j)) then
          go to 18
        else
          j=j-1
        endif
      enddo

18    j=j-1
      do jm=ib,j
        x1=x(jm-1)
        x2=x(jm)
        x3=x(jm+1)
        term1=y(jm-1)/((x1-x2)*(x1-x3))
        term2=y(jm)/((x2-x1)*(x2-x3))
        term3=y(jm+1)/((x3-x1)*(x3-x2))
        a=term1+term2+term3
        b=-(x2+x3)*term1-(x1+x3)*term2-(x1+x2)*term3
        c=x2*x3*term1+x1*x3*term2+x1*x2*term3
        if(jm.eq.ib) then
          ca=a
          cb=b
          cc=c
        else
          ca=0.5d0*(a+ca)
          cb=0.5d0*(b+cb)
          cc=0.5d0*(c+cc)
        endif
        syu=x(jm)
        sum=sum+ca*(syu**3-syl**3)/3.d0
     *         +cb*0.5d0*(syu**2-syl**2)+cc*(syu-syl)
        ca=a

        cb=b
        cc=c
        syl=syu
      enddo

      avint=sum+ca*(xup**3-syl**3)/3.d0
     *         +cb*0.5d0*(xup**2-syl**2)+cc*(xup-syl)
      return
      end
*******************************************************************************
*                                                                             *
*  This program computes the integral of the function FUN(Q) in the interval  *
*  QMIN < Q < QMAX with the GAUSS-LAGUERRE quadrature using 15, 20, 24 or 48  *
*  points.                                                                    *
*                                                                             *
*      N : number of points (15,20,24,48).                                    *
*  ALPHA : velocity parameter.                                                *
*                                                                             *
*******************************************************************************
      function galag(n,fun,alpha,qmin)
      implicit none
      integer n,i
      real*8 galag,fun,alpha,qmin,sum
      real*8 x15(15), x20(20), x24(24), x48(48)
      real*8 w15(15), w20(20), w24(24), w48(48)
      external fun

      data x15/0.0933078120,0.4926917403,1.2155954120,2.2699495262,
     *         3.6676227217,5.4253366274,7.5659162266,10.1202285680,
     *         13.1302824821,16.6544077083,20.7764788994,25.6238942267,
     *         31.4075191697,38.5306833064,48.0260855726/
      data w15/0.2395781703,0.5601008427,0.8870082629,1.22366440215,
     *         1.57444872163,1.94475197653,2.34150205664,2.77404192683,
     *         3.25564334640,3.80631171423,4.45847775384,5.27001778443,
     *         6.35956346973,8.03178763212,11.5277721009/

      data x20/0.070539890,0.37212682,0.91658210,1.7073065,2.7491993,
     *         4.0489253,5.6151750,7.4590175,9.5943929,12.038803,
     *         14.814293,17.948896,21.478788,25.451703,29.932555,
     *         35.013434,40.833057,47.619994,55.810796,66.524417/
      data w20/0.18108006,0.42255677,0.66690955,0.91535237,1.1695397,
     *         1.4313550,1.7029811,1.9870159,2.2866358,2.6058347,
     *         2.9497837,3.3253958,3.7422555,4.2142367,4.7625185,
     *         5.4217260,6.2540124,7.3873144,9.1513287,12.893389/

      data x24/0.059019852,0.31123915,0.76609691,1.4255976,2.2925621,
     *         3.3707743,4.6650837,6.1815351,7.9275392,9.9120980,
     *         12.146103,14.642732,17.417993,20.491460,23.887330,
     *         27.635937,31.776041,36.358406,41.451720,47.153106,
     *         53.608575,61.058531,69.962240,81.498279/
      data w24/0.15149441,0.35325658,0.55678456,0.76268532,0.97187263,
     *         1.1853579,1.4042656,1.6298686,1.8636351,2.1072912,
     *         2.3629059,2.6330088,2.9207576,3.2301851,3.5665734,
     *         3.9370438,4.3515312,4.8244819,5.3780221,6.0484178,
     *         6.9008984,8.0699652,9.9027933,13.820532/

      data x48/0.029811236,0.15710799,0.38626504,0.71757469,1.1513938,
     *         1.6881858,2.3285270,3.0731109,3.9227524,4.8783934,
     *         5.9411081,7.1121105,8.3927626,9.7845832,11.289259,
     *         12.908658,14.644841,16.500081,18.476882,20.577999,
     *         22.806462,25.165612,27.659128,30.291071,33.065931,
     *         35.988681,39.064849,42.300590,45.702792,49.279186,
     *         53.038498,56.990625,61.146865,65.520207,70.125706,
     *         74.980978,80.106857,85.528311,91.275708,97.386668,
     *         103.90883,110.90422,118.45643,126.68343,135.76259,
     *         145.98643,157.91561,172.99633/
      data w48/0.076509177,0.17817675,0.28018214,0.38249636,0.48521906,
     *         0.58846091,0.69233697,0.79696649,0.90247346,1.0089874,
     *         1.1166442,1.2255875,1.3359696,1.4479526,1.5617104,
     *         1.6774302,1.7953145,1.9155834,2.0384773,2.1642601,
     *         2.2932231,2.4256893,2.5620186,2.7026147,2.8479326,
     *         2.9984884,3.1548711,3.3177578,3.4879319,3.6663079,
     *         3.8539614,4.0521703,4.2624683,4.4867171,4.7272065,
     *         4.9867940,5.2691058,5.5788333,5.9221856,6.3076046,
     *         6.7469484,7.2575591,7.8661312,8.6166144,9.5883932,
     *         10.945845,13.115122,17.849139/

      sum=0.d0
      if(n.le.15) then
        do 15 i=1,15
15      sum=sum+w15(i)*fun(x15(i)/alpha+qmin)
        galag=sum/alpha
        return
      else if(n.le.20) then
        do 20 i=1,20
20      sum=sum+w20(i)*fun(x20(i)/alpha+qmin)
        galag=sum/alpha
        return
      else if(n.le.24) then
        do 24 i=1,24
24      sum=sum+w24(i)*fun(x24(i)/alpha+qmin)
        galag=sum/alpha
        return
      else
        do 48 i=1,48
48      sum=sum+w48(i)*fun(x48(i)/alpha+qmin)
        galag=sum/alpha
        return
      end if
      end
c*********
      subroutine poch(n,cx,cpoch)
*
*  Routine  to calculate  the POCHAMMER symbol
*  (A)n = A(A+1)(A+2)...(A+n-1) , (A)o = 1
*  Author: PabloF 3/2/94    Version 1.1
*  New version extended to complex arguments
*
      integer n
      complex*16 cx,cpoch
      integer i,nmax
      parameter(nmax=100)
      complex*16 cvpo(0:nmax),cnil,cone
      parameter(cnil=(0.d0,1.d0),cone=(1.d0,1.d0))
      if(n.eq.0) then
        cpoch=cone
        return
      else
        if(cx.eq.cnil) then
	  cpoch=cnil
	  return
	else
	  cvpo(0)=cone
          do i=0,n-1
            cvpo(i+1)=cvpo(i)*(cx+dcmplx(dfloat(i),0.d0))
          enddo
          cpoch=cvpo(n)
          return
	end if
      end if
      end
