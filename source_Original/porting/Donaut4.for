      subroutine donaut(QP,Z1,AP,EP,Z2,AC,RHO,EPM,iflag,E1,S1,E2,S2,
     s      istp,iprt,ilgn)
      real*8 zp8,zt8,E8
	integer ninis
	dimension sec(34),cor(48),secs(48),seci(48),SeSE(66),StSE(12)
	common/seceff/sec,cor,secs
	common/SecSE/SeSE,StSE
	common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4
	common/tol/ep0,ep1,erel,erabs
	common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn
	common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
	common/don/zp,E,zt
	common/corr/ibin
	character kk*2,ll*2,kc*2, nomf*24
c_______________________________________________________________________
c
c	table of cross sections (units 10-20cm2):
c	indices			secs			indices			sec
c
c		1			MeC1s	SEIK		1			C1s
c		2			MeC2s	SEIK		2			C2s
c		3			MeC2p	SEIK		3			C2p
c		4			MeC3s	SEIK		4			C3s
c		5			MeC3p	SEIK		5			C3p
c		6			MeC3d	SEIK		6			C3d
c		7			MeC4	SEIK		7			C4
c
c		8			ReC1s				8			Ion1s
c		9			ReC2s				9			Ion2s
c		10			ReC2p				10			Ion2p
c		11			ReC3s				11			Ion3s
c		12			ReC3p				12			Ion3p
c		13			ReC3d				13			Ion3d
c		14			ReC4				14			Ion4
c
c		15			Ion1s	CDW	EIS		15			Se1s2s
c		16			Ion2s	CDW	EIS		16			Se1s2p
c		17			Ion2p	CDW	EIS		17			Se1s3s
c		18			Ion3s	Scaling		18			Se1s3p
c		19			Ion3p	Scaling		19			Se1s3d
c		20			Ion3d	Scaling		20			Se1s4
c		21			Ion4	Scaling	
c
c     Scaling means Ion(3l,Z) = Ion(1s,Z/3) and ion(4,Z) = Ion(2p,Z/2)
c										
c		22			Se1s5	SE  		21			Se2s3s
c		23			Se2s5	SE  		22			Se2s3p
c		24			Se2p5	SE  		23			Se2s3d
c		25			Se3s5	SE          24			Se2s4
c		26			Se3p5	SE		
c		27			Se3d5	Scaling SE	25			Se2p3s
c		28			Se456	Scaling SE	26			Se2p3p
c										27			Se2p3d
c		29			Se1s2s	SE          28			Se2p4
c		30			Se1s2p	SE			
c		31			Se1s3s	SE			29			Se3s4
c		32			Se1s3p	SE			30			Se3p4
c		33			Se1s3d	SE	        31			Se3d4
c		34			Se1s4	SE			
c										
c		35			Se2s3s	SE			32			Se2s2p
c		36			Se2s3p	SE	        33			Se3s3p
c		37			Se2s3d	SE          34			Se3p3d
c		38			Se2s4	SE          
c
c		39			Se2p3s	SE
c		40			Se2p3p	SE
c		41			Se2p3d	SE
c		42			Se2p4	SE
c
c		43			Se3s4	SE
c		44			Se3p4	SE
c		45			Se3d4	Scaling SE
c
c		46			Se2s2p
c		47			Se3s3p
c		48			Se3p3d
c_________________________________________________________________________________
c...... ibin=0 empirical saturation correction for born .....
c...... ibin=1 binding correction included in born (not recommended) .....
c...... ibin=2 no empirical correction and no binding correction
c...... correction are only for PWBA (BORN1) cross sections
c...... as calculated in PION, SEXI, SEX2 and SEX3,
c...... are only used for the evolution of cross sections with effective charge,
c...... and should not make much a change 
c
c...... Excitation to n=5 (secs(i), i=22,28) is not added to ionisation and set 
c...... to 0
c
c...... CDW is not used for capture cross sections 
c________________________________________________________________________________
	ibin=0
	iflag=1
	iter=0
	y1s=0.
	y2s=0.
	y2p=0.
	y3s=0.
	y3p=0.
	y3d=0.
	y4=0.
	yl=0.
	ym=0.
	ymp=0.
	ymm=0.
	ykm=0.
	yl1m=0.
	yl2m=0.
	ym1m=0.
	ym2m=0.
	tetak=1.
	tetal1=1.
	tetal2=1.
	tetam1=1.
	tetam2=1.
	tetan=1.
	nomf = 'files\etadon.etacha'
	open(UNIT=13,ERR=1003,STATUS='old',FILE=nomf)
	read(13,11,err=1000) QP,ZP,AP,zt,AC
11	format(4(1x,f4.0),1x,f6.2/)
	read(13,12,err=1000) E,RHO,EPM
12	format(2(1x,F8.3),1x,F10.3,/)
	read(13,15,err=1000) E1,S1,E2,S2,ISTP
15	format(4(1x,F8.3),1xI2,/)
	read(13,16,err=1000) ep0,ep1,erel,erabs
16	format(2(1x,g10.3),2(1x,g11.4),/)
	read(13,17,err=1000) iprt,ilgn
17	format(2(1x,I2),/)
	CLOSE(UNIT=13)
	write(6,*)'      Data used in previous calculation :'
	write(6,*)' '
500	write(6,*)'      General data :'
	write(6,*)' '
	write (6,20)zp,qp,ap,e,zt,ac,rho
20	format(2x26hPROJECTILE: atomic number=,f4.0,2x11hincident ch,
     s	5harge=,f4.0,2x12hatomic mass=f4.0,/,12x16hincident energy=,
     s	f8.3,6h MeV/u,/,6x22hTARGET: atomic number=,f4.0,2x6hatomic,
     s	6h mass=,f4.0,2x8hdensity=,f6.3,6h g/cm3)
	write (6,21)epm,ep0,ep1,erabs,erel
21	format(2x35hmaximum target thickness (mg/cm2) =,f10.3,/,
     s	2x23hminimum step (µg/cm2) =,g10.3,2x13hmaximum step=,g10.3,
     s	/,2x35hnumerical uncertainties on output :,/,2x9habsolute=,
     s	g11.4,2x9hrelative=,g11.4,/)
	write (6,*) ' Want to change any of these value? (No/y)'
	write (6,*) ' (Type ''y'' to change'
	write (6,*) '  Type ''return'' not to change) '
	read (5,'(A)') KK
c...........................................
	if ((KK.eq.'Y').or.(KK.eq.'y')) then
501	  write (6,22) zp
22	format(2x26hProjectile atomic number =,f4.0," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    iter=0
	    write (6,*) ' New value?'
	    read (5,*) Zp
	  end if
	  write (6,23) qp
23	format(2x19hProjectile charge =,f4.0," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) QP
	  end if
	  write (6,24) Ap
24	format(2x17hProjectile mass =,f4.0," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    iter=0
	    write (6,*) ' New value? '
	    read (5,*) AP
	  end if
	  write (6,25) E
25	format(2x25hIncident energy (MeV/u) =,f8.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    iter=0
	    write (6,*) ' New value? '
	    read (5,*) E
	  end if
	  write (6,26) Zt
26	format(2x22hTarget atomic number =,f4.0," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    iter=0
	    write (6,*) ' New value? '
	    read (5,*) zt
	  end if
	  write (6,27) Ac
27	format(2x20hTarget atomic mass =,f4.0," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    iter=0
	    write (6,*) ' New value? '
	    read (5,*) Ac
	  end if
	  write (6,28) Rho
28	format(2x24hTarget density (g/cm3) =,f6.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) RHO
	  end if
	  write (6,29) Epm
29	format(2x35hmaximum target thickness (mg/cm2) =,f10.3,
     s	" Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) EPM
	  end if
	  write (6,30) Ep0
30	format(2x23hminimum step (µg/cm2) =,g10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Ep0
	  end if
	  write (6,31) Ep1
31	format(2x23hmaximum step (µg/cm2) =,g10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Ep1
	  end if
	write(6,32) Erabs,Erel
32	format(2x46hMaximun numerical uncertainties on populations,
     s	28h P(i) (=erel*(Max(P,erabs))),/,2x9habsolute=,
     s	g11.4,2x9hrelative=,g11.4," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value for erabs? '
	    read (5,*) ERABS
	    write (6,*) ' New value for erel? '
	    read (5,*) EREL
	  end if
	  write(6,*)' New values :'
	  write(6,*)' '
	  write (6,20)zp,qp,ap,e,zt,ac,rho
	  write (6,21)epm,ep0,ep1,erabs,erel
	  write (6,*) ' Want to change again any of these value?',
     s	  ' (No/y) '
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    go to 501
	  end if
	end if
c.......................................................................
	write(6,*)' Do you want this program to modify cross '
	write(6,*)'sections when the projectile ion loss'
	write(6,*)' energy in thick targets? (No/y) '
	  read (5,'(A)') KK
c.................................................
	if ((KK.eq.'Y').or.(KK.eq.'y')) then
	istp=1
	write(6,*)' This program needs two values of energy',
     s	' in order to perform this task'
	write(6,*)' (You can choose one value close to the initial',
     s' energy and one close to '
	write(6,*)' the final one)'
	write (6,*) ' Previous values :'
cc.........
	write (6,33)E1,S1
	write (6,34)E2,S2
33	format(2x12hE1 (MeV/u) =,f8.3,2x17hS1 (MeV/mg/cm2) =,f8.3)
34	format(2x12hE2 (MeV/u) =,f8.3,2x17hS2 (MeV/mg/cm2) =,f8.3)
	write (6,*) ' Change these values ? (No/y) '
	read (5,'(A)') KK
cc.........
	if ((KK.eq.'Y').or.(KK.eq.'y')) then
502	  write (6,35)E1,S1
35	  format(2x12hE1 (MeV/u) =,f8.3,2x17hS1 (MeV/mg/cm2) =,f8.3,
     s	  " Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value for E1 (MeV/u)? '
	    Write (6,*)' (If your input for E1 is 0 there will be NO',
     s	  ' stopping correction)'
	    read (5,*) E10
	    if (E10.eq.0.) then
	    istp=0
	    go to 503
	    else
	    E1=E10
		call zstop(zp,E1,zt,S1)
	    end if
		write (6,*) ' Calculated value : '
		write (6,351)E1,S1
351	  format(2x12hE1 (MeV/u) =,f8.3,2x17hS1 (MeV/mg/cm2) =,f8.3)
	    write (6,*) ' Change value for S1? (No/y)'
		read (5,'(A)') KK
		if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value for S1 (MeV/mg/cm2)? '
	    read (5,*) S1
		end if
	  end if
	  write (6,36)E2,S2
36	  format(2x12hE2 (MeV/u) =,f8.3,2x17hS2 (MeV/mg/cm2) =,f8.3,
     s	  " Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value for E2 (MeV/u)? '
	    read (5,*) E2
		call zstop(zp,E2,zt,S2)
		write (6,*) ' Calculated value : '
		write (6,361)E2,S2
361	  format(2x12hE2 (MeV/u) =,f8.3,2x17hS2 (MeV/mg/cm2) =,f8.3)
	    write (6,*) ' Change value for S2? (No/y)'
		read (5,'(A)') KK
		if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value for S2 (MeV/mg/cm2)? '
	    read (5,*) S2
		end if
	  end if
503	  write(6,*)' new values for stopping power correction:'
	  if (istp.eq.0) then
	  write(6,*)' no stopping power correction ...'
	  else
	  write (6,33)E1,S1
	  write (6,34)E2,S2
	  end if
	  write (6,*) ' Change these values ? (No/y) '
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	  go to 502
	  end if
	end if
	else
	  istp=0
	end if
c---------------------------------------------------
c	CROSS SECTIONS 
c---------------------------------------------------
	write (6,*) ' '
	write (6,*) ' please, wait......'
	write (6,*) ' Hydrogenic (reference) cross sections',
     s	' calculation ....'
	if (iter.eq.0) then
	  do 50 I=1,48
	    cor(I)=1.
50	  continue
	  call CSEC(E)
	  do 51 I=1,48
	    seci(I)=secs(I)
51	  continue
cc
c     Preliminary calculation of scaling factors for excitation cross sections
c     of d and f states in SE approximation @ v=zp
c
      Ev=E
      bet=zp/137.036
      E=931.5*(1/(1-bet*bet)**0.5-1)
      write (6,*)' E=',E
	call sex3(e3s4,e3p4,e3d4,e3s5,e3p5,e3d5,e4s5,e4p5,e45,e456)
	R3d4=(e3d4)/(e3s4+e3p4)
	R3d5=(e3d5)/(e3s5+e3p5)
c     e456 is not corrected for 1/n**3 law anymore...
	R45=e456/(e4s5+e4p5)
	E=Ev
cc
	write (6,*) ' '
	write(6,*)' CDWEIS Ionisation Cross sections'
	ninis=10
	call tceis(tot,ninis)
	secs(15)=tot*1.e20
	ninis=20
	call tceis(tot,ninis)
	secs(16)=tot*1.e20
	ninis=21
	call tceis(tot,ninis)
	secs(17)=tot*1.e20
c
      zp3=zp
      zp=zp3/3.
	ninis=10
	call tceis(tot,ninis)
	secs(18)=tot*1.e20
	secs(19)=tot*1.e20
	secs(20)=tot*1.e20
	zp=zp3
c	
      zp2=zp
      zp=zp2/2.
	ninis=21
	call tceis(tot,ninis)
	secs(21)=tot*1.e20
	zp=zp2
c	
	cor(15)=secs(15)/seci(15)
	cor(16)=secs(16)/seci(16)
	cor(17)=secs(17)/seci(17)
	cor(18)=secs(18)/seci(18)
	cor(19)=secs(19)/seci(19)
	cor(20)=secs(20)/seci(20)
	cor(21)=secs(21)/seci(21)     	
c      	  
	write(6,*)' SE Excitation Cross sections'
	zp8=dble(zp)
	zt8=dble(zt)
	E8=dble(E)
	call SEnlm(zp8,zt8,E8)
c
c     Excitation to n>=5 set to 0 in the present version
c	
c     1s - 5
	secs(22)=StSE(6)*1.e20*0.
c     2s,2p - 5	
	secs(23)=StSE(7)*1.e20*0.
	secs(24)=StSE(8)*1.e20*0.
c     3s,3p,3d - 5	
	secs(25)=StSE(9)*1.e20*0.
	secs(26)=StSE(10)*1.e20*0.
	secs(27)=R3d5*(secs(25)+secs(26))
c     R45 does not include correction for the sum on 1/n**3	anymore
	secs(28)=R45*(StSE(11)+StSE(12))*1.e20*0.
c	 
	secs(29)=SeSE(1)*1.e20
	secs(30)=SeSE(2)*1.e20
	secs(31)=SeSE(3)*1.e20
	secs(32)=SeSE(4)*1.e20
	secs(33)=SeSE(5)*1.e20
	secs(34)=StSE(1)*1.e20
	secs(35)=SeSE(15)*1.e20
	secs(36)=SeSE(16)*1.e20
	secs(37)=SeSE(17)*1.e20
	secs(38)=StSE(2)*1.e20
	secs(39)=SeSE(27)*1.e20
	secs(40)=SeSE(28)*1.e20
	secs(41)=SeSE(29)*1.e20
	secs(42)=StSE(3)*1.e20
	secs(43)=StSE(4)*1.e20
	secs(44)=StSE(5)*1.e20
	secs(45)=R3d4*(secs(43)+secs(44))
c
c     Excitation to n=5 and 6 set to 0	
c
	cor(22)=0.
	cor(23)=0.
	cor(24)=0.
	cor(25)=0.
	cor(26)=0.
	cor(27)=0.
	cor(28)=0.
c	
	cor(29)=secs(29)/seci(29)
	cor(30)=secs(30)/seci(30)
	cor(31)=secs(31)/seci(31)
	cor(32)=secs(32)/seci(32)
	cor(33)=secs(33)/seci(33)
	cor(34)=secs(34)/seci(34)
	cor(35)=secs(35)/seci(35)
	cor(36)=secs(36)/seci(36)
	cor(37)=secs(37)/seci(37)
	cor(38)=secs(38)/seci(38)
	cor(39)=secs(39)/seci(39)
	cor(40)=secs(40)/seci(40)
	cor(41)=secs(41)/seci(41)
	cor(42)=secs(42)/seci(42)
	cor(43)=secs(43)/seci(43)
	cor(44)=secs(44)/seci(44)
	cor(45)=secs(45)/seci(45)
c
c     note that secs(22) -> secs(28) are set to 0...
c	  
	sec(8)=secs(15)+secs(22)
	sec(9)=secs(16)+secs(23)
	sec(10)=secs(17)+secs(24)
	sec(11)=secs(18)+secs(25)
	sec(12)=secs(19)+secs(26)
	sec(13)=secs(20)+secs(27)
	sec(14)=secs(21)+secs(28)
c	
	sec(15)=secs(29)
	sec(16)=secs(30)
	sec(17)=secs(31)
	sec(18)=secs(32)
	sec(19)=secs(33)
	sec(20)=secs(34)
	sec(21)=secs(35)
	sec(22)=secs(36)
	sec(23)=secs(37)
	sec(24)=secs(38)
	sec(25)=secs(39)
	sec(26)=secs(40)
	sec(27)=secs(41)
	sec(28)=secs(42)
	sec(29)=secs(43)
	sec(30)=secs(44)
	sec(31)=secs(45)
	  iter=1
	end if
	write (6,*) ' '
	write (6,40) (sec(I),I=1,33)
40	FORMAT(2x49hCapture cross sections into fully stripped projec,
     s	8htile ...,/,2x43hIonization and excitation cross sections of,
     s	28h one electron projectile ...,/,2x22h(corrected for "satura,
     s	36htion" and screening + antiscreening),/,2x14h(all cross sec,
     s	21htions in 10E-20 cmý ),/,2x11hCAPTURE TO:,5x12h(capture inc,
     s	17hludes MEC + REC ),/,21x4h1s =,G10.3,2x4h2s =,G10.3,2x4h2p =,
     s	G10.3,/,21x4h3s =,G10.3,2x4h3p =,G10.3,2x4h3d =,G10.3,/,
     s	21x7h(n=4) =,G10.3,/,
     s   2x57hIONIZATION OF:  (ionization includes excitation to n > 4),
     s   /,21x4h1s =,G10.3,2x4h2s =,G10.3,2x4h2p =,G10.3,/,21x4h3s =,
     s	G10.3,2x4h3p =,G10.3,2x4h3d =,G10.3,/,21x7h(n=4) =,G10.3,/,
     s	2x17h1s EXCITATION TO:,2x4h2s =,
     s	G10.3,2x4h2p =,G10.3,/,21x4h3s =,G10.3,2x4h3p =,G10.3,2x4h3d =,
     s	G10.3,/,21x7h(n=4) =,G10.3,/,
     s	2x17h2s EXCITATION TO:,2x4h3s =,G10.3,2x4h3p =,
     s	G10.3,2x4h3d =,G10.3,/,21x7h(n=4) =,G10.3,/,
     s	2x17h2p EXCITATION TO:,2x4h3s =,G10.3,
     s	2x4h3p =,G10.3,2x4h3d =,G10.3,/,21x7h(n=4) =,G10.3,/,
     s	2x23hEXCITATION TO n=4 from:,2x4h3s =,G10.3,2x4h3p =,
     s	G10.3,2x4h3d =,G10.3,/,
     s	2x22hINTRASHELL EXCITATION:,
     s	2x10h2s to 2p =,G10.3,2x10h3s to 3p =,G10.3)
	Write(6,*)'  During charge state calculation, cross sections',
     s	' will be periodicaly adjusted as a function of mean charge',
     s	' state (if necessary) and as a function of energy  (if',
     s	' desired). If however you are not pleased with some of the',
     s	' above values,you may now enter your own ones',
     s	'(fully stripped or one electron ion cross sections)',
     s	'   The program will in this case scale his calculat',
     s	'ions to yours.'
	Write(6,*)'     Do you want to enter such corrected'
	Write(6,*)' values ? (No/y) '
	read (5,'(A)') KK
c ...................... debut des modifs ..........................
	if ((KK.eq.'Y').or.(KK.eq.'y')) then
504	  write (6,*) ' Enter your cross sections in 10E-20 cm2'
	  Write(6,*)'  CAPTURE to substates of fully stripped ions :'
	Write(6,*)' Want to change ? (No/y) '
	read (5,'(A)') KC
	if ((KC.eq.'Y').or.(KC.eq.'y')) then
	  write (6,101) secs(1)
101	  Format(1x18h1s capture (MEC) =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	   write (6,*) ' New value? '
	   read (5,*) C1s
	   cor(1)=C1s/seci(1)
	  end if
	  write (6,102) secs(2)
102	  Format(1x18h2s capture (MEC) =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) C2s
	    cor(2)=C2s/seci(2)
	  end if
	  write (6,103) secs(3)
103	  Format(1x18h2p capture (MEC) =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) C2p
	    cor(3)=C2p/seci(3)
	  end if
	  write (6,104) secs(4)
104	  Format(1x18h3s capture (MEC) =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) C3s
	    cor(4)=C3s/seci(4)
	  end if
	  write (6,105) secs(5)
105	  Format(1x18h3p capture (MEC) =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) C3p
	    cor(5)=C3p/seci(5)
	  end if
	  write (6,106) secs(6)
106	  Format(1x18h3d capture (MEC) =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) C3d
	    cor(6)=C3d/seci(6)
	  end if
	  write (6,107) secs(7)
107	  Format(1x19hn=4 capture (MEC) =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	   write (6,*) ' New value? '
	   read (5,*) C4
	   cor(7)=C4/seci(7)
	  end if
	  write (6,108) secs(8)
108	  Format(1x8h1s REC =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	   write (6,*) ' New value? '
	   read (5,*) SRK
	   cor(8)=SRK/seci(8)
	  end if
	  write (6,109) secs(9)
109	  Format(1x8h2s REC =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) SRls
	    cor(9)=SRls/seci(9)
	  end if
	  write (6,110) secs(10)
110	  Format(1x8h2p REC =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) SRlp
	    cor(10)=SRlp/seci(10)
	  end if
	  write (6,111) secs(11)
111	  Format(1x8h3s REC =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) SR3s
	    cor(11)=SR3s/seci(11)
	  end if
	  write (6,112) secs(12)
112	  Format(1x8h3p REC =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) SR3p
	    cor(12)=SR3p/seci(12)
	  end if
	  write (6,113) secs(13)
113	  Format(1x8h3d REC =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) SR3d
	    cor(13)=SR3d/seci(13)
	  end if
	  write (6,114) secs(14)
114	  Format(1x9hn=4 REC =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) SR4
	    cor(14)=SR4/seci(14)
	  end if
	end if
	  Write(6,*)'  IONIZATION of hydrogenlike ions :'
	Write(6,*)' Want to change ? (No/y) '
	read (5,'(A)') KC
	if ((KC.eq.'Y').or.(KC.eq.'y')) then
	  write (6,115) secs(15)
115	  Format(1x15h1s Ionization =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Dk
	    cor(15)=Dk/seci(15)
	  end if
	  write (6,116) secs(16)
116	  Format(1x15h2s Ionization =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Dl1
	    cor(16)=Dl1/seci(16)
	  end if
	  write (6,117) secs(17)
117	  Format(1x15h2p Ionization =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Dl2
	    cor(17)=Dl2/seci(17)
	  end if
	  write (6,118) secs(18)
118	  Format(1x15h3s Ionization =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Dm1
	    cor(18)=Dm1/seci(18)
	  end if
	  write (6,119) secs(19)
119	  Format(1x15h3p Ionization =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Dm2
	    cor(19)=Dm2/seci(19)
	  end if
	  write (6,120) secs(20)
120	  Format(1x15h3d Ionization =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Dm3
	    cor(20)=Dm3/seci(20)
	  end if
	  write (6,121) secs(21)
121	  Format(1x16hn=4 Ionization =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Dn
	    cor(21)=Dn/seci(21)
	  end if
	end if
	  write (6,*)' excitation to n>4 to be added to ionization :'
	Write(6,*)' Want to change ? (No/y) '
	read (5,'(A)') KC
	if ((KC.eq.'Y').or.(KC.eq.'y')) then
	  write (6,122) secs(22)
122	  Format(1x22h1s -> n>4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) es5
	    cor(22)=es5/seci(22)
	  end if
	  write (6,123) secs(23)
123	  Format(1x22h2s -> n>4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e2s5
	    cor(23)=e2s5/seci(23)
	  end if
	  write (6,124) secs(24)
124	  Format(1x22h2p -> n>4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e2p5
	    cor(24)=e2p5/seci(24)
	  end if
	  write (6,125) secs(25)
125	  Format(1x22h3s -> n>4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e3s5
	    cor(25)=e3s5/seci(25)
	  end if
	  write (6,126) secs(26)
126	  Format(1x22h3p -> n>4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e3p5
	    cor(26)=e3p5/seci(26)
	  end if
	  write (6,127) secs(27)
127	  Format(1x22h3d -> n>4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e3d5
	    cor(27)=e3d5/seci(27)
	  end if
	  write (6,128) secs(28)
128	  Format(1x23hn=4 -> n>4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e45
	    cor(28)=e45/seci(28)
	  end if
	end if
	  Write(6,*)'  EXCITATION of hydrogenlike ions :'
	  write(6,*)' Excitation from 1s :'
	Write(6,*)' Want to change ? (No/y) '
	read (5,'(A)') KC
	if ((KC.eq.'Y').or.(KC.eq.'y')) then
	  write (6,129) secs(29)
129	  Format(1x21h1s -> 2s excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) E2s
	    cor(29)=E2s/seci(29)
	  end if
	  write (6,130) secs(30)
130	  Format(1x21h1s -> 2p excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) E2p
	    cor(30)=E2p/seci(30)
	  end if
	  write (6,131) secs(31)
131	  Format(1x21h1s -> 3s excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) E3s
	    cor(31)=E3s/seci(31)
	  end if
	  write (6,132) secs(32)
132	  Format(1x21h1s -> 3p excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) E3p
	    cor(32)=E3p/seci(32)
	  end if
	  write (6,133) secs(33)
133	  Format(1x21h1s -> 3d excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) E3d
	    cor(33)=E3d/seci(33)
	  end if
	  write (6,134) secs(34)
134	  Format(1x22h1s -> n=4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) es4
	    cor(34)=es4/seci(34)
	  end if
	end if
	  write(6,*)' Excitation from n=2 :'
	Write(6,*)' Want to change ? (No/y) '
	read (5,'(A)') KC
	if ((KC.eq.'Y').or.(KC.eq.'y')) then
	  write (6,135) secs(35)
135	  Format(1x21h2s -> 3s excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) s3s
	    cor(35)=s3s/seci(35)
	  end if
	  write (6,136) secs(36)
136	  Format(1x21h2s -> 3p excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) s3p
	    cor(36)=s3p/seci(36)
	  end if
	  write (6,137) secs(37)
137	  Format(1x21h2s -> 3d excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) s3d
	    cor(37)=s3d/seci(37)
	  end if
	  write (6,138) secs(38)
138	  Format(1x22h2s -> n=4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e2s4
	    cor(38)=e2s4/seci(38)
	  end if
	  write (6,139) secs(39)
139	  Format(1x21h2p -> 3s excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) p3s
	    cor(39)=p3s/seci(39)
	  end if
	  write (6,140) secs(40)
140	  Format(1x21h2p -> 3p excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) p3p
	    cor(40)=p3p/seci(40)
	  end if
	  write (6,141) secs(41)
141	  Format(1x21h2p -> 3d excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) p3d
	    cor(41)=p3d/seci(41)
	  end if
	  write (6,142) secs(42)
142	  Format(1x22h2p -> n=4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e2p4
	    cor(42)=e2p4/seci(42)
	  end if
	end if
	  write(6,*)' Excitation from n=3 :'
	Write(6,*)' Want to change ? (No/y) '
	read (5,'(A)') KC
	if ((KC.eq.'Y').or.(KC.eq.'y')) then
	  write (6,143) secs(43)
143	  Format(1x22h3s -> n=4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e3s4
	    cor(43)=e3s4/seci(43)
	  end if
	  write (6,144) secs(44)
144	  Format(1x22h3p -> n=4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e3p4
	    cor(44)=e3p4/seci(44)
	  end if
	  write (6,145) secs(45)
145	  Format(1x22h3d -> n=4 excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) e3d4
	    cor(45)=e3d4/seci(45)
	  end if
	end if
	  Write(6,*)'  INTRASHELL EXCITATION of hydrogenlike ions :'
	Write(6,*)' Want to change ? (No/y) '
	read (5,'(A)') KC
	if ((KC.eq.'Y').or.(KC.eq.'y')) then
	  write (6,146) secs(46)
146	  Format(1x21h2s -> 2p excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Esp
	    cor(46)=Esp/seci(46)
	  end if
	  write (6,147) secs(47)
147	  Format(1x21h3s -> 3p excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Esp3
	    cor(47)=Esp3/seci(47)
	  end if
	  write (6,148) secs(48)
148	  Format(1x21h3p -> 3d excitation =,G10.3," Change? (No/y) ")
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    write (6,*) ' New value? '
	    read (5,*) Epd3
	    cor(48)=Epd3/seci(48)
	  end if
	end if
c..........................................
	  do 52 I=1,48
	  secs(I)=cor(I)*seci(I)
52	  continue
	  do 53 I=1,7
c	add MEC and REC
	  sec(i)=secs(i)+secs(i+7)
c	add ionisation and nl -> 5 excitation
	  sec(i+7)=secs(i+14)+secs(i+21)
53	  continue
c	shift index for other processes (excitation)
	do 54 i=1,20
	  sec(i+14)=secs(i+28)
54    continue
c..........................................
	  write (6,*) ' '
	  write (6,40) (sec(I),I=1,33)
	  write (6,*) ' Want to change again any of these value? '
	  write (6,*) '(No/y) '
	  read (5,'(A)') KK
	  if ((KK.eq.'Y').or.(KK.eq.'y')) then
	    go to 504
	  end if
	end if
c ...................... fin des modifs ..........................
	write (6,*) ' Want to recheck all input values? (No/y) '
	read (5,'(A)') KK
	if ((KK.eq.'Y').or.(KK.eq.'y')) then
	  go to 500
	end if
	open(UNIT=13,ERR=1001,STATUS='old',FILE='files\etadon.etacha')
	write(13,11,err=1001) QP,ZP,AP,ZT,AC
	write(13,12,err=1001) E,RHO,EPM
	write(13,15,err=1001) E1,S1,E2,S2,ISTP
	write(13,16,err=1001) ep0,ep1,erel,erabs
	write(13,17,err=1001) iprt,ilgn
	close(UNIT=13)
	z1=zp
	ep=e
	z2=zt
	Epm=1000.*EPM
	go to 1002
1003	write(6,*) 'Problem in opening data file... check', nomf
	iflag=2
	goto 1002
1000	write(6,*) 'Problem in reading data... check', nomf
	iflag=3
	goto 1002
1001	write(6,*) 'Problem in writing data... check', nomf
	iflag=2
1002	continue
	return
	end
cc      ************** fin de donnees ***************************
	subroutine CSEC(EF)
	dimension sec(34),cor(48),secs(48)
	common/seceff/sec,cor,secs
	common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4
	common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn
	common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn
	common/don/zp,E,zt
	common/corr/ibin
	QM=Zp-(y1s+yl+ym+yn)
	E=EF
	VC=(1-(1+E/931.5)**(-2))**(0.5)
	vu=vc*137.036
	write(6,2000)E,QM
2000	format(33h Cross sections calculation at E=,g12.5,7h(MeV/u),
     s	3x8h, Qmean=,f7.3)
	write(6,*)'ionization'
	call PION(dk,dl1,dl2,dm1,dm2,dn)
	write(6,*)'capture'
	call seik(ck,cl,cm,cn,ct,cs)
	sig14=ck+cl+cm+cn
	if (sig14.le.cs) then
	sig530=cs-sig14
	else
c
c     ct is the sum of seik from 5 to 30
c	
	sig530=ct
	end if
	call rec(srk,srls,srlp)
	srm=8./27.*(srls+srlp)
	stm=cm+srm
	srn=1./8.*(srls+srlp)
c	MEC
	secs(1)=ck
	secs(2)=.25*cl
	secs(3)=.75*cl
	secs(4)=cm*1./9.
	secs(5)=cm*3./9.
	secs(6)=cm*5./9.
c************************	
	secs(7)=cn
c	if (zp.gt.60.) then
c	secs(7)=secs(7)+sig530
c	end if
c************************	
c	REC
	secs(8)=srk
	secs(9)=srls
	secs(10)=srlp
	secs(11)=srls*8./27.
	secs(12)=srlp*8./27.
	secs(13)=srlp*8./45.
	secs(14)=srn
c	EXCITATION
	write(6,*)'excitation 1s-nl'
	call sexi(e2s,e2p,e3s,e3p,e3d,es4,es5)
	write(6,*)'excitation 2l-nl'
	call sex2(s3s,s3p,s3d,p3s,p3p,p3d,e2s4,e2p4,e2s5,e2p5)
	write(6,*)'excitation 3l-nl'
	call sex3(e3s4,e3p4,e3d4,e3s5,e3p5,e3d5,e4s5,e4p5,e45,e456)
	write(6,*)'excitation nl-nl'
	call snl(esp,esp3,epd3)
	if ((ibin.eq.1).or.(ibin.eq.2)) then
cc.......... binding correction included (1) or no correction (2) ..........
	secs(15)=Dk
	secs(16)=Dl1
	secs(17)=Dl2
	secs(18)=Dm1
	secs(19)=Dm1
	secs(20)=Dm2
	secs(21)=Dn
cc excitation to n=5 will be added to ionization
c     set to 0
	secs(22)=es5*0.
c	(es5 includes a sum up to infinity)
c	secs(23)=e2s5*3.05
c	secs(24)=e2p5*3.05
	secs(23)=e2s5*0.
	secs(24)=e2p5*0.
c	(theoretical factor in 1/n**3 = 3.049358 for n=5)
c	secs(25)=e3s5*2.5
c	secs(26)=e3p5*2.5
c	secs(27)=e3d5*2.5
	secs(25)=e3s5*0.
	secs(26)=e3p5*0.
	secs(27)=e3d5*0.
c	(3 ->5 => 2.5 instead of 3.05)
	secs(28)=e456*0.
c	(e456=e45+e46 instead of 
c	e456=e45+2.5*e46 -> 2.5 au lieu de 3.5)
c.....................
cc excitation from 1s
	secs(29)=e2s
	secs(30)=e2p
	secs(31)=e3s
	secs(32)=e3p
	secs(33)=e3d
	secs(34)=es4
cc excitation from n=2
	secs(35)=s3s
	secs(36)=s3p
	secs(37)=s3d
	secs(38)=e2s4
	secs(39)=p3s
	secs(40)=p3p
	secs(41)=p3d
	secs(42)=e2p4
cc excitation from n=3 to n=4
	secs(43)=e3s4
	secs(44)=e3p4
	secs(45)=e3d4
c...........................
	else
cc....... empirical saturation correction (09/09/94) ...........
cc    not used anymore since a long time ...
	cse=1.06
	csi=0.735
	zs1=(10.96*vu**2*exp(0.111*zp)/(zp)**1.946)**0.5
	zs2=(10.96*vu**2*exp(0.111*zp/1.5)/(zp/2.)**1.946)**0.5
	zs3=(10.96*vu**2*exp(0.111*zp/1.5)/(zp/3.)**1.946)**0.5
	sate1=(zs1/zt)**2*exp(-cse*zt/vu**2.1)
	sate1=sate1*(1.-exp(-(zt/zs1)**2.5))**0.8
	sate2=(zs2/zt)**2*exp(-cse*zt/vu**2.1)
	sate2=sate2*(1.-exp(-(zt/zs2)**2.5))**0.8
	sate3=(zs3/zt)**2*exp(-cse*zt/vu**2.1)
	sate3=sate3*(1.-exp(-(zt/zs3)**2.5))**0.8
	sati1=(zs1/zt)**2*exp(-csi*zt/vu**2.1)
	sati1=sati1*(1.-exp(-(zt/zs1)**2.5))**0.8
	sati2=(zs2/zt)**2*exp(-csi*zt/vu**2.1)
	sati2=sati2*(1.-exp(-(zt/zs2)**2.5))**0.8
	sati3=(zs3/zt)**2*exp(-csi*zt/vu**2.1)
	sati3=sati3*(1.-exp(-(zt/zs3)**2.5))**0.8
cc......................................................
	secs(15)=Dk*sati1
	secs(16)=Dl1*sati2
	secs(17)=Dl2*sati2
	secs(18)=Dm1*sati3
	secs(19)=Dm1*sati3
	secs(20)=Dm2*sati3
	secs(21)=Dn
cc excitation to n=5 to be added to ionization
c     set to 0
	secs(22)=es5*sate1*0.
c	(es5 includes a sum up to infinity)
	secs(23)=e2s5*3.05*sate2*0.
	secs(24)=e2p5*3.05*sate2*0.
c	(theoretical factor in 1/n**3 = 3.049358 for n=5)
	secs(25)=e3s5*2.5*sate3*0.
	secs(26)=e3p5*2.5*sate3*0.
	secs(27)=e3d5*2.5*sate3*0.
c	(3 ->5 => 2.5 in place of 3.05)
	secs(28)=e456*0.
c	(e456=e45+2.5*e46 -> 2.5 in place of 3.5)
cc excitation from 1s
	secs(29)=e2s*sate1
	secs(30)=e2p*sate1
	secs(31)=e3s*sate1
	secs(32)=e3p*sate1
	secs(33)=e3d*sate1
	secs(34)=es4*sate1
cc excitation a partir de n=2
	secs(35)=s3s*sate2
	secs(36)=s3p*sate2
	secs(37)=s3d*sate2
	secs(38)=e2s4*sate2
	secs(39)=p3s*sate2
	secs(40)=p3p*sate2
	secs(41)=p3d*sate2
	secs(42)=e2p4*sate2
cc excitation a partir de n=3 vers n=4
	secs(43)=e3s4*sate3
	secs(44)=e3p4*sate3
	secs(45)=e3d4*sate3
cc
	end if
cc intracouche
	secs(46)=esp
	secs(47)=esp3
	secs(48)=epd3
c......................
	do 2001 I=1,48
	  secs(I)=cor(I)*secs(I)
2001	  continue
	do 2002 I=1,7
c  capture totale (1 a 7)
	  sec(i)=secs(i)+secs(i+7)
c  perte totale (8 a 14)
	  sec(i+7)=secs(i+14)+secs(i+21)
2002	continue
c  excitation (15 a 34)
	do 2003 i=1,20
	  sec(i+14)=secs(i+28)
2003  continue
	return
	end
cc	************** fin de Csec ***************************