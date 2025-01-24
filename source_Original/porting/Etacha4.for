	program etacha
cc*****************************************************************************
cc***     This version of ETACHA is for ions with up to 60 électrons        ***
cc***     First PC version of ETACHA was 18/12/91                           ***
cc***     First distributed PC version of ETACHA was 07/96                  ***
cc*****************************************************************************
cc
cc   05/2015 version going up to nmax=4 with partially correlated states
cc   Last corrections 09/2015 
cc   
cc   NOTES :    1) COA = 1-(1.75*qmin/vu)**2
cc			  2) CDW calculations are not used any more, SEIK instead
cc			  3) CDWEIS +SC & ASC calculations for ionisation 
cc			  4) SE + SC & ASC calculations for excitation
cc              5) Excitation to n >= 5 not added to ionization anymore    // Oleg may be make it optional?
cc                 (See Donaut4.for)
cc
cc   Contact: J.P.ROZET  eMail: rozet@insp.jussieu.fr
cc            D. Vernhet eMail: dominique.vernhet@insp.jussieu.fr
cc            E. Lamour  eMail: lamour@insp.jussieu.fr
cc
cc	   Calculates charge states of swift ions and their evolution as a fonction
cc	of traversed target thickness for a number of initial electrons =< 60.
cc	In a first step, only 1s,2s et 2p states are considered (as in oldest ETACHA),
cc	and 3s,3p and 3d state populations independently estimated.
cc	This program then solves first a set of 84 coupled differential equations : 
cc    63 are for the evolution of the 63 "correlated states of the type Y(i,j,k),
cc	where i,j and k stand for the number of 2p,2s and 1s electrons
cc	(i between 0 et 6, j and k between 0 and 2), and 21 are for the evolution
cc	of 3s (3 equations), 3p (7 equations) and 3d (11 equations) substates.
cc	Differential equations are integrated using a Runge-Kutta type method and associated subroutines
cc    (see EQDIFF.for, an improvement of Adams code - a variable order predictor-corrector method)
cc	
cc	In the actual version, the evolution of n<4 states is calculated in an improved way, 
cc	by considering the  11*19=209 Y(n12,n3) type states, where n12 is the electron number in n=1 and 2,
cc    and n3 the electron number in n=3. For this calculation, actual (averadged) cross sections,
cc    functions of n=n12+n3 are used (see SecMean.for).
cc	This corresponds to a total of 84+209=293 equations (this is the ETACHA3 version).
cc
cc	A final evolution is calculated by considering first in an "indépendant" way the n=4 shell (33 equations),
cc    and then finally  the 29*33 Y(n123,n4)type states, where n123 is the total electron number in n<4 and n4 is
cc    the electron number in n=4.
cc
cc__________________________________________________________________
cc		We then have eventually 293+33+29*33=293+990=1283 equations.
cc__________________________________________________________________
cc    
cc    Files to be included for building the exe file (15 files):
cc          -Etacha4.for (present file)
cc          -Auger4.for
cc          -Donaut4.for
cc          -EQDIF.for
cc          -F4.for
cc          -INTG.for
cc          -Pion.for (PWBA ionisation cross sections)
cc          -SecMean4.for
cc          -SEIK.for (capture cross sections in the SE approximation + REC)
cc          -senlm.for (excitation cross sections in teh SE approximation)
cc          -Sex2.for (PWBA excitation cross sections from n=2 and n=3)
cc          -Sexi.for (PWBA excitation cross sections from n=1)
cc          -Snl.for (PWBA intrashell excitation cross sections)
cc          -Tceis.for (CDW-EIS ionisation cross sections for n=1 and 2)
cc          -Zstop.for (Ziegler's SRIM derived stopping power)
cc    Files to be included as data files (4 files)
cc          -Etadon.dat   ==> etadon.etacha
cc          -SCOEFgas.95, SCOEF.95A, SCOEF.95B
cc    Output files :
cc          -ETA09.prn, ETA1019.prn, (ETA2029.prn, ETA3039.prn, ETA4049.prn and ETA5059.prn if usefull)
cc          -ETAPied.prn, POPMean.prn, seceff.prn
cc
cc
cc	***********************************************************
	dimension Y(1284),YX(1284),WORK(1680000),IWORK(1350),PR(62),
     s	P1s(2),P2s(2),P2p(2),PRF(62)
c     LENW=N*N+17*N+204=1670688
c     LENIW=N+21=1305     
	character nomf*24
	common /etats/ICO1(63),ICO2(294),ICO3(1284),NUM1(733),NUM2(1911),
     s	NUM3(3329)
	common/seceff/sec(34),cor(48),secs(48)
	common/secKLM/C12(209),D12(209),RAD2(209),AKL2(209),PA2(209),
     s AKM3(209),RAD3(209),ALM3(209),C3(209),D3(209),E3(209),DE3(209),
     s PA3(209),PA23(209)
	common/sec1234/P14(957),C13(957),D13(957),C4M(957),D4M(957),
     s	E4M(957),DE4M(957),PA4(957),AKLM4(957),PA123(957)
	common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4
	common/tol/ep0,ep1,erel,erabs
	common/aug/AKLL,AKLM,ALMM,AM4
	common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn
	common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m
	common/popi/n02s,n02p,n03s,n03p,n03d,nl0,nm0
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn
	common/proj/zp
      external F
cc    ***********************************************************
cc	initializations and data input procedures
cc	***********************************************************
	ilgn=64
	iprt=1
	NEQ=1283
	iflag=1
	icor=0
	call ini(ICO1,ICO2,ICO3,NUM1,NUM2,NUM3)
 20	call donaut(QP,Z1,AP,EP,Z2,AC,RHO,EPM,iflag,E1,S1,E2,S2,
     s	istp,iprt,ilgn)
	if (iflag.eq.2) then
	write (6,*) 'iflag=2 !! Input data file problem'
	goto 5000
	endif
	ilp=ilgn-55
	do 21 i=1,38
	  if (cor(i).ne.1.) then
	    icor=1
	  end if
 21	continue
	if (icor.eq.1) then
	ilp=ilp-2
	end if
	pas0=20.*ep0
	zp=z1
	ip=zp
	e=ep
	zc=z2
	VC=(1-(1+E/931.5)**(-2))**(0.5)
	VU=137.036*VC
	ok=omk(zp)
	ol=oml(zp)
cc....... slowing down coefficients ..........
	if(istp.eq.0) then
	  tc=1.e9
	  as=0.
	  bs=e
cc.......!!!
	else
	  as=(s2-s1)/((e2-e1)*ap)
	  bs=(s2*e1-s1*e2)/(s2-s1)
	  stp=s1+(e-E1)*(s2-s1)/(e2-e1)
	  if(stp.lt.0.) then
	  write(6,*)' '
	  write(6,*)' '
	  write(6,*)' stp=',stp
	  STOP' problem with stopping power values, end of calculation'
	  end if
	  stp1=log10(ap*e/stp)
	  astp=10.**aint(stp1)
	  if (stp1.lt.0.) then
	    astp=astp/10.
	  end if
	  sstp=ap*e/(stp*astp)
	  if (sstp.le.1.5) then
	    tc=2.5*astp
	  else if (sstp.le.3.5) then
	    tc=5.*astp
	  else if (sstp.lt.7.) then
	    tc=10.*astp
	  else
	    tc=25.*astp
	  end if
	end if
	if (zp.le.6) then
	  ddq=8.
	  if (erel.gt.1.e-4) then
	    ddq=ddq/2.
	  end if
	else if (zp.le.20) then
	  ddq=4.
	  if (erel.gt.1.e-4) then
	    ddq=ddq/2.
	  end if
	else if (zp.le.40) then
	  ddq=2.
	  if (erel.gt.1.e-4) then
	    ddq=ddq/2.
	  end if
	else
	  ddq=1.
	end if
	write (6,*) 'initialization of populations'
cc	***** number of electrons in shells of neutral atom ***********
	n02s=0
	n02p=0
	n03s=0
	n03p=0
	n03d=0
	n04s=0
	n04p=0
	n04d=0
	n04f=0
	if (IP.le.2) then
	continue
	elseif (IP.le.4) then
	n02s=IP-2
	elseif (IP.le.10) then
	n02s=2
	n02p=IP-4
	elseif (IP.le.12) then
	n02s=2
	n02p=6
	n03s=IP-10
	elseif (IP.le.18) then
	n02s=2
	n02p=6
	n03s=2
	n03p=IP-12
	elseif (IP.le.20) then
	n02s=2
	n02p=6
	n03s=2
	n03p=6
	n04s=IP-18
	elseif (IP.le.28) then
	n02s=2
	n02p=6
	n03s=2
	n03p=6
	n04s=2
	n03d=IP-20
	elseif (IP.le.36) then
	n02s=2
	n02p=6
	n03s=2
	n03p=6
	n04s=2
	n03d=10
	n04p=IP-30
	elseif (IP.le.38) then
	n02s=2
	n02p=6
	n03s=2
	n03p=6
	n04s=2
	n03d=10
	n04p=6
	n05s=IP-36
	elseif (IP.le.48) then
	n02s=2
	n02p=6
	n03s=2
	n03p=6
	n04s=2
	n03d=10
	n04p=6
	n05s=2
	n04d=IP-38
	elseif (IP.le.54) then
	n02s=2
	n02p=6
	n03s=2
	n03p=6
	n04s=2
	n03d=10
	n04p=6
	n05s=2
	n04d=10
	n05p=IP-48
	elseif (IP.le.56) then
	n02s=2
	n02p=6
	n03s=2
	n03p=6
	n04s=2
	n03d=10
	n04p=6
	n05s=2
	n04d=10
	n05p=6
	n06s=IP-54
	else
	n02s=2
	n02p=6
	n03s=2
	n03p=6
	n04s=2
	n03d=10
	n04p=6
	n05s=2
	n04d=10
	n05p=6
	n06s=2
	n04f=IP-56
	endif
	nl0=n02s+n02p
	nm0=n03s+n03p+n03d
	nn0=n04s+n04p+n04d+n04f
cc	***** Y(n) initialization *****************************
	do 1 n=1,1283
	Y(n)=0.
  1	continue
cc	***** number of electrons in shells of incident ion ***********
	NEL=Zp-Qp
	n1s=0
	n2s=0
	n2p=0
	n3s=0
	n3p=0
	n3d=0
	n4s=0
	n4p=0
	n4d=0
	n4f=0
	if (NEL.le.2) then
	n1s=NEL
	elseif (NEL.le.4) then
	n1s=2
	n2s=NEL-2
	elseif (NEL.le.10) then
	n1s=2
	n2s=2
	n2p=NEL-4
	elseif (NEL.le.12) then
	n1s=2
	n2s=2
	n2p=6
	n3s=NEL-10
	elseif (NEL.le.18) then
	n1s=2
	n2s=2
	n2p=6
	n3s=2
	n3p=NEL-12
	elseif (NEL.le.28) then
	n1s=2
	n2s=2
	n2p=6
	n3s=2
	n3p=6
	n3d=NEL-18
	elseif (NEL.le.30) then
	n1s=2
	n2s=2
	n2p=6
	n3s=2
	n3p=6
	n3d=10
	n4s=NEL-28
	elseif (NEL.le.36) then
	n1s=2
	n2s=2
	n2p=6
	n3s=2
	n3p=6
	n3d=10
	n4s=2
	n4p=NEL-30
	elseif (NEL.le.46) then
	n1s=2
	n2s=2
	n2p=6
	n3s=2
	n3p=6
	n3d=10
	n4s=2
	n4p=6
	n4d=NEL-36
	elseif (NEL.le.60) then
	n1s=2
	n2s=2
	n2p=6
	n3s=2
	n3p=6
	n3d=10
	n4s=2
	n4p=6
	n4d=10
	n4f=NEL-46
	else
	write (6,*) ' '
	write (6,*) '**************************************************'
	write (6,*) ' too much electrons ...'
	write (6,*) ' try again with ZP-Q < 61'
	write (6,*) '**************************************************'
	write (6,*) ' '
	go to 20
	endif
	nkl=n1s+n2s+n2p
	nm=n3s+n3p+n3d
	nklm=nkl+nm
	n4l=n4s+n4p+n4d+n4f
	Ic1=n2p*100+n2s*10+n1s
	Ic2=nm*100+nkl
	n1=num(Ic1)
	n2=numP(Ic2)
	write (6,*) 'initialization of populations for actual charge'
	Y(n1)=1.
	Y(64+n3s)=1.
	Y(67+n3p)=1.
	Y(74+n3d)=1.
	Y(n2)=1.
	Y(294+n4l)=1.
	n1234=327+nklm+29*n4l
	Y(n1234)=1.
cc	******************************************************************
cc	*** radiative and Auger "cross sections" (units 10e-20 cm²) ***
	conv=6.023E-3/AC
	RAD=3.461E-6*(Zk*tetak)**4*AC/(RHO*VC)/2.
	RAD3s=.034891E-6*(Zk*tetak)**4*AC/(RHO*VC)/6.
	RAD3p=1.0301E-6*(Zk*tetak)**4*AC/(RHO*VC)/2.
	RAD3d=.3544E-6*(Zk*tetak)**4*AC/(RHO*VC)/6.
	RAD4=0.0454E-6*(Zk*tetak)**4*AC/(RHO*VC)/18.
	if (n02p.ge.1) then
	AKLL=rad*(1./ok-1.)*n02p/(nl0*(nl0-1))
	else
	AKLL=7.6E-2*AC/(RHO*VC)
	endif
	if (nm0.ge.1) then
	AKLM=rad3p*0.882*(1./ok-1.)*n02p/(nl0*nm0)
	else
	AKLM=0.3*AC/(RHO*VC)
	endif
	if (nm0.ge.2) then
	ALMM=(rad3s*n03s+0.118*rad3p*n03p+rad3d*n03d)*(1./ol-1.)
     s /(nm0*(nm0-1))
	else
	ALMM=0.05*AC/(RHO*VC)
	endif
	AM4=3.*ALMM
cc      ***************************************************************
cc      Open and write headings of output files
cc      ***************************************************************
      call datetime(iyr,imon,iday,ihour,imin,isec,imil)
cc ...................................................................
200      format(2x,i2.2,1h/,i2.2,1h/,i4,21x,7hEtacha4,13x,
     s      3h(C),2x,9hINSP-ASUR,2x,11hJPR 05/2015)
202      format(89('*'))
203      format(118('*'))
204      format(2x,"PROJECTILE: atomic number=",f4.0,2x,"incident ch",
     s   "harge=",f4.0,2x,"atomic mass=",f4.0,/,12x,"incident energy=",
     s   f8.3," MeV/u",2x,"velocity=",f8.3," (au)",/,
     s   6x,"TARGET: atomic number=",f4.0,2x,"atomic",
     s   " mass=",f4.0,2x,"density=",f6.3," g/cm3")
205      format(2x,"relative error=",g11.4,3x,"absolute=",g11.4)
206      format(2x,"No stopping power correction",/)
207      format(2x,"Energy loss: S1=",f8.3," MeV/mg/cm² at E1=",f8.3,
     s "MeV/u",/,10x,"and: S2=",f8.3," MeV/mg/cm² at E2=",f8.3," MeV/u")   
208	FORMAT(2x,"Reference ('hydrogenic') cross sections in 10E-20 cm²:",
     s    /,"  CAPTURE (MEC+REC) TO:  ",
     s    " 1s =",G10.3,"  2s =",G10.3,"  2p =",G10.3,/,25x,
     s    " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
     s    "  (n=4) =",G10.3,/,
     s	"  IONIZATION OF:",9x,
     s    " 1s =",G10.3,"  2s =",G10.3,"  2p =",G10.3,/,25x,
     s    " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
     s    "  (n=4) =",G10.3,/,
     s	2x,"1s EXCITATION TO:",
     s    23x,            "2s =",G10.3,"  2p =",G10.3,/,25x,
     s    " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
     s	"  (n=4) =",G10.3,/,
     s	2x,"2s EXCITATION TO:",6x,
     s	" 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
     s	"  (n=4) =",G10.3,/,
     s    2x,"2p EXCITATION TO:",6x,
     s    " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
     s    "  (n=4) =",G10.3,/,
     s	2x,"EXCITATION To n=4 from:",
     s	" 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,/,
     s	2x,"INTRASHELL EXCITATION:",
     s	"  2s to 2p =",G10.3,2x,"3s to 3p =",G10.3)
309   FORMAT (\,"  Binding Energy correction factors =",6(2x,F8.4))
310   FORMAT (\,"  Effective charges =",6(2x,F8.4))    
311   FORMAT ("  Binding energies  =",6(2x,F8.1)) 
cc ....................................................................
210      format('Tar.thick.    0e-       1e-       2e-       3e-',
     s  '       4e-       5e-       6e-       7e-       8e-       9e-',
     s  '     Eout')
211      format(1x,"(ug/cm²) ",10(3x,"(",i2,"+)  "))
215      format(1x,"(ug/cm²) ",10(3x,"(",i2,"+)  "))
212      format('Tar.thick.    10e-      11e-      12e-      13e-',
     s  '      14e-      15e-      16e-      17e-      18e-      19e-')
213      format('Tar.thick.    20e-      21e-      22e-      23e-',
     s  '      24e-      25e-      26e-      27e-      28e-      29e-')
214      format('Tar.thick.    30e-      31e-      32e-      33e-',
     s  '      34e-      35e-      36e-      37e-      38e-      39e-')
218      format('Tar.thick.    40e-      41e-      42e-      43e-',
     s  '      44e-      45e-      46e-      47e-      48e-      49e-')
219      format('Tar.thick.    50e-      51e-      52e-      53e-',
     s  '      54e-      55e-      56e-      57e-      58e-      59e-')
216      format('T (ug/cm²)  bare    1s     2s     2p    1s²',
     s  '    1s2s   1s2p  1s²2s  1s²2p   tot')
217      format('T (ug/cm²)  y1s       y2s      y2p      ym',
     s  '      yn      Qm      Qm in   Qm out    PTOT')
cc ...................................................................
	nomf='results\ETA0009.PRN'
	open (10,file=nomf)
	write (10,200) imon,iday,iyr
	write (10,202)
	write (10,204) zp,qp,ap,e,vu,zc,ac,rho
	write (10,205) erel,erabs
	if (istp.eq.0) then
	    write (10,206)
	else
	    write (10,207) S1,E1,S2,E2
	end if
	write (10,208) (sec(I),I=1,33)
	write (10,203)
	write (10,210)
	write (10,211)ip,ip-1,ip-2,ip-3,ip-4,ip-5,ip-6,ip-7,ip-8,ip-9
	write (10,203)
c
	nomf='results\ETA1019.PRN'
	open (11,file=nomf)
	write (11,200) imon,iday,iyr
	write (11,202)
	write (11,204) zp,qp,ap,e,vu,zc,ac,rho
	write (11,205) erel,erabs
	if (istp.eq.0) then
	    write (11,206)
	else
	  write (11,207) S1,E1,S2,E2
	end if
	write (11,208) (sec(I),I=1,33)
	write (11,203)
	write (11,212)
	if (ip.ge.19) then
	  write (11,211)ip-10,ip-11,ip-12,ip-13,ip-14,ip-15,ip-16,ip-17,ip-18,
     s ip-19
	else
	  write (11,*) ' '
	endif
	write (11,203)
c
	if (ip.ge.20) then
	nomf='results\ETA2029.PRN'
	open (12,file=nomf)
	write (12,200) imon,iday,iyr
	write (12,202)
	write (12,204) zp,qp,ap,e,vu,zc,ac,rho
	write (12,205) erel,erabs
	if (istp.eq.0) then
	    write (12,206)
	else
	    write (12,207) S1,E1,S2,E2
	end if
	write (12,208) (sec(I),I=1,33)
	write (12,203)
	write (12,213)
	if (ip.ge.29) then
	  write (12,211)ip-20,ip-21,ip-22,ip-23,ip-24,ip-25,ip-26,ip-27,ip-28,
     s ip-29
	else
	  write (12,*) ' '
	endif
	write (12,203)
	end if
c
	if (ip.ge.30) then
	nomf='results\ETA3039.PRN'
	open (13,file=nomf)
	write (13,200) imon,iday,iyr
	write (13,202)
	write (13,204) zp,qp,ap,e,vu,zc,ac,rho
	write (13,205) erel,erabs
	if (istp.eq.0) then
	    write (13,206)
	else
	    write (13,207) S1,E1,S2,E2
	end if
	write (13,208) (sec(I),I=1,33)
	write (13,203)
	write (13,214)
	if (ip.ge.39) then
	  write (13,215)ip-30,ip-31,ip-32,ip-33,ip-34,ip-35,ip-36,ip-37,ip-38,
     s ip-39
	else
	  write (13,*) '   '
	endif
	write (13,203)
	end if
c	
	if (ip.ge.40) then
	nomf='results\ETA4049.PRN'
	open (20,file=nomf)
	write (20,200) imon,iday,iyr
	write (20,202)
	write (20,204) zp,qp,ap,e,vu,zc,ac,rho
	write (20,205) erel,erabs
	if (istp.eq.0) then
	    write (20,206)
	else
	    write (20,207) S1,E1,S2,E2
	end if
	write (20,208) (sec(I),I=1,33)
	write (20,203)
	write (20,218)
	if (ip.ge.49) then
	  write (20,215)ip-40,ip-41,ip-42,ip-43,ip-44,ip-45,ip-46,ip-47,ip-48,
     s ip-49
	else
	  write (20,*) '   '
	endif
	write (20,203)
	end if
c	
	if (ip.ge.50) then
	nomf='results\ETA5059.PRN'
	open (21,file=nomf)
	write (21,200) imon,iday,iyr
	write (21,202)
	write (21,204) zp,qp,ap,e,vu,zc,ac,rho
	write (21,205) erel,erabs
	if (istp.eq.0) then
	    write (21,206)
	else
	    write (21,207) S1,E1,S2,E2
	end if
	write (21,208) (sec(I),I=1,33)
	write (21,203)
	write (21,219)
	if (ip.ge.59) then
	  write (21,215)ip-50,ip-51,ip-52,ip-53,ip-54,ip-55,ip-56,ip-57,ip-58,
     s ip-59
	else
	  write (21,*) '   '
	endif
	write (21,203)
	end if
c
	nomf='results\ETAPIED.PRN'
	open (14,file=nomf)
	write (14,200) imon,iday,iyr
	write (14,202)
	write (14,204) zp,qp,ap,e,vu,zc,ac,rho
	write (14,205) erel,erabs
	if (istp.eq.0) then
	  write (14,206)
	else
	  write (14,207) S1,E1,S2,E2
	end if
	write (14,208) (sec(I),I=1,33)
	write (14,203)
	write (14,*) ' '
	write (14,216)
	write (14,203)
c
	nomf='results\POPMEAN.PRN'
	open (15,file=nomf)
	write (15,200) imon,iday,iyr
	write (15,202)
	write (15,204) zp,qp,ap,e,vu,zc,ac,rho
	write (15,205) erel,erabs
	if (istp.eq.0) then
	  write (15,206)
	else
	  write (15,207) S1,E1,S2,E2
	end if
	write (15,208) (sec(I),I=1,33)
	write (15,202)
	write (15,*) '  '
	write (15,217)
	write (15,203)
cc
	nomf='results\seceff.PRN'
	open (16,file=nomf)
	write (16,200) imon,iday,iyr
	write (16,202)
	write (16,204) zp,qp,ap,e,vu,zc,ac,rho
	write (16,205) erel,erabs
	if (istp.eq.0) then
	  write (16,206)
	else
	  write (16,207) S1,E1,S2,E2
	end if
	write (16,208) (sec(I),I=1,33)
	write (16,309) tetak,tetal1,tetal2,tetam1,tetam2,tetan
	write (16,310) zk,zl1,zl2,zm1,zm2,zn
	write (16,311) Bk,Bl1,Bl2,Bm1,Bm2,Bn
	write (16,202)
cc
cc	*******************************************************
cc	calculates cross sections (incident ion) in reduced units (µg/cm²)-1
cc	*******************************************************
	call POPMEAN(Y)
	QM=Zp-(y1s+yl+ym+yn)
	call CSEC(E)
	write (16,*) 'QM=',QM,'E=',E
	write (16,208) (sec(I),I=1,33)
	write (16,309) tetak,tetal1,tetal2,tetam1,tetam2,tetan
	write (16,310) zk,zl1,zl2,zm1,zm2,zn
	write (16,311) Bk,Bl1,Bl2,Bm1,Bm2,Bn
	write (16,203)
220	format(1x,f8.3,1x,f7.3,6(1x,G10.3))
	do 2 I=1,34
	Sec(I)=Sec(I)*conv
2	continue
cc
cc    cross sections per electron and per hole
cc
	Sec(1)=sec(1)/2.
	sec(2)=sec(2)/2.
	sec(3)=sec(3)/6.
	sec(4)=sec(4)/2.
	sec(5)=sec(5)/6.
	sec(6)=sec(6)/10.
	sec(7)=sec(7)/32.
	sec(15)=sec(15)/2.
	sec(16)=sec(16)/6.
	sec(17)=sec(17)/2.
	sec(18)=sec(18)/6.
	sec(19)=sec(19)/10.
	sec(20)=sec(20)/32.
	sec(21)=sec(21)/2.
	sec(22)=sec(22)/6.
	sec(23)=sec(23)/10.
	sec(24)=sec(24)/32.
	sec(25)=sec(25)/2.
	sec(26)=sec(26)/6.
	sec(27)=sec(27)/10.
	sec(28)=sec(28)/32.
	sec(29)=sec(29)/32.
	sec(30)=sec(30)/32.
	sec(31)=sec(31)/32.
	sec(32)=sec(32)/6.
	sec(33)=sec(33)/6.
	sec(34)=sec(34)/10.
cc	
	Rad=rad*conv
	Rad3s=Rad3s*conv
	Rad3p1=Rad3p*0.882*conv
	Rad3p2=Rad3p*0.118*conv
	Rad3d=Rad3d*conv
	Rad4=Rad4*conv
	rad3s=rad3s+sec(25)
	rad3p1=rad3p1+sec(18)
	rad3p2=rad3p2+sec(22)
	rad3d=rad3d+sec(27)
	AKLL=AKLL*conv
	AKLM=AKLM*conv
	ALMM=ALMM*conv
	AM4=AM4*conv
cc	***********************************************
cc	actual cross sections
cc	***********************************************
	T=0.
	Call SecMean(Y,T)
cc	*******************************************************
cc	start integration
cc	*******************************************************
	call datetime(iyr,imon,iday,ih0,imi0,isec0,it0)
	write (6,160) ih0,imi0,isec0
160	format (2x,12hbegin time =,2x,i2,3h h ,i2,4h mn ,i2,2h s)
	t0=(ih0*3600+imi0*60+isec0+it0/1000.)
	iboucle=0
	MINT=3
	NROOT=0
	T=0.
	TOUT=ep0
	toutc=0.
	QM0=QM
100	dtout=(tout-toutc)/tc
	DQM=ABS((QM-QM0)*ddq)
	if((dtout.ge.1.).or.(DQM.ge.1.)) then
	  dt=tout-toutc
	  call chgt(e,zp,ac,rho,as,bs,dt,ef,qm)
	  e=ef
	  call popmean(Y)
	  QM=Zp-(y1s+yl+ym+yn)
	  Call SecMean(Y,T)
	  toutc=tout
	  QM0=QM
	  call datetime(iyr,imon,iday,ih1,imi1,isec1,it1)
	  t1=(ih1*3600+imi1*60+isec1+it1/1000.)
	  ttot=t1-t0
	  iht=ttot/3600
	  imt=(ttot-3600*iht)/60
	  tst=(ttot-3600*iht-60*imt)
	  write (6,105) iht,imt,tst
	end if
105	format (2x23hrestart integration ...,14h(elapsed time=,
     s	2x,i2,3h h ,i2,4h mn ,f5.2,3h s))
	MSTATE=1
	EPS=erel
	EWT=erabs
	LENW=1670000
	LENIW=1310
	X=T
	nequ=0
	YPMax=0.
cc	****************************************************
	call EQDIF(NEQ,T,Y,F,TOUT,MSTATE,NROOT,EPS,EWT,MINT,WORK,
     s	LENW,IWORK,LENIW,F)
cc	*************** IDID=1 as long as T <= TOUT **********
	if (MSTATE.GT.2)THEN
	call message(mstate)
	go to 5000
	end if
	NSTEP=IWORK(3)
	iboucle=NSTEP
	AVGORD=WORK(3)
	do 101 IN=1,1283
	if (Y(IN).lt.0.) then
	Y(IN)=0.
	end if
	YX(IN)=Y(IN)*100.
 101	continue
cc	**************** end of one step calculation ************
	X=T
cc	***********************************************
cc	charge state probabilities calculation
cc	***********************************************
	call popmean(Y)
	QM=Zp-(y1s+yl+ym+yn)
c     charge states before autoionization
      Do 500 N=1,61
	PRF(N)=0.
500	Continue
	Do 600 N1=1,29
	Do 600 N4=1,33
	I=N1-1
	N=N4-1
	INM=100*N+I
	NN=NUMPP(INM)
	Nel=I+N4
	PRF(Nel)=PRF(Nel)+Y(NN)
600	Continue
      Qin=0.
      PTF=0.
	do 700 M=1,61
	RM=Real(M)-1.
	PTF=PTF+PRF(M)
	Qin=Qin+(Zp-RM)*PRF(M)
700	continue
cc	***********************************************
cc	re-calculate cross sections for actual charge state
cc	***********************************************
	Call SecMean(Y,T)
cc-------- K shell autoionization and K+L  populations --------
	Call Auger(Y,Zp,PR,YTOT)
	PT=0.
	QF=0.
	do 9 M=1,62
	RM=Real(M)-1.
	PR(M)=100.*PR(M)
	PT=PT+PR(M)
	QF=QF+(Zp-RM)*PR(M)/100.
9	continue
	do 50 K=1,2
	P1s(K)=0
 50	continue
	do 51 K=1,2
	do 51 J=0,2
	do 51 I=0,6
	L=100*I+10*J+K
	N=NUM(L)
	P1s(K)=P1s(K)+Y(N)
 51	continue
	do 52 J=1,2
	P2s(J)=0
 52	continue
	do 53 J=1,2
	do 53 K=0,2
	do 53 I=0,6
	L=100*I+10*J+K
	N=NUM(L)
	P2s(J)=P2s(J)+Y(N)
 53	continue
	do 54 I=1,2
	P2p(I)=0
 54	continue
	do 55 I=1,2
	do 55 K=0,2
	do 55 J=0,2
	L=100*I+10*J+K
	N=NUM(L)
	P2p(I)=P2p(I)+Y(N)
 55	continue
cc	********* prints probabilities at each step *****
c	write (10,125) T,(PR(I),I=1,10),e,PT
c 125	format(F10.3,F9.5,9(1x,F9.5),2x,2(1x,g10.5))
	write (10,125) T,(PR(I),I=1,10),e
	write (11,25) T,(PR(I),I=11,20)
	if (ip.ge.20) then
	write (12,25) T,(PR(I),I=21,30)
	end if
	if (ip.ge.30) then
	write (13,26) T,(PR(I),I=31,40)
	end if
	if (ip.ge.40) then
	write (20,26) T,(PR(I),I=41,50)
	end if
	if (ip.ge.50) then
	write (21,26) T,(PR(I),I=51,60)
	end if
	tots=yx(1)+yx(2)+yx(4)+yx(10)+yx(3)+yx(5)+yx(11)+yx(6)+yx(12)
	write (14,27) T,yx(1),yx(2),yx(4),yx(10),yx(3),yx(5),yx(11),
     s	yx(6),yx(12),tots
c	
	write (15,28) T,y1s,y2s,y2p,ym,yn,Qm,Qin,QF,PTF	
25	format(F10.3,F9.5,9(1x,F9.5))
125	format(F10.3,F9.5,9(1x,F9.5),3x,g12.5)
26	format(F9.2,F9.5,9(1x,F9.5))
27	format(F9.2,10(1x,F6.2))
28	format(F9.2,3(1x,F7.5),6(1x,F8.5))
	T=TOUT
	write(6,29) T,iboucle,nequ
29	format(2x,19hachieved thickness=,f9.2,4x,17h(iteration # for ,
     s	10hthis step=,i5,1x,i5,1h))
	iboucle=0
cc********* calculation goes on as long as T <= EPM *****
	if (TOUT.lt.EPM) then
		if (TOUT.lt.pas0)then
cc*************** (pas0=20*ep0) ******************
		  TOUT=TOUT+ep0
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.le.0.2) then
		  TOUT=TOUT+ep0
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.1.4999) then
		  if (ep1.gt.0.05) then
		   TOUT=TOUT+0.05
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.1.9999) then
		  if (ep1.gt.0.1) then
		   TOUT=TOUT+0.1
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.4.9999) then
		  if (ep1.gt.0.2) then
		   TOUT=TOUT+0.2
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.14.999) then
		  if (ep1.gt.0.5) then
		   TOUT=TOUT+0.5
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.19.99) then
		  if (ep1.gt.1.) then
		   TOUT=TOUT+1.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.49.99) then
		  if (ep1.gt.2.) then
		   TOUT=TOUT+2.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.149.99) then
		  if (ep1.gt.5.) then
		   TOUT=TOUT+5.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.199.99) then
		  if (ep1.gt.10.) then
		   TOUT=TOUT+10.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.499.99) then
		  if (ep1.gt.20.) then
		   TOUT=TOUT+20.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.1499.9) then
		  if (ep1.gt.50.) then
		   TOUT=TOUT+50.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.1999.) then
		  if (ep1.gt.100.) then
		   TOUT=TOUT+100.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.4999.) then
		  if (ep1.gt.200.) then
		   TOUT=TOUT+200.
		  else
		   TOUT=TOUT+ep1
		  end if
		   if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.14999.) then
		  if (ep1.gt.500.) then
		   TOUT=TOUT+500.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.19999.) then
		  if (ep1.gt.1000.) then
		   TOUT=TOUT+1000.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else if (TOUT.lt.49999.) then
		  if (ep1.gt.2000.) then
		   TOUT=TOUT+2000.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		else
		  if (ep1.gt.5000.) then
		   TOUT=TOUT+5000.
		  else
		   TOUT=TOUT+ep1
		  end if
		  if (TOUT.gt.EPM) TOUT=EPM
		  iboucle=0.
		  goto 100
		endif
	endif
cc ************* end of output files **************************
      call datetime(iyr,imon,iday,ih1,imi1,isec1,it1)
	write (6,161) ih1,imi1,isec1
	t1=(ih1*3600+imi1*60+isec1+it1/1000.)
	ttot=t1-t0
	iht=ttot/3600
	imt=(ttot-3600*iht)/60
	tst=(ttot-3600*iht-60*imt)
	write (6,162) iht,imt,tst
161	format ("     end time = ",i2," h ",i2," mn ",i2," s")
162	format (" elapsed time = ",i2," h ",i2," mn ",f5.2," s")
	if(tc.lt.epm) then
	  dt=tout-toutc
	  if (dt.gt.0.) then
	    ef=e+(bs-e)*(1.-exp(-0.001*tc*as))
	    e=ef
	  end if
	end if
	write (6,163) e
163	format (35x,"final energy : ",g11.4,"(MeV/u)")
	write (10,203)
	write (10,210)
	if (ip.ge.9) then
	write (10,211)ip,ip-1,ip-2,ip-3,ip-4,ip-5,ip-6,ip-7,ip-8,ip-9
	end if
	write (10,203)
	write (10,160) ih0,imi0,isec0
	write (10,161) ih1,imi1,isec1
	write (10,163) e	
	close(unit=10)
c	
	write (11,203)
	write (11,212)
	if (ip.ge.19) then
	  write (11,211)ip-10,ip-11,ip-12,ip-13,ip-14,ip-15,ip-16,ip-17,ip-18,
     s ip-19
	endif
	write (11,203)
	write (11,160) ih0,imi0,isec0
	write (11,161) ih1,imi1,isec1
	write (11,163) e
	close(unit=11)
c
	if (ip.ge.20) then
	 write (12,203)
	 write (12,213)
	 if (ip.ge.29) then
	 write (12,211)ip-20,ip-21,ip-22,ip-23,ip-24,ip-25,ip-26,ip-27,ip-28,
     s ip-29
       endif
	  write (12,203)
	  write (12,160) ih0,imi0,isec0
	  write (12,161) ih1,imi1,isec1
	  write (12,163) e
	  close(unit=12)
	end if
c	
	if (ip.ge.30) then
	  write (13,203)
	  write (13,214)
	  if (ip.ge.39) then
	  write (13,215)ip-30,ip-31,ip-32,ip-33,ip-34,ip-35,ip-36,ip-37,ip-38,
     s  ip-39
	  endif
	  write (13,203)
	  write (13,160) ih0,imi0,isec0
	  write (13,161) ih1,imi1,isec1
	  write (13,163) e
	  close(unit=13)
	end if
c	
	if (ip.ge.40) then
	  write (20,203)
	  write (20,218)
	  if (ip.ge.49) then
	  write (20,215)ip-40,ip-41,ip-42,ip-43,ip-44,ip-45,ip-46,ip-47,ip-48,
     s  ip-49
	  endif
	  write (20,203)
	  write (20,160) ih0,imi0,isec0
	  write (20,161) ih1,imi1,isec1
	  write (20,163) e
	  close(unit=20)
	end if
c
      if (ip.ge.50) then
	  write (21,203)
	  write (21,219)
	  if (ip.ge.59) then
	  write (21,215)ip-50,ip-51,ip-52,ip-53,ip-54,ip-55,ip-56,ip-57,ip-58,
     s ip-59
	  endif
	  write (21,203)
	  write (21,160) ih0,imi0,isec0
	  write (21,161) ih1,imi1,isec1
	  write (21,163) e
	  close(unit=21)
	end if
c	
	write (14,203)
	write (14,216)
	write (14,203)
	write (14,160) ih0,imi0,isec0
	write (14,161) ih1,imi1,isec1
	write (14,163) e
c	
	do n3=327,1284
	if (Y(n3).ge.0.01) then
	write(14,*) 'Code=',Ico3(n3), 'pop=',Y(n3)
	end if
	end do
	close(unit=14)
c	
	write (15,203)
	write (15,217)
	write (15,203)
	write (15,160) ih0,imi0,isec0
	write (15,161) ih1,imi1,isec1
	write (15,163) e
	close(unit=15)
c
	close(unit=16)
	write(6,*) 'output data in files:'
	write(6,*) '   0 to  9 e- charge states    in ETA09.PRN'
	write(6,*) '  10 to 19 e- charge states    in ETA1019.PRN'
	if (ip.ge.20) then
	write(6,*) '  20 to 29 e- charge states    in ETA2029.PRN'
	end if
	if (ip.ge.30) then
	write(6,*) '  30 to 39 e- charge states    in ETA3039.PRN'
	end if
	if (ip.ge.40) then
	write(6,*) '  40 to 49 e- charge states    in ETA4049.PRN'
	end if
	if (ip.ge.50) then
	write(6,*) '  50 to 59 e- charge states    in ETA5059.PRN'
	end if
	write(6,*) ' bare,1s,2s,2p,1sý,1s2s,1s2p,1sý2s,1sý2p ions '
	write(6,*) ' and sum of these             in ETAPIED.PRN'
	write(6,*) ' mean 1s,2s,2p,3s,3p and 3d populations'
	write(6,*) '                              in POPMEAN.PRN'
	write(6,*) '  '
	write(6,*)'WARNING! Next calculation will overwrite these files'
	write(6,*)'Consider saving or renaming these results ! '
5000	continue
	END
cc    ***************************************************
cc	*********** END of etacha *************************
cc    ***************************************************
	subroutine chgt(e,zp,ac,rho,as,bs,tc,ef,qm)
	common/seceff/sec(34),cor(48),secs(48)
	common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4
	common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,
     s tetan
      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn
	common/aug/AKLL,AKLM,ALMM,AM4
	common/popi/n02s,n02p,n03s,n03p,n03d,nl0,nm0
	ef=e+(bs-e)*(1.-exp(-0.001*tc*as))
	e=ef
	call CSEC(ef)
220	format(1x,f8.3,1x,f7.3,6(1x,G10.3))
203	format(89('*'))
308	FORMAT(2x,"Reference ('hydrogenic') cross sections in 10E-20 cm²:",
     s    /,"  CAPTURE (MEC+REC) TO:  ",
     s    " 1s =",G10.3,"  2s =",G10.3,"  2p =",G10.3,/,25x,
     s    " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
     s    "  (n=4) =",G10.3,/,
     s	"  IONIZATION OF:",9x,
     s    " 1s =",G10.3,"  2s =",G10.3,"  2p =",G10.3,/,25x,
     s    " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
     s    "  (n=4) =",G10.3,/,
     s	2x,"1s EXCITATION TO:",
     s    23x,            "2s =",G10.3,"  2p =",G10.3,/,25x,
     s    " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
     s	"  (n=4) =",G10.3,/,
     s	2x,"2s EXCITATION TO:",6x,
     s	" 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
     s	"  (n=4) =",G10.3,/,
     s    2x,"2p EXCITATION TO:",6x,
     s    " 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,
     s    "  (n=4) =",G10.3,/,
     s	2x,"EXCITATION To n=4 from:",
     s	" 3s =",G10.3,"  3p =",G10.3,"  3d =",G10.3,/,
     s	2x,"INTRASHELL EXCITATION:",
     s	"  2s to 2p =",G10.3,2x,"3s to 3p =",G10.3)
309   FORMAT (\,"  Binding Energy correction factors =",6(2x,F8.4))
310   FORMAT (\,"  Effective charges =",6(2x,F8.4))    
311   FORMAT ("  Binding energies  =",6(2x,F8.1)) 
	write (16,*) 'QM=',QM,'E=',E
	write (16,308) (sec(I),I=1,33)
	write (16,309) tetak,tetal1,tetal2,tetam1,tetam2,tetan
	write (16,310) zk,zl1,zl2,zm1,zm2,zn
	write (16,311) Bk,Bl1,Bl2,Bm1,Bm2,Bn
	write (16,203)
	dum=qm
	e=ef
	conv=6.023E-3/AC
	write(6,*)'Radiative and Auger yields'
	VC=(1-(1+E/931.5)**(-2))**(0.5)
	RAD=3.461E-6*(Zk*tetak)**4*AC/(RHO*VC)/2.
	RAD3s=.034891E-6*(Zk*tetak)**4*AC/(RHO*VC)/6.
	RAD3p=1.0301E-6*(Zk*tetak)**4*AC/(RHO*VC)/2.
	RAD3d=.3544E-6*(Zk*tetak)**4*AC/(RHO*VC)/6.
	RAD4=0.0454E-6*(Zk*tetak)**4*AC/(RHO*VC)/18.
	if (n02p.ge.1) then
	AKLL=rad*(1./ok-1.)*n02p/(nl0*(nl0-1))
	else
	AKLL=7.6E-2*AC/(RHO*VC)
	endif
	if (nm0.ge.1) then
	AKLM=rad3p*0.882*(1./ok-1.)*n02p/(nl0*nm0)
	else
	AKLM=0.3*AC/(RHO*VC)
	endif
	if (nm0.ge.2) then
	ALMM=(rad3s*n03s+0.118*rad3p*n03p+rad3d*n03d)*(1./ol-1.)
     s /(nm0*(nm0-1))
	else
	ALMM=0.05*AC/(RHO*VC)
	endif
	AM4=3.*ALMM
	do 1 I=1,34
	Sec(I)=Sec(I)*conv
1	continue
	Sec(1)=sec(1)/2.
	sec(2)=sec(2)/2.
	sec(3)=sec(3)/6.
	Sec(4)=sec(4)/2.
	sec(5)=sec(5)/6.
	sec(6)=sec(6)/10.
	sec(7)=sec(7)/32.
	sec(15)=sec(15)/2.
	sec(16)=sec(16)/6.
	sec(17)=sec(17)/2.
	sec(18)=sec(18)/6.
	sec(19)=sec(19)/10.
	sec(20)=sec(20)/32.
	sec(21)=sec(21)/2.
	sec(22)=sec(22)/6.
	sec(23)=sec(23)/10.
	sec(24)=sec(24)/32.
	sec(25)=sec(25)/2.
	sec(26)=sec(26)/6.
	sec(27)=sec(27)/10.
	sec(28)=sec(28)/32.
	sec(29)=sec(29)/32.
	sec(30)=sec(30)/32.
	sec(31)=sec(31)/32.
	sec(32)=sec(32)/6.
	sec(33)=sec(33)/6.
	sec(34)=sec(34)/10.
	Rad=rad*conv
	Rad3s=Rad3s*conv
	Rad3p1=Rad3p*0.882*conv
	Rad3p2=Rad3p*0.118*conv
	Rad3d=Rad3d*conv
	Rad4=Rad4*conv
	rad3s=rad3s+sec(25)
	rad3p1=rad3p1+sec(18)
	rad3p2=rad3p2+sec(22)
	rad3d=rad3d+sec(27)
	AKLL=AKLL*conv
	AKLM=AKLM*conv
	ALMM=ALMM*conv
	AM4=AM4*conv
	return
	end
cc	***************************************************
	subroutine popmean(u)
	dimension u(1283)
	common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn
	common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m
cc	mean number of electrons in each subshell
	y3s=U(65)+2.*U(66)
	y3p=U(68)+2.*U(69)+3.*U(70)+4.*U(71)+5.*U(72)+6.*U(73)
	y3d=U(75)+2.*U(76)+3.*U(77)+4.*U(78)+5.*U(79)+6.*U(80)
	y3d=y3d+7.*U(81)+8.*U(82)+9.*U(83)+10.*U(84)
	ym=y3s+y3p+y3d
	yn=0.
	ynn=0.
	do 10 n=1,32
	rn=real(n)
	yn=yn+rn*U(294+n)
	if (n.ge.2) then
	ynn=ynn+rn*(rn-1.)*U(294+n)
	end if
10	continue
c -------- see F.for (ALMM) ------------------------------
	ymp=0.
	ymm=0.
	do 1 I=0,2
	do 1 J=0,6
	do 1 K=0,10
	N=I+J+K
	ymp=ymp+U(64+I)*U(67+J)*U(74+K)*N
	if (N.ge.2) then
	ymm=ymm+U(64+I)*U(67+J)*U(74+K)*N*(N-1)
	endif
1	continue
	ym1m=0.
	do 2 I=0,2
	do 2 J=0,6
	N=I+J-1
	if (N.ge.2) then
	ym1m=ym1m+U(64+I)*U(67+J)*(N-1)
	endif
2	continue
	ym2m=U(76)+2.*U(77)+3.*U(78)+4.*U(79)+5.*U(80)
	ym2m=ym2m+6.*U(81)+7.*U(82)+8.*U(83)+9.*U(84)
c------------------------------------------------------
c.....mean values for 1s, 2s, 2p ......
	y1s=0.
	y2s=0.
	y2p=0.
	do 3 N=1,63
	I=II(N)
	J=JJ(I,N)
	K=KK(I,J,N)
	y1s=y1s+K*U(N)
	y2s=y2s+J*U(N)
	y2p=y2p+I*U(N)
3	continue
	yl=y2s+y2p
c.....mean values for 1s², 2s², 2p² (in case is needed) ......
	ykm=0.
	yl1m=0.
	yl2m=0.
	do 4 N=1,63
	I=II(N)
	J=JJ(I,N)
	K=KK(I,J,N)
	if(K.ge.2) then
	ykm=ykm+U(N)
	end if
	if(J.ge.2) then
	yl1m=yl1m+U(N)
	end if
	if(I.ge.2) then
	yl2m=yl2m+(I-1)*U(N)
	end if
4	continue
	return
	end
cc*********************************************************************
	function omk(z)
	dimension ok(92)
	data (ok(i),i=1,92) / 0.,0.,0.001,0.0012,
     s	.170E-02,.280E-02,.520E-02,.830E-02,.130E-01,.180E-01,.230E-01,
     s	.300E-01,.390E-01,.500E-01,.630E-01,.780E-01,.970E-01,.118,
     s	.140	 ,.163	  ,.188    ,.214    ,.243    ,.275    ,.308,
     s	.340	 ,.373	  ,.406    ,.440    ,.474    ,.507    ,.535,
     s	.562	 ,.589	  ,.618    ,.643    ,.667    ,.690    ,.710,
     s	.730	 ,.747	  ,.765    ,.780    ,.794    ,.808    ,.820,
     s	.831	 ,.843	  ,.853    ,.862    ,.870    ,.877    ,.884,
     s	.891	 ,.897	  ,.902    ,.907    ,.912    ,.917    ,.921,
     s	.925	 ,.929	  ,.932    ,.935    ,.938    ,.941    ,.944,
     s	.947	 ,.949	  ,.951    ,.953    ,.955    ,.957    ,.958,
     s	.959	 ,.961	  ,.962    ,.963    ,.964    ,.965    ,.966,
     s	.967	 ,.968	  ,.968    ,.969    ,.969    ,.970    ,.970,
     s	.971	 ,.971	  ,.972	   ,.972 /
	i=Z
	omk=ok(i)
	return
	end
cc*********************************************************************
	function oml(z)
	dimension ol(92)
	data (ol(i),i=1,92) / 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     s	.120E-02,.750E-03,.380E-03,.310E-03,.260E-03,.240E-03,.220E-03,
     s  .270E-03,.330E-03,.840E-03,.150E-02,.260E-02,.370E-02,.500E-02,
     s  .630E-02,.770E-02,.930E-02,.110E-01,.120E-01,.130E-01,.150E-01,
     s  .160E-01,.180E-01,.200E-01,.220E-01,.240E-01,.260E-01,.280E-01,
     s  .310E-01,.340E-01,.370E-01,.400E-01,.430E-01,.460E-01,.490E-01,
     s  .520E-01,.560E-01,.600E-01,.640E-01,.690E-01,.740E-01,.790E-01,
     s  .850E-01,.910E-01,.970E-01,.104    ,.111    ,.118    ,.125    ,
     s	.132	,.139	 ,.147	  ,.155    ,.164    ,.174    ,.182    ,
     s  .192    ,.201    ,.210    ,.220    ,.231    ,.243    ,.255    ,
     s  .268    ,.281    ,.294    ,.306    ,.320    ,.333    ,.347    ,
     s  .360    ,.373    ,.386    ,.399    ,.411    ,.424    ,.437    ,
     s	.450	,.463	 ,.476	  ,.489 /
	i=Z
	oml=ol(i)
	return
	end
cc*********************************************************************
	subroutine message(idid)
	write (6,*)'problem integrating with QDIFF'
	if (idid.eq.3) then
	write (6,*)' more than 1000 iterations have been atempted ...'
	write (6,*)' try reducing maximum step size'
	end if
	if (idid.eq.4) then
	write (6,*)' you are probably asking too much accuracy ...'
	write (6,*)' try again with larger uncertainties'
	end if
	if (idid.GT.4) then
	write (6,*)' for some reason, the problem is very stiff and',
     s	' cannot be solved with the present integration routine'
	end if
	return
	end
cc*********************************************************************
      subroutine datetime(iyr,imon,iday,ihour,imin,isec,imil)
      integer Datime(8)
      Character*12 Dum(3)
      Call date_and_time(Dum(1),Dum(2),Dum(3),Datime)
      iyr=datime(1)
      imon=datime(2)
      iday=datime(3)
      ihour=datime(5)
      imin=datime(6)
      isec=datime(7)
      imil=datime(8)
      return
      end