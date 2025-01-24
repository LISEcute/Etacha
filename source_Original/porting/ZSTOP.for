	subroutine zstop(z1,ea1,z2,stotal)
C
C	THIS PROGRAM IS A FORTRAN VERSION OF ZIEGLER'S BASIC PROGRAM STOP
C     ea1 is energy in MeV/u
C     stotal is return value of total stopping power in MeV/mg/cm2
C
C
	DIMENSION SCOEF(93,54),SCOEFG(24,11)
	CHARACTER*3 SOLGAS
	COMMON/A/ SCOEF,SCOEFG
	COMMON/B/ ATDENS
c
	iz1=int(z1+.5)
	iz2=int(z2+.5)
C
	CALL STOPCO
C
	IF(IZ1.LT.1.OR.IZ1.GT.92)THEN
	STOP ' ILLEGAL VALUE OF Z1 - PROGRAM TERMINATED.'
	END IF
	IF(IZ2.LT.1.OR.IZ2.GT.92)THEN
	STOP ' ILLEGAL VALUE OF Z2 - PROGRAM TERMINATED.'
	END IF
	IF(EA1.LE.0)THEN
	STOP ' ILLEGAL VALUE OF ION ENERGY - PROGRAM TERMINATED.'
	END IF
C
	ISOGAS=0
	IGASNO=0
	DO 135 I=1,8
	IF(IZ2.EQ.SCOEFG(I,1))THEN
	ISOGAS=1
	IGASNO=I
	END IF
135	CONTINUE
	SOLGAS='SOL'
	IF(ISOGAS.EQ.1)SOLGAS='GAS'
C
	RM1=SCOEF(IZ1,3)
	RM2=SCOEF(IZ2,4)
	IF(ISOGAS.EQ.0)THEN
	RHO=SCOEF(IZ2,5)
	ATDENS=SCOEF(IZ2,6)
	ELSE
	RHO=SCOEFG(IGASNO,2)
	ATDENS=SCOEFG(IGASNO,3)
	END IF
C
	IF(ISOGAS.EQ.0)THEN
	ATDENS=ATDENS*RHO/SCOEF(IZ2,5)*SCOEF(IZ2,4)/RM2
	ELSE
	ATDENS=ATDENS*RHO/SCOEF(IGASNO,2)*SCOEF(IZ2,4)/RM2
	END IF
	ENERGY=EA1*RM1*1000.
C
	CORR=.06022/RM2
C
	E0KEV=ENERGY
	SEE=E0KEV
C
	CALL STOP96(IZ1,RM1,IZ2,SEE,ISOGAS)
C
C	Convert (eV-A2) to Correct Units and put into SE
	SE=SEE*CORR
C
C	Get proton stopping for Effective Charge
C
	IZ=1
	RM=1
	SEE=E0KEV/RM1
	CALL STOP96(IZ,RM,IZ2,SEE,ISOGAS)
	IF(SEE.LT.1.E-8)THEN
	WRITE(*,*)' Z1 = ',IZ1,' Z2 = ',IZ2,' SEE = ',SEE
	STOP
     &	' SEE (SUBROUTINE STOP96) LESS THEN 1.E-8 - PROGRAM TERMINATED'
	END IF
C
	EFFCHG=(SE/(SEE*CORR))**0.5
C
C	Calculate ZBL universal nuclear stopping powers.
C	Epsilon is the reduced energy of the ion/target combination.
C
	E=E0KEV
	EPSIL=32.53*RM2*E/(IZ1*IZ2*(RM1+RM2)*
     &	(FLOAT(IZ1)**.23+FLOAT(IZ2)**.23))
	IF(EPSIL.GE.30)THEN
	SN=LOG(EPSIL)/(2*EPSIL)
	ELSE
	A=(.01321*EPSIL**.21226)+(.19593*EPSIL**.5)
	SN=.5*LOG(1+1.1383*EPSIL)/(EPSIL+A)
	END IF
C
C	convert from LSS reduced units to eV-cm2/1E15
	SN=SN*IZ1*IZ2*RM1*8.462/((RM1+RM2)*
     &	(FLOAT(IZ1)**.23+FLOAT(IZ2)**.23))
	SN=SN*CORR*10
C	(Convert (eV-cm2) -> (eV-A2) -> Correct Units)
C
	STOTAL=SE+SN
C
	return
	END
C
C
C
	SUBROUTINE STOP96(IZ1,RM1,IZ2,SE,ISOGAS)
C
C	Output is in variable SE in units of eV/A2.
C	NUCLEAR STOPPING IS NOT CALCULATED.
C	If ISOGAS=0, then stopping in solids, ISOGAS=1 for gases.
C
C	ษอออออออออออออออออออออออออออออออออออออออออออออออออออออออออป
C	ฬออออออออออออออออออป array: SCOEF.95 ษออออออออออออออออออออน
C	ฬออ> Columns 1-8 : ศอออออออออออออออออผ                    บ
C	บ    (1=Atomic Number),(2=Mass of MAI),(3=Weight of MAI)  บ
C	บ    (4=Natural Mass),(5=Target density: g/cm3),          บ
C	บ    (6=Target Density: atoms/cm3),(7=Fermi Velocity)     บ
C	บ    (8=Sublimation Energy)                               บ
C	ฬออ> Columns  9-16 : Proton stopping coefficients.        บ
C	ฬออ> Columns 17-20 : Proton stopping > 10 MeV.            บ
C	ฬออ> Columns 21-39 : Ion Lambda values (width of ion)     บ
C	ฬออ> Columns 40-54 : Tgt Fermi Velocity Correction (<5%)  บ
C	บ    ** MAI = Most Abundant Isotope                       บ
C	ฬอออออออออออออออออออออออออออออออออออออออออออออออออออออออออน
C	บ      Variable definitions for STOP95                    บ
C	ฬออ>   IZ1 = ION ATOMIC NUMBER                            บ
C	ฬออ>   MM1= ION ATOMIC MASS                               บ
C	ฬออ>   RM1 = ION ATOMIC WEIGHT (AMU)                      บ
C	ฬออ>   IZ2 = TARGET ATOMIC NUMBER                         บ
C	ฬออ>   RM2 = TARGET ATOMIC WEIGHT (AMU)                   บ
C	ฬออ>   RHODENSITY = TARGET DENSITY (G/CM3)                บ
C	ฬออ>   ATDENS=TARGET DENSITY (ATOMS/CM3)                  บ
C	ฬออ>   VFERMI = (FERMI VELOCITY OF SOLID) / V0            บ
C	ฬออ>   SE = CALCULATED ELECTRONIC STOPPING (EV-A2)        บ
C	ศอออออออออออออออออออออออออออออออออออออออออออออออออออออออออผ
C
	DIMENSION SCOEF(93,54),SCOEFG(24,11),IGAS(8)
	COMMON/A/ SCOEF,SCOEFG
	COMMON/B/ ATDENS


C	IF(RM1.EQ.0)RM1=SCOEF(IZ1,3)
C	(Weight (amu) of most abundant isotope.)
C	ATDENS=SCOEF(IZ2,6)
C	(Atomic density of target: atoms/cm3)
	VFERMI=SCOEF(IZ2,7)
C	(VFermi of Target)

	IF(ISOGAS.EQ.1)THEN
C	(Determine if GAS request is valid)
	IGAS(1)=1
	IGAS(2)=2
	IGAS(3)=7
	IGAS(4)=8
	IGAS(5)=10
	IGAS(6)=18
	IGAS(7)=36
	IGAS(8)=54
	IGASNO=0
	DO 110 I=1,8
	IF(IZ2.EQ.IGAS(I))THEN
	IGASNO=I
	GO TO 120
	END IF
110	CONTINUE
120	IF(IGASNO.EQ.0)ISOGAS=0
C	(Reset for Solid target.)
	END IF
C
	EKEV=SE/RM1
C	(Energy passed in array SE.)
C
	E=EKEV
C
C	(Branch for Stopping of : IZ1=1, IZ1=2 OR IZ1>2)
C
        IF(IZ1.EQ.1)THEN
C	(Proton Stopping Powers)
	IF(ISOGAS.EQ.1)THEN
	CALL PSTOPG(IZ1,IZ2,E,SE,IGASNO)
	ELSE
	CALL PSTOP(IZ2,E,SE)
	END IF
	END IF
C
	IF(IZ1.EQ.2)THEN
C	(HELIUM Stopping Powers)
C	(VELOCITY PROPORTIONAL STOPPING BELOW KEV/AMU HE0.)
	HE0=1
C	(Units = keV/amu)
	IF(HE0.GT.E)THEN
	HE=HE0
	ELSE
	HE=E
	END IF
	B=LOG(HE)
	A=.2865+B*(.1266+B*(-.001429+B*(.02402+B*(-.01135+B*.001475))))
	IF(A.GT.30)A=30
	HEH=1.-EXP(-A)
C	(ADD IZ1**3 EFFECT TO HE/H STOPPING POWER RATIO HEH.)
	IF(HE.LT.1)HE=1
	A=(1.+(.007+.00005*IZ2)*EXP(-(7.6-LOG(HE))**2))
	HEH=HEH*A*A
	IF(ISOGAS.EQ.1)THEN
	CALL PSTOPG(IZ1,IZ2,HE,SP,IGASNO)
	ELSE
	CALL PSTOP(IZ2,HE,SP)
	END IF
	SE=SP*HEH*4.
	IF(E.GT.HE0)GO TO 310
C	(CALC. HE VELOCITY PROPORTIONAL STOPPING)
	SE=SE*SQRT(E/HE0)
	END IF
C
C
        IF(IZ1.GT.2)THEN
C	(HEAVY ION ELECTRONIC STOPPING POWERS)
C	(USE VELOCITY STOPPING FOR (YRMIN=VR/FLOAT(IZ1)**.67) <= 0.13)
C	(OR FOR VR <= 1.0)
	YRMIN=0.13
	VRMIN=1.0
	V=SQRT(E/25)/VFERMI
C	(Relative Velocity)
	IF(V.LT.1)THEN
	VR=(3*VFERMI/4)*(1+V*V*(2./3-V*V/15.))
	ELSE
	VR=V*VFERMI*(1+1./5/V/V)
	END IF
C	(SET YR=MAXIMUM OF (VR/FLOAT(IZ1)**.67),(VRMIN/FLOAT(IZ1)**.67) OR YRMIN.)
	YR=VR/FLOAT(IZ1)**.6667
	IF(YR.LT.YRMIN)YR=YRMIN
	A=VRMIN/FLOAT(IZ1)**.6667
	IF(YR.LT.A)YR=A
	A=-.803*YR**0.3+1.3167*YR**0.6+.38157*YR+.008983*YR*YR
	IF(A.GT.50)A=50
C	(Prevents Underflow)
	Q=1-EXP(-A)
	IF(Q.LT.0)Q=0
	IF(Q.GT.1)Q=1
C	(Q = IONIZATION LEVEL OF THE ION AT VELOCITY YR.)
C	(NOW WE CONVERT IONIZATION LEVEL TO EFFECTIVE CHARGE.)
C
C	(Screening Distance of Ion (Lambda in B.& K.))
C	(Lambda is in SCOEF(IZ2,22-39))
C	(Interpolation values in SCOEF(93,22-39))
C
	DO 130 J=22,39
C	(Find Q Interpolation)
	IF(Q.LE.SCOEF(93,J))GO TO 140
C	(in SCOEF(93,22-39))
130	CONTINUE
140	J=J-1
	IF(J.LT.22)J=22
	IF(J.GT.38)J=38
	RLAMB0=SCOEF(IZ1,J)
	RLAMB1=(Q-SCOEF(93,J))*(SCOEF(IZ1,J+1)-
     &	SCOEF(IZ1,J))/(SCOEF(93,J+1)-SCOEF(93,J))
	EL=(RLAMB0+RLAMB1)/FLOAT(IZ1)**.33333
C
	ZETA0=Q+(1./(2.*VFERMI**2))*(1.-Q)*LOG(1+(4*EL*VFERMI/1.919)**2)
C	(ADD IZ1**3 EFFECT AS SHOWN IN REF. 779.)
	A=LOG(E)
	IF(A.LT.0)A=0
	ZETA=ZETA0*(1.+(1./IZ1**2)*(.08+.0015*IZ2)*EXP(-(7.6-A)**2))
	A=VRMIN/FLOAT(IZ1)**.6667
	IF(A.LT.YRMIN)A=YRMIN
	IF(YR.LE.A)GO TO 145
	IF(ISOGAS.EQ.1)THEN
	CALL PSTOPG(IZ1,IZ2,E,SP,IGASNO)
	ELSE
	CALL PSTOP(IZ2,E,SP)
	END IF
C
	SE=SP*(ZETA*IZ1)**2
C	(Add Fermi Velocity Correction - 1995)
C	(VFCORR is in SCOEF(IZ2,41-54))
C	(Interpolation values in SCOEF(93,41-54))
	EION=E
	IF(EION.GT.9999)EION=9999
C	(Not valid >1E4 keV/amu)
	DO 150 J=41,53
C	(Find E Interpolation)
	IF(EION.LT.SCOEF(93,J))GO TO 160
C	(in SCOEF(93,41-54))
150	CONTINUE
160	J=J-1
	IF(J.LT.41)J=41
	IF(J.GT.53)J=53
	VFCOR0=SCOEF(IZ2,J)
	VFCOR1=(EION-SCOEF(93,J))*(SCOEF(IZ2,J+1)-
     &	SCOEF(IZ2,J))/(SCOEF(93,J+1)-SCOEF(93,J))
	SE=SE*(VFCOR0+VFCOR1)
	GO TO 310
C
C
C	(CALCULATE VELOCITY STOPPING FOR YR LESS THAN YRMIN.)
C
145	A=YRMIN*FLOAT(IZ1)**.6667
	IF(VRMIN.LT.A)VRMIN=A
	A=VRMIN**2-0.8*VFERMI**2
	IF(A.LT.0)A=0
	VMIN=.5*(VRMIN+SQRT(A))
	EEE=25*VMIN**2
	IF(ISOGAS.EQ.1)THEN
	CALL PSTOPG(IZ1,IZ2,EEE,SP,IGASNO)
	ELSE
	CALL PSTOP(IZ2,EEE,SP)
	END IF
C
C	(Add Fermi Velocity Correction to Low Energy value)
C	(VFCORR is in SCOEF(IZ2,41-54))
C	(Interpolation values in SCOEF(93,41-54))
	EION=EEE
	IF(EION.GT.9999)EION=9999
C	(Not valid >1E4 keV/amu)
	DO 170 J=41,53
C	(Find E Interpolation)
	IF(EION.LT.SCOEF(93,J))GO TO 180
C	(in SCOEF(93,41-54))
170	CONTINUE
180	J=J-1
	IF(J.LT.41)J=41
	IF(J.GT.53)J=53
	VFCOR0=SCOEF(IZ2,J)
	VFCOR1=(EION-SCOEF(93,J))*(SCOEF(IZ2,J+1)-
     &	SCOEF(IZ2,J))/(SCOEF(93,J+1)-SCOEF(93,J))
	SP=SP*(VFCOR0+VFCOR1)
C
C	Following corrects for low-energy stopping, where little data exists.
C	Traditionally, this is velocity-proportional stopping, however for
C	light ions, light targets and semiconductors, a large variation exists.
C
C	(Note: HIPOWR down = Se up      LAMBDA down = Se down)
C
	HIPOWR=.47
C	(TRIM-88 used 0.50)
	IF(IZ1.EQ.3)THEN
	HIPOWR= 0.55
	GO TO 290
	END IF
	IF(IZ2.LT.7)THEN
	HIPOWR= 0.375
	GO TO 290
	END IF
C	(Following compensates for semiconductor band-gap)
	IF(IZ1.LT.18.AND.(IZ2.EQ.14.OR.IZ2.EQ.32))HIPOWR=0.375
C
290	SE=(SP*(ZETA*IZ1)**2)*(E/EEE)**HIPOWR
C
	END IF
C
C
310	SE=SE*10
C	(This converts Stopping in eV/(1E15-cm2) to eV-A2)
	RETURN
	END
C
C
C
	SUBROUTINE PSTOP(IZ2,E,SP)
C
C	(CALCULATES PROTON ELECTRONIC STOPPING POWERS IN SOLIDS)
C
	DIMENSION SCOEF(93,54),SCOEFG(24,11)
	COMMON/A/ SCOEF,SCOEFG
C
	IF(E.GT.1.E4)THEN
C	(High Energy Stopping (6/87))
	X=LOG(E)/E
	SP=SCOEF(IZ2,17)+(SCOEF(IZ2,18)*X)+
     &	(SCOEF(IZ2,19)*X*X)+(SCOEF(IZ2,20)/X)
	RETURN
	END IF
C
C	(VELOCITY PROPORTIONAL STOPPING BELOW VELOCITY PE0.)
	PE0=10
	PE=PE0
	IF(PE.LT.E)PE=E
	SL=(SCOEF(IZ2,9)*PE**SCOEF(IZ2,10))+
     &	SCOEF(IZ2,11)*PE**SCOEF(IZ2,12)
	SH=SCOEF(IZ2,13)/PE**SCOEF(IZ2,14)*
     &	LOG((SCOEF(IZ2,15)/PE)+SCOEF(IZ2,16)*PE)
	SP=SL*SH/(SL+SH)
	IF(E.GT.PE0)RETURN
C
C	(PPOWER IS THE POWER OF VELOCITY STOPPING BELOW PE0.)
C
	PPOWER=0.45
C	(Low Energy Stopping: S(E)**0.45)
	IF(IZ2.LT.7)PPOWER=PPOWER-0.1
C	(Z2=3-6 has low S(E)**0.35)
	SP=SP*(E/PE0)**PPOWER
	RETURN
C
	END
C
C
C
	SUBROUTINE PSTOPG(IZ1,IZ2,E,SP,IGASNO)
C
C 	Calculates proton Stopping Powers in GASES
C
	DIMENSION SCOEF(93,54),SCOEFG(24,11)
	COMMON/A/ SCOEF,SCOEFG
C
	IF(E.GT.1.E4)THEN
C	(High Energy Stopping is Normal)
	X=LOG(E)/E
	SP=SCOEF(IZ2,17)+(SCOEF(IZ2,18)*X)+
     &	(SCOEF(IZ2,19)*X*X)+(SCOEF(IZ2,20)/X)
	RETURN
	END IF
C
	II=IGASNO
C	(IGASNO is from 1 to 8)
	IF(IZ1.GT.1)II=IGASNO+8
	IF(IZ1.GT.2)II=IGASNO+16
C
C	(VELOCITY PROPORTIONAL STOPPING BELOW VELOCITY PE0.)
C
	PE0=10
	PE=PE0
	IF(PE.LT.E)PE=E
	SL=(SCOEFG(II,04)*PE**SCOEFG(II,05))+
     &	SCOEFG(II,06)*PE**SCOEFG(II,07)
	SH=SCOEFG(II,08)/PE**SCOEFG(II,09)*LOG((SCOEFG(II,10)/PE)+
     &	SCOEFG(II,11)*PE)
	SP=SL*SH/(SL+SH)
	IF(E.GT.PE0)RETURN
C
C	(PPOWER IS THE POWER OF VELOCITY STOPPING BELOW PE0.)
C
	PPOWER=0.45
C	(Low Energy Stopping: S(E)**0.45)
	IF(IZ2.LT.7)PPOWER=0.35
C	(Z2=3-6 has low S(E)**0.35)
	SP=SP*(E/PE0)**PPOWER
	RETURN
C
	END
C
C
C
	SUBROUTINE STOPCO
C
C	Gets array SCOEF.95A+SCOEF.95B (93,54)
C	Gets array SCOEFGAS.95 (24,11)
C	Contains all Data and Coefficients
C
C     ษออออออออออออออออหหออออออออออออหหอออออออออออออออออออออออป
C     ฬออออออออออออออออสผ  SCOEF.95A ศสอออออออออออออออออออออออน
C     ฬออ> Columns 1-8 :                                      บ
C     บ    (1=Atomic Number),(2=Mass of MAI),(3=Weight of MAI)บ
C     บ    (4=Natural Mass),(5=Target density: g/cm3),        บ
C     บ    (6=Target Density: atoms/cm3),(7=Fermi Velocity)   บ
C     บ    (8=Sublimation Energy)                             บ
C     ฬออ> Columns  9-16 : Proton stopping coefficients.      บ
C     บ                 ษหอออออออออออหป                       บ
C     ฬอออออออออออออออออสผ SCOEF.95B ศสอออออออออออออออออออออออน
C     ฬออ> Columns 1 - 4 : Proton stopping > 10 MeV.          บ
C     ฬออ> Columns  5-13 : Ion Lambda values (width of ion)   บ
C     ฬออ> Columns 14-28 : Tgt Fermi Velocity Correction (<5%)บ
C     ศอออออออออออออออออออออออออออออออออออออออออออออออออออออออผ
C
	DIMENSION  SCOEF(93,54), SCOEFG(24,11)
	COMMON/A/ SCOEF,SCOEFG
C
	OPEN(95,FILE='SCOEF.95A',STATUS='OLD')
	READ(95,'(/)')
	DO 110 I=1,93
110	READ(95,*)(SCOEF(I,J),J=1,16)
	CLOSE(95)
C
	OPEN(95,FILE='SCOEF.95B',STATUS='OLD')
	READ(95,'(/)')
	DO 120 I=1,93
120	READ(95,*)(SCOEF(I,J),J=17,54)
	CLOSE(95)
C
	OPEN(95,FILE='SCOEFGAS.95',STATUS='OLD')
	READ(95,'(/)')
	DO 130 I=1,24
130	READ(95,*)(SCOEFG(I,J),J=1,11)
	CLOSE(95)
C
	RETURN
	END
c
c