	 SUBROUTINE INTG(F,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
     1   IER,LIMIT,LENW,LAST,IWORK,WORK)
cc*********************************************************************
cc	Routine d'int�gration
cc*********************************************************************
      REAL ABSERR,BOUND,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,INF,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(LIMIT),WORK(LENW)
C
      EXTERNAL F
C
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0E+00
      ABSERR = 0.0E+00
      IF(LIMIT.LT.1.OR.LENW.LT.LIMIT*4) GO TO 10
C
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
C
      CALL INTGE(F,BOUND,INF,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,
     1  NEVAL,IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF(IER.NE.0) CALL xerreur( 'ABNORMAL RETURN FROM  INTG',
     1  26,IER,LVL)
      RETURN
      END
cc
      SUBROUTINE EA(NEWFLG,SVALUE,LIMEXP,RESULT,ABSERR,EPSTAB,IERR)
      REAL ABSERR,DELTA1,DELTA2,DELTA3,EPRN,EPSTAB(*),
     1   ERROR,ERR1,ERR2,ERR3,E0,E1,E2,E3,RELPR,RES,RESULT,
     2	 RES3LA(3),r2mach,SS,SVALUE,TOL1,TOL2,TOL3
      INTEGER I,IB,IB2,IE,IERR,IN,K1,K2,K3,LIMEXP,N,NEWELM,NUM,NRES
      LOGICAL NEWFLG
C
      IF(LIMEXP.LT.3) THEN
        IERR = 1
	CALL xerreur('LIMEXP IS LESS THAN 3',21,1,1)
        GO TO 110
      ENDIF
      IERR = 0
      RES3LA(1)=EPSTAB(LIMEXP+5)
      RES3LA(2)=EPSTAB(LIMEXP+6)
      RES3LA(3)=EPSTAB(LIMEXP+7)
      RESULT=SVALUE
      IF(NEWFLG) THEN
        N=1
        NRES=0
        NEWFLG=.FALSE.
        EPSTAB(N)=SVALUE
        ABSERR=ABS(RESULT)
        GO TO 100
      ELSE
        N=INT(EPSTAB(LIMEXP+3))
        NRES=INT(EPSTAB(LIMEXP+4))
        IF(N.EQ.2) THEN
          EPSTAB(N)=SVALUE
          ABSERR=.6E+01*ABS(RESULT-EPSTAB(1))
          GO TO 100
        ENDIF
      ENDIF
      EPSTAB(N)=SVALUE
      RELPR=r2mach(4)
      EPRN=1.0E+01*RELPR
      EPSTAB(N+2)=EPSTAB(N)
      NEWELM=(N-1)/2
      NUM=N
      K1=N
      DO 40 I=1,NEWELM
        K2=K1-1
        K3=K1-2
        RES=EPSTAB(K1+2)
        E0=EPSTAB(K3)
        E1=EPSTAB(K2)
        E2=RES
        DELTA2=E2-E1
        ERR2=ABS(DELTA2)
        TOL2=MAX(ABS(E2),ABS(E1))*RELPR
        DELTA3=E1-E0
        ERR3=ABS(DELTA3)
        TOL3=MAX(ABS(E1),ABS(E0))*RELPR
        IF(ERR2.GT.TOL2.OR.ERR3.GT.TOL3) GO TO 10
C
        RESULT=RES
        ABSERR=ERR2+ERR3
        GO TO 50
   10   IF(I.NE.1) THEN
          E3=EPSTAB(K1)
          EPSTAB(K1)=E1
          DELTA1=E1-E3
          ERR1=ABS(DELTA1)
          TOL1=MAX(ABS(E1),ABS(E3))*RELPR
C
          IF(ERR1.LE.TOL1.OR.ERR2.LE.TOL2.OR.ERR3.LE.TOL3) GO TO 20
          SS=0.1E+01/DELTA1+0.1E+01/DELTA2-0.1E+01/DELTA3
        ELSE
          EPSTAB(K1)=E1
          IF(ERR2.LE.TOL2.OR.ERR3.LE.TOL3) GO TO 20
          SS=0.1E+01/DELTA2-0.1E+01/DELTA3
        ENDIF
C
        IF(ABS(SS*E1).GT.0.1E-03) GO TO 30
   20   N=I+I-1
        IF(NRES.EQ.0) THEN
          ABSERR=ERR2+ERR3
          RESULT=RES
        ELSE IF(NRES.EQ.1) THEN
          RESULT=RES3LA(1)
        ELSE IF(NRES.EQ.2) THEN
          RESULT=RES3LA(2)
        ELSE
          RESULT=RES3LA(3)
        ENDIF
        GO TO 50
C
   30   RES=E1+0.1E+01/SS
        EPSTAB(K1)=RES
        K1=K1-2
        IF(NRES.EQ.0) THEN
          ABSERR=ERR2+ABS(RES-E2)+ERR3
          RESULT=RES
          GO TO 40
        ELSE IF(NRES.EQ.1) THEN
          ERROR=.6E+01*(ABS(RES-RES3LA(1)))
        ELSE IF(NRES.EQ.2) THEN
          ERROR=.2E+01*(ABS(RES-RES3LA(2))+ABS(RES-RES3LA(1)))
        ELSE
          ERROR=ABS(RES-RES3LA(3))+ABS(RES-RES3LA(2))
     1          +ABS(RES-RES3LA(1))
        ENDIF
        IF(ERROR.GT.1.0E+01*ABSERR) GO TO 40
        ABSERR=ERROR
        RESULT=RES
   40 CONTINUE
C
        IF(NRES.EQ.1) THEN
          ABSERR=.6E+01*(ABS(RESULT-RES3LA(1)))
        ELSE IF(NRES.EQ.2) THEN
          ABSERR=.2E+01*ABS(RESULT-RES3LA(2))+ABS(RESULT-RES3LA(1))
        ELSE IF(NRES.GT.2) THEN
          ABSERR=ABS(RESULT-RES3LA(3))+ABS(RESULT-RES3LA(2))
     1          +ABS(RESULT-RES3LA(1))
        ENDIF
C
   50 IF(N.EQ.LIMEXP) N=2*(LIMEXP/2)-1
      IB=1
      IF((NUM/2)*2.EQ.NUM) IB=2
      IE=NEWELM+1
      DO 60 I=1,IE
        IB2=IB+2
        EPSTAB(IB)=EPSTAB(IB2)
        IB=IB2
   60 CONTINUE
      IF(NUM.EQ.N) GO TO 80
      IN=NUM-N+1
      DO 70 I=1,N
        EPSTAB(I)=EPSTAB(IN)
        IN=IN+1
   70 CONTINUE
C
   80 IF(NRES.EQ.0) THEN
        RES3LA(1)=RESULT
      ELSE IF(NRES.EQ.1) THEN
        RES3LA(2)=RESULT
      ELSE IF(NRES.EQ.2) THEN
        RES3LA(3)=RESULT
      ELSE
        RES3LA(1)=RES3LA(2)
        RES3LA(2)=RES3LA(3)
        RES3LA(3)=RESULT
      ENDIF
   90 ABSERR=MAX(ABSERR,EPRN*ABS(RESULT))
      NRES=NRES+1
  100 N=N+1
*     IF(N.LE.3) ABSERR = r2mach(2) * (0.1D-03)
      EPSTAB(LIMEXP+3)=REAL(N)
      EPSTAB(LIMEXP+4)=REAL(NRES)
      EPSTAB(LIMEXP+5)=RES3LA(1)
      EPSTAB(LIMEXP+6)=RES3LA(2)
      EPSTAB(LIMEXP+7)=RES3LA(3)
  110 RETURN
      END
cc
      SUBROUTINE INTGE(F,BOUND,INF,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,
     1   NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
C
      REAL ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,
     1  A2,BLIST,BOUN,BOUND,B1,B2,CORREC,DEFABS,DEFAB1,DEFAB2,
     2  DRES,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND,
     3  ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,ERTEST,F,RESABS,
     4	RESEPS,RESULT,RLIST,RLIST2,r2mach,SMALL,UFLOW
      INTEGER ID,IER,IERR,IERRO,INF,IORD,IROFF1,IROFF2,IROFF3,JUPBND,K,
     1  KSGN,KTMIN,LAST,LIMEXP,LIMIT,MAXERR,NEVAL,NRMAX
      LOGICAL EXTRAP,LERR,NEWFLG,NOEXT
C
      PARAMETER (LIMEXP = 50)
C
      DIMENSION ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),IORD(LIMIT),
     1  RLIST(LIMIT),RLIST2(LIMEXP+7)
C
      EXTERNAL F
C
       EPMACH = r2mach(4)
C
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0E+00
      ABSERR = 0.0E+00
      ALIST(1) = 0.0E+00
      BLIST(1) = 0.1E+01
      RLIST(1) = 0.0E+00
      ELIST(1) = 0.0E+00
      IORD(1) = 0
      NEWFLG = .TRUE.
      IF(EPSABS.LE.0.0E+00.AND.EPSREL.LT.AMAX1(0.5E+02*EPMACH,0.5E-14))
     1  IER = 6
      IF(IER.EQ.6) GO TO 999
C
      BOUN = BOUND
      IF(INF.EQ.2) BOUN = 0.0E+00
      CALL QK15I(F,BOUN,INF,0.0E+00,0.1E+01,RESULT,ABSERR,
     1  DEFABS,RESABS)
C
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
      DRES = ABS(RESULT)
      ERRBND = AMAX1(EPSABS,EPSREL*DRES)
      IF(ABSERR.LE.1.0E+02*EPMACH*DEFABS.AND.ABSERR.GT.
     1  ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS).OR.
     1  ABSERR.EQ.0.0E+00) GO TO 130
C
      UFLOW = r2mach(1) * 0.2E+01
      LERR = .FALSE.
      CALL EA(NEWFLG,RESULT,LIMEXP,RESEPS,ABSEPS,RLIST2,IERR)
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      NRMAX = 1
      KTMIN = 0
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      IERRO = 0
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      KSGN = -1
      IF(DRES.GE.(0.1E+01-0.5E+02*EPMACH)*DEFABS) KSGN = 1
C
      DO 90 LAST = 2,LIMIT
C
        A1 = ALIST(MAXERR)
        B1 = 0.5E+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL QK15I(F,BOUN,INF,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        CALL QK15I(F,BOUN,INF,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
C
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2)GO TO 15
        IF(ABS(RLIST(MAXERR)-AREA12).GT.0.1E-04*ABS(AREA12)
     1  .OR.ERRO12.LT.0.99E+00*ERRMAX) GO TO 10
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   10   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   15   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = AMAX1(EPSABS,EPSREL*ABS(AREA))
C
        IF(IROFF1+IROFF2.GE.10.OR.IROFF3.GE.20) IER = 2
        IF(IROFF2.GE.5) IERRO = 3
C
        IF(LAST.EQ.LIMIT) IER = 1
C
        IF(AMAX1(ABS(A1),ABS(B2)).LE.(0.1E+01+0.1E+03*EPMACH)*
     1  (ABS(A2)+0.1E+04*UFLOW)) IER = 4
C
        IF(ERROR2.GT.ERROR1) GO TO 20
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 30
   20   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
   30   CALL QPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
        IF(ERRSUM.LE.ERRBND) GO TO 115
        IF(IER.NE.0) GO TO 100
        IF(LAST.EQ.2) GO TO 80
        IF(NOEXT) GO TO 90
        ERLARG = ERLARG-ERLAST
        IF(ABS(B1-A1).GT.SMALL) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 40
C
        IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
        EXTRAP = .TRUE.
        NRMAX = 2
   40   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 60
C
        ID = NRMAX
        JUPBND = LAST
        IF(LAST.GT.(2+LIMIT/2)) JUPBND = LIMIT+3-LAST
        DO 50 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
          IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
          NRMAX = NRMAX+1
   50   CONTINUE
C
   60   CALL EA(NEWFLG,AREA,LIMEXP,RESEPS,ABSEPS,RLIST2,IERR)
        KTMIN = KTMIN+1
        IF((KTMIN.GT.5).AND.(ABSERR.LT.0.1E-02*ERRSUM).AND.(LERR))
     1     IER = 5
        IF((ABSEPS.GE.ABSERR).AND.(LERR)) GO TO 70
        KTMIN = 0
        ABSERR = ABSEPS
        LERR = .TRUE.
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = AMAX1(EPSABS,EPSREL*ABS(RESEPS))
        IF((ABSERR.LE.ERTEST).AND.(LERR)) GO TO 100
C
   70   IF(RLIST2(LIMEXP+3).EQ.1) NOEXT = .TRUE.
        IF(IER.EQ.5) GO TO 100
        MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        SMALL = SMALL*0.5E+00
        ERLARG = ERRSUM
        GO TO 90
   80   SMALL = 0.375E+00
        ERLARG = ERRSUM
        ERTEST = ERRBND
        CALL EA(NEWFLG,AREA,LIMEXP,RESEPS,ABSEPS,RLIST2,IERR)
   90 CONTINUE
C
  100 IF(.NOT.LERR) GO TO 115
      IF((IER+IERRO).EQ.0) GO TO 110
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0E+00.AND.AREA.NE.0.0E+00)GO TO 105
      IF(ABSERR.GT.ERRSUM)GO TO 115
      IF(AREA.EQ.0.0E+00) GO TO 130
      GO TO 110
  105 IF(ABSERR/ABS(RESULT).GT.ERRSUM/ABS(AREA))GO TO 115
C
  110 IF(KSGN.EQ.(-1).AND.AMAX1(ABS(RESULT),ABS(AREA)).LE.
     1 DEFABS*0.1E-01) GO TO 130
      IF(0.1E-01.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1E+03.
     1OR.ERRSUM.GT.ABS(AREA)) IER = 6
      GO TO 130
C
  115 RESULT = 0.0E+00
      DO 120 K = 1,LAST
        RESULT = RESULT+RLIST(K)
  120 CONTINUE
      ABSERR = ERRSUM
  130 NEVAL = 30*LAST-15
      IF(INF.EQ.2) NEVAL = 2*NEVAL
      IF(IER.GT.2) IER=IER-1
  999 RETURN
      END
cc
      SUBROUTINE QK15I(F,BOUN,INF,A,B,RESULT,ABSERR,RESABS,RESASC)
C
      REAL A,ABSC,ABSC1,ABSC2,ABSERR,B,BOUN,CENTR,
     1	DINF,r2mach,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,
     2  FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,TABSC1,TABSC2,
     3  UFLOW,WG,WGK,XGK
      INTEGER INF,J,MIN0
      EXTERNAL F
C
      DIMENSION FV1(7),FV2(7),XGK(8),WGK(8),WG(8)
C
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),
     1  XGK(8)/
     2     0.9914553711208126E+00,     0.9491079123427585E+00,
     3     0.8648644233597691E+00,     0.7415311855993944E+00,
     4     0.5860872354676911E+00,     0.4058451513773972E+00,
     5     0.2077849550078985E+00,     0.0000000000000000E+00/
C
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),
     1  WGK(8)/
     2     0.2293532201052922E-01,     0.6309209262997855E-01,
     3     0.1047900103222502E+00,     0.1406532597155259E+00,
     4     0.1690047266392679E+00,     0.1903505780647854E+00,
     5     0.2044329400752989E+00,     0.2094821410847278E+00/
C
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8)/
     1     0.0000000000000000E+00,     0.1294849661688697E+00,
     2     0.0000000000000000E+00,     0.2797053914892767E+00,
     3     0.0000000000000000E+00,     0.3818300505051189E+00,
     4     0.0000000000000000E+00,     0.4179591836734694E+00/
C
      EPMACH = r2mach(4)
      UFLOW = r2mach(1)
      DINF = MIN0(1,INF)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      TABSC1 = BOUN+DINF*(0.1E+01-CENTR)/CENTR
      FVAL1 = F(TABSC1)
      IF(INF.EQ.2) FVAL1 = FVAL1+F(-TABSC1)
      FC = (FVAL1/CENTR)/CENTR
C
      RESG = WG(8)*FC
      RESK = WGK(8)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,7
        ABSC = HLGTH*XGK(J)
        ABSC1 = CENTR-ABSC
        ABSC2 = CENTR+ABSC
        TABSC1 = BOUN+DINF*(0.1E+01-ABSC1)/ABSC1
        TABSC2 = BOUN+DINF*(0.1E+01-ABSC2)/ABSC2
        FVAL1 = F(TABSC1)
        FVAL2 = F(TABSC2)
        IF(INF.EQ.2) FVAL1 = FVAL1+F(-TABSC1)
        IF(INF.EQ.2) FVAL2 = FVAL2+F(-TABSC2)
        FVAL1 = (FVAL1/ABSC1)/ABSC1
        FVAL2 = (FVAL2/ABSC2)/ABSC2
        FV1(J) = FVAL1
        FV2(J) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(J)*FSUM
        RESABS = RESABS+WGK(J)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(8)*ABS(FC-RESKH)
      DO 20 J=1,7
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESASC = RESASC*HLGTH
      RESABS = RESABS*HLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.E0) ABSERR = RESASC*
     1 AMIN1(0.1E+01,(0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = AMAX1
     1 ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
cc
      SUBROUTINE QPSRT(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
C
      REAL ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
     1  NRMAX
      DIMENSION ELIST(LAST),IORD(LAST)
C
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
C
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
C
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
C
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
C
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
C
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
C
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
C
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
C
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END
c*********************************************************************
      REAL FUNCTION r2mach(I)
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      REAL RMACH(5) 
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C-------------------------------------------------------------
C     MACHINE CONSTANTS FOR THE IBM PC FAMILY (D. KAHANER NBS)
C-------------------------------------------------------------
      DATA RMACH/1.18E-38,3.40E+38,0.595E-07,1.19E-07,0.30102999566/
C-------------------------------------------------------------------
      IF((I.LT.1).OR.(I.GT.5))
     s CALL xerreur ( 'r2mach -- I OUT OF BOUNDS',25,1,2)
      r2mach = RMACH(I)
      RETURN
      END
c**********************************************************************
      SUBROUTINE xerreur(MESSG,NMESSG,NERR,LEVEL)
      CHARACTER*(*) MESSG
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
cc
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
      CHARACTER*(*) MESSG
      write(*,*) nmessg,nerr
      WRITE(*,*) MESSG
      IF(NI.EQ.2)THEN
        WRITE(*,*) I1,I2
      ELSEIF(NI.EQ.1) THEN
        WRITE(*,*) I1
      ENDIF
      IF(NR.EQ.2) THEN
        WRITE(*,*) R1,R2
      ELSEIF(NR.EQ.1) THEN
        WRITE(*,*) R1
      ENDIF
      IF(ABS(LEVEL).LT.2)RETURN
      STOP 
      END
c*********************************************************************
      INTEGER FUNCTION ISAMAX(N,SX,INCX)
      REAL SX(*),SMAX,XMAG
      ISAMAX = 0
      IF(N.LE.0) RETURN
      ISAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
C
      SMAX = ABS(SX(1))
      NS = N*INCX
      II = 1
          DO 10 I=1,NS,INCX
          XMAG = ABS(SX(I))
          IF(XMAG.LE.SMAX) GO TO 5
          ISAMAX = II
          SMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
C
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N 
         XMAG = ABS(SX(I))
         IF(XMAG.LE.SMAX) GO TO 30
         ISAMAX = I 
         SMAX = XMAG
   30 CONTINUE
      RETURN
      END
cc
      REAL FUNCTION SASUM(N,SX,INCX)
C
      REAL SX(*)
      SASUM = 0.0E0 
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
      NS = N*INCX
          DO 10 I=1,NS,INCX
          SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      RETURN
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        SASUM = SASUM + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        SASUM = SASUM + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))
     1  + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
      RETURN
      END
cc
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
C
      REAL SX(*),SY(*),SA
      IF(N.LE.0.OR.SA.EQ.0.E0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
      DO 10 I = 1,N 
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END
cc
      SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
C
      REAL SX(*),SY(*)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
      DO 10 I = 1,N 
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        SY(I) = SX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I + 1) = SX(I + 1) 
        SY(I + 2) = SX(I + 2) 
        SY(I + 3) = SX(I + 3) 
        SY(I + 4) = SX(I + 4) 
        SY(I + 5) = SX(I + 5) 
        SY(I + 6) = SX(I + 6) 
   50 CONTINUE
      RETURN
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SX(I)
   70     CONTINUE
      RETURN
      END
cc
      REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
C
      REAL SX(*),SY(*)
      SDOT = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
    5 CONTINUE
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
      DO 10 I = 1,N 
        SDOT = SDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        SDOT = SDOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SDOT = SDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +
     1   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
      RETURN
C
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        SDOT = SDOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END
cc
      REAL FUNCTION SNRM2(N,SX,INCX)
      INTEGER          NEXT
      REAL   SX(*),  CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO, ONE
      DATA   ZERO, ONE /0.0E0, 1.0E0/
C
      DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      IF(N .GT. 0) GO TO 10
         SNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX 
C
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( ABS(SX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
   50 IF( SX(I) .EQ. ZERO) GO TO 200
      IF( ABS(SX(I)) .GT. CUTLO) GO TO 85
C
      ASSIGN 70 TO NEXT
      GO TO 105
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / SX(I)) / SX(I)
  105 XMAX = ABS(SX(I))
      GO TO 115
C
   70 IF( ABS(SX(I)) .GT. CUTLO ) GO TO 75
C
  110 IF( ABS(SX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / SX(I))**2
         XMAX = ABS(SX(I))
         GO TO 200
C
  115 SUM = SUM + (SX(I)/XMAX)**2
      GO TO 200
C
   75 SUM = (SUM * XMAX) * XMAX
C
   85 HITEST = CUTHI/FLOAT( N )
C
      DO 95 J =I,NN,INCX
      IF(ABS(SX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + SX(J)**2 
      SNRM2 = SQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
      SNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END
cc
      SUBROUTINE SSCAL(N,SA,SX,INCX)
C
      REAL SA,SX(*) 
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
      NS = N*INCX
          DO 10 I = 1,NS,INCX 
          SX(I) = SA*SX(I)
   10     CONTINUE
      RETURN
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END
cc
      SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
C
      REAL SX(*),SY(*),STEMP1,STEMP2,STEMP3
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1 
      DO 10 I = 1,N 
        STEMP1 = SX(IX)
        SX(IX) = SY(IY)
        SY(IY) = STEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40 
      DO 30 I = 1,M 
        STEMP1 = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        STEMP1 = SX(I)
        STEMP2 = SX(I+1)
        STEMP3 = SX(I+2)
        SX(I) = SY(I)
        SX(I+1) = SY(I+1)
        SX(I+2) = SY(I+2)
        SY(I) = STEMP1
        SY(I+1) = STEMP2
        SY(I+2) = STEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
C
      NS = N*INCX
        DO 70 I=1,NS,INCX
        STEMP1 = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP1
   70   CONTINUE
      RETURN
      END
c*********************************************************************
      REAL FUNCTION UNI()
      PARAMETER(
     *    CSAVE=362436./16777216.  ,
     *    CD=7654321./16777216.,
     *    CM=16777213./16777216.  )
C                            2**24=16777216
      REAL U(17),S,T,USTART,C,UNIB
      INTEGER I,J,II,JJ,K,KK,I1,J1,K1,L1,M1,ISEED
C      
      SAVE U,I,J,K,C
C
      DATA U/
     *0.8668672834288,  0.3697986366357,  0.8008968294805,
     *0.4173889774680,  0.8254561579836,  0.9640965269077,
     *0.4508667414265,  0.6451309529668,  0.1645456024730,
     *0.2787901807898,  0.06761531340295, 0.9663226330820,
     *0.01963343943798, 0.02947398211399, 0.1636231515294,
     *0.3976343250467,  0.2631008574685/
      DATA I,J,K,C/17,5,24,CSAVE/
C 
      UNI = U(I)-U(J)
      IF(UNI.LT.0.0)UNI = UNI+1.0
      U(I) = UNI
      I = I-1
      IF(I.EQ.0)I = 17
      J = J-1
      IF(J.EQ.0)J = 17
C
      C = C-CD
      IF(C.LT.0.0) C=C+CM
C
      UNI = UNI-C
      IF(UNI.LT.0.0)UNI = UNI+1.0
      RETURN
C
      ENTRY USTART(ISEED)
C
        I1 = MOD(ABS(ISEED),177)+1
        J1 = MOD(ABS(ISEED),167)+1
        K1 = MOD(ABS(ISEED),157)+1
        L1 = MOD(ABS(ISEED),147)+1
C 
        DO 2 II = 1,17
          S = 0.0
          T = 0.5
C
          DO 3 JJ = 1,K
                  M1 = MOD(MOD(I1*J1,179)*K1,179)
                  I1 = J1
                  J1 = K1
                  K1 = M1
                  L1 = MOD(53*L1+1,169)
                  IF(MOD(L1*M1,64).GE.32)S=S+T
    3             T = .5*T
    2   U(II) = S
        USTART = FLOAT(ISEED)
        RETURN
C
      ENTRY UNIB(KK)
        IF(KK.LE.24)THEN
             K=24
        ELSE
             K=KK
        ENDIF
        UNIB=FLOAT(K)
      END
cc*********************************************************************
