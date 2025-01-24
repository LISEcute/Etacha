	SUBROUTINE EQDIF (N,T,Y,F,TOUT,MSTATE,NROOT,EPS,EWT,MINT,WORK,
     s	 LENW,IWORK,LENIW,G)
cc*********************************************************************
cc	routine de resolution des equations differentielles
cc*********************************************************************
      EXTERNAL F, G
      REAL EPS, EWT, EWTCOM(1), G, HMAX, T, TOUT,
     8     WORK(*), Y(*)
      INTEGER IWORK(*)
      CHARACTER MSG*81
      PARAMETER(IMPL = 0, MXSTEP = 1000)
      IF (MINT .LT. 1 .OR. MINT .GT. 3) THEN
	WRITE(MSG, '(''EQDIF1FE Illegal input.  Improper value for '',
     8  ''the integration method flag,'', I8)') MINT
	CALL xerreur(MSG(1:81), 81, 21, 2)
        RETURN
      END IF
      IF (MSTATE .GE. 0) THEN
        NSTATE = MSTATE
        NTASK = 1
      ELSE
        NSTATE = - MSTATE
        NTASK = 3
      END IF
      EWTCOM(1) = EWT
      IF (EWT .NE. 0.E0) THEN
        IERROR = 3
      ELSE
        IERROR = 2
      END IF
      IF (MINT .EQ. 1) THEN
        MITER = 0
        MXORD = 12
      ELSE IF (MINT .EQ. 2) THEN
        MITER = 2
        MXORD = 5
      ELSE IF (MINT .EQ. 3) THEN
        MITER = 2
        MXORD = 12
      END IF
      HMAX = 2.E0*ABS(TOUT - T)
      CALL EQDIF3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWTCOM,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENW, IWORK, LENIW, F, F, NDE, MXSTEP, G, F)
      IF (MSTATE .GE. 0) THEN
        MSTATE = NSTATE
      ELSE
        MSTATE = - NSTATE
      END IF
      END
c***********************************************************************      
      SUBROUTINE SDCOR (DFDY,EL,FA,H,IMPL,IPVT,MATDIM,MITER,ML,MU,N,
     8   NDE,NQ,T,USERS,Y,YH,YWT,EVALFA,SAVE1,SAVE2,A,D,JSTATE)
C***ROUTINES CALLED  SGESL,SGBSL,SNRM2
      REAL A(MATDIM,*), D, DFDY(MATDIM,*), EL(13,12), H,
     8     SAVE1(*), SAVE2(*), SNRM2, T, Y(*), YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL EVALFA
      IF (MITER .EQ. 0) THEN
        DO 100 I = 1,N
 100      SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/YWT(I)
        D = SNRM2(N, SAVE1, 1)/SQRT(REAL(N))
        DO 105 I = 1,N
 105      SAVE1(I) = H*SAVE2(I) - YH(I,2)
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (IMPL .EQ. 0) THEN
          DO 130 I = 1,N
 130        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 150 I = 1,N
 150        SAVE2(I) = H*SAVE2(I)
          DO 160 J = 1,N
            DO 160 I = 1,N
 160          SAVE2(I) = SAVE2(I) - A(I,J)*(YH(J,2) + SAVE1(J))
        ELSE IF (IMPL .EQ. 2) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 180 I = 1,N
 180        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        END IF
        CALL SGESL (DFDY, MATDIM, N, IPVT, SAVE2, 0)
        DO 200 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 200      SAVE2(I) = SAVE2(I)/YWT(I)
        D = SNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (IMPL .EQ. 0) THEN
          DO 230 I = 1,N
 230        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 250 I = 1,N
 250        SAVE2(I) = H*SAVE2(I)
          MW = ML + 1 + MU
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 260 I = I1,I2
              I3 = I + J - MW
 260          SAVE2(I3) = SAVE2(I3) - A(I,J)*(YH(J,2) + SAVE1(J))
        ELSE IF (IMPL .EQ. 2) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 280 I = 1,N
 280        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        END IF
        CALL SGBSL (DFDY, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
        DO 300 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 300      SAVE2(I) = SAVE2(I)/YWT(I)
        D = SNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 2
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        IF (N .EQ. 0) THEN
          JSTATE = 10
          RETURN
        END IF
        DO 320 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 320      SAVE2(I) = SAVE2(I)/YWT(I)
        D = SNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      END IF
      END
c************************************************************      
      SUBROUTINE SDCST (MAXORD,MINT,ISWFLG,EL,TQ)
      REAL EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)
      FACTRL(1) = 1.E0
      DO 10 I = 2,MAXORD
 10     FACTRL(I) = REAL(I)*FACTRL(I-1)
      IF (MINT .EQ. 1) THEN
        GAMMA(1) = 1.E0
        DO 40 I = 1,MAXORD+1
          SUM = 0.E0
          DO 30 J = 1,I
 30         SUM = SUM - GAMMA(J)/REAL(I-J+2)
 40       GAMMA(I+1) = SUM
        EL(1,1) = 1.E0
        EL(2,1) = 1.E0
        EL(2,2) = 1.E0
        EL(3,2) = 1.E0
        DO 60 J = 3,MAXORD
          EL(2,J) = FACTRL(J-1)
          DO 50 I = 3,J
 50         EL(I,J) = REAL(J-1)*EL(I,J-1) + EL(I-1,J-1)
 60       EL(J+1,J) = 1.E0
        DO 80 J = 2,MAXORD
          EL(1,J) = EL(1,J-1) + GAMMA(J)
          EL(2,J) = 1.E0
          DO 80 I = 3,J+1
 80         EL(I,J) = EL(I,J)/(REAL(I-1)*FACTRL(J-1))
        DO 100 J = 1,MAXORD
          TQ(1,J) = -1.E0/(FACTRL(J)*GAMMA(J))
          TQ(2,J) = -1.E0/GAMMA(J+1)
 100      TQ(3,J) = -1.E0/GAMMA(J+2)
      ELSE IF (MINT .EQ. 2) THEN
        EL(1,1) = 1.E0
        EL(2,1) = 1.E0
        DO 130 J = 2,MAXORD
          EL(1,J) = FACTRL(J)
          DO 120 I = 2,J
 120        EL(I,J) = REAL(J)*EL(I,J-1) + EL(I-1,J-1)
 130      EL(J+1,J) = 1.E0
        SUM = 1.E0
        DO 150 J = 2,MAXORD
          SUM = SUM + 1.E0/REAL(J)
          DO 150 I = 1,J+1
 150        EL(I,J) = EL(I,J)/(FACTRL(J)*SUM)
        DO 170 J = 1,MAXORD
          IF (J .GT. 1) TQ(1,J) = 1.E0/FACTRL(J-1)
          TQ(2,J) = REAL(J+1)/EL(1,J)
 170      TQ(3,J) = REAL(J+2)/EL(1,J)
      END IF
      IF (ISWFLG .EQ. 3) THEN
        MXRD = MIN(MAXORD, 5)
        IF (MINT .EQ. 2) THEN
          GAMMA(1) = 1.E0
          DO 190 I = 1,MXRD
            SUM = 0.E0
            DO 180 J = 1,I
 180          SUM = SUM - GAMMA(J)/REAL(I-J+2)
 190        GAMMA(I+1) = SUM
        END IF
        SUM = 1.E0
        DO 200 I = 2,MXRD
          SUM = SUM + 1.E0/REAL(I)
 200      EL(1+I,1) = -REAL(I+1)*SUM*GAMMA(I+1)
      END IF
      END
c*********************************************************************      
      SUBROUTINE SDNTL (EPS,F,FA,HMAX,HOLD,IMPL,JTASK,MATDIM,MAXORD,
     8   MINT,MITER,ML,MU,N,NDE,SAVE1,T,UROUND,USERS,Y,YWT,H,MNTOLD,
     8   MTROLD,NFE,RC,YH,A,CONVRG,EL,FAC,IER,IPVT,NQ,NWAIT,RH,RMAX,
     8   SAVE2,TQ,TREND,ISWFLG,JSTATE)
      REAL A(MATDIM,*), EL(13,12), EPS, FAC(*), H, HMAX,
     8     HOLD, OLDL0, RC, RH, RMAX, RMINIT, SAVE1(*), SAVE2(*), SMAX,
     8     SMIN, SNRM2, SUM, SUM0, T, TQ(3,12), TREND, UROUND, Y(*),
     8     YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL CONVRG, IER
      PARAMETER(RMINIT = 10000.E0)
      IER = .FALSE.
      IF (JTASK .GE. 0) THEN
        IF (JTASK .EQ. 0) THEN
          CALL SDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RMAX = RMINIT
        END IF
        RC = 0.E0
        CONVRG = .FALSE.
        TREND = 1.E0
        NQ = 1
        NWAIT = 3
        CALL F (N, T, Y, SAVE2)
        IF (N .EQ. 0) THEN
          JSTATE = 6
          RETURN
        END IF
        NFE = NFE + 1
        IF (IMPL .NE. 0) THEN
          IF (MITER .EQ. 3) THEN
            IFLAG = 0
            CALL USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N,
     8                  NDE, IFLAG)
            IF (N .EQ. 0) THEN
              JSTATE = 10
              RETURN
            END IF
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL SGEFA (A, MATDIM, N, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL SGESL (A, MATDIM, N, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL SGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL SGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
            DO 150 I = 1,NDE
              IF (A(I,1) .EQ. 0.E0) THEN
                IER = .TRUE.
                RETURN
              ELSE
                SAVE2(I) = SAVE2(I)/A(I,1)
              END IF
 150          CONTINUE
            DO 155 I = NDE+1,N
 155          A(I,1) = 0.E0
          END IF
        END IF
        DO 170 I = 1,NDE
 170      SAVE1(I) = SAVE2(I)/YWT(I)
        SUM = SNRM2(NDE, SAVE1, 1)
        SUM0 = 1.E0/MAX(1.E0, ABS(T))
        SMAX = MAX(SUM0, SUM)
        SMIN = MIN(SUM0, SUM)
        SUM = SMAX*SQRT(1.E0 + (SMIN/SMAX)**2)/SQRT(REAL(NDE))
        H = SIGN(MIN(2.E0*EPS/SUM, ABS(H)), H)
        DO 180 I = 1,N
 180      YH(I,2) = H*SAVE2(I)
        IF (MITER .EQ. 2 .OR. MITER .EQ. 5 .OR. ISWFLG .EQ. 3) THEN
          DO 20 I = 1,N
 20         FAC(I) = SQRT(UROUND)
        END IF
      ELSE
        IF (MITER .NE. MTROLD) THEN
          MTROLD = MITER
          RC = 0.E0
          CONVRG = .FALSE.
        END IF
        IF (MINT .NE. MNTOLD) THEN
          MNTOLD = MINT
          OLDL0 = EL(1,NQ)
          CALL SDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RC = RC*EL(1,NQ)/OLDL0
          NWAIT = NQ + 2
        END IF
        IF (H .NE. HOLD) THEN
          NWAIT = NQ + 2
          RH = H/HOLD
          CALL SDSCL (HMAX, N, NQ, RMAX,  HOLD, RC, RH, YH)
        END IF
      END IF
      END
      SUBROUTINE SDNTP (H,K,N,NQ,T,TOUT,YH,Y)
      REAL FACTOR, H, R, T, TOUT, Y(*), YH(N,*)
      IF (K .EQ. 0) THEN
        DO 10 I = 1,N
 10       Y(I) = YH(I,NQ+1)
        R = ((TOUT - T)/H)
        DO 20 JJ = 1,NQ
          J = NQ + 1 - JJ
          DO 20 I = 1,N
 20         Y(I) = YH(I,J) + R*Y(I)
      ELSE
        KUSED = MIN(K, NQ)
        FACTOR = 1.E0
        DO 40 KK = 1,KUSED
 40       FACTOR = FACTOR*REAL(NQ+1-KK)
        DO 50 I = 1,N
 50       Y(I) = FACTOR*YH(I,NQ+1)
        DO 80 JJ = KUSED+1,NQ
          J = K + 1 + NQ - JJ
          FACTOR = 1.E0
          DO 60 KK = 1,KUSED
 60         FACTOR = FACTOR*REAL(J-KK)
          DO 70 I = 1,N
 70         Y(I) = FACTOR*YH(I,J) + R*Y(I)
 80       CONTINUE
        DO 100 I = 1,N
 100      Y(I) = Y(I)*H**(-KUSED)
      END IF
      END
c********************************************************      
      SUBROUTINE SDPSC (KSGN,N,NQ,YH)
      REAL YH(N,*)
      IF (KSGN .GT. 0) THEN
        DO 10 J1 = 1,NQ
          DO 10 J2 = J1,NQ
            J = NQ - J2 + J1
            DO 10 I = 1,N
 10           YH(I,J) = YH(I,J) + YH(I,J+1)
      ELSE
        DO 30 J1 = 1,NQ
          DO 30 J2 = J1,NQ
            J = NQ - J2 + J1
            DO 30 I = 1,N
 30           YH(I,J) = YH(I,J) - YH(I,J+1)
      END IF
      END
c**********************************************************************        
      SUBROUTINE SDPST (EL,F,FA,H,IMPL,JACOBN,MATDIM,MITER,ML,MU,N,NDE,
     8   NQ,SAVE2,T,USERS,Y,YH,YWT,UROUND,NFE,NJE,A,DFDY,FAC,IER,IPVT,
     8   SAVE1,ISWFLG,BND,JSTATE)
C***ROUTINES CALLED  SGEFA,SGBFA,SNRM2
      REAL A(MATDIM,*), BL, BND, BP, BR, BU, DFDY(MATDIM,*),
     8     DFDYMX, DIFF, DY, EL(13,12), FAC(*), FACMAX, FACMIN, FACTOR,
     8     H, SAVE1(*), SAVE2(*), SCALE, SNRM2, T, UROUND, Y(*),
     8     YH(N,*), YJ, YS, YWT(*)
      INTEGER IPVT(*)
      LOGICAL IER
      PARAMETER(FACMAX = .5E0)
      NJE = NJE + 1
      IER = .FALSE.
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (MITER .EQ. 1) THEN
          CALL JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          IF (ISWFLG .EQ. 3) BND = SNRM2(N*N, DFDY, 1)
          FACTOR = -EL(1,NQ)*H
          DO 110 J = 1,N
            DO 110 I = 1,N
 110          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 2) THEN
          BR = UROUND**(.875E0)
          BL = UROUND**(.75E0)
          BU = UROUND**(.25E0)
          BP = UROUND**(-.15E0)
          FACMIN = UROUND**(.78E0)
          DO 170 J = 1,N
            YS = MAX(ABS(YWT(J)), ABS(Y(J)))
 120        DY = FAC(J)*YS
            IF (DY .EQ. 0.E0) THEN
              IF (FAC(J) .LT. FACMAX) THEN
                FAC(J) = MIN(100.E0*FAC(J), FACMAX)
                GO TO 120
              ELSE
                DY = YS
              END IF
            END IF
            IF (NQ .EQ. 1) THEN
              DY = SIGN(DY, SAVE2(J))
            ELSE
              DY = SIGN(DY, YH(J,3))
            END IF
            DY = (Y(J) + DY) - Y(J)
            YJ = Y(J)
            Y(J) = Y(J) + DY
            CALL F (N, T, Y, SAVE1)
            IF (N .EQ. 0) THEN
              JSTATE = 6
              RETURN
            END IF
            Y(J) = YJ
            FACTOR = -EL(1,NQ)*H/DY
            DO 140 I = 1,N
 140          DFDY(I,J) = (SAVE1(I) - SAVE2(I))*FACTOR
C                                                                 Step 1
            DIFF = ABS(SAVE2(1) - SAVE1(1))
            IMAX = 1
            DO 150 I = 2,N
              IF (ABS(SAVE2(I) - SAVE1(I)) .GT. DIFF) THEN
                IMAX = I
                DIFF = ABS(SAVE2(I) - SAVE1(I))
              END IF
 150          CONTINUE
C                                                                 Step 2
            IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT. 0.E0) THEN
              SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
C                                                                 Step 3
              IF (DIFF .GT. BU*SCALE) THEN
                FAC(J) = MAX(FACMIN, FAC(J)*.1E0)
              ELSE IF (BR*SCALE .LE. DIFF .AND. DIFF .LE. BL*SCALE) THEN
                FAC(J) = MIN(FAC(J)*10.E0, FACMAX)
C                                                                 Step 4
              ELSE IF (DIFF .LT. BR*SCALE) THEN
                FAC(J) = MIN(BP*FAC(J), FACMAX)
              END IF
            END IF
 170        CONTINUE
          IF (ISWFLG .EQ. 3) BND = SNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)
          NFE = NFE + N
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 190 I = 1,N
 190        DFDY(I,I) = DFDY(I,I) + 1.E0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 210 J = 1,N
            DO 210 I = 1,N
 210          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 230 I = 1,NDE
 230        DFDY(I,I) = DFDY(I,I) + A(I,1)
        END IF
        CALL SGEFA (DFDY, MATDIM, N, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (MITER .EQ. 4) THEN
          CALL JACOBN (N, T, Y, DFDY(ML+1,1), MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          FACTOR = -EL(1,NQ)*H
          MW = ML + MU + 1
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 260 I = I1,I2
 260          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 5) THEN
          BR = UROUND**(.875E0)
          BL = UROUND**(.75E0)
          BU = UROUND**(.25E0)
          BP = UROUND**(-.15E0)
          FACMIN = UROUND**(.78E0)
          MW = ML + MU + 1
          J2 = MIN(MW, N)
          DO 340 J = 1,J2
            DO 290 K = J,N,MW
              YS = MAX(ABS(YWT(K)), ABS(Y(K)))
 280          DY = FAC(K)*YS
              IF (DY .EQ. 0.E0) THEN
                IF (FAC(K) .LT. FACMAX) THEN
                  FAC(K) = MIN(100.E0*FAC(K), FACMAX)
                  GO TO 280
                ELSE
                  DY = YS
                END IF
              END IF
              IF (NQ .EQ. 1) THEN
                DY = SIGN(DY, SAVE2(K))
              ELSE
                DY = SIGN(DY, YH(K,3))
              END IF
              DY = (Y(K) + DY) - Y(K)
              DFDY(MW,K) = Y(K)
 290          Y(K) = Y(K) + DY
            CALL F (N, T, Y, SAVE1)
            IF (N .EQ. 0) THEN
              JSTATE = 6
              RETURN
            END IF
            DO 330 K = J,N,MW
              Y(K) = DFDY(MW,K)
              YS = MAX(ABS(YWT(K)), ABS(Y(K)))
              DY = FAC(K)*YS
              IF (DY .EQ. 0.E0) DY = YS
              IF (NQ .EQ. 1) THEN
                DY = SIGN(DY, SAVE2(K))
              ELSE
                DY = SIGN(DY, YH(K,3))
              END IF
              DY = (Y(K) + DY) - Y(K)
              FACTOR = -EL(1,NQ)*H/DY
              I1 = MAX(ML+1, MW+1-K)
              I2 = MIN(MW+N-K, MW+ML)
              DO 300 I = I1,I2
                I3 = K + I - MW
 300            DFDY(I,K) = FACTOR*(SAVE1(I3) - SAVE2(I3))
C                                                                 Step 1
              IMAX = MAX(1, K - MU)
              DIFF = ABS(SAVE2(IMAX) - SAVE1(IMAX))
              I1 = IMAX
              I2 = MIN(K + ML, N)
              DO 310 I = I1+1,I2
                IF (ABS(SAVE2(I) - SAVE1(I)) .GT. DIFF) THEN
                  IMAX = I
                  DIFF = ABS(SAVE2(I) - SAVE1(I))
                END IF
 310            CONTINUE
C                                                                 Step 2
              IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT.0.E0) THEN
                SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
C                                                                 Step 3
                IF (DIFF .GT. BU*SCALE) THEN
                  FAC(K) = MAX(FACMIN, FAC(K)*.1E0)
                ELSE IF (BR*SCALE .LE.DIFF .AND. DIFF .LE.BL*SCALE) THEN
                  FAC(K) = MIN(FAC(K)*10.E0, FACMAX)
C                                                                 Step 4
                ELSE IF (DIFF .LT. BR*SCALE) THEN
                  FAC(K) = MIN(BP*FAC(K), FACMAX)
                END IF
              END IF
 330          CONTINUE
 340        CONTINUE
          NFE = NFE + J2
        END IF
        IF (ISWFLG .EQ. 3) THEN
          DFDYMX = 0.E0
          DO 345 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 345 I = I1,I2
 345          DFDYMX = MAX(DFDYMX, ABS(DFDY(I,J)))
          BND = 0.E0
          IF (DFDYMX .NE. 0.E0) THEN
            DO 350 J = 1,N
              I1 = MAX(ML+1, MW+1-J)
              I2 = MIN(MW+N-J, MW+ML)
              DO 350 I = I1,I2
 350            BND = BND + (DFDY(I,J)/DFDYMX)**2
            BND = DFDYMX*SQRT(BND)/(-EL(1,NQ)*H)
          END IF
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 360 J = 1,N
 360        DFDY(MW,J) = DFDY(MW,J) + 1.E0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 380 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 380 I = I1,I2
 380          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 400 J = 1,NDE
 400        DFDY(MW,J) =  DFDY(MW,J) + A(J,1)
        END IF
        CALL SGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 1
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        IF (N .EQ. 0) THEN
          JSTATE = 10
          RETURN
        END IF
      END IF
      END
c***********************************************************************
      SUBROUTINE SDSCL (HMAX,N,NQ,RMAX,H,RC,RH,YH)
      REAL H, HMAX, RC, RH, RMAX, R1, YH(N,*)
      IF (H .LT. 1.E0) THEN
        RH = MIN(ABS(H)*RH, ABS(H)*RMAX, HMAX)/ABS(H)
      ELSE
        RH = MIN(RH, RMAX, HMAX/ABS(H))
      END IF
      R1 = 1.E0
      DO 10 J = 1,NQ
        R1 = R1*RH
        DO 10 I = 1,N
 10       YH(I,J+1) = YH(I,J+1)*R1
      H = H*RH
      RC = RC*RH
      END
c***********************************************************************      
      SUBROUTINE SDSTP (EPS,F,FA,HMAX,IMPL,JACOBN,MATDIM,MAXORD,MINT,
     8   MITER,ML,MU,N,NDE,YWT,UROUND,USERS,AVGH,AVGORD,H,HUSED,JTASK,
     8   MNTOLD,MTROLD,NFE,NJE,NQUSED,NSTEP,T,Y,YH,A,CONVRG,DFDY,EL,FAC,
     8   HOLD,IPVT,JSTATE,NQ,NWAIT,RC,RMAX,SAVE1,SAVE2,TQ,TREND,ISWFLG,
     8   MTRSV,MXRDSV)
C***ROUTINES CALLED  SDNTL,SDPST,SDCOR,SDPSC,SDSCL,SNRM2
      EXTERNAL F, JACOBN, FA, USERS
      REAL A(MATDIM,*), AVGH, AVGORD, BIAS1, BIAS2, BIAS3,
     8     BND, CTEST, D, DENOM, DFDY(MATDIM,*), D1, EL(13,12), EPS,
     8     ERDN, ERUP, ETEST, FAC(*), H, HMAX, HN, HOLD, HS, HUSED,
     8     NUMER, RC, RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL, RMNORM,
     8     SAVE1(*), SAVE2(*), SNRM2, T, TOLD, TQ(3,12), TREND, TRSHLD,
     8     UROUND, Y(*), YH(N,*), YWT(*), Y0NRM
      INTEGER IPVT(*)
      LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH
      PARAMETER(BIAS1 = 1.3E0, BIAS2 = 1.2E0, BIAS3 = 1.4E0, MXFAIL = 3,
     8          MXITER = 3, MXTRY = 50, RCTEST = .3E0, RMFAIL = 2.E0,
     8          RMNORM = 10.E0, TRSHLD = 1.E0)
      DATA IER /.FALSE./
      NSV = N
      BND = 0.E0
      SWITCH = .FALSE.
      NTRY = 0
      TOLD = T
      NFAIL = 0
      IF (JTASK .LE. 0) THEN
        CALL SDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC,
     8              YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH,
     8              RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
        IF (N .EQ. 0) GO TO 440
        IF (H .EQ. 0.E0) GO TO 400
        IF (IER) GO TO 420
      END IF
 100  NTRY = NTRY + 1
      IF (NTRY .GT. MXTRY) GO TO 410
      T = T + H
      CALL SDPSC (1, N, NQ,  YH)
      EVALJC = ((ABS(RC - 1.E0) .GT. RCTEST) .AND. (MITER .NE. 0))
      EVALFA = .NOT. EVALJC
C
 110  ITER = 0
      DO 115 I = 1,N
 115    Y(I) = YH(I,1)
      CALL F (N, T, Y, SAVE2)
      IF (N .EQ. 0) THEN
        JSTATE = 6
        GO TO 430
      END IF
      NFE = NFE + 1
      IF (EVALJC .OR. IER) THEN
        CALL SDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,
     8              MU, N, NDE, NQ, SAVE2, T, USERS, Y, YH, YWT, UROUND,
     8              NFE, NJE,  A, DFDY, FAC, IER, IPVT, SAVE1, ISWFLG,
     8              BND, JSTATE)
        IF (N .EQ. 0) GO TO 430
        IF (IER) GO TO 160
        CONVRG = .FALSE.
        RC = 1.E0
      END IF
      DO 125 I = 1,N
 125    SAVE1(I) = 0.E0
C
 130  CALL SDCOR (DFDY, EL, FA, H, IMPL, IPVT, MATDIM, MITER, ML,
     8            MU, N, NDE, NQ, T, USERS, Y, YH, YWT,  EVALFA, SAVE1,
     8            SAVE2,  A, D, JSTATE)
        IF (N .EQ. 0) GO TO 430
      IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
        IF (ITER .EQ. 0) THEN
          NUMER = SNRM2(N, SAVE1, 1)
          DO 132 I = 1,N
 132        DFDY(1,I) = SAVE1(I)
          Y0NRM = SNRM2(N, YH, 1)
        ELSE
          DENOM = NUMER
          DO 134 I = 1,N
 134        DFDY(1,I) = SAVE1(I) - DFDY(1,I)
          NUMER = SNRM2(N, DFDY, MATDIM)
          IF (EL(1,NQ)*NUMER .LE. 100.E0*UROUND*Y0NRM) THEN
            IF (RMAX .EQ. RMFAIL) THEN
              SWITCH = .TRUE.
              GO TO 170
            END IF
          END IF
          DO 136 I = 1,N
 136        DFDY(1,I) = SAVE1(I)
          IF (DENOM .NE. 0.E0)
     8    BND = MAX(BND, NUMER/(DENOM*ABS(H)*EL(1,NQ)))
        END IF
      END IF
      IF (ITER .GT. 0) TREND = MAX(.9E0*TREND, D/D1)
      D1 = D
      CTEST = MIN(2.E0*TREND, 1.E0)*D
      IF (CTEST .LE. EPS) GO TO 170
      ITER = ITER + 1
      IF (ITER .LT. MXITER) THEN
        DO 140 I = 1,N
 140      Y(I) = YH(I,1) + EL(1,NQ)*SAVE1(I)
        CALL F (N, T, Y, SAVE2)
        IF (N .EQ. 0) THEN
          JSTATE = 6
          GO TO 430
        END IF
        NFE = NFE + 1
        GO TO 130
      END IF
c
      IF (CONVRG) THEN
        EVALJC = .TRUE.
        EVALFA = .FALSE.
        GO TO 110
      END IF
 160  T = TOLD
      CALL SDPSC (-1, N, NQ,  YH)
      NWAIT = NQ + 2
      IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
      IF (ITER .EQ. 0) THEN
        RH = .3E0
      ELSE
        RH = .9E0*(EPS/CTEST)**(.2E0)
      END IF
      IF (RH*H .EQ. 0.E0) GO TO 400
      CALL SDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
      GO TO 100
C
 170  CONVRG = (MITER .NE. 0)
      DO 180 I = 1,NDE
 180    SAVE2(I) = SAVE1(I)/YWT(I)
      ETEST = SNRM2(NDE, SAVE2, 1)/(TQ(2,NQ)*SQRT(REAL(NDE)))
C
      IF (ETEST .GT. EPS) THEN
        T = TOLD
        CALL SDPSC (-1, N, NQ,  YH)
        NFAIL = NFAIL + 1
        IF (NFAIL .LT. MXFAIL) THEN
          IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
          RH2 = 1.E0/(BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))
          IF (NQ .GT. 1) THEN
            DO 190 I = 1,NDE
 190          SAVE2(I) = YH(I,NQ+1)/YWT(I)
            ERDN = SNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(REAL(NDE)))
            RH1 = 1.E0/MAX(1.E0, BIAS1*(ERDN/EPS)**(1.E0/REAL(NQ)))
            IF (RH2 .LT. RH1) THEN
              NQ = NQ - 1
              RC = RC*EL(1,NQ)/EL(1,NQ+1)
              RH = RH1
            ELSE
              RH = RH2
            END IF
          ELSE
            RH = RH2
          END IF
          NWAIT = NQ + 2
          IF (RH*H .EQ. 0.E0) GO TO 400
          CALL SDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
          GO TO 100
        END IF
C
        NFAIL = 0
        JTASK = 2
        DO 215 I = 1,N
 215      Y(I) = YH(I,1)
        CALL SDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC,
     8              YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH,
     8              RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
        RMAX = RMNORM
        IF (N .EQ. 0) GO TO 440
        IF (H .EQ. 0.E0) GO TO 400
        IF (IER) GO TO 420
        GO TO 100
      END IF
C
      NSTEP = NSTEP + 1
      HUSED = H
      NQUSED = NQ
      AVGH = (REAL(NSTEP-1)*AVGH + H)/REAL(NSTEP)
      AVGORD = (REAL(NSTEP-1)*AVGORD + REAL(NQ))/REAL(NSTEP)
      DO 230 J = 1,NQ+1
        DO 230 I = 1,N
 230      YH(I,J) = YH(I,J) + EL(J,NQ)*SAVE1(I)
      DO 235 I = 1,N
 235    Y(I) = YH(I,1)
C                                          If ISWFLG is 3, consider
C                                          changing integration methods.
      IF (ISWFLG .EQ. 3) THEN
        IF (BND .NE. 0.E0) THEN
          IF (MINT .EQ. 1 .AND. NQ .LE. 5) THEN
            HN = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.E0/REAL(NQ+1)))
            HN = MIN(HN, 1.E0/(2.E0*EL(1,NQ)*BND))
            HS = ABS(H)/MAX(UROUND,
     8      (ETEST/(EPS*EL(NQ+1,1)))**(1.E0/REAL(NQ+1)))
            IF (HS .GT. 1.2E0*HN) THEN
              MINT = 2
              MNTOLD = MINT
              MITER = MTRSV
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 5)
              RC = 0.E0
              RMAX = RMNORM
              TREND = 1.E0
              CALL SDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          ELSE IF (MINT .EQ. 2) THEN
            HS = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.E0/REAL(NQ+1)))
            HN = ABS(H)/MAX(UROUND,
     8      (ETEST*EL(NQ+1,1)/EPS)**(1.E0/REAL(NQ+1)))
            HN = MIN(HN, 1.E0/(2.E0*EL(1,NQ)*BND))
            IF (HN .GE. HS) THEN
              MINT = 1
              MNTOLD = MINT
              MITER = 0
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 12)
              RMAX = RMNORM
              TREND = 1.E0
              CONVRG = .FALSE.
              CALL SDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          END IF
        END IF
      END IF
      IF (SWITCH) THEN
        MINT = 2
        MNTOLD = MINT
        MITER = MTRSV
        MTROLD = MITER
        MAXORD = MIN(MXRDSV, 5)
        NQ = MIN(NQ, MAXORD)
        RC = 0.E0
        RMAX = RMNORM
        TREND = 1.E0
        CALL SDCST (MAXORD, MINT, ISWFLG, EL, TQ)
        NWAIT = NQ + 2
      END IF
C
      IF (JTASK .EQ. 0 .OR. JTASK .EQ. 2) THEN
        RH = 1.E0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))
        IF (RH.GT.TRSHLD) CALL SDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
      ELSE IF (NWAIT .GT. 1) THEN
        NWAIT = NWAIT - 1
        IF (NWAIT .EQ. 1 .AND. NQ .LT. MAXORD) THEN
          DO 250 I = 1,NDE
 250        YH(I,MAXORD+1) = SAVE1(I)
        END IF
C
      ELSE
        IF (NQ .EQ. 1) THEN
          RH1 = 0.E0
        ELSE
          DO 270 I = 1,NDE
 270        SAVE2(I) = YH(I,NQ+1)/YWT(I)
          ERDN = SNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(REAL(NDE)))
          RH1 = 1.E0/MAX(UROUND, BIAS1*(ERDN/EPS)**(1.E0/REAL(NQ)))
        END IF
        RH2 = 1.E0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))
        IF (NQ .EQ. MAXORD) THEN
          RH3 = 0.E0
        ELSE
          DO 290 I = 1,NDE
 290        SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/YWT(I)
          ERUP = SNRM2(NDE, SAVE2, 1)/(TQ(3,NQ)*SQRT(REAL(NDE)))
          RH3 = 1.E0/MAX(UROUND, BIAS3*(ERUP/EPS)**(1.E0/REAL(NQ+2)))
        END IF
        IF (RH1 .GT. RH2 .AND. RH1 .GE. RH3) THEN
          RH = RH1
          IF (RH .LE. TRSHLD) GO TO 380
          NQ = NQ - 1
          RC = RC*EL(1,NQ)/EL(1,NQ+1)
        ELSE IF (RH2 .GE. RH1 .AND. RH2 .GE. RH3) THEN
          RH = RH2
          IF (RH .LE. TRSHLD) GO TO 380
        ELSE
          RH = RH3
          IF (RH .LE. TRSHLD) GO TO 380
          DO 360 I = 1,N
 360        YH(I,NQ+2) = SAVE1(I)*EL(NQ+1,NQ)/REAL(NQ+1)
          NQ = NQ + 1
          RC = RC*EL(1,NQ)/EL(1,NQ-1)
        END IF
        IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
          IF (BND.NE.0.E0) RH = MIN(RH, 1.E0/(2.E0*EL(1,NQ)*BND*ABS(H)))
        END IF
        CALL SDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
        RMAX = RMNORM
 380    NWAIT = NQ + 2
      END IF
C
      JSTATE = 1
      HOLD = H
      RETURN
C
 400  JSTATE = 2
      HOLD = H
      DO 405 I = 1,N
 405    Y(I) = YH(I,1)
      RETURN
C
 410  JSTATE = 3
      HOLD = H
      RETURN
C
 420  JSTATE = 4
      HOLD = H
      RETURN
C
 430  T = TOLD
      CALL SDPSC (-1, NSV, NQ,  YH)
      DO 435 I = 1,NSV
 435    Y(I) = YH(I,1)
 440  HOLD = H
      RETURN
      END
c**********************************************************************      
      SUBROUTINE SDZRO (AE,F,H,N,NQ,IROOT,RE,T,YH,UROUND,B,C,FB,FC,Y)
C***ROUTINES CALLED  SDNTP
      REAL A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC,
     8     H, P, Q, RE, RW, T, TOL, UROUND, Y(*), YH(N,*)
      ER = 4.E0*UROUND
      RW = MAX(RE, ER)
      IC = 0
      ACBS = ABS(B - C)
      A = C
      FA = FC
      KOUNT = 0
C
 10   IF (ABS(FC) .LT. ABS(FB)) THEN
        A = B
        FA = FB
        B = C
        FB = FC
        C = A
        FC = FA
      END IF
      CMB = 0.5E0*(C - B)
      ACMB = ABS(CMB)
      TOL = RW*ABS(B) + AE
C
      IF (ACMB .LE. TOL) RETURN
      IF (KOUNT .GT. 50) RETURN
C
      P = (B - A)*FB
      Q = FA - FB
      IF (P .LT. 0.E0) THEN
        P = -P
        Q = -Q
      END IF
C
      A = B
      FA = FB
      IC = IC + 1
      IF (IC .GE. 4) THEN
        IF (8.E0*ACMB .GE. ACBS) THEN
C
          B = 0.5E0*(C + B)
          GO TO 20
        END IF
        IC = 0
      END IF
      ACBS = ACMB
C
      IF (P .LE. ABS(Q)*TOL) THEN
        B = B + SIGN(TOL, CMB)
      ELSE IF (P .LT. CMB*Q) THEN
        B = B + P/Q
      ELSE
        B = 0.5E0*(C + B)
      END IF
 20   CALL SDNTP (H, 0, N, NQ, T, B, YH,  Y)
      FB = F(N, B, Y, IROOT)
      IF (N .EQ. 0) RETURN
      IF (FB .EQ. 0.E0) RETURN
      KOUNT = KOUNT + 1
C
      IF (SIGN(1.0E0, FB) .EQ. SIGN(1.0E0, FC)) THEN
        C = A
        FC = FA
      END IF
      GO TO 10
      END
c***********************************************************************      
      SUBROUTINE EQDIF3 (N,T,Y,F,NSTATE,TOUT,NTASK,NROOT,EPS,EWT,IERROR,
     8   MINT,MITER,IMPL,ML,MU,MXORD,HMAX,WORK,LENW,IWORK,LENIW,JACOBN,
     8   FA,NDE,MXSTEP,G,USERS)
C
      EXTERNAL F, JACOBN, FA, G, USERS
      REAL AE, BIG, EPS, EWT(*), G, GLAST, H, HMAX, HSIGN,
     8	   NROUND, RE, r2mach, SIZE, SNRM2, SUM, T, TLAST, TOUT, TROOT,
     8     UROUND, WORK(*), Y(*)
      INTEGER IWORK(*)
      LOGICAL CONVRG
      CHARACTER MSG*205
      PARAMETER(NROUND = 20.E0)
      PARAMETER(IAVGH = 1, IHUSED = 2, IAVGRD = 3,
     8          IEL = 4, IH = 160, IHMAX = 161, IHOLD = 162,
     8          IHSIGN = 163, IRC = 164, IRMAX = 165, IT = 166,
     8          ITOUT = 167, ITQ = 168, ITREND = 204, IYH = 205,
     8          INDMXR = 1, INQUSD = 2, INSTEP = 3, INFE = 4, INJE = 5,
     8          INROOT = 6, ICNVRG = 7, IJROOT = 8, IJTASK = 9,
     8          IMNTLD = 10, IMTRLD = 11, INQ = 12, INRTLD = 13,
     8          INDTRT = 14, INWAIT = 15, IMNT = 16, IMTRSV = 17,
     8          IMTR = 18, IMXRDS = 19, IMXORD = 20, INDPRT = 21,
     8          INDPVT = 22)
      NPAR = N
      UROUND = r2mach (4)
      IF (NROOT .NE. 0) THEN
	AE = r2mach(1)
        RE = UROUND
      END IF
      IF (EPS .LT. 0.E0) THEN
	WRITE(MSG, '(''EQDIF36FE Illegal input.  EPS,'', E16.8,
     8  '', is negative.'')') EPS
	CALL xerreur(MSG(1:60), 60, 6, 2)
        RETURN
      END IF
      IF (N .LE. 0) THEN
	WRITE(MSG, '(''EQDIF37FE Illegal input.  Number of equations,'',
     8  I8, '', is not positive.'')') N
	CALL xerreur(MSG(1:72), 72, 7, 2)
        RETURN
      END IF
      IF (MXORD .LE. 0) THEN
	WRITE(MSG, '(''EQDIF314FE Illegal input.  Maximum order,'', I8,
     8  '', is not positive.'')') MXORD
	CALL xerreur(MSG(1:67), 67, 14, 2)
        RETURN
      END IF
      IF ((MINT .LT. 1 .OR. MINT .GT. 3) .OR. (MINT .EQ. 3 .AND.
     8  (MITER .EQ. 0 .OR. MITER .EQ. 3 .OR. IMPL .NE. 0))
     8  .OR. (MITER .LT. 0 .OR. MITER .GT. 5) .OR.
     8  (IMPL .NE. 0 .AND. IMPL .NE. 1 .AND. IMPL .NE. 2) .OR.
     8  ((IMPL .EQ. 1 .OR. IMPL .EQ. 2) .AND. MITER .EQ. 0) .OR.
     8  (IMPL .EQ. 2 .AND. MINT .EQ. 1) .OR.
     8  (NSTATE .LT. 1 .OR. NSTATE .GT. 10)) THEN
	WRITE(MSG, '(''EQDIF39FE Illegal input.  Improper value for '',
     8  ''NSTATE(MSTATE), MINT, MITER or IMPL.'')')
	CALL xerreur(MSG(1:81), 81, 9, 2)
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        LIWCHK = INDPVT - 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2 .OR. MITER .EQ. 4 .OR.
     8  MITER .EQ. 5) THEN
        LIWCHK = INDPVT + N - 1
      END IF
      IF (LENIW .LT. LIWCHK) THEN
	WRITE(MSG, '(''EQDIF310FE Illegal input.  Insufficient '',
     8  ''storage allocated for the IWORK array.  Based on the '')')
        WRITE(MSG(94:), '(''value of the input parameters involved, '',
     8  ''the required storage is'', I8)') LIWCHK
	CALL xerreur(MSG(1:164), 164, 10, 2)
        RETURN
      END IF
      IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
        MAXORD = MIN(MXORD, 12)
      ELSE IF (MINT .EQ. 2) THEN
        MAXORD = MIN(MXORD, 5)
      END IF
      IDFDY = IYH + (MAXORD + 1)*N
C
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3)  THEN
        IYWT = IDFDY
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2)  THEN
        IYWT = IDFDY + N*N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5)  THEN
        IYWT = IDFDY + (2*ML + MU + 1)*N
      END IF
      ISAVE1 = IYWT + N
      ISAVE2 = ISAVE1 + N
      IGNOW = ISAVE2 + N
      ITROOT = IGNOW + NROOT
      IFAC = ITROOT + NROOT
      IF (MITER .EQ. 2 .OR. MITER .EQ. 5 .OR. MINT .EQ. 3) THEN
        IA = IFAC + N
      ELSE
        IA = IFAC
      END IF
      IF (IMPL .EQ. 0 .OR. MITER .EQ. 3) THEN
        LENCHK = IA - 1
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
        LENCHK = IA - 1 + N*N
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 4 .OR. MITER .EQ. 5)) THEN
        LENCHK = IA - 1 + (2*ML + MU + 1)*N
      ELSE IF (IMPL .EQ. 2 .AND. MITER .NE. 3) THEN
        LENCHK = IA - 1 + N
      END IF
      IF (LENW .LT. LENCHK) THEN
	WRITE(MSG, '(''EQDIF38FE Illegal input.  Insufficient '',
     8  ''storage allocated for the WORK array.  Based on the '')')
        WRITE(MSG(92:), '(''value of the input parameters involved, '',
     8  ''the required storage is'', I8)') LENCHK
	CALL xerreur(MSG(1:162), 162, 8, 2)
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        MATDIM = 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        MATDIM = N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        MATDIM = 2*ML + MU + 1
      END IF
      IF (IMPL .EQ. 0 .OR. IMPL .EQ. 1) THEN
        NDECOM = N
      ELSE IF (IMPL .EQ. 2) THEN
        NDECOM = NDE
      END IF
      IF (NSTATE .EQ. 1) THEN
C
        IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
          IWORK(IMXORD) = MIN(MXORD, 12)
        ELSE IF (MINT .EQ. 2) THEN
          IWORK(IMXORD) = MIN(MXORD, 5)
        END IF
        IWORK(IMXRDS) = MXORD
        IF (MINT .EQ. 1 .OR. MINT .EQ. 2) THEN
          IWORK(IMNT) = MINT
          IWORK(IMTR) = MITER
          IWORK(IMNTLD) = MINT
          IWORK(IMTRLD) = MITER
        ELSE IF (MINT .EQ. 3) THEN
          IWORK(IMNT) = 1
          IWORK(IMTR) = 0
          IWORK(IMNTLD) = IWORK(IMNT)
          IWORK(IMTRLD) = IWORK(IMTR)
          IWORK(IMTRSV) = MITER
        END IF
        WORK(IHMAX) = HMAX
        H = (TOUT - T)*(1.E0 - 4.E0*UROUND)
        H = SIGN(MIN(ABS(H), HMAX), H)
        WORK(IH) = H
        HSIGN = SIGN(1.E0, H)
        WORK(IHSIGN) = HSIGN
        IWORK(IJTASK) = 0
        WORK(IAVGH) = 0.E0
        WORK(IHUSED) = 0.E0
        WORK(IAVGRD) = 0.E0
        IWORK(INDMXR) = 0
        IWORK(INQUSD) = 0
        IWORK(INSTEP) = 0
        IWORK(INFE) = 0
        IWORK(INJE) = 0
        IWORK(INROOT) = 0
        WORK(IT) = T
        IWORK(ICNVRG) = 0
        IWORK(INDPRT) = 0
C
        DO 30 I = 1,N
          JYH = I + IYH - 1
 30       WORK(JYH) = Y(I)
        IF (T .EQ. TOUT) RETURN
        GO TO 180
      END IF
C
      IF (IWORK(ICNVRG) .EQ. 1) THEN
        CONVRG = .TRUE.
      ELSE
        CONVRG = .FALSE.
      END IF
      T = WORK(IT)
      H = WORK(IH)
      HSIGN = WORK(IHSIGN)
      IF (IWORK(IJTASK) .EQ. 0) GO TO 180
C
      IF (NROOT .NE. 0) THEN
        JROOT = IWORK(IJROOT)
        IF (JROOT .GT. 0) THEN
          IF (NSTATE .NE. 5) THEN
            IF (TOUT*HSIGN .GE. WORK(ITOUT)*HSIGN) THEN
              TROOT = WORK(ITOUT)
              CALL SDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
              T = TROOT
              NSTATE = 5
              GO TO 580
            END IF
          ELSE
            TROOT = T
            IROOT = 0
            DO 50 I = 1,IWORK(INRTLD)
              JTROOT = IWORK(INDTRT) + I - 1
              IF (WORK(JTROOT)*HSIGN .LE. TROOT*HSIGN) THEN
C
                IF (WORK(JTROOT) .EQ. WORK(ITOUT) .AND.
     8          I .GT. IWORK(INROOT)) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                  GO TO 60
                END IF
                IF (WORK(JTROOT)*HSIGN .GT. WORK(ITOUT)*HSIGN) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                END IF
              END IF
 50           CONTINUE
 60         IWORK(INROOT) = IROOT
            WORK(ITOUT) = TROOT
            IWORK(IJROOT) = NTASK
            IF (NTASK .EQ. 1) THEN
              IF (IROOT .EQ. 0) THEN
                IWORK(IJROOT) = 0
              ELSE
                IF (TOUT*HSIGN .GE. TROOT*HSIGN) THEN
                  CALL SDNTP(H, 0, N, IWORK(INQ), T, TROOT,WORK(IYH),Y)
                  NSTATE = 5
                  T = TROOT
                  GO TO 580
                END IF
              END IF
            ELSE IF (NTASK .EQ. 2 .OR. NTASK .EQ. 3) THEN
C
              IF (IROOT .EQ. 0 .OR. (TOUT*HSIGN .LT. TROOT*HSIGN)) THEN
                IWORK(IJROOT) = 0
              ELSE
                CALL SDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH), Y)
                NSTATE = 5
                T = TROOT
                GO TO 580
              END IF
            END IF
          END IF
        END IF
      END IF
C
      IF (NTASK .EQ. 1) THEN
        NSTATE = 2
        IF (T*HSIGN .GE. TOUT*HSIGN) THEN
          CALL SDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
      ELSE IF (NTASK .EQ. 2) THEN
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
	  WRITE(MSG, '(''EQDIF32WRN With NTASK='', I1, '' on input, '',
     8    ''T,'', E16.8, '', was beyond TOUT,'', E16.8, ''.  Solution'',
     8    '' obtained by interpolation.'')') NTASK, T, TOUT
	  CALL xerreur(MSG(1:124), 124, 2, 0)
          CALL SDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          NSTATE = 2
          GO TO 580
        END IF
C
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          NSTATE = 2
          GO TO 560
        END IF
        IF (NSTATE .EQ. 5) THEN
          NSTATE = 2
          GO TO 560
        END IF
        NSTATE = 2
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.E0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        NSTATE = 2
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
	  WRITE(MSG, '(''EQDIF32WRN With NTASK='', I1, '' on input, '',
     8    ''T,'', E16.8, '', was beyond TOUT,'', E16.8, ''.  Solution'',
     8    '' obtained by interpolation.'')') NTASK, T, TOUT
	  CALL xerreur(MSG(1:124), 124, 2, 0)
          CALL SDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          GO TO 560
        END IF
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.E0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      END IF
C
      IF ((MINT .NE. IWORK(IMNTLD) .OR. MITER .NE. IWORK(IMTRLD)) .AND.
     8  MINT .NE. 3 .AND. IWORK(IMNTLD) .NE. 3) IWORK(IJTASK) = -1
      IF (HMAX .NE. WORK(IHMAX)) THEN
        H = SIGN(MIN(ABS(H), HMAX), H)
        IF (H .NE. WORK(IH)) THEN
          IWORK(IJTASK) = -1
          WORK(IH) = H
        END IF
        WORK(IHMAX) = HMAX
      END IF
C
 180  NSTEPL = IWORK(INSTEP)
      DO 190 I = 1,N
        JYH = IYH + I - 1
 190    Y(I) = WORK(JYH)
      IF (NROOT .NE. 0) THEN
        DO 200 I = 1,NROOT
          JGNOW = IGNOW + I - 1
          WORK(JGNOW) = G (NPAR, T, Y, I)
          IF (NPAR .EQ. 0) THEN
            IWORK(INROOT) = I
            NSTATE = 7
            RETURN
          END IF
 200     CONTINUE
      END IF
      IF (IERROR .EQ. 1) THEN
        DO 230 I = 1,N
          JYWT = I + IYWT - 1
 230      WORK(JYWT) = 1.E0
        GO TO 410
      ELSE IF (IERROR .EQ. 5) THEN
        DO 250 I = 1,N
          JYWT = I + IYWT - 1
 250      WORK(JYWT) = EWT(I)
        GO TO 410
      END IF
C                                       Reset YWT array.  Looping point.
 260  IF (IERROR .EQ. 2) THEN
        DO 280 I = 1,N
          IF (Y(I) .EQ. 0.E0) GO TO 290
          JYWT = I + IYWT - 1
 280      WORK(JYWT) = ABS(Y(I))
        GO TO 410
 290    IF (IWORK(IJTASK) .EQ. 0) THEN
          CALL F (NPAR, T, Y, WORK(ISAVE2))
          IF (NPAR .EQ. 0) THEN
            NSTATE = 6
            RETURN
          END IF
          IWORK(INFE) = IWORK(INFE) + 1
          IF (MITER .EQ. 3 .AND. IMPL .NE. 0) THEN
            IFLAG = 0
            CALL USERS(Y, WORK(IYH), WORK(IYWT), WORK(ISAVE1),
     8                 WORK(ISAVE2), T, H, WORK(IEL), IMPL, NPAR,
     8                 NDECOM, IFLAG)
            IF (NPAR .EQ. 0) THEN
              NSTATE = 10
              RETURN
            END IF
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL SGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL SGESL(WORK(IA),MATDIM,N,IWORK(INDPVT),WORK(ISAVE2),0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              JAML = IA + ML
              CALL FA (NPAR, T, Y, WORK(JAML), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL SGBFA (WORK(IA),MATDIM,N,ML,MU,IWORK(INDPVT),INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL SGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
            IF (NPAR .EQ. 0) THEN
              NSTATE = 9
              RETURN
            END IF
            DO 340 I = 1,NDECOM
              JA = I + IA - 1
              JSAVE2 = I + ISAVE2 - 1
              IF (WORK(JA) .EQ. 0.E0) GO TO 690
 340          WORK(JSAVE2) = WORK(JSAVE2)/WORK(JA)
          END IF
        END IF
        DO 360 J = I,N
          JYWT = J + IYWT - 1
          IF (Y(J) .NE. 0.E0) THEN
            WORK(JYWT) = ABS(Y(J))
          ELSE
            IF (IWORK(IJTASK) .EQ. 0) THEN
              JSAVE2 = J + ISAVE2 - 1
              WORK(JYWT) = ABS(H*WORK(JSAVE2))
            ELSE
              JHYP = J + IYH + N - 1
              WORK(JYWT) = ABS(WORK(JHYP))
            END IF
          END IF
          IF (WORK(JYWT) .EQ. 0.E0) WORK(JYWT) = UROUND
 360      CONTINUE
      ELSE IF (IERROR .EQ. 3) THEN
        DO 380 I = 1,N
          JYWT = I + IYWT - 1
 380      WORK(JYWT) = MAX(EWT(1), ABS(Y(I)))
      ELSE IF (IERROR .EQ. 4) THEN
        DO 400 I = 1,N
          JYWT = I + IYWT - 1
 400      WORK(JYWT) = MAX(EWT(I), ABS(Y(I)))
      END IF
C
 410  DO 420 I = 1,N
        JYWT = I + IYWT - 1
        JSAVE2 = I + ISAVE2 - 1
 420    WORK(JSAVE2) = Y(I)/WORK(JYWT)
      SUM = SNRM2(N, WORK(ISAVE2), 1)/SQRT(REAL(N))
      IF (EPS .LT. SUM*UROUND) THEN
        EPS = SUM*UROUND*(1.E0 + 10.E0*UROUND)
	WRITE(MSG, '(''EQDIF34REC At T,'', E16.8, '', the requested '',
     8  ''accuracy, EPS, was not obtainable with the machine '',
     8  ''precision.  EPS has been increased to'')') T
        WRITE(MSG(137:), '(E16.8)') EPS
	CALL xerreur(MSG(1:152), 152, 4, 1)
        NSTATE = 4
        GO TO 560
      END IF
      IF (ABS(H) .GE. UROUND*ABS(T)) THEN
        IWORK(INDPRT) = 0
      ELSE IF (IWORK(INDPRT) .EQ. 0) THEN
c	WRITE(MSG, '(''EQDIF35WRN At T,'', E16.8, '', the step size,'',
c    8  E16.8, '', is smaller than the roundoff level of T.  '')') T, H
c        WRITE(MSG(109:), '(''This may occur if there is an abrupt '',
c     8  ''change in the right hand side of the differential '',
c     8  ''equations.'')')
c	CALL xerreur(MSG(1:205), 205, 5, 0)
        IWORK(INDPRT) = 1
      END IF
      IF (NTASK.NE.2) THEN
        IF ((IWORK(INSTEP)-NSTEPL) .GT. MXSTEP) THEN
	  WRITE(MSG, '(''EQDIF33WRN At T,'', E16.8, '', '', I8,
     8    '' steps have been taken without reaching TOUT,'', E16.8)')
     8    T, MXSTEP, TOUT
	  CALL xerreur(MSG(1:103), 103, 3, 0)
          NSTATE = 3
          GO TO 560
        END IF
      END IF
C
      CALL SDSTP (EPS, F, FA, WORK(IHMAX), IMPL, JACOBN, MATDIM,
     8            IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML, MU, NPAR,
     8           NDECOM, WORK(IYWT), UROUND, USERS,  WORK(IAVGH),
     8           WORK(IAVGRD), WORK(IH), WORK(IHUSED), IWORK(IJTASK),
     8           IWORK(IMNTLD), IWORK(IMTRLD), IWORK(INFE), IWORK(INJE),
     8            IWORK(INQUSD), IWORK(INSTEP), WORK(IT), Y, WORK(IYH),
     8            WORK(IA), CONVRG, WORK(IDFDY), WORK(IEL), WORK(IFAC),
     8            WORK(IHOLD), IWORK(INDPVT), JSTATE, IWORK(INQ),
     8            IWORK(INWAIT), WORK(IRC), WORK(IRMAX), WORK(ISAVE1),
     8            WORK(ISAVE2), WORK(ITQ), WORK(ITREND), MINT,
     8            IWORK(IMTRSV), IWORK(IMXRDS))
      T = WORK(IT)
      H = WORK(IH)
      GO TO (470, 670, 680, 690, 690, 660, 660, 660, 660, 660), JSTATE
 470  IWORK(IJTASK) = 1
      IF (NROOT .NE. 0) THEN
        IROOT = 0
        DO 500 I = 1,NROOT
          JTROOT = ITROOT + I - 1
          JGNOW = IGNOW + I - 1
          GLAST = WORK(JGNOW)
          WORK(JGNOW) = G (NPAR, T, Y, I)
          IF (NPAR .EQ. 0) THEN
            IWORK(INROOT) = I
            NSTATE = 7
            RETURN
          END IF
          IF (GLAST*WORK(JGNOW) .GT. 0.E0) THEN
            WORK(JTROOT) = T + H
          ELSE
            IF (WORK(JGNOW) .EQ. 0.E0) THEN
              WORK(JTROOT) = T
              IROOT = I
            ELSE
              IF (GLAST .EQ. 0.E0) THEN
                WORK(JTROOT) = T + H
              ELSE
                IF (ABS(WORK(IHUSED)) .GE. UROUND*ABS(T)) THEN
                  TLAST = T - WORK(IHUSED)
                  IROOT = I
                  TROOT = T
                  CALL SDZRO (AE, G, H, NPAR, IWORK(INQ), IROOT, RE, T,
     8                        WORK(IYH), UROUND,  TROOT, TLAST,
     8                        WORK(JGNOW), GLAST,  Y)
                  DO 480 J = 1,N
  480               Y(J) = WORK(IYH + J -1)             
                  IF (NPAR .EQ. 0) THEN
                    IWORK(INROOT) = I
                    NSTATE = 7
                    RETURN
                  END IF
                  WORK(JTROOT) = TROOT
                ELSE
                  WORK(JTROOT) = T
                  IROOT = I
                END IF
              END IF
            END IF
          END IF
 500      CONTINUE
        IF (IROOT .EQ. 0) THEN
          IWORK(IJROOT) = 0
        ELSE
          IWORK(IJROOT) = NTASK
          IWORK(INRTLD) = NROOT
          IWORK(INDTRT) = ITROOT
          TROOT = T + H
          DO 510 I = 1,NROOT
            JTROOT = ITROOT + I - 1
            IF (WORK(JTROOT)*HSIGN .LT. TROOT*HSIGN) THEN
              TROOT = WORK(JTROOT)
              IROOT = I
            END IF
 510        CONTINUE
          IWORK(INROOT) = IROOT
          WORK(ITOUT) = TROOT
          IF (TROOT*HSIGN .LE. TOUT*HSIGN) THEN
            CALL SDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
            NSTATE = 5
            T = TROOT
            GO TO 580
          END IF
        END IF
      END IF
C
      NSTATE = 2
      IF (NTASK .EQ. 1) THEN
        IF (T*HSIGN .LT. TOUT*HSIGN) GO TO 260
        CALL SDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
        T = TOUT
        GO TO 580
C
      ELSE IF (NTASK .EQ. 2) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.E0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.E0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
          GO TO 260
        END IF
      END IF
C
 560  DO 570 I = 1,N
        JYH = I + IYH - 1
 570    Y(I) = WORK(JYH)
 580  IF (CONVRG) THEN
        IWORK(ICNVRG) = 1
      ELSE
        IWORK(ICNVRG) = 0
      END IF
      IF (IWORK(IJTASK) .EQ. 0) RETURN
      BIG = 0.E0
      IMXERR = 1
      IWORK(INDMXR) = IMXERR
      DO  590 I = 1,N
C                                            SIZE = ABS(ERROR(I)/YWT(I))
        JYWT = I + IYWT - 1
        JERROR = I + ISAVE1 - 1
        SIZE = ABS(WORK(JERROR)/WORK(JYWT))
        IF (BIG .LT. SIZE) THEN
          BIG = SIZE
          IMXERR = I
          IWORK(INDMXR) = IMXERR
        END IF
 590    CONTINUE
      RETURN
C
 660  NSTATE = JSTATE
      RETURN
C
 670  WRITE(MSG, '(''EQDIF311FE At T,'', E16.8, '', the attempted '',
     8  ''step size has gone to zero.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')') T
      CALL xerreur(MSG(1:129), 129, 11, 2)
      RETURN
C
 680  WRITE(MSG, '(''EQDIF312FE At T,'', E16.8, '', the step size has'',
     8  '' been reduced about 50 times without advancing the '')') T
      WRITE(MSG(103:), '(''solution.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')')
      CALL xerreur(MSG(1:165), 165, 12, 2)
      RETURN
C
 690  WRITE(MSG, '(''EQDIF313FE At T,'', E16.8, '', while solving'',
     8  '' A*YDOT = F, A is singular.'')') T
      CALL xerreur(MSG(1:74), 74, 13, 2)
      RETURN
      END
c***********************************************************************      
      SUBROUTINE SGBFA(ABD,LDA,N,ML,MU,IPVT,INFO) 
C
      INTEGER LDA,N,ML,MU,IPVT(*),INFO
      REAL ABD(LDA,*)
C
      REAL T
      INTEGER I,ISAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1 
C
      M = ML + MU + 1
      INFO = 0
C
      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0E0 
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
C
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0E0
   40       CONTINUE
   50    CONTINUE
C
         LM = MIN0(ML,N-K)
         L = ISAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
C
         IF (ABD(L,K) .EQ. 0.0E0) GO TO 100
C
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
C
            T = -1.0E0/ABD(M,K)
            CALL SSCAL(LM,T,ABD(M+1,K),1)
C
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU 
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL SAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0E0) INFO = N 
      RETURN
      END
c***********************************************************       
      SUBROUTINE SGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
C
      INTEGER LDA,N,ML,MU,IPVT(*),JOB
      REAL ABD(LDA,*),B(*)
C
      REAL SDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN0(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K) 
                  B(K) = T
   10          CONTINUE
               CALL SAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL SAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = SDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1 
               K = N - KB
               LM = MIN0(ML,N-K)
               B(K) = B(K) + SDOT(LM,ABD(M+1,K),1,B(K+1),1) 
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K) 
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
C****************************************************************
      SUBROUTINE SGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(*),INFO
      REAL A(LDA,*) 
C
      REAL T
      INTEGER ISAMAX,J,K,KP1,L,NM1
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
         L = ISAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
         IF (A(L,K) .EQ. 0.0E0) GO TO 40
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
            T = -1.0E0/A(K,K) 
            CALL SSCAL(N-K,T,A(K+1,K),1)
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL SAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
c***********************************************************************       
      SUBROUTINE SGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(*),JOB
      REAL A(LDA,*),B(*)
C
      REAL SDOT,T
      INTEGER K,KB,L,NM1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL SAXPY(N-K,T,A(K+1,K),1,B(K+1),1) 
   20    CONTINUE
   30    CONTINUE
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL SAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
         DO 60 K = 1, N
            T = SDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + SDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
