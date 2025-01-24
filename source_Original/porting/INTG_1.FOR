	 SUBROUTINE INTG(F,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
     1   IER,LIMIT,LENW,LAST,IWORK,WORK)
cc*********************************************************************
cc	Routine d'int‚gration
cc*********************************************************************
        N=1
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
      END
cc*********************************************************************
