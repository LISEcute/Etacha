//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
 void INTG[F][BOUND][INF][EPSABS][EPSREL][RESULT][ABSERR][NEVAL][
             IER][LIMIT][LENW][LAST][IWORK][WORK]; {
//c*********************************************************************
//c    Routine d'int‚gration
//c*********************************************************************
        N=1;

L20:  M = ((N)%(6));
       if( M  ==  0 ) goto L40 ;
      for(/*L30*/ I=1; I<=M; I++) {
        SASUM = SASUM + fabs(SX[I]);
L30:  } //continue;
       if( N  <  6 ) return;
L40:  MP1 = M + 1;
      for(/*L50*/ I=MP1; I<=N; I+6) {
        SASUM = SASUM + fabs(SX[I]) + fabs(SX[I + 1]) + ABS[SX[I + 2]]
            + fabs(SX[I + 3]) + fabs(SX[I + 4]) + ABS[SX[I + 5]];
L50:  } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SAXPY(N,SA,SX,INCX,SY,INCY) {

      double SX[*],SY[*],SA;
       if(N <= 0 || SA == 0.E0) return;
       if(INCX == INCY)  if(INCX-1) 5,20,60;
L5:   } //continue;

      IX = 1;
      IY = 1;
       if(INCX < 0)IX = (-N+1)*INCX + 1 ;
       if(INCY < 0)IY = (-N+1)*INCY + 1 ;
      for(/*L10*/ I=1; I<=N; I++) {
        SY[IY] = SY[IY] + SA*SX[IX];
        IX = IX + INCX;
        IY = IY + INCY;
L10:  } //continue;
      return;

L20:  M = ((N)%(4));
       if( M  ==  0 ) goto L40 ;
      }
//c*********************************************************************
