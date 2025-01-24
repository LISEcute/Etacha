//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
 void INTG[F][BOUND][INF][EPSABS][EPSREL][RESULT][ABSERR][NEVAL][
             IER][LIMIT][LENW][LAST][IWORK][WORK]; {
//c*********************************************************************
//c    Routine d'int‚gration
//c*********************************************************************
      double ABSERR,BOUND,EPSABS,EPSREL,F,RESULT,WORK;
      int IER,INF,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL;

      double IWORK[LIMIT],WORK[LENW];

      EXTERNAL F;

      IER = 6;
      NEVAL = 0;
      LAST = 0;
      RESULT = 0.0E+00;
      ABSERR = 0.0E+00;
       if(LIMIT < 1 || LENW < LIMIT*4) goto L10;

      L1 = LIMIT+1;
      L2 = LIMIT+L1;
      L3 = LIMIT+L2;

      CALL INTGE[F][BOUND][INF][EPSABS][EPSREL][LIMIT][RESULT][ABSERR][
            NEVAL][IER][WORK[1]][WORK[L1]][WORK[L2]][WORK[L3]][IWORK][LAST];

      LVL = 0;
L10:   if(IER == 6) LVL = 1;
       if(IER != 0) CALL xerreur[ 'ABNORMAL return FROM  INTG'][
            26][IER][LVL];
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void EA(NEWFLG,SVALUE,LIMEXP,RESULT,ABSERR,EPSTAB,IERR) {
      double ABSERR,DELTA1,DELTA2,DELTA3,EPRN,EPSTAB[*],
             ERROR,ERR1,ERR2,ERR3,E0,E1,E2,E3,RELPR,RES,RESULT,
                 RES3LA[3],r2mach,SS,SVALUE,TOL1,TOL2,TOL3;
      int I,IB,IB2,IE,IERR,IN,K1,K2,K3,LIMEXP,N,NEWELM,NUM,NRES;
      LOGICAL NEWFLG;

       if(LIMEXP < 3) {
        IERR = 1;
      CALL xerreur['LIMEXP IS LESS THAN 3'][21][1][1];
        goto L110;
      }
      IERR = 0;
      RES3LA[1]=EPSTAB[LIMEXP+5];
      RES3LA[2]=EPSTAB[LIMEXP+6];
      RES3LA[3]=EPSTAB[LIMEXP+7];
      RESULT=SVALUE;
       if(NEWFLG) {
        N=1;
        NRES=0;
        NEWFLG=.FALSE.;
        EPSTAB[N]=SVALUE;
        ABSERR=fabs(RESULT);
        goto L100;
      }
 else {
        N=INT[EPSTAB[LIMEXP+3]];
        NRES=INT[EPSTAB[LIMEXP+4]];
         if(N == 2) {
          EPSTAB[N]=SVALUE;
          ABSERR=.6E+01*fabs(RESULT-EPSTAB[1]);
          goto L100;
        }
      }
      EPSTAB[N]=SVALUE;
      RELPR=r2mach[4];
      EPRN=1.0E+01*RELPR;
      EPSTAB[N+2]=EPSTAB[N];
      NEWELM=(N-1)/2;
      NUM=N;
      K1=N;
      for(/*L40*/ I=1; I<=NEWELM; I++) {
        K2=K1-1;
        K3=K1-2;
        RES=EPSTAB[K1+2];
        E0=EPSTAB[K3];
        E1=EPSTAB[K2];
        E2=RES;
        DELTA2=E2-E1;
        ERR2=fabs(DELTA2);
        TOL2=MAX[fabs(E2)][ABS[E1]]*RELPR;
        DELTA3=E1-E0;
        ERR3=fabs(DELTA3);
        TOL3=MAX[fabs(E1)][ABS[E0]]*RELPR;
         if(ERR2 > TOL2 || ERR3 > TOL3) goto L10;

        RESULT=RES;
        ABSERR=ERR2+ERR3;
        goto L50;
L10:     if(I != 1) {
          E3=EPSTAB[K1];
          EPSTAB[K1]=E1;
          DELTA1=E1-E3;
          ERR1=fabs(DELTA1);
          TOL1=MAX[fabs(E1)][ABS[E3]]*RELPR;

           if(ERR1 <= TOL1 || ERR2 <= TOL2 || ERR3 <= TOL3) goto L20;
          SS=0.1E+01/DELTA1+0.1E+01/DELTA2-0.1E+01/DELTA3;
        }
 else {
          EPSTAB[K1]=E1;
           if(ERR2 <= TOL2 || ERR3 <= TOL3) goto L20;
          SS=0.1E+01/DELTA2-0.1E+01/DELTA3;
        }

         if(fabs(SS*E1) > 0.1E-03) goto L30;
L20:    N=I+I-1;
         if(NRES == 0) {
          ABSERR=ERR2+ERR3;
          RESULT=RES;
        }
 else if(NRES == 1) {
          RESULT=RES3LA[1];
        }
 else if(NRES == 2) {
          RESULT=RES3LA[2];
        }
 else {
          RESULT=RES3LA[3];
        }
        goto L50;

L30:    RES=E1+0.1E+01/SS;
        EPSTAB[K1]=RES;
        K1=K1-2;
         if(NRES == 0) {
          ABSERR=ERR2+fabs(RES-E2)+ERR3;
          RESULT=RES;
          goto L40;
        }
 else if(NRES == 1) {
          ERROR=.6E+01*(fabs(RES-RES3LA[1]));
        }
 else if(NRES == 2) {
          ERROR=.2E+01*(fabs(RES-RES3LA[2])+ABS[RES-RES3LA[1]]);
        }
 else {
          ERROR=fabs(RES-RES3LA[3])+ABS[RES-RES3LA[2]]
                    +fabs(RES-RES3LA[1]);
        }
         if(ERROR > 1.0E+01*ABSERR) goto L40;
        ABSERR=ERROR;
        RESULT=RES;
L40:  } //continue;

         if(NRES == 1) {
          ABSERR=.6E+01*(fabs(RESULT-RES3LA[1]));
        }
 else if(NRES == 2) {
          ABSERR=.2E+01*fabs(RESULT-RES3LA[2])+ABS[RESULT-RES3LA[1]];
        }
 else if(NRES > 2) {
          ABSERR=fabs(RESULT-RES3LA[3])+ABS[RESULT-RES3LA[2]]
                    +fabs(RESULT-RES3LA[1]);
        }

L50:   if(N == LIMEXP) N=2*(LIMEXP/2)-1;
      IB=1;
       if((NUM/2)*2 == NUM) IB=2;
      IE=NEWELM+1;
      for(/*L60*/ I=1; I<=IE; I++) {
        IB2=IB+2;
        EPSTAB[IB]=EPSTAB[IB2];
        IB=IB2;
L60:  } //continue;
       if(NUM == N) goto L80;
      IN=NUM-N+1;
      for(/*L70*/ I=1; I<=N; I++) {
        EPSTAB[I]=EPSTAB[IN];
        IN=IN+1;
L70:  } //continue;

L80:   if(NRES == 0) {
        RES3LA[1]=RESULT;
      }
 else if(NRES == 1) {
        RES3LA[2]=RESULT;
      }
 else if(NRES == 2) {
        RES3LA[3]=RESULT;
      }
 else {
        RES3LA[1]=RES3LA[2];
        RES3LA[2]=RES3LA[3];
        RES3LA[3]=RESULT;
      }
L90:  ABSERR=MAX[ABSERR][EPRN*fabs(RESULT)];
      NRES=NRES+1;
L100: N=N+1;
*      if(N <= 3) ABSERR = r2mach[2] * (0.1D-03);
      EPSTAB[LIMEXP+3]=double[N];
      EPSTAB[LIMEXP+4]=double[NRES];
      EPSTAB[LIMEXP+5]=RES3LA[1];
      EPSTAB[LIMEXP+6]=RES3LA[2];
      EPSTAB[LIMEXP+7]=RES3LA[3];
L110: return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void INTGE(F,BOUND,INF,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,
             NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST); {

      double ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,
            A2,BLIST,BOUN,BOUND,B1,B2,CORREC,DEFABS,DEFAB1,DEFAB2,
            DRES,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND,
            ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,ERTEST,F,RESABS,
                RESEPS,RESULT,RLIST,RLIST2,r2mach,SMALL,UFLOW;
      int ID,IER,IERR,IERRO,INF,IORD,IROFF1,IROFF2,IROFF3,JUPBND,K,
            KSGN,KTMIN,LAST,LIMEXP,LIMIT,MAXERR,NEVAL,NRMAX;
      LOGICAL EXTRAP,LERR,NEWFLG,NOEXT;

      PARAMETER [LIMEXP = 50];

      double ALIST[LIMIT],BLIST[LIMIT],ELIST[LIMIT],IORD[LIMIT],
            RLIST[LIMIT],RLIST2[LIMEXP+7];

      EXTERNAL F;

       EPMACH = r2mach[4];

      IER = 0;
      NEVAL = 0;
      LAST = 0;
      RESULT = 0.0E+00;
      ABSERR = 0.0E+00;
      ALIST[1] = 0.0E+00;
      BLIST[1] = 0.1E+01;
      RLIST[1] = 0.0E+00;
      ELIST[1] = 0.0E+00;
      IORD[1] = 0;
      NEWFLG = .TRUE.;
       if(EPSABS <= 0.0E+00 && EPSREL < max(0.5E+02*EPMACH,0.5E-14))
            IER = 6;
       if(IER == 6) goto L999;

      BOUN = BOUND;
       if(INF == 2) BOUN = 0.0E+00;
      CALL QK15I[F][BOUN][INF][0.0E+00][0.1E+01][RESULT][ABSERR][
            DEFABS][RESABS];

      LAST = 1;
      RLIST[1] = RESULT;
      ELIST[1] = ABSERR;
      IORD[1] = 1;
      DRES = fabs(RESULT);
      ERRBND = max(EPSABS,EPSREL*DRES);
       if(ABSERR <= 1.0E+02*EPMACH*DEFABS && ABSERR > 
            ERRBND) IER = 2;
       if(LIMIT == 1) IER = 1;
       if(IER != 0 || [ABSERR <= ERRBND && ABSERR != RESABS] || 
            ABSERR == 0.0E+00) goto L130;

      UFLOW = r2mach[1] * 0.2E+01;
      LERR = .FALSE.;
      CALL EA[NEWFLG][RESULT][LIMEXP][RESEPS][ABSEPS][RLIST2][IERR];
      ERRMAX = ABSERR;
      MAXERR = 1;
      AREA = RESULT;
      ERRSUM = ABSERR;
      NRMAX = 1;
      KTMIN = 0;
      EXTRAP = .FALSE.;
      NOEXT = .FALSE.;
      IERRO = 0;
      IROFF1 = 0;
      IROFF2 = 0;
      IROFF3 = 0;
      KSGN = -1;
       if(DRES >= (0.1E+01-0.5E+02*EPMACH)*DEFABS) KSGN = 1;

      for(/*L90*/ LAST=2; LAST<=LIMIT; LAST++) {

        A1 = ALIST[MAXERR];
        B1 = 0.5E+00*(ALIST[MAXERR]+BLIST[MAXERR]);
        A2 = B1;
        B2 = BLIST[MAXERR];
        ERLAST = ERRMAX;
        CALL QK15I[F][BOUN][INF][A1][B1][AREA1][ERROR1][RESABS][DEFAB1];
        CALL QK15I[F][BOUN][INF][A2][B2][AREA2][ERROR2][RESABS][DEFAB2];

        AREA12 = AREA1+AREA2;
        ERRO12 = ERROR1+ERROR2;
        ERRSUM = ERRSUM+ERRO12-ERRMAX;
        AREA = AREA+AREA12-RLIST[MAXERR];
         if(DEFAB1 == ERROR1 || DEFAB2 == ERROR2)goto L15;
         if(fabs(RLIST[MAXERR]-AREA12) > 0.1E-04*ABS[AREA12]
             || ERRO12 < 0.99E+00*ERRMAX) goto L10;
         if(EXTRAP) IROFF2 = IROFF2+1;
         if(.NOT.EXTRAP) IROFF1 = IROFF1+1;
L10:     if(LAST > 10 && ERRO12 > ERRMAX) IROFF3 = IROFF3+1;
L15:    RLIST[MAXERR] = AREA1;
        RLIST[LAST] = AREA2;
        ERRBND = max(EPSABS,EPSREL*fabs(AREA));

         if(IROFF1+IROFF2 >= 10 || IROFF3 >= 20) IER = 2;
         if(IROFF2 >= 5) IERRO = 3;

         if(LAST == LIMIT) IER = 1;

         if(max(fabs(A1),ABS[B2]) <= (0.1E+01+0.1E+03*EPMACH)*
            [fabs(A2)+0.1E+04*UFLOW]) IER = 4;

         if(ERROR2 > ERROR1) goto L20;
        ALIST[LAST] = A2;
        BLIST[MAXERR] = B1;
        BLIST[LAST] = B2;
        ELIST[MAXERR] = ERROR1;
        ELIST[LAST] = ERROR2;
        goto L30;
L20:    ALIST[MAXERR] = A2;
        ALIST[LAST] = A1;
        BLIST[LAST] = B1;
        RLIST[MAXERR] = AREA2;
        RLIST[LAST] = AREA1;
        ELIST[MAXERR] = ERROR2;
        ELIST[LAST] = ERROR1;

L30:    CALL QPSRT[LIMIT][LAST][MAXERR][ERRMAX][ELIST][IORD][NRMAX];
         if(ERRSUM <= ERRBND) goto L115;
         if(IER != 0) goto L100;
         if(LAST == 2) goto L80;
         if(NOEXT) goto L90;
        ERLARG = ERLARG-ERLAST;
         if(fabs(B1-A1) > SMALL) ERLARG = ERLARG+ERRO12;
         if(EXTRAP) goto L40;

         if(fabs(BLIST[MAXERR]-ALIST[MAXERR]) > SMALL) goto L90;
        EXTRAP = .TRUE.;
        NRMAX = 2;
L40:     if(IERRO == 3 || ERLARG <= ERTEST) goto L60;

        ID = NRMAX;
        JUPBND = LAST;
         if(LAST > [2+LIMIT/2]) JUPBND = LIMIT+3-LAST;
        for(/*L50*/ K=ID; K<=JUPBND; K++) {
          MAXERR = IORD[NRMAX];
          ERRMAX = ELIST[MAXERR];
           if(fabs(BLIST[MAXERR]-ALIST[MAXERR]) > SMALL) goto L90;
          NRMAX = NRMAX+1;
L50:    } //continue;

L60:    CALL EA[NEWFLG][AREA][LIMEXP][RESEPS][ABSEPS][RLIST2][IERR];
        KTMIN = KTMIN+1;
         if((KTMIN > 5) && [ABSERR < 0.1E-02*ERRSUM] && [LERR])
               IER = 5;
         if((ABSEPS >= ABSERR) && [LERR]) goto L70;
        KTMIN = 0;
        ABSERR = ABSEPS;
        LERR = .TRUE.;
        RESULT = RESEPS;
        CORREC = ERLARG;
        ERTEST = max(EPSABS,EPSREL*fabs(RESEPS));
         if((ABSERR <= ERTEST) && [LERR]) goto L100;

L70:     if(RLIST2[LIMEXP+3] == 1) NOEXT = .TRUE.;
         if(IER == 5) goto L100;
        MAXERR = IORD[1];
        ERRMAX = ELIST[MAXERR];
        NRMAX = 1;
        EXTRAP = .FALSE.;
        SMALL = SMALL*0.5E+00;
        ERLARG = ERRSUM;
        goto L90;
L80:    SMALL = 0.375E+00;
        ERLARG = ERRSUM;
        ERTEST = ERRBND;
        CALL EA[NEWFLG][AREA][LIMEXP][RESEPS][ABSEPS][RLIST2][IERR];
L90:  } //continue;

L100:  if(.NOT.LERR) goto L115;
       if((IER+IERRO) == 0) goto L110;
       if(IERRO == 3) ABSERR = ABSERR+CORREC;
       if(IER == 0) IER = 3;
       if(RESULT != 0.0E+00 && AREA != 0.0E+00)goto L105;
       if(ABSERR > ERRSUM)goto L115;
       if(AREA == 0.0E+00) goto L130;
      goto L110;
L105:  if(ABSERR/fabs(RESULT) > ERRSUM/ABS[AREA])goto L115;

L110:  if(KSGN == (-1) && max(fabs(RESULT),ABS[AREA]) <= 
           DEFABS*0.1E-01) goto L130;
       if(0.1E-01 > [RESULT/AREA] || [RESULT/AREA] > 0.1E+03.
          OR.ERRSUM > fabs(AREA)) IER = 6;
      goto L130;

L115: RESULT = 0.0E+00;
      for(/*L120*/ K=1; K<=LAST; K++) {
        RESULT = RESULT+RLIST[K];
L120: } //continue;
      ABSERR = ERRSUM;
L130: NEVAL = 30*LAST-15;
       if(INF == 2) NEVAL = 2*NEVAL;
       if(IER > 2) IER=IER-1;
L999: return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void QK15I(F,BOUN,INF,A,B,RESULT,ABSERR,RESABS,RESASC) {

      double A,ABSC,ABSC1,ABSC2,ABSERR,B,BOUN,CENTR,
                DINF,r2mach,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,
            FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,TABSC1,TABSC2,
            UFLOW,WG,WGK,XGK;
      int INF,J,MIN0;
      EXTERNAL F;

      double FV1[7],FV2[7],XGK[8],WGK[8],WG[8];

      double /*data*/ XGK[1],XGK[2],XGK[3],XGK[4],XGK[5],XGK[6],XGK[7],
            XGK[8]/
               0.9914553711208126E+00,     0.9491079123427585E+00,
               0.8648644233597691E+00,     0.7415311855993944E+00,
               0.5860872354676911E+00,     0.4058451513773972E+00,
               0.2077849550078985E+00,     0.0000000000000000E+00/;

      double /*data*/ WGK[1],WGK[2],WGK[3],WGK[4],WGK[5],WGK[6],WGK[7],
            WGK[8]/
               0.2293532201052922E-01,     0.6309209262997855E-01,
               0.1047900103222502E+00,     0.1406532597155259E+00,
               0.1690047266392679E+00,     0.1903505780647854E+00,
               0.2044329400752989E+00,     0.2094821410847278E+00/;

      double /*data*/ WG[1],WG[2],WG[3],WG[4],WG[5],WG[6],WG[7],WG[8]/
               0.0000000000000000E+00,     0.1294849661688697E+00,
               0.0000000000000000E+00,     0.2797053914892767E+00,
               0.0000000000000000E+00,     0.3818300505051189E+00,
               0.0000000000000000E+00,     0.4179591836734694E+00/;

      EPMACH = r2mach[4];
      UFLOW = r2mach[1];
      DINF = MIN0[1][INF];

      CENTR = 0.5E+00*(A+B);
      HLGTH = 0.5E+00*(B-A);
      TABSC1 = BOUN+DINF*(0.1E+01-CENTR)/CENTR;
      FVAL1 = F[TABSC1];
       if(INF == 2) FVAL1 = FVAL1+F[-TABSC1];
      FC = (FVAL1/CENTR)/CENTR;

      RESG = WG[8]*FC;
      RESK = WGK[8]*FC;
      RESABS = fabs(RESK);
      for(/*L10*/ J=1; J<=7; J++) {
        ABSC = HLGTH*XGK[J];
        ABSC1 = CENTR-ABSC;
        ABSC2 = CENTR+ABSC;
        TABSC1 = BOUN+DINF*(0.1E+01-ABSC1)/ABSC1;
        TABSC2 = BOUN+DINF*(0.1E+01-ABSC2)/ABSC2;
        FVAL1 = F[TABSC1];
        FVAL2 = F[TABSC2];
         if(INF == 2) FVAL1 = FVAL1+F[-TABSC1];
         if(INF == 2) FVAL2 = FVAL2+F[-TABSC2];
        FVAL1 = (FVAL1/ABSC1)/ABSC1;
        FVAL2 = (FVAL2/ABSC2)/ABSC2;
        FV1[J] = FVAL1;
        FV2[J] = FVAL2;
        FSUM = FVAL1+FVAL2;
        RESG = RESG+WG[J]*FSUM;
        RESK = RESK+WGK[J]*FSUM;
        RESABS = RESABS+WGK[J]*(fabs(FVAL1)+ABS[FVAL2]);
L10:  } //continue;
      RESKH = RESK*0.5E+00;
      RESASC = WGK[8]*fabs(FC-RESKH);
      for(/*L20*/ J=1; J<=7; J++) {
        RESASC = RESASC+WGK[J]*(fabs(FV1[J]-RESKH)+ABS[FV2[J]-RESKH]);
L20:  } //continue;
      RESULT = RESK*HLGTH;
      RESASC = RESASC*HLGTH;
      RESABS = RESABS*HLGTH;
      ABSERR = fabs((RESK-RESG)*HLGTH);
       if(RESASC != 0.0E+00 && ABSERR != 0.E0) ABSERR = RESASC*
           min(0.1E+01,pow((0.2E+03*ABSERR/RESASC),1.5E)+00);
       if(RESABS > UFLOW/(0.5E+02*EPMACH)) ABSERR = max
           [(EPMACH*0.5E+02]*RESABS,ABSERR);
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void QPSRT(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX) {

      double ELIST,ERMAX,ERRMAX,ERRMIN;
      int I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
            NRMAX;
      double ELIST[LAST],IORD[LAST];

       if(LAST > 2) goto L10;
      IORD[1] = 1;
      IORD[2] = 2;
      goto L90;

L10:  ERRMAX = ELIST[MAXERR];
       if(NRMAX == 1) goto L30;
      IDO = NRMAX-1;
      for(/*L20*/ I=1; I<=IDO; I++) {
        ISUCC = IORD[NRMAX-1];

         if(ERRMAX <= ELIST[ISUCC]) goto L30;
        IORD[NRMAX] = ISUCC;
        NRMAX = NRMAX-1;
L20:     } //continue;

L30:  JUPBN = LAST;
       if(LAST > [LIMIT/2+2]) JUPBN = LIMIT+3-LAST;
      ERRMIN = ELIST[LAST];

      JBND = JUPBN-1;
      IBEG = NRMAX+1;
       if(IBEG > JBND) goto L50;
      for(/*L40*/ I=IBEG; I<=JBND; I++) {
        ISUCC = IORD[I];

         if(ERRMAX >= ELIST[ISUCC]) goto L60;
        IORD[I-1] = ISUCC;
L40:  } //continue;
L50:  IORD[JBND] = MAXERR;
      IORD[JUPBN] = LAST;
      goto L90;

L60:  IORD[I-1] = MAXERR;
      K = JBND;
      for(/*L70*/ J=I; J<=JBND; J++) {
        ISUCC = IORD[K];

         if(ERRMIN < ELIST[ISUCC]) goto L80;
        IORD[K+1] = ISUCC;
        K = K-1;
L70:  } //continue;
      IORD[I] = LAST;
      goto L90;
L80:  IORD[K+1] = LAST;

L90:  MAXERR = IORD[NRMAX];
      ERMAX = ELIST[MAXERR];
      return;
      }
//*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function r2mach[I] {
      int SMALL[2];
      int LARGE[2];
      int RIGHT[2];
      int DIVER[2];
      int LOG10[2];
      double RMACH[5] ;
      EQUIVALENCE [RMACH[1]][SMALL[1]];
      EQUIVALENCE [RMACH[2]][LARGE[1]];
      EQUIVALENCE [RMACH[3]][RIGHT[1]];
      EQUIVALENCE [RMACH[4]][DIVER[1]];
      EQUIVALENCE [RMACH[5]][LOG10[1]];
//-------------------------------------------------------------
//     MACHINE CONSTANTS FOR THE IBM PC FAMILY (D. KAHANER NBS)
//-------------------------------------------------------------
      double /*data*/ RMACH/1.18E-38,3.40E+38,0.595E-07,1.19E-07,0.30102999566/;
//-------------------------------------------------------------------
       if((I < 1) || [I > 5])
           CALL xerreur [ 'r2mach -- I OUT OF BOUNDS'][25][1][2];
      r2mach = RMACH[I];
      return;
      }
//**********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void xerreur(MESSG,NMESSG,NERR,LEVEL) {
      char*(*) MESSG;
      CALL XERRWV[MESSG][NMESSG][NERR][LEVEL][0][0][0][0][0.][0.];
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2) {
      char*(*) MESSG;
      printf(*,*) nmessg,nerr;
      printf(*,*) MESSG;
       if(NI == 2){
        printf(*,*) I1,I2;
      }
 else if(NI == 1) {
        printf(*,*) I1;
      }
       if(NR == 2) {
        printf(*,*) R1,R2;
      }
 else if(NR == 1) {
        printf(*,*) R1;
      }
       if(fabs(LEVEL) < 2)return;
      return ;
      }
//*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int function ISAMAX[N][SX][INCX] {
      double SX[*],SMAX,XMAG;
      ISAMAX = 0;
       if(N <= 0) return;
      ISAMAX = 1;
       if(N <= 1)return;
       if(INCX == 1)goto L20;

      SMAX = fabs(SX[1]);
      NS = N*INCX;
      II = 1;
          for(/*L10*/ I=1; I<=NS; I+INCX) {
          XMAG = fabs(SX[I]);
           if(XMAG <= SMAX) goto L5;
          ISAMAX = II;
          SMAX = XMAG;
L5:       II = II + 1;
L10:      } //continue;
      return;

L20:  SMAX = fabs(SX[1]);
      for(/*L30*/ I=2; I<=N; I++) {
         XMAG = fabs(SX[I]);
          if(XMAG <= SMAX) goto L30;
         ISAMAX = I ;
         SMAX = XMAG;
L30:  } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function SASUM[N][SX][INCX] {

      double SX[*];
      SASUM = 0.0E0 ;
       if(N <= 0)return;
       if(INCX == 1)goto L20;

      NS = N*INCX;
          for(/*L10*/ I=1; I<=NS; I+INCX) {
          SASUM = SASUM + fabs(SX[I]);
L10:      } //continue;
      return;

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
      for(/*L30*/ I=1; I<=M; I++) {
        SY[I] = SY[I] + SA*SX[I];
L30:  } //continue;
       if( N  <  4 ) return;
L40:  MP1 = M + 1;
      for(/*L50*/ I=MP1; I<=N; I+4) {
        SY[I] = SY[I] + SA*SX[I];
        SY[I + 1] = SY[I + 1] + SA*SX[I + 1];
        SY[I + 2] = SY[I + 2] + SA*SX[I + 2];
        SY[I + 3] = SY[I + 3] + SA*SX[I + 3];
L50:  } //continue;
      return;

L60:  } //continue;
      NS = N*INCX;
          for(/*L70*/ I=1; I<=NS; I+INCX) {
          SY[I] = SA*SX[I] + SY[I];
L70:      } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SCOPY(N,SX,INCX,SY,INCY) {

      double SX[*],SY[*];
       if(N <= 0)return;
       if(INCX == INCY)  if(INCX-1) 5,20,60;
L5:   } //continue;

      IX = 1;
      IY = 1;
       if(INCX < 0)IX = (-N+1)*INCX + 1 ;
       if(INCY < 0)IY = (-N+1)*INCY + 1 ;
      for(/*L10*/ I=1; I<=N; I++) {
        SY[IY] = SX[IX];
        IX = IX + INCX;
        IY = IY + INCY;
L10:  } //continue;
      return;

L20:  M = ((N)%(7));
       if( M  ==  0 ) goto L40 ;
      for(/*L30*/ I=1; I<=M; I++) {
        SY[I] = SX[I];
L30:  } //continue;
       if( N  <  7 ) return;
L40:  MP1 = M + 1;
      for(/*L50*/ I=MP1; I<=N; I+7) {
        SY[I] = SX[I];
        SY[I + 1] = SX[I + 1] ;
        SY[I + 2] = SX[I + 2] ;
        SY[I + 3] = SX[I + 3] ;
        SY[I + 4] = SX[I + 4] ;
        SY[I + 5] = SX[I + 5] ;
        SY[I + 6] = SX[I + 6] ;
L50:  } //continue;
      return;

L60:  } //continue;
      NS = N*INCX;
          for(/*L70*/ I=1; I<=NS; I+INCX) {
          SY[I] = SX[I];
L70:      } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function SDOT[N][SX][INCX][SY][INCY] {

      double SX[*],SY[*];
      SDOT = 0.0E0;
       if(N <= 0)return;
       if(INCX == INCY) IF[INCX-1]5,20,60;
L5:   } //continue;

      IX = 1;
      IY = 1;
       if(INCX < 0)IX = (-N+1)*INCX + 1 ;
       if(INCY < 0)IY = (-N+1)*INCY + 1 ;
      for(/*L10*/ I=1; I<=N; I++) {
        SDOT = SDOT + SX[IX]*SY[IY];
        IX = IX + INCX;
        IY = IY + INCY;
L10:  } //continue;
      return;

L20:  M = ((N)%(5));
       if( M  ==  0 ) goto L40 ;
      for(/*L30*/ I=1; I<=M; I++) {
        SDOT = SDOT + SX[I]*SY[I];
L30:  } //continue;
       if( N  <  5 ) return;
L40:  MP1 = M + 1;
      for(/*L50*/ I=MP1; I<=N; I+5) {
        SDOT = SDOT + SX[I]*SY[I] + SX[I + 1]*SY[I + 1] +
             SX[I + 2]*SY[I + 2] + SX[I + 3]*SY[I + 3] + SX[I + 4]*SY[I + 4];
L50:  } //continue;
      return;

L60:  } //continue;
      NS=N*INCX;
      for(/*L70*/ I=1; I<=NS; I+INCX) {
        SDOT = SDOT + SX[I]*SY[I];
L70:    } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function SNRM2[N][SX][INCX] {
      int          NEXT;
      double   SX[*],  CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO, ONE;
      double /*data*/   ZERO, ONE /0.0E0, 1.0E0/;

      double /*data*/ CUTLO, CUTHI / 4.441E-16,  1.304E19 /;
       if(N  >  0) goto L10;
         SNRM2  = ZERO;
         goto L300;

L10:  ASSIGN 30 TO NEXT;
      SUM = ZERO;
      NN = N * INCX ;

      I = 1;
L20:     goto LNEXT,(30, 50, 70, 110);
L30:   if( fabs(SX[I])  >  CUTLO) goto L85;
      ASSIGN 50 TO NEXT;
      XMAX = ZERO;

L50:   if( SX[I]  ==  ZERO) goto L200;
       if( fabs(SX[I])  >  CUTLO) goto L85;

      ASSIGN 70 TO NEXT;
      goto L105;

L100: I = J;
      ASSIGN 110 TO NEXT;
      SUM = (SUM / SX[I]) / SX[I];
L105: XMAX = fabs(SX[I]);
      goto L115;

L70:   if( fabs(SX[I])  >  CUTLO ) goto L75;

L110:  if( fabs(SX[I])  <=  XMAX ) goto L115;
         SUM = ONE + SUM *pow( (XMAX / SX[I]),2);
         XMAX = fabs(SX[I]);
         goto L200;

L115: SUM = SUM +pow( (SX[I]/XMAX),2);
      goto L200;

L75:  SUM = (SUM * XMAX) * XMAX;

L85:  HITEST = CUTHI/FLOAT[ N ];

      for(/*L95*/ J=I; J<=NN; J+INCX) {
       if(fabs(SX[J])  >=  HITEST) goto L100;
L95:     SUM = SUM +pow( SX[J],2 );
      SNRM2 = sqrt( SUM );
      goto L300;

L200: } //continue;
      I = I + INCX;
       if ( I  <=  NN ) goto L20;

      SNRM2 = XMAX * sqrt(SUM);
L300: } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SSCAL(N,SA,SX,INCX) {

      double SA,SX[*] ;
       if(N <= 0)return;
       if(INCX == 1)goto L20;

      NS = N*INCX;
          for(/*L10*/ I=1; I<=NS; I+INCX) {
          SX[I] = SA*SX[I];
L10:      } //continue;
      return;

L20:  M = ((N)%(5));
       if( M  ==  0 ) goto L40 ;
      for(/*L30*/ I=1; I<=M; I++) {
        SX[I] = SA*SX[I];
L30:  } //continue;
       if( N  <  5 ) return;
L40:  MP1 = M + 1;
      for(/*L50*/ I=MP1; I<=N; I+5) {
        SX[I] = SA*SX[I];
        SX[I + 1] = SA*SX[I + 1];
        SX[I + 2] = SA*SX[I + 2];
        SX[I + 3] = SA*SX[I + 3];
        SX[I + 4] = SA*SX[I + 4];
L50:  } //continue;
      return;
      }
//c
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SSWAP(N,SX,INCX,SY,INCY) {

      double SX[*],SY[*],STEMP1,STEMP2,STEMP3;
       if(N <= 0)return;
       if(INCX == INCY)  if(INCX-1) 5,20,60;
L5:   } //continue;

      IX = 1;
      IY = 1;
       if(INCX < 0)IX = (-N+1)*INCX + 1 ;
       if(INCY < 0)IY = (-N+1)*INCY + 1 ;
      for(/*L10*/ I=1; I<=N; I++) {
        STEMP1 = SX[IX];
        SX[IX] = SY[IY];
        SY[IY] = STEMP1;
        IX = IX + INCX;
        IY = IY + INCY;
L10:  } //continue;
      return;

L20:  M = ((N)%(3));
       if( M  ==  0 ) goto L40 ;
      for(/*L30*/ I=1; I<=M; I++) {
        STEMP1 = SX[I];
        SX[I] = SY[I];
        SY[I] = STEMP1;
L30:  } //continue;
       if( N  <  3 ) return;
L40:  MP1 = M + 1;
      for(/*L50*/ I=MP1; I<=N; I+3) {
        STEMP1 = SX[I];
        STEMP2 = SX[I+1];
        STEMP3 = SX[I+2];
        SX[I] = SY[I];
        SX[I+1] = SY[I+1];
        SX[I+2] = SY[I+2];
        SY[I] = STEMP1;
        SY[I+1] = STEMP2;
        SY[I+2] = STEMP3;
L50:  } //continue;
      return;
L60:  } //continue;

      NS = N*INCX;
        for(/*L70*/ I=1; I<=NS; I+INCX) {
        STEMP1 = SX[I];
        SX[I] = SY[I];
        SY[I] = STEMP1;
L70:    } //continue;
      return;
      }
//*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double function UNI[] {
      PARAMETER[
              CSAVE=362436./16777216.  ][
              CD=7654321./16777216.][
              CM=16777213./16777216.  ];
//                            2**24=16777216
      double U[17],S,T,USTART,C,UNIB;
      int I,J,II,JJ,K,KK,I1,J1,K1,L1,M1,ISEED;
// 
      SAVE U,I,J,K,C;

      double /*data*/ U/
          0.8668672834288,  0.3697986366357,  0.8008968294805,
          0.4173889774680,  0.8254561579836,  0.9640965269077,
          0.4508667414265,  0.6451309529668,  0.1645456024730,
          0.2787901807898,  0.06761531340295, 0.9663226330820,
          0.01963343943798, 0.02947398211399, 0.1636231515294,
          0.3976343250467,  0.2631008574685/;
      double /*data*/ I,J,K,C/17,5,24,CSAVE/;
// 
      UNI = U[I]-U[J];
       if(UNI < 0.0)UNI = UNI+1.0;
      U[I] = UNI;
      I = I-1;
       if(I == 0)I = 17;
      J = J-1;
       if(J == 0)J = 17;

      C = C-CD;
       if(C < 0.0) C=C+CM;

      UNI = UNI-C;
       if(UNI < 0.0)UNI = UNI+1.0;
      return;

      ENTRY USTART[ISEED];

        I1 = ((fabs(ISEED))%(177))+1;
        J1 = ((fabs(ISEED))%(167))+1;
        K1 = ((fabs(ISEED))%(157))+1;
        L1 = ((fabs(ISEED))%(147))+1;
// 
        for(/*L2*/ II=1; II<=17; II++) {
          S = 0.0;
          T = 0.5;

          for(/*L3*/ JJ=1; JJ<=K; JJ++) {
                  M1 = ((mod(I1*J1,179)*K1)%(179));
                  I1 = J1;
                  J1 = K1;
                  K1 = M1;
                  L1 = ((53*L1+1)%(169));
                   if(((L1*M1)%(64)) >= 32)S=S+T;
L3:               T = .5*T;
L2:     U[II] = S;
        USTART = FLOAT[ISEED];
        return;

      ENTRY UNIB[KK];
         if(KK <= 24){
             K=24;
        }
 else {
             K=KK;
        }
        UNIB=FLOAT[K];
      }
//c*********************************************************************
