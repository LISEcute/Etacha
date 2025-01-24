//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void EQDIF (N,T,Y,F,TOUT,MSTATE,NROOT,EPS,EWT,MINT,WORK,
                 LENW,IWORK,LENIW,G); {
//c*********************************************************************
//c    routine de resolution des equations differentielles
//c*********************************************************************
      EXTERNAL F, G;
      double EPS, EWT, EWTCOM[1], G, HMAX, T, TOUT,
               WORK[*], Y[*];
      int IWORK[*];
      char MSG*81;
      PARAMETER[IMPL = 0][ MXSTEP = 1000];
       if (MINT  <  1  ||  MINT  >  3) {
      printf(MSG, '[''EQDIF1FE Illegal input.  Improper value for ''][
            ''the integration method flag][''][ I8]') MINT;
      CALL xerreur[MSG[1:81]][ 81][ 21][ 2];
        return;
      }
       if (MSTATE  >=  0) {
        NSTATE = MSTATE;
        NTASK = 1;
      }
 else {
        NSTATE = - MSTATE;
        NTASK = 3;
      }
      EWTCOM[1] = EWT;
       if (EWT  !=  0.E0) {
        IERROR = 3;
      }
 else {
        IERROR = 2;
      }
       if (MINT  ==  1) {
        MITER = 0;
        MXORD = 12;
      }
 else if (MINT  ==  2) {
        MITER = 2;
        MXORD = 5;
      }
 else if (MINT  ==  3) {
        MITER = 2;
        MXORD = 12;
      }
      HMAX = 2.E0*fabs(TOUT - T);
      CALL EQDIF3 [N][ T][ Y][ F][ NSTATE][ TOUT][ NTASK][ NROOT][ EPS][ EWTCOM][
                       IERROR][ MINT][ MITER][ IMPL][ ML][ MU][ MXORD][ HMAX][ WORK][
                       LENW][ IWORK][ LENIW][ F][ F][ NDE][ MXSTEP][ G][ F];
       if (MSTATE  >=  0) {
        MSTATE = NSTATE;
      }
 else {
        MSTATE = - NSTATE;
      }
      }
//*********************************************************************** 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SDCOR (DFDY,EL,FA,H,IMPL,IPVT,MATDIM,MITER,ML,MU,N,
             NDE,NQ,T,USERS,Y,YH,YWT,EVALFA,SAVE1,SAVE2,A,D,JSTATE); {
//***ROUTINES CALLED  SGESL,SGBSL,SNRM2
      double A[MATDIM][*], D, DFDY[MATDIM][*], EL[13][12], H,
               SAVE1[*], SAVE2[*], SNRM2, T, Y[*], YH[N][*], YWT[*];
      int IPVT[*];
      LOGICAL EVALFA;
       if (MITER  ==  0) {
        for(/*L100*/ I=1; I<=N; I++) {
L100:     SAVE1[I] = (H*SAVE2[I] - YH[I][2] - SAVE1[I])/YWT[I];
        D = SNRM2[N][ SAVE1][ 1]/sqrt(double[N]);
        for(/*L105*/ I=1; I<=N; I++) {
L105:     SAVE1[I] = H*SAVE2[I] - YH[I][2];
      }
 else if (MITER  ==  1  ||  MITER  ==  2) {
         if (IMPL  ==  0) {
          for(/*L130*/ I=1; I<=N; I++) {
L130:       SAVE2[I] = H*SAVE2[I] - YH[I][2] - SAVE1[I];
        }
 else if (IMPL  ==  1) {
           if (EVALFA) {
            CALL FA [N][ T][ Y][ A][ MATDIM][ ML][ MU][ NDE];
             if (N  ==  0) {
              JSTATE = 9;
              return;
            }
          }
 else {
            EVALFA = .TRUE.;
          }
          for(/*L150*/ I=1; I<=N; I++) {
L150:       SAVE2[I] = H*SAVE2[I];
          for(/*L160*/ J=1; J<=N; J++) {
            for(/*L160*/ I=1; I<=N; I++) {
L160:         SAVE2[I] = SAVE2[I] - A[I][J]*(YH[J][2] + SAVE1[J]);
        }
 else if (IMPL  ==  2) {
           if (EVALFA) {
            CALL FA [N][ T][ Y][ A][ MATDIM][ ML][ MU][ NDE];
             if (N  ==  0) {
              JSTATE = 9;
              return;
            }
          }
 else {
            EVALFA = .TRUE.;
          }
          for(/*L180*/ I=1; I<=N; I++) {
L180:       SAVE2[I] = H*SAVE2[I] - A[I][1]*(YH[I][2] + SAVE1[I]);
        }
        CALL SGESL [DFDY][ MATDIM][ N][ IPVT][ SAVE2][ 0];
        for(/*L200*/ I=1; I<=N; I++) {
          SAVE1[I] = SAVE1[I] + SAVE2[I];
L200:     SAVE2[I] = SAVE2[I]/YWT[I];
        D = SNRM2[N][ SAVE2][ 1]/sqrt(double[N]);
      }
 else if (MITER  ==  4  ||  MITER  ==  5) {
         if (IMPL  ==  0) {
          for(/*L230*/ I=1; I<=N; I++) {
L230:       SAVE2[I] = H*SAVE2[I] - YH[I][2] - SAVE1[I];
        }
 else if (IMPL  ==  1) {
           if (EVALFA) {
            CALL FA [N][ T][ Y][ A[ML+1][1]][ MATDIM][ ML][ MU][ NDE];
             if (N  ==  0) {
              JSTATE = 9;
              return;
            }
          }
 else {
            EVALFA = .TRUE.;
          }
          for(/*L250*/ I=1; I<=N; I++) {
L250:       SAVE2[I] = H*SAVE2[I];
          MW = ML + 1 + MU;
          for(/*L260*/ J=1; J<=N; J++) {
            I1 = MAX[ML+1][ MW+1-J];
            I2 = MIN[MW+N-J][ MW+ML];
            for(/*L260*/ I=I1; I<=I2; I++) {
              I3 = I + J - MW;
L260:         SAVE2[I3] = SAVE2[I3] - A[I][J]*(YH[J][2] + SAVE1[J]);
        }
 else if (IMPL  ==  2) {
           if (EVALFA) {
            CALL FA [N][ T][ Y][ A][ MATDIM][ ML][ MU][ NDE];
             if (N  ==  0) {
              JSTATE = 9;
              return;
            }
          }
 else {
            EVALFA = .TRUE.;
          }
          for(/*L280*/ I=1; I<=N; I++) {
L280:       SAVE2[I] = H*SAVE2[I] - A[I][1]*(YH[I][2] + SAVE1[I]);
        }
        CALL SGBSL [DFDY][ MATDIM][ N][ ML][ MU][ IPVT][ SAVE2][ 0];
        for(/*L300*/ I=1; I<=N; I++) {
          SAVE1[I] = SAVE1[I] + SAVE2[I];
L300:     SAVE2[I] = SAVE2[I]/YWT[I];
        D = SNRM2[N][ SAVE2][ 1]/sqrt(double[N]);
      }
 else if (MITER  ==  3) {
        IFLAG = 2;
        CALL USERS [Y][ YH[1][2]][ YWT][ SAVE1][ SAVE2][ T][ H][ EL[1][NQ]][ IMPL][
                        N][ NDE][ IFLAG];
         if (N  ==  0) {
          JSTATE = 10;
          return;
        }
        for(/*L320*/ I=1; I<=N; I++) {
          SAVE1[I] = SAVE1[I] + SAVE2[I];
L320:     SAVE2[I] = SAVE2[I]/YWT[I];
        D = SNRM2[N][ SAVE2][ 1]/sqrt(double[N]);
      }
      }
//************************************************************ 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SDCST (MAXORD,MINT,ISWFLG,EL,TQ) {
      double EL[13][12], FACTRL[12], GAMMA[14], SUM, TQ[3][12];
      FACTRL[1] = 1.E0;
      for(/*L10*/ I=2; I<=MAXORD; I++) {
L10:    FACTRL[I] = double[I]*FACTRL[I-1];
       if (MINT  ==  1) {
        GAMMA[1] = 1.E0;
        for(/*L40*/ I=1; I<=MAXORD+1; I++) {
          SUM = 0.E0;
          for(/*L30*/ J=1; J<=I; J++) {
L30:        SUM = SUM - GAMMA[J]/double[I-J+2];
L40:      GAMMA[I+1] = SUM;
        EL[1][1] = 1.E0;
        EL[2][1] = 1.E0;
        EL[2][2] = 1.E0;
        EL[3][2] = 1.E0;
        for(/*L60*/ J=3; J<=MAXORD; J++) {
          EL[2][J] = FACTRL[J-1];
          for(/*L50*/ I=3; I<=J; I++) {
L50:        EL[I][J] = double[J-1]*EL[I][J-1] + EL[I-1][J-1];
L60:      EL[J+1][J] = 1.E0;
        for(/*L80*/ J=2; J<=MAXORD; J++) {
          EL[1][J] = EL[1][J-1] + GAMMA[J];
          EL[2][J] = 1.E0;
          for(/*L80*/ I=3; I<=J+1; I++) {
L80:        EL[I][J] = EL[I][J]/(double[I-1]*FACTRL[J-1]);
        for(/*L100*/ J=1; J<=MAXORD; J++) {
          TQ[1][J] = -1.E0/(FACTRL[J]*GAMMA[J]);
          TQ[2][J] = -1.E0/GAMMA[J+1];
L100:     TQ[3][J] = -1.E0/GAMMA[J+2];
      }
 else if (MINT  ==  2) {
        EL[1][1] = 1.E0;
        EL[2][1] = 1.E0;
        for(/*L130*/ J=2; J<=MAXORD; J++) {
          EL[1][J] = FACTRL[J];
          for(/*L120*/ I=2; I<=J; I++) {
L120:       EL[I][J] = double[J]*EL[I][J-1] + EL[I-1][J-1];
L130:     EL[J+1][J] = 1.E0;
        SUM = 1.E0;
        for(/*L150*/ J=2; J<=MAXORD; J++) {
          SUM = SUM + 1.E0/double[J];
          for(/*L150*/ I=1; I<=J+1; I++) {
L150:       EL[I][J] = EL[I][J]/(FACTRL[J]*SUM);
        for(/*L170*/ J=1; J<=MAXORD; J++) {
           if (J  >  1) TQ[1][J] = 1.E0/FACTRL[J-1];
          TQ[2][J] = double[J+1]/EL[1][J];
L170:     TQ[3][J] = double[J+2]/EL[1][J];
      }
       if (ISWFLG  ==  3) {
        MXRD = MIN[MAXORD][ 5];
         if (MINT  ==  2) {
          GAMMA[1] = 1.E0;
          for(/*L190*/ I=1; I<=MXRD; I++) {
            SUM = 0.E0;
            for(/*L180*/ J=1; J<=I; J++) {
L180:         SUM = SUM - GAMMA[J]/double[I-J+2];
L190:       GAMMA[I+1] = SUM;
        }
        SUM = 1.E0;
        for(/*L200*/ I=2; I<=MXRD; I++) {
          SUM = SUM + 1.E0/double[I];
L200:     EL[1+I][1] = -double[I+1]*SUM*GAMMA[I+1];
      }
      }
//********************************************************************* 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SDNTL (EPS,F,FA,HMAX,HOLD,IMPL,JTASK,MATDIM,MAXORD,
             MINT,MITER,ML,MU,N,NDE,SAVE1,T,UROUND,USERS,Y,YWT,H,MNTOLD,
             MTROLD,NFE,RC,YH,A,CONVRG,EL,FAC,IER,IPVT,NQ,NWAIT,RH,RMAX,
             SAVE2,TQ,TREND,ISWFLG,JSTATE); {
      double A[MATDIM][*], EL[13][12], EPS, FAC[*], H, HMAX,
               HOLD, OLDL0, RC, RH, RMAX, RMINIT, SAVE1[*], SAVE2[*], SMAX,
               SMIN, SNRM2, SUM, SUM0, T, TQ[3][12], TREND, UROUND, Y[*],
               YH[N][*], YWT[*];
      int IPVT[*];
      LOGICAL CONVRG, IER;
      PARAMETER[RMINIT = 10000.E0];
      IER = .FALSE.;
       if (JTASK  >=  0) {
         if (JTASK  ==  0) {
          CALL SDCST [MAXORD][ MINT][ ISWFLG][  EL][ TQ];
          RMAX = RMINIT;
        }
        RC = 0.E0;
        CONVRG = .FALSE.;
        TREND = 1.E0;
        NQ = 1;
        NWAIT = 3;
        CALL F [N][ T][ Y][ SAVE2];
         if (N  ==  0) {
          JSTATE = 6;
          return;
        }
        NFE = NFE + 1;
         if (IMPL  !=  0) {
           if (MITER  ==  3) {
            IFLAG = 0;
            CALL USERS [Y][ YH][ YWT][ SAVE1][ SAVE2][ T][ H][ EL][ IMPL][ N][
                            NDE][ IFLAG];
             if (N  ==  0) {
              JSTATE = 10;
              return;
            }
          }
 else if (IMPL  ==  1) {
             if (MITER  ==  1  ||  MITER  ==  2) {
              CALL FA [N][ T][ Y][ A][ MATDIM][ ML][ MU][ NDE];
               if (N  ==  0) {
                JSTATE = 9;
                return;
              }
              CALL SGEFA [A][ MATDIM][ N][ IPVT][ INFO];
               if (INFO  !=  0) {
                IER = .TRUE.;
                return;
              }
              CALL SGESL [A][ MATDIM][ N][ IPVT][ SAVE2][ 0];
            }
 else if (MITER  ==  4  ||  MITER  ==  5) {
              CALL FA [N][ T][ Y][ A[ML+1][1]][ MATDIM][ ML][ MU][ NDE];
               if (N  ==  0) {
                JSTATE = 9;
                return;
              }
              CALL SGBFA [A][ MATDIM][ N][ ML][ MU][ IPVT][ INFO];
               if (INFO  !=  0) {
                IER = .TRUE.;
                return;
              }
              CALL SGBSL [A][ MATDIM][ N][ ML][ MU][ IPVT][ SAVE2][ 0];
            }
          }
 else if (IMPL  ==  2) {
            CALL FA [N][ T][ Y][ A][ MATDIM][ ML][ MU][ NDE];
             if (N  ==  0) {
              JSTATE = 9;
              return;
            }
            for(/*L150*/ I=1; I<=NDE; I++) {
               if (A[I][1]  ==  0.E0) {
                IER = .TRUE.;
                return;
              }
 else {
                SAVE2[I] = SAVE2[I]/A[I][1];
              }
L150:         } //continue;
            for(/*L155*/ I=NDE+1; I<=N; I++) {
L155:         A[I][1] = 0.E0;
          }
        }
        for(/*L170*/ I=1; I<=NDE; I++) {
L170:     SAVE1[I] = SAVE2[I]/YWT[I];
        SUM = SNRM2[NDE][ SAVE1][ 1];
        SUM0 = 1.E0/MAX[1.E0][ fabs(T)];
        SMAX = MAX[SUM0][ SUM];
        SMIN = MIN[SUM0][ SUM];
        SUM = SMAX*sqrt(1.E0 +pow( (SMIN/SMAX),2))/sqrt(double[NDE]);
        H = SIGN[MIN[2.E0*EPS/SUM][ fabs(H)]][ H];
        for(/*L180*/ I=1; I<=N; I++) {
L180:     YH[I][2] = H*SAVE2[I];
         if (MITER  ==  2  ||  MITER  ==  5  ||  ISWFLG  ==  3) {
          for(/*L20*/ I=1; I<=N; I++) {
L20:        FAC[I] = sqrt(UROUND);
        }
      }
 else {
         if (MITER  !=  MTROLD) {
          MTROLD = MITER;
          RC = 0.E0;
          CONVRG = .FALSE.;
        }
         if (MINT  !=  MNTOLD) {
          MNTOLD = MINT;
          OLDL0 = EL[1][NQ];
          CALL SDCST [MAXORD][ MINT][ ISWFLG][  EL][ TQ];
          RC = RC*EL[1][NQ]/OLDL0;
          NWAIT = NQ + 2;
        }
         if (H  !=  HOLD) {
          NWAIT = NQ + 2;
          RH = H/HOLD;
          CALL SDSCL [HMAX][ N][ NQ][ RMAX][  HOLD][ RC][ RH][ YH];
        }
      }
      }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SDNTP (H,K,N,NQ,T,TOUT,YH,Y) {
      double FACTOR, H, R, T, TOUT, Y[*], YH[N][*];
       if (K  ==  0) {
        for(/*L10*/ I=1; I<=N; I++) {
L10:      Y[I] = YH[I][NQ+1];
        R = ((TOUT - T)/H);
        for(/*L20*/ JJ=1; JJ<=NQ; JJ++) {
          J = NQ + 1 - JJ;
          for(/*L20*/ I=1; I<=N; I++) {
L20:        Y[I] = YH[I][J] + R*Y[I];
      }
 else {
        KUSED = MIN[K][ NQ];
        FACTOR = 1.E0;
        for(/*L40*/ KK=1; KK<=KUSED; KK++) {
L40:      FACTOR = FACTOR*double[NQ+1-KK];
        for(/*L50*/ I=1; I<=N; I++) {
L50:      Y[I] = FACTOR*YH[I][NQ+1];
        for(/*L80*/ JJ=KUSED+1; JJ<=NQ; JJ++) {
          J = K + 1 + NQ - JJ;
          FACTOR = 1.E0;
          for(/*L60*/ KK=1; KK<=KUSED; KK++) {
L60:        FACTOR = FACTOR*double[J-KK];
          for(/*L70*/ I=1; I<=N; I++) {
L70:        Y[I] = FACTOR*YH[I][J] + R*Y[I];
L80:      } //continue;
        for(/*L100*/ I=1; I<=N; I++) {
L100:     Y[I] = Y[I]*pow((,-KUSED);)
      }
      }
//******************************************************** 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SDPSC (KSGN,N,NQ,YH) {
      double YH[N][*];
       if (KSGN  >  0) {
        for(/*L10*/ J1=1; J1<=NQ; J1++) {
          for(/*L10*/ J2=J1; J2<=NQ; J2++) {
            J = NQ - J2 + J1;
            for(/*L10*/ I=1; I<=N; I++) {
L10:          YH[I][J] = YH[I][J] + YH[I][J+1];
      }
 else {
        for(/*L30*/ J1=1; J1<=NQ; J1++) {
          for(/*L30*/ J2=J1; J2<=NQ; J2++) {
            J = NQ - J2 + J1;
            for(/*L30*/ I=1; I<=N; I++) {
L30:          YH[I][J] = YH[I][J] - YH[I][J+1];
      }
      }
//********************************************************************** 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SDPST (EL,F,FA,H,IMPL,JACOBN,MATDIM,MITER,ML,MU,N,NDE,
             NQ,SAVE2,T,USERS,Y,YH,YWT,UROUND,NFE,NJE,A,DFDY,FAC,IER,IPVT,
             SAVE1,ISWFLG,BND,JSTATE); {
//***ROUTINES CALLED  SGEFA,SGBFA,SNRM2
      double A[MATDIM][*], BL, BND, BP, BR, BU, DFDY[MATDIM][*],
               DFDYMX, DIFF, DY, EL[13][12], FAC[*], FACMAX, FACMIN, FACTOR,
               H, SAVE1[*], SAVE2[*], SCALE, SNRM2, T, UROUND, Y[*],
               YH[N][*], YJ, YS, YWT[*];
      int IPVT[*];
      LOGICAL IER;
      PARAMETER[FACMAX = .5E0];
      NJE = NJE + 1;
      IER = .FALSE.;
       if (MITER  ==  1  ||  MITER  ==  2) {
         if (MITER  ==  1) {
          CALL JACOBN [N][ T][ Y][ DFDY][ MATDIM][ ML][ MU];
           if (N  ==  0) {
            JSTATE = 8;
            return;
          }
           if (ISWFLG  ==  3) BND = SNRM2[N*N][ DFDY][ 1];
          FACTOR = -EL[1][NQ]*H;
          for(/*L110*/ J=1; J<=N; J++) {
            for(/*L110*/ I=1; I<=N; I++) {
L110:         DFDY[I][J] = FACTOR*DFDY[I][J];
        }
 else if (MITER  ==  2) {
          BR =pow( UROU(.,875E0);
)          BL =pow( UROU(.,75E0);
)          BU =pow( UROU(.,25E0);
)          BP =pow( UROU(-,.15E0);
)          FACMIN =pow( UROU(.,78E0);
)          for(/*L170*/ J=1; J<=N; J++) {
            YS = MAX[fabs(YWT[J])][ ABS[Y[J]]];
L120:       DY = FAC[J]*YS;
             if (DY  ==  0.E0) {
               if (FAC[J]  <  FACMAX) {
                FAC[J] = MIN[100.E0*FAC[J]][ FACMAX];
                goto L120;
              }
 else {
                DY = YS;
              }
            }
             if (NQ  ==  1) {
              DY = SIGN[DY][ SAVE2[J]];
            }
 else {
              DY = SIGN[DY][ YH[J][3]];
            }
            DY = (Y[J] + DY) - Y[J];
            YJ = Y[J];
            Y[J] = Y[J] + DY;
            CALL F [N][ T][ Y][ SAVE1];
             if (N  ==  0) {
              JSTATE = 6;
              return;
            }
            Y[J] = YJ;
            FACTOR = -EL[1][NQ]*H/DY;
            for(/*L140*/ I=1; I<=N; I++) {
L140:         DFDY[I][J] = (SAVE1[I] - SAVE2[I])*FACTOR;
//                                                                 Step 1
            DIFF = fabs(SAVE2[1] - SAVE1[1]);
            IMAX = 1;
            for(/*L150*/ I=2; I<=N; I++) {
               if (fabs(SAVE2[I] - SAVE1[I])  >  DIFF) {
                IMAX = I;
                DIFF = fabs(SAVE2[I] - SAVE1[I]);
              }
L150:         } //continue;
//                                                                 Step 2
             if (MIN[fabs(SAVE2[IMAX])][ ABS[SAVE1[IMAX]]]  >  0.E0) {
              SCALE = MAX[fabs(SAVE2[IMAX])][ ABS[SAVE1[IMAX]]];
//                                                                 Step 3
               if (DIFF  >  BU*SCALE) {
                FAC[J] = MAX[FACMIN][ FAC[J]*.1E0];
              }
 else if (BR*SCALE  <=  DIFF  &&  DIFF  <=  BL*SCALE) {
                FAC[J] = MIN[FAC[J]*10.E0][ FACMAX];
//                                                                 Step 4
              }
 else if (DIFF  <  BR*SCALE) {
                FAC[J] = MIN[BP*FAC[J]][ FACMAX];
              }
            }
L170:       } //continue;
           if (ISWFLG  ==  3) BND = SNRM2[N*N][ DFDY][ 1]/(-EL[1][NQ]*H);
          NFE = NFE + N;
        }
         if (IMPL  ==  0) {
          for(/*L190*/ I=1; I<=N; I++) {
L190:       DFDY[I][I] = DFDY[I][I] + 1.E0;
        }
 else if (IMPL  ==  1) {
          CALL FA [N][ T][ Y][ A][ MATDIM][ ML][ MU][ NDE];
           if (N  ==  0) {
            JSTATE = 9;
            return;
          }
          for(/*L210*/ J=1; J<=N; J++) {
            for(/*L210*/ I=1; I<=N; I++) {
L210:         DFDY[I][J] = DFDY[I][J] + A[I][J];
        }
 else if (IMPL  ==  2) {
          CALL FA [N][ T][ Y][ A][ MATDIM][ ML][ MU][ NDE];
           if (N  ==  0) {
            JSTATE = 9;
            return;
          }
          for(/*L230*/ I=1; I<=NDE; I++) {
L230:       DFDY[I][I] = DFDY[I][I] + A[I][1];
        }
        CALL SGEFA [DFDY][ MATDIM][ N][ IPVT][ INFO];
         if (INFO  !=  0) IER = .TRUE.;
      }
 else if (MITER  ==  4  ||  MITER  ==  5) {
         if (MITER  ==  4) {
          CALL JACOBN [N][ T][ Y][ DFDY[ML+1][1]][ MATDIM][ ML][ MU];
           if (N  ==  0) {
            JSTATE = 8;
            return;
          }
          FACTOR = -EL[1][NQ]*H;
          MW = ML + MU + 1;
          for(/*L260*/ J=1; J<=N; J++) {
            I1 = MAX[ML+1][ MW+1-J];
            I2 = MIN[MW+N-J][ MW+ML];
            for(/*L260*/ I=I1; I<=I2; I++) {
L260:         DFDY[I][J] = FACTOR*DFDY[I][J];
        }
 else if (MITER  ==  5) {
          BR =pow( UROU(.,875E0);
)          BL =pow( UROU(.,75E0);
)          BU =pow( UROU(.,25E0);
)          BP =pow( UROU(-,.15E0);
)          FACMIN =pow( UROU(.,78E0);
)          MW = ML + MU + 1;
          J2 = MIN[MW][ N];
          for(/*L340*/ J=1; J<=J2; J++) {
            for(/*L290*/ K=J; K<=N; K+MW) {
              YS = MAX[fabs(YWT[K])][ ABS[Y[K]]];
L280:         DY = FAC[K]*YS;
               if (DY  ==  0.E0) {
                 if (FAC[K]  <  FACMAX) {
                  FAC[K] = MIN[100.E0*FAC[K]][ FACMAX];
                  goto L280;
                }
 else {
                  DY = YS;
                }
              }
               if (NQ  ==  1) {
                DY = SIGN[DY][ SAVE2[K]];
              }
 else {
                DY = SIGN[DY][ YH[K][3]];
              }
              DY = (Y[K] + DY) - Y[K];
              DFDY[MW][K] = Y[K];
L290:         Y[K] = Y[K] + DY;
            CALL F [N][ T][ Y][ SAVE1];
             if (N  ==  0) {
              JSTATE = 6;
              return;
            }
            for(/*L330*/ K=J; K<=N; K+MW) {
              Y[K] = DFDY[MW][K];
              YS = MAX[fabs(YWT[K])][ ABS[Y[K]]];
              DY = FAC[K]*YS;
               if (DY  ==  0.E0) DY = YS;
               if (NQ  ==  1) {
                DY = SIGN[DY][ SAVE2[K]];
              }
 else {
                DY = SIGN[DY][ YH[K][3]];
              }
              DY = (Y[K] + DY) - Y[K];
              FACTOR = -EL[1][NQ]*H/DY;
              I1 = MAX[ML+1][ MW+1-K];
              I2 = MIN[MW+N-K][ MW+ML];
              for(/*L300*/ I=I1; I<=I2; I++) {
                I3 = K + I - MW;
L300:           DFDY[I][K] = FACTOR*(SAVE1[I3] - SAVE2[I3]);
//                                                                 Step 1
              IMAX = MAX[1][ K - MU];
              DIFF = fabs(SAVE2[IMAX] - SAVE1[IMAX]);
              I1 = IMAX;
              I2 = MIN[K + ML][ N];
              for(/*L310*/ I=I1+1; I<=I2; I++) {
                 if (fabs(SAVE2[I] - SAVE1[I])  >  DIFF) {
                  IMAX = I;
                  DIFF = fabs(SAVE2[I] - SAVE1[I]);
                }
L310:           } //continue;
//                                                                 Step 2
               if (MIN[fabs(SAVE2[IMAX])][ ABS[SAVE1[IMAX]]]  > 0.E0) {
                SCALE = MAX[fabs(SAVE2[IMAX])][ ABS[SAVE1[IMAX]]];
//                                                                 Step 3
                 if (DIFF  >  BU*SCALE) {
                  FAC[K] = MAX[FACMIN][ FAC[K]*.1E0];
                }
 else if (BR*SCALE  <= DIFF  &&  DIFF  <= BL*SCALE) {
                  FAC[K] = MIN[FAC[K]*10.E0][ FACMAX];
//                                                                 Step 4
                }
 else if (DIFF  <  BR*SCALE) {
                  FAC[K] = MIN[BP*FAC[K]][ FACMAX];
                }
              }
L330:         } //continue;
L340:       } //continue;
          NFE = NFE + J2;
        }
         if (ISWFLG  ==  3) {
          DFDYMX = 0.E0;
          for(/*L345*/ J=1; J<=N; J++) {
            I1 = MAX[ML+1][ MW+1-J];
            I2 = MIN[MW+N-J][ MW+ML];
            for(/*L345*/ I=I1; I<=I2; I++) {
L345:         DFDYMX = MAX[DFDYMX][ fabs(DFDY[I][J])];
          BND = 0.E0;
           if (DFDYMX  !=  0.E0) {
            for(/*L350*/ J=1; J<=N; J++) {
              I1 = MAX[ML+1][ MW+1-J];
              I2 = MIN[MW+N-J][ MW+ML];
              for(/*L350*/ I=I1; I<=I2; I++) {
L350:           BND = BND +pow( (DFDY[I][J]/DFDYM2;,
)            BND = DFDYMX*sqrt(BND)/(-EL[1][NQ]*H);
          }
        }
         if (IMPL  ==  0) {
          for(/*L360*/ J=1; J<=N; J++) {
L360:       DFDY[MW][J] = DFDY[MW][J] + 1.E0;
        }
 else if (IMPL  ==  1) {
          CALL FA [N][ T][ Y][ A[ML+1][1]][ MATDIM][ ML][ MU][ NDE];
           if (N  ==  0) {
            JSTATE = 9;
            return;
          }
          for(/*L380*/ J=1; J<=N; J++) {
            I1 = MAX[ML+1][ MW+1-J];
            I2 = MIN[MW+N-J][ MW+ML];
            for(/*L380*/ I=I1; I<=I2; I++) {
L380:         DFDY[I][J] = DFDY[I][J] + A[I][J];
        }
 else if (IMPL  ==  2) {
          CALL FA [N][ T][ Y][ A][ MATDIM][ ML][ MU][ NDE];
           if (N  ==  0) {
            JSTATE = 9;
            return;
          }
          for(/*L400*/ J=1; J<=NDE; J++) {
L400:       DFDY[MW][J] =  DFDY[MW][J] + A[J][1];
        }
        CALL SGBFA [DFDY][ MATDIM][ N][ ML][ MU][ IPVT][ INFO];
         if (INFO  !=  0) IER = .TRUE.;
      }
 else if (MITER  ==  3) {
        IFLAG = 1;
        CALL USERS [Y][ YH[1][2]][ YWT][ SAVE1][ SAVE2][ T][ H][ EL[1][NQ]][ IMPL][
                        N][ NDE][ IFLAG];
         if (N  ==  0) {
          JSTATE = 10;
          return;
        }
      }
      }
//***********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SDSCL (HMAX,N,NQ,RMAX,H,RC,RH,YH) {
      double H, HMAX, RC, RH, RMAX, R1, YH[N][*];
       if (H  <  1.E0) {
        RH = MIN[fabs(H)*RH][ ABS[H]*RMAX][ HMAX]/ABS[H];
      }
 else {
        RH = MIN[RH][ RMAX][ HMAX/fabs(H)];
      }
      R1 = 1.E0;
      for(/*L10*/ J=1; J<=NQ; J++) {
        R1 = R1*RH;
        for(/*L10*/ I=1; I<=N; I++) {
L10:      YH[I][J+1] = YH[I][J+1]*R1;
      H = H*RH;
      RC = RC*RH;
      }
//*********************************************************************** 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SDSTP (EPS,F,FA,HMAX,IMPL,JACOBN,MATDIM,MAXORD,MINT,
             MITER,ML,MU,N,NDE,YWT,UROUND,USERS,AVGH,AVGORD,H,HUSED,JTASK,
             MNTOLD,MTROLD,NFE,NJE,NQUSED,NSTEP,T,Y,YH,A,CONVRG,DFDY,EL,FAC,
             HOLD,IPVT,JSTATE,NQ,NWAIT,RC,RMAX,SAVE1,SAVE2,TQ,TREND,ISWFLG,
             MTRSV,MXRDSV); {
//***ROUTINES CALLED  SDNTL,SDPST,SDCOR,SDPSC,SDSCL,SNRM2
      EXTERNAL F, JACOBN, FA, USERS;
      double A[MATDIM][*], AVGH, AVGORD, BIAS1, BIAS2, BIAS3,
               BND, CTEST, D, DENOM, DFDY[MATDIM][*], D1, EL[13][12], EPS,
               ERDN, ERUP, ETEST, FAC[*], H, HMAX, HN, HOLD, HS, HUSED,
               NUMER, RC, RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL, RMNORM,
               SAVE1[*], SAVE2[*], SNRM2, T, TOLD, TQ[3][12], TREND, TRSHLD,
               UROUND, Y[*], YH[N][*], YWT[*], Y0NRM;
      int IPVT[*];
      LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH;
      PARAMETER[BIAS1 = 1.3E0][ BIAS2 = 1.2E0][ BIAS3 = 1.4E0][ MXFAIL = 3][
                    MXITER = 3][ MXTRY = 50][ RCTEST = .3E0][ RMFAIL = 2.E0][
                    RMNORM = 10.E0][ TRSHLD = 1.E0];
      double /*data*/ IER /.FALSE./;
      NSV = N;
      BND = 0.E0;
      SWITCH = .FALSE.;
      NTRY = 0;
      TOLD = T;
      NFAIL = 0;
       if (JTASK  <=  0) {
        CALL SDNTL [EPS][ F][ FA][ HMAX][ HOLD][ IMPL][ JTASK][ MATDIM][
                        MAXORD][ MINT][ MITER][ ML][ MU][ N][ NDE][ SAVE1][ T][
                        UROUND][ USERS][ Y][ YWT][  H][ MNTOLD][ MTROLD][ NFE][ RC][
                        YH][  A][ CONVRG][ EL][ FAC][ IER][ IPVT][ NQ][ NWAIT][ RH][
                        RMAX][ SAVE2][ TQ][ TREND][ ISWFLG][ JSTATE];
         if (N  ==  0) goto L440;
         if (H  ==  0.E0) goto L400;
         if (IER) goto L420;
      }
L100: NTRY = NTRY + 1;
       if (NTRY  >  MXTRY) goto L410;
      T = T + H;
      CALL SDPSC [1][ N][ NQ][  YH];
      EVALJC = ((fabs(RC - 1.E0)  >  RCTEST)  &&  [MITER  !=  0]);
      EVALFA = .NOT. EVALJC;

L110: ITER = 0;
      for(/*L115*/ I=1; I<=N; I++) {
L115:   Y[I] = YH[I][1];
      CALL F [N][ T][ Y][ SAVE2];
       if (N  ==  0) {
        JSTATE = 6;
        goto L430;
      }
      NFE = NFE + 1;
       if (EVALJC  ||  IER) {
        CALL SDPST [EL][ F][ FA][ H][ IMPL][ JACOBN][ MATDIM][ MITER][ ML][
                        MU][ N][ NDE][ NQ][ SAVE2][ T][ USERS][ Y][ YH][ YWT][ UROUND][
                        NFE][ NJE][  A][ DFDY][ FAC][ IER][ IPVT][ SAVE1][ ISWFLG][
                        BND][ JSTATE];
         if (N  ==  0) goto L430;
         if (IER) goto L160;
        CONVRG = .FALSE.;
        RC = 1.E0;
      }
      for(/*L125*/ I=1; I<=N; I++) {
L125:   SAVE1[I] = 0.E0;

L130: CALL SDCOR [DFDY][ EL][ FA][ H][ IMPL][ IPVT][ MATDIM][ MITER][ ML][
                      MU][ N][ NDE][ NQ][ T][ USERS][ Y][ YH][ YWT][  EVALFA][ SAVE1][
                      SAVE2][  A][ D][ JSTATE];
         if (N  ==  0) goto L430;
       if (ISWFLG  ==  3  &&  MINT  ==  1) {
         if (ITER  ==  0) {
          NUMER = SNRM2[N][ SAVE1][ 1];
          for(/*L132*/ I=1; I<=N; I++) {
L132:       DFDY[1][I] = SAVE1[I];
          Y0NRM = SNRM2[N][ YH][ 1];
        }
 else {
          DENOM = NUMER;
          for(/*L134*/ I=1; I<=N; I++) {
L134:       DFDY[1][I] = SAVE1[I] - DFDY[1][I];
          NUMER = SNRM2[N][ DFDY][ MATDIM];
           if (EL[1][NQ]*NUMER  <=  100.E0*UROUND*Y0NRM) {
             if (RMAX  ==  RMFAIL) {
              SWITCH = .TRUE.;
              goto L170;
            }
          }
          for(/*L136*/ I=1; I<=N; I++) {
L136:       DFDY[1][I] = SAVE1[I];
           if (DENOM  !=  0.E0)
              BND = MAX[BND][ NUMER/(DENOM*fabs(H)*EL[1][NQ]]);
        }
      }
       if (ITER  >  0) TREND = MAX[.9E0*TREND][ D/D1];
      D1 = D;
      CTEST = MIN[2.E0*TREND][ 1.E0]*D;
       if (CTEST  <=  EPS) goto L170;
      ITER = ITER + 1;
       if (ITER  <  MXITER) {
        for(/*L140*/ I=1; I<=N; I++) {
L140:     Y[I] = YH[I][1] + EL[1][NQ]*SAVE1[I];
        CALL F [N][ T][ Y][ SAVE2];
         if (N  ==  0) {
          JSTATE = 6;
          goto L430;
        }
        NFE = NFE + 1;
        goto L130;
      }

       if (CONVRG) {
        EVALJC = .TRUE.;
        EVALFA = .FALSE.;
        goto L110;
      }
L160: T = TOLD;
      CALL SDPSC [-1][ N][ NQ][  YH];
      NWAIT = NQ + 2;
       if (JTASK  !=  0  &&  JTASK  !=  2) RMAX = RMFAIL;
       if (ITER  ==  0) {
        RH = .3E0;
      }
 else {
        RH = .9E0*pow((EPS/CTES(.,2E0);
)      }
       if (RH*H  ==  0.E0) goto L400;
      CALL SDSCL [HMAX][ N][ NQ][ RMAX][  H][ RC][ RH][ YH];
      goto L100;

L170: CONVRG = (MITER  !=  0);
      for(/*L180*/ I=1; I<=NDE; I++) {
L180:   SAVE2[I] = SAVE1[I]/YWT[I];
      ETEST = SNRM2[NDE][ SAVE2][ 1]/(TQ[2][NQ]*sqrt(double[NDE]));

       if (ETEST  >  EPS) {
        T = TOLD;
        CALL SDPSC [-1][ N][ NQ][  YH];
        NFAIL = NFAIL + 1;
         if (NFAIL  <  MXFAIL) {
           if (JTASK  !=  0  &&  JTASK  !=  2) RMAX = RMFAIL;
          RH2 = 1.E0/(BIAS2*pow((ETEST/EP(1,.E0/double[NQ+1]));)
           if (NQ  >  1) {
            for(/*L190*/ I=1; I<=NDE; I++) {
L190:         SAVE2[I] = YH[I][NQ+1]/YWT[I];
            ERDN = SNRM2[NDE][ SAVE2][ 1]/(TQ[1][NQ]*sqrt(double[NDE]));
            RH1 = 1.E0/MAX[1.E0][ BIAS1*(ERDN/pow(EP(1,.E0/double[NQ]));)
             if (RH2  <  RH1) {
              NQ = NQ - 1;
              RC = RC*EL[1][NQ]/EL[1][NQ+1];
              RH = RH1;
            }
 else {
              RH = RH2;
            }
          }
 else {
            RH = RH2;
          }
          NWAIT = NQ + 2;
           if (RH*H  ==  0.E0) goto L400;
          CALL SDSCL [HMAX][ N][ NQ][ RMAX][  H][ RC][ RH][ YH];
          goto L100;
        }

        NFAIL = 0;
        JTASK = 2;
        for(/*L215*/ I=1; I<=N; I++) {
L215:     Y[I] = YH[I][1];
        CALL SDNTL [EPS][ F][ FA][ HMAX][ HOLD][ IMPL][ JTASK][ MATDIM][
                        MAXORD][ MINT][ MITER][ ML][ MU][ N][ NDE][ SAVE1][ T][
                        UROUND][ USERS][ Y][ YWT][  H][ MNTOLD][ MTROLD][ NFE][ RC][
                        YH][  A][ CONVRG][ EL][ FAC][ IER][ IPVT][ NQ][ NWAIT][ RH][
                        RMAX][ SAVE2][ TQ][ TREND][ ISWFLG][ JSTATE];
        RMAX = RMNORM;
         if (N  ==  0) goto L440;
         if (H  ==  0.E0) goto L400;
         if (IER) goto L420;
        goto L100;
      }

      NSTEP = NSTEP + 1;
      HUSED = H;
      NQUSED = NQ;
      AVGH = (double[NSTEP-1]*AVGH + H)/REAL[NSTEP];
      AVGORD = (double[NSTEP-1]*AVGORD + REAL[NQ])/REAL[NSTEP];
      for(/*L230*/ J=1; J<=NQ+1; J++) {
        for(/*L230*/ I=1; I<=N; I++) {
L230:     YH[I][J] = YH[I][J] + EL[J][NQ]*SAVE1[I];
      for(/*L235*/ I=1; I<=N; I++) {
L235:   Y[I] = YH[I][1];
//                                          If ISWFLG is 3, consider
//                                          changing integration methods.
       if (ISWFLG  ==  3) {
         if (BND  !=  0.E0) {
           if (MINT  ==  1  &&  NQ  <=  5) {
            HN = fabs(H)/MAX[UROUND][ (ETEST/pow(EP(1,.E0/double[NQ+1]));)
            HN = MIN[HN][ 1.E0/(2.E0*EL[1][NQ]*BND]);
            HS = fabs(H)/MAX[UROUND][
                [ETEST/(EPS*EL[NQ+pow(1][1](1,.E0/double[NQ+1]));)
             if (HS  >  1.2E0*HN) {
              MINT = 2;
              MNTOLD = MINT;
              MITER = MTRSV;
              MTROLD = MITER;
              MAXORD = MIN[MXRDSV][ 5];
              RC = 0.E0;
              RMAX = RMNORM;
              TREND = 1.E0;
              CALL SDCST [MAXORD][ MINT][ ISWFLG][ EL][ TQ];
              NWAIT = NQ + 2;
            }
          }
 else if (MINT  ==  2) {
            HS = fabs(H)/MAX[UROUND][ (ETEST/pow(EP(1,.E0/double[NQ+1]));)
            HN = fabs(H)/MAX[UROUND][
                [ETEST*EL[NQ+1][1]/pow(EP(1,.E0/double[NQ+1]]);
)            HN = MIN[HN][ 1.E0/(2.E0*EL[1][NQ]*BND]);
             if (HN  >=  HS) {
              MINT = 1;
              MNTOLD = MINT;
              MITER = 0;
              MTROLD = MITER;
              MAXORD = MIN[MXRDSV][ 12];
              RMAX = RMNORM;
              TREND = 1.E0;
              CONVRG = .FALSE.;
              CALL SDCST [MAXORD][ MINT][ ISWFLG][ EL][ TQ];
              NWAIT = NQ + 2;
            }
          }
        }
      }
       if (SWITCH) {
        MINT = 2;
        MNTOLD = MINT;
        MITER = MTRSV;
        MTROLD = MITER;
        MAXORD = MIN[MXRDSV][ 5];
        NQ = MIN[NQ][ MAXORD];
        RC = 0.E0;
        RMAX = RMNORM;
        TREND = 1.E0;
        CALL SDCST [MAXORD][ MINT][ ISWFLG][ EL][ TQ];
        NWAIT = NQ + 2;
      }

       if (JTASK  ==  0  ||  JTASK  ==  2) {
        RH = 1.E0/MAX[UROUND][ BIAS2*(ETEST/pow(EP(1,.E0/double[NQ+1]));)
         if (RH > TRSHLD) CALL SDSCL [HMAX][ N][ NQ][ RMAX][ H][ RC][ RH][ YH];
      }
 else if (NWAIT  >  1) {
        NWAIT = NWAIT - 1;
         if (NWAIT  ==  1  &&  NQ  <  MAXORD) {
          for(/*L250*/ I=1; I<=NDE; I++) {
L250:       YH[I][MAXORD+1] = SAVE1[I];
        }

      }
 else {
         if (NQ  ==  1) {
          RH1 = 0.E0;
        }
 else {
          for(/*L270*/ I=1; I<=NDE; I++) {
L270:       SAVE2[I] = YH[I][NQ+1]/YWT[I];
          ERDN = SNRM2[NDE][ SAVE2][ 1]/(TQ[1][NQ]*sqrt(double[NDE]));
          RH1 = 1.E0/MAX[UROUND][ BIAS1*(ERDN/pow(EP(1,.E0/double[NQ]));)
        }
        RH2 = 1.E0/MAX[UROUND][ BIAS2*(ETEST/pow(EP(1,.E0/double[NQ+1]));)
         if (NQ  ==  MAXORD) {
          RH3 = 0.E0;
        }
 else {
          for(/*L290*/ I=1; I<=NDE; I++) {
L290:       SAVE2[I] = (SAVE1[I] - YH[I][MAXORD+1])/YWT[I];
          ERUP = SNRM2[NDE][ SAVE2][ 1]/(TQ[3][NQ]*sqrt(double[NDE]));
          RH3 = 1.E0/MAX[UROUND][ BIAS3*(ERUP/pow(EP(1,.E0/double[NQ+2]));)
        }
         if (RH1  >  RH2  &&  RH1  >=  RH3) {
          RH = RH1;
           if (RH  <=  TRSHLD) goto L380;
          NQ = NQ - 1;
          RC = RC*EL[1][NQ]/EL[1][NQ+1];
        }
 else if (RH2  >=  RH1  &&  RH2  >=  RH3) {
          RH = RH2;
           if (RH  <=  TRSHLD) goto L380;
        }
 else {
          RH = RH3;
           if (RH  <=  TRSHLD) goto L380;
          for(/*L360*/ I=1; I<=N; I++) {
L360:       YH[I][NQ+2] = SAVE1[I]*EL[NQ+1][NQ]/double[NQ+1];
          NQ = NQ + 1;
          RC = RC*EL[1][NQ]/EL[1][NQ-1];
        }
         if (ISWFLG  ==  3  &&  MINT  ==  1) {
           if (BND != 0.E0) RH = MIN[RH][ 1.E0/(2.E0*EL[1][NQ]*BND*fabs(H)]);
        }
        CALL SDSCL [HMAX][ N][ NQ][ RMAX][  H][ RC][ RH][ YH];
        RMAX = RMNORM;
L380:   NWAIT = NQ + 2;
      }

      JSTATE = 1;
      HOLD = H;
      return;

L400: JSTATE = 2;
      HOLD = H;
      for(/*L405*/ I=1; I<=N; I++) {
L405:   Y[I] = YH[I][1];
      return;

L410: JSTATE = 3;
      HOLD = H;
      return;

L420: JSTATE = 4;
      HOLD = H;
      return;

L430: T = TOLD;
      CALL SDPSC [-1][ NSV][ NQ][  YH];
      for(/*L435*/ I=1; I<=NSV; I++) {
L435:   Y[I] = YH[I][1];
L440: HOLD = H;
      return;
      }
//********************************************************************** 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SDZRO (AE,F,H,N,NQ,IROOT,RE,T,YH,UROUND,B,C,FB,FC,Y) {
//***ROUTINES CALLED  SDNTP
      double A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC,
               H, P, Q, RE, RW, T, TOL, UROUND, Y[*], YH[N][*];
      ER = 4.E0*UROUND;
      RW = MAX[RE][ ER];
      IC = 0;
      ACBS = fabs(B - C);
      A = C;
      FA = FC;
      KOUNT = 0;

L10:   if (fabs(FC)  <  ABS[FB]) {
        A = B;
        FA = FB;
        B = C;
        FB = FC;
        C = A;
        FC = FA;
      }
      CMB = 0.5E0*(C - B);
      ACMB = fabs(CMB);
      TOL = RW*fabs(B) + AE;

       if (ACMB  <=  TOL) return;
       if (KOUNT  >  50) return;

      P = (B - A)*FB;
      Q = FA - FB;
       if (P  <  0.E0) {
        P = -P;
        Q = -Q;
      }

      A = B;
      FA = FB;
      IC = IC + 1;
       if (IC  >=  4) {
         if (8.E0*ACMB  >=  ACBS) {

          B = 0.5E0*(C + B);
          goto L20;
        }
        IC = 0;
      }
      ACBS = ACMB;

       if (P  <=  fabs(Q)*TOL) {
        B = B + SIGN[TOL][ CMB];
      }
 else if (P  <  CMB*Q) {
        B = B + P/Q;
      }
 else {
        B = 0.5E0*(C + B);
      }
L20:  CALL SDNTP [H][ 0][ N][ NQ][ T][ B][ YH][  Y];
      FB = F[N][ B][ Y][ IROOT];
       if (N  ==  0) return;
       if (FB  ==  0.E0) return;
      KOUNT = KOUNT + 1;

       if (SIGN[1.0E0][ FB]  ==  SIGN[1.0E0][ FC]) {
        C = A;
        FC = FA;
      }
      goto L10;
      }
//*********************************************************************** 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void EQDIF3 (N,T,Y,F,NSTATE,TOUT,NTASK,NROOT,EPS,EWT,IERROR,
             MINT,MITER,IMPL,ML,MU,MXORD,HMAX,WORK,LENW,IWORK,LENIW,JACOBN,
             FA,NDE,MXSTEP,G,USERS); {

      EXTERNAL F, JACOBN, FA, G, USERS;
      double AE, BIG, EPS, EWT[*], G, GLAST, H, HMAX, HSIGN,
                   NROUND, RE, r2mach, SIZE, SNRM2, SUM, T, TLAST, TOUT, TROOT,
               UROUND, WORK[*], Y[*];
      int IWORK[*];
      LOGICAL CONVRG;
      char MSG*205;
      PARAMETER[NROUND = 20.E0];
      PARAMETER[IAVGH = 1][ IHUSED = 2][ IAVGRD = 3][
                    IEL = 4][ IH = 160][ IHMAX = 161][ IHOLD = 162][
                    IHSIGN = 163][ IRC = 164][ IRMAX = 165][ IT = 166][
                    ITOUT = 167][ ITQ = 168][ ITREND = 204][ IYH = 205][
                    INDMXR = 1][ INQUSD = 2][ INSTEP = 3][ INFE = 4][ INJE = 5][
                    INROOT = 6][ ICNVRG = 7][ IJROOT = 8][ IJTASK = 9][
                    IMNTLD = 10][ IMTRLD = 11][ INQ = 12][ INRTLD = 13][
                    INDTRT = 14][ INWAIT = 15][ IMNT = 16][ IMTRSV = 17][
                    IMTR = 18][ IMXRDS = 19][ IMXORD = 20][ INDPRT = 21][
                    INDPVT = 22];
      NPAR = N;
      UROUND = r2mach [4];
       if (NROOT  !=  0) {
      AE = r2mach[1];
        RE = UROUND;
      }
       if (EPS  <  0.E0) {
      printf(MSG, '[''EQDIF36FE Illegal input.  EPS][''][ E16.8][
            ''][ is negative.'']') EPS;
      CALL xerreur[MSG[1:60]][ 60][ 6][ 2];
        return;
      }
       if (N  <=  0) {
      printf(MSG, '[''EQDIF37FE Illegal input.  Number of equations][''][
            I8][ ''][ is not positive.'']') N;
      CALL xerreur[MSG[1:72]][ 72][ 7][ 2];
        return;
      }
       if (MXORD  <=  0) {
      printf(MSG, '[''EQDIF314FE Illegal input.  Maximum order][''][ I8][
            ''][ is not positive.'']') MXORD;
      CALL xerreur[MSG[1:67]][ 67][ 14][ 2];
        return;
      }
       if ((MINT  <  1  ||  MINT  >  3)  ||  [MINT  ==  3  && 
            [MITER  ==  0  ||  MITER  ==  3  ||  IMPL  !=  0]]
             ||  [MITER  <  0  ||  MITER  >  5]  || 
            [IMPL  !=  0  &&  IMPL  !=  1  &&  IMPL  !=  2]  || 
            [(IMPL  ==  1  ||  IMPL  ==  2]  &&  MITER  ==  0)  || 
            [IMPL  ==  2  &&  MINT  ==  1]  || 
            [NSTATE  <  1  ||  NSTATE  >  10]) {
      printf(MSG, '[''EQDIF39FE Illegal input.  Improper value for ''][
            ''NSTATE[MSTATE]][ MINT][ MITER or IMPL.'']');
      CALL xerreur[MSG[1:81]][ 81][ 9][ 2];
        return;
      }
       if (MITER  ==  0  ||  MITER  ==  3) {
        LIWCHK = INDPVT - 1;
      }
 else if (MITER  ==  1  ||  MITER  ==  2  ||  MITER  ==  4  || 
            MITER  ==  5) {
        LIWCHK = INDPVT + N - 1;
      }
       if (LENIW  <  LIWCHK) {
      printf(MSG, '[''EQDIF310FE Illegal input.  Insufficient ''][
            ''storage allocated for the IWORK array.  Based on the '']');
        printf(MSG[94:], '[''value of the input parameters involved][ ''][
            ''the required storage is''][ I8]') LIWCHK;
      CALL xerreur[MSG[1:164]][ 164][ 10][ 2];
        return;
      }
       if (MINT  ==  1  ||  MINT  ==  3) {
        MAXORD = MIN[MXORD][ 12];
      }
 else if (MINT  ==  2) {
        MAXORD = MIN[MXORD][ 5];
      }
      IDFDY = IYH + (MAXORD + 1)*N;

       if (MITER  ==  0  ||  MITER  ==  3)  {
        IYWT = IDFDY;
      }
 else if (MITER  ==  1  ||  MITER  ==  2)  {
        IYWT = IDFDY + N*N;
      }
 else if (MITER  ==  4  ||  MITER  ==  5)  {
        IYWT = IDFDY + (2*ML + MU + 1)*N;
      }
      ISAVE1 = IYWT + N;
      ISAVE2 = ISAVE1 + N;
      IGNOW = ISAVE2 + N;
      ITROOT = IGNOW + NROOT;
      IFAC = ITROOT + NROOT;
       if (MITER  ==  2  ||  MITER  ==  5  ||  MINT  ==  3) {
        IA = IFAC + N;
      }
 else {
        IA = IFAC;
      }
       if (IMPL  ==  0  ||  MITER  ==  3) {
        LENCHK = IA - 1;
      }
 else if (IMPL  ==  1  &&  [MITER  ==  1  ||  MITER  ==  2]) {
        LENCHK = IA - 1 + N*N;
      }
 else if (IMPL  ==  1  &&  [MITER  ==  4  ||  MITER  ==  5]) {
        LENCHK = IA - 1 + (2*ML + MU + 1)*N;
      }
 else if (IMPL  ==  2  &&  MITER  !=  3) {
        LENCHK = IA - 1 + N;
      }
       if (LENW  <  LENCHK) {
      printf(MSG, '[''EQDIF38FE Illegal input.  Insufficient ''][
            ''storage allocated for the WORK array.  Based on the '']');
        printf(MSG[92:], '[''value of the input parameters involved][ ''][
            ''the required storage is''][ I8]') LENCHK;
      CALL xerreur[MSG[1:162]][ 162][ 8][ 2];
        return;
      }
       if (MITER  ==  0  ||  MITER  ==  3) {
        MATDIM = 1;
      }
 else if (MITER  ==  1  ||  MITER  ==  2) {
        MATDIM = N;
      }
 else if (MITER  ==  4  ||  MITER  ==  5) {
        MATDIM = 2*ML + MU + 1;
      }
       if (IMPL  ==  0  ||  IMPL  ==  1) {
        NDECOM = N;
      }
 else if (IMPL  ==  2) {
        NDECOM = NDE;
      }
       if (NSTATE  ==  1) {

         if (MINT  ==  1  ||  MINT  ==  3) {
          IWORK[IMXORD] = MIN[MXORD][ 12];
        }
 else if (MINT  ==  2) {
          IWORK[IMXORD] = MIN[MXORD][ 5];
        }
        IWORK[IMXRDS] = MXORD;
         if (MINT  ==  1  ||  MINT  ==  2) {
          IWORK[IMNT] = MINT;
          IWORK[IMTR] = MITER;
          IWORK[IMNTLD] = MINT;
          IWORK[IMTRLD] = MITER;
        }
 else if (MINT  ==  3) {
          IWORK[IMNT] = 1;
          IWORK[IMTR] = 0;
          IWORK[IMNTLD] = IWORK[IMNT];
          IWORK[IMTRLD] = IWORK[IMTR];
          IWORK[IMTRSV] = MITER;
        }
        WORK[IHMAX] = HMAX;
        H = (TOUT - T)*(1.E0 - 4.E0*UROUND);
        H = SIGN[MIN[fabs(H)][ HMAX]][ H];
        WORK[IH] = H;
        HSIGN = SIGN[1.E0][ H];
        WORK[IHSIGN] = HSIGN;
        IWORK[IJTASK] = 0;
        WORK[IAVGH] = 0.E0;
        WORK[IHUSED] = 0.E0;
        WORK[IAVGRD] = 0.E0;
        IWORK[INDMXR] = 0;
        IWORK[INQUSD] = 0;
        IWORK[INSTEP] = 0;
        IWORK[INFE] = 0;
        IWORK[INJE] = 0;
        IWORK[INROOT] = 0;
        WORK[IT] = T;
        IWORK[ICNVRG] = 0;
        IWORK[INDPRT] = 0;

        for(/*L30*/ I=1; I<=N; I++) {
          JYH = I + IYH - 1;
L30:      WORK[JYH] = Y[I];
         if (T  ==  TOUT) return;
        goto L180;
      }

       if (IWORK[ICNVRG]  ==  1) {
        CONVRG = .TRUE.;
      }
 else {
        CONVRG = .FALSE.;
      }
      T = WORK[IT];
      H = WORK[IH];
      HSIGN = WORK[IHSIGN];
       if (IWORK[IJTASK]  ==  0) goto L180;

       if (NROOT  !=  0) {
        JROOT = IWORK[IJROOT];
         if (JROOT  >  0) {
           if (NSTATE  !=  5) {
             if (TOUT*HSIGN  >=  WORK[ITOUT]*HSIGN) {
              TROOT = WORK[ITOUT];
              CALL SDNTP[H][ 0][ N][ IWORK[INQ]][ T][ TROOT][ WORK[IYH]][  Y];
              T = TROOT;
              NSTATE = 5;
              goto L580;
            }
          }
 else {
            TROOT = T;
            IROOT = 0;
            for(/*L50*/ I=; I<=IWORK[INRTLD; I++]) {
              JTROOT = IWORK[INDTRT] + I - 1;
               if (WORK[JTROOT]*HSIGN  <=  TROOT*HSIGN) {

                 if (WORK[JTROOT]  ==  WORK[ITOUT]  && 
                    I  >  IWORK[INROOT]) {
                  IROOT = I;
                  TROOT = WORK[JTROOT];
                  goto L60;
                }
                 if (WORK[JTROOT]*HSIGN  >  WORK[ITOUT]*HSIGN) {
                  IROOT = I;
                  TROOT = WORK[JTROOT];
                }
              }
L50:          } //continue;
L60:        IWORK[INROOT] = IROOT;
            WORK[ITOUT] = TROOT;
            IWORK[IJROOT] = NTASK;
             if (NTASK  ==  1) {
               if (IROOT  ==  0) {
                IWORK[IJROOT] = 0;
              }
 else {
                 if (TOUT*HSIGN  >=  TROOT*HSIGN) {
                  CALL SDNTP[H][ 0][ N][ IWORK[INQ]][ T][ TROOT][WORK[IYH]][Y];
                  NSTATE = 5;
                  T = TROOT;
                  goto L580;
                }
              }
            }
 else if (NTASK  ==  2  ||  NTASK  ==  3) {

               if (IROOT  ==  0  ||  [TOUT*HSIGN  <  TROOT*HSIGN]) {
                IWORK[IJROOT] = 0;
              }
 else {
                CALL SDNTP[H][ 0][ N][ IWORK[INQ]][ T][ TROOT][ WORK[IYH]][ Y];
                NSTATE = 5;
                T = TROOT;
                goto L580;
              }
            }
          }
        }
      }

       if (NTASK  ==  1) {
        NSTATE = 2;
         if (T*HSIGN  >=  TOUT*HSIGN) {
          CALL SDNTP [H][ 0][ N][ IWORK[INQ]][ T][ TOUT][ WORK[IYH]][  Y];
          T = TOUT;
          goto L580;
        }
      }
 else if (NTASK  ==  2) {
         if (T*HSIGN  >  TOUT*HSIGN) {
        printf(MSG, '[''EQDIF32WRN With NTASK=''][ I1][ '' on input][ ''][
              ''T][''][ E16.8][ ''][ was beyond TOUT][''][ E16.8][ ''.  Solution''][
              '' obtained by interpolation.'']') NTASK, T, TOUT;
        CALL xerreur[MSG[1:124]][ 124][ 2][ 0];
          CALL SDNTP [H][ 0][ N][ IWORK[INQ]][ T][ TOUT][ WORK[IYH]][  Y];
          T = TOUT;
          NSTATE = 2;
          goto L580;
        }

         if (fabs(TOUT - T) <= NROUND*UROUND*MAX[ABS[T]][ ABS[TOUT]]) {
          T = TOUT;
          NSTATE = 2;
          goto L560;
        }
         if (NSTATE  ==  5) {
          NSTATE = 2;
          goto L560;
        }
        NSTATE = 2;
         if ((T + H)*HSIGN  >  TOUT*HSIGN) {
          H = TOUT - T;
           if ((T + H)*HSIGN  >  TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND);
          WORK[IH] = H;
           if (H  ==  0.E0) goto L670;
          IWORK[IJTASK] = -1;
        }
      }
 else if (NTASK  ==  3) {
        NSTATE = 2;
         if (T*HSIGN  >  TOUT*HSIGN) {
        printf(MSG, '[''EQDIF32WRN With NTASK=''][ I1][ '' on input][ ''][
              ''T][''][ E16.8][ ''][ was beyond TOUT][''][ E16.8][ ''.  Solution''][
              '' obtained by interpolation.'']') NTASK, T, TOUT;
        CALL xerreur[MSG[1:124]][ 124][ 2][ 0];
          CALL SDNTP [H][ 0][ N][ IWORK[INQ]][ T][ TOUT][ WORK[IYH]][  Y];
          T = TOUT;
          goto L580;
        }
         if (fabs(TOUT - T) <= NROUND*UROUND*MAX[ABS[T]][ ABS[TOUT]]) {
          T = TOUT;
          goto L560;
        }
         if ((T + H)*HSIGN  >  TOUT*HSIGN) {
          H = TOUT - T;
           if ((T + H)*HSIGN  >  TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND);
          WORK[IH] = H;
           if (H  ==  0.E0) goto L670;
          IWORK[IJTASK] = -1;
        }
      }

       if ((MINT  !=  IWORK[IMNTLD]  ||  MITER  !=  IWORK[IMTRLD])  && 
            MINT  !=  3  &&  IWORK[IMNTLD]  !=  3) IWORK[IJTASK] = -1;
       if (HMAX  !=  WORK[IHMAX]) {
        H = SIGN[MIN[fabs(H)][ HMAX]][ H];
         if (H  !=  WORK[IH]) {
          IWORK[IJTASK] = -1;
          WORK[IH] = H;
        }
        WORK[IHMAX] = HMAX;
      }

L180: NSTEPL = IWORK[INSTEP];
      for(/*L190*/ I=1; I<=N; I++) {
        JYH = IYH + I - 1;
L190:   Y[I] = WORK[JYH];
       if (NROOT  !=  0) {
        for(/*L200*/ I=1; I<=NROOT; I++) {
          JGNOW = IGNOW + I - 1;
          WORK[JGNOW] = G [NPAR][ T][ Y][ I];
           if (NPAR  ==  0) {
            IWORK[INROOT] = I;
            NSTATE = 7;
            return;
          }
L200:    } //continue;
      }
       if (IERROR  ==  1) {
        for(/*L230*/ I=1; I<=N; I++) {
          JYWT = I + IYWT - 1;
L230:     WORK[JYWT] = 1.E0;
        goto L410;
      }
 else if (IERROR  ==  5) {
        for(/*L250*/ I=1; I<=N; I++) {
          JYWT = I + IYWT - 1;
L250:     WORK[JYWT] = EWT[I];
        goto L410;
      }
//                                       Reset YWT array.  Looping point.
L260:  if (IERROR  ==  2) {
        for(/*L280*/ I=1; I<=N; I++) {
           if (Y[I]  ==  0.E0) goto L290;
          JYWT = I + IYWT - 1;
L280:     WORK[JYWT] = fabs(Y[I]);
        goto L410;
L290:    if (IWORK[IJTASK]  ==  0) {
          CALL F [NPAR][ T][ Y][ WORK[ISAVE2]];
           if (NPAR  ==  0) {
            NSTATE = 6;
            return;
          }
          IWORK[INFE] = IWORK[INFE] + 1;
           if (MITER  ==  3  &&  IMPL  !=  0) {
            IFLAG = 0;
            CALL USERS[Y][ WORK[IYH]][ WORK[IYWT]][ WORK[ISAVE1]][
                           WORK[ISAVE2]][ T][ H][ WORK[IEL]][ IMPL][ NPAR][
                           NDECOM][ IFLAG];
             if (NPAR  ==  0) {
              NSTATE = 10;
              return;
            }
          }
 else if (IMPL  ==  1) {
             if (MITER  ==  1  ||  MITER  ==  2) {
              CALL FA [NPAR][ T][ Y][ WORK[IA]][ MATDIM][ ML][ MU][ NDECOM];
               if (NPAR  ==  0) {
                NSTATE = 9;
                return;
              }
              CALL SGEFA [WORK[IA]][ MATDIM][ N][ IWORK[INDPVT]][ INFO];
               if (INFO  !=  0) goto L690;
              CALL SGESL[WORK[IA]][MATDIM][N][IWORK[INDPVT]][WORK[ISAVE2]][0];
            }
 else if (MITER  ==  4  ||  MITER  ==  5) {
              JAML = IA + ML;
              CALL FA [NPAR][ T][ Y][ WORK[JAML]][ MATDIM][ ML][ MU][ NDECOM];
               if (NPAR  ==  0) {
                NSTATE = 9;
                return;
              }
              CALL SGBFA [WORK[IA]][MATDIM][N][ML][MU][IWORK[INDPVT]][INFO];
               if (INFO  !=  0) goto L690;
              CALL SGBSL [WORK[IA]][ MATDIM][ N][ ML][ MU][ IWORK[INDPVT]][
                              WORK[ISAVE2]][ 0];
            }
          }
 else if (IMPL  ==  2) {
            CALL FA [NPAR][ T][ Y][ WORK[IA]][ MATDIM][ ML][ MU][ NDECOM];
             if (NPAR  ==  0) {
              NSTATE = 9;
              return;
            }
            for(/*L340*/ I=1; I<=NDECOM; I++) {
              JA = I + IA - 1;
              JSAVE2 = I + ISAVE2 - 1;
               if (WORK[JA]  ==  0.E0) goto L690;
L340:         WORK[JSAVE2] = WORK[JSAVE2]/WORK[JA];
          }
        }
        for(/*L360*/ J=I; J<=N; J++) {
          JYWT = J + IYWT - 1;
           if (Y[J]  !=  0.E0) {
            WORK[JYWT] = fabs(Y[J]);
          }
 else {
             if (IWORK[IJTASK]  ==  0) {
              JSAVE2 = J + ISAVE2 - 1;
              WORK[JYWT] = fabs(H*WORK[JSAVE2]);
            }
 else {
              JHYP = J + IYH + N - 1;
              WORK[JYWT] = fabs(WORK[JHYP]);
            }
          }
           if (WORK[JYWT]  ==  0.E0) WORK[JYWT] = UROUND;
L360:     } //continue;
      }
 else if (IERROR  ==  3) {
        for(/*L380*/ I=1; I<=N; I++) {
          JYWT = I + IYWT - 1;
L380:     WORK[JYWT] = MAX[EWT[1]][ fabs(Y[I])];
      }
 else if (IERROR  ==  4) {
        for(/*L400*/ I=1; I<=N; I++) {
          JYWT = I + IYWT - 1;
L400:     WORK[JYWT] = MAX[EWT[I]][ fabs(Y[I])];
      }

L410: for(/*L420*/ I=1; I<=N; I++) {
        JYWT = I + IYWT - 1;
        JSAVE2 = I + ISAVE2 - 1;
L420:   WORK[JSAVE2] = Y[I]/WORK[JYWT];
      SUM = SNRM2[N][ WORK[ISAVE2]][ 1]/sqrt(double[N]);
       if (EPS  <  SUM*UROUND) {
        EPS = SUM*UROUND*(1.E0 + 10.E0*UROUND);
      printf(MSG, '[''EQDIF34REC At T][''][ E16.8][ ''][ the requested ''][
            ''accuracy][ EPS][ was not obtainable with the machine ''][
            ''precision.  EPS has been increased to'']') T;
        printf(MSG[137:], '[E16.8]') EPS;
      CALL xerreur[MSG[1:152]][ 152][ 4][ 1];
        NSTATE = 4;
        goto L560;
      }
       if (fabs(H)  >=  UROUND*ABS[T]) {
        IWORK[INDPRT] = 0;
      }
 else if (IWORK[INDPRT]  ==  0) {
//     WRITE(MSG, '(''EQDIF35WRN At T,'', E16.8, '', the step size,'',
//    8  E16.8, '', is smaller than the roundoff level of T.  '')') T, H
//        WRITE(MSG(109:), '(''This may occur if there is an abrupt '',
//     8  ''change in the right hand side of the differential '',
//     8  ''equations.'')')
//     CALL xerreur(MSG(1:205), 205, 5, 0)
        IWORK[INDPRT] = 1;
      }
       if (NTASK != 2) {
         if ((IWORK[INSTEP]-NSTEPL)  >  MXSTEP) {
        printf(MSG, '[''EQDIF33WRN At T][''][ E16.8][ ''][ ''][ I8][
              '' steps have been taken without reaching TOUT][''][ E16.8]')
              T, MXSTEP, TOUT;
        CALL xerreur[MSG[1:103]][ 103][ 3][ 0];
          NSTATE = 3;
          goto L560;
        }
      }

      CALL SDSTP [EPS][ F][ FA][ WORK[IHMAX]][ IMPL][ JACOBN][ MATDIM][
                      IWORK[IMXORD]][ IWORK[IMNT]][ IWORK[IMTR]][ ML][ MU][ NPAR][
                     NDECOM][ WORK[IYWT]][ UROUND][ USERS][  WORK[IAVGH]][
                     WORK[IAVGRD]][ WORK[IH]][ WORK[IHUSED]][ IWORK[IJTASK]][
                     IWORK[IMNTLD]][ IWORK[IMTRLD]][ IWORK[INFE]][ IWORK[INJE]][
                      IWORK[INQUSD]][ IWORK[INSTEP]][ WORK[IT]][ Y][ WORK[IYH]][
                      WORK[IA]][ CONVRG][ WORK[IDFDY]][ WORK[IEL]][ WORK[IFAC]][
                      WORK[IHOLD]][ IWORK[INDPVT]][ JSTATE][ IWORK[INQ]][
                      IWORK[INWAIT]][ WORK[IRC]][ WORK[IRMAX]][ WORK[ISAVE1]][
                      WORK[ISAVE2]][ WORK[ITQ]][ WORK[ITREND]][ MINT][
                      IWORK[IMTRSV]][ IWORK[IMXRDS]];
      T = WORK[IT];
      H = WORK[IH];
      goto L(470, 670, 680, 690, 690, 660, 660, 660, 660, 660), JSTATE;
L470: IWORK[IJTASK] = 1;
       if (NROOT  !=  0) {
        IROOT = 0;
        for(/*L500*/ I=1; I<=NROOT; I++) {
          JTROOT = ITROOT + I - 1;
          JGNOW = IGNOW + I - 1;
          GLAST = WORK[JGNOW];
          WORK[JGNOW] = G [NPAR][ T][ Y][ I];
           if (NPAR  ==  0) {
            IWORK[INROOT] = I;
            NSTATE = 7;
            return;
          }
           if (GLAST*WORK[JGNOW]  >  0.E0) {
            WORK[JTROOT] = T + H;
          }
 else {
             if (WORK[JGNOW]  ==  0.E0) {
              WORK[JTROOT] = T;
              IROOT = I;
            }
 else {
               if (GLAST  ==  0.E0) {
                WORK[JTROOT] = T + H;
              }
 else {
                 if (fabs(WORK[IHUSED])  >=  UROUND*ABS[T]) {
                  TLAST = T - WORK[IHUSED];
                  IROOT = I;
                  TROOT = T;
                  CALL SDZRO [AE][ G][ H][ NPAR][ IWORK[INQ]][ IROOT][ RE][ T][
                                  WORK[IYH]][ UROUND][  TROOT][ TLAST][
                                  WORK[JGNOW]][ GLAST][  Y];
                  for(/*L480*/ J=1; J<=N; J++) {
L480:               Y[J] = WORK[IYH + J -1] ;
                   if (NPAR  ==  0) {
                    IWORK[INROOT] = I;
                    NSTATE = 7;
                    return;
                  }
                  WORK[JTROOT] = TROOT;
                }
 else {
                  WORK[JTROOT] = T;
                  IROOT = I;
                }
              }
            }
          }
L500:     } //continue;
         if (IROOT  ==  0) {
          IWORK[IJROOT] = 0;
        }
 else {
          IWORK[IJROOT] = NTASK;
          IWORK[INRTLD] = NROOT;
          IWORK[INDTRT] = ITROOT;
          TROOT = T + H;
          for(/*L510*/ I=1; I<=NROOT; I++) {
            JTROOT = ITROOT + I - 1;
             if (WORK[JTROOT]*HSIGN  <  TROOT*HSIGN) {
              TROOT = WORK[JTROOT];
              IROOT = I;
            }
L510:       } //continue;
          IWORK[INROOT] = IROOT;
          WORK[ITOUT] = TROOT;
           if (TROOT*HSIGN  <=  TOUT*HSIGN) {
            CALL SDNTP [H][ 0][ N][ IWORK[INQ]][ T][ TROOT][ WORK[IYH]][  Y];
            NSTATE = 5;
            T = TROOT;
            goto L580;
          }
        }
      }

      NSTATE = 2;
       if (NTASK  ==  1) {
         if (T*HSIGN  <  TOUT*HSIGN) goto L260;
        CALL SDNTP [H][ 0][ N][ IWORK[INQ]][ T][ TOUT][ WORK[IYH]][  Y];
        T = TOUT;
        goto L580;

      }
 else if (NTASK  ==  2) {
         if (fabs(TOUT - T) <= NROUND*UROUND*MAX[ABS[T]][ ABS[TOUT]]) {
          T = TOUT;
        }
 else {
           if ((T + H)*HSIGN  >  TOUT*HSIGN) {
            H = TOUT - T;
             if ((T + H)*HSIGN > TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND);
            WORK[IH] = H;
             if (H  ==  0.E0) goto L670;
            IWORK[IJTASK] = -1;
          }
        }
      }
 else if (NTASK  ==  3) {
         if (fabs(TOUT - T) <= NROUND*UROUND*MAX[ABS[T]][ ABS[TOUT]]) {
          T = TOUT;
        }
 else {
           if ((T + H)*HSIGN  >  TOUT*HSIGN) {
            H = TOUT - T;
             if ((T + H)*HSIGN > TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND);
            WORK[IH] = H;
             if (H  ==  0.E0) goto L670;
            IWORK[IJTASK] = -1;
          }
          goto L260;
        }
      }

L560: for(/*L570*/ I=1; I<=N; I++) {
        JYH = I + IYH - 1;
L570:   Y[I] = WORK[JYH];
L580:  if (CONVRG) {
        IWORK[ICNVRG] = 1;
      }
 else {
        IWORK[ICNVRG] = 0;
      }
       if (IWORK[IJTASK]  ==  0) return;
      BIG = 0.E0;
      IMXERR = 1;
      IWORK[INDMXR] = IMXERR;
      for(/*L590*/ I=1; I<=N; I++) {
//                                            SIZE = ABS(ERROR(I)/YWT(I))
        JYWT = I + IYWT - 1;
        JERROR = I + ISAVE1 - 1;
        SIZE = fabs(WORK[JERROR]/WORK[JYWT]);
         if (BIG  <  SIZE) {
          BIG = SIZE;
          IMXERR = I;
          IWORK[INDMXR] = IMXERR;
        }
L590:   } //continue;
      return;

L660: NSTATE = JSTATE;
      return;

L670: printf(MSG, '[''EQDIF311FE At T][''][ E16.8][ ''][ the attempted ''][
            ''step size has gone to zero.  Often this occurs  if the ''][
            ''problem setup is incorrect.'']') T;
      CALL xerreur[MSG[1:129]][ 129][ 11][ 2];
      return;

L680: printf(MSG, '[''EQDIF312FE At T][''][ E16.8][ ''][ the step size has''][
            '' been reduced about 50 times without advancing the '']') T;
      printf(MSG[103:], '[''solution.  Often this occurs  if the ''][
            ''problem setup is incorrect.'']');
      CALL xerreur[MSG[1:165]][ 165][ 12][ 2];
      return;

L690: printf(MSG, '[''EQDIF313FE At T][''][ E16.8][ ''][ while solving''][
            '' A*YDOT = F][ A is singular.'']') T;
      CALL xerreur[MSG[1:74]][ 74][ 13][ 2];
      return;
      }
//*********************************************************************** 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)  {

      int LDA,N,ML,MU,IPVT[*],INFO;
      double ABD[LDA][*];

      double T;
      int I,ISAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1 ;

      M = ML + MU + 1;
      INFO = 0;

      J0 = MU + 2;
      J1 = MIN0[N][M] - 1;
       if (J1  <  J0) goto L30;
      for(/*L20*/ JZ=J0; JZ<=J1; JZ++) {
         I0 = M + 1 - JZ;
         for(/*L10*/ I=I0; I<=ML; I++) {
            ABD[I][JZ] = 0.0E0 ;
L10:     } //continue;
L20:  } //continue;
L30:  } //continue;
      JZ = J1;
      JU = 0;

      NM1 = N - 1;
       if (NM1  <  1) goto L130;
      for(/*L120*/ K=1; K<=NM1; K++) {
         KP1 = K + 1;

         JZ = JZ + 1;
          if (JZ  >  N) goto L50;
          if (ML  <  1) goto L50;
            for(/*L40*/ I=1; I<=ML; I++) {
               ABD[I][JZ] = 0.0E0;
L40:        } //continue;
L50:     } //continue;

         LM = MIN0[ML][N-K];
         L = ISAMAX[LM+1][ABD[M][K]][1] + M - 1;
         IPVT[K] = L + K - M;

          if (ABD[L][K]  ==  0.0E0) goto L100;

             if (L  ==  M) goto L60;
               T = ABD[L][K];
               ABD[L][K] = ABD[M][K];
               ABD[M][K] = T;
L60:        } //continue;

            T = -1.0E0/ABD[M][K];
            CALL SSCAL[LM][T][ABD[M+1][K]][1];

            JU = MIN0[MAX0[JU][MU+IPVT[K]]][N];
            MM = M;
             if (JU  <  KP1) goto L90;
            for(/*L80*/ J=KP1; J<=JU; J++) {
               L = L - 1;
               MM = MM - 1;
               T = ABD[L][J];
                if (L  ==  MM) goto L70;
                  ABD[L][J] = ABD[MM][J];
                  ABD[MM][J] = T;
L70:           } //continue;
               CALL SAXPY[LM][T][ABD[M+1][K]][1][ABD[MM+1][J]][1];
L80:        } //continue;
L90:        } //continue;
         goto L110;
L100:    } //continue;
            INFO = K;
L110:    } //continue;
L120: } //continue;
L130: } //continue;
      IPVT[N] = N;
       if (ABD[M][N]  ==  0.0E0) INFO = N ;
      return;
      }
//*********************************************************** 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB) {

      int LDA,N,ML,MU,IPVT[*],JOB;
      double ABD[LDA][*],B[*];

      double SDOT,T;
      int K,KB,L,LA,LB,LM,M,NM1;
      M = MU + ML + 1;
      NM1 = N - 1;
       if (JOB  !=  0) goto L50;

          if (ML  ==  0) goto L30;
          if (NM1  <  1) goto L30;
            for(/*L20*/ K=1; K<=NM1; K++) {
               LM = MIN0[ML][N-K];
               L = IPVT[K];
               T = B[L];
                if (L  ==  K) goto L10;
                  B[L] = B[K] ;
                  B[K] = T;
L10:           } //continue;
               CALL SAXPY[LM][T][ABD[M+1][K]][1][B[K+1]][1];
L20:        } //continue;
L30:     } //continue;

         for(/*L40*/ KB=1; KB<=N; KB++) {
            K = N + 1 - KB;
            B[K] = B[K]/ABD[M][K];
            LM = MIN0[K][M] - 1;
            LA = M - LM;
            LB = K - LM;
            T = -B[K];
            CALL SAXPY[LM][T][ABD[LA][K]][1][B[LB]][1];
L40:     } //continue;
      goto L100;
L50:  } //continue;

         for(/*L60*/ K=1; K<=N; K++) {
            LM = MIN0[K][M] - 1;
            LA = M - LM;
            LB = K - LM;
            T = SDOT[LM][ABD[LA][K]][1][B[LB]][1];
            B[K] = (B[K] - T)/ABD[M][K];
L60:     } //continue;

          if (ML  ==  0) goto L90;
          if (NM1  <  1) goto L90;
            for(/*L80*/ KB=1; KB<=NM1; KB++) {
               K = N - KB;
               LM = MIN0[ML][N-K];
               B[K] = B[K] + SDOT[LM][ABD[M+1][K]][1][B[K+1]][1] ;
               L = IPVT[K];
                if (L  ==  K) goto L70;
                  T = B[L];
                  B[L] = B[K] ;
                  B[K] = T;
L70:           } //continue;
L80:        } //continue;
L90:     } //continue;
L100: } //continue;
      return;
      }
//****************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SGEFA(A,LDA,N,IPVT,INFO) {
      int LDA,N,IPVT[*],INFO;
      double A[LDA][*] ;

      double T;
      int ISAMAX,J,K,KP1,L,NM1;

      INFO = 0;
      NM1 = N - 1;
       if (NM1  <  1) goto L70;
      for(/*L60*/ K=1; K<=NM1; K++) {
         KP1 = K + 1;

         L = ISAMAX[N-K+1][A[K][K]][1] + K - 1;
         IPVT[K] = L;

          if (A[L][K]  ==  0.0E0) goto L40;

             if (L  ==  K) goto L10;
               T = A[L][K];
               A[L][K] = A[K][K];
               A[K][K] = T;
L10:        } //continue;

            T = -1.0E0/A[K][K] ;
            CALL SSCAL[N-K][T][A[K+1][K]][1];

            for(/*L30*/ J=KP1; J<=N; J++) {
               T = A[L][J];
                if (L  ==  K) goto L20;
                  A[L][J] = A[K][J];
                  A[K][J] = T;
L20:           } //continue;
               CALL SAXPY[N-K][T][A[K+1][K]][1][A[K+1][J]][1];
L30:        } //continue;
         goto L50;
L40:     } //continue;
            INFO = K;
L50:     } //continue;
L60:  } //continue;
L70:  } //continue;
      IPVT[N] = N;
       if (A[N][N]  ==  0.0E0) INFO = N;
      return;
      }
//*********************************************************************** 
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SGESL(A,LDA,N,IPVT,B,JOB) {
      int LDA,N,IPVT[*],JOB;
      double A[LDA][*],B[*];

      double SDOT,T;
      int K,KB,L,NM1;
      NM1 = N - 1;
       if (JOB  !=  0) goto L50;

          if (NM1  <  1) goto L30;
         for(/*L20*/ K=1; K<=NM1; K++) {
            L = IPVT[K];
            T = B[L];
             if (L  ==  K) goto L10;
               B[L] = B[K];
               B[K] = T;
L10:        } //continue;
            CALL SAXPY[N-K][T][A[K+1][K]][1][B[K+1]][1] ;
L20:     } //continue;
L30:     } //continue;

         for(/*L40*/ KB=1; KB<=N; KB++) {
            K = N + 1 - KB;
            B[K] = B[K]/A[K][K];
            T = -B[K];
            CALL SAXPY[K-1][T][A[1][K]][1][B[1]][1];
L40:     } //continue;
      goto L100;
L50:  } //continue;

         for(/*L60*/ K=1; K<=N; K++) {
            T = SDOT[K-1][A[1][K]][1][B[1]][1];
            B[K] = (B[K] - T)/A[K][K];
L60:     } //continue;

          if (NM1  <  1) goto L90;
         for(/*L80*/ KB=1; KB<=NM1; KB++) {
            K = N - KB;
            B[K] = B[K] + SDOT[N-K][A[K+1][K]][1][B[K+1]][1];
            L = IPVT[K];
             if (L  ==  K) goto L70;
               T = B[L];
               B[L] = B[K];
               B[K] = T;
L70:        } //continue;
L80:     } //continue;
L90:     } //continue;
L100: } //continue;
      return;
      }
