#include "../win/e_myextern.h"
#include "../win/e_Constant.h"
#include <stdio.h>
#include <stdlib.h>
#include <QtMath>
#include <QDebug>


void INTG( double (*F)(double x), double BOUND, int INF,
                 double EPSABS, double EPSREL, double &RESULT, double &ABSERR,
                 int &NEVAL, int &IER, int LIMIT, int LENW, int &LAST,
                 int *IWORK, double *WORK);

void INTGE(     double (*F)(double x), double BOUND, int INF,
                double EPSABS, double EPSREL,
                int LIMIT,
                double &RESULT, double &ABSERR,
                int &NEVAL, int &IER,
                double *ALIST, double *BLIST, double *RLIST, double *ELIST,
                int *IORD,
                int &LAST
                );

void EA(bool &NEWFLG, double SVALUE,int LIMEXP, double &RESULT, double &ABSERR, double *EPSTAB, int &IERR);

void QK15I( double (*F)(double x), double BOUND, int INF,
        double A, double B,
        double &RESULT, double &ABSERR,
        double &RESABS, double &RESASC);

void QPSRT(int LIMIT, int LAST, int &MAXERR, double &ERMAX, double *ELIST, int *IORD, int &NRMAX);

double r2mach(int I);
void xerreur(char *MESSG, int NMESSG, int NERR, int LEVEL);
void XERRWV (char *MESSG, int NMESSG, int NERR, int LEVEL, int NI,int I1,int I2, int NR,double R1, double R2);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void INTG( double (*F)(double x), double BOUND, int INF, double EPSABS, double EPSREL, double &RESULT, double &ABSERR,
                 int &NEVAL, int &IER, int LIMIT, int LENW, int &LAST, int *IWORK, double *WORK)

 //[F][BOUND][INF][EPSABS][EPSREL][RESULT][ABSERR][NEVAL][
//             IER][LIMIT][LENW][LAST][IWORK][WORK]

{
//c*********************************************************************
//c    Routine d'int‚gration
//c*********************************************************************
      int LVL,L1,L2,L3;

      IER = 6;
      NEVAL = 0;
      LAST = 0;
      RESULT = 0.0E+00;
      ABSERR = 0.0E+00;
      if(LIMIT < 1 || LENW < LIMIT*4) goto L10;

      L1 = LIMIT+1;
      L2 = LIMIT+L1;
      L3 = LIMIT+L2;

//      CALL INTGE[F][BOUND][INF][EPSABS][EPSREL][LIMIT][RESULT][ABSERR][
//            NEVAL][IER][WORK[1]][WORK[L1]][WORK[L2]][WORK[L3]][IWORK][LAST];

         INTGE(
                F, BOUND, INF,
                EPSABS, EPSREL,
                LIMIT,
                RESULT, ABSERR,
                NEVAL, IER,
                &WORK[0], &WORK[L1], &WORK[L2], &WORK[L3],
                IWORK,
                LAST
                );

       LVL = 0;
L10:   if(IER == 6) LVL = 1;
       if(IER != 0)  qDebug() << "ABNORMAL return FROM  INTG" << 26 << IER << LVL;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void EA(bool &NEWFLG, double SVALUE,int LIMEXP, double &RESULT, double &ABSERR, double *EPSTAB, int &IERR)
{
double RES3LA[4];

double EPRN, RELPR, RES;
double DELTA1,DELTA2,DELTA3,ERR1,ERR2,ERR3,ERROR,E0,E1,E2,E3,SS,TOL1,TOL2,TOL3;

int IB,IB2, IE, IN, K1, K2, K3, N, NEWELM, NRES, NUM;

/*      double  ERROR, SS;
      int I,IE,IERR,IN,K1,K2,K3,LIMEXP,N,NEWELM,NUM,NRES;
*/

    if(LIMEXP < 3)
        {
        IERR = 1;
        qDebug() << "LIMEXP IS LESS THAN 3" << 21 << 1 << 1;
        return;
        }

      IERR = 0;
      RES3LA[1] = EPSTAB[LIMEXP+5];
      RES3LA[2] = EPSTAB[LIMEXP+6];
      RES3LA[3] = EPSTAB[LIMEXP+7];
      RESULT    =  SVALUE;

if(NEWFLG)
        {
        N=1;
        NRES=0;
        NEWFLG=false;
        EPSTAB[N]=SVALUE;
        ABSERR=fabs(RESULT);
        goto L100;
        }
 else   {
        N    =int(EPSTAB[LIMEXP+3]);
        NRES =int(EPSTAB[LIMEXP+4]);

        if(N == 2)
                {
                EPSTAB[N]=SVALUE;
                ABSERR=.6E+01*fabs(RESULT-EPSTAB[1]);
                goto L100;
                }
        }
      EPSTAB[N]=SVALUE;
      RELPR = r2mach(4);
      EPRN  = 10.*RELPR;
      EPSTAB[N+2]=EPSTAB[N];
      NEWELM = (N-1)/2;
      NUM=N;
      K1=N;


///---------------------------------------
      for(int I=1; I<=NEWELM; I++)  //L40
        {
        K2 = K1-1;
        K3 = K1-2;
        RES = EPSTAB[K1+2];
        E0  = EPSTAB[K3];
        E1  = EPSTAB[K2];
        E2  = RES;
        DELTA2 = E2-E1;
        DELTA3 = E1-E0;
        ERR2 = fabs(DELTA2);
        ERR3 = fabs(DELTA3);
        TOL2 = qMax(fabs(E2),fabs(E1))*RELPR;
        TOL3 = qMax(fabs(E1),fabs(E0))*RELPR;

        if(ERR2 > TOL2 || ERR3 > TOL3) goto L10;

        RESULT=RES;
        ABSERR=ERR2+ERR3;
        goto L50;

L10:
        if(I != 1)
          {
          E3 = EPSTAB[K1];
          EPSTAB[K1]=E1;
          DELTA1 = E1-E3;
          ERR1   = fabs(DELTA1);
          TOL1   = qMax(fabs(E1),fabs(E3))*RELPR;

          if(ERR1 <= TOL1 || ERR2 <= TOL2 || ERR3 <= TOL3) goto L20;
          SS=1./DELTA1+1./DELTA2-1./DELTA3;
          }
         else
          {
          EPSTAB[K1]=E1;
          if(ERR2 <= TOL2 || ERR3 <= TOL3) goto L20;
          SS=1./DELTA2-1./DELTA3;
          }

       if(fabs(SS*E1) > 0.1E-03) goto L30;

L20:   N=I+I-1;
              if(NRES == 0)
                 {
                 ABSERR=ERR2+ERR3;
                 RESULT=RES;
                 }
        else if(NRES == 1) RESULT=RES3LA[1];
        else if(NRES == 2) RESULT=RES3LA[2];
        else               RESULT=RES3LA[3];

        goto L50;

L30:    RES = E1 + 1./SS;
        EPSTAB[K1]=RES;
        K1=K1-2;

        if(NRES == 0)
                {
                ABSERR=ERR2+fabs(RES-E2)+ERR3;
                RESULT=RES;
                continue;
                }
         else if(NRES == 1) ERROR=6.*(fabs(RES-RES3LA[1]));
         else if(NRES == 2) ERROR=2.*(fabs(RES-RES3LA[2])+fabs(RES-RES3LA[1]));
         else               ERROR=    fabs(RES-RES3LA[3])+fabs(RES-RES3LA[2])+fabs(RES-RES3LA[1]);

      if(ERROR > 10.*ABSERR) continue;

      ABSERR=ERROR;
      RESULT=RES;

      } //continue;   L40:
//----------------------------

      if(NRES == 1)   ABSERR=6.*fabs(RESULT-RES3LA[1]);
 else if(NRES == 2)   ABSERR=2.*fabs(RESULT-RES3LA[2]) + fabs(RESULT-RES3LA[1]);
 else if(NRES > 2)    ABSERR=   fabs(RESULT-RES3LA[3]) + fabs(RESULT-RES3LA[2])  + fabs(RESULT-RES3LA[1]);


L50:   if(N == LIMEXP) N=2*(LIMEXP/2)-1;

      IB=1;                             // ODD
      if((NUM/2)*2 == NUM) IB=2;        // EVEN

      IE = NEWELM+1;

      for(int  I=1; I<=IE; I++)  //*L60
        {
        IB2=IB+2;
        EPSTAB[IB]=EPSTAB[IB2];
        IB=IB2;
        }


      if(NUM == N) goto L80;
      IN = NUM-N+1;

      for(int  I=1; I<=N; I++)  //L70
        {
        EPSTAB[I]=EPSTAB[IN];
        IN++;
        }


L80:
       if(NRES == 0)    RES3LA[1]=RESULT;
 else  if(NRES == 1)    RES3LA[2]=RESULT;
 else  if(NRES == 2)    RES3LA[3]=RESULT;
 else {
        RES3LA[1] =RES3LA[2];
        RES3LA[2] =RES3LA[3];
        RES3LA[3] =RESULT;
      }

//L90:
      ABSERR=qMax(ABSERR,EPRN*fabs(RESULT));
      NRES++;

L100: N++;
//  commented originally : *      if(N <= 3) ABSERR = r2mach[2] * (0.1D-03);
      EPSTAB[LIMEXP+3]=double(N);
      EPSTAB[LIMEXP+4]=double(NRES);
      EPSTAB[LIMEXP+5]=RES3LA[1];
      EPSTAB[LIMEXP+6]=RES3LA[2];
      EPSTAB[LIMEXP+7]=RES3LA[3];
//L110:

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void INTGE(
                double (*F)(double x), double BOUND, int INF,
                double EPSABS, double EPSREL,
                int LIMIT,
                double &RESULT, double &ABSERR,
                int &NEVAL, int &IER,
                double *ALIST, double *BLIST, double *RLIST, double *ELIST,
                int *IORD,
                int &LAST
                )
{

    double DEFABS,RESABS, BOUN, DRES,ERRBND,ERRSUM, UFLOW, CORREC=0;


//    double ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,
//            A2,BLIST,BOUN,BOUND,B1,B2,CORREC,DEFABS,DEFAB1,DEFAB2,
//            DRES,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND,
//            ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,ERTEST,F,RESABS,
//            RESEPS,RESULT,RLIST,RLIST2,r2mach,SMALL,UFLOW;

//    int ID,IER,IERR,IERRO,INF,IORD,IROFF1,IROFF2,IROFF3,JUPBND,K,
//            KSGN,KTMIN,LAST,LIMEXP,LIMIT,MAXERR,NEVAL,NRMAX;

int ID, KTMIN, NRMAX, IERRO, IROFF1, IROFF2, IROFF3, JUPBND, KSGN, MAXERR;
double ERRMAX,AREA;
double SMALL = 0.375; //  Oleg init
double ERTEST = 7.37e-3;    //  Oleg init
double ERLARG = 0;          //  Oleg init

#define LIMEXP 50
#define LIMEXP2 (LIMEXP+8)

bool EXTRAP,LERR,NEWFLG,NOEXT;

double  RESEPS=0, ABSEPS=0;
int IERR=0;

double RLIST2a[LIMEXP2];

double EPMACH = r2mach(4);

      IER = 0;
      NEVAL = 0;
      LAST = 0;
      RESULT = 0.;
      ABSERR = 0.;
      ALIST[1] = 0.;
      BLIST[1] = 1;
      RLIST[1] = 0;
      ELIST[1] = 0;
      IORD[1] = 0;
      NEWFLG = true;

     if(EPSABS <= 0 && EPSREL < qMax(50.*EPMACH,5E-15)) IER = 6;

     if(IER == 6) return; //goto L999;

      BOUN = BOUND;
      if(INF == 2) BOUN = 0;

        QK15I(F, BOUN, INF,0,1,RESULT,ABSERR,DEFABS,RESABS);


      LAST = 1;
      RLIST[1] = RESULT;
      ELIST[1] = ABSERR;
      IORD[1] = 1;
      DRES = fabs(RESULT);
      ERRBND = qMax(EPSABS,EPSREL*DRES);

      if(ABSERR <= 100.*EPMACH*DEFABS  &&  ABSERR >  ERRBND) IER = 2;

      if(LIMIT == 1) IER = 1;
      if(IER != 0 || (ABSERR <= ERRBND &&  ABSERR != RESABS) ||  ABSERR == 0) goto L130;

      UFLOW = r2mach(1) * 2.;
      LERR = false;


//      EA[NEWFLG][RESULT][LIMEXP][RESEPS][ABSEPS][RLIST2][IERR];
//  EA(bool &NEWFLG, double SVALUE,int LIMEXP, double &RESULT, double &ABSERR, double *EPSTAB, int &IERR)
    EA(NEWFLG, RESULT, LIMEXP, RESEPS, ABSEPS, RLIST2a, IERR);


      ERRMAX = ABSERR;
      MAXERR = 1;
      AREA   = RESULT;
      ERRSUM = ABSERR;
      NRMAX  = 1;
      KTMIN  = 0;
      EXTRAP = false;
      NOEXT  = false;
      IERRO  = 0;
      IROFF1 = 0;
      IROFF2 = 0;
      IROFF3 = 0;
      KSGN   = -1;

      if(DRES >= (1.-50*EPMACH)*DEFABS) KSGN = 1;

    //======================================================  L90-for-begin
      for( LAST=2; LAST<=LIMIT; LAST++)  //L90
        {

        double A1 = ALIST[MAXERR];
        double B1 = 0.5E+00*(ALIST[MAXERR]+BLIST[MAXERR]);
        double A2 = B1;
        double B2 = BLIST[MAXERR];
        double ERLAST = ERRMAX;

        double AREA1,AREA2,ERROR1,ERROR2,DEFAB1,DEFAB2;
        QK15I(F,BOUN,INF,A1,B1,AREA1,ERROR1,RESABS,DEFAB1);
        QK15I(F,BOUN,INF,A2,B2,AREA2,ERROR2,RESABS,DEFAB2);


        double AREA12 = AREA1+AREA2;
        double ERRO12 = ERROR1+ERROR2;
        ERRSUM = ERRSUM+ERRO12-ERRMAX;
        AREA = AREA+AREA12-RLIST[MAXERR];
         if(DEFAB1 == ERROR1 || DEFAB2 == ERROR2) goto L15;
         if(fabs(RLIST[MAXERR]-AREA12) > 1.e-3*fabs(AREA12) || ERRO12 < 0.99*ERRMAX) goto L10;
         if(EXTRAP) IROFF2++;
         if(!EXTRAP) IROFF1++;

L10:     if(LAST > 10 && ERRO12 > ERRMAX) IROFF3++;

L15:    RLIST[MAXERR] = AREA1;
        RLIST[LAST]   = AREA2;

        ERRBND = qMax(EPSABS,EPSREL*fabs(AREA));

         if(IROFF1+IROFF2 >= 10 || IROFF3 >= 20) IER = 2;
         if(IROFF2 >= 5) IERRO = 3;

         if(LAST == LIMIT) IER = 1;

         if(qMax(fabs(A1),fabs(B2)) <= (1.+100.*EPMACH)* (fabs(A2)+1000.*UFLOW)) IER = 4;

         if(ERROR2 > ERROR1) goto L20;

        ALIST[LAST]   = A2;
        BLIST[MAXERR] = B1;
        BLIST[LAST]   = B2;
        ELIST[MAXERR] = ERROR1;
        ELIST[LAST]   = ERROR2;
        goto L30;

L20:    ALIST[MAXERR] = A2;
        ALIST[LAST]   = A1;
        BLIST[LAST]   = B1;
        RLIST[MAXERR] = AREA2;
        RLIST[LAST]   = AREA1;
        ELIST[MAXERR] = ERROR2;
        ELIST[LAST]   = ERROR1;

L30:
        //QPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX);
        //void QPSRT(int LIMIT, int LAST, int &MAXERR, double &ERMAX, double *ELIST, double *IORD, int &NRMAX)
        QPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX);

         if(ERRSUM <= ERRBND) goto L115;
         if(IER != 0)         goto L100;
         if(LAST == 2)        goto L80;
         if(NOEXT)            continue; //goto L90;

        ERLARG -= ERLAST;

         if(fabs(B1-A1) > SMALL) ERLARG += ERRO12;
         if(EXTRAP) goto L40;

         if(fabs(BLIST[MAXERR]-ALIST[MAXERR]) > SMALL) continue; //goto L90;

        EXTRAP = true;
        NRMAX = 2;

L40:     if(IERRO == 3 || ERLARG <= ERTEST) goto L60;

        ID = NRMAX;
        JUPBND = LAST;
        if(LAST > (2+LIMIT/2)) JUPBND = LIMIT+3-LAST;

        for(int K=ID; K<=JUPBND; K++)   //L50
                {
                MAXERR = IORD[NRMAX];
                ERRMAX = ELIST[MAXERR];
                if(fabs(BLIST[MAXERR]-ALIST[MAXERR]) > SMALL) continue; //goto L90;
                NRMAX++;
                }

L60:    EA(NEWFLG,AREA,LIMEXP,RESEPS,ABSEPS,RLIST2a,IERR);

        KTMIN++;

         if(KTMIN > 5 &&  ABSERR < 0.001*ERRSUM && LERR) IER = 5;

         if( ABSEPS >= ABSERR && LERR) goto L70;

        KTMIN = 0;
        ABSERR = ABSEPS;
        LERR = true;
        RESULT = RESEPS;
        CORREC = ERLARG;

        ERTEST = qMax(EPSABS,EPSREL*fabs(RESEPS));

        if( ABSERR <= ERTEST && LERR) goto L100;

L70:     if(RLIST2a[LIMEXP+3] == 1) NOEXT = true;
         if(IER == 5) goto L100;

        MAXERR = IORD[1];
        ERRMAX = ELIST[MAXERR];
        NRMAX = 1;
        EXTRAP = false;
        SMALL *= 0.5;
        ERLARG = ERRSUM;
        continue; //goto L90;

L80:    SMALL = 0.375;
        ERLARG = ERRSUM;
        ERTEST = ERRBND;
        EA(NEWFLG,AREA,LIMEXP,RESEPS,ABSEPS,RLIST2a,IERR);

//L90:
      }
      //======================================================  L90-for-end
//----------------


L100:  if(!LERR) goto L115;
       if((IER+IERRO) == 0) goto L110;
       if(IERRO == 3)       ABSERR += CORREC;
       if(IER == 0)         IER = 3;

       if(RESULT != 0 && AREA != 0)     goto L105;
       if(ABSERR > ERRSUM)              goto L115;
       if(AREA == 0)                    goto L130;

       goto L110;

L105:  if(ABSERR/fabs(RESULT) > ERRSUM/fabs(AREA)) goto L115;

L110:  if(KSGN == (-1) && qMax(fabs(RESULT),fabs(AREA)) <= DEFABS*0.01) goto L130;

       if(0.01 > RESULT/AREA  ||  RESULT/AREA > 100. || ERRSUM > fabs(AREA)) IER = 6;
       goto L130;

L115:
        RESULT = 0;
        for(int K=1; K<=LAST; K++)  RESULT += RLIST[K];


      ABSERR = ERRSUM;

L130:  NEVAL = 30*LAST-15;

       if(INF == 2) NEVAL *= 2;
       if(IER > 2)  IER--;

//L999:
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void QK15I( double (*F)(double x), double BOUND, int INF,
        double A, double B,
        double &RESULT, double &ABSERR,
        double &RESABS, double &RESASC)
{
/*
      double A,ABSC,ABSC1,ABSC2,ABSERR,B,BOUN,CENTR,
                DINF,r2mach,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,
            FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,TABSC1,TABSC2,
            UFLOW,WG,WGK,XGK;
      int INF,J,MIN0;
      EXTERNAL F;
          */
      double FV1[8],FV2[8];

      double  XGK[9] = {
                0,
               0.9914553711208126E+00,     0.9491079123427585E+00,
               0.8648644233597691E+00,     0.7415311855993944E+00,
               0.5860872354676911E+00,     0.4058451513773972E+00,
               0.2077849550078985E+00,     0.0000000000000000E+00};

      double WGK[9] = {
               0,
               0.2293532201052922E-01,     0.6309209262997855E-01,
               0.1047900103222502E+00,     0.1406532597155259E+00,
               0.1690047266392679E+00,     0.1903505780647854E+00,
               0.2044329400752989E+00,     0.2094821410847278E+00};

      double  WG[9] = {
                0,
               0.0000000000000000E+00,     0.1294849661688697E+00,
               0.0000000000000000E+00,     0.2797053914892767E+00,
               0.0000000000000000E+00,     0.3818300505051189E+00,
               0.0000000000000000E+00,     0.4179591836734694E+00};

      double EPMACH = r2mach(4);
      double UFLOW  = r2mach(1);

      double DINF = qMin(1,INF);

      double CENTR = 0.5*(A+B);
      double HLGTH = 0.5*(B-A);
      double TABSC1 = BOUND+DINF*(1.-CENTR)/CENTR;
      double FVAL1 = F(TABSC1);

      if(INF == 2) FVAL1 += F(-TABSC1);

      double FC = (FVAL1/CENTR)/CENTR;

      double RESG = WG[8] *FC;
      double RESK = WGK[8]*FC;
      RESABS = fabs(RESK);

      for(int J=1; J<=7; J++)
        {
        double ABSC  = HLGTH*XGK[J];
        double ABSC1 = CENTR-ABSC;
        double ABSC2 = CENTR+ABSC;

        double TABSC1 = BOUND+DINF*(1.-ABSC1)/ABSC1;
        double TABSC2 = BOUND+DINF*(1.-ABSC2)/ABSC2;
        double FVAL1 = F(TABSC1);
        double FVAL2 = F(TABSC2);

         if(INF == 2)   {
                        FVAL1 = FVAL1+F(-TABSC1);
                        FVAL2 = FVAL2+F(-TABSC2);
                        }

        FVAL1 = (FVAL1/ABSC1)/ABSC1;
        FVAL2 = (FVAL2/ABSC2)/ABSC2;
        FV1[J] = FVAL1;
        FV2[J] = FVAL2;
        double FSUM = FVAL1+FVAL2;
        RESG   = RESG + WG [J]*FSUM;
        RESK   = RESK + WGK[J]*FSUM;
        RESABS = RESABS+WGK[J]*(fabs(FVAL1)+fabs(FVAL2));
        }

      double RESKH  = RESK*0.5;
      RESASC = WGK[8]*fabs(FC-RESKH);

      for(int J=1; J<=7; J++) {
                RESASC += WGK[J]*(fabs(FV1[J]-RESKH)+fabs(FV2[J]-RESKH));
          }

      RESULT = RESK*HLGTH;
      RESASC *= HLGTH;
      RESABS *= HLGTH;

      ABSERR = fabs((RESK-RESG)*HLGTH);

       if(RESASC != 0.0 && ABSERR != 0.)
                        ABSERR = RESASC * qMin(1.,pow((200.*ABSERR/RESASC),1.5));

       if(RESABS > UFLOW/(50.*EPMACH) )
            ABSERR = qMax((EPMACH*50)*RESABS,ABSERR);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void QPSRT(int LIMIT, int LAST, int &MAXERR, double &ERMAX, double *ELIST, int *IORD, int &NRMAX)
{

    double ERRMAX,ERRMIN;
    int IBEG, IDO, JBND, JUPBN;
    int I, K; // I from cycle -- bad issue

      if(LAST <= 2)
               {
              IORD[1] = 1;
              IORD[2] = 2;
              goto L91;
              }

      ERRMAX = ELIST[MAXERR];

      if(NRMAX != 1)
                {
                IDO = NRMAX-1;

                for(I=1; I<=IDO; I++)  //L20
                        {
                        int ISUCC = IORD[NRMAX-1];

                        if(ERRMAX <= ELIST[ISUCC]) break;
                        IORD[NRMAX] = ISUCC;
                        NRMAX--;
                        }
                }


      JUPBN = LAST;

      if(LAST > LIMIT/2+2) JUPBN = LIMIT+3-LAST;
      ERRMIN = ELIST[LAST];

      JBND = JUPBN-1;
      IBEG = NRMAX+1;

      if(IBEG <= JBND)
           {
              for(I=IBEG; I<=JBND; I++) //*L40
                {
                int ISUCC = IORD[I];
                if(ERRMAX >= ELIST[ISUCC]) goto L61;
                IORD[I-1] = ISUCC;
                }
           }

      IORD[JBND]  = MAXERR;
      IORD[JUPBN] = LAST;
      goto L91;

L61:  IORD[I-1] = MAXERR;   ///  I from Cycle. It's not good

      K = JBND;

      for(int J=I; J<=JBND; J++)   //L70
        {
        int ISUCC = IORD[K];

        if(ERRMIN < ELIST[ISUCC]) {IORD[K+1] = LAST; goto L91;}

        IORD[K+1] = ISUCC;
        K--;
        }

      IORD[I] = LAST;

L91:  MAXERR = IORD[NRMAX];
      ERMAX  = ELIST[MAXERR];
}
//*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

double r2mach(int I)
{
//-------------------------------------------------------------
//     MACHINE CONSTANTS FOR THE IBM PC FAMILY (D. KAHANER NBS)
//-------------------------------------------------------------
      double RMACH[6] = {0,
                       1.18E-38,3.40E+38,0.595E-07,1.19E-07,0.30102999566};


/* Oleg Commented ---   should be checked later
      int SMALL[3];
      int LARGE[3];
      int RIGHT[3];
      int DIVER[3];
      int LOG10[3];

      EQUIVALENCE [RMACH[1]][SMALL[1]];
      EQUIVALENCE [RMACH[2]][LARGE[1]];
      EQUIVALENCE [RMACH[3]][RIGHT[1]];
      EQUIVALENCE [RMACH[4]][DIVER[1]];
      EQUIVALENCE [RMACH[5]][LOG10[1]];
*/
//-------------------------------------------------------------------

       if( I<1 || I>5)
                qDebug() << "r2mach -- I OUT OF BOUNDS" << 25 << 1 << 2;

      double r2mach = RMACH[I];
      return r2mach;
      }

//**********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void xerreur(char *MESSG,int NMESSG,int NERR,int LEVEL)    {XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.);}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void XERRWV(char *MESSG, int NMESSG, int NERR, int LEVEL, int NI,int I1,int I2, int NR,double R1, double R2)
{

      printf("NMESSG=%d;   NERR=%d",NMESSG,NERR);
      printf("message : %s",MESSG);

      if(NI == 2) printf("I1=%d;  I2=%d",I1,I2);
 else if(NI == 1) printf("I1=%d  ",I1);

        if(NR == 2)         printf("R1=%g;  R2=%g", R1,R2);
 else   if(NR == 1)         printf("R1=%g  ", R1);

if(fabs(LEVEL) < 2)return;
exit(2);
}
//*********************************************************************

