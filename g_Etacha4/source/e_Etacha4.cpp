#include "e_Etacha4.h"

#include <QDateTime>
#include <QtMath>
#include <stdio.h>
#include <QKeyEvent>
#include <QApplication>
#include <QMessageBox>
#include <QDebug>

#include "../source/e_Etacha4.h"
//#include "../win/e_mainwindow.h"
#include "e_Declare.h"
#include "e_AtomicShell.h"
#include "w_Stuff/liseStrcpyOS.h"
/*
//----------------------------------------------  utils
extern double Velocity_au(double E);
extern double E_to_Beta(double E);*/
extern double pow2(double par);
extern double pow2I(double par);
extern QString ElapsedTime(qint64 msec);
extern double EnergyResidueSP(double Zp, double Ap,  double Zt, double At, double E, double thick );
extern FILE *mfopen(const QString& filename, const char* operand);
/*
extern bool test_ESC(void);
//----------------------------------------------  F4
*/
extern void INI_F4();
extern int DONAUT(QWidget *w);
extern int f_num(int NCO);
extern int f_numP(int NCO);
extern int f_numPP(int NCO);


extern void ode
(
  void f ( double t, double y[], double yp[] ),
  int neqn,
  double y[],
  double &t,
  double tout,
  double relerr,
  double abserr,
  int &iflag,
  double work[],
  int iwork[]
);


extern int r8_rkf45
(
void f ( double t, double y[], double yp[] ), int neqn,
  double y[], double yp[], double *t, double tout, double *relerr,
  double abserr, int flag );


//void EqDif(int N,double  T,double *Y,void (*F)(int,double,double *,double*),double TOUT,int &MSTATE,
//                int NROOT,double EPS,double EWT,int MINT,
//                double *WORK,int LENW,int *IWORK, int LENIW,void (*G)(int,double,double *,double*));

//extern void F(int NEQ, double X, double *U,double *UP);
extern void F(double X, double *U, double *UP);

//------------------------------

extern void SecMean(double *Y);
extern double Auger(double *Y, double Zp, double *PR);
double PopMean(double* u, double *pr);

FILE* CreateEtaTxtFile(const char *filename, const QString &LinitialDir, const QString &LfileName, QDateTime  &s_time, int option);
int PrintBindingEnergies(FILE *f);
int PrintCrossSections(FILE *f);
int PrintCurrentDistribution(FILE *fil, int option, double  T, double EE);
int PrintCurrentExcelDistribution(FILE *fil, double  T, double EE);
int PrintFinalDistribution  (FILE *fil, int option, QDateTime &Time3, qint64 dif, double EE);

// LISE external
extern double gAb,gZb,gQb,gEnergy;
extern double gAt,gZt;
extern double gThick, gMinStepTarget,gMaxStepTarget;
extern double gO_EnergyL1,gO_EnergyL2,gO_dEdX1,gO_dEdX2, gZb, gAb,gQb,gZt, gAt,gEnergy, gDensity,
ProjectileVelocity,gUncertainRel,gUncertainAbs;

extern char BufRtf[1000];
extern atom_shell Zshell;
extern atom_shell Qshell;
extern int EtachaVersion;
extern int DifEqModel; //  0 - ODE, 1 - RKB45
extern int UseEloss;
extern const char *eta_filenames[];
extern bool GlobalBreak;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void CalcOlegSum(double *v);
void NormOlegSum(double *V);
double sumM[8];
extern int ODEsteps;

double QM, QF, dQF, currentThick;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
ETACHA::ETACHA(QWidget *parentI)
{
parent=parentI;

connect(this, SIGNAL(appendText(int,int)),  parent, SLOT(CM_appendText(int,int)));
connect(this, SIGNAL(appendShell(int,int)), parent, SLOT(CM_appendShell(int,int)));
connect(this, SIGNAL(updateStatusBar()), parent, SLOT(CM_updateStatusBar()));
connect(this, SIGNAL(updateGraph()), parent, SLOT(CM_updateGraph()));

};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ETACHA::init()
{
Zp = gZb;
Zt = gZt;
Ep  = gEnergy;
z1  = Zp;
z2  = Zt;

eRel  = gUncertainRel;
eAbs  = gUncertainAbs;

/*strcpy(BufRtf, " text test");
emit appendText(fsoBold,clRed);
emit appendShell(fsoItalic,clBlue);
*/
INI_F4();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

int ETACHA::Etacha(const QString &LinitialDir, const QString &LfileName, bool DebugMode)
{
//c*****************************************************************************
//c***     This version of ETACHA is for ions with up to 60 electrons        ***
//c***     First PC version of ETACHA was 18/12/91                           ***
//c***     First distributed PC version of ETACHA was 07/96                  ***
//c*****************************************************************************
//c
//c   05/2015 version going up to nmax=4 with partially correlated states
//c   Last corrections 09/2015
//c
//c   NOTES :    1) COA = 1-(1.75*qmin/vu)**2
//c              2) CDW calculations are not used any more, SEIK instead
//c              3) CDWEIS +SC & ASC calculations for ionisation
//c              4) SE + SC & ASC calculations for excitation
//c              5) Excitation to n >= 5 not added to ionization anymore
//c                 (See Donaut4.for)
//c
//c   Contact: J.P.ROZET  eMail: rozet@insp.jussieu.fr
//c            D. Vernhet eMail: dominique.vernhet@insp.jussieu.fr
//c            E. Lamour  eMail: lamour@insp.jussieu.fr
//c
//c       Calculates charge states of swift ions and their evolution as a fonction
//c    of traversed target thickness for a number of initial electrons =< 60.
//c    In a first step, only 1s,2s et 2p states are considered (as in oldest ETACHA),
//c    and 3s,3p and 3d state populations independently estimated.
//c    This program then solves first a set of 84 coupled differential equations :
//c    63 are for the evolution of the 63 "correlated states of the type Y(i,j,k),
//c    where i,j and k stand for the number of 2p,2s and 1s electrons
//c    (i between 0 et 6, j and k between 0 and 2), and 21 are for the evolution
//c    of 3s (3 equations), 3p (7 equations) and 3d (11 equations) substates.
//c    Differential equations are integrated using a Runge-Kutta type method and associated subroutines
//c    (see EQDIFF.for, an improvement of Adams code - a variable order predictor-corrector method)
//c
//c    In the actual version, the evolution of n<4 states is calculated in an improved way,
//c    by considering the  11*19=209 Y(n12,n3) type states, where n12 is the electron number in n=1 and 2,
//c    and n3 the electron number in n=3. For this calculation, actual (averadged) cross sections,
//c    functions of n=n12+n3 are used (see SecMean.for).
//c    This corresponds to a total of 84+209=293 equations (this is the ETACHA3 version).
//c
//c    A final evolution is calculated by considering first in an "independant" way the n=4 shell (33 equations),
//c    and then finally  the 29*33 Y(n123,n4)type states, where n123 is the total electron number in n<4 and n4 is
//c    the electron number in n=4.
//c
//c__________________________________________________________________
//c          We then have eventually 293+33+29*33=293+990=1283 equations.
//c__________________________________________________________________
//c
//c    Files to be included for building the exe file (15 files):
//c          -Etacha4.for (present file)
//c          -Auger4.for
//c          -Donaut4.for
//c          -EQDIF.for                                     1
//c          -F4.for
//c          -INTG.for
//c          -Pion.for (PWBA ionisation cross sections)
//c          -SecMean4.for
//c          -SEIK.for (capture cross sections in the SE approximation + REC)
//c          -senlm.for (excitation cross sections in teh SE approximation)
//c          -Sex2.for (PWBA excitation cross sections from n=2 and n=3)
//c          -Sexi.for (PWBA excitation cross sections from n=1)
//c          -Snl.for (PWBA intrashell excitation cross sections)
//c          -Tceis.for (CDW-EIS ionisation cross sections for n=1 and 2)
//c          -Zstop.for (Ziegler's SRIM derived stopping power)
//c    Files to be included as data files (4 files)
//c          -Etadon.dat   ==> etadon.etacha
//c          -SCOEFgas.95, SCOEF.95A, SCOEF.95B
//c    Output files :
//c          -ETA09.txt, ETA1019.txt, (ETA2029.txt, ETA3039.txt, ETA4049.txt and ETA5059.txt if usefull)
//c          -ETAPied.txt, POPMean.txt, seceff.txt
//c
//c
//c    ***********************************************************
    size_t len;
    char *bufFile;
    qint64 dif;
    QString dtime;
    QString ResultName;

const char *ShellYieldFormat =
        "<font size=-1><font color=\"blue\">T: %02.3f &gt;&nbsp; "
        "</font><font color=\"purple\"><b>k</b></font>: %4.2f &nbsp; "
        "<font color=\"DarkMagenta\"><b>2s</b></font>: %4.2f &nbsp; "
        "<font color=\"navy\"><b>2p</b></font>: %4.2f &nbsp; "
        "<font color=\"green\"><b>m</b></font>: %4.2f &nbsp; "
        "<font color=\"SeaGreen\"><b>n</b></font>: %4.2f &nbsp; "
        "<font color=\"gray\">&Sigma;e<sup>-</sup>: %5.3f</font></font>";
//c    ***********************************************************
//c    initializations and data input procedures
//c    ***********************************************************
//TMessage Message;
strcpy(BufRtf,"The Main module starts to calculate"); emit appendText(fsoUnderline,clBlue);

double PTF,/*YTOT,*/ tots, Sum7=0;    // YTOT is not used
double PT,EPM;
double Qin, /*dQin,*/ dQ, expectedT;
int iqmax; double vqmax;
double Ap  = gAb;
double Qp  = gQb;
double At  = gAt;
double Rho = gDensity;

int ReturnFlag =1;

      ep0 = gMinStepTarget;
      ep1 = gMaxStepTarget;
      EPM = gThick*1000.;
      ODEsteps = 0;
// not used int        ilgn=64;
// not used more int   iprt=1;

int NEQ;

switch(EtachaVersion)
               {
               case etacha_v23 : NEQ = 84; break;
               case etacha_v3  : NEQ = 293; break;
               case etacha_v34 : NEQ = 326; break;
               case etacha_v4  :
               case etacha_v45 : NEQ = 1283; break;
               };

//if(DifEqModel==1) NEQ++;  // temporary Oleg

// not used int int      icor=0;


//  not used:  ilp = ilgn-55;


// not used

//for(int  i=1; i<=38; i++)
//         if (Gcor[i] != 1.)  {icor=1; break;}


//  not used:  if (icor == 1) ilp -= 2;


double pas0=20.*ep0;
//Zp=z1;

double EE=Ep;
//Zt=z2;
//Qt not used double VC = sqrt(1.-pow2I((1.+EE/931.5)));
//double VU= 137.036 * VC;

double tc, as,  bs, stp;
//double stp1,astp,sstp;

//c....... slowing down coefficients ..........
if(UseEloss == 0)
        {
        tc=1.e9;
        as=0.;
        bs=EE;
        }
 else   {

        double d_dedx = gO_dEdX2    - gO_dEdX1;
        double d_ener = gO_EnergyL2 - gO_EnergyL1;
        if(d_dedx == 0) d_dedx = 1;
        if(d_ener == 0) d_ener = 100;

        as  = d_dedx/(d_ener*Ap);
        bs  = (gO_dEdX2*gO_EnergyL1-gO_dEdX1*gO_EnergyL2)/d_dedx;
        stp =  gO_dEdX1 +  (EE-gO_EnergyL1)*d_dedx / d_ener;

        if(stp < 0.) {
                sprintf(BufRtf,"    stp=%f", stp);                                                       emit appendText(fsoBold,clRed);
                sprintf(BufRtf,"<<<--  problem with stopping power values, end of calculation <<<--  "); emit appendText(fsoBold,clRed);
                return -2;
                }

        double stp1 = log10(Ap*EE/stp);
        double astp = pow(10.,int(stp1));
        if (stp1 < 0.) astp/=10.;


        double sstp = Ap*EE/(stp*astp);
              if (sstp <= 1.5)  tc=2.5*astp;
         else if (sstp <= 3.5)  tc=5.0*astp;
         else if (sstp < 7.)    tc=10.*astp;
         else                   tc=25.*astp;
         }

double ddq = 1;

      if (Zp <=  6)  {  ddq=8.;   if (eRel > 1e-4) ddq/=2.;  }
 else if (Zp <= 20)  {  ddq=4.;   if (eRel > 1e-4) ddq/=2.;  }
 else if (Zp <= 40)  {  ddq=2.;   if (eRel > 1e-4) ddq/=2.;  }

//------------------------------------------------------
sprintf(BufRtf,"  initialization of populations");
emit appendText(fsoItalic,clBlack);

                 //c    ***** Y(n) initialization *****************************
for(int n=0; n<=1284; n++)  { Y[n]=0.;   YX[n]=0.; }

                //c    ***** number of electrons in shells of incident ion ***********
int NEL=Zp-Qp;

if (NEL > 60) {
               sprintf(BufRtf,"<br>"
              "**************************************************"
              " too much electrons ..."
              " try again with ZP-Q < 61"
              "**************************************************<br>");
              emit appendText(fsoBold,clRed);
              return -3;
              }

//Qt not used int nkl  = Qshell.nKL;
//Qt not used int nm   = Qshell.nM;
int nklm = Qshell.nKLM;
int n4l  = Qshell.nN;
int n1234 = 327+nklm+29*n4l;

int Ic1=Qshell.n_2p()*100 + Qshell.n_2s()*10+Qshell.n_1s();
int Ic2=Qshell.nM *100 + Qshell.nKL;

int n1=f_num (Ic1);
int n2=f_numP(Ic2);


sprintf(BufRtf,"initialization of populations for actual charge<br>");

Y[n1]=1.;
Y[64+Qshell.n_3s()]=1.;
Y[67+Qshell.n_3p()]=1.;
Y[74+Qshell.n_3d()]=1.;

if(EtachaVersion >= etacha_v3 ) Y[n2]=1.;
if(EtachaVersion >= etacha_v34) Y[294+n4l]=1.;
if(EtachaVersion >= etacha_v4 ) Y[n1234]  =1.;


//c      ***************************************************************
//c      Open and write headings of output files
//c      ***************************************************************
QDateTime Time1,Time2,Time3;
QDateTime Time0 = QDateTime::currentDateTime();

FILE *f09=0,*f19=0, *f29=0, *f39=0, *f49=0, *f59=0, *fPied=0, *fMean=0, *fSEff=0, *fExcel=0;

            f09   = CreateEtaTxtFile(eta_filenames[en_f09], LinitialDir, LfileName, Time0, 0);     //10
            f19   = CreateEtaTxtFile(eta_filenames[en_f19], LinitialDir, LfileName, Time0, 1);     //11
if(Zp>=20)  f29   = CreateEtaTxtFile(eta_filenames[en_f29], LinitialDir, LfileName, Time0, 2);     //12
if(Zp>=30)  f39   = CreateEtaTxtFile(eta_filenames[en_f39], LinitialDir, LfileName, Time0, 3);     //13
if(Zp>=40)  f49   = CreateEtaTxtFile(eta_filenames[en_f49], LinitialDir, LfileName, Time0, 4);     //20
if(Zp>=50) {f59   = CreateEtaTxtFile(eta_filenames[en_f59], LinitialDir, LfileName, Time0, 5);}     //21
            fPied = CreateEtaTxtFile(eta_filenames[en_fPied], LinitialDir, LfileName, Time0, 6);     //14
            fMean = CreateEtaTxtFile(eta_filenames[en_fMean], LinitialDir, LfileName, Time0, 7);     //15
            fSEff = CreateEtaTxtFile(eta_filenames[en_fSEff], LinitialDir, LfileName, Time0, 8);     //16
            fExcel = CreateEtaTxtFile(eta_filenames[en_fExcel], LinitialDir, LfileName, Time0, 9);     //16

//c    *******************************************************
//c    calculates cross sections (incident ion) in reduced units (ug/cm2)-1
//c    *******************************************************

//c    ******************************************************************
//c    *** radiative and Auger "cross sections" (units 10e-20 cm²) ***

//Oleg  because it was twice --->        PrintBindingEnergies(fSEff);

      QM = PopMean(Y, PR);
      CHGT(EE,Zp,At,Rho,as,EE,0,QM,fSEff);   // bs = EE and  dt=0 -- to omit energy loss for the 1st step,
      SecMean(Y);

      sprintf(BufRtf,"\"T\" - thicknenss in mg/cm<sup>2</sup><br>");     emit appendShell(fsoItalic,clOlive);
      sprintf(BufRtf,ShellYieldFormat,0,y1s,y2s,y2p,yM,yN,y1s+y2s+y2p+yM+yN);
      emit appendShell(0,clBlack);


      Time1 = QDateTime::currentDateTime();
      QString time1text=Time1.time().toString("hh:mm:ss");
      sprintf(BufRtf,"       begin time = %s<br>", time1text.toStdString().c_str()); emit appendText(fsoItalic,clBlack);


      int iboucle=0;
                                //not used more       int MINT=3;
                                //not used more       int NROOT=0;
      double T=0.;
      double tout=ep0;
      double toutc=0.;
      double QM0=QM;
      int Ocounter=0;
//c    *******************************************************
//c    start integration
//c    *******************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw  etacha cycle start
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw  etacha cycle start
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw  etacha cycle start
      int MSTATE=1;

      double EPS=eRel;
      double EWT=eAbs;
      int Flag3_counter=0;
      int Flag4_counter=0;

      if(fMean)  fprintf(fMean,"\n%9.2f %7.4f %7.4f %7.4f %8.4f %8.4f %8.4f",
                        T,y1s,y2s,y2p,yM,yN,QM) ;
//      if(fMean)  fprintf(fMean,"\n%9.2f %7.4f %7.4f %7.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f",
//                        T,y1s,y2s,y2p,yM,yN,QM,Qin,QF,PTF) ;

L100:

      double dtout=(tout-toutc)/tc;
      double DQM=fabs((QM-QM0)*ddq);

      if(dtout >= 1. || DQM >= 1.)
//      if(false)
                {
                double dt=tout-toutc;

                QM = PopMean(Y, PR);

                EE=CHGT(EE,Zp,At,Rho,as,bs,dt,QM,fSEff);

                if(EE<=0)
                        {
                        sprintf(BufRtf," &nbsp; &nbsp; &nbsp; &nbsp;     &lt;&lt;&lt;&lt;  Energy < 0 &gt;&gt;&gt;&gt;");
                        emit appendText(fsoBold,clRed);
                        ReturnFlag=-3;
                        goto L5000;
                        }

                SecMean(Y);

                toutc=tout;
                QM0=QM;
                Time2=QDateTime::currentDateTime();
                dif = Time1.msecsTo(Time2);
                dtime = ElapsedTime(dif);
                sprintf(BufRtf,"&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
                "<font size=-1> CHGT:  %s</font>",dtime.toStdString().c_str());
                emit appendText(fsoItalic,clBlack);
                }

                              // oleg       never used more       int LENW=1670000;
                              // oleg       never used more       int LENIW=1310;
                              // oleg    X- never used       double X=T;
                              // oleg       never used      int nequ=0;

                              /// temporary test  start
                              //    F(0, 0, Y, YX);
                              /// temporary test  end

                                //not used      YPMax=0.;


                        // Oleg
                        //      EqDif(NEQ,T,Y,F,tout,MSTATE,
                        //            NROOT,EPS,EWT,MINT,
                        //            WORK,LENW,IWORK,LENIW,F);
                        //

CalcOlegSum(Y);
MSTATE=1;

expectedT  = tout;
if(DifEqModel==0)
                {
                //    * NEQN, the number of equations to be integrated;
                //    * Y(1:NEQN), the vector of initial conditions;
                //    * T, the starting point of integration;
                //    * TOUT, the point at which a solution is desired;
                // ---Input/output, double &T, the current value of the independent    variable.
                // ---Input, double TOUT, the desired value of T on output.
                //    * RELERR, ABSERR, the relative and absolute local error tolerances;
                //    * IFLAG, an indicator to initialize the code.  Normal input
                //      is +1.  The user should set IFLAG = -1 only if it is
                //      impossible to continue the integration beyond TOUT.
                //     The subroutine integrates from T to TOUT.
                ode(F, NEQ, &Y[1], T, expectedT, EPS, EWT, MSTATE, WORK, IWORK);
//                ode(F, NEQ, Y, T, tout, EPS, EWT, MSTATE, WORK, IWORK);
                }
else            {
                //     Typically the subroutine is used to integrate from T to TOUT but it
                //    can be used as a one-step integrator to advance the solution a
                //    single step in the direction of TOUT.


                MSTATE = r8_rkf45 (F, NEQ, &Y[1], YX, &T, expectedT, &EPS, EWT, MSTATE);
// Oleg  temp               MSTATE = r8_rkf45 (F, NEQ, Y, YX, &T, tout, &EPS, EWT, MSTATE);
//                int r8_rkf45 ( void f ( double t, double y[], double yp[] ), int neqn,
//                                  double y[], double yp[], double *t, double tout, double *relerr,
//                                  double abserr, int flag )
                //    * initialize the parameters:
                //      NEQN, Y(1:NEQN), T, TOUT, RELERR, ABSERR, FLAG.
                //      In particular, T should initially be the starting point for integration,
                //      Y should be the value of the initial conditions, and FLAG should
                //      normally be +1.
                //
                //    Normally, the user only sets the value of FLAG before the first call, and
                //    thereafter, the program manages the value.  On the first call, FLAG should
                //    normally be +1 (or -1 for single step mode.)  On normal return, FLAG will
                //    have been reset by the program to the value of 2 (or -2 in single
                //    step mode), and the user can continue to call the routine with that
                //    value of FLAG.
                //
                //    (When the input magnitude of FLAG is 1, this indicates to the program
                //    that it is necessary to do some initialization work.  An input magnitude
                //    of 2 lets the program know that that initialization can be skipped,
                //    and that useful information was computed earlier.)
                //
                //    The routine returns with all the information needed to continue
                //    the integration.  If the integration reached TOUT, the user need only
                //    define a new TOUT and call again.  In the one-step integrator
                //    mode, returning with FLAG = -2, the user must keep in mind that
                //    each step taken is in the direction of the current TOUT.  Upon
                //    reaching TOUT, indicated by the output value of FLAG switching to 2,
                //    the user must define a new TOUT and reset FLAG to -2 to continue
                //    in the one-step integrator mode.
                //

                }

Ocounter++;
if(DebugMode)
    {
    sprintf(BufRtf,"debug &gt;&gt; number of steps = %d &nbsp; "
                   "counter = %d; &nbsp; MSTATE = %d" , ODEsteps, Ocounter, MSTATE);
    emit appendText(0,clGray);
    }

sprintf(BufRtf,"number of steps = %d", ODEsteps); emit updateStatusBar();

if(T !=  expectedT )
                {
                sprintf(BufRtf,"Treal = %.2f    Texpected = %.2f", T,expectedT);
                emit appendText(fsoItalic,clRed);
                }


//c    *************** IDID=1 as long as T <= tout **********
       if (MSTATE > 2) {
                        message(MSTATE);

                        if( (MSTATE==3 && Flag3_counter==0)  || MSTATE!=3)  message(MSTATE);

                               if(MSTATE==3) { Flag3_counter++; /*MSTATE=2;*/}
                       else    if(MSTATE==4) { Flag4_counter++; /*MSTATE=2;*/}
                       else         {ReturnFlag=-1; goto L5000; }
                       }

CalcOlegSum(Y); Sum7 = sumM[7];

iboucle = IWORK[3];       // not used AVGORD  = WORK[3];

for(int I=0; I<=1283; I++)  if(Y[I] < 0.)  Y[I]=0.;

CalcOlegSum(Y);  NormOlegSum(Y);

for(int I=0; I<=1283; I++)  YX[I]=Y[I]*100.;


                       //c    **************** end of one step calculation ************
                       // oleg    X- never used   X=T;
//c    ***********************************************
//c    charge state probabilities calculation
//c    ***********************************************

      QM = PopMean(Y, PR);

//*  Oleg 09/21/21 v.4.3.14
         PTF=Qin=0; //dQin=0.;

//---------------------------
      if(EtachaVersion >= etacha_v4)                   //     charge states before autoionization
        {
          for(int N=0; N<=61; N++)  PRF[N]=0.;         // 123&4(KLM&N) 28+32+1

          for(int N1=1; N1<=29; N1++)                  // it is only for Etacha4              29 = 28 electron orbitals + 1
            for(int N4=1; N4<=33; N4++)                // it is only for Etacha4              33 = 4s(2) + 4p(6) + 4d(10) + 4f(14) + 1
              {
                int I   = N1-1;
                int N   = N4-1;
                int INM = 100*N+I;
                int NN  = f_numPP(INM);
                int Nel   = I+N4;                // number of electrons
                PRF[Nel] += Y[NN];
              }

          for(int M=1; M<=61; M++)        { double RM=M-1.;    PTF += PRF[M];     dQ = Zp-RM;  Qin  +=         dQ*PRF[M]; }
     //     for(int M=1; M<=61; M++)        { double RM=M-1.;                       dQ = Zp-RM; dQin += pow2(Qin-dQ)*PRF[M]; }
     //     dQin = sqrt(dQin);
        }
      //---------------------------

                                //c    ***********************************************
                                //c    re-calculate cross sections for actual charge state
                                //c    ***********************************************

      SecMean(Y);
// put debug lines
                                //c-------- K shell autoionization and K+L  populations --------
      /*YTOT = */Auger(Y,Zp,PR);
// should be again compared to  original code


//---------------------------


       PT=QF=dQF=0.;

      for(int M=1; M<=62; M++)  {double RM=M-1.;  dQ = Zp-RM;   QF +=         dQ *PR[M]; }
      for(int M=1; M<=62; M++)  {double RM=M-1.;  dQ = Zp-RM;  dQF += pow2(QF-dQ)*PR[M];  PR[M] *=100.;   PT += PR[M];}
      dQF = sqrt(dQF);

//---------------------------



//  only for debug purpose -- start     ------------------
      for(int K=0; K<=2; K++) P1s[K]=0;
      for(int J=0; J<=2; J++) P2s[J]=0;
      for(int I=0; I<=6; I++) P2p[I]=0;

      for(int K=0; K<=2; K++)
        for(int J=0; J<=2; J++)
          for(int I=0; I<=6; I++)
                  {
                  int L = 100*I+10*J+K;
                  int N = f_num(L);
                  P1s[K] += Y[N];
                  P2s[J] += Y[N];
                  P2p[I] += Y[N];
                  }
//  only for debug purpose -- stop ---------------------------------


//c    ********* prints probabilities at each step *****

                          PrintCurrentDistribution(f09, 0, T, EE);
                          PrintCurrentDistribution(f19, 1, T, EE);
      if (Zp >= 20)       PrintCurrentDistribution(f29, 2, T, EE);
      if (Zp >= 30)       PrintCurrentDistribution(f39, 3, T, EE);
      if (Zp >= 40)       PrintCurrentDistribution(f49, 4, T, EE);
      if (Zp >= 50)       PrintCurrentDistribution(f59, 5, T, EE);

     PrintCurrentExcelDistribution(fExcel,T,EE);


      tots=YX[1]+YX[2]+YX[4]+YX[10]+YX[3]+YX[5]+YX[11]+YX[6]+YX[12];

 if(fPied)  fprintf(fPied,"\n%9.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f",
                T,YX[1],YX[2],YX[4],YX[10],YX[3],YX[5],YX[11], YX[6],YX[12],tots);

 if(fMean)  fprintf(fMean,"\n%9.2f %7.4f %7.4f %7.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f",
                        T,y1s,y2s,y2p,yM,yN,QM,Qin,QF,PT) ;


currentThick=T/1000;

iqmax=0; vqmax=0;
for(int M=1; M<=62; M++)  if(vqmax<PR[M]) {vqmax=PR[M]; iqmax=M;}
iqmax = Zp - (iqmax-1);

//sprintf(BufRtf,"achieved thickness=%.3f mg/cm2   iterations for this step=%d", T/1000,iboucle);
sprintf(BufRtf,"achieved T=%.3f mg/cm<sup>2</sup>  &nbsp;  iter=%d &nbsp; &nbsp; "
                "<font color=\"green\">q<sub> prob</sub>=<b>%.2f</b></font> &nbsp; "
                "<font color=\"blue\">&lt;q&gt;=%.2f(%.2f)</font> &nbsp; "
                "q<sub> max</sub>=%i &nbsp; &nbsp; E=%.3f &nbsp; <font size=-1>&delta;&Sigma;=%.3f</font>",
                T/1000,iboucle,QM,QF,dQF,iqmax,EE,Sum7);
emit appendText(0,clBlack);

sprintf(BufRtf,ShellYieldFormat, T/1000.,y1s,y2s,y2p,yM,yN,y1s+y2s+y2p+yM+yN);


emit appendShell(0,clBlack);
emit updateGraph();

//Qt not used       iboucle=0;

//-------------------------------------------
// check for Escape
//
//

if(GlobalBreak)
        {
        int flagVer=QMessageBox::question(parent,"The \"Cancel\" button has been pressed", "Do you want to cancel calculations?");
        if(flagVer==QMessageBox::Yes)
                {
                sprintf(BufRtf," &nbsp; &nbsp; &nbsp; &nbsp;     &lt;&lt;&lt;&lt;  User's break!! &gt;&gt;&gt;&gt;");
                emit appendText(fsoBold,clRed);
                ReturnFlag=-2;
                goto L5000;
                }
        else  GlobalBreak=false;
        }

//c********* calculation goes on as long as T <= EPM *****
if (T < EPM)
      {                         //c*************** (pas0=20*ep0) ******************
      if (T < pas0  ) { tout = T + ep0;               if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T <= 0.2  ) { tout = T + ep0;               if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 1.4999) { tout = T + qMin(ep1,0.05);    if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 1.9999) { tout = T + qMin(ep1,0.1);     if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 4.9999) { tout = T + qMin(ep1,0.2);     if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 14.999) { tout = T + qMin(ep1,0.5);     if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 19.99 ) { tout = T + qMin(ep1,1.);      if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 49.99 ) { tout = T + qMin(ep1,2.);      if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 149.99) { tout = T + qMin(ep1,5.);      if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 199.99) { tout = T + qMin(ep1,10.);     if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 499.99) { tout = T + qMin(ep1,20.);     if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 1499.9) { tout = T + qMin(ep1,50.);     if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 1999. ) { tout = T + qMin(ep1,100.);    if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 4999. ) { tout = T + qMin(ep1,200.);    if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 14999.) { tout = T + qMin(ep1,500.);    if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 19999.) { tout = T + qMin(ep1,1000.);   if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else if (T < 49999.) { tout = T + qMin(ep1,2000.);   if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
 else                 { tout = T + qMin(ep1,5000.);   if (tout > EPM) tout=EPM;  /* iboucle=0.; */  goto L100; }
      }

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw  etacha cycle end
 //c ************* end of output files **************************

Time3=QDateTime::currentDateTime();
dif = Time1.msecsTo(Time3);
dtime = ElapsedTime(dif);

      sprintf(BufRtf,Star89);                                              emit appendText(0,clBlack);
      sprintf(BufRtf," &nbsp; Finished at &nbsp; %s",
              Time3.time().toString().toStdString().c_str());              emit appendText(fsoNone,clBlack);
      sprintf(BufRtf," &nbsp; %s",dtime.toStdString().c_str());            emit appendText(fsoNone,clBlack);
      sprintf(BufRtf,Star89);                                              emit appendText(0,clBlack);

if(tc < EPM) {
         double dt=tout-toutc;
         if (dt > 0. && UseEloss) {
                  double Ef=EE+(bs-EE)*(1.-exp(-0.001*tc*as));
                  EE=Ef;
                }
      }

sprintf(BufRtf," &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  "
                "&nbsp;  &nbsp;  Final energy : %.3f (MeV/u)", EE);  emit appendText(fsoItalic,clRed);


PrintFinalDistribution(f09, 0, Time3, dif, EE);
PrintFinalDistribution(f19, 1, Time3, dif, EE);
PrintFinalDistribution(f29, 2, Time3, dif, EE);
PrintFinalDistribution(f39, 3, Time3, dif, EE);
PrintFinalDistribution(f49, 4, Time3, dif, EE);
PrintFinalDistribution(f59, 5, Time3, dif, EE);

PrintFinalDistribution(fPied, 6, Time3, dif, EE);
PrintFinalDistribution(fMean, 7, Time3, dif, EE);
PrintFinalDistribution(fSEff, 8, Time3, dif, EE);

ResultName = LinitialDir+LfileName;
len = ResultName.size()+1;
bufFile = new char[len];
strcpyL(bufFile, len, ResultName.toStdString().c_str());


      sprintf(BufRtf,"<br>output data in files:");                                         emit appendText(fsoItalic,clBlue);
                       sprintf(BufRtf,"<font size=-1>  00 to 09 EE- charge states  in &nbsp; %s_%s</font>",bufFile,eta_filenames[en_f09]);     emit appendText(0,clNavy);
                       sprintf(BufRtf,"<font size=-1>  10 to 19 EE- charge states  in &nbsp; %s_%s</font>",bufFile,eta_filenames[en_f19]);     emit appendText(0,clNavy);
       if (Zp >= 20) { sprintf(BufRtf,"<font size=-1>  20 to 29 EE- charge states  in &nbsp; %s_%s</font>",bufFile,eta_filenames[en_f29]);     emit appendText(0,clNavy);   }
       if (Zp >= 30) { sprintf(BufRtf,"<font size=-1>  30 to 39 EE- charge states  in &nbsp; %s_%s</font>",bufFile,eta_filenames[en_f39]);     emit appendText(0,clNavy);   }
       if (Zp >= 40) { sprintf(BufRtf,"<font size=-1>  40 to 49 EE- charge states  in &nbsp; %s_%s</font>",bufFile,eta_filenames[en_f49]);     emit appendText(0,clNavy);   }
       if (Zp >= 50) { sprintf(BufRtf,"<font size=-1>  50 to 59 EE- charge states  in &nbsp; %s_%s</font>",bufFile,eta_filenames[en_f59]);     emit appendText(0,clNavy);   }

      sprintf(BufRtf,"<font size=-1> bare,1s,2s,2p,1s2,1s2s,1s2p,1s2 2s,1s2+2p ions and sum of these in %s_ETAPIED.txt</font>",bufFile); emit appendText(0,clNavy);
      sprintf(BufRtf,"<font size=-1> mean 1s,2s,2p,3s,3p and 3d populations in %s_POPMEAN.txt</font>",bufFile); emit appendText(0,clNavy);
      sprintf(BufRtf,"<br>WARNING! Next calculation will overwrite these files. Consider saving or renaming these results !"); emit appendText(fsoItalic,clOlive);


L5000:


if(f09) {fclose(f09); f09=nullptr;}
if(f19) {fclose(f19); f19=nullptr;}
if(f29) {fclose(f29); f29=nullptr;}
if(f39) {fclose(f39); f39=nullptr;}
if(f49) {fclose(f49); f49=nullptr;}
if(f59) {fclose(f59); f59=nullptr;}
if(fPied) {fclose(fPied); fPied=nullptr;}
if(fMean) {fclose(fMean); fMean=nullptr;}
if(fSEff) {fclose(fSEff); fSEff=nullptr;}
if(fExcel){fclose(fExcel); fExcel=nullptr;}


// Last message

if(ReturnFlag>0) {
        sprintf(BufRtf,"<BR> FINAL achieved &gt;&gt; T=%.3f mg/cm<sup>2</sup>  &nbsp; <b>&lt;Q&gt;=%.3f &nbsp;"
        " dQ=%.3f</b> &nbsp; E=%.3f &nbsp;  dSum=%.3f", T/1000,QF,dQF,EE,Sum7);
        emit appendText(0,clBlue);
        }

emit updateGraph();
return ReturnFlag;
}
//c    ***************************************************
//c    *********** END of etacha *************************
//c    ***************************************************

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double ETACHA::CHGT(double EE, double Zp, double At,double Rho,
            double as, double bs, double tc,
            double QM, FILE *f)
{

       EE += (bs-EE)*(1.-exp(-0.001*tc*as));
//double E2 = EnergyResidueSP(gZb, gAb, gZt, gAt, EE, tc*0.001);
//qDebug() << "energy " << EE << E2;

      CSEC(EE,2);

//      sprintf(BufRtf,"CHGT => QM=%.2f  EE=%.3f",QM,EE);   emit appendShell(0,clBlack);
      if(f) fprintf(f,"QM=%.4f     EE=%.4f",QM,EE);

       PrintCrossSections(f);
       PrintBindingEnergies(f);

// not used       double dum=QM;

       int iZp=Zp;
       double ok = omk[iZp];
       double ol = oml[iZp];

       double conv = 6.0221408e-3 / At;                   // number of atoms per ug/cm2 / 10e20

      double VC = sqrt(1.-pow2I((1.+EE/931.5)));

      double tem = pow((zk*tetak),4)*At/(Rho*VC);
       Rad   = 3.461E-6    *tem /2.;
       Rad3s = 0.034891E-6 *tem /6.;
double Rad3p = 1.0301E-6   *tem /2.;
       Rad3d = 0.3544E-6   *tem /6.;
       Rad4  = 0.0454E-6   *tem /18.;

       int n02p = Zshell.n_2p();
       int nl0  = Zshell.nL;
       int nm0  = Zshell.nM;
       int n03s = Zshell.n_3s();
       int n03p = Zshell.n_3p();
       int n03d = Zshell.n_3d();

       if (n02p >= 1)  AKLL = Rad*(1./ok-1.)*n02p/(nl0*(nl0-1));
       else            AKLL = 7.6E-2*At/(Rho*VC);

       if (nm0 >= 1)   AKLM = Rad3p*0.882*(1./ok-1.)*n02p/(nl0*nm0);
       else            AKLM = 0.3*At/(Rho*VC);

       if (nm0 >= 2)   ALMM=(Rad3s*n03s+0.118*Rad3p*n03p+Rad3d*n03d)*(1./ol-1.)/(nm0*(nm0-1));
       else            ALMM=0.05*At/(Rho*VC);


      AM4=3.*ALMM;

      for(int I=1; I<=34; I++) GSEC[I] *= conv;    // probabilty for 1 ug/cm2 , where  CS in 1e-20 cm2


      for(int I=1; I<=34; I++) GSEC[I] /= Coef_Electron_Holes[I];

      Rad     *= conv;
      Rad3s   *= conv;
      Rad3p1  = Rad3p*0.882*conv;
      Rad3p2  = Rad3p*0.118*conv;
      Rad3d   *= conv;
      Rad4    *= conv;

      Rad3s   += GSEC[25];
      Rad3p1  += GSEC[18];
      Rad3p2  += GSEC[22];
      Rad3d   += GSEC[27];

      AKLL *= conv;
      AKLM *= conv;
      ALMM *= conv;
      AM4  *= conv;

      return EE;
      }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

double PopMean(double *U, double *PR)
{
extern int f_II(int num);
extern int f_JJ(int I, int num);
extern int f_KK(int I,int J,int num);

      y3s = U[65]+2.*U[66];                                          // U[64] --> 0 electrons @ 3s
      y3p = U[68]+2.*U[69]+3.*U[70]+4.*U[71]+5.*U[72]+6.*U[73];      // U[67] --> 0 electrons @ 3p
      y3d = U[75]+ 2.*U[76]+ 3.*U[77]+ 4.*U[78]+  5.*U[79]+       // U[74] --> 0 electrons @ 3d
          6.*U[80]+ 7.*U[81]+ 8.*U[82]+ 9.*U[83]+ 10.*U[84];

      yM=y3s+y3p+y3d;

      yN=0.;
      yNn=0.;
      for(int n=1; n<=32; n++)
         {
         double rn = n;
         double t = rn*U[294+n];

                     yN  += t;
         if (n >= 2) yNn += t*(rn-1.);
         }
// -------- see F.for (ALMM) ------------------------------
      yMp=0.;
      yMm=0.;

      for(int I=0; I<=2; I++)               // 3s
         for(int J=0; J<=6; J++)            // 3p
           for(int K=0; K<=10; K++)         // 3d
              {
              int N=I+J+K;                           // n of e- @ 3
              double t = U[64+I]*U[67+J]*U[74+K]*N;

                          yMp += t;                 //
              if (N >= 2) yMm += t*(N-1);
              }

      yM1m=0.;

      for(int I=0; I<=2; I++)                // 3s
        for(int J=0; J<=6; J++)              // 3p
                {
                int N=I+J-1;
                if(N >= 2) yM1m += U[64+I]*U[67+J]*(N-1);
                }

      yM2m=U[76]+2.*U[77]+3.*U[78]+4.*U[79]+5.*U[80]+6.*U[81]+7.*U[82]+8.*U[83]+9.*U[84];  //only 3d

//------------------------------------------------------
//.....mean values for 1s, 2s, 2p ......
      y1s=0.;
      y2s=0.;
      y2p=0.;


      for(int N=1; N<=63; N++)
          {
          int I=f_II(N);       // 2p
          int J=f_JJ(I,N);     // 2s
          int K=f_KK(I,J,N);   // 1s

          y1s += K*U[N];
          y2s += J*U[N];
          y2p += I*U[N];
          }

      yL=y2s+y2p;

//.....mean values for 1s², 2s², 2p² (in case is needed) ......
      yKm  = 0.;
      yL1m = 0.;
      yL2m = 0.;

      for(int N=1; N<=63; N++)
        {
        int I = f_II(N);
        int J = f_JJ(I,N);
        int K = f_KK(I,J,N);

        if(K >= 2) yKm  += U[N];           // 1s2
        if(J >= 2) yL1m += U[N];           // 2s2
        if(I >= 2) yL2m += (I-1)*U[N];     // 2p2
        }

//***********************************************************  Oleg - PR for old versions  START
      if(EtachaVersion < etacha_v4)                   //  modified  09/21/2021
        {
          for(int N=0; N<=62; N++)   PR[N]=0.;
          //  check e_SecMean4
          for(int N12=1; N12<=63; N12++)
            for(int L=64; L<=66; L++)           // 3s
              for(int LL=67; LL<=73; LL++)       // 3p
                for(int LLL=74; LLL<=84; LLL++)  // 3d
                  {
                    int I = f_II(N12);      //  2p
                    int J = f_JJ(I,N12);    //  2s
                    int K = f_KK(I,J,N12);  //  1s

                    int ne12  = I+J+K              ;
                    int ne3   = L+LL+LLL-205       ;
                    int ne    = ne12 + ne3 + 1     ;

                    if(EtachaVersion >= etacha_v4) ne+=2;      //  modified  09/21/2021

                    double PT = Y[N12]*Y[L]*Y[LL]*Y[LLL];
                    PR[ne] += PT;
                  }
        }
//***********************************************************  Oleg - PR for old versions  START

double QM = Zp-(y1s+yL+yM+yN);   //  mmore probable
return QM;
}
//c*********************************************************************
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void ETACHA::message(int idid)
{
sprintf(BufRtf,"\nMSTATE=%d problem integrating with QDIFF\n", idid);
emit appendText(fsoUnderline,clRed );
if(DifEqModel==0)
{
     if (idid == 3)
       sprintf(BufRtf," 3, integration did not reach TOUT because the error tolerances were too small. \n"
                      "But RELERR and ABSERR were increased appropriately for continuing");

else if (idid == 4)
       sprintf(BufRtf," 4, integration did not reach TOUT because more than 500 steps were taken");

else if (idid == 5)
       sprintf(BufRtf,"5. integration did not reach TOUT because the equations appear to be stiff\n"
                          "try to use max step value and decrease absolute uncertainty");

else if (idid == 6)
       sprintf(BufRtf,"6, invalid input parameters (fatal error)");

//if (idid == 3) {
//      sprintf(BufRtf," more than 1000 iterations have been atempted ...\n"
//                     " try reducing maximum step size\n");
//      }
//else  if (idid == 4) {
//      sprintf(BufRtf," you are probably asking too much accuracy ...\n"
//                     " try again with larger uncertainties\n");
//      }
//if (idid > 4) {
//      sprintf(BufRtf," for some reason, the problem is very stiff and\n"
//                     " cannot be solved with the present integration routine\n");
//      }
}
else
{
     if (idid == 3)
      sprintf(BufRtf,"3, integration was not completed because the input value of RELERR, the\n"
         "relative error tolerance, was too small.  RELERR has been increased\n"
         "appropriately for continuing.  If the user accepts the output value of\n"
         "RELERR, then simply reset FLAG to 2 and continue.\n");
else if (idid == 4)
      sprintf(BufRtf,"4, integration was not completed because more than MAXNFE derivative\n"
         "evaluations were needed.  This is approximately (MAXNFE/6) steps.\n"
         "The user may continue by simply calling again.  The function counter\n"
         "will be reset to 0, and another MAXNFE function evaluations are allowed.\n");
else if (idid == 5)
      sprintf(BufRtf,"5, integration was not completed because the solution vanished,\n"
         "making a pure relative error test impossible.  The user must use\n"
         "a non-zero ABSERR to continue.  Using the one-step integration mode\n"
         "for one step is a good way to proceed.\n");
else if (idid == 6)
      sprintf(BufRtf,"6, integration was not completed because the requested accuracy\n"
         "could not be achieved, even using the smallest allowable stepsize.\n"
         "The user must increase the error tolerances ABSERR or RELERR before\n"
         "continuing.  It is also necessary to reset FLAG to 2 (or -2 when\n"
         "the one-step integration mode is being used).  The occurrence of\n"
         "FLAG = 6 indicates a trouble spot.  The solution is changing\n"
         "rapidly, or a singularity may be present.  It often is inadvisable\n"
         "to continue.\n");
else if (idid == 7)
      sprintf(BufRtf,"7, it is likely that this routine is inefficient for solving\n"
         "this problem.  Too much output is restricting the natural stepsize\n"
         "choice.  The user should use the one-step integration mode with\n"
         "the stepsize determined by the code.  If the user insists upon\n"
         "continuing the integration, reset FLAG to 2 before calling\n"
         "again.  Otherwise, execution will be terminated.\n");
else if (idid == 8)
      sprintf(BufRtf,"8, invalid input parameters, indicates one of the following:\n"
         "NEQN <= 0;\n"
         "T = TOUT and |FLAG| /= 1;\n"
         "RELERR < 0 or ABSERR < 0;\n"
         "FLAG == 0  or FLAG < -2 or 8 < FLAG.\n");
}

emit appendText(fsoItalic,clNavy );
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

FILE* CreateEtaTxtFile(const char *filename, const QString &LinitialDir, const QString &LfileName, QDateTime  &s_time, int option)
{

QString LocalName = LinitialDir + LfileName + "_" + filename;
QString mdate=s_time.toString(Qt::ISODate);
mdate += "    Etacha4 (GUI)      (C) INSP-ASUR JPR  + MSU  03/2017     (file \"";
mdate += LfileName;
mdate += "\")\n";

FILE *f = mfopen(LocalName,"wt");
if(!f) return nullptr;

if(option < 9)  // < en_fExcel
      {
      fprintf(f,"%s", mdate.toStdString().c_str());

      fputs(Star89,f);
      fprintf(f,"PROJECTILE: atomic number=%4.0f  incident charge=%4.0f  atomic mass=%4.0f\n"
                "    TARGET: atomic number=%4.0f      atomic mass=%4.0f      density=%6.3f g/cm3\n"
                "          incident energy=%8.3f MeV/u  velocity=%8.3f (au)\n",
             gZb,gQb,gAb, gZt,gAt,gDensity,gEnergy,ProjectileVelocity);
      fprintf(f,"           relative error=%.4g absolute error=%.4g\n",gUncertainRel,gUncertainAbs);

       if (UseEloss == 0) fprintf(f," No stopping power correction\n");
       else               fprintf(f,"  Energy loss: S1=%8.3f MeV/mg/cm2 at E1=%8.3f MeV/u \n"
                                    "          and: S2=%8.3f MeV/mg/cm2 at E2=%8.3f MeV/u \n",
                                        gO_dEdX1,gO_EnergyL1,gO_dEdX2,gO_EnergyL2);

      PrintCrossSections(f);
      fputs(Star118,f);
      }

int iZp = gZb;

if(option <=5)
        {
        fprintf(f,"Tar.thick.    ");
        for(int i=0; i<10; i++) fprintf(f,"%2de-      ",i+option*10);
        fprintf(f,"E_out\n");

        fprintf(f," (ug/cm2)    ");
        for(int i=0; i<10; i++) fprintf(f,"(%2de+)    ",iZp-i-option*10);
        fprintf(f," MeV/u\n");
        }
else if(option==6)
        {
        fprintf(f,"\n"
                "T (ug/cm2)  bare    1s     2s     2p    1s²    1s2s   1s2p  1s²2s  1s²2p   tot\n");
        }
else if(option==7)
        {
        fprintf(f,"\n"
                "T (ug/cm2)  y1s       y2s      y2p      ym      yn      Qm      Qm in   Qm out    PTOT\n");
        }
//----------------------------------
if(option<=7) fputs(Star118,f);

        //-------------
if(option==9)   //en_fExcel
        {
        fprintf(f,"Thick(ug)\tEnergy");
        int L = qMin(60.,gZb);
        for(int i=0; i<L; i++)    fprintf(f,"\t q%d   ",int(gZb)-i);
        fprintf(f,"\n");
        }
        //-------------



return f;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,tetan;
//      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn;
//      double o_Bk, o_Bl1, o_Bl2, o_Bm1, o_Bm2, o_Bn;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*void ETACHA::keyPressEvent(QKeyEvent *event)
{
if(event->key() == Qt::Key_Escape) cancelStatus = true;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
bool ETACHA::test_ESC()
{
qApp->processEvents();
return   cancelStatus;
}*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

int PrintBindingEnergies(FILE *fil)
{
if(!fil) return -1;

      fprintf(fil,"\n  Binding Energy");
      fprintf(fil,"\n  correction factors = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f",tetak,tetal1,tetal2,tetam1,tetam2,tetan);
      fprintf(fil,"\n   Effective charges = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f",zk,zl1,zl2,zm1,zm2,zn);
      fprintf(fil,"\n    Binding energies = %8.1f  %8.1f  %8.1f  %8.1f  %8.1f  %8.1f",o_Bk,o_Bl1,o_Bl2,o_Bm1,o_Bm2,o_Bn);
      fprintf(fil,"\n");

fputs(Star89,fil);
return 1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int PrintCrossSections(FILE *fil)
{
if(!fil) return -1;

fprintf(fil,
        "\nReference ('hydrogenic') cross sections in 10E-20 cm2:\n"
        "  CAPTURE (MEC+REC) TO:   1s =%10.3e  2s =%10.3e  2p =%10.3e\n"
        "                          3s =%10.3e  3p =%10.3e  3d =%10.3e   (n=4) =%10.3e\n"
        "  IONIZATION OF:          1s =%10.3e  2s =%10.3e  2p =%10.3e\n"
        "                          3s =%10.3e  3p =%10.3e  3d =%10.3e   (n=4) =%10.3e\n"
        "1s EXCITATION TO:         2s =%10.3e  2p =%10.3e\n"
        "                          3s =%10.3e  3p =%10.3e  3d =%10.3e   (n=4) =%10.3e\n"
        "2s EXCITATION TO:         3s =%10.3e  3p =%10.3e  3d =%10.3e   (n=4) =%10.3e\n"
        "2p EXCITATION TO:         3s =%10.3e  3p =%10.3e  3d =%10.3e   (n=4) =%10.3e\n"
        "EXCITATION To n=4 from:   3s =%10.3e  3p =%10.3e  3d =%10.3e\n"
        "INTRASHELL EXCITATION:    2s to 2p =%10.3e      3s to 3p =%10.3e\n",
                GSEC[ 1],GSEC[ 2],GSEC[ 3],GSEC[ 4],GSEC[ 5],GSEC[ 6],GSEC[ 7],GSEC[ 8],GSEC[ 9],GSEC[10],
                GSEC[11],GSEC[12],GSEC[13],GSEC[14],GSEC[15],GSEC[16],GSEC[17],GSEC[18],GSEC[19],GSEC[20],
                GSEC[21],GSEC[22],GSEC[23],GSEC[24],GSEC[25],GSEC[26],GSEC[27],GSEC[28],GSEC[29],GSEC[30],
                GSEC[31],GSEC[32],GSEC[33]);

return 1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int PrintCurrentDistribution(FILE *fil, int option, double  T, double EE)
{
if(!fil) return -1;

sprintf(BufRtf,"\n%9.2f ",T);
for(int I=1+option*10; I<=10+option*10; I++)
                        sprintf(&BufRtf[strlen(BufRtf)],"%9.5f ",PR[I]);

sprintf(&BufRtf[strlen(BufRtf)]," %9.3f ",EE);

fputs(BufRtf,fil);
return 1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int PrintCurrentExcelDistribution(FILE *fil, double  T, double EE)
{
if(!fil) return -1;

int L = qMin(60.,gZb);

fprintf(fil,"%.3e\t%.5e",T,EE);
        for(int i=1; i<=L; i++)    fprintf(fil,"\t%.3e",PR[i]);

fprintf(fil,"\n");

return 1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int PrintFinalDistribution  (FILE *f, int option, QDateTime &Time3, qint64 dif, double EE)
{
if(!f) {return -1;}

       fprintf(f,"\n");
       fputs(Star118,f);

if(option <=5)
        {
        fprintf(f,"Tar.thick.    ");
        for(int i=0; i<10; i++) fprintf(f,"%2de-      ",i+option*10);
        fprintf(f,"E_out\n");

        fprintf(f," (ug/cm2)    ");
        for(int i=0; i<10; i++) fprintf(f,"(%2de+)    ",int(gZb)-i-option*10);
        fprintf(f," MeV/u\n");
        }
else if(option==6)
      {
      fprintf(f,
                "T (ug/cm2)  bare    1s     2s     2p    1s²    1s2s   1s2p  1s²2s  1s²2p   tot\n");
      fputs(Star118,f);
      for(int n3=327; n3<=1284; n3++)
              if (Y[n3] >= 0.01)
                fprintf(f,"Code=%3d    pop=%.3e\n",ICO3[n3],Y[n3]);
                ////c    ICO3=100*n4+n123
      }
else if(option==7)
      {
      fprintf(f,
                "T (ug/cm2)  y1s       y2s      y2p      ym      yn      Qm      Qm in   Qm out    PTOT\n");
       }
fputs(Star118,f);

QString elapsed = ElapsedTime(dif);
QString time1text=Time3.time().toString("hh:mm:ss");
fprintf(f,"\n        End time is  %s",time1text.toStdString().c_str());
fprintf(f,"\n    %s",elapsed.toStdString().c_str());
fprintf(f,"\n             Final energy : %.3f (MeV/u)", EE);

return 1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void CalcOlegSum(double *V)
{
static int counter=0;

      for(int i=0; i<8; i++) sumM[i]=0;

      for(int I=1;   I<=63;   I++) sumM[0] += V[I];
      for(int I=64;  I<=66;   I++) sumM[1] += V[I];
      for(int I=67;  I<=73;   I++) sumM[2] += V[I];
      for(int I=74;  I<=84;   I++) sumM[3] += V[I];
      for(int I=85;  I<=293;  I++) sumM[4] += V[I];
      for(int I=294; I<=326;  I++) sumM[5] += V[I];
      for(int I=327; I<=1283; I++) sumM[6] += V[I];

      for(int i=0; i<7; i++)
                        if(i!=4)
                                sumM[7] += (sumM[i]-1.);

counter++;

}
//--------------------OPTIMIZED CODE-----------------------------------------------------------------------------
/*
void CalcOlegSum(double *V) {
    static int counter = 0;

    constexpr int ranges[7][2] = {
        {1, 63},    // sumM[0]
        {64, 66},   // sumM[1]
        {67, 73},   // sumM[2]
        {74, 84},   // sumM[3]
        {85, 293},  // sumM[4]
        {294, 326}, // sumM[5]
        {327, 1283} // sumM[6]
    };

    double localSum[7] = {0};  // Temporary storage to avoid modifying sumM in each iteration

    for (int i = 1; i <= 1283; ++i) {
        if (i <= ranges[0][1])       localSum[0] += V[i];
        else if (i <= ranges[1][1])  localSum[1] += V[i];
        else if (i <= ranges[2][1])  localSum[2] += V[i];
        else if (i <= ranges[3][1])  localSum[3] += V[i];
        else if (i <= ranges[4][1])  localSum[4] += V[i];
        else if (i <= ranges[5][1])  localSum[5] += V[i];
        else                         localSum[6] += V[i];
    }

    // Store results in sumM
    for (int i = 0; i < 7; i++) sumM[i] = localSum[i];

    // Compute sumM[7] efficiently
    sumM[7] = -1.0 * localSum[4]; // Subtract sumM[4] from the total
    for (int i = 0; i < 7; i++) if (i != 4) sumM[7] += localSum[i];

    counter++;
}
*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void NormOlegSum(double *)
{
return;   //due to  bugs found by Toshi Sumikama 09/28/2021

/*
if(sumM[0]>1e-5)      for(int I=1;   I<=63;   I++)  Y[I]/=sumM[0];
//if(sumM[1]>1e-5)      for(int I=64;  I<=66;   I++)  Y[I]/=sumM[1];
//if(sumM[2]>1e-5)      for(int I=67;  I<=73;   I++)  Y[I]/=sumM[2];
//if(sumM[3]>1e-5)      for(int I=74;  I<=84;   I++)  Y[I]/=sumM[3];
//if(sumM[4]>1e-5)      for(int I=85;  I<=293;  I++)  Y[I]/=sumM[4];
//if(sumM[5]>1e-5)      for(int I=294; I<=326;  I++)  Y[I]/=sumM[5];
if(sumM[6]>1e-5)      for(int I=327; I<=1283; I++)  Y[I]/=sumM[6];
CalcOlegSum(V);*/
}



