#include "e_Etacha4.h"
#include "../win/e_myextern.h"
#include "../win/e_Constant.h"
#include <math.h>
//#include "MyMessages.h"
//#include <stdio.h>
//#define max(a,b) (a > b ? a : b)




extern double pow2(double par);
extern double pow2I(double par);

extern double Velocity_au(double E);
extern double Beta_to_E(double Beta);

extern void PION(double Ecurr);
extern void SEIK(double Ecurr, double &ck, double &cl, double &cm, double &cn, double &ct, double &cs);
extern void REC(double Ecurr,double &srk,double &srls,double &srlp);

extern void sexi(double Ecurr, double &e2s,double &e2p,double &e3s,double &e3p,double &e3d,double &es4,double &es5);
extern void sex2(double Ecurr, double &s3s,double &s3p,double &s3d,double &p3s,double &p3p,double &p3d,double &e2s4,double &e2p4,double &e2s5,double &e2p5);
extern void sex3(double Ecurr, double &e3s4,double &e3p4,double &e3d4,double &e3s5,double &e3p5,double &e3d5,double &e4s5,double &e4p5,double &e45,double &e456);
extern void snl(double Ecurr, double &esp,double &esp3,double &epd3);

/*extern void Tceis(double &tot, int netat,  TForm *tobject);
extern void


*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int ETACHA::DONAUT()
//void donaut(QP,Z1,AP,EP,Z2,AC,RHO,EPM,iflag,E1,S1,E2,S2,istp,iprt,ilgn)

{
int iflag;  //  subroutine string parameters

//_______________________________________________________________________

//     table of cross sections (units 10-20cm2):
//     indices                 secs              indices                 GSEC

//           1                 MeC1s SEIK              1                 C1s
//           2                 MeC2s SEIK              2                 C2s
//           3                 MeC2p SEIK              3                 C2p
//           4                 MeC3s SEIK              4                 C3s
//           5                 MeC3p SEIK              5                 C3p
//           6                 MeC3d SEIK              6                 C3d
//           7                 MeC4  SEIK              7                 C4

//           8                 ReC1s                   8                 Ion1s          15+22  (Ion1s CDW EIS + Se1s5 SE)
//           9                 ReC2s                   9                 Ion2s          16+23
//           10                ReC2p                   10                Ion2p          17+24
//           11                ReC3s                   11                Ion3s
//           12                ReC3p                   12                Ion3p
//           13                ReC3d                   13                Ion3d
//           14                ReC4                    14                Ion4           21+28


//           15                Ion1s CDW   EIS         15                Se1s2s
//           16                Ion2s CDW   EIS         16                Se1s2p
//           17                Ion2p CDW   EIS         17                Se1s3s
//           18                Ion3s Scaling           18                Se1s3p
//           19                Ion3p Scaling           19                Se1s3d
//           20                Ion3d Scaling           20                Se1s4
//           21                Ion4  Scaling

//     Scaling means Ion(3l,Z) = Ion(1s,Z/3) and ion(4,Z) = Ion(2p,Z/2)
//
//           22                Se1s5 SE                21                Se2s3s
//           23                Se2s5 SE                22                Se2s3p
//           24                Se2p5 SE                23                Se2s3d
//           25                Se3s5 SE                24                Se2s4
//           26                Se3p5 SE
//           27                Se3d5 Scaling SE        25                Se2p3s
//           28                Se456 Scaling SE        26                Se2p3p
//                                                     27                Se2p3d
//           29                Se1s2s      SE          28                Se2p4
//           30                Se1s2p      SE
//           31                Se1s3s      SE            29              Se3s4
//           32                Se1s3p      SE            30              Se3p4
//           33                Se1s3d      SE            31              Se3d4
//StSE[1]    34                Se1s4       SE
//
//SeSE[15]   35                Se2s3s      SE            32             Se2s2p
//SeSE[16]   36                Se2s3p      SE            33             Se3s3p
//SeSE[17]   37                Se2s3d      SE            34             Se3p3d
//StSE[2 ]   38                Se2s4       SE

//SeSE[27]   39                Se2p3s      SE
//SeSE[28]   40                Se2p3p      SE
//SeSE[29]   41                Se2p3d      SE

//StSE[3 ]   42                Se2p4 SE
//StSE[4 ]   43                Se3s4 SE
//StSE[5 ]   44                Se3p4 SE
//           45                Se3d4 Scaling SE                               ==> R3d4*(Gsecs[43]+Gsecs[44]);

//esp        46                Se2s2p
//esp3       47                Se3s3p
//epd3       48                Se3p3d


//__________secs_____________________________________________________Gsecs________________
//...... ibin=0 empirical saturation correction for born .....
//
//...... ibin=1 binding correction included in born (not recommended) .....
//...... ibin=2 no empirical correction and no binding correction
//...... correction are only for PWBA (BORN1) cross sections
//...... as calculated in PION, SEXI, SEX2 and SEX3,
//...... are only used for the evolution of cross sections with effective charge,
//...... and should not make much a change

//...... Excitation to n=5 (secs(i), i=22,28) is not added to ionisation and set
//...... to 0

//...... CDW is not used for capture cross sections
//________________________________________________________________________________
int      iter=0;
//double   y4=0.;

      ibin = ibinParameter;  // Oleg

      iflag=1;
      y1s=0.;
      yKm=0.;

      y2s=0.;
      y2p=0.;
      yL=0.;
      yL1m=0.;
      yL2m=0.;

      y3s=0.;
      y3p=0.;
      y3d=0.;
      yM=0.;
      yMp=0.;
      yMm=0.;
      yM1m=0.;
      yM2m=0.;

      tetak=1.;
      tetal1=1.;
      tetal2=1.;
      tetam1=1.;
      tetam2=1.;
      tetan=1.;

      // make here initialization of local etacha variables


//      Qp = gQb;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


//---------------------------------------------------
//     CROSS SECTIONS
//---------------------------------------------------
// Hydrogenic (reference) cross sections calculation ....'
strcpy(BufRtf,"Hydrogenic (reference) cross sections calculation ....");
emit appendText(fsoBold,clBlack);


///---------------------------------------------- if(iter==0) starts
if (iter == 0)
        {
        for(int I=0; I<=48; I++) Gcor[I]=1.;

        CSEC(Ep,0);

        for(int I=1; I<=48; I++) Gseci[I]=Gsecs[I];   /// Oleg : Initial

                //     Preliminary calculation of scaling factors for excitation cross sections
                //     of d and f states in SE approximation @ v=zp

      double Esave = Ep;
      double beta = Zp / 137.036;
      Ep = Beta_to_E(beta);

      strcpy(BufRtf,"excitation 3l&rarr;4l"); emit appendText(fsoItalic,clBlack );
      sprintf(BufRtf,"E=%.3f in front of sex3",Ep); emit appendText(fsoBold,clBlack );

      double e3s4=0,e3p4=0,e3d4=0,e3s5=0,e3p5=0,e3d5=0,e4s5=0,e4p5=0,e45=0,e456=0;
      sex3(Ep,e3s4,e3p4,e3d4,e3s5,e3p5,e3d5,e4s5,e4p5,e45,e456);

      double  R3d4 = e3d4 / qMax(1e-50, (e3s4+e3p4));
      double  R3d5 = e3d5 / qMax(1e-50, (e3s5+e3p5));

                //     e456 is not corrected for 1/n**3 law anymore...
      double R45=e456/(e4s5+e4p5);
      Ep = Esave;


      //------------------------------------- ionization model ==0 start
  if(IonizationModel==0)
      {

      strcpy(BufRtf,"<br>CDWEIS Ionisation Cross sections"); emit appendText(fsoUnderline, clNavy );
      strcpy(BufRtf," Total cross sections for single ionization of Hydrogenic or Helium<br>"
                    " targets by bare ion impact using the CDW-EIS (RHF) approximation<br>"); emit appendText(0, clGreen );

      //ninis=10;     // 1s - initial state
      strcpy(BufRtf,"Hydrogenic Target - 1s Initial State"); emit appendText(fsoUnderline, clBlue);emit appendShell(fsoUnderline, clBlue);
      double tot;
      Tceis(Ep,tot,10);
      Gsecs[15]=tot*1.e20;

      //ninis=20;     // 2s - initial state
      strcpy(BufRtf,"Hydrogenic Target - 2s Initial State"); emit appendText(fsoUnderline, clBlue);emit appendShell(fsoUnderline, clBlue);
      Tceis(Ep,tot,20);
      Gsecs[16]=tot*1.e20;

      //ninis=21;    // 2p - initial state
      strcpy(BufRtf,"Hydrogenic Target - 2p Initial State"); emit appendText(fsoUnderline, clBlue);emit appendShell(fsoUnderline, clBlue);
      Tceis(Ep,tot,21);
      Gsecs[17]=tot*1.e20;

      strcpy(BufRtf,"Hydrogenic Target - 1s Initial State; Z/3"); emit appendText(fsoUnderline, clBlue);emit appendShell(fsoUnderline, clBlue);
      double zp3=Zp;  // 3s,3p,3d - initial state SCALE
      Zp=zp3/3.;
      //ninis=10;
      Tceis(Ep,tot,10);
      Gsecs[18]=Gsecs[19]=Gsecs[20]=tot*1.e20;
      Zp=zp3;

      if(EtachaVersion >= etacha_v34)
                {
                strcpy(BufRtf,"Hydrogenic Target - 2p Initial State; Z/3"); emit appendText(fsoUnderline, clBlue);emit appendShell(fsoUnderline, clBlue);
                double zp2=Zp;
                Zp=zp2/2.;
                Tceis(Ep,tot,21);
                Gsecs[21]=tot*1.e20;
                Zp=zp2;
                }
      }
      //------------------------------------- ionization model ==0 end

      for(int i=15; i<=21; i++) Gcor[i]= Gseci[i] > 0 ? Gsecs[i]/Gseci[i] : 1;


      //------------------------------------- excitation model ==0 start
if(ExcitationModel==0)
      {
      strcpy(BufRtf,"<br>SE Excitation Cross sections"); emit appendText(fsoUnderline,clNavy);emit appendShell(fsoBold,clRed);
      strcpy(BufRtf," total cross sections of monoelectronic atoms excitation<BR>"
                    " for colision with nude projectiles in the SE (Symmetric-Eikonal) aproximation"); emit appendText(fsoItalic,0);

      double zp8 = Zp;
      double zt8 = Zt;
      double E8  = Ep;
      SEnlm(zp8,zt8,E8);                       //   version 23 gives the same result

                        //     Excitation to n>=5 set to 0 in the present version
        if(EtachaVersion == etacha_v45)
              {                          //     1s - 5
              Gsecs[22]=StSE[6]*1.e20;
                                        //     2s,2p - 5
              Gsecs[23]=StSE[7]*1.e20;
              Gsecs[24]=StSE[8]*1.e20;
                                        //     3s,3p,3d - 5
              Gsecs[25]=StSE[9]*1.e20;
              Gsecs[26]=StSE[10]*1.e20;

              Gsecs[27]=R3d5*(Gsecs[25]+Gsecs[26]);
                                //     R45 does not include correction for the sum on 1/n**3 anymore
              Gsecs[28]=R45*(StSE[11]+StSE[12])*1.e20;

              for(int i=22; i<=28; i++) Gcor[i]= Gseci[i] > 0 ? Gsecs[i]/Gseci[i] : 1;
              }
        else  {
                                //     Excitation to n=5 and 6 set to 0
              for(int i=22; i<=28; i++) Gcor[i]=1;
              }


      Gsecs[29]=SeSE[1 ]*1.e20;
      Gsecs[30]=SeSE[2 ]*1.e20;
      Gsecs[31]=SeSE[3 ]*1.e20;
      Gsecs[32]=SeSE[4 ]*1.e20;
      Gsecs[33]=SeSE[5 ]*1.e20;

      Gsecs[34]=StSE[1 ]*1.e20;
      Gsecs[38]=StSE[2 ]*1.e20;

      if(EtachaVersion >= etacha_v34)
              {
              Gsecs[35]=SeSE[15]*1.e20;    // other indexation SeSE for etacha23, even function is the same!       ... somewhere is wrong!!!
              Gsecs[36]=SeSE[16]*1.e20;
              Gsecs[37]=SeSE[17]*1.e20;
              Gsecs[39]=SeSE[27]*1.e20;
              Gsecs[40]=SeSE[28]*1.e20;
              Gsecs[41]=SeSE[29]*1.e20;
              }
      else    {
              Gsecs[35]=SeSE[10]*1.e20;    // other indexation  for etacha23         10
              Gsecs[36]=SeSE[11]*1.e20;
              Gsecs[37]=SeSE[12]*1.e20;
              Gsecs[39]=SeSE[17]*1.e20;
              Gsecs[40]=SeSE[18]*1.e20;
              Gsecs[41]=SeSE[19]*1.e20;
              }

      Gsecs[42]=StSE[3 ]*1.e20;
      Gsecs[43]=StSE[4 ]*1.e20;
      Gsecs[44]=StSE[5 ]*1.e20;
      Gsecs[45]=R3d4*(Gsecs[43]+Gsecs[44]);
      }
      //------------------------------------- excitation model ==0 stop


      for(int i=29; i<=45; i++) Gcor[i]= Gseci[i] > 0 ? Gsecs[i]/Gseci[i] : 1;

                //     note that Gsecs[22) -> Gsecs[28) are possibly set to 0...

      GSEC[8] =Gsecs[15]+Gsecs[22];
      GSEC[9] =Gsecs[16]+Gsecs[23];
      GSEC[10]=Gsecs[17]+Gsecs[24];
      GSEC[11]=Gsecs[18]+Gsecs[25];
      GSEC[12]=Gsecs[19]+Gsecs[26];
      GSEC[13]=Gsecs[20]+Gsecs[27];
      GSEC[14]=Gsecs[21]+Gsecs[28];

      for(int i=15; i<=31; i++)  GSEC[i]=Gsecs[i+14];

      iter=1;
      }
///---------------------------------------------- if(iter==0) stop

      for(int I=1; I<=48; I++) Gsecs[I]= Gcor[I]*Gseci[I];

        for(int i=1; i<=7; i++)
                {
                GSEC[i  ]=Gsecs[i   ]+Gsecs[i+7 ]; //     add MEC and REC
                GSEC[i+7]=Gsecs[i+14]+Gsecs[i+21]; //     add ionisation and nl -> 5 excitation
                }

                                //     shift index for other processes (excitation)
      for(int i=15; i<=34; i++)  GSEC[i]=Gsecs[i+14];


//************************************************************************


if(EtachaVersion <= etacha_v3)
                {
                GSEC[ 8] += GSEC[20];  GSEC[20]=0;                       // Se1s4
                GSEC[ 9] += GSEC[24];  GSEC[24]=0;                       // Se2s4
                GSEC[10] += GSEC[28];  GSEC[28]=0;                       // Se2p4
                GSEC[11] += GSEC[29];  GSEC[29]=0;                       // Se3s4
                GSEC[12] += GSEC[30];  GSEC[30]=0;                       // Se3p4
                GSEC[13] += GSEC[31];  GSEC[31]=0;                       // Se3d4
                }

      // ...................... fin des modifs ..........................

        ///Oleg,  write data to file!
return  iflag;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void ETACHA::CSEC(double Ecurr, int iStep)
{

//      double GSEC[34),Gcor[48),Gsecs[48);
//      common/seceff/GSEC,Gcor,Gsecs;
//      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
//      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
//      common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m;
//      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,tetan;
//      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn;
//      common/don/zp,E,zt;
//      common/corr/ibin;
//   iStep  -  0 from donaut, 1 - from


for(int i=0; i<49; i++) Gsecs[i]=0;   // Oleg

//OLEG double QM = Zp-(y1s+yl+ym+yn);
double E=Ecurr;
double VC = sqrt(1.-pow2I(1+E/931.5));         // beta
double vu=VC*137.036;
//double vu = Velocity_au(Ecurr);
bool Flag34 = EtachaVersion >= etacha_v34;


//-------------------------------------------------------   ionization
if(iStep==0) { strcpy(BufRtf,"Ionization (PWBA)"); emit appendText(fsoItalic,clBlack );}
PION(Ecurr);                  // SECTIONS EFFICACES D'IONISATION PWBA

double Dk  = o_TCk;
double Dl1 = o_TCl1;
double Dl2 = o_TCl2;
double Dm1 = o_TCm1;
double Dm2 = o_TCm2;
double Dn  = o_TCn;

    //  it is used later for Gsecs[15-21] with correctios
//---------------------------------------------------------------
if(iStep==0) {strcpy(BufRtf,"Capture (MEC)"); emit appendText(fsoItalic,clBlack );}

double ck=0,cl=0,cm=0,cn=0,ct=0,cs=0;   //Calcule les sections efficaces de capture non radiative
                                        //dans l'approximation eikonale
SEIK(E,ck,cl,cm,cn,ct,cs);

double sig14=ck+cl+cm+cn;
//Oleg double sig530;

//     ct is the sum of seik from 5 to 30
//Oleg if (sig14 <= cs)   sig530=cs-sig14;
//Oleg else               sig530=ct;

//     MEC
     Gsecs[1]=ck;
      Gsecs[2]=0.25*cl;
      Gsecs[3]=0.75*cl;
      Gsecs[4]=cm*1./9.;
      Gsecs[5]=cm*3./9.;
      Gsecs[6]=cm*5./9.;
      Gsecs[7]=cn;
// originally commented     if (Zp > 60.) Gsecs[7] += sig530;

//------------------------------------------------------------------------------------------
//************************
//     REC
double srk=0,srls=0,srlp=0;             // Calcule les sections efficaces de REC 1s,2s et 2p
                                        // dans l'approximation type Bethe et Salpeter

if(iStep==0) {strcpy(BufRtf,"Capture (REC)"); emit appendText(fsoItalic,clBlack );}
REC(E,srk,srls,srlp);


      double srm = 8./27.*(srls+srlp);
//Oleg       double stm = cm+srm;
      double srn = 1./8.*(srls+srlp);

      Gsecs[8] =srk;
      Gsecs[9] =srls;
      Gsecs[10]=srlp;
      Gsecs[11]=srls*8./27.;
      Gsecs[12]=srlp*8./27.;
      Gsecs[13]=srlp*8./45.;

if(Flag34)  Gsecs[14]=srn;
else        Gsecs[14] = 0;    

//------------------------------------------------------------------------------------------
//     EXCITATION

double e2s=0,e2p=0,e3s=0,e3p=0,e3d=0,es4=0, es5=0;
double s3s=0,s3p=0,s3d=0,p3s=0,p3p=0,p3d=0,e2s4=0,e2p4=0,e2s5=0,e2p5=0;
double e3s4=0,e3p4=0,e3d4=0,e3s5=0,e3p5=0,e3d5=0,e4s5=0,e4p5=0,e45=0,e456=0;
double esp=0,esp3=0,epd3=0;

 //Calcule les sections efficaces d'excitation 1s-2s, 1s-2p et   1s-n dans l'approximation PWBA
if(iStep==0) {strcpy(BufRtf,"excitation 1s&rarr;nl"); emit appendText(fsoItalic,clBlack );}
sexi(E,e2s,e2p,e3s,e3p,e3d,es4,es5);

if(iStep==0) {strcpy(BufRtf,"excitation 2l&rarr;nl"); emit appendText(fsoItalic,clBlack );}
sex2(E,s3s,s3p,s3d,p3s,p3p,p3d,e2s4,e2p4,e2s5,e2p5);            // d'excitation 2l-3l' et 2l-4l'

if(iStep==0) {strcpy(BufRtf,"excitation 3l&rarr;nl"); emit appendText(fsoItalic,clBlack );}
sex3(E,e3s4,e3p4,e3d4,e3s5,e3p5,e3d5,e4s5,e4p5,e45,e456);       // d'excitation 3l-4l'

                                //Calcule les sections efficaces d'excitation ns-np  dans    l'approximation PWBA
if(iStep==0) {strcpy(BufRtf,"excitation ns&rarr;np"); emit appendText(fsoItalic,clBlack );}
snl(E,esp,esp3,epd3);                                           // ns-np

//-------------------------------------  if start for "ibin"
if(ibin==1|| ibin==2)
        {
        //c.......... binding correction included (1) or no correction (2) ..........
      Gsecs[15]=Dk;
      Gsecs[16]=Dl1;
      Gsecs[17]=Dl2;
      Gsecs[18]=Dm1;
      Gsecs[19]=Dm1;
      Gsecs[20]=Dm2;
      Gsecs[21]=Dn;

         //c excitation to n=5 will be added to ionization
        if(EtachaVersion == etacha_v45)
              {
              Gsecs[22]=es5;
        //     (es5 includes a sum up to infinity)
             Gsecs[23]=e2s5*3.05;
             Gsecs[24]=e2p5*3.05;
        //     (theoretical factor in 1/n**3 = 3.049358 for n=5)
             Gsecs[25]=e3s5*2.5;
             Gsecs[26]=e3p5*2.5;
             Gsecs[27]=e3d5*2.5;
        //     (3 ->5 => 2.5 instead of 3.05)
              Gsecs[28]=e456*2.5;
        //     (e456=e45+e46 instead of
        //     e456=e45+2.5*e46 -> 2.5 au lieu de 3.5)
        }
        else {  Gsecs[22] = Gsecs[23] = Gsecs[24] = Gsecs[25] = Gsecs[26] = Gsecs[27] = Gsecs[28] = 0.;}

        //c excitation from 1s
        Gsecs[29]=e2s;
        Gsecs[30]=e2p;
        Gsecs[31]=e3s;
        Gsecs[32]=e3p;
        Gsecs[33]=e3d;
        Gsecs[34]=es4;
        //c excitation from n=2
        Gsecs[35]=s3s;
        Gsecs[36]=s3p;
        Gsecs[37]=s3d;
        Gsecs[38]=e2s4;
        Gsecs[39]=p3s;
        Gsecs[40]=p3p;
        Gsecs[41]=p3d;
        Gsecs[42]=e2p4;
        //c excitation from n=3 to n=4
        Gsecs[43]=e3s4;
        Gsecs[44]=e3p4;
        Gsecs[45]=e3d4;
        //...........................
        }
//-------------------------------------  else for "ibin"
 else   {    // ibin==0
        //c....... empirical saturation correction (09/09/94) ...........
        //c    not used anymore since a long time ...
      double cse = 1.06;
      double csi = 0.735;

      double loc1 = 10.96*pow2(vu);
      double zs1 = sqrt((loc1*exp(0.111*Zp    )/pow((Zp   ),1.946)));
      double zs2 = sqrt((loc1*exp(0.111*Zp/1.5)/pow((Zp/2.),1.946)));
      double zs3 = sqrt((loc1*exp(0.111*Zp/1.5)/pow((Zp/3.),1.946)));


      double loc2e = exp(-cse*Zt/pow(vu,2.1));
      double loc2i = exp(-csi*Zt/pow(vu,2.1));

      double loc31 = pow2(zs1/Zt)* pow((1.-exp(-pow((Zt/zs1),2.5))),0.8);
      double loc32 = pow2(zs2/Zt)* pow((1.-exp(-pow((Zt/zs2),2.5))),0.8);
      double loc33 = pow2(zs3/Zt)* pow((1.-exp(-pow((Zt/zs3),2.5))),0.8);

      double sate1 = loc2e * loc31;
      double sate2 = loc2e * loc32;
      double sate3 = loc2e * loc33;
      double sati1 = loc2i * loc31;
      double sati2 = loc2i * loc32;
      double sati3 = loc2i * loc33;
        //c......................................................
      Gsecs[15] = Dk *sati1;
      Gsecs[16] = Dl1*sati2;
      Gsecs[17] = Dl2*sati2;
      Gsecs[18] = Dm1*sati3;
      Gsecs[19] = Dm1*sati3;
      Gsecs[20] = Dm2*sati3;
      Gsecs[21] = Dn;

      if(EtachaVersion == etacha_v45)
              {         //c excitation to n=5 to be added to ionization
              Gsecs[22] = es5*sate1;
                        //     (es5 includes a sum up to infinity)
              Gsecs[23] = e2s5*3.05*sate2;
              Gsecs[24] = e2p5*3.05*sate2;
                //     (theoretical factor in 1/n**3 = 3.049358 for n=5)
              Gsecs[25] = e3s5*2.5*sate3;
              Gsecs[26] = e3p5*2.5*sate3;
              Gsecs[27] = e3d5*2.5*sate3;
                //     (3 ->5 => 2.5 in place of 3.05)
              Gsecs[28] = e456*2.5;
                //     (e456=e45+2.5*e46 -> 2.5 in place of 3.5)
              }
         else {
              Gsecs[22] = Gsecs[23] = Gsecs[24] = Gsecs[25] = Gsecs[26] = Gsecs[27] = Gsecs[28] = 0.;
              }

        //c excitation from 1s
      Gsecs[29] = e2s*sate1;
      Gsecs[30] = e2p*sate1;
      Gsecs[31] = e3s*sate1;
      Gsecs[32] = e3p*sate1;
      Gsecs[33] = e3d*sate1;
      Gsecs[34] = es4*sate1;
        //c excitation a partir de n=2
      Gsecs[35] = s3s*sate2;
      Gsecs[36] = s3p*sate2;
      Gsecs[37] = s3d*sate2;
      Gsecs[38] = e2s4*sate2;
      Gsecs[39] = p3s*sate2;
      Gsecs[40] = p3p*sate2;
      Gsecs[41] = p3d*sate2;
      Gsecs[42] = e2p4*sate2;
        //c excitation a partir de n=3 vers n=4
      Gsecs[43] = e3s4*sate3;
      Gsecs[44] = e3p4*sate3;
      Gsecs[45] = e3d4*sate3;
      }
//-------------------------------------  end if for "ibin"

                                        // intracouche
      Gsecs[46] = esp;
      Gsecs[47] = esp3;
      Gsecs[48] = epd3;
                                        //......................

      for(int I=1; I<=48; I++)   Gsecs[I] *= Gcor[I];

      for(int i=1; i<=7;  i++)   GSEC[i] = Gsecs[i]   + Gsecs[i+7];         //  capture totale (1 a 7)
      for(int i=8; i<=14; i++)   GSEC[i] = Gsecs[i+7] + Gsecs[i+14];        //  perte totale (8 a 14)

      for(int i=15; i<=34; i++)  GSEC[i] = Gsecs[i+14];                     //  excitation (15 a 34)


      if(EtachaVersion <= etacha_v3)
                {
                GSEC[ 8] += GSEC[20];  GSEC[20]=0;                       // Se1s4
                GSEC[ 9] += GSEC[24];  GSEC[24]=0;                       // Se2s4
                GSEC[10] += GSEC[28];  GSEC[28]=0;                       // Se2p4
                GSEC[11] += GSEC[29];  GSEC[29]=0;                       // Se3s4
                GSEC[12] += GSEC[30];  GSEC[30]=0;                       // Se3p4
                GSEC[13] += GSEC[31];  GSEC[31]=0;                       // Se3d4
                }

}


