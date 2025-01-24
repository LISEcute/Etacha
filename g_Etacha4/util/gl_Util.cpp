#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <QString>
#include "L_Loss/dll_eloss.h"
#include "g_Etacha4/win/e_Constant.h"

using namespace std;


extern double  EnergyResidueFromRange(double Zp, double Ap,  double Zt, double At, double E, double thick);
QString numberKillZero(double v, char ch,  int k, int m=0); // only for 'f'

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
const char *symb=
"n H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZn"
"GaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeCsBaLaCePrNdPm"
"SmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaU "
"NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtB0B1B2B3B4B5B6B7B8B9C0C1C2C3C"
"4C5C6C7C8C9D0??";


double pow_int(double par, int power);
long double pow_int(long double par, int power);
int pow_m1(int power);

double pow1(double par);
double pow2(double par);
double pow3(double par);
double pow4(double par);
double pow5(double par);
double pow6(double par);
double pow8(double par);
double pow10(double par);

complex<double> pow2(complex<double> par);
double pow2I(double par);

char *GetNextSymbol(char *s);
char *GetNextDelimeter(char *s);
char *GetIntFromString(int &V, char *s);
char *GetDoubleFromString(double &V, char *s);
double mzsqrt(double X);
char* eos(char *b);
const char* ElementName(int IZ);
int GetZfromName(char *el_init);

//double GetPrivateProfileDouble(char* lpAppName, char* lpKeyName,  double def_value, char* ini_file);
//unsigned long GetPrivateProfileULong(char* lpAppName, char* lpKeyName,  unsigned long def_value, char* ini_file);

double EnergyResidue(double Zp, double Ap,  double Zt, double At, double E, double thick );
double E_to_Beta(double E);
double E_to_Gamma(double E);
double Beta_to_Gamma(double beta);
double Beta_to_E(double Beta);
double Gamma_to_E(double gamma);
double Velocity_au(double E);
double ANINT_O (double x);
//void DateTime(Struct_Datetime &s);
QString ElapsedTime(qint64 msec);
double fsign(double a, double b);
//bool test_ESC(void);
double StoppingPower( int Zp, int Ap, int Zt, int At, double E);
double EnergyResidueRng(double Zp, double Ap,  double Zt, double At, double E, double thick );
double EnergyResidueSP(double Zp, double Ap,  double Zt, double At, double E, double thick );

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
long double pow_int(long double par, int power)
{
long double S=1;

for(int i=1; i<=power; i++) S*=par;
return S;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow_int(double par, int power)
{
double S=1;

for(int i=1; i<=power; i++) S*=par;
return S;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int pow_m1(int power)
{
power = abs(power);
int k = power%2;

return k==1 ? -1 : 1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
complex<double> pow2(complex<double> par)
{
return par*par;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow1(double par)
{
return par;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow2(double par)
{
return par*par;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow2I(double par)
{
if(par==0) return -777;
return 1./pow2(par);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow3(double par)
{
return par*par*par;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow4(double par)
{
double par2 = par*par;
return par2*par2;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow5(double par)
{
return pow4(par)*par;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow6(double par)
{
double p2 = pow2(par);
double p4 = p2*p2;
return p2*p4;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow8(double par)
{
double p2 = pow2(par);
double p4 = p2*p2;
return p4*p4;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow10(double par)
{
double p2 = pow2(par);
double p4 = p2*p2;
return p4*p4*p2;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

//---------------------------------------------------------------------------
char *GetNextSymbol(char *s)
{
char c;

while(true) {
        c=*s;
        if(c=='\0' || c=='\n' || c=='\r') {return NULL;}

        if( !(c==' ' || c=='\t' || c==',')) break;
        s++;
        }
return s;
}
//---------------------------------------------------------------------------
char *GetNextDelimeter(char *s)
{
char c;

while(true) {
        c=*s;
        if(c=='\0' || c=='\n' || c=='\r')   {return NULL;}
        if( c==' ' || c=='\t' || c==',') break;
        s++;
        }
return s;
}

//---------------------------------------------------------------------------
char *GetIntFromString(int &V, char *s)
{
s=GetNextSymbol(s); if(!s) return NULL;
V=atoi(s);
return GetNextDelimeter(s);
}
//---------------------------------------------------------------------------
char *GetDoubleFromString(double &V, char *s)
{
s=GetNextSymbol(s); if(!s) return NULL;
V=atof(s);
return GetNextDelimeter(s);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double mzsqrt(double X) {
      if(X <0) X=0;
      return sqrt(X);
      };
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char* eos(char *b)
{ char *s; s=b;

while(1)  {
   if(*s==0 || *s=='\n' || *s==';' ) break;
   s++;
   }
return s;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
const char* ElementName(int IZ) {


if(IZ<0 || IZ>130) IZ=131;
return &symb[IZ*2];
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int GetZfromName(char *el_init)
 {
int i;

char *s, *seos;                                     // format: ii cc  ii+\0
s=el_init;
seos=eos(s);

while(*s <= '9') { if(s==seos) return 131;      s++;   }

if( *s     <= 'z' && *s       >= 'a' ) 	  *s  =  char(*s     - ('a'-'A'));
if( *(s+1) <= 'Z' && *(s+1)   >= 'A' ) *(s+1) =  char(*(s+1) + ('a'-'A'));


//int shift=2;

 for(i=0; i<=130; i++)							// from 2 symbols
 		if(strncmp(s, symb+2*i, 2) == 0) break;

 if(i > 130 )  {
  //	shift=1;
 	for(i=0; i<=130; i++)                                 // from 1 symbol
      	if(  *(symb+2*i  ) == *s   &&
                 *(symb+2*i+1) == ' ' 	) break;
      }

return i;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double EnergyResidueRng(double Zp, double Ap,  double Zt, double At, double E, double thick )
{
double Er =  EnergyResidueFromRange(Zp, Ap,  Zt, At, E, thick);
return Er;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double EnergyResidueSP(double Zp, double Ap,  double Zt, double At, double E, double thick )
{
if(thick <=0) return E;

int N_Steps=1000;
if(thick<10) N_Steps/=10;
if(thick< 1) N_Steps/=10;

double t = thick/double(N_Steps);

for(int i=0; i<N_Steps; i++)
        {
        double SE = StoppingPower(Zp, Ap, Zt, At, E);
        double Eloss = SE * t / Ap;
        E -= Eloss;
        if(E<=0) return 0;
        }
return E;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double StoppingPower( int Zp, int Ap, int Zt, int At, double E)
{
if(At<1)  At=1;
double Lz = ZieglerStopping(Zp,E,Zt,0)*(0.6022/At);
double Ln = NuclearStopping(Zp,Ap,Zt,At,E);
double d = Lz+Ln;

return d;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double E_to_Beta(double E)
{
double gamma=E_to_Gamma(E);
if(gamma <= 1. ) return 0;
else            return sqrt(1.-1./gamma/gamma);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double E_to_Gamma(double E)
{
  if(E<=0) return 1;
  else    return E / amu_MeV + 1.;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW==
double Beta_to_E(double beta)
{
if(beta<=0 || beta>=1) return 0;
double gamma=Beta_to_Gamma(beta);
return Gamma_to_E(gamma);
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW==
double Beta_to_Gamma(double beta)   {
if(beta<=0)return 1;
double beta2 = beta*beta;
double k = 1.- beta2;
#define minv 1e-40
if(fabs(k)< minv)  k=minv;
return sqrt(1./k);
 };
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW==
double Gamma_to_E(double gamma)
{ return (gamma-1.) * amu_MeV; };
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double Velocity_au(double E)
{
double beta =  E_to_Beta(E);
double vau = beta / FineStructConst;
return vau;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double ANINT_O (double x)  // Oleg
{
int k  =  x + 0.5;     //ANINT(A [, KIND]) rounds its argument to the nearest whole number.
return double(k);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*void DateTime(Struct_Datetime &s)
{
struct date s_date;
struct time s_time;
getdate(&s_date);
gettime(&s_time);

s.iyr  = s_date.da_year;
s.imon = s_date.da_mon;
s.iday = s_date.da_day;
s.ihour = s_time.ti_hour;
s.imin  = s_time.ti_min;
s.isec  = s_time.ti_sec;
s.ihund = s_time.ti_hund;

double mons = s.iyr * 12  + s.imon;
double days =  mons * 365./12. + s.iday;
double hous =  days * 24.      + s.ihour;
double mins =  hous * 60.      + s.imin;
     s.secs =  mins * 60.      + s.isec  + s.ihund/100.;
}*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QString ElapsedTime(qint64 msec)
{
QString buff;
    double dsec = msec/1000.;
    int sec    =  dsec;  // ms
    int hours  =  sec/3600;
    int mins    = (sec - hours*3600 )/60;
    sec %= 60;

QString result = buff.asprintf("Elapsed time is  %02d:%02d:%02d (or %.3f sec)",
                                hours, mins, sec, dsec);
return result;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double fsign(double a, double b)
{
//The kind of the return value is that of A and B. If B >= 0 then the result is ABS(A), else it is -ABS(A).

if(b>=0) return  fabs(a);
else     return -fabs(a);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
