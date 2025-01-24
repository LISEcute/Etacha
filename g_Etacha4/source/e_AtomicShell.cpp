//#pragma hdrstop
//#include <stdio.h>

#include "e_AtomicShell.h"

#define min(a,b) (a < b ? a : b)
#define max(a,b) (a > b ? a : b)

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
atom_shell::atom_shell()
{input_zq(1);}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
atom_shell::atom_shell(int iz, int iq)
{input_zq(iz,iq);}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int atom_shell::input_zq(int iz, int iq)
{
iz = min(110,iz);
iz = max(1  ,iz);
iq = max(0,  iq);
iq = min(iq, iz);

z  = iz;
ne = z-iq;

n1s=
n2s=n2p=
n3s=n3p=n3d=
n4s=n4p=n4d=n4f=
n5s=n5p=n5d=n5f=
n6s=n6p=n6d=
n7s=n7p=0;


//--------------------------------
        if (z <= 2)  n1s = z;
      else           n1s = 2;


      if (z <= 4)   {
                     n2s=z-2;
                     }
 else if (z <= 10)   {
                      n2s=2;
                      n2p=z-4;
                      }
 else if (z <= 12)  {
                      n2s=2;
                      n2p=6;
                      n3s=z-10;
                      }
 else if (z <= 18)  {
                      n2s=2;
                      n2p=6;
                      n3s=2;
                      n3p=z-12;
                      }
 else if (z <= 20)  {
                      n2s=2;
                      n2p=6;
                      n3s=2;
                      n3p=6;
                      n4s=z-18;
                      }
 else if (z <= 28)  {
                      n2s=2;
                      n2p=6;
                      n3s=2;
                      n3p=6;
                      n3d=z-20;
                      n4s=2;
                      }
 else if (z <= 36) {
                      n2s=2;
                      n2p=6;
                      n3s=2;
                      n3p=6;
                      n3d=10;
                      n4s=2;
                      n4p=z-30;
                      }
 else if (z <= 38) {
                      n2s=2;
                      n2p=6;
                      n3s=2;
                      n3p=6;
                      n4s=2;
                      n3d=10;
                      n4p=6;
                      n5s=z-36;
                      }
 else if (z <= 48) {
                      n2s=2;
                      n2p=6;
                      n3s=2;
                      n3p=6;
                      n4s=2;
                      n3d=10;
                      n4p=6;
                      n4d=z-38;
                      n5s=2;
                      }
 else if (z <= 54)    {
                      n2s=2;
                      n2p=6;
                      n3s=2;
                      n3p=6;
                      n4s=2;
                      n3d=10;
                      n4p=6;
                      n4d=10;
                      n5s=2;
                      n5p=z-48;
                      }
 else                 {
                      // [Xe] - configuration  up to [Hg]
                      n2s=2;
                      n2p=6;
                      n3s=2;
                      n3p=6;
                      n4s=2;
                      n3d=10;
                      n4p=6;
                      n4d=10;
                      n5s=2;
                      n5p=6;

                               if (z <= 56) {                   n6s =z-54;  }  //+2
                      else     if (z == 57) {           n5d =2; n6s =1;     }  //+3
                      else     if (z == 58) {n4f =1;    n5d =1; n6s =2;     }  //+4
                      else     if (z == 59) {n4f =3;    n5d =0; n6s =2;     }  //+5
                      else     if (z == 60) {n4f =4;    n5d =0; n6s =2;     }  //+6
                      else     if (z == 61) {n4f =5;    n5d =0; n6s =2;     }  //+7
                      else     if (z == 62) {n4f =6;    n5d =0; n6s =2;     }  //+8
                      else     if (z == 63) {n4f =7;    n5d =0; n6s =2;     }  //+9
                      else     if (z == 64) {n4f =7;    n5d =1; n6s =2;     }  //+10
                      else     if (z == 65) {n4f =9;    n5d =0; n6s =2;     }  //+11
                      else     if (z == 66) {n4f =10;   n5d =0; n6s =2;     }  //+12
                      else     if (z == 67) {n4f =11;   n5d =0; n6s =2;     }  //+13
                      else     if (z == 68) {n4f =12;   n5d =0; n6s =2;     }  //+14
                      else     if (z == 69) {n4f =13;   n5d =0; n6s =2;     }  //+15
                      else     if (z == 70) {n4f =14;   n5d =0; n6s =2;     }  //+16
                      else     if (z == 71) {n4f =14;   n5d =1; n6s =2;     }  //+17
                      else     if (z == 72) {n4f =14;   n5d =2; n6s =2;     }  //+18
                      else     if (z == 73) {n4f =14;   n5d =3; n6s =2;     }  //+19
                      else     if (z == 74) {n4f =14;   n5d =4; n6s =2;     }  //+20
                      else     if (z == 75) {n4f =14;   n5d =5; n6s =2;     }  //+21
                      else     if (z == 76) {n4f =14;   n5d =6; n6s =2;     }  //+22
                      else     if (z == 77) {n4f =14;   n5d =7; n6s =2;     }  //+23
                      else     if (z == 78) {n4f =14;   n5d =9;  n6s =1;    }  //+24
                      else     if (z == 79) {n4f =14;   n5d =10; n6s =1;    }  //+25
                      else
                             {
                             // [Hg] - configuration  up to [Rn]
                                             n4f =14;   n5d =10; n6s =2;       //+26

                              if (z <= 86) {  n6p = z-80;}  //
                              else
                                   {
                                   n6p=6;

                                   //[Rn] - configuration
                                        if (z == 87){                 n7s=1; }  //+1
                                   else if (z == 88){                 n7s=2; }  //+2
                                   else if (z == 89){        n6d=1;   n7s=2; }  //+3
                                   else if (z == 90){        n6d=2;   n7s=2; }  //+4
                                   else if (z == 91){n5f=2;  n6d=1;   n7s=2; }  //+5
                                   else if (z == 92){n5f=3;  n6d=1;   n7s=2; }  //+6
                                   else if (z == 93){n5f=4;  n6d=1;   n7s=2; }  //+7
                                   else if (z == 94){n5f=6;  n6d=0;   n7s=2; }  //+8
                                   else if (z == 95){n5f=7;  n6d=0;   n7s=2; }  //+9
                                   else if (z == 96){n5f=7;  n6d=1;   n7s=2; }  //+10
                                   else if (z == 97){n5f=9;  n6d=0;   n7s=2; }  //+11
                                   else if (z == 98){n5f=10; n6d=0;   n7s=2; }  //+12
                                   else if (z == 99){n5f=11; n6d=0;   n7s=2; }  //+13
                                   else if (z ==100){n5f=12; n6d=0;   n7s=2; }  //+14
                                   else if (z ==101){n5f=13; n6d=0;   n7s=2; }  //+15
                                   else if (z ==102){n5f=14; n6d=0;   n7s=2; }  //+16
                                   else if (z ==103){n5f=14; n6d=0;   n7s=2;  n7p=1; }  //+17
                                   else if (z ==104){n5f=14; n6d=2;   n7s=2;  n7p=0; }  //+18
                                   else if (z ==105){n5f=14; n6d=3;   n7s=2;  n7p=0; }  //+19
                                   else if (z ==106){n5f=14; n6d=4;   n7s=2;  n7p=0; }  //+20
                                   else if (z ==107){n5f=14; n6d=5;   n7s=2;  n7p=0; }  //+21
                                   else if (z ==108){n5f=14; n6d=6;   n7s=2;  n7p=0; }  //+22
                                   else if (z ==109){n5f=14; n6d=7;   n7s=2;  n7p=0; }  //+23
                                   else if (z ==110){n5f=14; n6d=8;   n7s=2;  n7p=0; }  //+24
                                   else if (z ==111){n5f=14; n6d=9;   n7s=2;  n7p=0; }  //+25
                                   else if (z ==112){n5f=14; n6d=10;  n7s=2;  n7p=0; }  //+26
                                   else
                                        {
                                        n5f=14; n6d=10;  n7s=2;  n7p=z-112;
                                        }
                                   }
                             }
                      }

//https://en.wikipedia.org/wiki/Electron_shell


for(int ishell=0; ishell<7; ishell++)
        {
        ShellMap[ishell] =0;
        for(int isubshell=0; isubshell<5; isubshell++)
                                Map[ishell][isubshell]=0;
        }

Map[0][0]=n1s;
Map[1][0]=n2s; Map[1][1]=n2p;
Map[2][0]=n3s; Map[2][1]=n3p; Map[2][2]=n3d;
Map[3][0]=n4s; Map[3][1]=n4p; Map[3][2]=n4d; Map[3][3]=n4f;
Map[4][0]=n5s; Map[4][1]=n5p; Map[4][2]=n5d; Map[4][3]=n5f; //+n5g
Map[5][0]=n6s; Map[5][1]=n6p; Map[5][2]=n6d; ; //+n5g
Map[6][0]=n7s; Map[6][1]=n7p;


int Sum=0;
int MaxShellIndex = -1;



for(int ishell=0; ishell<7; ishell++)
        for(int isubshell=0; isubshell<5; isubshell++)
                                {
                                if(Sum >= ne) Map[ishell][isubshell]=0;
                                else {
                                      int  k = Map[ishell][isubshell];
                                      if(k>0)
                                           {
                                           int k2 = ne-Sum;
                                           if(k2 < k)  {
                                               k = k2;
                                               Map[ishell][isubshell]=k;
                                               }
                                           Sum += k;
                                           ShellMap[ishell] += k;
                                           int NewIndex = 10*ishell+isubshell;
                                           if( MaxShellIndex < NewIndex) MaxShellIndex = NewIndex;
                                           }
                                   }
                                }

//if(Sum!= iq) return -1; // for debugging
if(MaxShellIndex<0)  MaxShellIndex=0;   // 7/11/2021
nshell   = MaxShellIndex/10;
subshell = MaxShellIndex%10;


nK=ShellMap[0];
nL=ShellMap[1];
nM=ShellMap[2];
nN=ShellMap[3];
nO=ShellMap[4];
nP=ShellMap[5];
nQ=ShellMap[6];


nKL   = nK  + nL;
nKLM  = nKL + nM;
nNOPQ = nN+nP+nO+nQ;

const char subshell_name[] = {'s', 'p', 'd', 'f', 'g'};
QString buff;

LastOrbit = buff.asprintf("%1d %c %1d",nshell+1,subshell_name[subshell],Map[nshell][subshell]);

return nshell+1;
};
