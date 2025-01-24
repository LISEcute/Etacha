#ifndef ATOMIC_SHELL_H
#define ATOMIC_SHELL_H
#include <QString>

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class atom_shell
{
public:

atom_shell();
atom_shell  (int iz, int iq=0);
int input_zq(int iz, int iq=0);

int nK,nL,nM,nN,nO,nP,nQ,nKL,nKLM,nNOPQ;
int Map[7][5];    //shell =7, subshell=5
int ShellMap[7];  //shell =7, subshell=5


int nshell;
int subshell;

QString LastOrbit;

int n_1s() {return Map[0][0];};

int n_2s() {return Map[1][0];};
int n_2p() {return Map[1][1];};

int n_3s() {return Map[2][0];};
int n_3p() {return Map[2][1];};
int n_3d() {return Map[2][2];};

int n_4s() {return Map[3][0];};
int n_4p() {return Map[3][1];};
int n_4d() {return Map[3][2];};
int n_4f() {return Map[3][3];};

int n_5s() {return Map[4][0];};
int n_5p() {return Map[4][1];};
int n_5d() {return Map[4][2];};
int n_5f() {return Map[4][3];};

int n_6s() {return Map[5][0];};
int n_6p() {return Map[5][1];};
int n_6d() {return Map[5][2];};

int n_7s() {return Map[6][0];};
int n_7p() {return Map[6][1];};

private :

int ne;
int z;

int n1s;

int n2s;
int n2p;

int n3s;
int n3p;
int n3d;

int n4s;
int n4p;
int n4d;
int n4f;

int n5s;
int n5p;
int n5d;
int n5f;

int n6s;
int n6p;
int n6d;

int n7s;
int n7p;
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#endif
