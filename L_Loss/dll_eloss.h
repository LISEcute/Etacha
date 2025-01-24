//----------------------------------------------------------------------------
// ObjecQWidgets
// Copyright (c) 1991, 1995 by Borland International, All Rights Reserved
//----------------------------------------------------------------------------
double  SimpleEnergyResidue(double Zp, double Ap,  double Zt, double At, double E, double thick);
double  ZieglerStopping(int Zp,  double E, int Zt, int SOLIDGAS);
double  NuclearStopping(int zp, double Mp,  int zt, double Mt, double e);
double  HubertStopping(int Zp, double Mp, int Zt, double Mt, double e, int SOLIDGAS);
int  GetRangeBinMass(int Z);
double  GetRangeBinR(int Z, int Zmat);

