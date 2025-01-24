//      double Y[1284],YX[1284],WORK[1680000],IWORK[1350],PR[62], P1s[2],P2s[2],P2p[2],PRF[62];
//extern double Y[1285],YX[1285],WORK[1680001],PR[63], P1s[3],P2s[3],P2p[3],PRF[63];
//extern int IWORK[1351];
extern double Y[1285],YX[1285],WORK[27050],PR[63], P1s[3],P2s[3],P2p[7],PRF[63];  // bug found by Toshi Sumikama 09/28/2021
extern int IWORK[7];

//      common /etats/ICO1[63],ICO2[294],ICO3[1284],NUM1[733],NUM2[1911],NUM3[3329];
extern int ICO1[64],ICO2[295],ICO3[1285],NUM1[734],NUM2[1912],NUM3[3330];

//      common/seceff/sec,cor,secs;
//double sec[34],cor[48],secs[48],seci[48];
extern double GSEC[35],Gcor[49],Gsecs[49],Gseci[49];

//      common/SecSE/SeSE,StSE;
//double SeSE[66],StSE[12];
extern double SeSE[67],StSE[13];

//      common/secKLM/C12[209],D12[209],RAD2[209],AKL2[209],PA2[209],AKM3[209],RAD3[209],ALM3[209],C3[209],D3[209],E3[209],DE3[209], PA3[209],PA23[209];
extern double C12[210],D12[210],RAD2[210],AKL2[210],PA2[210],AKM3[210],RAD3[210],ALM3[210],C3[210],D3[210],E3[210],DE3[210], PA3[210],PA23[210];

//      common/sec1234/P14[957],C13[957],D13[957],C4M[957],D4M[957],E4M[957],DE4M[957],PA4[957],AKLM4[957],PA123[957];
extern double P14[958],C13[958],D13[958],C4M[958],D4M[958],E4M[958],DE4M[958],PA4[958],AKLM4[958],PA123[958];

//      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
extern double Rad,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;

//      common/tol/ep0,ep1,erel,erabs;
extern double ep0,ep1,eRel,eAbs;

//      common/aug/AKLL,AKLM,ALMM,AM4;
extern double AKLL,AKLM,ALMM,AM4;

//      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
//      common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m;
extern double y1s, y2s,y2p,yL, y3s,y3p,y3d,yM, yMp,yMm, yNn,yN;
extern double yKm, yL1m,yL2m,  yM1m,yM2m;

//      common/popi/n02s,n02p,n03s,n03p,n03d,nl0,nm0;
//extern int n02s,n02p,n03s,n03p,n03d,nl0,nm0;

//      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
extern double epo,z1,z2,etak,etal1,etal2,etan;

//      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
extern double etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;

//      common/par3/coak,coal1,coal2,coam1,coam2,coan;
extern double coak,coal1,coal2,coam1,coam2,coan;

//      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,tetan;
extern double zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,tetan;


//      common/proj/zp;
//      common/don/zp,Ep,zt;
extern double Zp,Ep,Zt;

//      common/corr/ibin;
extern int ibin;


//     common/parf1s/W;
extern double   parf1s_W;

//      common/parf2s/W1;
extern double   parf2s_W1;

//      common/parf2p/W2;
extern double parf2p_W2;

//      common/parfm1/W1M;
extern double parfm1_W1M;

//      common/parfm2/W2M;
extern double parfm2_W2M;

//      common/parfn/W4;
extern double parfn_W4;

//      common/pion_Oleg/
extern double o_TCk,o_TCl1,o_TCl2,o_TCm1,o_TCm2,o_TCn;
//      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn;
extern double o_Bk,o_Bl1,o_Bl2,o_Bm1,o_Bm2,o_Bn;


//      common /nn/ n_rn,n_coa,n_scf,n_ll;
extern double n_rn,n_coa,n_scf;
extern int n_ll;

//	common /nsp/ n
extern int nsp_n;

extern int EtachaVersion;
extern int ibinParameter;
extern int IonizationModel;
extern int ExcitationModel;
