//      double Y[1284],YX[1284],WORK[1680000],IWORK[1350],PR[62], P1s[2],P2s[2],P2p[2],PRF[62];
//double Y[1285],YX[1285],WORK[1680001],PR[63], P1s[3],P2s[3],P2p[3],PRF[63];
//int IWORK[1351];
double Y[1285],YX[1285],WORK[27050],PR[63], P1s[3],P2s[3],P2p[7],PRF[63];
int IWORK[7];

//      common /etats/ICO1[63],ICO2[294],ICO3[1284],NUM1[733],NUM2[1911],NUM3[3329];
int ICO1[64],ICO2[295],ICO3[1285],NUM1[734],NUM2[1912],NUM3[3330];

//      common/seceff/sec,cor,secs;
//double sec[34],cor[48],secs[48],seci[48];
double GSEC[35],Gcor[49],Gsecs[49],Gseci[49];

//      common/SecSE/SeSE,StSE;
//double SeSE[66],StSE[12];
double SeSE[67],StSE[13];


//      common/secKLM/C12[209],D12[209],RAD2[209],AKL2[209],PA2[209],
//           AKM3[209],RAD3[209],ALM3[209],C3[209],D3[209],E3[209],DE3[209],
//           PA3[209],PA23[209];

double C12[210],D12[210],RAD2[210],AKL2[210],PA2[210],
           AKM3[210],RAD3[210],ALM3[210],C3[210],D3[210],E3[210],DE3[210],
           PA3[210],PA23[210];

//      common/sec1234/P14[957],C13[957],D13[957],C4M[957],D4M[957],
//                E4M[957],DE4M[957],PA4[957],AKLM4[957],PA123[957];

double P14[958],C13[958],D13[958],C4M[958],D4M[958],
                E4M[958],DE4M[958],PA4[958],AKLM4[958],PA123[958];


//      common/secrad/Rad,ok,ol,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;
double Rad,Rad3s,Rad3p1,Rad3p2,Rad3d,Rad4;

//      common/tol/ep0,ep1,erel,erabs;
double ep0,ep1,eRel,eAbs;

//      common/aug/AKLL,AKLM,ALMM,AM4;
double AKLL,AKLM,ALMM,AM4;

//      common/mean/y1s,y2s,y2p,y3s,y3p,y3d,yl,ym,ymp,ymm,ynn,yn;
//      common/qmoy/ykm,yl1m,yl2m,ym1m,ym2m;
double y1s, y2s, y2p,yL, y3s,y3p,y3d,yM, yMp,yMm, yNn,yN;
double yKm, yL1m, yL2m, yM1m, yM2m;

//      common/popi/n02s,n02p,n03s,n03p,n03d,nl0,nm0;
//int n02s,n02p,n03s,n03p,n03d,nl0,nm0;

//      common/par1/epo,z1,z2,etak,etal1,etal2,etan;
double epo,z1,z2,etak,etal1,etal2,etan;

//      common/par2/etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;
double etam1,etam2,scfk,scfl1,scfl2,scfm1,scfm2,scfn;

//      common/par3/coak,coal1,coal2,coam1,coam2,coan;
double coak,coal1,coal2,coam1,coam2,coan;

//      common/par4/zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,tetan;
double zk,zl1,zl2,zm1,zm2,zn,tetak,tetal1,tetal2,tetam1,tetam2,tetan;

//      common/proj/zp;
//      common/don/zp,E,zt;
double Zp,Ep,Zt;

//      common/corr/ibin;
int ibin;

///      common/parf1s/W;
double   parf1s_W;

//      common/parf2s/W1;
double   parf2s_W1;

//      common/parf2p/W2;
double parf2p_W2;

//      common/parfm1/W1M;
double parfm1_W1M;

//      common/parfm2/W2M;
double parfm2_W2M;

//      common/parfn/W4;
double parfn_W4;


//      common /nn/ n_rn,n_coa,n_scf,n_ll;
double n_rn,n_coa,n_scf;
int n_ll;

int nsp_n;

//      common/bind/Bk,Bl1,Bl2,Bm1,Bm2,Bn;
      double o_Bk, o_Bl1, o_Bl2, o_Bm1, o_Bm2, o_Bn;
//      common/pion_Oleg/
      double o_TCk,o_TCl1,o_TCl2,o_TCm1,o_TCm2,o_TCn;

      //      external F;
// ?????? external F;

//      char nomf*24;
char nomf[25];

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
const double omk[93] = { 0.,
                0.,0.,0.001,0.0012,
                .170E-02,.280E-02,.520E-02,.830E-02,.130E-01,.180E-01,.230E-01,
                .300E-01,.390E-01,.500E-01,.630E-01,.780E-01,.970E-01,.118,
                .140   ,.163        ,.188    ,.214    ,.243    ,.275    ,.308,
                .340   ,.373        ,.406    ,.440    ,.474    ,.507    ,.535,
                .562   ,.589        ,.618    ,.643    ,.667    ,.690    ,.710,
                .730   ,.747        ,.765    ,.780    ,.794    ,.808    ,.820,
                .831   ,.843        ,.853    ,.862    ,.870    ,.877    ,.884,
                .891   ,.897        ,.902    ,.907    ,.912    ,.917    ,.921,
                .925   ,.929        ,.932    ,.935    ,.938    ,.941    ,.944,
                .947   ,.949        ,.951    ,.953    ,.955    ,.957    ,.958,
                .959   ,.961        ,.962    ,.963    ,.964    ,.965    ,.966,
                .967   ,.968        ,.968    ,.969    ,.969    ,.970    ,.970,
                .971   ,.971        ,.972        ,.972  };

const double oml[93] = { 0.,
                 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                .120E-02,.750E-03,.380E-03,.310E-03,.260E-03,.240E-03,.220E-03,
            .270E-03,.330E-03,.840E-03,.150E-02,.260E-02,.370E-02,.500E-02,
            .630E-02,.770E-02,.930E-02,.110E-01,.120E-01,.130E-01,.150E-01,
            .160E-01,.180E-01,.200E-01,.220E-01,.240E-01,.260E-01,.280E-01,
            .310E-01,.340E-01,.370E-01,.400E-01,.430E-01,.460E-01,.490E-01,
            .520E-01,.560E-01,.600E-01,.640E-01,.690E-01,.740E-01,.790E-01,
            .850E-01,.910E-01,.970E-01,.104    ,.111    ,.118    ,.125    ,
                .132  ,.139  ,.147        ,.155    ,.164    ,.174    ,.182    ,
            .192    ,.201    ,.210    ,.220    ,.231    ,.243    ,.255    ,
            .268    ,.281    ,.294    ,.306    ,.320    ,.333    ,.347    ,
            .360    ,.373    ,.386    ,.399    ,.411    ,.424    ,.437    ,
                .450  ,.463  ,.476        ,.489 };

const double Coef_Electron_Holes[35] = {1., // 0 - not used in C
2.,                     // 1
2.,6.,2.,6.,10.,32.,     //2-6
1.,1.,1.,1.,1., 1.,1.,   //8-14
2.,6.,2.,6.,10.,32.,      //15-20
2.,6.,10.,32.,2.,6.,10.,32.,32.,32.,32.,6.,6.,10.  // 21-34
};
                
//c*********************************************************************
//L202:    format(89('*'));
const char *Star89 =
"*****************************************************************************************\n";
//L203:    format(118('*'));
const char *Star118 =
"********************************************************************************************************\n";

/*

//c ....................................................................
//L211:    format(1x,"(ug/cm²) ",10(3x,"(",i2,"+)  "));
//L215:    format(1x,"(ug/cm²) ",10(3x,"(",i2,"+)  "));

//L216:    format('T (ug/cm²)  bare    1s     2s     2p    1s²',
            '    1s2s   1s2p  1s²2s  1s²2p   tot');
//L217:    format('T (ug/cm²)  y1s       y2s      y2p      ym',
            '      yn      Qm      Qm in   Qm out    PTOT');




//L210:    format('Tar.thick.    0e-       1e-       2e-       3e-       4e-       5e-       6e-       7e-       8e-       9e-',
            '     Eout');
//L212:    format('Tar.thick.    10e-      11e-      12e-      13e-      14e-      15e-      16e-      17e-      18e-      19e-');
//L213:    format('Tar.thick.    20e-      21e-      22e-      23e-      24e-      25e-      26e-      27e-      28e-      29e-');
//L214:    format('Tar.thick.    30e-      31e-      32e-      33e-',
            '      34e-      35e-      36e-      37e-      38e-      39e-');
//L218:    format('Tar.thick.    40e-      41e-      42e-      43e-',
            '      44e-      45e-      46e-      47e-      48e-      49e-');
//L219:    format('Tar.thick.    50e-      51e-      52e-      53e-',
            '      54e-      55e-      56e-      57e-      58e-      59e-');

//L220: format(1x,f8.3,1x,f7.3,6(1x,G10.3));
*/

