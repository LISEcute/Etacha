#define ETACHAG_version "Version 4.5.2"
#define ETACHAG_date    "04-APR-2025"


// 4.L 2.122   23-JUN-2016     correction for IonCS  v.23 & 3
// 4.L 2.123   23-JUN-2016     correction for   v.23 & 3   Se2s3s and others between shells  (BTW : somewhere (or 23 or 4) is wrong )
// 4.L 2.124   23-JUN-2016     option UseEnergyLoss
// 4.L 2.125   24-JUN-2016     PR-array, adaption of old versions
// 4.L 2.126   28-JUN-2016     Pentium compilation, double word alignment, EqDif.cpp compilation
// 4.L 2.127   19-JUL-2016     Correction of MainFrame, default Options
// 4.L 2.128   19-JUL-2016     Correction non-use energy loss option
// 4.L 2.129   19-JUL-2016     Correction for go_EnergyL* values
// 4.L 2.130   19-JUL-2016     Energy loss based on Ziegler

// 4.L 2.131   29-MAR-2017     Misspelling in word "Version" (DOS)
// 4.L 2.132   30-MAR-2017     Output files body contains the name of main file
// 4.L 2.133   30-MAR-2017     Output files name contains the name of main file
// 4.L 2.134   30-MAR-2017     Modifications for cell "Beam-q"
// 4.L 2.135   30-MAR-2017     Modifications for cell "Beam-A"
// 4.L 2.136   31-MAR-2017     Output file for EXCEL
// 4.L 2.137   19-APR-2019     Change of default parameters, Zero cross section at beginning

// Qt porting
// 4.3 1       01-APR-2021     Ksenia's work
// 4.3 2       29-APR-2021     Ksenia's work
// 4.3 3       30-APR-2021     Oleg's update
// 4.3 4       23-JUN-2021     Global update
// 4.3 5       11-JUL-2021     Global update
// 4.3 6       13-JUL-2021     Global update
// 4.3.7       26-JUL-2021     Global update
// 4.3.8       02-AUG-2021     Global update -- relase version to include in LISE package
// 4.3.9       09-AUG-2021     Graphs have been completed
// 4.3.10      10-AUG-2021     benchmarks
// 4.3.11      11-AUG-2021     updates for liseStrcpyOS.h
// 4.3.12      11-AUG-2021     modifications for debug messages and icons
// 4.3.13      13-SEP-2021     q-max value in logs
// 4.3.14      21-SEP-2021     Modification in PopMean  (it was disagreement for etacha4 option)
// 4.3.15      21-SEP-2021     Printing and plotting Qmean and Qprob
// 4.3.16      22-SEP-2021     Shell evolution -- QChartArea instead QSpline

// 4.3.17      29-SEP-2021     4 bugs were fixed by Toshi Sumikama (RIKEN)
// 4.3.18      29-SEP-2021     function "NormOlegSum" became VOID
// 4.4.01      29-SEP-2021     no more BETA sign!!  Some initial default values were changed
// 4.4.02      30-SEP-2021     Energy value bug in PION was fixed by Toshi Sumikama (RIKEN)
// 4.4.03      30-SEP-2021     Similar modifications are necessary for
//                             e_SEIK.cpp, e_Sex2.cpp, e_Sexi.cpp, e_Snl.cpp, e_Tceis.cpp.  Toshi Sumikama (RIKEN)
// 4.4.04      25-OCT-2021     ETACHA page, links have been changed
// 4.4.05       12/03/2021     modified for non-latin paths
// 4.4.06       12/04/21       read/write modification for const QString&
// 4.4.07       12/29/21       compatability with Linux
// 4.4.08       01/24/22       set filename at reading

// 4.4.09       05/30/22       Plot option in etacha file  // HW's requests
// 4.4.10       05/30/22       Correction with density value input from file
// 4.4.11       05/30/22       batch mode with command string argument  -r

// 4.4.12       07/13/22       Migrating to Qt63
// 4.4.13       07/27/22       updating title for graph window
// 4.4.14       10/16/22       Evolution plot start
// 4.4.15       10/16/22       close graph windows at exit
// 4.4.16       10/16/22       Evolution plot full completion
// 4.4.17       03/18/23       More digits for sum in the Shell window
// 4.4.18       07/01/23       modificications, compiled with MSVC

// 4.5.0        04/25/25       Arjus's optimization
// 4.5.1        04/25/25       ETACHA site haas been changed in About dialog
// 4.5.2        04/25/25       the status bar shows ETACHA version at Wecome
// 4.5.3        04/25/25       Ar0j1, Arj2

