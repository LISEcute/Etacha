#include "e_mainwindow.h"
#include "ui_e_mainwindow.h"

#include <QFileDialog>
#include <QFile>
#include <QDir>
#include <QString>
#include <QScreen>
#include <QMessageBox>
#include <QSignalMapper>
#include <QTextStream>
#include <QKeyEvent>
#include <QDebug>
#include <QElapsedTimer>
#include <QDateTime>
#include <QTimer>
#include <QFontDatabase>


#include "e_Constant.h"
#include "e_ftype.h"
#include "../source/e_TargetData.h"
#include "../source/e_AtomicShell.h"
#include "../source/e_Etacha4.h"
#include "w_Stuff/liseStrcpyOS.h"

#include "e_CSresults.h"


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWW
extern QString numberKillZero(double v, char ch,  int k, int m=0); // only for 'f'
extern QString numberKillExpZero(double v, int k, bool strong=true);
extern int EraseExtention(QString &fileName);
extern QString getShortFileName(const QString &source, int opt=0);
extern void ELEMENT(double &Z, double &A, char *CZ, int IOPT, int &IRC );
//extern double EnergyResidueRng(double Zp, double Ap,  double Zt, double At, double E, double thick );
extern double EnergyResidueSP(double Zp, double Ap,  double Zt, double At, double E, double thick );
extern double StoppingPower( int Zp, int Ap, int Zt, int At, double E);
extern double Velocity_au(double E);
extern double pow2(double par);
extern QString ElapsedTime(qint64 msec);

extern void fastAppend(QTextEdit *editWidget, const QString& message);

//---------------------------------------------
extern QString ResultsDir;
extern QString localPATH;
extern int fontsizeGlobal;
extern const char *FileNameAbsent;
extern QString FileArg;

//char DecimalPoint='.';

//--------------------------------  ETACHA
double gAb=207, gZb=82, gQb=64, gEnergy=28.9;
double gAt=12,gZt=6,gThick=1,gDensity=2.25;
double gMinStepTarget=5,gMaxStepTarget=250;
double gO_EnergyL1=28.9,gO_EnergyL2=27.9;
double gO_dEdX1=78.819,gO_dEdX2=76.998;
double gUncertainAbs=1e-12,gUncertainRel=1e-5;

double ProjectileVelocity;     /// au
double ElectronVelocity;       /// m/s

int EtachaVersion = etacha_v3;
int ibinParameter = 0;
int IonizationModel = 0;      /// 1 - PWBA
int ExcitationModel = 0;      /// 1 - PWBA
int ShowItermediate = 1;
int ShowResult = 0;
int DifEqModel = 0 ;          ///  0 - ODE, 1 - RKB45

int UseEloss = 1;
int RunningFlag = -1;
bool GlobalBreak=false;
extern bool RunBatch;
//--------------------------------  ETACHA   }
/*                  //Original comments
enum enum_file {
en_log,
en_f09,
en_f19,
en_f29,
en_f39,
en_f49,
en_f59,
en_fPied,
en_fMean,
en_fSEff,
en_fExcel
};
*/

const char *eta_filenames[] = {
  "EtaLog.txt",
  "Eta0009.txt",
  "Eta1019.txt",
  "Eta2029.txt",
  "Eta3039.txt",
  "Eta4049.txt",
  "Eta5059.txt",
  "EtaPied.txt",
  "PopMean.txt",
  "SecEff.txt",
  "forExcel.txt" };

double gEnergyF;

//extern int e_etacha(QWidget *tobject, const QString *LinitialDir, const char *LfileName);

ETACHA *etacha=nullptr;
char BufRtf[1000];
atom_shell Zshell;
atom_shell Qshell;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  FlagPermit = false;
  ZPedit_permit=true;
  ZTedit_permit=true;

  Modified=false;
  etacha = new ETACHA(this);
  //-----------------------------------------------------------
  setWindowFlags( Qt::Window | Qt::CustomizeWindowHint |
                  Qt::WindowTitleHint | Qt::WindowSystemMenuHint |
                  Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint );

  ui->setupUi(this);
  //-----------------------------------------------------------
  QButtonGroup *bg;
  //================================================
  bg_Version = new QButtonGroup(0);     bg = bg_Version;     bg->setExclusive(true);
  bg->addButton(ui->rB_23);
  bg->addButton(ui->rB_3);
  bg->addButton(ui->rB_34);
  bg->addButton(ui->rB_4);
  bg->addButton(ui->rB_45);

  int k=0;
  const QList<QAbstractButton *> buttons = bg->buttons();
  for (QAbstractButton *myButton : buttons) bg->setId(myButton,k++);

  connect(bg_Version, QOverload<int>::of(&QButtonGroup::idClicked),
          this, &MainWindow::bg_Version_clicked);

  ui->rB_34->setEnabled(false);
  //================================================
  bg_ibin = new QButtonGroup(0);     bg = bg_ibin;     bg->setExclusive(true);
  bg->addButton(ui->rB_corEmpirical);
  bg->addButton(ui->rB_corBinding);
  bg->addButton(ui->rB_corNo);

  k=0;
  const QList<QAbstractButton *> buttons2 = bg->buttons();
  for (QAbstractButton *myButton : buttons2) bg->setId(myButton,k++);
  connect(bg, QOverload<int>::of(&QButtonGroup::idClicked),
          this, &MainWindow::bg_ibin_clicked);
  //================================================

  bg_Ionization = new QButtonGroup(0);     bg = bg_Ionization;     bg->setExclusive(true);
  bg->addButton(ui->rb_IC);
  bg->addButton(ui->rb_IP);

  k=0;
  const QList<QAbstractButton *> buttons3 = bg->buttons();
  for (QAbstractButton *myButton : buttons3) bg->setId(myButton,k++);
  connect(bg, QOverload<int>::of(&QButtonGroup::idClicked),
          this, &MainWindow::bg_Ionization_clicked);

  //================================================

  bg_Excitation = new QButtonGroup(0);     bg = bg_Excitation;     bg->setExclusive(true);
  bg->addButton(ui->rb_ES);
  bg->addButton(ui->rb_EP);

  k=0;
  const QList<QAbstractButton *> buttons4 = bg->buttons();
  for (QAbstractButton *myButton : buttons4) bg->setId(myButton,k++);
  connect(bg, QOverload<int>::of(&QButtonGroup::idClicked),
          this, &MainWindow::bg_Excitation_clicked);

  //================================================

  bg_Integration = new QButtonGroup(0);     bg = bg_Integration;     bg->setExclusive(true);
  bg->addButton(ui->rb_ODE);
  bg->addButton(ui->rb_RKF45);
  bg->addButton(ui->rb_EM);
  bg->addButton(ui->rb_EQ);
  k=0;
  const QList<QAbstractButton *> buttons5 = bg->buttons();
  for (QAbstractButton *myButton : buttons5) bg->setId(myButton,k++);
  connect(bg, QOverload<int>::of(&QButtonGroup::idClicked),
          this, &MainWindow::bg_Integration_clicked);

  ui->rb_EM->setEnabled(true);
  ui->rb_EQ->setEnabled(true);
  //---------------------------------------------------------------------------

  cb_Show = ui->cb_ShowResults;
  connect(cb_Show, SIGNAL(currentIndexChanged(int)), this, SLOT(cb_Show_clicked()));

  Plots=true;
  PlotEvolution=true;
  ui->check_Plots->setChecked(Plots);
  ui->check_EvolutionPlot->setChecked(PlotEvolution);
  connect(ui->check_Plots, SIGNAL(clicked()), this, SLOT(plots_clicked()));
  connect(ui->check_EvolutionPlot, SIGNAL(clicked()), this, SLOT(plotsEvolution_clicked()));

  connect(ui->label_Construction, SIGNAL(clicked()), this, SLOT(cmLinkConstruction()));

  //    ui->check_Debug->setEnabled(false);
  //-------

  //ffffffffffffffffffffffffffffffffffffffffffffffffff    font family start
  const QStringList myFontNameList = {"Consolas", "Lucida Console", "Monaco", "Courier New",  "Courier"};


  const QStringList fontFamilies = QFontDatabase::families();

  for ( const auto& item : myFontNameList)
    {
      bool res = fontFamilies.contains(item);
      if(res)
        {
          //  qDebug() << item;
          //font.setFamily(item);
          font = QFont(item,fontsizeGlobal);
          // QFont font2(item,fontsizeGlobal-2);
          ui->text_Main->setCurrentFont(font);
          ui->text_Shells->setCurrentFont(font);
          break;
        }
    }
  //ffffffffffffffffffffffffffffffffffffffffffffffffff    font family stop
  ui->pB_DOC_version->setHidden(true);

  FlagPermit = true;

  FFileName  = FileArg.size()>0? FileArg : localPATH + FileNameAbsent;
  if(FileArg.size()>0) readFile(FFileName);
  else SetPage();

  if(RunBatch)
    QTimer::singleShot(100, this, &MainWindow::on_SB_Run_clicked);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
MainWindow::~MainWindow()
{
  delete ui;
  delete etacha;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::closeEvent(QCloseEvent *e)
{
  for (auto *item : m_charts)
      {
      if(item)
        item->close();
      }

  m_charts.clear();

  e->accept();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::SetPage(bool fromFile)
{ 
  double gAb_save=gAb;
  double gAt_save=gAt;

  ZPedit_permit=false;
  ui->edit_BeamZ->setText(numberKillZero(gZb, 'f', 0));
  ui->edit_BeamQ->setText(numberKillZero(gQb, 'f', 0));
  ZPedit_permit=true;

  makeBeamZ();

  gAb=gAb_save;
  ui->edit_BeamA ->setText(numberKillZero(gAb,'f',2));
  ui->edit_Energy->setText(numberKillZero(gEnergy, 'f', 3));

  //-----------------------------------------
  ZTedit_permit=false;
  ui->edit_TargZ->setText(numberKillZero(gZt, 'f', 0));
  ui->edit_Thick->setText(numberKillZero(gThick, 'f', 4));
  ui->edit_Density->setText(numberKillZero(gDensity, 'f', 4));
  ZTedit_permit=true;

  makeTargetZ(fromFile);

  gAt=gAt_save;
  ui->edit_TargA->setText(numberKillZero(gAt, 'f', 2));
  makeTargetThick();

  ui->edit_MinStep->setText(numberKillZero(gMinStepTarget, 'f', 3));
  ui->edit_MaxStep->setText(numberKillZero(gMaxStepTarget, 'f', 3));

  ui->editUncertainAbs->setText(numberKillExpZero(gUncertainAbs,2));
  ui->editUncertainRel->setText(numberKillExpZero(gUncertainRel,2));

  emit bg_Version     ->idClicked(EtachaVersion);
  emit bg_ibin        ->idClicked(ibinParameter);
  emit bg_Ionization  ->idClicked(IonizationModel);
  emit bg_Excitation  ->idClicked(ExcitationModel);
  emit bg_Integration ->idClicked(DifEqModel);

  cb_Show->setCurrentIndex(ShowResult);

  ui->check_CS->setChecked(ShowItermediate);
  ui->check_Plots->setChecked(Plots);
  ui->check_EvolutionPlot->setChecked(PlotEvolution);
  ui->check_EnergyLoss->setChecked(UseEloss);
  ShowEnergyLossCells();


  ui->text_Main->clear();
  ui->text_Shells->clear();
  QString firstPage =
      "<center><font color=\"blue\"><br><br>This version of ETACHA is for ions with up to 60 electrons </font><BR><br>"
"<i><b>Ready to modify input parameters and run</b></i></center>";
  QString welc = "Welcome to ETACHA (";
  welc += QString(ETACHAG_version) + " from " + QString(ETACHAG_date) + ")";
  statusBar()->showMessage(welc);

  ui->text_Main->setHtml(firstPage);

  Modified=false;
  setFileName();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void MainWindow::ReadPage()
{
  const char *cc = "Warning: wrong initial settings";

  FlagReadPage = false;

  if(gEnergy <0.1 || gEnergy > 1000) {
      QMessageBox::information(this,cc, "Energy should be more than 0.1 MeV/u\n"
                                          "and less than 1000 MeV/u!");
      return;
    }

  if(gEnergyF <=0 && UseEloss)
    {
      QMessageBox::information(this, "Zero final energy is 0!","Check initial energy and target thickness!");
      return;
    }

  if(gZb< 1 || gZb >96) {
      QMessageBox::information(this,cc,"Check Z of projectile");
      return;
    }

  gQb = ui->edit_BeamQ->text().toDouble();
  if(gQb<0 || gQb>gZb ) {
      QMessageBox::information(this,cc,"Check Q of projectile");
      return;
    }

  if(gThick<=0){
      QMessageBox::information(this, cc, "Thickness should be positive!");
      return;
    }

  FlagReadPage=true;

  gAb             = ui->edit_BeamA->text().toDouble();
  gDensity        = ui->edit_Density->text().toDouble();
  gMinStepTarget  = ui->edit_MinStep->text().toDouble();
  gMaxStepTarget  = ui->edit_MaxStep->text().toDouble();
  gUncertainAbs   = ui->editUncertainAbs->text().toDouble();
  gUncertainRel   = ui->editUncertainRel->text().toDouble();
  ShowResult = cb_Show->currentIndex();
  ShowItermediate = ui->check_CS->isChecked();
  DebugMode = ui->check_Debug->isChecked();
  //===================================================================
  /*bool flagStep=false;
if(flagStep){QMessageBox::information(this, cc,"Number of Steps > 0 and < 100!",);}
*/
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_Energy_textChanged(const QString &arg)
{
  gEnergy = arg.toDouble();
  Modified=true; setFileName();
  CalculateEnergy();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_BeamZ_textChanged(const QString &arg)
{
  if(!ZPedit_permit) return;
  if(arg.length()==0) return;
  gZb = arg.toInt();
  makeBeamZ();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::makeBeamZ()
{
  strcpyL(CZ,3,"**");
  int IRC1;

  if(gZb>0 && gZb<=97)
    {
      ELEMENT(gZb, gAb, CZ, 3, IRC1 );

      if(IRC1 == -1)
        QMessageBox::information(this, "Error!","Z <-> Element conversion");
    }

  bool ZPedit_permit_save = ZPedit_permit;
  ZPedit_permit=false;

  ui->edit_BeamEl->setText(CZ);
  ui->edit_BeamA->setText(numberKillZero(gAb,'f',2));

  ZPedit_permit=ZPedit_permit_save;

  makeBeamQ();
  Modified=true; setFileName();
  CalculateEnergy();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_BeamQ_textChanged(const QString &arg)
{
  if(!ZPedit_permit) return;
  if(arg.length()==0) return;
  gQb = arg.toDouble();
  makeBeamQ();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_BeamA_textChanged(const QString &arg)
{   gAb = arg.toDouble(); CalculateEnergy();}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_BeamEl_textChanged(const QString &arg)
{
  if(!ZPedit_permit) return;

  if(arg.length()==0) return;

  strcpyL(CZ, 3, arg.toStdString().c_str());
  if(strlen(CZ)==1) strcatL(CZ,3," ");

  int IRC1;

  ELEMENT(gZb, gAb, CZ, 4, IRC1 );

  ZPedit_permit=false;
  if(IRC1==0){
      ui->edit_BeamZ->setText(numberKillZero(gZb, 'f', 0));
      ui->edit_BeamA->setText(numberKillZero(gAb, 'f', 2));
    }
  else    {
      ui->edit_BeamZ->text()  = "**";
      ui->edit_BeamA->text()  = "**";
    }

  ZPedit_permit=true;
  makeBeamQ();
  Modified=true; setFileName();
  CalculateEnergy();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::makeBeamQ()
{
  if(gQb<0 || gQb>gZb) {
      gQb = gZb;
      ui->edit_BeamQ->setText(numberKillZero(gQb, 'f' ,0));
    }

  Qshell.input_zq(gZb,gQb);
  Zshell.input_zq(gZb,  0);
  ui->label_Zorbital->setText(Zshell.LastOrbit);
  ui->label_Qorbital->setText(Qshell.LastOrbit);
  PertrubationParameters();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_TargZ_textChanged(const QString &arg)
{
  if(!ZTedit_permit) return;

  if(arg.length()==0) return;

  gZt = arg.toInt();
  makeTargetZ();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::makeTargetZ(bool fromFile)
{
  strcpy(CZ,"**");
  int IRC1;
  int gZt_den = gZt;

  if(gZt>0 && gZt<=97) {
      ELEMENT(gZt, gAt, CZ, 1, IRC1 );

      if(IRC1 == -1)
        QMessageBox::information(this, "Error!","Z <-> Element conversion");
    }
  else gZt_den=4;

  ZTedit_permit=false;
  ui->edit_TargEl->setText(CZ);
  ui->edit_TargA->setText(numberKillZero(gAt, 'f', 2));

  if(!fromFile)
    {
      gDensity = el_densite[gZt_den];
      ui->edit_Density->setText(numberKillZero(gDensity, 'f', 4));
    }

  ZTedit_permit=true;
  Modified=true; setFileName();
  makeTargetThick();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_TargEl_textChanged(const QString &arg)
{
  if(!ZTedit_permit) return;

  if(arg.length()==0) return;

  strcpyL(CZ,3,arg.toStdString().c_str());
  if(strlen(CZ)==1) strcatL(CZ,3," ");

  int IRC1;

  int gZt_save=gZt; int gAt_save=gAt;

  ELEMENT(gZt, gAt, CZ, 2, IRC1 );

  if(gZt<=0 && gZt>97)
    {
      gZt = gZt_save;  gAt = gAt_save;
      IRC1=-1;
    }
  //------------------------------------------
  ZTedit_permit=false;
  if(IRC1==0){
      ui->edit_TargZ->setText(numberKillZero(gZt, 'f', 0));
      ui->edit_TargA->setText(numberKillZero(gAt, 'f',  3));         //check it, was ("0.###",gAt);
      gDensity = el_densite[int(gZt)];
      ui->edit_Density->setText(numberKillZero(gDensity, 'f' ,4));
    }
  else    {
      ui->edit_TargZ->setText("**");
      ui->edit_TargA->setText("**");
      ui->edit_Density->setText("**");
    }
  ZTedit_permit=true;

  Modified=true; setFileName();
  makeTargetThick();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_TargA_textChanged(const QString &arg)
{
  if(!ZTedit_permit) return;
  gAt = arg.toDouble();
  makeTargetThick();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::makeTargetThick(void)
{
  if(!ZTedit_permit) return;
  double TTHICK_A = gThick /(gAt*U_FC);             /// atoms/cm2;
  ui->label_TargetAtoms->setText(numberKillExpZero(TTHICK_A,3));
  CalculateEnergy();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::CalculateEnergy()
{
  //gEnergyF = EnergyResidue(gZb, gAb, gZt, gAt, gEnergy, gThick);    //original comments

  if(UseEloss) {
      gEnergyF  = EnergyResidueSP(gZb, gAb, gZt, gAt, gEnergy, gThick);
      //            double gEnergyF2 = EnergyResidueRng(gZb, gAb, gZt, gAt, gEnergy, gThick);
      //            gEnergyF = qMax (gEnergyF2,gEnergyF);
    }
  else         gEnergyF = gEnergy;

  // test 207Pb 28.90 target-C 50 mg/cm2 //
  // LISE EL = 4     ER=5.6MeV/u
  // LISE EL = 1     ER=6.9MeV/u
  // LISE EL = 0     ER=3.5MeV/u
  // using SP = 6.91
  // using Range = 0

  //--------------------------------------------------
  ProjectileVelocity = Velocity_au(gEnergy);  /// in a.u.
  ElectronVelocity =  gZb;                    /// in au

  if(gEnergyF>0)
    {
      gO_EnergyL1=gEnergy;
      gO_EnergyL2=gEnergyF;

      gO_dEdX1 = StoppingPower( gZb, gAb, gZt, gAt, gEnergy);
      gO_dEdX2 = StoppingPower( gZb, gAb, gZt, gAt, gEnergyF);
      ui->label_EnergyFinal->setText(numberKillZero(gEnergyF, 'f', 3));

      ui->label_Sinit ->setText(numberKillZero(gO_dEdX1, 'f', 3));
      ui->label_Sfinal->setText(numberKillZero(gO_dEdX2, 'f', 3));

      ui->label_Vp->setText(numberKillZero(ProjectileVelocity, 'f', 3));
    }
  else    {
      ui->label_EnergyFinal->setText("0");
      ui->label_Sinit->setText(" ");
      ui->label_Sfinal->setText("0");
      ui->label_Vp->setText("-");
    }

  PertrubationParameters();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::PertrubationParameters()
{
  double t = qMax(1e-5,ProjectileVelocity*gZb);
  K_perturbation = gZt * ElectronVelocity/ t;
  K_perturbationN = K_perturbation / pow2(Qshell.nshell+1);

  ui->label_Kp_at_n->setText(tr("Kp (n=%1) =").arg(Qshell.nshell+1));

  if(gEnergyF>0)
    {
      ui->label_Kp1->setText(numberKillZero(K_perturbation , 'f', 3));
      ui->label_Kpn->setText(numberKillZero(K_perturbationN, 'f', 3));
    }
  else    {
      ui->label_Kp1->setText("-");
      ui->label_Kpn->setText("-");
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_check_CS_clicked()
{   ShowItermediate = ui->check_CS->isChecked();}
//wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
/*void MainWindow::buttonGroupClicked(QAbstractButton *button) // keep it
{
    EtachaVersion = bg_Version->id(button);

       const QList<QAbstractButton *> buttons = bg_Version->buttons();
       for (QAbstractButton *myButton : buttons) {
           if (myButton != button)  button->setChecked(false);
       }
}*/
//---------------------------------------------------------------------------
void MainWindow::bg_Version_clicked(int k)
{
  EtachaVersion = k; Modified=true; setFileName();
  k=0;
  const QList<QAbstractButton *> buttons = bg_Version->buttons();
  for (QAbstractButton *myButton : buttons)
    myButton->setChecked(k++ == EtachaVersion);
}
//---------------------------------------------------------------------------
void MainWindow::bg_ibin_clicked(int k)
{
  ibinParameter = k; Modified=true; setFileName();
  k=0;
  const QList<QAbstractButton *> buttons = bg_ibin->buttons();
  for (QAbstractButton *myButton : buttons)
    myButton->setChecked(k++ == ibinParameter);
}
//---------------------------------------------------------------------------
void MainWindow::bg_Ionization_clicked(int k)
{
  IonizationModel = k; Modified=true; setFileName();
  k=0;
  const QList<QAbstractButton *> buttons = bg_Ionization->buttons();
  for (QAbstractButton *myButton : buttons)
    myButton->setChecked(k++ == IonizationModel);
}
//---------------------------------------------------------------------------
void MainWindow::bg_Excitation_clicked(int k)
{
  ExcitationModel = k; Modified=true; setFileName();
  k=0;
  const QList<QAbstractButton *> buttons = bg_Excitation->buttons();
  for (QAbstractButton *myButton : buttons)
    myButton->setChecked(k++ == ExcitationModel);
}
//---------------------------------------------------------------------------
void MainWindow::bg_Integration_clicked(int k)
{
  DifEqModel = k; Modified=true; setFileName();
  k=0;
  const QList<QAbstractButton *> buttons = bg_Integration->buttons();
  for (QAbstractButton *myButton : buttons)
    myButton->setChecked(k++ == DifEqModel);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::cb_Show_clicked()
{
  ShowResult = cb_Show->currentIndex();
  if(RunningFlag!=0) return ;   ///  -1 start; 1 - running; 0 - after calc
  Results();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_check_EnergyLoss_clicked()
{
  UseEloss = ui->check_EnergyLoss->isChecked();
  ShowEnergyLossCells();
  Modified=true; setFileName();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::ShowEnergyLossCells()
{
  CalculateEnergy();
  ui->label_Sinit->setVisible(UseEloss);
  ui->label_Sfinal->setVisible(UseEloss);
  ui->label_Stopping_stat->setVisible(UseEloss);

  ui->label_EnergyFinal->setVisible(UseEloss);
  ui->label_Ei_stat->setVisible(UseEloss);
  ui->label_Ef_stat->setVisible(UseEloss);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::keyPressEvent(QKeyEvent *e)
{
  if(e->key() == Qt::Key_Escape) GlobalBreak=true;
  //QMessageBox::information(this,"escape", "pressed");
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_Density_textChanged(const QString &arg1)
{
  gDensity = arg1.toDouble();
  Modified=true; setFileName();
  makeTargetThick();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_edit_Thick_textChanged(const QString &arg1)
{
  gThick = arg1.toDouble();
  Modified=true; setFileName();
  makeTargetThick();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW



//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWW                                                     WWWWWWWWWWWWWWW
//WWWWWWW             R U N      H E R E                      WWWWWWWWWWWWWWW
//WWWWWWW                                                     WWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_pb_Run2_clicked() {on_SB_Run_clicked();}
void MainWindow::on_actionRun_ETACHA_triggered(){on_SB_Run_clicked();}
//----------------------------------------------------------------------
void MainWindow::on_SB_Run_clicked()
{
  static int att = 0;
  static int attE = 0;

  ui->text_Main->clear();
  ui->text_Shells->clear();

  QString result =
      "<center><font color=\"blue\"><br><br><b>The \"Run\" button has been clicked</b></font><BR><br>"
"<i>please, wait .... </i></center>";


  ui->text_Main->setHtml(result);

  statusBar()->showMessage("Cross-sections part");

  //ui->text_Shells->setText("One moment please .....");

  ReadPage();

  //double er=gEnergy;

  if(FlagReadPage)
    {
      if(gEnergyF<=0 && UseEloss)
        {
          statusBar()->showMessage("Ready");
          ui->text_Main->clear();
          ui->text_Shells->clear();
          QMessageBox::information(this,"Zero final energy!",
                                   "Check initial energy and target thickness!");
          return;
        }

      QDateTime timeSr = QDateTime::currentDateTime();

      result = timeSr.toString();
      ui->text_Main->setHtml(result);

      statusBar()->showMessage("One moment please .....");
      RunningFlag=1;
      etacha->init();
      fastAppend(ui->text_Main, "<p> F4 initialization <br> Donaut initialization </p>");

      //-----------------------------------  init   OLD -> NEW  end
      QElapsedTimer timer;
      timer.start();

      etacha->DONAUT();

      result = "<p>Donaut finished. <br>" + ElapsedTime(timer.elapsed()) + "</p>"; // msec
      fastAppend(ui->text_Main, result);

      if(ShowItermediate) {
          (new FormResults())->exec();
        }

      if(Plots)
        {
          currentGraph = new e_graph();
          setPlotGeometry(currentGraph);
          QString mTitle="ETACHA4 charts #" + QString::number(++att);
          mTitle += " :  " + ui->edit_BeamA->text()+ ui->edit_BeamEl->text() + " (" + ui->edit_Energy->text() +
                    " MeV/u) + " + ui->edit_TargEl->text();
          currentGraph->setWindowTitle(mTitle);
          currentGraph->show();
          m_charts.append(currentGraph);
      }
      else currentGraph=nullptr;

      if(PlotEvolution)
        {
          currentEvolutionGraph = new e_graphEvolution();
          setPlotGeometryE(currentEvolutionGraph);
          QString mTitle="ETACHA4 Evolution plot #" + QString::number(++attE);
          mTitle += " :  " + ui->edit_BeamA->text()+ ui->edit_BeamEl->text() + " (" + ui->edit_Energy->text() +
                    " MeV/u) + " + ui->edit_TargEl->text();
          currentEvolutionGraph->setWindowTitle(mTitle);
          currentEvolutionGraph->show();
          m_charts.append(currentEvolutionGraph);
        }
      else currentEvolutionGraph=nullptr;


      FShortFileName = getShortFileName(FFileName,0);
      GlobalBreak=false;
      //////////////////////////////////////////
      etacha->Etacha(ResultsDir, FShortFileName, DebugMode);
      //////////////////////////////////////////

      RunningFlag=0;
      statusBar()->showMessage("Ready");

      FFileNameOut = ResultsDir + FShortFileName + "_" + eta_filenames[en_log];
      QFile f(FFileNameOut);
      f.open(QIODevice::WriteOnly);
      if(!f.isOpen()){ QMessageBox::information(this, "File writing error", FFileNameOut); }
      else {
          QTextStream stream(&f);
          stream <<   ui->text_Main->toHtml();
          f.close();
        }
    }

  if(currentEvolutionGraph) currentEvolutionGraph->squeeze();
  if(ShowResult!=en_log) Results();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::Results()
{

  ui->text_Main->clear();
  ui->text_Main->setCurrentFont(font);

  const int reindex[] = {en_log, en_fMean, en_fSEff, en_fPied,
                         en_f09,en_f19,en_f29,en_f39,en_f49,en_f59};

  int index = reindex[ShowResult];

  if(ShowResult > 5  &&
     ((gZb<20 && index >= en_f29) ||
      (gZb<30 && index >= en_f39) ||
      (gZb<40 && index >= en_f49) ||
      (gZb<50 && index >= en_f59))
     )
    {
      QMessageBox::information(this, "Check output choice",
                               "Z of the projectile is less than"
                                     " e- state of requested file!");
    }
  else  {
      QString LocalName = ResultsDir + FShortFileName + "_" + eta_filenames[index];
      QFile f(LocalName);
      f.open(QIODevice::ReadOnly);
      if(!f.isOpen()){ QMessageBox::information(this, "File reading error", LocalName); }
      else {
          QTextStream stream(&f);
          QString text;
          while(!stream.atEnd())
            {
              QString byg;
              byg= stream.readLine();
              text += byg + '\n';
            }

          if(index==en_log) ui->text_Main->setHtml(text);
          else              ui->text_Main->setText(text);
          f.close();
        }
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::CM_appendText(int style, int color)
{
  appendTextEdit(ui->text_Main,style,color);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::CM_appendShell(int style, int color)
{
  appendTextEdit(ui->text_Shells,style,color);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::CM_updateStatusBar()
{
  statusBar()->showMessage(BufRtf);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::plots_clicked()
{
  Plots = ui->check_Plots->isChecked();
  //if(drawPlots();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::plotsEvolution_clicked()
{
  PlotEvolution = ui->check_EvolutionPlot->isChecked();
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::appendTextEdit(QTextEdit *qte, int style, int color)
{
  QString message;
  message = BufRtf;
  const char *colorCh[] = {"red", "blue", "navy", "olive", "gray", "green"};

  switch(style)
    {
    case  fsoBold : message.insert(0,"<b>");
      message.append("</b>");
      break;
    case  fsoItalic : message.insert(0,"<i>");
      message.append("</i>");
      break;
    case  fsoUnderline : message.insert(0,"<u>");
      message.append("</u>");
      break;
    }

  if(color > clBlack && color <= clGreen)
    {
      QString fc="<font color=\"";
      fc += colorCh[color-1];  fc += "\">";
      message.insert(0,fc);
      message.append("</font>");
    }

  fastAppend(qte, message);
  qApp->processEvents();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::setPlotGeometry(e_graph *graph)
{
  QSize scr = QGuiApplication::primaryScreen()->availableSize();
  int w=scr.width()/10;
  int h=scr.height()/10;

  int W = 5*w;
  int H = 5*h;
  QRect rect(4.8*w,0.4*h,W,H);
  graph->setGeometry(rect);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::setPlotGeometryE(e_graphEvolution *graph)
{
  QSize scr = QGuiApplication::primaryScreen()->availableSize();
  int w=scr.width()/10;
  int h=scr.height()/10;

  int W = 5*w;
  int H = 4.2*h;
  QRect rect(4.8*w,5.75*h,W,H);
  graph->setGeometry(rect);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::CM_updateGraph()
{
  if(Plots)  currentGraph->updatePlots();
  if(PlotEvolution)  currentEvolutionGraph->updatePlot();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::cmLinkConstruction()
{
  QString link = "http://lise.nscl.msu.edu/paper/etacha/constructionNotes.pdf";
  QDesktopServices::openUrl(QUrl(link));
}
