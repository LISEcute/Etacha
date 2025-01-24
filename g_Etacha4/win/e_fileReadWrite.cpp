#include "e_mainwindow.h"
#include "ui_e_mainwindow.h"
#include "e_myextern.h"
#include "e_Constant.h"

#include <QApplication>
#include <QDir>
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <iostream>
#include <QStandardPaths>
#include <QSettings>

const char *FileNameAbsent = "/eUntitled.etacha";
const char *LISEini="/lisepp.ini";
QString *LiseppPathForData;

QString LISErootPATH;
QString MyDocCompPATH;
QString InitialDir;
QString ResultsDir;
QString localPATH;

//QString FFileNameCS;
//QString FFileNameHtml;

extern double gAb,gZb,gQb,gEnergy;
extern double gAt,gZt,gThick,gDensity,gMinStepTarget,gMaxStepTarget;
extern double gO_EnergyL1,gO_EnergyL2,gO_dEdX1,gO_dEdX2;
extern double gUncertainAbs,gUncertainRel;

extern int EtachaVersion,ibinParameter,IonizationModel,ExcitationModel,ShowItermediate,DifEqModel,UseEloss;

extern int fontsizeGlobal;
//extern int useHighDpiScaling;
extern FILE *mfopen(const QString& filename, const char* operand);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//void getDPIscaling()
//{
//      MyDocCompPATH = QStandardPaths::standardLocations(QStandardPaths::DocumentsLocation).constFirst();
//      QString FN1 = MyDocCompPATH + "/LISEcute";  FN1 += LISEini;
//      QSettings myLiseIni1(FN1,QSettings::IniFormat);
//      useHighDpiScaling  = myLiseIni1.value("font/scaling",  useHighDpiScaling).toInt();
//}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void getInitialDir(void)
{
//---------------------------------------------------------  paths begin

    MyDocCompPATH = QStandardPaths::standardLocations(QStandardPaths::DocumentsLocation).constFirst();

    QString FileCheck(LISErootPATH); FileCheck += LISEini;
    FILE *fcheck=mfopen(FileCheck, "at");
    int work_in_LISEroot_main = 0;

    if(fcheck) {                    // work in root directory
          fclose(fcheck);
          QSettings myLiseIni0(FileCheck,QSettings::IniFormat);
          work_in_LISEroot_main = myLiseIni0.value("Version/WorkInROOT",0).toInt();
          if(work_in_LISEroot_main) localPATH = LISErootPATH;
          }

    if(work_in_LISEroot_main==0)
          {
          localPATH = MyDocCompPATH;
          localPATH += "/LISEcute";
          }

    //--------------------------------------------------------- lise.ini  begin
     QString FN1=localPATH;  FN1 += LISEini;
      QSettings myLiseIni1(FN1,QSettings::IniFormat);
      fontsizeGlobal     = myLiseIni1.value("font/size",     fontsizeGlobal).toInt();
//      useHighDpiScaling  = myLiseIni1.value("font/scaling",  useHighDpiScaling).toInt();
     //--------------------------------------------------------- lise.ini  end
ResultsDir = localPATH + "/results/";
localPATH += "/files";

 QDir pathDir(localPATH);
 if(!pathDir.exists()) pathDir.mkdir(localPATH);
//---------------------------------------------------------  paths end
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::setFileName(const QString& FileName)
{
    if(FileName!=nullptr) FFileName = FileName;

    QFileInfo fI(FFileName);
    QString windowName = fI.baseName();
    if(Modified) windowName += " *";
    if(!windowName.contains(&FileNameAbsent[1]) || Modified)  this->setWindowTitle("ETACHA4 - " + windowName);
    else                                                      this->setWindowTitle("ETACHA4");
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionOpen_triggered()
{
QString file = QFileDialog::getOpenFileName(this,tr("Open ETACHA's input file"),
                    FFileName, tr("ETACHA files (*.etacha)"));

if(file.size()<=0) return;

readFile(file);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_pB_Open_clicked() {on_actionOpen_triggered();}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::readFile(const QString& file)
{
 QFile f(file);

    if(!f.open(QIODevice::ReadOnly | QIODevice::Text))
      {
        QMessageBox MB;
        MB.setWindowTitle("Warning!");
        MB.setText(tr("Could not open file.\n%1").arg(file));
        MB.setIcon(QMessageBox::Warning);
        MB.exec();
        SetPage(); //necessary if from command string
        return;
        }
    f.close();

    QSettings myOpt(file,QSettings::IniFormat);


    myOpt.beginGroup("Projectile");
         gAb  = myOpt.value("A",  238).toDouble();
         gZb  = myOpt.value("Z",  92).toDouble();
         gQb  = myOpt.value("q",  60).toDouble();
         gEnergy = myOpt.value("Energy",   30).toDouble();
    myOpt.endGroup();


    myOpt.beginGroup("Target");
       gAt = myOpt.value("A",  12).toDouble();
       gZt = myOpt.value("Z",  6).toDouble();
       gThick = myOpt.value("Thick",  1).toDouble();
       gDensity = myOpt.value("Density",  2.25).toDouble();
       gMinStepTarget = myOpt.value("MinStep",  1).toDouble();
       gMaxStepTarget = myOpt.value("MaxStep",  10).toDouble();
    myOpt.endGroup();

    myOpt.beginGroup("EnergyLoss");
        gO_EnergyL1 = myOpt.value("EnergyL1").toDouble();
        gO_EnergyL2 = myOpt.value("EnergyL2").toDouble();
        gO_dEdX1 = myOpt.value("dEdX1").toDouble();
        gO_dEdX2 = myOpt.value("dEdX2").toDouble();
    myOpt.endGroup();


    myOpt.beginGroup("Options");
        EtachaVersion= myOpt.value("Version",3).toInt();
        DifEqModel= myOpt.value("ODE_Model",0).toInt();
        ibinParameter= myOpt.value("ibin",0).toInt();
        IonizationModel= myOpt.value("Ionization",1).toInt();
        ExcitationModel= myOpt.value("Excitation",1).toInt();
        ShowItermediate= myOpt.value("ShowInterm",0).toInt();
        Plots= myOpt.value("ShowPlots",1).toInt();
        PlotEvolution= myOpt.value("EvolutionPlot",1).toInt();
        UseEloss= myOpt.value("EnergyLoss",1).toInt();

        gUncertainAbs = myOpt.value("UncertainAbs",1e-3).toDouble();
        gUncertainRel = myOpt.value("UncertainRel",1e-3).toDouble();
    myOpt.endGroup();

    setFileName(file);

    //------------------------------------------------
SetPage(true);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionSave_triggered()
{
    QFile f(FFileName);

    if(!f.open(QIODevice::ReadWrite | QIODevice::Text)){
            QMessageBox MB;
            MB.setWindowTitle("Warning!");
            QString txt = tr("Could not Save file.\n%1").arg(FFileName);
            MB.setText(txt);
            MB.setIcon(QMessageBox::Warning);
            MB.exec();
            return;
            }
        else f.close();

        Modified=false;
        setFileName(FFileName);
        ReadPage();

        QSettings myOpt(FFileName,QSettings::IniFormat);

        myOpt.beginGroup("Projectile");
            myOpt.setValue("A",gAb);
            myOpt.setValue("Z",gZb);
            myOpt.setValue("q",gQb);
            myOpt.setValue("Energy",gEnergy);
        myOpt.endGroup();


        myOpt.beginGroup("Target");
            myOpt.setValue("A",gAt);
            myOpt.setValue("Z",gZt);
            myOpt.setValue("Thick",gThick);
            myOpt.setValue("Density",gDensity);
            myOpt.setValue("MinStep",gMinStepTarget);
            myOpt.setValue("MaxStep",gMaxStepTarget);
        myOpt.endGroup();

        myOpt.beginGroup("EnergyLoss");
             myOpt.setValue("EnergyL1",gO_EnergyL1);
             myOpt.setValue("EnergyL2",gO_EnergyL2);
             myOpt.setValue("dEdX1",gO_dEdX1);
             myOpt.setValue("dEdX2",gO_dEdX2);
        myOpt.endGroup();


        myOpt.beginGroup("Options");
            myOpt.setValue("Version",EtachaVersion);
            myOpt.setValue("ODE_Model",DifEqModel);
            myOpt.setValue("ibin",ibinParameter);
            myOpt.setValue("Ionization",IonizationModel);
            myOpt.setValue("Excitation",ExcitationModel);
            myOpt.setValue("ShowInterm",ShowItermediate);           
            myOpt.setValue("ShowPlots",int(Plots));
            myOpt.setValue("EvolutionPlot",int(PlotEvolution));
            myOpt.setValue("EnergyLoss",UseEloss);

             myOpt.setValue("UncertainAbs",gUncertainAbs);
             myOpt.setValue("UncertainRel",gUncertainRel);
        myOpt.endGroup();

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_pB_Save_clicked()
{
    if(FFileName.contains(&FileNameAbsent[1])) on_actionSave_As_triggered();
    else                                       on_actionSave_triggered();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionSave_As_triggered()
{
    QString file = QFileDialog::getSaveFileName(this,
                            "ETACHA: Save As",
                            FFileName,
                           "ETACHA files (*.etacha)");

    if(file.size()<=0) return;

    FFileName=file;
    on_actionSave_triggered();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

