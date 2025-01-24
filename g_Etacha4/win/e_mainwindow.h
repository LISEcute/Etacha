#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QDialog>
#include <QDesktopServices>
#include <QUrl>
#include <QPageSetupDialog>

#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QRadioButton>
#include <QButtonGroup>
#include <QCheckBox>
#include <QTextEdit>
#include <QComboBox>
#include "../win/e_Constant.h"
#include "e_graph.h"

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
namespace Ui {
class MainWindow;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    e_graph          *currentGraph=nullptr;
    e_graphEvolution *currentEvolutionGraph=nullptr;

protected:

    bool   FlagPermit;
    bool   Modified;
    bool   ZPedit_permit;
    bool   ZTedit_permit;
    bool   FlagReadPage;
    bool   DebugMode;
    bool   Plots;
    bool   PlotEvolution;
    char   CZ[10];

    double K_perturbation;
    double K_perturbationN;


    QString Buffer;


    QString FFileName;
    QString FFileNameCS;
    QString FShortFileName;
    QString FFileNameOut;

    QFont font;


    QButtonGroup *bg_Version;
    QButtonGroup *bg_ibin;
    QButtonGroup *bg_Ionization;
    QButtonGroup *bg_Excitation;
    QButtonGroup *bg_Integration;

    QComboBox *cb_Show;

private:
        void getInitialDir(void);
        void readFile(const QString& file);
        void CheckFileSave(void);

        void SetPage(bool fromFile=false);
        void ReadPage();
        void Results();

        void CalculateEnergy();
        void CalculateFraction(void);

        void GetInitialDir(void);
        void ShowEnergyLossCells();

      void PertrubationParameters();
      void setFileName(const QString& FileName=nullptr);
      void setPlotGeometry(e_graph*);
      void setPlotGeometryE(e_graphEvolution*);
    //  int e_etacha_init();

//------------------------------------------------
protected:
    void keyPressEvent(QKeyEvent *e) override;
    void closeEvent(QCloseEvent *event) override;


    void makeBeamZ();
    void makeBeamQ(void);
    void makeTargetZ(bool fromFile=false);
    void makeTargetThick(void);
    void appendTextEdit(QTextEdit *, int style, int color);

private slots:

    void bg_Version_clicked(int);
    void bg_ibin_clicked(int);
    void bg_Ionization_clicked(int);
    void bg_Excitation_clicked(int);
    void bg_Integration_clicked(int);

    void cb_Show_clicked();

    void on_pB_Open_clicked();
    void on_pB_Save_clicked();
    void on_pB_Print_clicked();
    void on_pB_About_clicked();
    void on_pB_DOC_version_clicked();

    void on_SB_Run_clicked();
    void on_pb_Run2_clicked();
    void on_actionRun_ETACHA_triggered();

    void on_actionOpen_triggered();
    void on_actionSave_triggered();
    void on_actionSave_As_triggered();
    void on_actionPrint_triggered();
    void on_actionPrint_Setup_triggered();
    void on_actionDOS_original_version_4_triggered();
    void on_actionAbout_triggered();
    void on_actionExit_triggered();

    void on_edit_Energy_textChanged(const QString &);
    void on_edit_BeamZ_textChanged(const QString &);
    void on_edit_BeamQ_textChanged(const QString &);
    void on_edit_BeamEl_textChanged(const QString &);
    void on_edit_BeamA_textChanged(const QString &);

    void on_check_EnergyLoss_clicked();

    void on_edit_TargZ_textChanged(const QString &);
    void on_edit_TargEl_textChanged(const QString &);
    void on_edit_TargA_textChanged(const QString &);
    void on_edit_Thick_textChanged(const QString &);
    void on_edit_Density_textChanged(const QString &arg1);

    void on_check_CS_clicked();
    void plots_clicked();
    void plotsEvolution_clicked();

    void CM_appendText(int style=fsoNone, int color=clBlack);
    void CM_appendShell(int style=fsoNone, int color=clBlack);
    void CM_updateStatusBar();
    void CM_updateGraph();
    void cmLinkConstruction();

    void on_actionETACHA4_in_LISE_triggered();
    void on_action1996_triggered();
    void on_action2009_triggered();
    void on_action2015_triggered();


private:

    QList<QWidget *> m_charts;
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H




