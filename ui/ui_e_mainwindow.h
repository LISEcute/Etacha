/********************************************************************************
** Form generated from reading UI file 'e_mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.15.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_E_MAINWINDOW_H
#define UI_E_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QIcon>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "w_Stuff/w_Label_clickable.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionOpen;
    QAction *actionSave;
    QAction *actionSave_As;
    QAction *actionPrint;
    QAction *actionExit;
    QAction *actionDOS_original_version_4;
    QAction *actionAbout;
    QAction *actionPrint_Setup;
    QAction *actionRun_ETACHA;
    QAction *actionCross_sections;
    QAction *actionETACHA_DOS_version;
    QAction *action2009;
    QAction *action2015;
    QAction *actionETACHA4_in_LISE;
    QAction *action1996;
    QWidget *centralWidget;
    QVBoxLayout *verticalLayout_10;
    QHBoxLayout *horizontalLayout_6;
    QVBoxLayout *verticalLayout_6;
    QFrame *frame2;
    QHBoxLayout *horizontalLayout_3;
    QPushButton *pB_Open;
    QPushButton *pB_Save;
    QPushButton *pB_Print;
    QPushButton *SB_Run;
    QSpacerItem *horizontalSpacer;
    QPushButton *pB_DOC_version;
    QPushButton *pB_About;
    QSpacerItem *horizontalSpacer_2;
    QGroupBox *gB_Projectile;
    QGridLayout *gridLayout_7;
    QFrame *frame3;
    QGridLayout *gridLayout_6;
    QLabel *label44;
    QLabel *label6;
    QLabel *label2;
    QLabel *label7;
    QLineEdit *edit_BeamA;
    QLineEdit *edit_BeamEl;
    QLineEdit *edit_BeamZ;
    QLineEdit *edit_BeamQ;
    QFrame *frame4;
    QGridLayout *gridLayout_5;
    QLabel *label_E_stat;
    QLabel *label_Stopping_stat;
    QLabel *label_Ei_stat;
    QLineEdit *edit_Energy;
    QLabel *label_Sinit;
    QLabel *label_Ef_stat;
    QLabel *label_EnergyFinal;
    QLabel *label_Sfinal;
    QCheckBox *check_EnergyLoss;
    QGroupBox *gB_1;
    QGridLayout *gridLayout_8;
    QLabel *label32;
    QLabel *label_Zorbital;
    QLabel *label33;
    QLabel *label_Qorbital;
    QGroupBox *gB_Target;
    QHBoxLayout *horizontalLayout_5;
    QFrame *frame5;
    QGridLayout *gridLayout_9;
    QLabel *label3;
    QLabel *labelElement;
    QLabel *label5;
    QLineEdit *edit_TargA;
    QLineEdit *edit_TargEl;
    QLineEdit *edit_TargZ;
    QFrame *frame6;
    QGridLayout *gridLayout_10;
    QLabel *label_TargetAtoms;
    QLabel *label39;
    QLabel *label_23;
    QLineEdit *edit_Thick;
    QLabel *label17;
    QLabel *label_24;
    QLineEdit *edit_Density;
    QLabel *label11;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *verticalLayout_5;
    ClickableLabel *label_Construction;
    QSpacerItem *horizontalSpacer_3;
    QPushButton *pb_Run2;
    QGroupBox *gB_ReactionChar;
    QGridLayout *gridLayout_4;
    QLabel *label1;
    QLabel *label_Kp1;
    QLabel *label_Kp_at_n;
    QLabel *label_Kpn;
    QLabel *label23;
    QLabel *label_vps;
    QLabel *label_Vp;
    QLabel *label_au;
    QLabel *label20;
    QVBoxLayout *verticalLayout;
    QGroupBox *gB_Version;
    QGridLayout *gridLayout_2;
    QRadioButton *rB_23;
    QLabel *label19;
    QRadioButton *rB_3;
    QLabel *label24;
    QLabel *label35;
    QRadioButton *rB_34;
    QLabel *label27;
    QLabel *label34;
    QRadioButton *rB_4;
    QLabel *label26;
    QLabel *label30;
    QRadioButton *rB_45;
    QLabel *label25;
    QLabel *label31;
    QHBoxLayout *horizontalLayout;
    QGroupBox *gB_Ionization;
    QVBoxLayout *verticalLayout_2;
    QRadioButton *rb_IC;
    QRadioButton *rb_IP;
    QGroupBox *gB_Excitation;
    QVBoxLayout *verticalLayout_3;
    QRadioButton *rb_ES;
    QRadioButton *rb_EP;
    QGroupBox *gB_ibin;
    QVBoxLayout *verticalLayout_4;
    QRadioButton *rB_corEmpirical;
    QRadioButton *rB_corNo;
    QRadioButton *rB_corBinding;
    QVBoxLayout *verticalLayout_7;
    QGroupBox *gB_Integration;
    QVBoxLayout *verticalLayout_8;
    QRadioButton *rb_ODE;
    QRadioButton *rb_RKF45;
    QRadioButton *rb_EM;
    QGroupBox *gB_NumUncertBox;
    QGridLayout *gridLayout_3;
    QLabel *label_NumberOfSteps;
    QLineEdit *editUncertainAbs;
    QLabel *label8;
    QLineEdit *editUncertainRel;
    QLabel *label_ThickStep;
    QLineEdit *edit_MinStep;
    QLabel *label12;
    QLabel *label_dZtarget;
    QLineEdit *edit_MaxStep;
    QLabel *label13;
    QGroupBox *gB_Show;
    QVBoxLayout *verticalLayout_9;
    QComboBox *cb_ShowResults;
    QCheckBox *check_CS;
    QCheckBox *check_Plots;
    QCheckBox *check_EvolutionPlot;
    QCheckBox *check_Debug;
    QHBoxLayout *horizontalLayout_4;
    QTextEdit *text_Main;
    QTextEdit *text_Shells;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuHelp;
    QMenu *menuTools;
    QStatusBar *statusBar;
    QButtonGroup *RadioG_Integration;
    QButtonGroup *RadioG_Ionization;
    QButtonGroup *RadioG_ibin;
    QButtonGroup *RadioG_Version;
    QButtonGroup *RadioG_Excitation;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1100, 800);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(MainWindow->sizePolicy().hasHeightForWidth());
        MainWindow->setSizePolicy(sizePolicy);
        MainWindow->setMinimumSize(QSize(1100, 700));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/etacha.ico"), QSize(), QIcon::Normal, QIcon::Off);
        MainWindow->setWindowIcon(icon);
        actionOpen = new QAction(MainWindow);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8("Icons/open.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionOpen->setIcon(icon1);
        actionSave = new QAction(MainWindow);
        actionSave->setObjectName(QString::fromUtf8("actionSave"));
        QIcon icon2;
        icon2.addFile(QString::fromUtf8("Icons/save.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSave->setIcon(icon2);
        actionSave_As = new QAction(MainWindow);
        actionSave_As->setObjectName(QString::fromUtf8("actionSave_As"));
        actionPrint = new QAction(MainWindow);
        actionPrint->setObjectName(QString::fromUtf8("actionPrint"));
        QIcon icon3;
        icon3.addFile(QString::fromUtf8("Icons/print.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionPrint->setIcon(icon3);
        actionExit = new QAction(MainWindow);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionDOS_original_version_4 = new QAction(MainWindow);
        actionDOS_original_version_4->setObjectName(QString::fromUtf8("actionDOS_original_version_4"));
        actionAbout = new QAction(MainWindow);
        actionAbout->setObjectName(QString::fromUtf8("actionAbout"));
        QIcon icon4;
        icon4.addFile(QString::fromUtf8("Icons/lisepp_small.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        actionAbout->setIcon(icon4);
        actionPrint_Setup = new QAction(MainWindow);
        actionPrint_Setup->setObjectName(QString::fromUtf8("actionPrint_Setup"));
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/print.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionPrint_Setup->setIcon(icon5);
        actionRun_ETACHA = new QAction(MainWindow);
        actionRun_ETACHA->setObjectName(QString::fromUtf8("actionRun_ETACHA"));
        actionCross_sections = new QAction(MainWindow);
        actionCross_sections->setObjectName(QString::fromUtf8("actionCross_sections"));
        actionETACHA_DOS_version = new QAction(MainWindow);
        actionETACHA_DOS_version->setObjectName(QString::fromUtf8("actionETACHA_DOS_version"));
        action2009 = new QAction(MainWindow);
        action2009->setObjectName(QString::fromUtf8("action2009"));
        action2015 = new QAction(MainWindow);
        action2015->setObjectName(QString::fromUtf8("action2015"));
        actionETACHA4_in_LISE = new QAction(MainWindow);
        actionETACHA4_in_LISE->setObjectName(QString::fromUtf8("actionETACHA4_in_LISE"));
        action1996 = new QAction(MainWindow);
        action1996->setObjectName(QString::fromUtf8("action1996"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        verticalLayout_10 = new QVBoxLayout(centralWidget);
        verticalLayout_10->setSpacing(6);
        verticalLayout_10->setContentsMargins(11, 11, 11, 11);
        verticalLayout_10->setObjectName(QString::fromUtf8("verticalLayout_10"));
        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(10);
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        verticalLayout_6 = new QVBoxLayout();
        verticalLayout_6->setSpacing(9);
        verticalLayout_6->setObjectName(QString::fromUtf8("verticalLayout_6"));
        frame2 = new QFrame(centralWidget);
        frame2->setObjectName(QString::fromUtf8("frame2"));
        frame2->setMinimumSize(QSize(0, 50));
        frame2->setMaximumSize(QSize(16777215, 50));
        frame2->setStyleSheet(QString::fromUtf8("background-color: rgb(192, 192, 192);"));
        horizontalLayout_3 = new QHBoxLayout(frame2);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        pB_Open = new QPushButton(frame2);
        pB_Open->setObjectName(QString::fromUtf8("pB_Open"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(pB_Open->sizePolicy().hasHeightForWidth());
        pB_Open->setSizePolicy(sizePolicy1);
        pB_Open->setMinimumSize(QSize(0, 31));
        pB_Open->setMaximumSize(QSize(16777215, 31));
        pB_Open->setStyleSheet(QString::fromUtf8("background-color: rgba(230, 230, 230, 230);"));
        QIcon icon6;
        icon6.addFile(QString::fromUtf8(":/open.png"), QSize(), QIcon::Normal, QIcon::Off);
        pB_Open->setIcon(icon6);
        pB_Open->setIconSize(QSize(32, 32));

        horizontalLayout_3->addWidget(pB_Open);

        pB_Save = new QPushButton(frame2);
        pB_Save->setObjectName(QString::fromUtf8("pB_Save"));
        sizePolicy1.setHeightForWidth(pB_Save->sizePolicy().hasHeightForWidth());
        pB_Save->setSizePolicy(sizePolicy1);
        pB_Save->setMinimumSize(QSize(0, 31));
        pB_Save->setMaximumSize(QSize(16777215, 31));
        pB_Save->setStyleSheet(QString::fromUtf8("background-color: rgba(230, 230, 230, 230);"));
        QIcon icon7;
        icon7.addFile(QString::fromUtf8(":/save.png"), QSize(), QIcon::Normal, QIcon::Off);
        pB_Save->setIcon(icon7);
        pB_Save->setIconSize(QSize(32, 32));

        horizontalLayout_3->addWidget(pB_Save);

        pB_Print = new QPushButton(frame2);
        pB_Print->setObjectName(QString::fromUtf8("pB_Print"));
        sizePolicy1.setHeightForWidth(pB_Print->sizePolicy().hasHeightForWidth());
        pB_Print->setSizePolicy(sizePolicy1);
        pB_Print->setMinimumSize(QSize(0, 31));
        pB_Print->setMaximumSize(QSize(16777215, 31));
        pB_Print->setStyleSheet(QString::fromUtf8("background-color: rgba(230, 230, 230, 230);"));
        pB_Print->setIcon(icon5);
        pB_Print->setIconSize(QSize(32, 32));

        horizontalLayout_3->addWidget(pB_Print);

        SB_Run = new QPushButton(frame2);
        SB_Run->setObjectName(QString::fromUtf8("SB_Run"));
        sizePolicy1.setHeightForWidth(SB_Run->sizePolicy().hasHeightForWidth());
        SB_Run->setSizePolicy(sizePolicy1);
        SB_Run->setMinimumSize(QSize(0, 31));
        SB_Run->setMaximumSize(QSize(16777215, 31));
        SB_Run->setStyleSheet(QString::fromUtf8("background-color: rgba(230, 230, 230, 230);"));
        QIcon icon8;
        icon8.addFile(QString::fromUtf8(":/trans.png"), QSize(), QIcon::Normal, QIcon::Off);
        SB_Run->setIcon(icon8);
        SB_Run->setIconSize(QSize(32, 32));

        horizontalLayout_3->addWidget(SB_Run);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Preferred, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer);

        pB_DOC_version = new QPushButton(frame2);
        pB_DOC_version->setObjectName(QString::fromUtf8("pB_DOC_version"));
        sizePolicy1.setHeightForWidth(pB_DOC_version->sizePolicy().hasHeightForWidth());
        pB_DOC_version->setSizePolicy(sizePolicy1);
        pB_DOC_version->setMinimumSize(QSize(90, 31));
        pB_DOC_version->setMaximumSize(QSize(16777215, 31));
        pB_DOC_version->setStyleSheet(QString::fromUtf8("background-color: rgba(230, 230, 230, 230);"));

        horizontalLayout_3->addWidget(pB_DOC_version);

        pB_About = new QPushButton(frame2);
        pB_About->setObjectName(QString::fromUtf8("pB_About"));
        sizePolicy1.setHeightForWidth(pB_About->sizePolicy().hasHeightForWidth());
        pB_About->setSizePolicy(sizePolicy1);
        pB_About->setMinimumSize(QSize(90, 31));
        pB_About->setMaximumSize(QSize(16777215, 31));
        pB_About->setStyleSheet(QString::fromUtf8("background-color: rgba(230, 230, 230, 230);"));
        QIcon icon9;
        icon9.addFile(QString::fromUtf8(":/lisepp_small.bmp"), QSize(), QIcon::Normal, QIcon::Off);
        pB_About->setIcon(icon9);
        pB_About->setIconSize(QSize(32, 32));

        horizontalLayout_3->addWidget(pB_About);

        horizontalSpacer_2 = new QSpacerItem(20, 20, QSizePolicy::Preferred, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_2);

        horizontalLayout_3->setStretch(5, 1);
        horizontalLayout_3->setStretch(6, 1);

        verticalLayout_6->addWidget(frame2);

        gB_Projectile = new QGroupBox(centralWidget);
        gB_Projectile->setObjectName(QString::fromUtf8("gB_Projectile"));
        sizePolicy.setHeightForWidth(gB_Projectile->sizePolicy().hasHeightForWidth());
        gB_Projectile->setSizePolicy(sizePolicy);
        gridLayout_7 = new QGridLayout(gB_Projectile);
        gridLayout_7->setSpacing(6);
        gridLayout_7->setContentsMargins(11, 11, 11, 11);
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        gridLayout_7->setHorizontalSpacing(12);
        gridLayout_7->setContentsMargins(-1, 6, -1, -1);
        frame3 = new QFrame(gB_Projectile);
        frame3->setObjectName(QString::fromUtf8("frame3"));
        frame3->setFrameShape(QFrame::NoFrame);
        frame3->setFrameShadow(QFrame::Sunken);
        gridLayout_6 = new QGridLayout(frame3);
        gridLayout_6->setSpacing(6);
        gridLayout_6->setContentsMargins(11, 11, 11, 11);
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        label44 = new QLabel(frame3);
        label44->setObjectName(QString::fromUtf8("label44"));
        label44->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(label44, 0, 0, 1, 1);

        label6 = new QLabel(frame3);
        label6->setObjectName(QString::fromUtf8("label6"));
        label6->setMinimumSize(QSize(55, 0));
        label6->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(label6, 0, 1, 1, 1);

        label2 = new QLabel(frame3);
        label2->setObjectName(QString::fromUtf8("label2"));
        label2->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(label2, 0, 2, 1, 1);

        label7 = new QLabel(frame3);
        label7->setObjectName(QString::fromUtf8("label7"));
        label7->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(label7, 0, 3, 1, 1);

        edit_BeamA = new QLineEdit(frame3);
        edit_BeamA->setObjectName(QString::fromUtf8("edit_BeamA"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(edit_BeamA->sizePolicy().hasHeightForWidth());
        edit_BeamA->setSizePolicy(sizePolicy2);
        edit_BeamA->setMinimumSize(QSize(40, 20));
        edit_BeamA->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(edit_BeamA, 1, 0, 1, 1);

        edit_BeamEl = new QLineEdit(frame3);
        edit_BeamEl->setObjectName(QString::fromUtf8("edit_BeamEl"));
        sizePolicy2.setHeightForWidth(edit_BeamEl->sizePolicy().hasHeightForWidth());
        edit_BeamEl->setSizePolicy(sizePolicy2);
        edit_BeamEl->setMinimumSize(QSize(45, 20));
        edit_BeamEl->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(edit_BeamEl, 1, 1, 1, 1);

        edit_BeamZ = new QLineEdit(frame3);
        edit_BeamZ->setObjectName(QString::fromUtf8("edit_BeamZ"));
        sizePolicy2.setHeightForWidth(edit_BeamZ->sizePolicy().hasHeightForWidth());
        edit_BeamZ->setSizePolicy(sizePolicy2);
        edit_BeamZ->setMinimumSize(QSize(40, 20));
        edit_BeamZ->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(edit_BeamZ, 1, 2, 1, 1);

        edit_BeamQ = new QLineEdit(frame3);
        edit_BeamQ->setObjectName(QString::fromUtf8("edit_BeamQ"));
        sizePolicy2.setHeightForWidth(edit_BeamQ->sizePolicy().hasHeightForWidth());
        edit_BeamQ->setSizePolicy(sizePolicy2);
        edit_BeamQ->setMinimumSize(QSize(40, 20));
        edit_BeamQ->setAlignment(Qt::AlignCenter);

        gridLayout_6->addWidget(edit_BeamQ, 1, 3, 1, 1);


        gridLayout_7->addWidget(frame3, 0, 0, 1, 1);

        frame4 = new QFrame(gB_Projectile);
        frame4->setObjectName(QString::fromUtf8("frame4"));
        frame4->setFrameShape(QFrame::Panel);
        frame4->setFrameShadow(QFrame::Sunken);
        gridLayout_5 = new QGridLayout(frame4);
        gridLayout_5->setSpacing(6);
        gridLayout_5->setContentsMargins(11, 11, 11, 11);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        gridLayout_5->setContentsMargins(-1, -1, -1, 7);
        label_E_stat = new QLabel(frame4);
        label_E_stat->setObjectName(QString::fromUtf8("label_E_stat"));
        sizePolicy2.setHeightForWidth(label_E_stat->sizePolicy().hasHeightForWidth());
        label_E_stat->setSizePolicy(sizePolicy2);
        label_E_stat->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(label_E_stat, 0, 1, 1, 1);

        label_Stopping_stat = new QLabel(frame4);
        label_Stopping_stat->setObjectName(QString::fromUtf8("label_Stopping_stat"));
        sizePolicy2.setHeightForWidth(label_Stopping_stat->sizePolicy().hasHeightForWidth());
        label_Stopping_stat->setSizePolicy(sizePolicy2);
        label_Stopping_stat->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(label_Stopping_stat, 0, 2, 1, 1);

        label_Ei_stat = new QLabel(frame4);
        label_Ei_stat->setObjectName(QString::fromUtf8("label_Ei_stat"));
        label_Ei_stat->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(label_Ei_stat, 1, 0, 1, 1);

        edit_Energy = new QLineEdit(frame4);
        edit_Energy->setObjectName(QString::fromUtf8("edit_Energy"));
        sizePolicy2.setHeightForWidth(edit_Energy->sizePolicy().hasHeightForWidth());
        edit_Energy->setSizePolicy(sizePolicy2);
        edit_Energy->setMinimumSize(QSize(70, 20));
        edit_Energy->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(edit_Energy, 1, 1, 1, 1);

        label_Sinit = new QLabel(frame4);
        label_Sinit->setObjectName(QString::fromUtf8("label_Sinit"));
        sizePolicy2.setHeightForWidth(label_Sinit->sizePolicy().hasHeightForWidth());
        label_Sinit->setSizePolicy(sizePolicy2);
        label_Sinit->setMinimumSize(QSize(70, 20));
        label_Sinit->setFrameShape(QFrame::Panel);
        label_Sinit->setFrameShadow(QFrame::Sunken);
        label_Sinit->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(label_Sinit, 1, 2, 1, 1);

        label_Ef_stat = new QLabel(frame4);
        label_Ef_stat->setObjectName(QString::fromUtf8("label_Ef_stat"));
        label_Ef_stat->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(label_Ef_stat, 2, 0, 1, 1);

        label_EnergyFinal = new QLabel(frame4);
        label_EnergyFinal->setObjectName(QString::fromUtf8("label_EnergyFinal"));
        sizePolicy2.setHeightForWidth(label_EnergyFinal->sizePolicy().hasHeightForWidth());
        label_EnergyFinal->setSizePolicy(sizePolicy2);
        label_EnergyFinal->setMinimumSize(QSize(70, 20));
        label_EnergyFinal->setFrameShape(QFrame::Panel);
        label_EnergyFinal->setFrameShadow(QFrame::Sunken);
        label_EnergyFinal->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(label_EnergyFinal, 2, 1, 1, 1);

        label_Sfinal = new QLabel(frame4);
        label_Sfinal->setObjectName(QString::fromUtf8("label_Sfinal"));
        sizePolicy2.setHeightForWidth(label_Sfinal->sizePolicy().hasHeightForWidth());
        label_Sfinal->setSizePolicy(sizePolicy2);
        label_Sfinal->setMinimumSize(QSize(70, 20));
        label_Sfinal->setFrameShape(QFrame::Panel);
        label_Sfinal->setFrameShadow(QFrame::Sunken);
        label_Sfinal->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(label_Sfinal, 2, 2, 1, 1);

        check_EnergyLoss = new QCheckBox(frame4);
        check_EnergyLoss->setObjectName(QString::fromUtf8("check_EnergyLoss"));
        sizePolicy2.setHeightForWidth(check_EnergyLoss->sizePolicy().hasHeightForWidth());
        check_EnergyLoss->setSizePolicy(sizePolicy2);

        gridLayout_5->addWidget(check_EnergyLoss, 3, 0, 1, 3, Qt::AlignHCenter);

        gridLayout_5->setColumnStretch(0, 1);
        gridLayout_5->setColumnStretch(1, 2);
        gridLayout_5->setColumnStretch(2, 2);

        gridLayout_7->addWidget(frame4, 0, 1, 2, 1);

        gB_1 = new QGroupBox(gB_Projectile);
        gB_1->setObjectName(QString::fromUtf8("gB_1"));
        gridLayout_8 = new QGridLayout(gB_1);
        gridLayout_8->setSpacing(6);
        gridLayout_8->setContentsMargins(11, 11, 11, 11);
        gridLayout_8->setObjectName(QString::fromUtf8("gridLayout_8"));
        label32 = new QLabel(gB_1);
        label32->setObjectName(QString::fromUtf8("label32"));
        label32->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_8->addWidget(label32, 0, 0, 1, 1);

        label_Zorbital = new QLabel(gB_1);
        label_Zorbital->setObjectName(QString::fromUtf8("label_Zorbital"));
        sizePolicy2.setHeightForWidth(label_Zorbital->sizePolicy().hasHeightForWidth());
        label_Zorbital->setSizePolicy(sizePolicy2);
        label_Zorbital->setMinimumSize(QSize(60, 20));
        label_Zorbital->setFrameShape(QFrame::Panel);
        label_Zorbital->setFrameShadow(QFrame::Sunken);
        label_Zorbital->setAlignment(Qt::AlignCenter);

        gridLayout_8->addWidget(label_Zorbital, 0, 1, 1, 1);

        label33 = new QLabel(gB_1);
        label33->setObjectName(QString::fromUtf8("label33"));
        label33->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_8->addWidget(label33, 1, 0, 1, 1);

        label_Qorbital = new QLabel(gB_1);
        label_Qorbital->setObjectName(QString::fromUtf8("label_Qorbital"));
        sizePolicy2.setHeightForWidth(label_Qorbital->sizePolicy().hasHeightForWidth());
        label_Qorbital->setSizePolicy(sizePolicy2);
        label_Qorbital->setMinimumSize(QSize(60, 20));
        label_Qorbital->setFrameShape(QFrame::Panel);
        label_Qorbital->setFrameShadow(QFrame::Sunken);
        label_Qorbital->setAlignment(Qt::AlignCenter);

        gridLayout_8->addWidget(label_Qorbital, 1, 1, 1, 1);


        gridLayout_7->addWidget(gB_1, 1, 0, 1, 1);

        gridLayout_7->setColumnStretch(0, 1);
        gridLayout_7->setColumnStretch(1, 1);

        verticalLayout_6->addWidget(gB_Projectile);

        gB_Target = new QGroupBox(centralWidget);
        gB_Target->setObjectName(QString::fromUtf8("gB_Target"));
        horizontalLayout_5 = new QHBoxLayout(gB_Target);
        horizontalLayout_5->setSpacing(12);
        horizontalLayout_5->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(-1, 5, -1, -1);
        frame5 = new QFrame(gB_Target);
        frame5->setObjectName(QString::fromUtf8("frame5"));
        frame5->setFrameShape(QFrame::NoFrame);
        frame5->setFrameShadow(QFrame::Sunken);
        gridLayout_9 = new QGridLayout(frame5);
        gridLayout_9->setSpacing(6);
        gridLayout_9->setContentsMargins(11, 11, 11, 11);
        gridLayout_9->setObjectName(QString::fromUtf8("gridLayout_9"));
        gridLayout_9->setContentsMargins(9, 7, -1, -1);
        label3 = new QLabel(frame5);
        label3->setObjectName(QString::fromUtf8("label3"));
        sizePolicy2.setHeightForWidth(label3->sizePolicy().hasHeightForWidth());
        label3->setSizePolicy(sizePolicy2);
        label3->setAlignment(Qt::AlignCenter);

        gridLayout_9->addWidget(label3, 0, 0, 1, 1);

        labelElement = new QLabel(frame5);
        labelElement->setObjectName(QString::fromUtf8("labelElement"));
        sizePolicy2.setHeightForWidth(labelElement->sizePolicy().hasHeightForWidth());
        labelElement->setSizePolicy(sizePolicy2);
        labelElement->setMinimumSize(QSize(55, 0));
        labelElement->setAlignment(Qt::AlignCenter);

        gridLayout_9->addWidget(labelElement, 0, 1, 1, 1);

        label5 = new QLabel(frame5);
        label5->setObjectName(QString::fromUtf8("label5"));
        sizePolicy2.setHeightForWidth(label5->sizePolicy().hasHeightForWidth());
        label5->setSizePolicy(sizePolicy2);
        label5->setAlignment(Qt::AlignCenter);

        gridLayout_9->addWidget(label5, 0, 2, 1, 1);

        edit_TargA = new QLineEdit(frame5);
        edit_TargA->setObjectName(QString::fromUtf8("edit_TargA"));
        sizePolicy2.setHeightForWidth(edit_TargA->sizePolicy().hasHeightForWidth());
        edit_TargA->setSizePolicy(sizePolicy2);
        edit_TargA->setMinimumSize(QSize(50, 20));
        edit_TargA->setAlignment(Qt::AlignCenter);

        gridLayout_9->addWidget(edit_TargA, 1, 0, 1, 1);

        edit_TargEl = new QLineEdit(frame5);
        edit_TargEl->setObjectName(QString::fromUtf8("edit_TargEl"));
        sizePolicy2.setHeightForWidth(edit_TargEl->sizePolicy().hasHeightForWidth());
        edit_TargEl->setSizePolicy(sizePolicy2);
        edit_TargEl->setMinimumSize(QSize(55, 20));
        edit_TargEl->setAlignment(Qt::AlignCenter);

        gridLayout_9->addWidget(edit_TargEl, 1, 1, 1, 1);

        edit_TargZ = new QLineEdit(frame5);
        edit_TargZ->setObjectName(QString::fromUtf8("edit_TargZ"));
        sizePolicy2.setHeightForWidth(edit_TargZ->sizePolicy().hasHeightForWidth());
        edit_TargZ->setSizePolicy(sizePolicy2);
        edit_TargZ->setMinimumSize(QSize(40, 20));
        edit_TargZ->setAlignment(Qt::AlignCenter);

        gridLayout_9->addWidget(edit_TargZ, 1, 2, 1, 1);


        horizontalLayout_5->addWidget(frame5);

        frame6 = new QFrame(gB_Target);
        frame6->setObjectName(QString::fromUtf8("frame6"));
        frame6->setFrameShape(QFrame::Panel);
        frame6->setFrameShadow(QFrame::Sunken);
        gridLayout_10 = new QGridLayout(frame6);
        gridLayout_10->setSpacing(6);
        gridLayout_10->setContentsMargins(11, 11, 11, 11);
        gridLayout_10->setObjectName(QString::fromUtf8("gridLayout_10"));
        gridLayout_10->setContentsMargins(-1, -1, -1, 7);
        label_TargetAtoms = new QLabel(frame6);
        label_TargetAtoms->setObjectName(QString::fromUtf8("label_TargetAtoms"));
        label_TargetAtoms->setMinimumSize(QSize(80, 20));
        label_TargetAtoms->setFrameShape(QFrame::Panel);
        label_TargetAtoms->setFrameShadow(QFrame::Sunken);
        label_TargetAtoms->setAlignment(Qt::AlignCenter);

        gridLayout_10->addWidget(label_TargetAtoms, 0, 1, 1, 1);

        label39 = new QLabel(frame6);
        label39->setObjectName(QString::fromUtf8("label39"));

        gridLayout_10->addWidget(label39, 0, 2, 1, 1);

        label_23 = new QLabel(frame6);
        label_23->setObjectName(QString::fromUtf8("label_23"));
        label_23->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_10->addWidget(label_23, 1, 0, 1, 1);

        edit_Thick = new QLineEdit(frame6);
        edit_Thick->setObjectName(QString::fromUtf8("edit_Thick"));
        edit_Thick->setMinimumSize(QSize(0, 20));
        edit_Thick->setAlignment(Qt::AlignCenter);

        gridLayout_10->addWidget(edit_Thick, 1, 1, 1, 1);

        label17 = new QLabel(frame6);
        label17->setObjectName(QString::fromUtf8("label17"));

        gridLayout_10->addWidget(label17, 1, 2, 1, 1);

        label_24 = new QLabel(frame6);
        label_24->setObjectName(QString::fromUtf8("label_24"));
        label_24->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_10->addWidget(label_24, 2, 0, 1, 1);

        edit_Density = new QLineEdit(frame6);
        edit_Density->setObjectName(QString::fromUtf8("edit_Density"));
        edit_Density->setMinimumSize(QSize(0, 20));
        edit_Density->setAlignment(Qt::AlignCenter);

        gridLayout_10->addWidget(edit_Density, 2, 1, 1, 1);

        label11 = new QLabel(frame6);
        label11->setObjectName(QString::fromUtf8("label11"));

        gridLayout_10->addWidget(label11, 2, 2, 1, 1);

        gridLayout_10->setColumnStretch(0, 2);
        gridLayout_10->setColumnStretch(1, 1);
        gridLayout_10->setColumnStretch(2, 2);

        horizontalLayout_5->addWidget(frame6);

        horizontalLayout_5->setStretch(0, 1);
        horizontalLayout_5->setStretch(1, 1);

        verticalLayout_6->addWidget(gB_Target);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(10);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        verticalLayout_5 = new QVBoxLayout();
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        verticalLayout_5->setContentsMargins(-1, -1, -1, 15);
        label_Construction = new ClickableLabel(centralWidget);
        label_Construction->setObjectName(QString::fromUtf8("label_Construction"));
        QFont font;
        font.setFamily(QString::fromUtf8("MS Shell Dlg 2"));
        font.setPointSize(9);
        font.setItalic(true);
        label_Construction->setFont(font);
        label_Construction->setStyleSheet(QString::fromUtf8("color: rgb(225, 152, 163);"));
        label_Construction->setAlignment(Qt::AlignCenter);

        verticalLayout_5->addWidget(label_Construction);

        horizontalSpacer_3 = new QSpacerItem(250, 1, QSizePolicy::Preferred, QSizePolicy::Minimum);

        verticalLayout_5->addItem(horizontalSpacer_3);

        pb_Run2 = new QPushButton(centralWidget);
        pb_Run2->setObjectName(QString::fromUtf8("pb_Run2"));
        sizePolicy2.setHeightForWidth(pb_Run2->sizePolicy().hasHeightForWidth());
        pb_Run2->setSizePolicy(sizePolicy2);
        pb_Run2->setMinimumSize(QSize(45, 25));
        pb_Run2->setMaximumSize(QSize(45, 16777215));
        pb_Run2->setIcon(icon8);
        pb_Run2->setIconSize(QSize(32, 32));

        verticalLayout_5->addWidget(pb_Run2, 0, Qt::AlignHCenter);


        horizontalLayout_2->addLayout(verticalLayout_5);

        gB_ReactionChar = new QGroupBox(centralWidget);
        gB_ReactionChar->setObjectName(QString::fromUtf8("gB_ReactionChar"));
        gridLayout_4 = new QGridLayout(gB_ReactionChar);
        gridLayout_4->setSpacing(6);
        gridLayout_4->setContentsMargins(11, 11, 11, 11);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        label1 = new QLabel(gB_ReactionChar);
        label1->setObjectName(QString::fromUtf8("label1"));
        label1->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(label1, 0, 1, 1, 1);

        label_Kp1 = new QLabel(gB_ReactionChar);
        label_Kp1->setObjectName(QString::fromUtf8("label_Kp1"));
        label_Kp1->setMinimumSize(QSize(80, 0));
        label_Kp1->setFrameShape(QFrame::Panel);
        label_Kp1->setFrameShadow(QFrame::Sunken);
        label_Kp1->setAlignment(Qt::AlignCenter);

        gridLayout_4->addWidget(label_Kp1, 0, 2, 1, 1);

        label_Kp_at_n = new QLabel(gB_ReactionChar);
        label_Kp_at_n->setObjectName(QString::fromUtf8("label_Kp_at_n"));
        label_Kp_at_n->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(label_Kp_at_n, 1, 1, 1, 1);

        label_Kpn = new QLabel(gB_ReactionChar);
        label_Kpn->setObjectName(QString::fromUtf8("label_Kpn"));
        label_Kpn->setMinimumSize(QSize(80, 0));
        label_Kpn->setFrameShape(QFrame::Panel);
        label_Kpn->setFrameShadow(QFrame::Sunken);
        label_Kpn->setAlignment(Qt::AlignCenter);

        gridLayout_4->addWidget(label_Kpn, 1, 2, 1, 1);

        label23 = new QLabel(gB_ReactionChar);
        label23->setObjectName(QString::fromUtf8("label23"));

        gridLayout_4->addWidget(label23, 2, 0, 1, 1);

        label_vps = new QLabel(gB_ReactionChar);
        label_vps->setObjectName(QString::fromUtf8("label_vps"));
        label_vps->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_4->addWidget(label_vps, 2, 1, 1, 1);

        label_Vp = new QLabel(gB_ReactionChar);
        label_Vp->setObjectName(QString::fromUtf8("label_Vp"));
        label_Vp->setMinimumSize(QSize(80, 0));
        label_Vp->setFrameShape(QFrame::Panel);
        label_Vp->setFrameShadow(QFrame::Sunken);
        label_Vp->setAlignment(Qt::AlignCenter);

        gridLayout_4->addWidget(label_Vp, 2, 2, 1, 1);

        label_au = new QLabel(gB_ReactionChar);
        label_au->setObjectName(QString::fromUtf8("label_au"));

        gridLayout_4->addWidget(label_au, 2, 3, 1, 1);

        label20 = new QLabel(gB_ReactionChar);
        label20->setObjectName(QString::fromUtf8("label20"));

        gridLayout_4->addWidget(label20, 0, 0, 2, 1);


        horizontalLayout_2->addWidget(gB_ReactionChar);

        horizontalLayout_2->setStretch(0, 1);
        horizontalLayout_2->setStretch(1, 1);

        verticalLayout_6->addLayout(horizontalLayout_2);

        verticalLayout_6->setStretch(1, 1);
        verticalLayout_6->setStretch(2, 1);
        verticalLayout_6->setStretch(3, 2);

        horizontalLayout_6->addLayout(verticalLayout_6);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(10);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        gB_Version = new QGroupBox(centralWidget);
        gB_Version->setObjectName(QString::fromUtf8("gB_Version"));
        gridLayout_2 = new QGridLayout(gB_Version);
        gridLayout_2->setSpacing(6);
        gridLayout_2->setContentsMargins(11, 11, 11, 11);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        rB_23 = new QRadioButton(gB_Version);
        RadioG_Version = new QButtonGroup(MainWindow);
        RadioG_Version->setObjectName(QString::fromUtf8("RadioG_Version"));
        RadioG_Version->addButton(rB_23);
        rB_23->setObjectName(QString::fromUtf8("rB_23"));
        QFont font1;
        font1.setFamily(QString::fromUtf8("MS Shell Dlg 2"));
        font1.setPointSize(9);
        font1.setBold(true);
        rB_23->setFont(font1);
        rB_23->setChecked(true);

        gridLayout_2->addWidget(rB_23, 0, 0, 1, 1);

        label19 = new QLabel(gB_Version);
        label19->setObjectName(QString::fromUtf8("label19"));

        gridLayout_2->addWidget(label19, 0, 1, 1, 2);

        rB_3 = new QRadioButton(gB_Version);
        RadioG_Version->addButton(rB_3);
        rB_3->setObjectName(QString::fromUtf8("rB_3"));
        rB_3->setFont(font1);

        gridLayout_2->addWidget(rB_3, 1, 0, 1, 1);

        label24 = new QLabel(gB_Version);
        label24->setObjectName(QString::fromUtf8("label24"));

        gridLayout_2->addWidget(label24, 1, 1, 1, 1);

        label35 = new QLabel(gB_Version);
        label35->setObjectName(QString::fromUtf8("label35"));
        label35->setFont(font);
        label35->setAlignment(Qt::AlignCenter);

        gridLayout_2->addWidget(label35, 1, 2, 1, 1);

        rB_34 = new QRadioButton(gB_Version);
        RadioG_Version->addButton(rB_34);
        rB_34->setObjectName(QString::fromUtf8("rB_34"));
        rB_34->setFont(font1);

        gridLayout_2->addWidget(rB_34, 2, 0, 1, 1);

        label27 = new QLabel(gB_Version);
        label27->setObjectName(QString::fromUtf8("label27"));

        gridLayout_2->addWidget(label27, 2, 1, 1, 1);

        label34 = new QLabel(gB_Version);
        label34->setObjectName(QString::fromUtf8("label34"));
        label34->setFont(font);
        label34->setAlignment(Qt::AlignCenter);

        gridLayout_2->addWidget(label34, 2, 2, 1, 1);

        rB_4 = new QRadioButton(gB_Version);
        RadioG_Version->addButton(rB_4);
        rB_4->setObjectName(QString::fromUtf8("rB_4"));
        rB_4->setFont(font1);

        gridLayout_2->addWidget(rB_4, 3, 0, 1, 1);

        label26 = new QLabel(gB_Version);
        label26->setObjectName(QString::fromUtf8("label26"));

        gridLayout_2->addWidget(label26, 3, 1, 1, 1);

        label30 = new QLabel(gB_Version);
        label30->setObjectName(QString::fromUtf8("label30"));
        label30->setFont(font);
        label30->setAlignment(Qt::AlignCenter);

        gridLayout_2->addWidget(label30, 3, 2, 1, 1);

        rB_45 = new QRadioButton(gB_Version);
        RadioG_Version->addButton(rB_45);
        rB_45->setObjectName(QString::fromUtf8("rB_45"));
        rB_45->setFont(font1);

        gridLayout_2->addWidget(rB_45, 4, 0, 1, 1);

        label25 = new QLabel(gB_Version);
        label25->setObjectName(QString::fromUtf8("label25"));

        gridLayout_2->addWidget(label25, 4, 1, 1, 1);

        label31 = new QLabel(gB_Version);
        label31->setObjectName(QString::fromUtf8("label31"));
        label31->setFont(font);
        label31->setAlignment(Qt::AlignCenter);

        gridLayout_2->addWidget(label31, 4, 2, 1, 1);

        gridLayout_2->setColumnStretch(0, 1);
        gridLayout_2->setColumnStretch(1, 2);
        gridLayout_2->setColumnStretch(2, 1);

        verticalLayout->addWidget(gB_Version);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(9);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        gB_Ionization = new QGroupBox(centralWidget);
        gB_Ionization->setObjectName(QString::fromUtf8("gB_Ionization"));
        verticalLayout_2 = new QVBoxLayout(gB_Ionization);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(14, -1, -1, -1);
        rb_IC = new QRadioButton(gB_Ionization);
        RadioG_Ionization = new QButtonGroup(MainWindow);
        RadioG_Ionization->setObjectName(QString::fromUtf8("RadioG_Ionization"));
        RadioG_Ionization->addButton(rb_IC);
        rb_IC->setObjectName(QString::fromUtf8("rb_IC"));
        rb_IC->setChecked(true);

        verticalLayout_2->addWidget(rb_IC);

        rb_IP = new QRadioButton(gB_Ionization);
        RadioG_Ionization->addButton(rb_IP);
        rb_IP->setObjectName(QString::fromUtf8("rb_IP"));

        verticalLayout_2->addWidget(rb_IP);


        horizontalLayout->addWidget(gB_Ionization);

        gB_Excitation = new QGroupBox(centralWidget);
        gB_Excitation->setObjectName(QString::fromUtf8("gB_Excitation"));
        verticalLayout_3 = new QVBoxLayout(gB_Excitation);
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setContentsMargins(11, 11, 11, 11);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        verticalLayout_3->setContentsMargins(14, -1, -1, -1);
        rb_ES = new QRadioButton(gB_Excitation);
        RadioG_Excitation = new QButtonGroup(MainWindow);
        RadioG_Excitation->setObjectName(QString::fromUtf8("RadioG_Excitation"));
        RadioG_Excitation->addButton(rb_ES);
        rb_ES->setObjectName(QString::fromUtf8("rb_ES"));
        rb_ES->setChecked(true);

        verticalLayout_3->addWidget(rb_ES);

        rb_EP = new QRadioButton(gB_Excitation);
        RadioG_Excitation->addButton(rb_EP);
        rb_EP->setObjectName(QString::fromUtf8("rb_EP"));

        verticalLayout_3->addWidget(rb_EP);


        horizontalLayout->addWidget(gB_Excitation);

        horizontalLayout->setStretch(0, 1);
        horizontalLayout->setStretch(1, 1);

        verticalLayout->addLayout(horizontalLayout);

        gB_ibin = new QGroupBox(centralWidget);
        gB_ibin->setObjectName(QString::fromUtf8("gB_ibin"));
        verticalLayout_4 = new QVBoxLayout(gB_ibin);
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setContentsMargins(11, 11, 11, 11);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        rB_corEmpirical = new QRadioButton(gB_ibin);
        RadioG_ibin = new QButtonGroup(MainWindow);
        RadioG_ibin->setObjectName(QString::fromUtf8("RadioG_ibin"));
        RadioG_ibin->addButton(rB_corEmpirical);
        rB_corEmpirical->setObjectName(QString::fromUtf8("rB_corEmpirical"));
        rB_corEmpirical->setChecked(true);

        verticalLayout_4->addWidget(rB_corEmpirical);

        rB_corNo = new QRadioButton(gB_ibin);
        RadioG_ibin->addButton(rB_corNo);
        rB_corNo->setObjectName(QString::fromUtf8("rB_corNo"));

        verticalLayout_4->addWidget(rB_corNo);

        rB_corBinding = new QRadioButton(gB_ibin);
        RadioG_ibin->addButton(rB_corBinding);
        rB_corBinding->setObjectName(QString::fromUtf8("rB_corBinding"));

        verticalLayout_4->addWidget(rB_corBinding);


        verticalLayout->addWidget(gB_ibin);

        verticalLayout->setStretch(0, 3);
        verticalLayout->setStretch(1, 2);
        verticalLayout->setStretch(2, 2);

        horizontalLayout_6->addLayout(verticalLayout);

        verticalLayout_7 = new QVBoxLayout();
        verticalLayout_7->setSpacing(6);
        verticalLayout_7->setObjectName(QString::fromUtf8("verticalLayout_7"));
        gB_Integration = new QGroupBox(centralWidget);
        gB_Integration->setObjectName(QString::fromUtf8("gB_Integration"));
        verticalLayout_8 = new QVBoxLayout(gB_Integration);
        verticalLayout_8->setSpacing(6);
        verticalLayout_8->setContentsMargins(11, 11, 11, 11);
        verticalLayout_8->setObjectName(QString::fromUtf8("verticalLayout_8"));
        rb_ODE = new QRadioButton(gB_Integration);
        RadioG_Integration = new QButtonGroup(MainWindow);
        RadioG_Integration->setObjectName(QString::fromUtf8("RadioG_Integration"));
        RadioG_Integration->addButton(rb_ODE);
        rb_ODE->setObjectName(QString::fromUtf8("rb_ODE"));
        rb_ODE->setChecked(true);

        verticalLayout_8->addWidget(rb_ODE);

        rb_RKF45 = new QRadioButton(gB_Integration);
        RadioG_Integration->addButton(rb_RKF45);
        rb_RKF45->setObjectName(QString::fromUtf8("rb_RKF45"));

        verticalLayout_8->addWidget(rb_RKF45);

        rb_EM = new QRadioButton(gB_Integration);
        rb_EM->setObjectName(QString::fromUtf8("rb_EM"));

        verticalLayout_8->addWidget(rb_EM);


        verticalLayout_7->addWidget(gB_Integration);

        gB_NumUncertBox = new QGroupBox(centralWidget);
        gB_NumUncertBox->setObjectName(QString::fromUtf8("gB_NumUncertBox"));
        gridLayout_3 = new QGridLayout(gB_NumUncertBox);
        gridLayout_3->setSpacing(6);
        gridLayout_3->setContentsMargins(11, 11, 11, 11);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        label_NumberOfSteps = new QLabel(gB_NumUncertBox);
        label_NumberOfSteps->setObjectName(QString::fromUtf8("label_NumberOfSteps"));
        label_NumberOfSteps->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_NumberOfSteps, 0, 0, 1, 1);

        editUncertainAbs = new QLineEdit(gB_NumUncertBox);
        editUncertainAbs->setObjectName(QString::fromUtf8("editUncertainAbs"));
        editUncertainAbs->setAlignment(Qt::AlignCenter);

        gridLayout_3->addWidget(editUncertainAbs, 0, 1, 1, 1);

        label8 = new QLabel(gB_NumUncertBox);
        label8->setObjectName(QString::fromUtf8("label8"));
        label8->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label8, 1, 0, 1, 1);

        editUncertainRel = new QLineEdit(gB_NumUncertBox);
        editUncertainRel->setObjectName(QString::fromUtf8("editUncertainRel"));
        editUncertainRel->setAlignment(Qt::AlignCenter);

        gridLayout_3->addWidget(editUncertainRel, 1, 1, 1, 1);

        label_ThickStep = new QLabel(gB_NumUncertBox);
        label_ThickStep->setObjectName(QString::fromUtf8("label_ThickStep"));
        label_ThickStep->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_ThickStep, 2, 0, 1, 1);

        edit_MinStep = new QLineEdit(gB_NumUncertBox);
        edit_MinStep->setObjectName(QString::fromUtf8("edit_MinStep"));
        edit_MinStep->setAlignment(Qt::AlignCenter);

        gridLayout_3->addWidget(edit_MinStep, 2, 1, 1, 1);

        label12 = new QLabel(gB_NumUncertBox);
        label12->setObjectName(QString::fromUtf8("label12"));

        gridLayout_3->addWidget(label12, 2, 2, 1, 1);

        label_dZtarget = new QLabel(gB_NumUncertBox);
        label_dZtarget->setObjectName(QString::fromUtf8("label_dZtarget"));
        label_dZtarget->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_dZtarget, 3, 0, 1, 1);

        edit_MaxStep = new QLineEdit(gB_NumUncertBox);
        edit_MaxStep->setObjectName(QString::fromUtf8("edit_MaxStep"));
        edit_MaxStep->setAlignment(Qt::AlignCenter);

        gridLayout_3->addWidget(edit_MaxStep, 3, 1, 1, 1);

        label13 = new QLabel(gB_NumUncertBox);
        label13->setObjectName(QString::fromUtf8("label13"));

        gridLayout_3->addWidget(label13, 3, 2, 1, 1);

        gridLayout_3->setColumnStretch(0, 4);
        gridLayout_3->setColumnStretch(1, 3);
        gridLayout_3->setColumnStretch(2, 1);

        verticalLayout_7->addWidget(gB_NumUncertBox);

        gB_Show = new QGroupBox(centralWidget);
        gB_Show->setObjectName(QString::fromUtf8("gB_Show"));
        verticalLayout_9 = new QVBoxLayout(gB_Show);
        verticalLayout_9->setSpacing(6);
        verticalLayout_9->setContentsMargins(11, 11, 11, 11);
        verticalLayout_9->setObjectName(QString::fromUtf8("verticalLayout_9"));
        verticalLayout_9->setContentsMargins(20, -1, 20, -1);
        cb_ShowResults = new QComboBox(gB_Show);
        cb_ShowResults->addItem(QString());
        cb_ShowResults->addItem(QString());
        cb_ShowResults->addItem(QString());
        cb_ShowResults->addItem(QString());
        cb_ShowResults->addItem(QString());
        cb_ShowResults->addItem(QString());
        cb_ShowResults->addItem(QString());
        cb_ShowResults->addItem(QString());
        cb_ShowResults->addItem(QString());
        cb_ShowResults->addItem(QString());
        cb_ShowResults->setObjectName(QString::fromUtf8("cb_ShowResults"));

        verticalLayout_9->addWidget(cb_ShowResults);

        check_CS = new QCheckBox(gB_Show);
        check_CS->setObjectName(QString::fromUtf8("check_CS"));

        verticalLayout_9->addWidget(check_CS);

        check_Plots = new QCheckBox(gB_Show);
        check_Plots->setObjectName(QString::fromUtf8("check_Plots"));

        verticalLayout_9->addWidget(check_Plots);

        check_EvolutionPlot = new QCheckBox(gB_Show);
        check_EvolutionPlot->setObjectName(QString::fromUtf8("check_EvolutionPlot"));

        verticalLayout_9->addWidget(check_EvolutionPlot);

        check_Debug = new QCheckBox(gB_Show);
        check_Debug->setObjectName(QString::fromUtf8("check_Debug"));

        verticalLayout_9->addWidget(check_Debug);


        verticalLayout_7->addWidget(gB_Show);


        horizontalLayout_6->addLayout(verticalLayout_7);

        horizontalLayout_6->setStretch(0, 3);
        horizontalLayout_6->setStretch(1, 2);
        horizontalLayout_6->setStretch(2, 1);

        verticalLayout_10->addLayout(horizontalLayout_6);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(9);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        text_Main = new QTextEdit(centralWidget);
        text_Main->setObjectName(QString::fromUtf8("text_Main"));
        sizePolicy.setHeightForWidth(text_Main->sizePolicy().hasHeightForWidth());
        text_Main->setSizePolicy(sizePolicy);

        horizontalLayout_4->addWidget(text_Main);

        text_Shells = new QTextEdit(centralWidget);
        text_Shells->setObjectName(QString::fromUtf8("text_Shells"));
        sizePolicy.setHeightForWidth(text_Shells->sizePolicy().hasHeightForWidth());
        text_Shells->setSizePolicy(sizePolicy);
        text_Shells->setMinimumSize(QSize(370, 0));

        horizontalLayout_4->addWidget(text_Shells);

        horizontalLayout_4->setStretch(0, 5);
        horizontalLayout_4->setStretch(1, 2);

        verticalLayout_10->addLayout(horizontalLayout_4);

        verticalLayout_10->setStretch(0, 1);
        verticalLayout_10->setStretch(1, 2);
        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1100, 22));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuHelp = new QMenu(menuBar);
        menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
        menuTools = new QMenu(menuBar);
        menuTools->setObjectName(QString::fromUtf8("menuTools"));
        MainWindow->setMenuBar(menuBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);
        QWidget::setTabOrder(pB_Save, pB_Print);
        QWidget::setTabOrder(pB_Print, SB_Run);
        QWidget::setTabOrder(SB_Run, pB_DOC_version);
        QWidget::setTabOrder(pB_DOC_version, pB_About);
        QWidget::setTabOrder(pB_About, edit_BeamA);
        QWidget::setTabOrder(edit_BeamA, edit_BeamEl);
        QWidget::setTabOrder(edit_BeamEl, edit_BeamZ);
        QWidget::setTabOrder(edit_BeamZ, edit_BeamQ);
        QWidget::setTabOrder(edit_BeamQ, edit_Energy);
        QWidget::setTabOrder(edit_Energy, check_EnergyLoss);
        QWidget::setTabOrder(check_EnergyLoss, edit_TargA);
        QWidget::setTabOrder(edit_TargA, edit_TargEl);
        QWidget::setTabOrder(edit_TargEl, edit_TargZ);
        QWidget::setTabOrder(edit_TargZ, edit_Thick);
        QWidget::setTabOrder(edit_Thick, edit_Density);
        QWidget::setTabOrder(edit_Density, rB_23);
        QWidget::setTabOrder(rB_23, rB_3);
        QWidget::setTabOrder(rB_3, rB_34);
        QWidget::setTabOrder(rB_34, rB_4);
        QWidget::setTabOrder(rB_4, rB_45);
        QWidget::setTabOrder(rB_45, rb_IC);
        QWidget::setTabOrder(rb_IC, rb_IP);
        QWidget::setTabOrder(rb_IP, rb_ES);
        QWidget::setTabOrder(rb_ES, rb_EP);
        QWidget::setTabOrder(rb_EP, rB_corEmpirical);
        QWidget::setTabOrder(rB_corEmpirical, rB_corNo);
        QWidget::setTabOrder(rB_corNo, rB_corBinding);
        QWidget::setTabOrder(rB_corBinding, rb_ODE);
        QWidget::setTabOrder(rb_ODE, rb_RKF45);
        QWidget::setTabOrder(rb_RKF45, editUncertainAbs);
        QWidget::setTabOrder(editUncertainAbs, editUncertainRel);
        QWidget::setTabOrder(editUncertainRel, edit_MinStep);
        QWidget::setTabOrder(edit_MinStep, edit_MaxStep);
        QWidget::setTabOrder(edit_MaxStep, cb_ShowResults);
        QWidget::setTabOrder(cb_ShowResults, check_CS);
        QWidget::setTabOrder(check_CS, pb_Run2);
        QWidget::setTabOrder(pb_Run2, text_Main);
        QWidget::setTabOrder(text_Main, text_Shells);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuTools->menuAction());
        menuBar->addAction(menuHelp->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addSeparator();
        menuFile->addAction(actionSave);
        menuFile->addAction(actionSave_As);
        menuFile->addSeparator();
        menuFile->addAction(actionPrint);
        menuFile->addAction(actionPrint_Setup);
        menuFile->addSeparator();
        menuFile->addAction(actionExit);
        menuHelp->addAction(action1996);
        menuHelp->addAction(action2009);
        menuHelp->addAction(action2015);
        menuHelp->addAction(actionETACHA4_in_LISE);
        menuHelp->addSeparator();
        menuHelp->addAction(actionAbout);
        menuTools->addAction(actionRun_ETACHA);
        menuTools->addAction(actionETACHA_DOS_version);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "ETACHA 4 (GUI)", nullptr));
        actionOpen->setText(QCoreApplication::translate("MainWindow", "Open...", nullptr));
        actionSave->setText(QCoreApplication::translate("MainWindow", "Save", nullptr));
        actionSave_As->setText(QCoreApplication::translate("MainWindow", "Save As...", nullptr));
        actionPrint->setText(QCoreApplication::translate("MainWindow", "Print...", nullptr));
        actionExit->setText(QCoreApplication::translate("MainWindow", "Exit", nullptr));
        actionDOS_original_version_4->setText(QCoreApplication::translate("MainWindow", "DOS (original) version 4", nullptr));
        actionAbout->setText(QCoreApplication::translate("MainWindow", "About...", nullptr));
        actionPrint_Setup->setText(QCoreApplication::translate("MainWindow", "Print Setup ...", nullptr));
        actionRun_ETACHA->setText(QCoreApplication::translate("MainWindow", "Run ETACHA", nullptr));
        actionCross_sections->setText(QCoreApplication::translate("MainWindow", "Cross sections", nullptr));
        actionETACHA_DOS_version->setText(QCoreApplication::translate("MainWindow", "ETACHA DOS (original) version 4", nullptr));
        action2009->setText(QCoreApplication::translate("MainWindow", "2009: J.P.Rozet, presentation @ S3WS", nullptr));
        action2015->setText(QCoreApplication::translate("MainWindow", "2015: E.Lamour et al., PRA 92, 042703", nullptr));
        actionETACHA4_in_LISE->setText(QCoreApplication::translate("MainWindow", "ETACHA page @ LISE", nullptr));
        action1996->setText(QCoreApplication::translate("MainWindow", "1996: J.P.Rozet et al., NIMB 107, 67", nullptr));
        pB_Open->setText(QString());
        pB_Save->setText(QString());
        pB_Print->setText(QString());
        SB_Run->setText(QString());
        pB_DOC_version->setText(QCoreApplication::translate("MainWindow", "DOS version", nullptr));
        pB_About->setText(QCoreApplication::translate("MainWindow", "  About", nullptr));
        gB_Projectile->setTitle(QCoreApplication::translate("MainWindow", "Projectile", nullptr));
        label44->setText(QCoreApplication::translate("MainWindow", "A", nullptr));
        label6->setText(QCoreApplication::translate("MainWindow", "Element", nullptr));
        label2->setText(QCoreApplication::translate("MainWindow", "Z", nullptr));
        label7->setText(QCoreApplication::translate("MainWindow", "Q", nullptr));
        edit_BeamA->setInputMask(QCoreApplication::translate("MainWindow", "999", nullptr));
        edit_BeamA->setText(QCoreApplication::translate("MainWindow", "208", nullptr));
        edit_BeamEl->setText(QCoreApplication::translate("MainWindow", "Pb", nullptr));
        edit_BeamZ->setInputMask(QCoreApplication::translate("MainWindow", "99", nullptr));
        edit_BeamZ->setText(QCoreApplication::translate("MainWindow", "82", nullptr));
        edit_BeamQ->setInputMask(QCoreApplication::translate("MainWindow", "99", nullptr));
        edit_BeamQ->setText(QCoreApplication::translate("MainWindow", "60", nullptr));
        label_E_stat->setText(QCoreApplication::translate("MainWindow", "Energy\n"
"(MeV/u)", nullptr));
        label_Stopping_stat->setText(QCoreApplication::translate("MainWindow", "Stopping<br>power<br><font size=-1>(MeV/mg/cm<sup>2</sup>)</font>", nullptr));
        label_Ei_stat->setText(QCoreApplication::translate("MainWindow", "Initial", nullptr));
        label_Sinit->setText(QCoreApplication::translate("MainWindow", "122.11", nullptr));
        label_Ef_stat->setText(QCoreApplication::translate("MainWindow", "Final", nullptr));
        label_EnergyFinal->setText(QCoreApplication::translate("MainWindow", "122.11", nullptr));
        label_Sfinal->setText(QCoreApplication::translate("MainWindow", "122.11", nullptr));
        check_EnergyLoss->setText(QCoreApplication::translate("MainWindow", "Use Energy Loss\n"
"Calculations", nullptr));
        gB_1->setTitle(QCoreApplication::translate("MainWindow", "Last orbital of", nullptr));
        label32->setText(QCoreApplication::translate("MainWindow", "Neutral atom =", nullptr));
        label_Zorbital->setText(QCoreApplication::translate("MainWindow", "122.11", nullptr));
        label33->setText(QCoreApplication::translate("MainWindow", "Ion in ground state =", nullptr));
        label_Qorbital->setText(QCoreApplication::translate("MainWindow", "122.11", nullptr));
        gB_Target->setTitle(QCoreApplication::translate("MainWindow", "Target", nullptr));
        label3->setText(QCoreApplication::translate("MainWindow", "A", nullptr));
        labelElement->setText(QCoreApplication::translate("MainWindow", "Element", nullptr));
        label5->setText(QCoreApplication::translate("MainWindow", "Z", nullptr));
        edit_TargA->setText(QCoreApplication::translate("MainWindow", "207", nullptr));
        edit_TargEl->setText(QCoreApplication::translate("MainWindow", "Pb", nullptr));
        edit_TargZ->setText(QCoreApplication::translate("MainWindow", "82", nullptr));
        label_TargetAtoms->setText(QCoreApplication::translate("MainWindow", "5.014e+19", nullptr));
        label39->setText(QCoreApplication::translate("MainWindow", "<font size=-1>atoms/cm<sup>2</sup></font>", nullptr));
        label_23->setText(QCoreApplication::translate("MainWindow", "Thickness =", nullptr));
        label17->setText(QCoreApplication::translate("MainWindow", "<font size=-1>mg/cm<sup>2</sup></font>", nullptr));
        label_24->setText(QCoreApplication::translate("MainWindow", "Density =", nullptr));
        label11->setText(QCoreApplication::translate("MainWindow", "<font size=-1>g/cm<sup>3</sup></font>", nullptr));
#if QT_CONFIG(tooltip)
        label_Construction->setToolTip(QCoreApplication::translate("MainWindow", "<html><head/><body><p><span style=\" color:#0000ff;\">link to notes</span></p></body></html>", nullptr));
#endif // QT_CONFIG(tooltip)
        label_Construction->setText(QCoreApplication::translate("MainWindow", "still under<br>construction !", nullptr));
        pb_Run2->setText(QString());
        gB_ReactionChar->setTitle(QCoreApplication::translate("MainWindow", "Reaction characteristics", nullptr));
        label1->setText(QCoreApplication::translate("MainWindow", "Kp (n=1) =", nullptr));
        label_Kp1->setText(QCoreApplication::translate("MainWindow", "1.23", nullptr));
        label_Kp_at_n->setText(QCoreApplication::translate("MainWindow", "Kp (n=3) =", nullptr));
        label_Kpn->setText(QCoreApplication::translate("MainWindow", "1.23", nullptr));
        label23->setText(QCoreApplication::translate("MainWindow", "projectile velocity", nullptr));
        label_vps->setText(QCoreApplication::translate("MainWindow", "V<sub>p</sub> =", nullptr));
        label_Vp->setText(QCoreApplication::translate("MainWindow", "1.23", nullptr));
        label_au->setText(QCoreApplication::translate("MainWindow", "au", nullptr));
        label20->setText(QCoreApplication::translate("MainWindow", "pertubation<br>parameter", nullptr));
        gB_Version->setTitle(QCoreApplication::translate("MainWindow", "Version", nullptr));
        rB_23->setText(QCoreApplication::translate("MainWindow", "v.23", nullptr));
        label19->setText(QCoreApplication::translate("MainWindow", " Y(1s,2s,2p),Y(3s),Y(3p),Y(3d)", nullptr));
        rB_3->setText(QCoreApplication::translate("MainWindow", "v.3", nullptr));
        label24->setText(QCoreApplication::translate("MainWindow", " + Y(12,3)", nullptr));
        label35->setText(QCoreApplication::translate("MainWindow", "fast, for high E", nullptr));
        rB_34->setText(QCoreApplication::translate("MainWindow", "v.34", nullptr));
        label27->setText(QCoreApplication::translate("MainWindow", " + Y(4)", nullptr));
        label34->setText(QCoreApplication::translate("MainWindow", "do not use", nullptr));
        rB_4->setText(QCoreApplication::translate("MainWindow", "v.4", nullptr));
        label26->setText(QCoreApplication::translate("MainWindow", " + Y(123, 4)", nullptr));
        label30->setText(QCoreApplication::translate("MainWindow", "default", nullptr));
        rB_45->setText(QCoreApplication::translate("MainWindow", "v.45", nullptr));
        label25->setText(QCoreApplication::translate("MainWindow", " + Y(5)", nullptr));
        label31->setText(QCoreApplication::translate("MainWindow", "beta", nullptr));
        gB_Ionization->setTitle(QCoreApplication::translate("MainWindow", "IONIZATION model", nullptr));
        rb_IC->setText(QCoreApplication::translate("MainWindow", "CDW-EIS\n"
" (default)", nullptr));
        rb_IP->setText(QCoreApplication::translate("MainWindow", "PWBA (fast)", nullptr));
        gB_Excitation->setTitle(QCoreApplication::translate("MainWindow", "EXCITATION model", nullptr));
        rb_ES->setText(QCoreApplication::translate("MainWindow", "Symmetric-Eikonal\n"
" (default)", nullptr));
        rb_EP->setText(QCoreApplication::translate("MainWindow", "PWBA (fast)", nullptr));
        gB_ibin->setTitle(QCoreApplication::translate("MainWindow", "Corrections for PWBA (parameter \"ibin\")", nullptr));
        rB_corEmpirical->setText(QCoreApplication::translate("MainWindow", "0: empirical saturation correction (default)", nullptr));
        rB_corNo->setText(QCoreApplication::translate("MainWindow", "1: binding correction included (not recommended)", nullptr));
        rB_corBinding->setText(QCoreApplication::translate("MainWindow", "2: no empirical correction and no binding correction", nullptr));
        gB_Integration->setTitle(QCoreApplication::translate("MainWindow", "Integration model", nullptr));
        rb_ODE->setText(QCoreApplication::translate("MainWindow", "ODE  ISBN: 0716704617\n"
"(ordinary differential equation solver)", nullptr));
        rb_RKF45->setText(QCoreApplication::translate("MainWindow", "RKF45\n"
"(Runge-Kutta-Fehlberg ODE solver)", nullptr));
        rb_EM->setText(QCoreApplication::translate("MainWindow", "Euler's method", nullptr));
        gB_NumUncertBox->setTitle(QCoreApplication::translate("MainWindow", "Steps && Numerical uncertainties", nullptr));
        label_NumberOfSteps->setText(QCoreApplication::translate("MainWindow", "Absolute =", nullptr));
        editUncertainAbs->setText(QCoreApplication::translate("MainWindow", "1e-3", nullptr));
        label8->setText(QCoreApplication::translate("MainWindow", "Relative =", nullptr));
        editUncertainRel->setText(QCoreApplication::translate("MainWindow", "1e-3", nullptr));
        label_ThickStep->setText(QCoreApplication::translate("MainWindow", "Minimum step =", nullptr));
        label12->setText(QCoreApplication::translate("MainWindow", "<font size=-1>&mu;g/cm<sup>2</sup></font>", nullptr));
        label_dZtarget->setText(QCoreApplication::translate("MainWindow", "Maximum step =", nullptr));
        label13->setText(QCoreApplication::translate("MainWindow", "<font size=-1>&mu;g/cm<sup>2</sup></font>", nullptr));
        gB_Show->setTitle(QCoreApplication::translate("MainWindow", "Show Results", nullptr));
        cb_ShowResults->setItemText(0, QCoreApplication::translate("MainWindow", "Event Logs", nullptr));
        cb_ShowResults->setItemText(1, QCoreApplication::translate("MainWindow", "Q mean", nullptr));
        cb_ShowResults->setItemText(2, QCoreApplication::translate("MainWindow", "Cross Sections", nullptr));
        cb_ShowResults->setItemText(3, QCoreApplication::translate("MainWindow", "populations 1s-2p", nullptr));
        cb_ShowResults->setItemText(4, QCoreApplication::translate("MainWindow", "e- states : 00-09", nullptr));
        cb_ShowResults->setItemText(5, QCoreApplication::translate("MainWindow", "e- states : 10-19", nullptr));
        cb_ShowResults->setItemText(6, QCoreApplication::translate("MainWindow", "e- states : 20-29", nullptr));
        cb_ShowResults->setItemText(7, QCoreApplication::translate("MainWindow", "e- states : 30-39", nullptr));
        cb_ShowResults->setItemText(8, QCoreApplication::translate("MainWindow", "e- states : 40-49", nullptr));
        cb_ShowResults->setItemText(9, QCoreApplication::translate("MainWindow", "e- states : 50-59", nullptr));

        check_CS->setText(QCoreApplication::translate("MainWindow", "Intermediate output of\n"
"cross sections", nullptr));
        check_Plots->setText(QCoreApplication::translate("MainWindow", "Plots (General set)", nullptr));
        check_EvolutionPlot->setText(QCoreApplication::translate("MainWindow", "Charge state evolution plot", nullptr));
        check_Debug->setText(QCoreApplication::translate("MainWindow", "Debug mode", nullptr));
        menuFile->setTitle(QCoreApplication::translate("MainWindow", "File", nullptr));
        menuHelp->setTitle(QCoreApplication::translate("MainWindow", "Help", nullptr));
        menuTools->setTitle(QCoreApplication::translate("MainWindow", "Execute", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_E_MAINWINDOW_H
