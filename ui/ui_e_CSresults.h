/********************************************************************************
** Form generated from reading UI file 'e_CSresults.ui'
**
** Created by: Qt User Interface Compiler version 5.15.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_E_CSRESULTS_H
#define UI_E_CSRESULTS_H

#include <QtCore/QVariant>
#include <QtGui/QIcon>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_FormResults
{
public:
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_4;
    QGroupBox *gB_g;
    QHBoxLayout *horizontalLayout_2;
    QTableWidget *tableGlobal;
    QGroupBox *gB_s;
    QHBoxLayout *horizontalLayout_3;
    QTableWidget *tableShell;
    QGroupBox *groupBox;
    QHBoxLayout *horizontalLayout;
    QGridLayout *gridLayout;
    QLabel *Label1;
    QLabel *Label_Ionization;
    QLabel *LabelVersion;
    QLabel *Label2;
    QLabel *Label_Excitation;
    QLabel *Label_Reaction;
    QSpacerItem *horizontalSpacer;
    QPushButton *pB_Results_Accept;
    QPushButton *pB_Cancel;
    QSpacerItem *horizontalSpacer_2;
    QLabel *label;

    void setupUi(QDialog *FormResults)
    {
        if (FormResults->objectName().isEmpty())
            FormResults->setObjectName(QString::fromUtf8("FormResults"));
        FormResults->resize(1005, 404);
        verticalLayout = new QVBoxLayout(FormResults);
        verticalLayout->setSpacing(12);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(15);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        gB_g = new QGroupBox(FormResults);
        gB_g->setObjectName(QString::fromUtf8("gB_g"));
        gB_g->setStyleSheet(QString::fromUtf8("background-color: rgb(170, 255, 255);"));
        horizontalLayout_2 = new QHBoxLayout(gB_g);
        horizontalLayout_2->setSpacing(9);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(-1, 12, -1, -1);
        tableGlobal = new QTableWidget(gB_g);
        tableGlobal->setObjectName(QString::fromUtf8("tableGlobal"));

        horizontalLayout_2->addWidget(tableGlobal);


        horizontalLayout_4->addWidget(gB_g);

        gB_s = new QGroupBox(FormResults);
        gB_s->setObjectName(QString::fromUtf8("gB_s"));
        gB_s->setStyleSheet(QString::fromUtf8("background-color: rgb(247, 255, 194);"));
        horizontalLayout_3 = new QHBoxLayout(gB_s);
        horizontalLayout_3->setSpacing(9);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(-1, 12, -1, -1);
        tableShell = new QTableWidget(gB_s);
        tableShell->setObjectName(QString::fromUtf8("tableShell"));

        horizontalLayout_3->addWidget(tableShell);


        horizontalLayout_4->addWidget(gB_s);

        horizontalLayout_4->setStretch(0, 3);
        horizontalLayout_4->setStretch(1, 5);

        verticalLayout->addLayout(horizontalLayout_4);

        groupBox = new QGroupBox(FormResults);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        horizontalLayout = new QHBoxLayout(groupBox);
        horizontalLayout->setSpacing(12);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(9, 9, 9, 9);
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        Label1 = new QLabel(groupBox);
        Label1->setObjectName(QString::fromUtf8("Label1"));
        Label1->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);

        gridLayout->addWidget(Label1, 0, 1, 1, 1);

        Label_Ionization = new QLabel(groupBox);
        Label_Ionization->setObjectName(QString::fromUtf8("Label_Ionization"));

        gridLayout->addWidget(Label_Ionization, 0, 2, 1, 1);

        LabelVersion = new QLabel(groupBox);
        LabelVersion->setObjectName(QString::fromUtf8("LabelVersion"));
        QFont font;
        font.setFamily(QString::fromUtf8("MS Shell Dlg 2"));
        font.setPointSize(10);
        LabelVersion->setFont(font);
        LabelVersion->setStyleSheet(QString::fromUtf8("color: rgb(15, 30, 165);"));
        LabelVersion->setAlignment(Qt::AlignCenter);

        gridLayout->addWidget(LabelVersion, 1, 0, 1, 1);

        Label2 = new QLabel(groupBox);
        Label2->setObjectName(QString::fromUtf8("Label2"));
        Label2->setLayoutDirection(Qt::LeftToRight);
        Label2->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);

        gridLayout->addWidget(Label2, 1, 1, 1, 1);

        Label_Excitation = new QLabel(groupBox);
        Label_Excitation->setObjectName(QString::fromUtf8("Label_Excitation"));

        gridLayout->addWidget(Label_Excitation, 1, 2, 1, 1);

        Label_Reaction = new QLabel(groupBox);
        Label_Reaction->setObjectName(QString::fromUtf8("Label_Reaction"));
        Label_Reaction->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);

        gridLayout->addWidget(Label_Reaction, 2, 1, 1, 2);

        gridLayout->setColumnStretch(0, 3);
        gridLayout->setColumnStretch(1, 1);
        gridLayout->setColumnStretch(2, 1);

        horizontalLayout->addLayout(gridLayout);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Preferred, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        pB_Results_Accept = new QPushButton(groupBox);
        pB_Results_Accept->setObjectName(QString::fromUtf8("pB_Results_Accept"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(pB_Results_Accept->sizePolicy().hasHeightForWidth());
        pB_Results_Accept->setSizePolicy(sizePolicy);
        pB_Results_Accept->setMinimumSize(QSize(0, 30));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/etacha.ico"), QSize(), QIcon::Normal, QIcon::Off);
        pB_Results_Accept->setIcon(icon);

        horizontalLayout->addWidget(pB_Results_Accept);

        pB_Cancel = new QPushButton(groupBox);
        pB_Cancel->setObjectName(QString::fromUtf8("pB_Cancel"));
        sizePolicy.setHeightForWidth(pB_Cancel->sizePolicy().hasHeightForWidth());
        pB_Cancel->setSizePolicy(sizePolicy);
        pB_Cancel->setMinimumSize(QSize(0, 30));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/trans.png"), QSize(), QIcon::Normal, QIcon::Off);
        pB_Cancel->setIcon(icon1);

        horizontalLayout->addWidget(pB_Cancel);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Preferred, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_2);

        label = new QLabel(groupBox);
        label->setObjectName(QString::fromUtf8("label"));
        label->setStyleSheet(QString::fromUtf8("color: rgb(110, 110, 160);"));
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout->addWidget(label);

        horizontalLayout->setStretch(0, 5);
        horizontalLayout->setStretch(1, 2);
        horizontalLayout->setStretch(2, 2);
        horizontalLayout->setStretch(3, 2);
        horizontalLayout->setStretch(4, 1);
        horizontalLayout->setStretch(5, 1);

        verticalLayout->addWidget(groupBox);

        verticalLayout->setStretch(0, 4);
        verticalLayout->setStretch(1, 1);
        QWidget::setTabOrder(tableGlobal, tableShell);
        QWidget::setTabOrder(tableShell, pB_Results_Accept);
        QWidget::setTabOrder(pB_Results_Accept, pB_Cancel);

        retranslateUi(FormResults);

        QMetaObject::connectSlotsByName(FormResults);
    } // setupUi

    void retranslateUi(QDialog *FormResults)
    {
        FormResults->setWindowTitle(QCoreApplication::translate("FormResults", "ETACHA cross sections", nullptr));
        gB_g->setTitle(QCoreApplication::translate("FormResults", "Capture and Ionization", nullptr));
        gB_s->setTitle(QCoreApplication::translate("FormResults", "Excitation", nullptr));
        Label1->setText(QCoreApplication::translate("FormResults", "Ionization:", nullptr));
        Label_Ionization->setText(QCoreApplication::translate("FormResults", "PWBA", nullptr));
        LabelVersion->setText(QCoreApplication::translate("FormResults", "Version", nullptr));
        Label2->setText(QCoreApplication::translate("FormResults", "Excitation:", nullptr));
        Label_Excitation->setText(QCoreApplication::translate("FormResults", "PWBA", nullptr));
        Label_Reaction->setText(QCoreApplication::translate("FormResults", "Reaction", nullptr));
        pB_Results_Accept->setText(QCoreApplication::translate("FormResults", "  Accept", nullptr));
        pB_Cancel->setText(QCoreApplication::translate("FormResults", "  Continue", nullptr));
        label->setText(QCoreApplication::translate("FormResults", "<font size=-1>Table cells can<br> be editted.<br><br>All cross sections<br> in 1e-20 cm<sup>2</sup></font>", nullptr));
    } // retranslateUi

};

namespace Ui {
    class FormResults: public Ui_FormResults {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_E_CSRESULTS_H
