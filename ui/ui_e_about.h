/********************************************************************************
** Form generated from reading UI file 'e_about.ui'
**
** Created by: Qt User Interface Compiler version 5.15.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_E_ABOUT_H
#define UI_E_ABOUT_H

#include <QtCore/QVariant>
#include <QtGui/QIcon>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include "w_Stuff/w_Label_clickable.h"

QT_BEGIN_NAMESPACE

class Ui_About
{
public:
    QGridLayout *gridLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *label_13;
    QLabel *label_11;
    QLabel *label_12;
    QLabel *label_8;
    QLabel *label_18;
    QLabel *label_version;
    QLabel *label_revision;
    ClickableLabel *label_mail;
    QFrame *frame_3;
    QGridLayout *gridLayout_3;
    ClickableLabel *label_etacha;
    QLabel *label_5;
    QLabel *label_6;
    ClickableLabel *label_LISE;
    QLabel *label_7;

    void setupUi(QDialog *About)
    {
        if (About->objectName().isEmpty())
            About->setObjectName(QString::fromUtf8("About"));
        About->resize(405, 504);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(About->sizePolicy().hasHeightForWidth());
        About->setSizePolicy(sizePolicy);
        QIcon icon;
        icon.addFile(QString::fromUtf8("icons/etacha32.png"), QSize(), QIcon::Normal, QIcon::Off);
        About->setWindowIcon(icon);
        About->setStyleSheet(QString::fromUtf8("background-color: rgb(85, 85, 225);"));
        gridLayout = new QGridLayout(About);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label_13 = new QLabel(About);
        label_13->setObjectName(QString::fromUtf8("label_13"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(label_13->sizePolicy().hasHeightForWidth());
        label_13->setSizePolicy(sizePolicy1);
        label_13->setMaximumSize(QSize(70, 70));
        QFont font;
        font.setPointSize(10);
        label_13->setFont(font);
        label_13->setAutoFillBackground(false);
        label_13->setStyleSheet(QString::fromUtf8(""));
        label_13->setPixmap(QPixmap(QString::fromUtf8(":/etacha.png")));
        label_13->setScaledContents(true);

        horizontalLayout->addWidget(label_13);

        label_11 = new QLabel(About);
        label_11->setObjectName(QString::fromUtf8("label_11"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(4);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(label_11->sizePolicy().hasHeightForWidth());
        label_11->setSizePolicy(sizePolicy2);
        QFont font1;
        font1.setFamily(QString::fromUtf8("Bookman Old Style"));
        font1.setPointSize(24);
        font1.setBold(true);
        font1.setItalic(false);
        label_11->setFont(font1);
        label_11->setStyleSheet(QString::fromUtf8("color: #fff;\n"
"margin-bottom:0;"));
        label_11->setAlignment(Qt::AlignCenter);

        horizontalLayout->addWidget(label_11);

        label_12 = new QLabel(About);
        label_12->setObjectName(QString::fromUtf8("label_12"));
        sizePolicy1.setHeightForWidth(label_12->sizePolicy().hasHeightForWidth());
        label_12->setSizePolicy(sizePolicy1);
        label_12->setMaximumSize(QSize(70, 70));
        label_12->setAutoFillBackground(false);
        label_12->setStyleSheet(QString::fromUtf8(""));
        label_12->setPixmap(QPixmap(QString::fromUtf8(":/etacha.png")));
        label_12->setScaledContents(true);
        label_12->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout->addWidget(label_12);

        horizontalLayout->setStretch(0, 1);
        horizontalLayout->setStretch(2, 1);

        gridLayout->addLayout(horizontalLayout, 0, 0, 1, 1);

        label_8 = new QLabel(About);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Minimum);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(label_8->sizePolicy().hasHeightForWidth());
        label_8->setSizePolicy(sizePolicy3);
        label_8->setMinimumSize(QSize(350, 0));
        QFont font2;
        font2.setFamily(QString::fromUtf8("Bookman Old Style"));
        font2.setPointSize(13);
        font2.setBold(true);
        font2.setItalic(false);
        label_8->setFont(font2);
        label_8->setStyleSheet(QString::fromUtf8("color: rgb(254, 255, 246);"));
        label_8->setAlignment(Qt::AlignCenter);

        gridLayout->addWidget(label_8, 1, 0, 1, 1);

        label_18 = new QLabel(About);
        label_18->setObjectName(QString::fromUtf8("label_18"));
        QFont font3;
        font3.setFamily(QString::fromUtf8("MS Shell Dlg 2"));
        font3.setPointSize(9);
        label_18->setFont(font3);
        label_18->setStyleSheet(QString::fromUtf8("color: rgb(222, 226, 197); background-color: rgb(0, 90, 255);"));
        label_18->setFrameShape(QFrame::StyledPanel);
        label_18->setAlignment(Qt::AlignCenter);
        label_18->setMargin(10);
        label_18->setIndent(5);

        gridLayout->addWidget(label_18, 2, 0, 1, 1);

        label_version = new QLabel(About);
        label_version->setObjectName(QString::fromUtf8("label_version"));
        QFont font4;
        font4.setFamily(QString::fromUtf8("MS Shell Dlg 2"));
        font4.setPointSize(10);
        font4.setBold(true);
        font4.setItalic(false);
        label_version->setFont(font4);
        label_version->setStyleSheet(QString::fromUtf8("color: rgb(240, 255, 200);"));
        label_version->setAlignment(Qt::AlignCenter);

        gridLayout->addWidget(label_version, 3, 0, 1, 1);

        label_revision = new QLabel(About);
        label_revision->setObjectName(QString::fromUtf8("label_revision"));
        label_revision->setFont(font4);
        label_revision->setStyleSheet(QString::fromUtf8("color: rgb(240, 255, 200);"));
        label_revision->setAlignment(Qt::AlignCenter);

        gridLayout->addWidget(label_revision, 4, 0, 1, 1);

        label_mail = new ClickableLabel(About);
        label_mail->setObjectName(QString::fromUtf8("label_mail"));
        label_mail->setCursor(QCursor(Qt::PointingHandCursor));
        label_mail->setStyleSheet(QString::fromUtf8("padding: 5px;\n"
"color: rgb(210, 210, 210);"));
        label_mail->setAlignment(Qt::AlignCenter);
        label_mail->setWordWrap(true);

        gridLayout->addWidget(label_mail, 5, 0, 1, 1);

        frame_3 = new QFrame(About);
        frame_3->setObjectName(QString::fromUtf8("frame_3"));
        frame_3->setStyleSheet(QString::fromUtf8("background-color: rgb(0, 90, 255);"));
        frame_3->setFrameShape(QFrame::StyledPanel);
        gridLayout_3 = new QGridLayout(frame_3);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_3->setHorizontalSpacing(9);
        gridLayout_3->setVerticalSpacing(12);
        gridLayout_3->setContentsMargins(8, 10, 8, 9);
        label_etacha = new ClickableLabel(frame_3);
        label_etacha->setObjectName(QString::fromUtf8("label_etacha"));
        QFont font5;
        font5.setFamily(QString::fromUtf8("Arial"));
        font5.setPointSize(9);
        label_etacha->setFont(font5);
        label_etacha->setCursor(QCursor(Qt::PointingHandCursor));
        label_etacha->setStyleSheet(QString::fromUtf8("color: rgb(222, 226, 197);"));
        label_etacha->setTextFormat(Qt::AutoText);
        label_etacha->setOpenExternalLinks(true);

        gridLayout_3->addWidget(label_etacha, 0, 1, 1, 2);

        label_5 = new QLabel(frame_3);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        QFont font6;
        font6.setPointSize(9);
        font6.setBold(true);
        label_5->setFont(font6);
        label_5->setStyleSheet(QString::fromUtf8("color: rgb(170, 170, 255);"));

        gridLayout_3->addWidget(label_5, 1, 0, 1, 1);

        label_6 = new QLabel(frame_3);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setFont(font6);
        label_6->setStyleSheet(QString::fromUtf8("color: rgb(170, 170, 255);"));

        gridLayout_3->addWidget(label_6, 0, 0, 1, 1);

        label_LISE = new ClickableLabel(frame_3);
        label_LISE->setObjectName(QString::fromUtf8("label_LISE"));
        label_LISE->setFont(font5);
        label_LISE->setCursor(QCursor(Qt::PointingHandCursor));
        label_LISE->setStyleSheet(QString::fromUtf8("color: rgb(222, 226, 197);"));
        label_LISE->setTextFormat(Qt::AutoText);
        label_LISE->setOpenExternalLinks(true);

        gridLayout_3->addWidget(label_LISE, 1, 1, 1, 1);

        label_7 = new QLabel(frame_3);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setPixmap(QPixmap(QString::fromUtf8(":/emblem_little_plus.png")));

        gridLayout_3->addWidget(label_7, 1, 2, 1, 1);


        gridLayout->addWidget(frame_3, 6, 0, 1, 1);


        retranslateUi(About);

        QMetaObject::connectSlotsByName(About);
    } // setupUi

    void retranslateUi(QDialog *About)
    {
        About->setWindowTitle(QCoreApplication::translate("About", "About  \"ETACHA\"", nullptr));
        label_13->setText(QString());
        label_11->setText(QCoreApplication::translate("About", "ETACHA", nullptr));
        label_12->setText(QString());
        label_8->setText(QCoreApplication::translate("About", "calculating charge state distributions", nullptr));
        label_18->setText(QCoreApplication::translate("About", "<html><head/><body><p>E.Lamour, P.D.Fainstein, M.Galassi,<br/>C.Prigent, C.A.Ramirez, R.D.Rivarola,<br/>J.-P.Rozet, M.Trassinelli, and D.Vernhet</p><hr/><p>PHYSICAL REVIEW A 92, 042703 (2015)</p></body></html>", nullptr));
        label_version->setText(QCoreApplication::translate("About", "Version 4.L 2.1300 beta", nullptr));
        label_revision->setText(QCoreApplication::translate("About", "Last revision 19-JUL-2016", nullptr));
#if QT_CONFIG(tooltip)
        label_mail->setToolTip(QCoreApplication::translate("About", "send mail to:", nullptr));
#endif // QT_CONFIG(tooltip)
        label_mail->setText(QCoreApplication::translate("About", "<html><head/><body><p>This program has been converted to C++ and ported to MS Windows GUI application within the LISE<span style=\" vertical-align:super;\">++</span> framework by O.B.Tarasov (FRIB/MSU)</p><p>This program has been ported to a cross-platform <br/>application using Qt by K.V. Tarasova</p><p>The GUI-version is currently maintained by O.B. Tarasov</p></body></html>", nullptr));
#if QT_CONFIG(tooltip)
        label_etacha->setToolTip(QCoreApplication::translate("About", "<html><head/><body><p><span style=\" color:#ffffff;\">Link to ETACHA4</span></p></body></html>", nullptr));
#endif // QT_CONFIG(tooltip)
        label_etacha->setText(QCoreApplication::translate("About", "https://lise.frib.msu.edu/etacha.html", nullptr));
        label_5->setText(QCoreApplication::translate("About", "LISE<sup>++</sup>", nullptr));
        label_6->setText(QCoreApplication::translate("About", "ETACHA", nullptr));
#if QT_CONFIG(tooltip)
        label_LISE->setToolTip(QCoreApplication::translate("About", "<html><head/><body><p><span style=\" color:#ffffff;\">link to LISE++</span></p></body></html>", nullptr));
#endif // QT_CONFIG(tooltip)
        label_LISE->setText(QCoreApplication::translate("About", "https://lise.frib.msu.edu", nullptr));
        label_7->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class About: public Ui_About {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_E_ABOUT_H
