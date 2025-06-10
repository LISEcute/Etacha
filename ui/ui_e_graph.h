/********************************************************************************
** Form generated from reading UI file 'e_graph.ui'
**
** Created by: Qt User Interface Compiler version 5.15.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_E_GRAPH_H
#define UI_E_GRAPH_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_e_graphForm
{
public:
    QGridLayout *gridLayout;
    QHBoxLayout *horizontalLayout;

    void setupUi(QWidget *e_graphForm)
    {
        if (e_graphForm->objectName().isEmpty())
            e_graphForm->setObjectName(QString::fromUtf8("e_graphForm"));
        e_graphForm->resize(684, 440);
        gridLayout = new QGridLayout(e_graphForm);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(7, 6, 5, 6);
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));

        gridLayout->addLayout(horizontalLayout, 0, 0, 1, 3);


        retranslateUi(e_graphForm);

        QMetaObject::connectSlotsByName(e_graphForm);
    } // setupUi

    void retranslateUi(QWidget *e_graphForm)
    {
        (void)e_graphForm;
    } // retranslateUi

};

namespace Ui {
    class e_graphForm: public Ui_e_graphForm {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_E_GRAPH_H
