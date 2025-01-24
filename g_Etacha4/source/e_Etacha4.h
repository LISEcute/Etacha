#ifndef ETACHA4_H
#define ETACHA4_H

#include <QWidget>
#include <QObject>
#include "../win/e_Constant.h"

extern char BufRtf[1000];

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class ETACHA : public QObject
{
Q_OBJECT

public: 
    explicit ETACHA(QWidget *parentI) ;

    ~ETACHA() {};


    void init();
    int  DONAUT();
    void Tceis(double Ecurr, double &tot, int netat);
    void CSEC(double Ecurr, int iStep);
    void SEnlm(double zp8, double zt8, double E8);

    int Etacha(const QString &LinitialDir, const QString &LfileName, bool DebugMode=false);
    double CHGT(double Ecurr, double Zp, double At,double Rho,
                double as, double bs, double tc,
                double QM, FILE *f);

    void message(int idid);
//    bool test_ESC();

protected:

//    void keyPressEvent(QKeyEvent *event);


private:
    QWidget *parent;
//    bool cancelStatus=false;

//public slots:

signals:
    void appendText(int style=fsoNone,  int color=clBlack);
    void appendShell(int style=fsoNone, int color=clBlack);
    void updateStatusBar();
    void updateGraph();

};
#endif
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
