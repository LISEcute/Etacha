#include "e_about.h"
#include "ui_e_about.h"

#include "e_ftype.h"
#include <QUrl>
#include <QDesktopServices>


extern void CmLISEsiteLink(QString str);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
About::About(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::About)
{
    ui->setupUi(this);
    this->setWindowFlags(Qt::WindowCloseButtonHint | Qt::WindowTitleHint);


    ui->label_version->setText(ETACHAG_version);
    ui->label_revision->setText(ETACHAG_date);

    connect(ui->label_mail, SIGNAL(clicked()), this, SLOT(cmMail()));
    connect(ui->label_LISE, SIGNAL(clicked()), this, SLOT(cmSite()));
    connect(ui->label_etacha, SIGNAL(clicked()), this, SLOT(cmLink()));

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
About::~About(){    delete ui;}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::cmMail()
{
    QString ss("mailto:tarasov@frib.msu.edu?subject=Etacha ");
    ss.append(ui->label_version->text());
    QDesktopServices::openUrl(QUrl(ss));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::cmSite() {    QDesktopServices::openUrl(QUrl(ui->label_LISE->text()));}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::cmLink()
{   QString link = "https://lise.frib.msu.edu/etacha.html";
    QDesktopServices::openUrl(QUrl(link));}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::cmEtacha() {    QDesktopServices::openUrl(QUrl(ui->label_etacha->text()));}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
