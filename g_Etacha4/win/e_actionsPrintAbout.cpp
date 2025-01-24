#include "e_mainwindow.h"
#include "ui_e_mainwindow.h"
#include "e_about.h"

#include <QPrinter>
#include <QPrintDialog>
#include <QProcess>
#include <QDebug>


extern QString LISErootPATH;
QString lise_url="http://lise.nscl.msu.edu/";
QString etacha_url=lise_url + "paper/etacha/";
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_pB_Print_clicked()
{
   QPrinter printer;
   printer.setPrinterName("desired printer name");
   QPrintDialog dialog(&printer,this);
   if(dialog.exec() == QDialog::Rejected) return;
   ui->text_Main->print(&printer);
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionPrint_triggered()
{
    QPrinter printer;
    printer.setPrinterName("desired printer name");
    QPrintDialog dialog(&printer,this);
    if(dialog.exec() == QDialog::Rejected) return;
    ui->text_Main->print(&printer);

    //if(PrintDialog->exec())
    // ui->text_Main->print(FFileNameOut.toStdString().c_str());

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionExit_triggered(){this->close();}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_pB_DOC_version_clicked()
{
    QString link = LISErootPATH + "/ETACHA4_DOS";
    QProcess *process = new QProcess();
    const QStringList arguments;
    process->start(link, arguments);
    //----------------------------------------------------------
    //qDebug() << "process log >>" << process->state()  << "Id="  << process->processId() << program;
    auto res= process->error();
    if(res!=QProcess::UnknownError) qDebug() << "ProcessError" << res;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionDOS_original_version_4_triggered()
{
on_pB_DOC_version_clicked();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionAbout_triggered()
{ ( new About())->show();}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_pB_About_clicked()
{
    About *ab = new About;
    ab->show();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionPrint_Setup_triggered()
{
//    PrinterSetupDialog1->exec();  ???
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionETACHA4_in_LISE_triggered()
{
    QString link = lise_url + "etacha.html";
    QDesktopServices::openUrl(QUrl(link));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_action2009_triggered()
{
    QString link = etacha_url+ "Rozet_S3WS09.pdf";
    QDesktopServices::openUrl(QUrl(link));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_action2015_triggered()
{
    QString link = etacha_url+ "Lamour_PRA15_ETACHA.pdf";
    QDesktopServices::openUrl(QUrl(link));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_action1996_triggered()
{
    QString link = etacha_url+ "Rozet_NIMB96_ETACHA.pdf";
    QDesktopServices::openUrl(QUrl(link));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

