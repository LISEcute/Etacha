#ifndef E_CSRESULTS_H
#define E_CSRESULTS_H

#include <QDialog>
#include <QTableWidget>

//=========================================================
namespace Ui {
class FormResults;
}

class FormResults : public QDialog
{
    Q_OBJECT

public:
    explicit FormResults(QWidget *parent = nullptr);
    ~FormResults();

private:

    int nGRows, nGCols;
    int nSRows, nSCols;

    QString buff;
    QTableWidget *tableG;
    QTableWidget *tableS;

    void fillTables();


private slots:
    void on_pB_Results_Accept_clicked();
    void on_pB_Cancel_clicked();

    void on_tableGlobal_cellChanged(int row, int column);
    void on_tableShell_cellChanged(int row, int column);

private:
    Ui::FormResults *ui;
};

#endif // E_CSRESULTS_H
//=========================================================
