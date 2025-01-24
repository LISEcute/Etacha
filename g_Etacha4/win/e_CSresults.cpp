#include "e_CSresults.h"
#include "ui_e_CSresults.h"

//#pragma hdrstop
//#pragma package(smart_init)
//#pragma resource "*.dfm"

#include "e_CSresults.h"
#include "e_Constant.h" //#include "Constants.h"
#include "e_myextern.h"
#include <QtMath>
#include <QString>
#include <QColor>


#define egNC 4
#define egNR 7

#define esNC 8
#define esNR 7

extern const char* ElementName(int IZ);
extern double gAb,gZb,gEnergy;
extern double gZt;

QStringList egShell={"1s","2s","2p","3s","3p","3d","n=4","n=5,6"};
const char *Versions[] = {"23","3","34","4","45"};
extern QString numberKillExpZero(double v, int k, bool strong=true);
/*
const int I2_Maska[6][6] = {
//2s,2p,3s,3p,3d,4n
 {15,16,17,18,19,20}, //1s
 { 0,32,21,22,23,24}, //2s
 { 0, 0,25,26,27,28}, //2p
 { 0, 0, 0,33, 0,29}, //3s
 { 0, 0, 0, 0,34,30}, //3p
 { 0, 0, 0, 0, 0,31}};//3d
*/
 const int E2_Maska[7][7] = {
//2s,2p,3s,3p,3d,4n,5n
 {29,30,31,32,33,34,22}, //1s
 { 0,46,35,36,37,38,23}, //2s
 { 0, 0,39,40,41,42,24}, //2p
 { 0, 0, 0,47, 0,43,25}, //3s
 { 0, 0, 0, 0,48,44,26}, //3p
 { 0, 0, 0, 0, 0,45,27}, //3d
 { 0, 0, 0, 0, 0, 0,28}  //n4
 };

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
FormResults::FormResults(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::FormResults)
{
//--------------------------------------------------------------------------------------
    this->setWindowFlags(Qt::WindowCloseButtonHint | Qt::WindowTitleHint | Qt::Dialog);
    ui->setupUi(this);
//--------------------------------------------------------------------------------------


    const char *chIon[] = {"CDW-EIS","PWBA"};
    const char *chExc[] = {"SE","PWBA"};

    //int total = 0;

    QString Buffer = "<font size=+1><b>Version ";
    Buffer += Versions[EtachaVersion];
    Buffer += "</font></b>";
    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWwww   Labels
    ui->LabelVersion->setText(Buffer);
    ui->Label_Ionization->setText(chIon[IonizationModel]);
    ui->Label_Excitation->setText(chExc[ExcitationModel]);

    char Cbeam[3];
    char Ctarg[3];
    strncpy(Cbeam,ElementName(gZb),2); int kb=2; if(Cbeam[1]==' ')kb=1; Cbeam[kb]=0;
    strncpy(Ctarg,ElementName(gZt),2); int kt=2; if(Ctarg[1]==' ')kt=1; Ctarg[kt]=0;

    Buffer = buff.asprintf("<sup>%.0f</sup>%s (%.1fMeV/u) + %s",gAb, Cbeam, gEnergy, Ctarg);
    ui->Label_Reaction->setText(Buffer);

    //========================================================================================
    tableG = ui->tableGlobal;

    const int egWidth[egNC] = {50,80,80,80};

    const QHeaderView::ResizeMode egMode[egNC]={QHeaderView::Stretch, QHeaderView::Stretch,
                       QHeaderView::Stretch,QHeaderView::Stretch};

    const QStringList egColName={"(sub)\nshell","MEC\n(capture)","REC\n(capture)","Ionization"};


    tableG->verticalHeader()->setMinimumSectionSize(28);
    tableG->verticalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    //tableS->verticalHeader()->setDefaultSectionSize(22);
    tableG->setRowCount(egNR);
    tableG->setColumnCount(egNC);
    tableG->setSelectionBehavior(QAbstractItemView::SelectItems);
    tableG->setSelectionMode    (QAbstractItemView::SingleSelection);

    QString stylesheetV = "::section{Background-color:rgb(190,245,245);border-radius:12px;text-align:center;}";
    tableG->verticalHeader()->setStyleSheet(stylesheetV);
    tableG->verticalHeader()->setMinimumWidth(15);
    tableG->verticalHeader()->setMinimumHeight(30);
    tableG->verticalHeader()->setDefaultAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    //tableG->verticalHeader()->set

    QString stylesheetH = "::section{Background-color:rgb(177,217,197);border-radius:12px;text-align:center;}";
    tableG->horizontalHeader()->setStyleSheet(stylesheetH);

//tableG->horizontalHeader()->set
    for(int i=0; i<egNC; i++)
       {
       tableG->horizontalHeader()->setSectionResizeMode(i, egMode[i]);
       tableG->setColumnWidth(i, egWidth[i]);
       tableG->setHorizontalHeaderLabels(egColName);
       }

    QTableWidgetItem* pCell;
    QBrush qb = QBrush(QColor(204,205,244));

    for(int i=0; i<egNR; i++)
       {
       pCell = new QTableWidgetItem;
       pCell->setBackground(qb);
       pCell->setText(egShell[i]);
       pCell->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
       pCell->setFlags( pCell->flags() & ~Qt::ItemIsEditable);

       tableG->setItem(i, 0, pCell);
       }

    //========================================================================================
    tableS = ui->tableShell;
    const int esWidth[esNC] = {50,70,70,70,70,70,70,70};
    const QStringList esShell={"From /\nTo","2s","2p","3s","3p","3d","n=4","n=5,6"};

    tableS->verticalHeader()->setMinimumSectionSize(28);
    tableS->verticalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    //tableS->verticalHeader()->setDefaultSectionSize(22);
    tableS->setRowCount(esNR);
    tableS->setColumnCount(esNC);
    tableS->setSelectionBehavior(QAbstractItemView::SelectItems);
    tableS->setSelectionMode    (QAbstractItemView::SingleSelection);

    stylesheetV = "::section{Background-color:rgb(237,247,167);border-radius:12px;text-align:center;}";
    tableS->verticalHeader()->setStyleSheet(stylesheetV);
    tableS->verticalHeader()->setMinimumWidth(15);
    tableS->verticalHeader()->setMinimumHeight(30);
    tableS->verticalHeader()->setDefaultAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    //tableS->verticalHeader()->set

    stylesheetH = "::section{Background-color:rgb(177,217,197);border-radius:12px;text-align:center;}";
    tableS->horizontalHeader()->setStyleSheet(stylesheetH);

//tableS->horizontalHeader()->set
    for(int i=0; i<esNC; i++)
       {
       tableS->horizontalHeader()->setSectionResizeMode(i, QHeaderView::Stretch);
       tableS->setColumnWidth(i, esWidth[i]);
       tableS->setHorizontalHeaderLabels(esShell);
       }


    QBrush qb2 = QBrush(QColor(204,205,244));

    for(int i=0; i<esNR; i++)
       {
       pCell = new QTableWidgetItem;
       pCell->setBackground(qb2);
       pCell->setText(egShell[i]);
       pCell->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
       pCell->setFlags( pCell->flags() & ~Qt::ItemIsEditable);

       tableS->setItem(i, 0, pCell);
       }

fillTables();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
FormResults::~FormResults(){    delete ui;}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void FormResults::on_pB_Cancel_clicked(){ reject();}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void FormResults::fillTables()
{
nGRows=egNR; nGCols=egNC;
nSRows=esNR; nSCols=esNC;

/*           R        R       C
etacha_v23,  6(3d)  5(3p)    6(3d)
etacha_v3,   6(3d)  5(3p)    6(3d)
etacha_v34,  7(n=4)  6(3d)   7(n=4)
etacha_v4,   7(n=4)  6(3d)   7(n=4)
etacha_v45   7(n=4)  7(n=4)  8(n=5,6)
*/

    if(EtachaVersion == etacha_v4 ||
       EtachaVersion == etacha_v34  )
                        {nGRows=7; nSRows=6; nSCols=7;}
else
    if(EtachaVersion == etacha_v3 ||
           EtachaVersion == etacha_v23  )
                            {nGRows=6; nSRows=5; nSCols=6;}


for(int i=nGRows; i<egNR; i++) tableG->setRowHidden(i, true);
for(int i=nSRows; i<esNR; i++) tableS->setRowHidden(i, true);
for(int i=nSCols; i<esNC; i++) tableS->setColumnHidden(i, true);

//=================================================================

QTableWidgetItem* pCell;
QBrush qb = QBrush(QColor(222,222,222));

//--------------------------------------------------------   edit 1   G

for(int row=0; row<7; row++)
   for(int col=0; col<3; col++)
        {
       pCell = new QTableWidgetItem;
       pCell->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
       pCell->setText(numberKillExpZero(Gsecs[col*7+row+1], 4));
       tableG->setItem(row, col+1, pCell);
       }

//--------------------------------------------------------   edi2 2  S

for(int row=0; row<7; row++)
   for(int col=0; col<7; col++)
        {
        pCell = new QTableWidgetItem;
        pCell->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);

        int Flag = E2_Maska[row][col];

        if(Flag > 0)  pCell->setText(numberKillExpZero(Gsecs[Flag],4));
        else          {
                      pCell->setFlags( pCell->flags() & ~Qt::ItemIsEditable);
                      pCell->setBackground(qb);
                      }

       tableS->setItem(row, col+1, pCell);
   }

ui->pB_Results_Accept->setEnabled(false);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void FormResults::on_pB_Results_Accept_clicked()
{
QTableWidgetItem* pCell;
double v;
//--------------------------------------------------------   edit G

for(int row=0; row<7; row++)
   for(int col=0; col<3; col++)
       {
       int index = col*7+row+1;
       pCell = tableG->item(row,col+1);
       v = fabs(pCell->text().toDouble());
       Gcor[index]= Gseci[index] > 0 ? v/Gseci[index] : 1;
       }

//--------------------------------------------------------   edit E

for(int row=0; row<7; row++)
   for(int col=0; col<7; col++)
        {
        int Flag = E2_Maska[row][col];
        if(Flag > 0)
              {
              pCell = tableS->item(row,col);
              v = fabs(pCell->text().toDouble());
              Gcor[Flag]= Gseci[Flag] > 0 ? v/Gseci[Flag] : 1;
              }
        }

//----------------------------
ui->pB_Results_Accept->setEnabled(false);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void FormResults::on_tableGlobal_cellChanged(int,int)
{ui->pB_Results_Accept->setEnabled(true);}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void FormResults::on_tableShell_cellChanged(int,int)
{ui->pB_Results_Accept->setEnabled(true);}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

