// chartthemes.pro

#include "e_graph.h"
#include "ui_e_graph.h"

#include <QtCharts/QChartView>
#include <QtCharts/QPieSeries>
#include <QtCharts/QPieSlice>
#include <QtCharts/QAbstractBarSeries>
#include <QtCharts/QPercentBarSeries>
#include <QtCharts/QStackedBarSeries>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QLineSeries>
#include <QtCharts/QSplineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QAreaSeries>
#include <QtCharts/QLegend>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QFormLayout>
//#include <QtWidgets/QComboBox>
//#include <QtWidgets/QSpinBox>
//#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QLabel>
#include <QtCore/QRandomGenerator>
#include <QtCharts/QBarCategoryAxis>
#include <QtWidgets/QApplication>
#include <QtCharts/QValueAxis>
#include <QtCharts/QLogValueAxis>

#include "e_myextern.h"
#include "e_Constant.h"

extern double gAb,gZb,gQb,gEnergy;
extern double gAt,gZt,gThick;
extern double QM, QF, dQF, currentThick;
const char *shellName[nSh] = {"1s","2s","2p","3s","3p","3d","n=4"};
const int   shellFull[nSh] = {2, 2, 6, 2, 6, 10, 32};
area_point initPoint;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

e_graph::e_graph(QWidget *parent) :
    QWidget(parent),
    m_ui(new Ui_e_graphForm)
{
    m_ui->setupUi(this);


    MeanChart = new QChartView(createMeanChart());  // mean value as function thickness
    m_ui->gridLayout->addWidget(MeanChart, 0, 0);
    m_charts << MeanChart;

    m_bars=qMin(8.,gZb);
    ChargeBarChart = new QChartView(createChargeBarChart(m_bars,findFirstBar(gQb))); // charge at current time
    m_ui->gridLayout->addWidget(ChargeBarChart, 1, 0);
    m_charts << ChargeBarChart;

    ShellChart = new QChartView(createShellChart());   //  shells
    m_ui->gridLayout->addWidget(ShellChart, 0, 1);
    m_charts << ShellChart;

    m_barsCS = 7;
    if(EtachaVersion == etacha_v3 || EtachaVersion == etacha_v23) m_barsCS=m_barsCS-1;
    CsChart = new QChartView(createCsChart(m_barsCS));   // cross sections
    m_ui->gridLayout->addWidget(CsChart, 1, 1);
    m_charts << CsChart;

    m_area.append(initPoint);

//--------------------------------------------------------------- properties start
    const auto charts = m_charts;
        if (!m_charts.isEmpty())
            {
            // int index=0;
            for (QChartView *chartView : charts) {
            //    if(index>0)chartView->chart()->setTheme(QChart::ChartThemeBlueNcs);
                chartView->setRenderHint(QPainter::Antialiasing, true);
                chartView->chart()->setAnimationOptions(QChart::NoAnimation);
            }
//    m_ui->animatedComboBox->addItem("Series Animations", QChart::SeriesAnimations);
            QPalette pal = window()->palette();

            pal.setColor(QPalette::Window, QRgb(0x018bba));
            pal.setColor(QPalette::WindowText, QRgb(0x404044));
            window()->setPalette(pal);
            }

  //--------------------------------------------------------------- properties end

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
e_graph::~e_graph(){    delete m_ui;}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void e_graph::updatePlots()
{
const QList<QAbstractSeries *> list =  MeanChart->chart()->series();
double err[4] = {0,dQF,-dQF,0};
double Q;

int index=0;
for (const QAbstractSeries *item : list)
    {

    if(index==3) Q=QM;
    else         Q=QF + err[qMin(index,2)];
    Q = qMax(0.,qMin(gZb, Q));
    ((QSplineSeries*)item)->append(currentThick,Q);
    index++;
    }

QValueAxis *axisY = qobject_cast<QValueAxis*>(MeanChart->chart()->axes(Qt::Vertical).constFirst());
double qmin = QM-dQF-0.3;
if(axisY->min() > qmin)
    if(qmin>0) axisY->setRange(qmin,gZb);

//==================================================================================
QChart *chart = ChargeBarChart->chart();
chart->removeAllSeries();
QBarCategoryAxis *axisX2 = qobject_cast<QBarCategoryAxis*>(chart->axes(Qt::Horizontal).constFirst());
axisX2->clear();

QBarSeries *series = new QBarSeries(chart);
QBarSet *set = new QBarSet("distribution");
QStringList categories;
int qlow = findFirstBar(QM);
static double procentMaxKeep=100;
double procentMax=0;

for (int i(0); i < m_bars; i++) {
    int qr = i+qlow;
    QString title = QString::number(qr) + "<sup>+</sup>";
    categories << title;
    double v = PR[int(gZb)-qr+1];
    procentMax = qMax(v,procentMax);
    *set << v;
    }

series->append(set);
chart->addSeries(series);

axisX2->append(categories);
axisY = qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).constFirst());
if(procentMax > procentMaxKeep+10 || procentMax < procentMaxKeep-10)
    {
    procentMaxKeep = procentMax;
    axisY->setRange(0,qMin(100.,procentMax+5.));
    }

series->attachAxis(axisY);
//==================================================================================


area_point T = {currentThick, {y1s, y2s, y2p, y3s, y3p, y3d, yN}};
m_area.append(T);

chart = ShellChart->chart();

chart->removeAllSeries();
//const QList<QAbstractSeries *> listShell =  chart->series();

QLineSeries *lowerSeries = nullptr;


//int nameIndex=0;
double h=0, v;
for (int shell(0); shell < nSh; shell++)
    {
    QLineSeries *upperSeries = new QLineSeries(chart);
//    QLineSeries *upperSeries =  ((QSplineSeries*)item)   QLineSeries *upperSeries = new QLineSeries(chart);

    int k=0;

    for(area_point &ap : m_area)
            {
            if (lowerSeries)
                {
                float y  = lowerSeries->at(k).y();
                v = y + ap.a[shell];
                h=qMax(h,v);
                upperSeries->append(QPointF(ap.x, v));
                }
            else {
                 upperSeries->append(QPointF(ap.x,ap.a[shell]));
                 }
            k++;
            }

    QAreaSeries *area = new QAreaSeries(upperSeries, lowerSeries);
    area->setName(shellName[shell]);
    chart->addSeries(area);
    lowerSeries = upperSeries;
    }


    chart->createDefaultAxes();
    chart->axes(Qt::Horizontal).constFirst()->setRange(0,gThick);
    chart->axes(Qt::Horizontal).constFirst()->setTitleText("Thickness &nbsp; <font size=-1>[mg/cm<sup>2</sup>]</font>");

    chart->axes(Qt::Vertical).constFirst()->setRange(0, int(h+0.9));
    axisY = qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).constFirst());
    Q_ASSERT(axisY);
    axisY->setLabelFormat("%.1f  ");

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QChart *e_graph::createMeanChart() const
{
    QChart *chart = new QChart();
    chart->setTitle("<b>&lt; q &gt;</b> &nbsp; - &nbsp; mean charge state value");
    chart->legend()->hide();

        QSplineSeries *series;
        series = new QSplineSeries(chart);
        QPen penMain(Qt::magenta);  penMain.setWidth(2);
        series->setPen(penMain);
        series->append(0,gQb);         series->setName("&lt;q&gt;");
        chart->addSeries(series);

        for(int k=0; k<2; k++)
            {
            series = new QSplineSeries(chart);
            QPen penError(Qt::blue);  penError.setWidth(1); penError.setStyle(Qt::DashLine);
            series->setPen(penError);
            series->append(0,gQb);
            chart->addSeries(series);
            }

        series = new QSplineSeries(chart);
        QPen penProb(Qt::darkGreen);  penProb.setWidth(1); penProb.setStyle(Qt::DotLine);
        series->setPen(penProb);
        series->append(0,gQb);         series->setName("prob");
        chart->addSeries(series);


    chart->createDefaultAxes();
    chart->axes(Qt::Horizontal).constFirst()->setRange(0,gThick);
    chart->axes(Qt::Horizontal).constFirst()->setTitleText("Thickness &nbsp; <font size=-1>[mg/cm<sup>2</sup>]</font>");

    chart->axes(Qt::Vertical).constFirst()->setRange( qMax(gQb-1.,0.), gZb);

    QValueAxis *axisY = qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).constFirst());
    Q_ASSERT(axisY);
    axisY->setLabelFormat("%.1f  ");

    return chart;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QChart *e_graph::createChargeBarChart(int valueCount, int qlow) const
{
    QChart *chart = new QChart();
    chart->setTitle("Charge distribution at current moment (%)");
    chart->legend()->hide();
    //chart->setAnimationOptions(QChart::AllAnimations);

    QBarSeries *series = new QBarSeries(chart);
    QBarSet *set = new QBarSet("distribution");
    QStringList categories;

    for (int i(0); i < valueCount; i++) {
        int qr = i+qlow;
        QString title = QString::number(qr) + "<sup>+</sup>";
        categories << title;
        if(qr==int(gQb)) *set << 100;
        else             *set << 0;
        }

    series->append(set);
    chart->addSeries(series);

    QBarCategoryAxis *axisX = new QBarCategoryAxis();
    axisX->append(categories);
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis();
    axisY->setRange(0,100);
    axisY->setLabelFormat("%.1f ");
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);

    return chart;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QChart *e_graph::createShellChart() const
{
  QChart *chart = new QChart();
  chart->setTitle("Shell population evolution");
  chart->legend()->setAlignment(Qt::AlignRight);

    int shellValue[nSh] = {0,0,0,0,0,0,0};

    int nE = gZb - gQb;
    int shell=0;

    while(nE>0)
        {
        if(shellValue[shell] >= shellFull[shell]) shell++;
        shellValue[shell]++;
        nE--;
        }

    initPoint.x=0;


//    int sum   = 0;

    for (int i(0); i < nSh; i++)
        {
        //sum +=
        //shellValue[i];
        initPoint.a[i]=shellValue[i];
        }
//----------------------------------------------
    QLineSeries *lowerSeries = nullptr;

    double h=0, v;
    for (int shell(0); shell < nSh; shell++)
        {
        QLineSeries *upperSeries = new QLineSeries(chart);

        if (lowerSeries)
                    {
                    float y  = lowerSeries->at(0).y();
                    v = y + shellValue[shell];
                    h=qMax(h,v);
                    upperSeries->append(QPointF(0, v));
                    }
                else {
                     upperSeries->append(QPointF(0,shellValue[shell]));
                     }

        QAreaSeries *area = new QAreaSeries(upperSeries, lowerSeries);
        area->setName(shellName[shell]);
        chart->addSeries(area);
        lowerSeries = upperSeries;
        }


        chart->createDefaultAxes();
        chart->axes(Qt::Horizontal).constFirst()->setRange(0,gThick);
        chart->axes(Qt::Horizontal).constFirst()->setTitleText("Thickness &nbsp; <font size=-1>[mg/cm<sup>2</sup>]</font>");

        chart->axes(Qt::Vertical).constFirst()->setRange(0, int(h+0.9));
        QValueAxis *axisY = qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).constFirst());
        Q_ASSERT(axisY);
        axisY->setLabelFormat("%.1f  ");

    return chart;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QChart *e_graph::createCsChart(int m_cat) const
{
    QChart *chart = new QChart();
    chart->setTitle("Initial Cross Sections (1e-20 cm<sup>2</sup>)");
    chart->legend()->setAlignment(Qt::AlignRight);


    extern QStringList egShell;
    QStringList processSets = {"MEC<i>apture<i>", "REC<i>apture</i>", "Ionization"};
    QStringList categories;

    QBarSeries *series = new QBarSeries();

// 3 sets ==> MEC, REC, Ionization
// categories  ==> nmber of points in set
// function of energy
double vmax=0, vmin=1e50;

for (int k(0); k < 3; k++)
    {
    QBarSet *set = new QBarSet(processSets[k]);
    for (int i(0); i < m_cat; i++)
        {
        if(k==0) categories << egShell[i];
        double cs = Gsecs[k*7+i+1];
        vmax = qMax(vmax,cs);
        vmin = qMin(vmin,cs);
        *set << cs;
        }
    series->append(set);
    }

chart->addSeries(series);

QBarCategoryAxis *axisX = new QBarCategoryAxis();
axisX->append(categories);
chart->addAxis(axisX, Qt::AlignBottom);
series->attachAxis(axisX);

QLogValueAxis *axisY = new QLogValueAxis();
axisY->setRange(vmin*0.5,vmax*1.1);
axisY->setLabelFormat("%.1e ");
chart->addAxis(axisY, Qt::AlignLeft);
series->attachAxis(axisY);

return chart;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int e_graph::findFirstBar(double q)
{
int qi = q + 0.5;
int qlow=qi;
int qupp=qi;

bool forward=false;
for(int i=0; i<m_bars; i++)
    {
            if(forward && qupp <= gZb ) { forward=false; qupp++; }
    else    if(qlow > 2 )               { forward=true;  qlow--; }
    }

return qlow;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

