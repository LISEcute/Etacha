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
//#include "e_Constant.h"

extern double gAb,gZb,gQb,gEnergy;
extern double gAt,gZt,gThick;
extern double QM, QF, dQF, currentThick;
/*const char *shellName[nSh] = {"1s","2s","2p","3s","3p","3d","n=4"};
const int   shellFull[nSh] = {2, 2, 6, 2, 6, 10, 32};
area_point initPoint;*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

e_graphEvolution::e_graphEvolution(QWidget *parent) :
  QWidget(parent),
  m_ui(new Ui_e_graphForm)
{
  m_ui->setupUi(this);

  for(int i=0; i<Number_Proton; i++) yield[i]=0;

  EvolutionChart = new QChartView(createEvolutionChart());
  m_ui->gridLayout->addWidget(EvolutionChart, 0, 0);

  m_charts << EvolutionChart;
  yield[int(gQb+0.5)]=100;


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
e_graphEvolution::~e_graphEvolution(){    delete m_ui;}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void e_graphEvolution::updatePlot()
{
  const QList<QAbstractSeries *> list =  EvolutionChart->chart()->series();

  int index=0;
  double v;

  for (const QAbstractSeries *item : list)
    {
      v = PR[index+1];
      yield[index] += v;
      ((QSplineSeries*)item)->append(currentThick,v);
      if(yield[index]<1) ((QSplineSeries*)item)->hide();
      else                  ((QSplineSeries*)item)->show();

      index++;
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void e_graphEvolution::squeeze()
{
  QList<QAbstractSeries *> list =  EvolutionChart->chart()->series();

  for (const QAbstractSeries *item : list)
    {
    QLineSeries *series = (QSplineSeries*)item;
    if(!series->isVisible())
      {
        EvolutionChart->chart()->removeSeries(series);
        delete series;
      }
    }

    list.squeeze();
    EvolutionChart->chart()->legend()->setAlignment(Qt::AlignRight);
    EvolutionChart->chart()->legend()->show();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QChart *e_graphEvolution::createEvolutionChart() const
{
  QChart *chart = new QChart();
  chart->setTitle("Charge state evolution (%)");
  chart->legend()->hide();

 int Z= gZb+0.5;
 int Q= gQb+0.5;

  for(int k=0; k<=qMin(Z,59); k++)
    {
//      QSplineSeries *series = new QSplineSeries(chart);
      QLineSeries *series = new QLineSeries(chart);
      series->setName("q="+QString::number(Z-k)+"<sup>+</sup>");

      int ran1= rand()%256;
      int ran2= rand()%256;
      int ran3= rand()%256;

      QPen penError(QColor(ran1,ran2,ran3));
      penError.setWidth(2);

      series->setPen(penError);
      series->append(0, k==Z-Q ? 100 : 0);
      chart->addSeries(series);

  //    series->attachAxis(axisX);
   //   series->attachAxis(axisY);
    }



  chart->createDefaultAxes();
  chart->axes(Qt::Horizontal).constFirst()->setRange(0,gThick);
  chart->axes(Qt::Horizontal).constFirst()->setTitleText("Thickness &nbsp; <font size=-1>[mg/cm<sup>2</sup>]</font>");

  chart->axes(Qt::Vertical).constFirst()->setRange( 0, 100);

  QValueAxis *axisY = qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).constFirst());
  Q_ASSERT(axisY);
  axisY->setLabelFormat("%.1f  ");

  return chart;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
