#ifndef E_GRAPH_H
#define E_GRAPH_H

#include <QtWidgets/QWidget>
#include <QtCharts/QChartGlobal>
#include <QtCharts/QAreaSeries>
#include "e_Constant.h"

QT_BEGIN_NAMESPACE
//class QComboBox;
//class QCheckBox;
class Ui_e_graphForm;
QT_END_NAMESPACE

//QT_CHARTS_BEGIN_NAMESPACE
class QChartView;
class QChart;
//QT_CHARTS_END_NAMESPACE

typedef QPair<QPointF, QString> Data;
typedef QList<Data> DataList;
typedef QList<DataList> DataTable;

#define nSh 7

struct area_point
{
    double x;
    double a[nSh];
//    area_point(double xi, float k[nSh])
//    {x=x;}
};


typedef QList<area_point> DataArea;

//-------------------------------------------------------------
//QT_CHARTS_USE_NAMESPACE

class e_graph: public QWidget
{
    Q_OBJECT
public:
    explicit e_graph(QWidget *parent = 0);
    ~e_graph();

    void updatePlots();


private:

    QChartView *MeanChart;
    QChartView *ChargeBarChart;
    QChartView *ShellChart;
    QChartView *CsChart;

    QChart *createMeanChart() const;
    QChart *createChargeBarChart(int valueCount, int qlow) const;
    QChart *createShellChart() const;
    QChart *createCsChart(int n) const;

 //   QAreaSeries *areaSeries[nSh];
 //   QLineSeries *lineSeries[nSh];

    int findFirstBar(double q);

private:
    int m_bars;
    int m_barsCS;

    QList<QChartView *> m_charts;
    DataArea m_area;

    Ui_e_graphForm *m_ui;
};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
class e_graphEvolution: public QWidget
{
    Q_OBJECT
public:
    explicit e_graphEvolution(QWidget *parent = 0);
    ~e_graphEvolution();

    void updatePlot();
    void squeeze();


private:

    QChartView *EvolutionChart;
    QChart *createEvolutionChart() const;

    double yield[Number_Proton];

 //   QAreaSeries *areaSeries[nSh];
 //   QLineSeries *lineSeries[nSh];

//    int findFirstBar(double q);

private:
//    int m_bars;
//    int m_barsCS;

    QList<QChartView *> m_charts;
  //  DataArea m_area;

    Ui_e_graphForm *m_ui;
};

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

#endif
