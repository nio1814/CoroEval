/* 
 * CoroEval
 *
 * An evaluation tool for coronary artery reconstructions.
 *
 * Copyright © 2014:
 *
 * Christoph Forman, Universität Erlangen-Nürnberg
 * christoph.forman@cs.fau.de
 *
 * Chris Schwemmer, Universität Erlangen-Nürnberg
 * chris.schwemmer@cs.fau.de
 *
 * Jens Wetzl, Universität Erlangen-Nürnberg
 * jens.wetzl@cs.fau.de
 *
 * CoroEval is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CoroEval is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CoroEval.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "PlotWidget.h"

#include "Interval.h"

#include <QPalette>
#include <QPen>

#include <qwt_plot_grid.h>
#include <qwt_plot_curve.h>
#include <qwt_painter.h>
#include <qwt_plot_canvas.h>
#include <qwt_plot_layout.h>
#include <qwt_scale_widget.h>
#include <qwt_round_scale_draw.h>
#include <qwt_scale_draw.h>
#include <qwt_plot_renderer.h>
#include <qwt_symbol.h>
#include <qwt_global.h>


PlotWidget::PlotWidget(QWidget* parent):
    QwtPlot(parent)
{

    m_curve     = 0;
    m_slicePos  = 0;
    m_lind80Pos = 0;
	m_lind50Pos = 0;
    m_lind20Pos = 0;
    m_rind80Pos = 0;
	m_rind50Pos = 0;
    m_rind20Pos = 0;
    m_lminPos   = 0;
    m_rminPos   = 0;
    m_maxPos    = 0;
	m_stenosisPos = 0;
	m_stenosisLeftPos = 0;
	m_stenosisRightPos = 0;
	m_referenceDiameter = 0;

    alignScales();

    // Insert grid
    m_grid = new QwtPlotGrid();
    m_grid->attach(this);
    m_grid->setVisible( true );

    // Axis
    setAxisScale(QwtPlot::xBottom, 0.0, 1.0);
    setAxisScale(QwtPlot::yLeft, 0.0, 1.0);

    m_axis = new QwtPlotCurve();
    m_axis->setPaintAttribute(QwtPlotCurve::ClipPolygons, true);
    m_axis->setRenderHint(QwtPlotCurve::RenderAntialiased, true);
    m_axis->attach(this);

	setWhiteOnBlack();

    std::vector<double> xAxisX, xAxisY;
    xAxisX.push_back(0);
    xAxisX.push_back(0);
    xAxisX.push_back(1);
    xAxisY.push_back(1);
    xAxisY.push_back(0);
    xAxisY.push_back(0);
    m_axis->setSamples(&xAxisX[0], &xAxisY[0], 3);
    m_axis->plot();

    canvas()->setAttribute(Qt::WA_PaintOnScreen, true);

    QwtPainter::setPolylineSplitting(false);
    replot();

	m_axisFont = axisWidget(QwtPlot::xBottom)->font();
}

PlotWidget::~PlotWidget()
{
    clear();
    delete m_axis;
}

// Current versions of Qwt have privatized the copy constructor...
static void copyQwtSymbol(const QwtSymbol *from, QwtSymbol *to)
{
    // Qwt privatized the copy constructor and operator= in
    // version 6.1.0 and also introduced some additional
    // member variables compared to 6.0.2, so neither version
    // of the code below works with both versions of Qwt,
    // hence the #if
#if QWT_VERSION < 0x060100
    *to = *from;
#else
    to->setCachePolicy(from->cachePolicy());
    to->setSize(from->size());
    to->setPinPoint(from->pinPoint());
    to->setPinPointEnabled(from->isPinPointEnabled());
    to->setBrush(from->brush());
    to->setPen(from->pen());
    to->setStyle(from->style());
    to->setPath(from->path());
    to->setPixmap(from->pixmap());
    to->setGraphic(from->graphic());
#endif
}

void PlotWidget::setBlackOnWhite()
{
	setStyleSheet("background-color: transparent;");
	setCanvasBackground(Qt::transparent);

    QPalette palette;
    palette.setColor(QPalette::WindowText, Qt::black);
	palette.setColor(QPalette::Text, Qt::black);
    palette.setColor(QPalette::Background, Qt::transparent);
	palette.setColor(QPalette::Base, Qt::black);
    axisWidget(QwtPlot::xBottom)->setPalette(palette);
    axisWidget(QwtPlot::yLeft)->setPalette(palette);

	QFont font("Helvetica", 16);
	axisWidget(QwtPlot::xBottom)->setFont(font);
	axisWidget(QwtPlot::yLeft)->setFont(font);
	
	font.setPointSize(18);
	QwtText title = axisTitle(QwtPlot::xBottom);
	title.setFont(font);
	setAxisTitle(QwtPlot::xBottom, title);

	title = axisTitle(QwtPlot::yLeft);
	title.setFont(font);
	setAxisTitle(QwtPlot::yLeft, title);

    canvas()->setPalette(palette);

	QPen pen_grid;
    pen_grid.setColor(QColor(0,0,0, 128));
    pen_grid.setStyle(Qt::DashLine);
	m_grid->setPen(pen_grid);

	m_axis->setPen(QPen(Qt::black));

	QPen pen(Qt::black);

	if (m_lind50Pos)
		m_lind50Pos->setPen(pen);
	if (m_rind50Pos)
		m_rind50Pos->setPen(pen);
	if (m_stenosisPos && m_maxPos) {
		m_stenosisPos->setPen(pen);
		m_maxPos->setPen(pen);
	}

	pen.setColor(Qt::darkYellow);
	if (m_stenosisLeftPos)
		m_stenosisLeftPos->setPen(pen);
	if (m_stenosisRightPos)
		m_stenosisRightPos->setPen(pen);

	if (m_curve) {
		QPen pen = m_curve->pen();
		pen.setWidth(3);
		m_curve->setPen(pen);
	}

	if (m_slicePos)
		m_slicePos->hide();

	if (m_maxPos && m_lminPos && m_rminPos && m_lind20Pos && m_rind20Pos && m_lind50Pos && m_rind50Pos && m_lind80Pos && m_rind80Pos) {
		QPen pen = m_curve->pen();
		pen.setWidth(4);
		m_curve->setPen(pen);			
			
		pen = m_maxPos->pen();
		pen.setWidth(3);
		m_maxPos->setPen(pen);

		pen = m_lminPos->pen();
		pen.setWidth(3);
		m_lminPos->setPen(pen);
		m_rminPos->setPen(pen);

		pen = m_lind20Pos->pen();
		pen.setWidth(3);
		m_lind20Pos->setPen(pen);
		m_rind20Pos->setPen(pen);
		m_lind80Pos->setPen(pen);
		m_rind80Pos->setPen(pen);

		pen = m_lind50Pos->pen();
		pen.setWidth(3);
		m_lind50Pos->setPen(pen);
		m_rind50Pos->setPen(pen);

        QwtSymbol* symbol = new QwtSymbol();
        copyQwtSymbol(m_lminPos->symbol(), symbol);
		symbol->setSize(10, 10);
		m_lminPos->setSymbol(symbol);
		
        symbol = new QwtSymbol();
        copyQwtSymbol(m_rminPos->symbol(), symbol);
		symbol->setSize(10, 10);
		m_rminPos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_lind20Pos->symbol(), symbol);
		symbol->setSize(10, 10);
		m_lind20Pos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_rind20Pos->symbol(), symbol);
		symbol->setSize(10, 10);
		m_rind20Pos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_lind50Pos->symbol(), symbol);
		symbol->setSize(10, 10);
		m_lind50Pos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_rind50Pos->symbol(), symbol);
		symbol->setSize(10, 10);
		m_rind50Pos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_lind80Pos->symbol(), symbol);
		symbol->setSize(10, 10);
		m_lind80Pos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_rind80Pos->symbol(), symbol);
		symbol->setSize(10, 10);
		m_rind80Pos->setSymbol(symbol);
	}

	if (m_maxPos && m_stenosisLeftPos && m_stenosisPos && m_stenosisRightPos && m_referenceDiameter) {
		QPen pen = m_curve->pen();
		pen.setWidth(4);
		m_curve->setPen(pen);

		pen = m_maxPos->pen();
		pen.setWidth(2);
		m_maxPos->setPen(pen);

		pen = m_stenosisLeftPos->pen();
		pen.setWidth(2);
		m_stenosisLeftPos->setPen(pen);
		m_stenosisRightPos->setPen(pen);

		pen = m_stenosisPos->pen();
		pen.setWidth(3);
		m_stenosisPos->setPen(pen);

		pen = m_referenceDiameter->pen();
		pen.setWidth(3);
		m_referenceDiameter->setPen(pen);
	}
}

void PlotWidget::setWhiteOnBlack(bool resetFonts)
{
	setStyleSheet("background-color: black;");
	setCanvasBackground(Qt::black);

    QPalette palette;
    palette.setColor(QPalette::WindowText, Qt::white);
    palette.setColor(QPalette::Text, Qt::white);
    palette.setColor(QPalette::Background, Qt::white);
    palette.setColor(QPalette::Base, Qt::white);
    axisWidget(QwtPlot::xBottom)->setPalette(palette);
    axisWidget(QwtPlot::yLeft)->setPalette(palette);

	if (resetFonts) {
		axisWidget(QwtPlot::xBottom)->setFont(m_axisFont);
		axisWidget(QwtPlot::yLeft)->setFont(m_axisFont);

		QwtText title = axisTitle(QwtPlot::xBottom);
		title.setFont(m_axisTitleFont);
		setAxisTitle(QwtPlot::xBottom, title);

		title = axisTitle(QwtPlot::yLeft);
		title.setFont(m_axisTitleFont);
		setAxisTitle(QwtPlot::yLeft, title);
	}

	canvas()->setPalette(palette);

	QPen pen_grid;
    pen_grid.setColor(QColor(255,255,255,50));
    pen_grid.setStyle(Qt::DashLine);
	m_grid->setPen(pen_grid);

	m_axis->setPen(QPen(Qt::white));

	QPen pen(Qt::white);

	if (m_lind50Pos)
		m_lind50Pos->setPen(pen);
	if (m_rind50Pos)
		m_rind50Pos->setPen(pen);
	if (m_stenosisPos && m_maxPos) {
		m_stenosisPos->setPen(pen);
		m_maxPos->setPen(pen);
	}

	pen.setColor(Qt::yellow);
	if (m_stenosisLeftPos)
		m_stenosisLeftPos->setPen(pen);
	if (m_stenosisRightPos)
		m_stenosisRightPos->setPen(pen);

	if (m_curve) {
		QPen pen = m_curve->pen();
		pen.setWidth(2);
		m_curve->setPen(pen);
	}

	if (m_slicePos)
		m_slicePos->show();

	if (m_maxPos && m_lminPos && m_rminPos && m_lind20Pos && m_rind20Pos && m_lind50Pos && m_rind50Pos && m_lind80Pos && m_rind80Pos) {
		QPen pen = m_maxPos->pen();
		pen.setWidth(0);
		m_maxPos->setPen(pen);

		pen = m_lminPos->pen();
		pen.setWidth(0);
		m_lminPos->setPen(pen);
		m_rminPos->setPen(pen);

		pen = m_lind20Pos->pen();
		pen.setWidth(0);
		m_lind20Pos->setPen(pen);
		m_rind20Pos->setPen(pen);
		m_lind80Pos->setPen(pen);
		m_rind80Pos->setPen(pen);

		pen = m_lind50Pos->pen();
		pen.setWidth(0);
		m_lind50Pos->setPen(pen);
		m_rind50Pos->setPen(pen);

        QwtSymbol* symbol = new QwtSymbol();
        copyQwtSymbol(m_lminPos->symbol(), symbol);
		symbol->setSize(5, 5);
		m_lminPos->setSymbol(symbol);
		
        symbol = new QwtSymbol();
        copyQwtSymbol(m_rminPos->symbol(), symbol);
		symbol->setSize(5, 5);
		m_rminPos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_lind20Pos->symbol(), symbol);
		symbol->setSize(5, 5);
		m_lind20Pos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_rind20Pos->symbol(), symbol);
		symbol->setSize(5, 5);
		m_rind20Pos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_lind50Pos->symbol(), symbol);
		symbol->setSize(5, 5);
		m_lind50Pos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_rind50Pos->symbol(), symbol);
		symbol->setSize(5, 5);
		m_rind50Pos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_lind80Pos->symbol(), symbol);
		symbol->setSize(5, 5);
		m_lind80Pos->setSymbol(symbol);

        symbol = new QwtSymbol();
        copyQwtSymbol(m_rind80Pos->symbol(), symbol);
		symbol->setSize(5, 5);
		m_rind80Pos->setSymbol(symbol);
	}

	if (m_maxPos && m_stenosisLeftPos && m_stenosisPos && m_stenosisRightPos && m_referenceDiameter) {
		QPen pen = m_maxPos->pen();
		pen.setWidth(1);
		m_maxPos->setPen(pen);

		pen = m_stenosisLeftPos->pen();
		pen.setWidth(1);
		m_stenosisLeftPos->setPen(pen);
		m_stenosisRightPos->setPen(pen);

		pen = m_stenosisPos->pen();
		pen.setWidth(2);
		m_stenosisPos->setPen(pen);

		pen = m_referenceDiameter->pen();
		pen.setWidth(1);
		m_referenceDiameter->setPen(pen);
	}
}

void PlotWidget::clear()
{
    if(m_curve)     { m_curve->detach();     delete m_curve;     m_curve = 0;}
    if(m_slicePos)  { m_slicePos->detach();  delete m_slicePos;  m_slicePos = 0;}
    if(m_lind80Pos) { m_lind80Pos->detach(); delete m_lind80Pos; m_lind80Pos = 0;}
	if(m_lind50Pos) { m_lind50Pos->detach(); delete m_lind50Pos; m_lind50Pos = 0;}
    if(m_lind20Pos) { m_lind20Pos->detach(); delete m_lind20Pos; m_lind20Pos = 0;}
    if(m_rind80Pos) { m_rind80Pos->detach(); delete m_rind80Pos; m_rind80Pos = 0;}
	if(m_rind50Pos) { m_rind50Pos->detach(); delete m_rind50Pos; m_rind50Pos = 0;}
    if(m_rind20Pos) { m_rind20Pos->detach(); delete m_rind20Pos; m_rind20Pos = 0;}
    if(m_lminPos)   { m_lminPos->detach();   delete m_lminPos;   m_lminPos = 0;}
    if(m_rminPos)   { m_rminPos->detach();   delete m_rminPos;   m_rminPos = 0;}
    if(m_maxPos)    { m_maxPos->detach();    delete m_maxPos;    m_maxPos = 0;}
	if(m_stenosisPos) { m_stenosisPos->detach();    delete m_stenosisPos;    m_stenosisPos = 0;}
	if(m_stenosisLeftPos) { m_stenosisLeftPos->detach();    delete m_stenosisLeftPos;    m_stenosisLeftPos = 0;}
	if(m_stenosisRightPos) { m_stenosisRightPos->detach();    delete m_stenosisRightPos;    m_stenosisRightPos = 0;}
	if(m_referenceDiameter) { m_referenceDiameter->detach();    delete m_referenceDiameter;    m_referenceDiameter = 0;}
}

void PlotWidget::setXTitle(const QString& title)
{
    setAxisTitle(QwtPlot::xBottom, title);
	m_axisTitleFont = axisTitle(QwtPlot::xBottom).font();
}

void PlotWidget::setYTitle(const QString& title)
{
    setAxisTitle(QwtPlot::yLeft, title);
}


void PlotWidget::plotVesselProfile(const std::vector<float>& values, const Interval& pos)
{
    clear();

    double center = (double)values.size() / 2.0;

    std::vector<double> x, y;
    for(size_t i=0; i<values.size(); i++)
    {
        x.push_back( (double)i - center );
        y.push_back( values.at(i) );
    }
    plotFunction(x, y);

    QPen pen;
    pen.setColor( QColor(0,255,255) );
    pen.setStyle( Qt::SolidLine );
    m_maxPos  = plotPoint(pos.getMaximumIndex() - center,
                         pos.getMaximumValue(),
                         pen, QwtSymbol::NoSymbol);

    pen.setStyle( Qt::DashLine );
    m_lminPos = plotPoint((double)pos.getMinimumIndex(Interval::Left) - center,
                         pos.getMinimumValue(Interval::Left),
                         pen);

    m_rminPos = plotPoint((double)pos.getMinimumIndex(Interval::Right) - center,
                         pos.getMinimumValue(Interval::Right),
                         pen);

    pen.setStyle( Qt::SolidLine );
    pen.setColor(Qt::green);
    m_lind20Pos = plotPoint((double)pos.get20Index(Interval::Left) - center,
                            pos.get20Value(Interval::Left),
                            pen);
	m_lind80Pos = plotPoint((double)pos.get80Index(Interval::Left) - center,
                            pos.get80Value(Interval::Left),
                            pen);


    m_rind20Pos = plotPoint((double)pos.get20Index(Interval::Right) - center,
                            pos.get20Value(Interval::Right),
                            pen);
    m_rind80Pos = plotPoint((double)pos.get80Index(Interval::Right) - center,
                            pos.get80Value(Interval::Right),
                            pen);

	pen.setColor(Qt::white);
	m_lind50Pos = plotPoint((double)pos.get50Index(Interval::Left) - center,
                            pos.get50Value(Interval::Left),
                            pen);
	m_rind50Pos = plotPoint((double)pos.get50Index(Interval::Right) - center,
                            pos.get50Value(Interval::Right),
                            pen);
}

void PlotWidget::plotSlicePosition(double pos)
{
    if( m_slicePos )
    {
        m_slicePos->detach();
        delete m_slicePos;
    }

    QPen pen;
    pen.setStyle( Qt::DashDotLine );
    pen.setColor( QColor(0,255,255) );
    m_slicePos = plotLine( pos, m_minY, pos, m_maxY, pen );

    replot();
}

void PlotWidget::plotStenosisInfo(double minPos, double maxPos, double maxVal,
				double lPos, double rPos, double m, double t) {
	if (m_maxPos) {
		m_maxPos->detach();
		delete m_maxPos;
	}

	if (m_stenosisPos) {
		m_stenosisPos->detach();
		delete m_stenosisPos;
	}

	if (m_stenosisLeftPos) {
		m_stenosisLeftPos->detach();
		delete m_stenosisLeftPos;
	}

	if (m_stenosisRightPos) {
		m_stenosisRightPos->detach();
		delete m_stenosisRightPos;
	}

	if (m_referenceDiameter) {
		m_referenceDiameter->detach();
		delete m_referenceDiameter;
	}

	QPen pen;

	pen.setColor(Qt::white);
	pen.setStyle(Qt::SolidLine);
	pen.setWidth(2);
	m_stenosisPos = plotLine(minPos, m_minY, minPos, m_maxY, pen);

	double y1 = m * m_minX * 10.0 + t;
	double y2 = m * m_maxX * 10.0 + t;
	pen.setWidth(1);
	pen.setColor(Qt::magenta);
	m_referenceDiameter = plotLine(m_minX, y1, m_maxX, y2, pen);
	
	pen.setStyle(Qt::DashLine);
	pen.setColor(Qt::white);
	m_maxPos = plotLine(maxPos, maxVal, m_maxX, maxVal, pen);

	pen.setColor(Qt::yellow);
	m_stenosisLeftPos = plotLine(lPos, m_minY, lPos, m_maxY, pen);
	m_stenosisRightPos = plotLine(rPos, m_minY, rPos, m_maxY, pen);

	replot();
}

QwtPlotCurve* PlotWidget::plotPoint(double x1, double y1, const QPen& pen, QwtSymbol::Style style )
{
    QwtPlotCurve* line  = new QwtPlotCurve();
    line->setPen( pen );
    line->setPaintAttribute(QwtPlotCurve::ClipPolygons, true);
    line->setRenderHint(QwtPlotCurve::RenderAntialiased, true);
    line->setStyle( QwtPlotCurve::Sticks );

    QwtSymbol* symbol = new QwtSymbol(style);
    symbol->setSize(5,5);
    symbol->setPen( pen );
    line->setSymbol( symbol );
    line->attach(this);

    std::vector<double> ptX, ptY;
    ptX.push_back(x1);
    ptY.push_back(y1);
    line->setSamples(&ptX[0], &ptY[0], 1);
    line->plot();

    return line;
}

QwtPlotCurve* PlotWidget::plotLine(double x1, double y1, double x2, double y2, const QPen& pen )
{
    QwtPlotCurve* line  = new QwtPlotCurve();
    line->setPen( pen );
    line->setPaintAttribute(QwtPlotCurve::ClipPolygons, true);
    line->setRenderHint(QwtPlotCurve::RenderAntialiased, true);
    line->attach(this);

    std::vector<double> ptX, ptY;
    ptX.push_back(x1);
    ptX.push_back(x2);
    ptY.push_back(y1);
    ptY.push_back(y2);
    line->setSamples(&ptX[0], &ptY[0], 2);
    line->plot();

    return line;

}

void PlotWidget::plotFunction(const  std::vector<double>& x, const std::vector<double>& y)
{
    int c = 0;
	
	if( x.size() != y.size() || x.empty() )
    {
        if(m_curve)
        {
            m_curve->detach();
            delete m_curve;
            m_curve = 0;
            replot();
        }
        return;
    }

    m_minX = x[0];
    m_maxX = m_minX;
    m_minY = y[0];
    m_maxY = m_minY;

    for(size_t i=0; i<x.size(); i++)
    {
        m_maxY = std::max(y[i], m_maxY);
        m_minY = std::min(y[i], m_minY);

        m_maxX = std::max(x[i], m_maxX);
        m_minX = std::min(x[i], m_minX);
    }

    double magn = m_maxY - m_minY;
    m_minY -= magn*0.05;
    m_maxY += magn*0.05;
    if(m_minY < 0.0) m_minY = 0.0;

    QPen pen_curve;
    pen_curve.setStyle( Qt::SolidLine );
    if(c==0)
        pen_curve.setColor(Qt::red);
    if(c==1)
        pen_curve.setColor(Qt::green);
    if(c==2)
        pen_curve.setColor(Qt::blue);
    pen_curve.setWidth(2);

    if(m_curve)
    {
        m_curve->detach();
        delete m_curve;
    }
    m_curve = new QwtPlotCurve();
    m_curve->setPaintAttribute(QwtPlotCurve::ClipPolygons, true);
    m_curve->setRenderHint(QwtPlotCurve::RenderAntialiased, true);
    m_curve->attach(this);
    m_curve->setSamples(&x[0], &y[0], x.size());
    m_curve->setPen(pen_curve);

    setAxisScale(QwtPlot::xBottom, m_minX, m_maxX);
    setAxisScale(QwtPlot::yLeft,   m_minY, m_maxY);

    std::vector<double> xAxisX, xAxisY;
    xAxisX.push_back(m_minX);
    xAxisX.push_back(m_minX);
    xAxisX.push_back(m_maxX);
    xAxisY.push_back(m_maxY);
    xAxisY.push_back(m_minY);
    xAxisY.push_back(m_minY);
    m_axis->setSamples(&xAxisX[0], &xAxisY[0], 3);
    m_axis->plot();

    m_curve->plot();
}

//
//  Set a plain canvas frame and align the scales to it
//
void PlotWidget::alignScales()
{
    // The code below shows how to align the scales to
    // the canvas frame, but is also a good example demonstrating
    // why the spreaded API needs polishing.

    for ( size_t i = 0; i < QwtPlot::axisCnt; i++ )
    {
        QwtScaleWidget *scaleWidget = (QwtScaleWidget *)axisWidget(i);
        if ( scaleWidget )
            scaleWidget->setMargin(0);

        QwtScaleDraw *scaleDraw = (QwtScaleDraw *)axisScaleDraw(i);
        if ( scaleDraw )
            scaleDraw->enableComponent(QwtAbstractScaleDraw::Backbone, false);
    }

    plotLayout()->setAlignCanvasToScales(true);
}

void PlotWidget::savePlot(const QString& filename)
{
	QwtPlotRenderer renderer;
	renderer.setDiscardFlag(QwtPlotRenderer::DiscardBackground, false);
	setBlackOnWhite();
	renderer.renderDocument(this, filename, QSizeF(300, 200), 85);
	setWhiteOnBlack(true);
}

