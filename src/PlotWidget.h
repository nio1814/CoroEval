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

#ifndef _PlotWidget_H_
#define _PlotWidget_H_

#include <qwt_plot.h>
#include <qwt_symbol.h>

#include <QList>
#include <QString>
#include <QFont>

#include <vector>


// Forward declarations
class QwtPlotGrid;
class QwtPlotCurve;
class Interval;
class QPen;


class PlotWidget: public QwtPlot
{
    Q_OBJECT

public:

    PlotWidget(QWidget* = NULL);
    ~PlotWidget();

    void plotFunction(const std::vector<double>& x, const std::vector<double>& y);
    void plotVesselProfile(const std::vector<float>& values, const Interval& pos);
    void plotSlicePosition(double pos);
	
	void plotStenosisInfo(double stenPos, double maxPos, double maxVal,
				double lPos, double rPos, double m, double t);

    void setXTitle(const QString& title);
    void setYTitle(const QString& title);

	void savePlot(const QString& filename);

private:
    void clear();
    void alignScales();
	void setWhiteOnBlack(bool resetFonts = false);
	void setBlackOnWhite();

    QwtPlotGrid*    m_grid;
    QwtPlotCurve*   m_curve;
    QwtPlotCurve*   m_axis;

    QwtPlotCurve*   m_slicePos;
    QwtPlotCurve*   m_lminPos;
    QwtPlotCurve*   m_lind80Pos;
	QwtPlotCurve*   m_lind50Pos;
    QwtPlotCurve*   m_lind20Pos;
    QwtPlotCurve*   m_maxPos;
    QwtPlotCurve*   m_rind80Pos;
	QwtPlotCurve*   m_rind50Pos;
    QwtPlotCurve*   m_rind20Pos;
    QwtPlotCurve*   m_rminPos;

	QwtPlotCurve*	m_stenosisPos;
	QwtPlotCurve*	m_stenosisLeftPos;
	QwtPlotCurve*	m_stenosisRightPos;
	QwtPlotCurve*	m_referenceDiameter;

    QwtPlotCurve* plotPoint(double x1, double y1, const QPen& pen, QwtSymbol::Style style = QwtSymbol::Ellipse);
    QwtPlotCurve* plotLine(double x1, double y1, double x2, double y2, const QPen& pen);

	QFont			m_axisFont;
	QFont			m_axisTitleFont;

    double          m_minX;
    double          m_maxX;
    double          m_minY;
    double          m_maxY;
};

#endif

