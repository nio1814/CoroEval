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

#ifndef SHARPNESSPLOT_H
#define SHARPNESSPLOT_H

#include "Wm3Vector3.h"
#include "Wm3Quaternion.h"

#include <QWidget>
#include <QString>

#include <vector>


// Forward declarations
class MeasurementPoint;
class Settings;
class PlotWidget;
class QLabel;
class QPushButton;


class SharpnessPlot : public QWidget
{
    Q_OBJECT
public:
    enum Mode{
        Sharpness,
        Diameter,
		Stenosis
    };

    explicit SharpnessPlot(QWidget* parent = 0);

    void setMode(SharpnessPlot::Mode m);
    void setMeasurements(const std::vector<MeasurementPoint*>& measurements);
	void setStenosisInfo(double minPos, double maxPos, double maxVal,
				double lPos, double rPos, double m, double t);
    void update();

signals:
    
public slots:
    void setPosition(size_t pos);

private:
	Settings&                      m_settings;
    SharpnessPlot::Mode            m_curMode;
    std::vector<MeasurementPoint*> m_measurements;
    PlotWidget*                    m_plot;
    QLabel*                        m_lSharpness;
	QPushButton*                   m_bExport;
	QPushButton*                   m_bSavePlot;
	std::vector<double>            m_lastX;
	std::vector<double>            m_lastY;
	QString                        m_lastDir;

private slots:
	void saveData();
	void savePlot();
};

#endif // SHARPNESSPLOT_H

