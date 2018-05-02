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

#ifndef NORMALPLOT_H
#define NORMALPLOT_H

#include <QDialog>
#include <QWidget>
#include <QString>

#define cimg_display 0
#include "CImg.h"

#include "Wm3Vector3.h"
#include "Wm3Quaternion.h"

#include "settings.h"


// Forward declarations
class QLabel;
class QSlider;
class QPushButton;
class PlotWidget;
class MeasurementPoint;


class NormalPlot : public QDialog
{
    Q_OBJECT
public:
    explicit NormalPlot(QWidget* parent = 0);
	~NormalPlot();
    void setMeasurements(const std::vector<MeasurementPoint*>& measurements);

signals:
    
public slots:
    void setPosition(size_t pos);
    void setCurrentNormal(int val);
	void saveData();
	void savePlot();

private:
	Settings&				m_settings;

    size_t                          m_pos;
    size_t                          m_curNormal;
    std::vector<MeasurementPoint*>  m_measurements;

    PlotWidget*             m_plot;
    QLabel*                 m_lSharpness;
    QSlider*                m_slider;
	QPushButton*			m_saveData;
	QPushButton*			m_savePlot;
	QString					m_lastDir;

    void update();
};

#endif // NORMALPLOT_H

