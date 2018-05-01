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

#ifndef EVALDIALOG_H
#define EVALDIALOG_H

#include <QDialog>
#include <QString>

#define cimg_display 0
#include "CImg.h"

#include "Wm3Vector3.h"


// Forward declarations
class QDoubleSpinBox;
class MeasurementPoint;
class ImageLabel;
class NormalPlot;


class EvalDialog : public QDialog
{
    Q_OBJECT
public:
	EvalDialog(char* type, std::vector<MeasurementPoint*> measurements, cimg_library::CImg<float>& cImage);

public slots:
	void imageClicked(Wm3::Vector3<double> pt);
	void changePosition(double val);
private:
	std::vector<MeasurementPoint*> m_measurements;
	ImageLabel*		m_label;
	NormalPlot*		m_normalPlot;
	QDoubleSpinBox* m_spinBox;
	QString			m_type;
};

#endif

