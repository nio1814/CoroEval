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

#ifndef STENOSISDIALOG_H
#define STENOSISDIALOG_H

#include <QDialog>

#include <vector>


// Forward declarations
class BSplineInterpolation;
class Data;
class MeasurementPoint;
class Settings;

namespace Ui {
class StenosisDialog;
}


class StenosisDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit StenosisDialog(QWidget* parent = 0);
    ~StenosisDialog();

	void reset();    
    void setInterpolator(BSplineInterpolation* bSpline);
    void setData(Data* d);
	void updateMax();
	void setMeasurements(const std::vector<MeasurementPoint*>& m);

public slots:
    void evaluate();

signals:

private:
	const Settings& m_settings;
    Ui::StenosisDialog* ui;

    Data* m_data;
    BSplineInterpolation* m_bSpline;
	std::vector<MeasurementPoint*> m_measurements;

	void fitLine(const std::vector<double>& x, const std::vector<double>& y, double& m, double& t);
};

#endif // StenosisDialog_H

