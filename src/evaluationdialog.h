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

#ifndef EVALUATIONDIALOG_H
#define EVALUATIONDIALOG_H

#include "settings.h"

#include <QDialog>


// Forward declarations
class BSplineInterpolation;
class Data;
class MeasurementPoint;

namespace Ui {
class EvaluationDialog;
}


class EvaluationDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit EvaluationDialog(QWidget* parent = 0);
    ~EvaluationDialog();

	void reset();    
    void setInterpolator(BSplineInterpolation* bSpline);
    void setData(Data* d);
	void setMeasurements(const std::vector<MeasurementPoint*>& m);
	void updateMax();

public slots:
    void evaluate();
    void evaluationStartChanged(size_t p);
    void evaluationEndChanged(size_t p);
    void evaluationStartChanged(double p);
    void evaluationEndChanged(double p);

private:
	const Settings& m_settings;
    Ui::EvaluationDialog* ui;

    Data* m_data;
    BSplineInterpolation* m_bSpline;

	std::vector<MeasurementPoint*> m_measurements;

    size_t                         m_startPos;
    std::vector<MeasurementPoint*> m_measurementsStart;

    size_t                         m_endPos;
    std::vector<MeasurementPoint*> m_measurementsEnd;

    std::vector<MeasurementPoint*> evaluationHelper(double start, double end);
    std::vector<MeasurementPoint*> evaluationHelper(double start, double end, double pos, size_t& vecPos);
};

#endif // EvaluationDialog_H

