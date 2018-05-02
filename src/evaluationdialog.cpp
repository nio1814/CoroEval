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

#include "evaluationdialog.h"
#include "ui_evaluationdialog.h"

#include "bsplineinterpolation.h"
#include "Data.h"
#include "measurementpoint.h"


EvaluationDialog::EvaluationDialog(QWidget *parent) :
    QDialog(parent),
    m_settings(Settings::getSettings()),
    ui(new Ui::EvaluationDialog)
{
    m_startPos = 1;
    m_endPos   = 1;

    ui->setupUi(this);
    setWindowTitle("Vessel Sharpness Evaluation");

	ui->sb_start->setSingleStep(m_settings.getSharpnessEvalInterval());
	ui->sb_end->setSingleStep(m_settings.getSharpnessEvalInterval());

    connect(ui->profileStart, SIGNAL(positionChanged(size_t)), this, SLOT(evaluationStartChanged(size_t)));
    connect(ui->profileEnd,   SIGNAL(positionChanged(size_t)), this, SLOT(evaluationEndChanged(size_t)));
    connect(ui->sb_start, SIGNAL(valueChanged(double)), this, SLOT(evaluationStartChanged(double)));
    connect(ui->sb_end, SIGNAL(valueChanged(double)),   this, SLOT(evaluationEndChanged(double)));

    connect(ui->pushButton, SIGNAL(clicked()), this, SLOT(evaluate()));
}

EvaluationDialog::~EvaluationDialog()
{
    delete ui;
}

void EvaluationDialog::reset()
{
	m_data = NULL;
	m_bSpline = NULL;
	ui->profileEnd->reset();
	ui->profileStart->reset();

	ui->sb_start->setMaximum(0);
	ui->sb_end->setMaximum(0);
}

void EvaluationDialog::setMeasurements(const std::vector<MeasurementPoint*>& m)
{
	m_measurements = m;
	evaluationStartChanged(ui->sb_start->value());
	evaluationEndChanged(ui->sb_end->value());
}

void EvaluationDialog::updateMax()
{
	if ((m_data && m_bSpline) && (m_bSpline->getLength() > 0)) {
		double max = m_bSpline->getLength() * m_data->getVoxelSize();
		if (!m_measurements.empty())
			max = m_measurements[m_measurements.size() - 1]->getPositionInSpline();

		ui->sb_end->setMaximum(max);
        ui->sb_start->setMaximum(max);

		ui->sb_end->setValue(max);

		ui->sb_start->setSingleStep(m_settings.getSharpnessEvalInterval());
		ui->sb_end->setSingleStep(m_settings.getSharpnessEvalInterval());
	} else {
		ui->sb_start->setMaximum(0);
		ui->sb_end->setMaximum(0);
	}
}

void EvaluationDialog::setInterpolator(BSplineInterpolation *bSpline)
{
    m_bSpline = bSpline;
}

void EvaluationDialog::setData(Data *d)
{
    m_data = d;
    ui->profileEnd->setData(d);
    ui->profileEnd->setProfileSize(15);
    ui->profileStart->setData(d);
    ui->profileStart->setProfileSize(15);
}

void EvaluationDialog::evaluate()
{
    if( !m_bSpline || !m_data )
        return;

    if(m_bSpline->getLength() < 0)
        return;

	if (m_measurements.empty())
		return;

	double start = ui->sb_start->value();
	double end   = ui->sb_end->value();
	std::vector<MeasurementPoint*> evalInterval = evaluationHelper(start, end);

    double mean_sharpness = 0;
    double mean_diameter  = 0;

    for(size_t i=0; i<evalInterval.size(); i++)
    {
        MeasurementPoint* mp = evalInterval.at(i);
        mean_sharpness += mp->getAverageSharpness();
        mean_diameter  += mp->getDiameter();
    }
    mean_sharpness /= double( evalInterval.size() );
    mean_diameter  /= double( evalInterval.size() );

    double std_diameter = 0;
	double std_sharpness  = 0;
    for(size_t i=0; i<evalInterval.size(); i++)
    {
        MeasurementPoint* mp = evalInterval.at(i);
        double value = mp->getDiameter() - mean_diameter;
        std_diameter   += (value*value);
		value = mp->getAverageSharpness() - mean_sharpness;
		std_sharpness  += (value*value);
    }
	if (evalInterval.size() > 1) {
		std_diameter = sqrt(std_diameter / double( evalInterval.size() - 1));
		std_sharpness = sqrt(std_sharpness / double(evalInterval.size() - 1));
	} else {
		std_diameter = sqrt(std_diameter / double( evalInterval.size() ));
		std_sharpness = sqrt(std_sharpness / double(evalInterval.size()));
	}

    ui->l_sharpness->setText( QString("%1 %2 %3").arg(mean_sharpness, 0, 'f', 3, '0').arg(QChar(177)).arg(std_sharpness,   0, 'f', 3, '0') );
    ui->l_diameter->setText(  QString("%1 %2 %3 mm").arg(mean_diameter,  0, 'f', 2, '0').arg(QChar(177)).arg(std_diameter, 0, 'f', 2, '0') );
	
	double len = evalInterval[evalInterval.size() - 1]->getPositionInSpline() - evalInterval[0]->getPositionInSpline();
	ui->l_length->setText(QString("%1 mm").arg(len,  0, 'f', 1, '0'));
}

std::vector<MeasurementPoint*> EvaluationDialog::evaluationHelper(double start, double end)
{
    size_t p;
    return evaluationHelper(start, end, start, p);
}

std::vector<MeasurementPoint*> EvaluationDialog::evaluationHelper(double start, double end, double pos, size_t &vecPos)
{
    std::vector<MeasurementPoint*> res;
	
	if( !m_bSpline || !m_data )
        return res;

    if(m_bSpline->getLength() < 0)
        return res;

	if (m_measurements.empty())
		return res;

	double length = m_bSpline->getLength() * m_data->getVoxelSize();
	double delta = m_settings.getSharpnessEvalInterval() / 2.0;
	
    start = std::max<double>(0.0, start);
    end   = std::min<double>(length, end);

	for(std::vector<MeasurementPoint*>::iterator i = m_measurements.begin();
		i != m_measurements.end();
		++i) {
		MeasurementPoint* mp = *i;

		double curPos = mp->getPositionInSpline(); 

		if (fabs(curPos - pos) < delta)
			vecPos = res.size();

		if (((fabs(curPos - start) < delta) || (curPos >= start)) &&
			((fabs(curPos - end) < delta) || (curPos <= end)))
			res.push_back(mp);
	}

	return res;
}

void EvaluationDialog::evaluationStartChanged(size_t p)
{
    m_startPos = p;
}

void EvaluationDialog::evaluationEndChanged(size_t p)
{
    m_endPos = p;
}

void EvaluationDialog::evaluationStartChanged(double p)
{
    if( !m_bSpline || !m_data )
        return;

    if( m_bSpline->getLength() < 0 )
        return;

    double start = p - 2.0;
    double end   = p + 2.0;

    size_t x;
    m_measurementsStart = evaluationHelper(start, end, p, x);
    ui->profileStart->setMeasurements(m_measurementsStart);
	if (!m_measurementsStart.empty())
		ui->profileStart->updateProfile(x);
}

void EvaluationDialog::evaluationEndChanged(double p)
{
    if( !m_bSpline || !m_data )
        return;
    if( m_bSpline->getLength() < 0 )
        return;

    double start = p - 2.0;
    double end   = p + 2.0;

    size_t x;
    m_measurementsEnd = evaluationHelper(start, end, p, x);
    ui->profileEnd->setMeasurements(m_measurementsEnd);
	if (!m_measurementsEnd.empty())
		ui->profileEnd->updateProfile(x);
}

