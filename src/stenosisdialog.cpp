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


/*
 * NOTE
 * This feature is untested and therefore not enabled in the software.
 * Do not use without verification!
 * NOTE
 */

#include "stenosisdialog.h"
#include "ui_stenosisdialog.h"

#include "os.h"
#include "bsplineinterpolation.h"
#include "measurementpoint.h"
#include "Data.h"
#include "settings.h"
#include "MPFilter.h"

#include <QMessageBox>

#include <cmath>
#include <algorithm>


#undef min
#undef max


StenosisDialog::StenosisDialog(QWidget* parent)
	:	m_settings(Settings::getSettings()),
		QDialog(parent),
		ui(new Ui::StenosisDialog)
{
    ui->setupUi(this);

	ui->plot->setMode(SharpnessPlot::Diameter);

    connect(ui->b_evaluate, SIGNAL(clicked()), this, SLOT(evaluate()));

	ui->sb_start->setSingleStep(m_settings.getSharpnessEvalInterval());
	ui->sb_end->setSingleStep(m_settings.getSharpnessEvalInterval());
	ui->sb_stenosisStart->setSingleStep(m_settings.getSharpnessEvalInterval());
	ui->sb_stenosisEnd->setSingleStep(m_settings.getSharpnessEvalInterval());
}

StenosisDialog::~StenosisDialog()
{
    delete ui;
}

void StenosisDialog::reset()
{
	m_data = NULL;
	m_bSpline = NULL;
}

void StenosisDialog::setInterpolator(BSplineInterpolation* bSpline)
{
    m_bSpline = bSpline;
}

void StenosisDialog::setData(Data* d)
{
    m_data = d;

	if ((m_data && m_bSpline) && (m_bSpline->getLength() > 0)) {
		ui->sb_end->setValue(m_bSpline->getLength() * m_data->getVoxelSize());
	}
}

void StenosisDialog::updateMax()
{
	if (m_data && m_bSpline && (m_bSpline->getLength() > 0)) {
		double max = m_bSpline->getLength() * m_data->getVoxelSize();
		if (!m_measurements.empty())
			max = m_measurements[m_measurements.size() - 1]->getPositionInSpline();
		
		ui->sb_end->setMaximum(max);
        ui->sb_start->setMaximum(max);
		ui->sb_stenosisEnd->setMaximum(max);
		ui->sb_stenosisStart->setMaximum(max);
		
		ui->sb_end->setValue(max);

		ui->sb_start->setSingleStep(m_settings.getSharpnessEvalInterval());
		ui->sb_end->setSingleStep(m_settings.getSharpnessEvalInterval());
		ui->sb_stenosisStart->setSingleStep(m_settings.getSharpnessEvalInterval());
		ui->sb_stenosisEnd->setSingleStep(m_settings.getSharpnessEvalInterval());
	}
}

void StenosisDialog::setMeasurements(const std::vector<MeasurementPoint*>& m)
{
	m_measurements = m;
	
	ui->plot->setMeasurements(m_measurements);
	ui->plot->update();
}

void StenosisDialog::evaluate()
{
    if( !m_bSpline || !m_data )
        return;

    if(m_bSpline->getLength() < 0)
        return;

	if (m_measurements.empty())
		return;

	double start = ui->sb_start->value();
	double end = ui->sb_end->value();
	double stenosisStart = ui->sb_stenosisStart->value();
	double stenosisEnd = ui->sb_stenosisEnd->value();

	if ((start >= end) || (stenosisStart >= stenosisEnd) || (start > stenosisStart) || (stenosisEnd > end)) {
		QMessageBox::critical(this, "Error", "Invalid evaluation range");
		return;
	}

	double delta = m_settings.getSharpnessEvalInterval() / 2.0;

	// Find start position
	size_t startPos = 0;
	while ((startPos < m_measurements.size()) && (fabs(m_measurements[startPos]->getPositionInSpline() - start) >= delta))
		++startPos;
	if (startPos == m_measurements.size()) {
		QMessageBox::critical(this, "Error", "Invalid evaluation range");
		return;
	}

	// Find stenosis start position
	size_t stenosisStartPos = startPos;
	while ((stenosisStartPos < m_measurements.size()) && (fabs(m_measurements[stenosisStartPos]->getPositionInSpline() - stenosisStart) >= delta))
		++stenosisStartPos;
	if (stenosisStartPos == m_measurements.size()) {
		QMessageBox::critical(this, "Error", "Invalid evaluation range");
		return;
	}

	// Find stenosis end position
	size_t stenosisEndPos = stenosisStartPos + 1;
	while ((stenosisEndPos < m_measurements.size()) && (fabs(m_measurements[stenosisEndPos]->getPositionInSpline() - stenosisEnd) >= delta))
		++stenosisEndPos;
	if (stenosisEndPos == m_measurements.size()) {
		QMessageBox::critical(this, "Error", "Invalid evaluation range");
		return;
	}

	// Find end position
	size_t endPos = stenosisEndPos;
	while ((endPos < m_measurements.size()) && (fabs(m_measurements[endPos]->getPositionInSpline() - end) >= delta))
		++endPos;
	if (endPos == m_measurements.size())
		--endPos;

	// Create actual measurements used for the evaluation
	std::vector<MeasurementPoint*> measurements;
	for (size_t i = startPos; i <= endPos; ++i)
		measurements.push_back(m_measurements[i]);

	stenosisStartPos -= startPos;
	stenosisEndPos -= startPos;
	endPos -= startPos;
	startPos = 0;

	// Find minimum
	size_t minPos = stenosisStartPos;
	double min = measurements[stenosisStartPos]->getDiameter();
	for (size_t pos = minPos + 1; pos <= stenosisEndPos; ++pos)
		if (measurements[pos]->getDiameter() < min) {
			minPos = pos;
			min = measurements[pos]->getDiameter();
		}

	ui->l_min->setText(QString("%1 mm").arg(min, 0, 'f', 2, '0'));
	ui->l_stenosisPos->setText(QString("%1 mm").arg(measurements[minPos]->getPositionInSpline(), 0, 'f', 1, '0'));
	ui->l_stenosisLength->setText(QString("%1 mm").arg(
					measurements[stenosisEndPos]->getPositionInSpline() - measurements[stenosisStartPos]->getPositionInSpline(),
					0, 'f', 1, '0'));
	
	// Evaluate general statistics
	size_t maxPos = startPos;
	double max = measurements[startPos]->getDiameter();
	for (size_t pos = startPos + 1; pos <= endPos; ++pos) {
		double v = measurements[pos]->getDiameter();
		if (v > max) {
			maxPos = pos;
			max = v;
		}
	}

	double mean = measurements[stenosisStartPos]->getDiameter();
	for (size_t pos = stenosisStartPos + 1; pos <= stenosisEndPos; ++pos) {
		double v = measurements[pos]->getDiameter();
		mean += v;
	}
	mean /= static_cast<double>(stenosisEndPos - stenosisStartPos + 1);

	double std = 0;
	for (size_t pos = stenosisStartPos; pos <= stenosisEndPos; ++pos) {
		double v = measurements[pos]->getDiameter() - mean;
		std += v*v;
	}
	std = std::sqrt(std / static_cast<double>(stenosisEndPos - stenosisStartPos));

	ui->l_max->setText(QString("%1 mm").arg(max, 0, 'f', 2, '0'));
	ui->l_mean->setText(QString("%1 mm").arg(mean, 0, 'f', 2, '0'));
	ui->l_std->setText(QString("%1 mm").arg(std, 0, 'f', 2, '0'));

	// Calculate reference diameter
	std::vector<double> x;
	std::vector<double> y;

	if (ui->cb_80fit->isChecked()) {
		std::vector<double> xOutside;
		std::vector<double> yOutside;
		
		for (size_t i = startPos; i <= endPos; ++i) {
			if ((i >= stenosisStartPos) && (i <= stenosisEndPos))
				continue;

			xOutside.push_back(measurements[i]->getPositionInSpline());
			yOutside.push_back(measurements[i]->getDiameter());
		}

		std::vector<double> ySort(yOutside);
		std::sort(ySort.begin(), ySort.end());

		/*size_t idx = static_cast<size_t>(std::floor(static_cast<float>(ySort.size()) * 0.8f));
		double maxVal = ySort[idx];
		
		for (size_t i = 0; i < xOutside.size(); ++i) {
			if (yOutside[i] <= maxVal) {
				x.push_back(xOutside[i]);
				y.push_back(yOutside[i]);
			}
		}*/

		size_t idxL = static_cast<size_t>(std::floor(static_cast<float>(ySort.size()) * 0.1f));
		size_t idxU = static_cast<size_t>(std::floor(static_cast<float>(ySort.size()) * 0.9f));
		double minVal = ySort[idxL];
		double maxVal = ySort[idxU];
		
		for (size_t i = 0; i < xOutside.size(); ++i) {
			if ((yOutside[i] >= minVal) && (yOutside[i] <= maxVal)) {
				x.push_back(xOutside[i]);
				y.push_back(yOutside[i]);
			}
		}
	} else {
		for (size_t i = startPos; i <= endPos; ++i) {
			if ((i >= stenosisStartPos) && (i <= stenosisEndPos))
				continue;

			x.push_back(measurements[i]->getPositionInSpline());
			y.push_back(measurements[i]->getDiameter());
		}
	}

	double m, t;
	fitLine(x, y, m, t);

	ui->l_ref_m->setText(QString("%1").arg(m, 0, 'f', 2, '0'));
	ui->l_ref_t->setText(QString("%1").arg(t, 0, 'f', 2, '0'));
	
	double ref = measurements[minPos]->getPositionInSpline() * m + t;
	ui->l_ref->setText(QString("%1 mm").arg(ref, 0, 'f', 2, '0'));

	// Calculate stenosis degree
	double degree = (1.0 - min / ref) * 100.0;
	ui->l_stenosisDegree->setText(QString("%1 %").arg(round(degree)));

	ui->plot->setMeasurements(measurements);
    ui->plot->update();

	ui->plot->setStenosisInfo(measurements[minPos]->getPositionInSpline(),
							measurements[maxPos]->getPositionInSpline(), max,
							measurements[stenosisStartPos]->getPositionInSpline(),
							measurements[stenosisEndPos]->getPositionInSpline(),
							m, t);
}

void StenosisDialog::fitLine(const std::vector<double>& x, const std::vector<double>& y, double& m, double& t) {
	if ((x.size() != y.size()) || x.size() < 2) {
		m = 0;
		t = 0;
		return;
	}

	double n = x.size();

	// Mean values
	double xMean = 0, yMean = 0;
	std::vector<double>::const_iterator i = x.begin();
	std::vector<double>::const_iterator j = y.begin();
	for (; i != x.end(); ++i, ++j) {
		xMean += *i;
		yMean += *j;
	}
	xMean /= n;
	yMean /= n;

	// Variance of x values and covariance of x and y
	double varX = 0, covXY = 0;
	i = x.begin();
	j = y.begin();
	for (; i != x.end(); ++i, ++j) {
		varX += (*i - xMean) * (*i - xMean);
		covXY += (*i - xMean) * (*j - yMean);
	}
	varX /= (n - 1);
	covXY /= (n - 1);

	// Linear fit
	m = covXY / varX;
	t = yMean - m * xMean;
}
