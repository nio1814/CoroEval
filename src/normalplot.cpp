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

#include "normalplot.h"

#include "PlotWidget.h"
#include "measurementpoint.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QSlider>
#include <QLabel>
#include <QPushButton>
#include <QString>
#include <QStringList>
#include <QFileDialog>
#include <QFileInfo>
#include <QImageWriter>

#include <fstream>
#include <iostream>


NormalPlot::NormalPlot(QWidget *parent) :
    QDialog(parent),
    m_settings(Settings::getSettings())
{
    m_curNormal = 0;

    this->setMinimumWidth(800);

    m_plot = new PlotWidget(this);
    m_plot->setXTitle("Voxel");
    m_plot->setYTitle("Intensity");
    m_plot->show();

    QVBoxLayout* layout = new QVBoxLayout(this);
    layout->addWidget(m_plot);

    QWidget* w = new QWidget();
    QHBoxLayout* layout2 = new QHBoxLayout(w);

    m_slider = new QSlider(Qt::Horizontal);
    m_slider->setMinimum(0);
    m_slider->setMaximum(m_settings.getNumNormalsAtPoint() - 1);
    m_slider->setTickInterval(1);
    m_slider->setValue(0);
    layout2->addWidget(m_slider);
    connect(m_slider, SIGNAL(valueChanged(int)), this, SLOT(setCurrentNormal(int)));

    m_lSharpness = new QLabel("Sharpness: 0.000");
    m_lSharpness->setAlignment( Qt::AlignRight );
    layout2->addWidget(m_lSharpness);

    layout->addWidget(w);

	QWidget* w2 = new QWidget();
	QHBoxLayout* layout3 = new QHBoxLayout(w2);

	m_saveData = new QPushButton("Export data...", this);
	m_savePlot = new QPushButton("Save plot...", this);
	connect(m_saveData, SIGNAL(clicked()), this, SLOT(saveData()));
	connect(m_savePlot, SIGNAL(clicked()), this, SLOT(savePlot()));
	layout3->addStretch();
	layout3->addWidget(m_saveData);
	layout3->addWidget(m_savePlot);

	layout->addWidget(w2);

	m_lastDir = m_settings.getDefaultDir();
}

NormalPlot::~NormalPlot() {
	if (m_lSharpness)
		delete m_lSharpness;
	if (m_slider)
		delete m_slider;
	if (m_plot)
		delete m_plot;
}

void NormalPlot::setMeasurements(const std::vector<MeasurementPoint*>& measurements)
{
    m_measurements = measurements;

	if (!m_measurements.empty()) 
	{
		setPosition( m_measurements.size()-1 );
		
		if( measurements.size() == 1 )
			m_slider->setHidden(true);
		else
			m_slider->setHidden(false);

		update();
	} else {
		m_pos = 0;
	}
}

void NormalPlot::setPosition(size_t pos)
{
    m_pos = pos;
    update();
}

void NormalPlot::setCurrentNormal(int val)
{
    m_curNormal = static_cast<size_t>(val);
    update();
}

void NormalPlot::update()
{
    if(m_pos >= m_measurements.size() )
        return;

    MeasurementPoint* p = m_measurements.at(m_pos);
    this->setWindowTitle( QString("Normal of Vessel at %1 mm").arg( p->getPositionInSpline(), 0, 'f', 2, '0') );

    m_lSharpness->setText( QString("Sharpness: %1").arg( p->getSharpness(m_curNormal), 0, 'f', 3, '0') );
    m_plot->plotVesselProfile(p->getNormal(m_curNormal), p->getCenterInterval(m_curNormal) );
    m_plot->replot();
}

void NormalPlot::saveData()
{
	if(m_pos >= m_measurements.size() )
        return;

	QString fileName = QFileDialog::getSaveFileName(this,
		"Save profile data", m_lastDir, tr("All Files (*)"));

	if (fileName.isEmpty())
		return;

	if (m_settings.getRememberLastDir()) {
		QFileInfo info(fileName);
		m_lastDir = info.path();
	}

	MeasurementPoint* p = m_measurements[m_pos];
	const std::vector<float>& values = p->getNormal(m_curNormal);

	double center = static_cast<double>(values.size()) / 2.0;

	std::vector<double> x, y;
	for(size_t i= 0 ; i < values.size(); i++) {
		x.push_back(static_cast<double>(i) - center);
		y.push_back(values[i]);
	}

	QByteArray fn = fileName.toLocal8Bit();
	std::ofstream file(fn.constData());
    if (!file.good()) {
		std::cerr << "Could not open file " << fileName.toStdString() << std::endl;
        return;
    }

	file.precision(6);

	std::vector<double>::const_iterator itX = x.begin();
	std::vector<double>::const_iterator itY = y.begin();

	for (; itX != x.end(); ++itX, ++itY)
		file << *itX << " " << *itY << std::endl;

	file.close();
}

void NormalPlot::savePlot()
{
	const QList<QByteArray> imageFormats = QImageWriter::supportedImageFormats();
	
	QStringList filter;
	filter += "PDF Documents (*.pdf)";
#ifndef QWT_NO_SVG
	filter += "SVG Documents (*.svg)";
#endif
	filter += "Postscript Documents (*.ps)";

	if (imageFormats.size() > 0) {
		QString imageFilter("Images (");
		for (int i = 0; i < imageFormats.size(); i++) {
			if (i > 0)
				imageFilter += " ";

			imageFilter += "*.";
			imageFilter += imageFormats[i];
		}
		imageFilter += ")";

		filter += imageFilter;
	}

	QString fileName = QFileDialog::getSaveFileName(this, "Save plot as", m_lastDir,
							filter.join(";;"));

	if (!fileName.isEmpty()) {
		if (m_settings.getRememberLastDir()) {
			QFileInfo info(fileName);
			m_lastDir = info.path();
		}

		m_plot->savePlot(fileName);
    }
}

