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

#include "sharpnessplot.h"

#include "PlotWidget.h"
#include "settings.h"
#include "measurementpoint.h"

#include <QString>
#include <QFileDialog>
#include <QFileInfo>
#include <QImageWriter>
#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QStringList>

#include <iostream>
#include <fstream>
#include <string>


SharpnessPlot::SharpnessPlot(QWidget* parent) :
	m_settings(Settings::getSettings()),
    QWidget(parent)
{
    m_plot = new PlotWidget(this);
    m_plot->setStyleSheet("background-color: black;");
    m_plot->setXTitle("cm");
    m_plot->setYTitle("Sharpness");
    m_plot->show();

    QVBoxLayout* layout = new QVBoxLayout(this);
    layout->addWidget(m_plot);

	QHBoxLayout* hLayout = new QHBoxLayout();
	layout->addLayout(hLayout);
	hLayout->addStretch();

	m_lSharpness = new QLabel("Sharpness: 0.000");
	m_lSharpness->setAlignment( Qt::AlignCenter );
    m_lSharpness->setStyleSheet("color:rgb(255,255,255)");
    hLayout->addWidget(m_lSharpness);

	m_bExport = new QPushButton("Export data...");
	hLayout->addWidget(m_bExport);
	connect(m_bExport, SIGNAL(clicked()), this, SLOT(saveData()));

	m_bSavePlot = new QPushButton("Save plot...");
	hLayout->addWidget(m_bSavePlot);
	connect(m_bSavePlot, SIGNAL(clicked()), this, SLOT(savePlot()));

	setMode(Sharpness);

	m_lastDir = m_settings.getDefaultDir();
}

void SharpnessPlot::setMode(SharpnessPlot::Mode m)
{
    m_curMode = m;
    switch(m_curMode)
    {
        case Sharpness:
            m_lSharpness->setText("Sharpness: 0.000");
            m_plot->setXTitle("Position [cm]");
            m_plot->setYTitle("Sharpness [1/mm]");
            break;
        case Stenosis:
            m_lSharpness->setText("Diameter: 0.00 mm");
            m_plot->setXTitle("Position [cm]");
            m_plot->setYTitle("Diameter [mm]");
            break;
		case Diameter:
			m_lSharpness->setText("");
			m_plot->setXTitle("Position [cm]");
            m_plot->setYTitle("Diameter [mm]");
			break;
    }
}

void SharpnessPlot::setMeasurements(const std::vector<MeasurementPoint*>& measurements)
{
    m_measurements = measurements;
    update();
}

void SharpnessPlot::setStenosisInfo(double minPos, double maxPos, double maxVal,
				double lPos, double rPos, double m, double t) {
	m_plot->plotStenosisInfo(minPos / 10.0, maxPos / 10.0, maxVal, lPos / 10.0, rPos / 10.0, m, t);
}

void SharpnessPlot::update()
{
    std::vector<double> x, y;

    for(size_t i=0; i<m_measurements.size(); i++)
    {
        MeasurementPoint* p = m_measurements.at(i);
        x.push_back( p->getPositionInSpline() / 10.0 );
    }

    if(m_curMode == Sharpness)
    {
        for(size_t i=0; i<m_measurements.size(); i++)
        {
            MeasurementPoint* p = m_measurements.at(i);
            y.push_back( p->getAverageSharpness() );     // Plot the Sharpness along the spline
        }
    } else if ((m_curMode == Diameter) || (m_curMode == Stenosis)) {
        for(size_t i=0; i<m_measurements.size(); i++)
        {
            MeasurementPoint* p = m_measurements.at(i);
            y.push_back( p->getDiameter() );    // Plot the Diameter  along the spline
        }
    }

    m_plot->plotFunction( x, y );
	m_plot->replot();

	m_lastX = x;
	m_lastY = y;
}

void SharpnessPlot::setPosition(size_t pos)
{
    if(pos >= m_measurements.size() )
    {
        setMode(m_curMode);
        return;
    }
    MeasurementPoint* mp = m_measurements.at(pos);
    switch(m_curMode)
    {
        case Sharpness:
            m_lSharpness->setText( QString("Sharpness: %1").arg( mp->getAverageSharpness(), 0, 'f', 3, '0') );
            break;
        case Diameter:
		case Stenosis:
            m_lSharpness->setText( QString("Diameter: %1 mm").arg( mp->getDiameter(), 0, 'f', 2, '0') );
            break;
    }

    m_lSharpness->setStyleSheet("color:rgb(0,255,255)");
    double splinePos = mp->getPositionInSpline() / 10.0;
    m_plot->plotSlicePosition( splinePos );
}

void SharpnessPlot::saveData() {
	QString msg;
	switch(m_curMode)
	{
		case Sharpness:
			msg = "Save sharpness measurements";
			break;
		case Diameter:
			msg = "Save diameter measurements";
			break;
	}

	QString fileName = QFileDialog::getSaveFileName(this,
		msg, m_lastDir, tr("All Files (*)"));

	if (fileName.isEmpty())
		return;

	if (m_settings.getRememberLastDir()) {
		QFileInfo info(fileName);
		m_lastDir = info.path();
	}

	QByteArray fn = fileName.toLocal8Bit();
	std::ofstream file(fn.constData());
    if (!file.good()) {
		std::cerr << "Could not open file " << fileName.toStdString() << std::endl;
        return;
    }

	file.precision(6);

	std::vector<double>::const_iterator itX = m_lastX.begin();
	std::vector<double>::const_iterator itY = m_lastY.begin();

	for (; itX != m_lastX.end(); ++itX, ++itY)
		file << *itX * 10.0 << " " << *itY << std::endl;

	file.close();
}

void SharpnessPlot::savePlot() {
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

