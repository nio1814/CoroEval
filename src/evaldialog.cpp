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

#include "evaldialog.h"

#include "measurementpoint.h"
#include "imagelabel.h"
#include "normalplot.h"
#include "Data.h"

#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QDoubleSpinBox>


EvalDialog::EvalDialog(char* type, std::vector<MeasurementPoint*> measurements, cimg_library::CImg<float>& cImage)
{
	m_type = QString( type );
	int width  = cImage.width();
	int height = cImage.height();

	m_measurements = measurements;

	float* image = cImage.data();
	Data* d = new Data(image);
	d->create( Data::COL, width, Data::LIN, height);
	d->setVoxelSize( 1.05 );
	
	QHBoxLayout* layout = new QHBoxLayout( this );
	
	QWidget* w = new QWidget(this);
	layout->addWidget(w);

	QVBoxLayout* layout2 = new QVBoxLayout( w );
	m_label = new ImageLabel();
	m_label->setInputData( d );
	m_label->setMode( ImageLabel::Segmentation );
	m_label->showGrid( true );
	layout2->addWidget(m_label);
	
	m_spinBox = new QDoubleSpinBox();
	m_spinBox->setSingleStep(0.5);
	m_spinBox->setMaximum( m_measurements.back()->getPositionInSpline() );
	layout2->addWidget(m_spinBox);

	float maxVal;
	float minVal = cImage.min_max(maxVal);
	float wcenter = (maxVal + minVal) / 2;
	m_label->updateWindowing( maxVal-minVal, wcenter );


	m_normalPlot = new NormalPlot();
	layout->addWidget( m_normalPlot );

	QObject::connect( m_label, SIGNAL( clickedAt(Wm3::Vector3<double>)), this, SLOT(imageClicked(Wm3::Vector3<double>) ) );
	QObject::connect( m_spinBox, SIGNAL( valueChanged(double)), this, SLOT( changePosition(double) ) );
	imageClicked( Wm3::Vector3<double>( 0, m_measurements.size() / 2, 0) );
}

void EvalDialog::imageClicked(Wm3::Vector3<double> pt)
{
	if( pt.Y() < 0 )
		return;
	if( pt.Y() > m_measurements.size()-1 )
		return;

	std::vector<MeasurementPoint*> meas;
	meas.push_back( m_measurements.at( pt.Y() ) );
	m_normalPlot->setMeasurements( meas );	
	m_label->updateGridPos( 0, pt.Y() );
	m_spinBox->setValue( m_measurements.at( pt.Y() )->getPositionInSpline() );
	this->setWindowTitle( QString("%1: Normal of Vessel at %2 mm").arg(m_type).arg( m_measurements.at( pt.Y() )->getPositionInSpline(), 0, 'f', 2, '0') );
}

void EvalDialog::changePosition(double val)
{
	Wm3::Vector3<double> pt( 0, val*2, 0 );
	emit( imageClicked( pt ));
}

