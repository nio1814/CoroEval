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

#include "imagewidget.h"
#include "ui_imagewidget.h"

#include "os.h"

#include "Data.h"

#include <QSlider>


ImageWidget::ImageWidget(QWidget* parent) :
    QWidget(parent),
    ui(new Ui::ImageWidget)
{
    m_data = 0;
    ui->setupUi(this);

    ui->imageLabel->setBackgroundRole(QPalette::Base);
    ui->imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    ui->imageLabel->setScaledContents(true);

    m_horizontalSlider = ui->horizontalSlider;
    m_verticalSlider   = ui->verticalSlider;

    m_horizontalSlider->setVisible(false);
    m_verticalSlider->setVisible(false);

    connect(m_horizontalSlider, SIGNAL(valueChanged(int)), this, SLOT(updateGridPos()) );
    connect(m_verticalSlider,   SIGNAL(valueChanged(int)), this, SLOT(updateGridPos()) );

    connect(ui->imageLabel, SIGNAL(clickedAt(const Wm3::Vector3<double>&)), this, SLOT(clickedAt(const Wm3::Vector3<double>&)));
    connect(ui->imageLabel, SIGNAL(sliceChanged(double)), this, SLOT(changeSlice(double)));

    connect(ui->imageLabel, SIGNAL(scaleChanged(double)), this, SLOT(sendScale(double)));
    connect(ui->imageLabel, SIGNAL(offsetChanged(const Wm3::Vector3<double>&)), this, SLOT(sendOffset(const Wm3::Vector3<double>&)));
}

ImageWidget::~ImageWidget()
{
    delete ui;
    delete m_data;
}

void ImageWidget::reset()
{
	if (m_data)
		delete m_data;
	m_data = NULL;
	ui->imageLabel->reset();
}

void ImageWidget::setMode(ImageLabel::Mode m)
{
    ui->imageLabel->setMode(m);
}

void ImageWidget::setInputData(Data* img, const Wm3::Matrix3<double>& rotMat)
{
    m_rotMat = rotMat;

    if( m_data )
		delete m_data;
    m_data = new Data();
    m_data->clone(img);
    m_data->setPar(m_data->getLen(Data::PAR) / 2);
    m_data->setOrientation(m_rotMat);

    int col = m_data->getLen(Data::COL);
    int lin = m_data->getLen(Data::LIN);

    this->setMinimumSize(col, lin);
    this->setMaximumSize(col, lin);

    ui->imageLabel->setInputData(m_data);
    ui->imageLabel->setMinimumSize( col, lin );
    ui->imageLabel->setMaximumSize( col, lin );

    ui->horizontalSlider->setMinimum(0);
    ui->horizontalSlider->setMaximum(col-1);

    ui->verticalSlider->setMinimum(0);
    ui->verticalSlider->setMaximum(lin-1);

    ui->verticalSlider->setValue( (lin-1)/2 );
    ui->horizontalSlider->setValue( (col-1)/2 );

    ui->scrollArea->setMinimumSize(col, lin);
    ui->scrollArea->setMaximumSize(col, lin);
}

QSlider* ImageWidget::getHorizontalSlider() const
{
    return m_horizontalSlider;
}

QSlider* ImageWidget::getVerticalSlider() const
{
    return m_verticalSlider;
}

void ImageWidget::changeSlice(const Wm3::Vector3<double>& pt)
{
    Wm3::Vector3<double> p = pt * m_rotMat;
    ui->horizontalSlider->setValue( round(p.X()) );
    ui->verticalSlider->setValue( round(p.Y()) );
    changeSlice( p.Z() );
}

void ImageWidget::changeSlice(int val)
{
    changeSlice( double(val) );
}

void ImageWidget::changeSlice(double val)
{
    m_data->setPar(val);
    ui->imageLabel->updateBaseImage();
    updateOverlay();

    emit sliceChanged( int(val) );
}

void ImageWidget::updateControlPoints(const std::vector< Wm3::Vector3<double> >& pts)
{
    m_cntPts.clear();
    for(size_t i=0; i<pts.size(); i++)
    {
        m_cntPts.push_back(pts.at(i) * m_rotMat);
    }
    updateOverlay();
}

void ImageWidget::updateSplinePoints(const std::vector< Wm3::Vector3<double> >& pts)
{
    m_splPts.clear();
    for(size_t i=0; i<pts.size(); i++)
    {
        m_splPts.push_back(pts.at(i) * m_rotMat);
    }
    updateOverlay();
}

void ImageWidget::updateOverlay() const
{
    ui->imageLabel->clearOverlayImage();

    QRgb color = qRgb(0, 255, 255);
    size_t size = m_splPts.size();
    for (size_t i=0; i<size; i++)
    {
        ui->imageLabel->drawPoint( m_splPts[i], 3, color);
    }

    color = qRgb(255, 0, 0);
    size = m_cntPts.size();
    for(size_t i=0; i<size; i++)
    {
        ui->imageLabel->drawPoint( m_cntPts[i], 3, color);
    }
    ui->imageLabel->updateOverlay();
}

void ImageWidget::updateGridPos() const
{
    ui->imageLabel->updateGridPos( m_horizontalSlider->value(), m_verticalSlider->value() );
}

void ImageWidget::updateWindowing(float wwidth, float wcenter) const
{
    ui->imageLabel->updateWindowing(wwidth, wcenter);
}

void ImageWidget::showGrid(bool val) const
{
    ui->imageLabel->showGrid(val);
}

void ImageWidget::clickedAt(const Wm3::Vector3<double>& pt) const
{
    if( pt.X() < 0.0 || pt.Y() < 0.0 || pt.Z() < 0.0 )
        return;
    if( pt.X() > m_data->getLen(Data::COL) || pt.Y() > m_data->getLen(Data::LIN) || pt.Z() > m_data->getLen(Data::PAR) )
        return;

    Wm3::Vector3<double> pos = pt * m_rotMat;
    emit clicked(pos);
}

void ImageWidget::changeOffset(const Wm3::Vector3<double>& offset) const
{
    Wm3::Vector3<double> pos = offset * m_rotMat;
    ui->imageLabel->setOffset(pos);
    updateOverlay();
}

void ImageWidget::sendOffset(const Wm3::Vector3<double>& offset) const
{
    updateOverlay();
    Wm3::Vector3<double> pos = offset * m_rotMat;
    emit offsetChanged(pos);
}

void ImageWidget::changeScale(double scale) const
{
    ui->imageLabel->setScale(scale);
    updateOverlay();
}

void ImageWidget::sendScale(double scale) const
{
    updateOverlay();
    emit scaleChanged(scale);
}

