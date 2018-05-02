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

#include "orthoviewer.h"
#include "ui_orthoWidget.h"

#include "Data.h"
#include "imagewidget.h"
#include "pointwidget.h"
#include "bsplineinterpolation.h"
#include "mipdialog.h"
#include "settings.h"

#include <QDebug>
#include <QTime>


OrthoViewer::OrthoViewer(QWidget* parent) :
    QWidget(parent),
    ui(new Ui::OrthoWidget)
{
    ui->setupUi(this);

    m_bSpline = 0;

	if (Settings::getSettings().getShowMIP()) {
		m_mipDialog = new MipDialog(this);
		connect(this, SIGNAL(setWindowing(float,float)),  m_mipDialog, SLOT(updateWindowing(float,float)));
	} else
		m_mipDialog = NULL;

    connect(ui->ww_slider, SIGNAL(sliderMoved(int)), this, SLOT( wWidthChanged(int)));
    connect(ui->wc_slider, SIGNAL(sliderMoved(int)), this, SLOT( wCenterChanged(int)));

    connect(ui->b_seg, SIGNAL(clicked()), this, SLOT(seg_clicked()));
    connect(ui->b_move, SIGNAL(clicked()), this, SLOT(move_clicked()));

    ImageWidget* coronalImg  = ui->w_coronal;
    ImageWidget* sagittalImg = ui->w_sagittal;
    ImageWidget* axialImg    = ui->w_axial;
    connect(ui->grid,                           SIGNAL(toggled(bool)),              coronalImg,                         SLOT(showGrid(bool)));
    connect(this,                               SIGNAL(setWindowing(float,float)),  coronalImg,                         SLOT(updateWindowing(float,float)));
    connect(coronalImg,                         SIGNAL(sliceChanged(int)),          axialImg->getHorizontalSlider(),    SLOT(setValue(int)));
    connect(coronalImg,                         SIGNAL(sliceChanged(int)),          sagittalImg->getHorizontalSlider(), SLOT(setValue(int)));

	connect(ui->grid,                           SIGNAL(toggled(bool)),              axialImg,                           SLOT(showGrid(bool)));
    connect(this,                               SIGNAL(setWindowing(float,float)),  axialImg,                           SLOT(updateWindowing(float,float)));
    connect(axialImg,                           SIGNAL(sliceChanged(int)),          sagittalImg->getVerticalSlider(),   SLOT(setValue(int)));
    connect(axialImg,                           SIGNAL(sliceChanged(int)),          coronalImg->getVerticalSlider(),    SLOT(setValue(int)));

    connect(ui->grid,                           SIGNAL(toggled(bool)),              sagittalImg,                        SLOT(showGrid(bool)));
    connect(this,                               SIGNAL(setWindowing(float,float)),  sagittalImg,                        SLOT(updateWindowing(float,float)));
    connect(sagittalImg,                        SIGNAL(sliceChanged(int)),          axialImg->getVerticalSlider(),      SLOT(setValue(int)));
    connect(sagittalImg,                        SIGNAL(sliceChanged(int)),          coronalImg->getHorizontalSlider(),  SLOT(setValue(int)));

    connect(ui->wc_slider,                      SIGNAL(valueChanged(int)),          this,                               SLOT(wCenterChanged(int)) );
    connect(ui->ww_slider,                      SIGNAL(valueChanged(int)),          this,                               SLOT(wWidthChanged(int)) );
    connect(ui->autoWindow,                     SIGNAL(toggled(bool)),              this,                               SLOT(autoWindowingChanged(bool)));

    // Receive Clicked Control-Points and send ControlPt-List / SplinePt-List to Overlay
    connect(coronalImg,  SIGNAL(clicked(const Wm3::Vector3<double>&)), this, SLOT(addPoint(const Wm3::Vector3<double>&)));
    connect(axialImg,    SIGNAL(clicked(const Wm3::Vector3<double>&)), this, SLOT(addPoint(const Wm3::Vector3<double>&)));
    connect(sagittalImg, SIGNAL(clicked(const Wm3::Vector3<double>&)), this, SLOT(addPoint(const Wm3::Vector3<double>&)));

    connect(this, SIGNAL(sendControlPoints(const std::vector< Wm3::Vector3<double> >&)), coronalImg,  SLOT(updateControlPoints(const std::vector< Wm3::Vector3<double> >&)));
    connect(this, SIGNAL(sendControlPoints(const std::vector< Wm3::Vector3<double> >&)), axialImg,    SLOT(updateControlPoints(const std::vector< Wm3::Vector3<double> >&)));
    connect(this, SIGNAL(sendControlPoints(const std::vector< Wm3::Vector3<double> >&)), sagittalImg, SLOT(updateControlPoints(const std::vector< Wm3::Vector3<double> >&)));

    connect(this, SIGNAL(sendSplinePoints(const std::vector< Wm3::Vector3<double> >&)), coronalImg,  SLOT(updateSplinePoints(const std::vector< Wm3::Vector3<double> >&)));
    connect(this, SIGNAL(sendSplinePoints(const std::vector< Wm3::Vector3<double> >&)), axialImg,    SLOT(updateSplinePoints(const std::vector< Wm3::Vector3<double> >&)));
    connect(this, SIGNAL(sendSplinePoints(const std::vector< Wm3::Vector3<double> >&)), sagittalImg, SLOT(updateSplinePoints(const std::vector< Wm3::Vector3<double> >&)));

    connect(this, SIGNAL(changeSlice(const Wm3::Vector3<double>&)), coronalImg, SLOT(changeSlice(const Wm3::Vector3<double>&)));
    connect(this, SIGNAL(changeSlice(const Wm3::Vector3<double>&)), axialImg, SLOT(changeSlice(const Wm3::Vector3<double>&)));
    connect(this, SIGNAL(changeSlice(const Wm3::Vector3<double>&)), sagittalImg, SLOT(changeSlice(const Wm3::Vector3<double>&)));


    connect(coronalImg, SIGNAL(scaleChanged(double)), sagittalImg, SLOT(changeScale(double)) );
    connect(coronalImg, SIGNAL(scaleChanged(double)), axialImg,    SLOT(changeScale(double)) );
    connect(sagittalImg, SIGNAL(scaleChanged(double)), coronalImg,    SLOT(changeScale(double)) );
    connect(sagittalImg, SIGNAL(scaleChanged(double)), axialImg,    SLOT(changeScale(double)) );
    connect(axialImg, SIGNAL(scaleChanged(double)), coronalImg,    SLOT(changeScale(double)) );
    connect(axialImg, SIGNAL(scaleChanged(double)), sagittalImg,    SLOT(changeScale(double)) );

    connect(coronalImg, SIGNAL(offsetChanged(const Wm3::Vector3<double>&)), sagittalImg, SLOT(changeOffset(const Wm3::Vector3<double>&)) );
    connect(coronalImg, SIGNAL(offsetChanged(const Wm3::Vector3<double>&)), axialImg,    SLOT(changeOffset(const Wm3::Vector3<double>&)) );
    connect(sagittalImg, SIGNAL(offsetChanged(const Wm3::Vector3<double>&)), coronalImg,    SLOT(changeOffset(const Wm3::Vector3<double>&)) );
    connect(sagittalImg, SIGNAL(offsetChanged(const Wm3::Vector3<double>&)), axialImg,    SLOT(changeOffset(const Wm3::Vector3<double>&)) );
    connect(axialImg, SIGNAL(offsetChanged(const Wm3::Vector3<double>&)), coronalImg,    SLOT(changeOffset(const Wm3::Vector3<double>&)) );
    connect(axialImg, SIGNAL(offsetChanged(const Wm3::Vector3<double>&)), sagittalImg,    SLOT(changeOffset(const Wm3::Vector3<double>&)) );

    ui->w_coronal->setStyleSheet("background-color: black;");
    ui->w_sagittal->setStyleSheet("background-color: black;");
    ui->w_axial->setStyleSheet("background-color: black;");

    m_pointWidget = new PointWidget(this);
    
    connect(m_pointWidget, SIGNAL( pointsChanged() ), this, SLOT( update()));

    move_clicked();
}

OrthoViewer::~OrthoViewer()
{
    delete ui;
}

void OrthoViewer::reset()
{
	m_pointWidget->reset();
	m_bSpline = NULL;
	ui->w_sagittal->reset();
	ui->w_axial->reset();
	ui->w_coronal->reset();
	m_data = NULL;
}

void OrthoViewer::setData(Data* data)
{
    m_data = data;

    if(m_bSpline)
		m_bSpline->reset();

    m_pointWidget->setData(data);

    float* dData = data->data0();
    m_minimum = dData[0];
    m_maximum = m_minimum;
    size_t size = data->size();
    for(size_t i=0; i<size; i++)
    {
        m_maximum = std::max<double>(m_maximum, dData[i]);
        m_minimum = std::min<double>(m_minimum, dData[i]);
    }

    qDebug() << "Prepare SagittalData";
    Wm3::Matrix3<double> matId = Wm3::Matrix3<double>::IDENTITY;
    ui->w_sagittal->setInputData(data, matId);

    qDebug() << "Prepare AxialData";
    Wm3::Matrix3<double> matRx(1,0,0, 0,0,1, 0,1,0);
    ui->w_axial->setInputData(data, matRx);

    qDebug() << "Prepare CoronalData";
    Wm3::Matrix3<double> matRy(0,0,1, 0,1,0, 1,0,0);
    ui->w_coronal->setInputData(data, matRy);

	// AutoScaling
    if( ui->autoWindow->isChecked() )
    {
        autoWindowingChanged(true);
    }

    showSlice( Wm3::Vector3<double>(data->getLen(Data::COL)/2.0,data->getLen(Data::LIN)/2.0,data->getLen(Data::PAR)/2.0) );

    this->setMinimumWidth( ui->w_coronal->width() + ui->w_sagittal->width() );
    this->setMaximumWidth( ui->w_coronal->width() + ui->w_sagittal->width() );
    emit setWindowing(m_wwidth, m_wcenter);
}

void OrthoViewer::setInterpolator(BSplineInterpolation* bSpline)
{
    m_bSpline = bSpline;
    m_pointWidget->setInterpolator(m_bSpline);
}

void OrthoViewer::wCenterChanged(int val)
{
   m_wcenter = (float)val / (float)ui->wc_slider->maximum() * (m_maximum-m_minimum);
   emit setWindowing(m_wwidth, m_wcenter);
}

void OrthoViewer::wWidthChanged(int val)
{
   m_wwidth = (float)val / (float)ui->ww_slider->maximum() * (m_maximum-m_minimum);
   emit setWindowing(m_wwidth, m_wcenter);
}

void OrthoViewer::autoWindowingChanged(bool val)
{
    if(val)
    {
        m_wwidth  = (m_maximum-m_minimum);
        m_wcenter = m_minimum + m_wwidth / 2.0;
        ui->ww_slider->setValue( ui->ww_slider->maximum() );
        ui->wc_slider->setValue( ui->wc_slider->maximum() / 2 );

        ui->ww_slider->setDisabled(true);
        ui->wc_slider->setDisabled(true);

        emit setWindowing(m_wwidth, m_wcenter);
    } else {
        ui->ww_slider->setDisabled(false);
        ui->wc_slider->setDisabled(false);
    }
}

void OrthoViewer::addPoint(const Wm3::Vector3<double>& pt)
{
    m_pointWidget->addControlPoint(pt);
    emit changeSlice(pt);

    update();
}

void OrthoViewer::update()
{
    if(!m_bSpline)
        return;

    emit sendControlPoints( m_pointWidget->getControlPoints() );
    m_bSpline->setSamples( m_pointWidget->getControlPoints() );

    const std::vector< Wm3::Vector3<double> >& splPts = m_bSpline->getInterpolatedPoints();

	if (Settings::getSettings().getShowMIP()) {
		if (m_mipDialog == NULL) {
			m_mipDialog = new MipDialog(this);
			connect(this, SIGNAL(setWindowing(float,float)),  m_mipDialog, SLOT(updateWindowing(float,float)));
		}

		m_mipDialog->setData(m_bSpline, m_data);
		m_mipDialog->show();
	} else if (m_mipDialog != NULL) {
			m_mipDialog->close();

			disconnect(m_mipDialog, SLOT(updateWindowing(float, float)));
			delete m_mipDialog;
			m_mipDialog = NULL;
	}

    emit sendSplinePoints( splPts );
    emit updateProfile();
}

void OrthoViewer::showSlice(const Wm3::Vector3<double>& pt)
{
    emit changeSlice(pt);
}

void OrthoViewer::seg_clicked() const
{
    ui->b_seg->setDown(true);
    ui->b_move->setDown(false);
    ui->w_coronal->setMode(ImageLabel::Segmentation);
    ui->w_sagittal->setMode(ImageLabel::Segmentation);
    ui->w_axial->setMode(ImageLabel::Segmentation);
}

void OrthoViewer::move_clicked() const
{
    ui->b_seg->setDown(false);
    ui->b_move->setDown(true);
    ui->w_coronal->setMode(ImageLabel::Resizing);
    ui->w_sagittal->setMode(ImageLabel::Resizing);
    ui->w_axial->setMode(ImageLabel::Resizing);
}

void OrthoViewer::showPointsDialog() {
	m_pointWidget->show();
	m_pointWidget->activateWindow();
	m_pointWidget->raise();
}

