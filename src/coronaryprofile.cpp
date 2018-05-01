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

#include "coronaryprofile.h"
#include "ui_coronaryprofile.h"

#include "Data.h"
#include "measurementpoint.h"
#include "settings.h"
#include "os.h"

#include <QMouseEvent>

#include <cstring>


CoronaryProfile::CoronaryProfile(QWidget *parent) :
	m_settings(Settings::getSettings()),
    QWidget(parent),
    ui(new Ui::CoronaryProfile)
{
    ui->setupUi(this);

    ui->imageLabel->setBackgroundRole(QPalette::Base);
    ui->imageLabel->setScaledContents(true);
    ui->imageLabel->setMode( ImageLabel::DisplayOnly );

	ui->profileLabel->setMode(ImageLabel::Segmentation);
	ui->profileLabel->showGrid(true);

    m_data          = 0;
    m_overlayImage  = 0;
    m_profile       = 0;
	m_lengthProfile = 0;
	m_refineMode = false;

    connect(ui->slider, SIGNAL(valueChanged(int)), this, SLOT(updateProfile(int)) );
	connect(ui->profileLabel, SIGNAL(clickedAt(Wm3::Vector3<double>)), this, SLOT(lengthProfileClicked(Wm3::Vector3<double>)));
	connect(ui->imageLabel, SIGNAL(clickedAt(Wm3::Vector3<double>)), this, SLOT(crossProfileClicked(Wm3::Vector3<double>)));
}

CoronaryProfile::~CoronaryProfile()
{
    delete ui;
}

void CoronaryProfile::reset()
{
	m_data = NULL;
	
	ui->imageLabel->reset();

	if (m_profile) {
		delete[] m_profile->data0();
		delete m_profile;
	}
	m_profile = NULL;

	ui->profileLabel->reset();

	if (m_lengthProfile) {
		delete m_lengthProfile;
		m_lengthProfile = NULL;
	}

	if (m_overlayImage)
		// This is deleted by imageLabel->reset()
		m_overlayImage = NULL;

	m_measurements.clear();

	update(true);
}

void CoronaryProfile::setData(Data* data)
{
    m_data = data;
}

void CoronaryProfile::setProfileSize(size_t size)
{
    m_length = size;

    if(m_profile)
    {
        delete [] m_profile->data0();
        delete m_profile;
    }

    size = static_cast<size_t>(ui->imageLabel->width());

    m_profile = new Data(new float[size*size]);
    m_profile->create(Data::COL, size, Data::LIN, size);
    m_profile->preset();
    ui->imageLabel->setInputData(m_profile);

    if( m_overlayImage )
        delete m_overlayImage;

    m_overlayImage = new QImage(size, size, QImage::Format_RGB32);
    ui->imageLabel->setOverlayImage(m_overlayImage);
}

void CoronaryProfile::setMeasurements(const std::vector<MeasurementPoint*>& m)
{
	m_measurements = m;

	update(true);
}

void CoronaryProfile::setMeasurements(const std::vector<MeasurementPoint*>& m, const cimg_library::CImg<float>& cImage, bool recache)
{
    m_measurements = m;

	ui->profileLabel->reset();
	if (m_lengthProfile) {
		delete m_lengthProfile;
		m_lengthProfile = NULL;
	}

	m_profileImg = cImage;

	float fac = static_cast<float>(ui->imageLabel->height()) / static_cast<float>(cImage.height());
	m_profileImg.resize(2 * cImage.width(), fac * cImage.height(), 1, 1, 3);

	int width  = m_profileImg.width();
	int height = m_profileImg.height();
	m_lengthProfile = new Data(m_profileImg.data());
	m_lengthProfile->create(Data::COL, width, Data::LIN, height);
	m_lengthProfile->setVoxelSize(m_data->getVoxelSize());
	
	ui->profileLabel->setInputData(m_lengthProfile);
	
	float maxVal;
	float minVal = cImage.min_max(maxVal);
	float wcenter = (maxVal + minVal) / 2.0;
	ui->profileLabel->updateWindowing(maxVal - minVal, wcenter);

	ui->profileLabel->setMinimumWidth(width);
	ui->profileLabel->setMinimumHeight(ui->imageLabel->height());
	ui->profileLabel->resize(width, ui->imageLabel->height());

	this->setMinimumWidth(ui->imageLabel->width() + ui->profileLabel->width());

	update(recache);
}

void CoronaryProfile::setRefineMode(bool on)
{
	if (on)
		ui->imageLabel->setMode(ImageLabel::Segmentation);
	else
		ui->imageLabel->setMode(ImageLabel::DisplayOnly);

	m_refineMode = on;
}

void CoronaryProfile::updateProfile(int pos)
{
    m_overlayImage->fill(0);
    if(pos >= static_cast<int>(m_measurements.size()))
        return;

    ui->slider->setValue(pos);
    m_pos = static_cast<size_t>(pos);
    drawPlane();
    drawOverlay();

    MeasurementPoint* p = m_measurements.at(m_pos);
    ui->l_splPos->setText( QString("%1 mm").arg(p->getPositionInSpline(), 0, 'f', 2, '0') );

	ui->profileLabel->updateGridPos(2 * pos, 0);

    emit showSlice( p->getPositionInVolume() );
    emit positionChanged(m_pos);
}

void CoronaryProfile::lengthProfileClicked(Wm3::Vector3<double> pt)
{
	if (pt.X() < 0)
		return;
	if (pt.X() >= 2 * m_measurements.size())
		return;

	updateProfile(pt.X() / 2.0);
}

void CoronaryProfile::crossProfileClicked(Wm3::Vector3<double> pt)
{
	if (pt.X() < 0)
		return;
	if (pt.Y() < 0)
		return;

	if (pt.X() >= ui->imageLabel->width())
		return;
	if (pt.Y() >= ui->imageLabel->height())
		return;

	if (m_length <= 0)
		return;

	if (m_pos >= m_measurements.size())
		return;

	Wm3::Vector3<double> centerPt(static_cast<double>(ui->imageLabel->width()) / 2.0,
								  static_cast<double>(ui->imageLabel->height()) /2.0,
								  0);
	double step = static_cast<double>(ui->imageLabel->width()) / static_cast<double>(m_length);

	Wm3::Vector3<double> p = (pt - centerPt) / step;
		
	MeasurementPoint* mp = m_measurements.at(m_pos);
	
	Wm3::Vector3<double> u, v;
	mp->getPlane(u, v);
	Wm3::Vector3<double> x = p.X() * u + p.Y() * v + mp->getPositionInVolume();
	
	mp->setPositionInVolume(x);

	redrawCurrentPlane();
	drawOverlay();
	emit measurementsEdit(m_pos);
}

void CoronaryProfile::update(bool recache)
{
    if(m_overlayImage)
        m_overlayImage->fill(0);

    if( m_measurements.empty() )
    {
		ui->slider->setDisabled(true);
        ui->l_splPos->setText( QString("%1 mm").arg(0, 0, 'f', 2, '0') );
        drawOverlay();

        if(m_profile)
            m_profile->preset();
        ui->imageLabel->updateBaseImage();

		if (m_lengthProfile)
			m_lengthProfile->preset();
		ui->profileLabel->updateBaseImage();
		ui->profileLabel->updateGridPos(0, 0);
		ui->profileLabel->clear();

		m_planeCache.resize(0);
		m_scaleCache.resize(0);

        return;
    }

    m_pos = m_measurements.size() - 1;
	
	ui->slider->setDisabled(false);
    ui->slider->setMaximum(m_pos);
    ui->slider->setMinimum( 0 );
    ui->slider->setValue( ui->slider->maximum() );

	if (recache)
		drawAllPlanes();

	drawPlane();
    drawOverlay();
}

void CoronaryProfile::drawAllPlanes()
{
	m_planeCache.resize(m_measurements.size());
	m_scaleCache.resize(m_measurements.size());

	for (size_t pos = 0; pos < m_measurements.size(); ++pos)
	{
		cimg_library::CImg<float>& img = m_planeCache[pos];
		MeasurementPoint* mp = m_measurements[pos];
		Wm3::Vector3<double> p = mp->getPositionInVolume();
		Wm3::Vector3<double> u, v;
		mp->getPlane(u, v);

		float minimum = m_data->getValue(p.X(), p.Y(), p.Z());
		float maximum = minimum;

		double step = static_cast<double>(m_length) / static_cast<double>(ui->imageLabel->width());

		double hwidth = static_cast<double>(ui->imageLabel->width()) / 2.0;
		
		size_t xSize = m_profile->getLen(Data::COL);
		size_t ySize = m_profile->getLen(Data::LIN);

		img.resize(xSize, ySize);
		img.fill(0);

		float* dProfile = img.data();

		for (double i = 0; i < ySize; ++i)
		{
			double y = (i - hwidth) * step;
			for (double j = 0; j < xSize; ++j)
			{
				double x = (j - hwidth) * step;
				Wm3::Vector3<double> curPt = p + x * u + y * v;
				float val = m_data->getValue(curPt);
				minimum = std::min<double>(val, minimum);
				maximum = std::max<double>(val, maximum);
				dProfile[static_cast<size_t>(round(i * static_cast<double>(xSize) + j))] = val;
			}
		}

		float scale = maximum - minimum;
		m_scaleCache[pos] = std::pair<float, float>(scale, minimum + scale / 2.0);
	}
}

void CoronaryProfile::redrawCurrentPlane()
{
	if (m_pos >= m_measurements.size())
		return;
	
	if ((m_planeCache.size() != m_measurements.size()) || (m_scaleCache.size() != m_measurements.size())) {
		drawAllPlanes();
		return;
	}

	cimg_library::CImg<float>& img = m_planeCache[m_pos];
	MeasurementPoint* mp = m_measurements[m_pos];
	Wm3::Vector3<double> p = mp->getPositionInVolume();
	Wm3::Vector3<double> u, v;
	mp->getPlane(u, v);

	float minimum = m_data->getValue(p.X(), p.Y(), p.Z());
	float maximum = minimum;

	double step = static_cast<double>(m_length) / static_cast<double>(ui->imageLabel->width());

	double hwidth = static_cast<double>(ui->imageLabel->width()) / 2.0;
		
	size_t xSize = m_profile->getLen(Data::COL);
	size_t ySize = m_profile->getLen(Data::LIN);

	img.resize(xSize, ySize);
	img.fill(0);

	float* dProfile = img.data();

	for (double i = 0; i < ySize; ++i)
	{
		double y = (i - hwidth) * step;
		for (double j = 0; j < xSize; ++j)
		{
			double x = (j - hwidth) * step;
			Wm3::Vector3<double> curPt = p + x * u + y * v;
			float val = m_data->getValue(curPt);
			minimum = std::min<double>(val, minimum);
			maximum = std::max<double>(val, maximum);
			dProfile[static_cast<size_t>(round(i * static_cast<double>(xSize) + j))] = val;
		}
	}

	float scale = maximum - minimum;
	m_scaleCache[m_pos] = std::pair<float, float>(scale, minimum + scale / 2.0);

	drawPlane();
}

void CoronaryProfile::drawPlane()
{
	if (m_planeCache.size() != m_measurements.size())
		drawAllPlanes();

	const float* img = m_planeCache[m_pos].data();
	float* dProfile = m_profile->data0();
	std::memcpy(dProfile, img, sizeof(float) * m_planeCache[m_pos].size());

	ui->imageLabel->updateWindowing(m_scaleCache[m_pos].first, m_scaleCache[m_pos].second);
}

inline void drawCross( QImage* image, QRgb color, const Wm3::Vector3<double>& p, int crossSize )
{
    if (p.X() < 0 || p.X() > image->width() ||
        p.Y() < 0 || p.Y() > image->height() )
        return;

    double hCrossSize = (double)crossSize/2.0;
    int x1, x2;
    x1 = round( p.X() - hCrossSize);
    x2 = x1 + crossSize;
    for(int x=x1; x<x2; x++)
        image->setPixel(x, p.Y(), color);

    x1 = round( p.Y() - hCrossSize );
    x2 = x1 + crossSize;
    for(int y=x1; y<x2; y++)
        image->setPixel(p.X(), y, color);
}

void CoronaryProfile::drawOverlay()
{
    if( m_overlayImage )
        m_overlayImage->fill(0);

    if( m_pos >= m_measurements.size() )
    {
        ui->imageLabel->updateBaseImage();
        ui->imageLabel->updateOverlay();
        return;
    }

    Wm3::Vector3<double> centerPt((double)ui->imageLabel->width()/2.0, (double)ui->imageLabel->height()/2.0, 0 );
    drawCross(m_overlayImage, qRgb(0, 255, 255), centerPt, 9);

    MeasurementPoint* mp = m_measurements.at(m_pos);
    double step = (double)ui->imageLabel->width() / (double)m_length;

    std::vector< Wm3::Vector3<double> > leftPts, rightPts;
    size_t numNormals = m_settings.getNumNormalsAtPoint();
    for(size_t j=0; j<numNormals; j++)
    {
        Wm3::Vector3<double> p1 = mp->get50PointInPlane(j, Interval::Left) * step + centerPt;
        drawCross(m_overlayImage, qRgb(0, 255, 0), p1, 3);
        leftPts.push_back(p1);

        Wm3::Vector3<double> p2 = mp->get50PointInPlane(j, Interval::Right) * step + centerPt;
        drawCross(m_overlayImage, qRgb(0, 255, 0), p2, 3);
        rightPts.push_back(p2);
    }

    Wm3::Vector3<double> newCenter = mp->getImprovedPositionInPlane() * step + centerPt;
    drawCross(m_overlayImage, qRgb(255, 0, 0), newCenter, 9);

    double radius = 0.5 * mp->getDiameter() * step / m_data->getVoxelSize();
    Wm3::Vector3<double> u,v;
    mp->getPlane(u,v);

    for(size_t i=0; i<rightPts.size(); i++)
    {
        Wm3::Vector3<double> sd = mp->getSamplingDirection(i);
        Wm3::Vector3<double> d =  Wm3::Vector3<double>(sd.Dot(u) / u.Dot(u), sd.Dot(v) / v.Dot(v), 0);

        Wm3::Vector3<double> lRadius = -1.0*radius*d + newCenter;
        Wm3::Vector3<double> rRadius =      radius*d + newCenter;

        if( !mp->isPointOk(i, Interval::Right) )
            drawCross(m_overlayImage, qRgb(255, 0, 0), rightPts.at(i), 3);

        if( !mp->isPointOk(i, Interval::Left) )
            drawCross(m_overlayImage, qRgb(255, 0, 0), leftPts.at(i), 3);

        drawCross(m_overlayImage, qRgb(255, 255, 255), lRadius, 3);
        drawCross(m_overlayImage, qRgb(255, 255, 255), rRadius, 3);
    }

    ui->imageLabel->updateOverlay();

    return;
}

