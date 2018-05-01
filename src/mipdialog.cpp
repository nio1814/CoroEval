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

#include "mipdialog.h"

#include <QSlider>
#include <QVBoxLayout>
#include <QTime>
#include <QLabel>

#include "bsplineinterpolation.h"
#include "Data.h"
#include "measurementpoint.h"


MipDialog::MipDialog(QWidget* parent):
    QDialog(parent)
{
    m_imageLabel = new QLabel();
    m_slider     = new QSlider( Qt::Horizontal );
    m_img        = 0;

    m_slider->setMinimum(0);
    m_slider->setMaximum(360);

    QVBoxLayout* layout = new QVBoxLayout(this);
    layout->addWidget(m_imageLabel);
    layout->addWidget(m_slider);

    setStyleSheet("background-color: black;");

    connect(m_slider, SIGNAL(valueChanged(int)), this, SLOT(updateAlpha(int)));
}

void MipDialog::setData(BSplineInterpolation* bSpline, Data* data)
{
    m_data    = data;
    m_bSpline = bSpline;

    if( m_img )
		delete m_img;

    int width  = m_data->getLen(Data::COL);
    int height = m_data->getLen(Data::LIN);
    m_img     = new cimg_library::CImg<float>(width, height);

    renderMip();
}

void MipDialog::updateAlpha(int val)
{
    m_alpha = double(val)/180.0 * Wm3::Mathd::PI;
    renderMip();
}

void MipDialog::updateWindowing(float wwidth, float wcenter)
{
    m_wwidth  = wwidth;
    m_wcenter = wcenter;
    renderMip();
}

void MipDialog::renderMip()
{
    if(!m_data || !m_bSpline || m_bSpline->getLength() < 0.0 )
        return;

    QTime t;
    t.start();

    int width  = m_data->getLen(Data::COL);
    int height = m_data->getLen(Data::LIN);

    Wm3::Vector3<double> vecOrigin( 0.5*m_data->getLen(Data::COL),
                                    0.5*m_data->getLen(Data::LIN),
                                    0.5*m_data->getLen(Data::PAR));
    Wm3::Quaternion<double> q(Wm3::Vector3<double>(0,1,0), m_alpha);

    m_img->fill(0);

    double currPos = 0.0;
    double incr    = 0.01;
    double thres   = 1.2;
    while(currPos < 0.99)
    {

        double l = (m_bSpline->pointAt(currPos) - m_bSpline->pointAt(currPos+incr)).Length();
        if( l < thres )
        {
            while(l < thres)
            {
                incr *= 2.0;
                l = (m_bSpline->pointAt(currPos) - m_bSpline->pointAt(currPos+incr)).Length();
            }
        } else {
            while( l > thres )
            {
                incr /= 2.0;
                l = (m_bSpline->pointAt(currPos) - m_bSpline->pointAt(currPos+incr)).Length();
            }
        }
        currPos += incr;

        MeasurementPoint mp(currPos, m_bSpline, m_data);
        const Wm3::Vector3<double>& pt = mp.getPositionInVolume();

        double radius = mp.getDiameter() * 0.5; // radius + 10% extra
        double wradius= radius * 0.75;
        double minX = pt.X() - radius - 1.;
        double maxX = pt.X() + radius + 1.;
        double minY = pt.Y() - radius - 1.;
        double maxY = pt.Y() + radius + 1.;
        double minZ = pt.Z() - radius - 1.;
        double maxZ = pt.Z() + radius + 1.;

        for(double z=minZ; z<maxZ; z += 1.0)
        {
            for(double y=minY; y<maxY; y += 1.0)
            {
                for(double x=minX; x<maxX; x += 1.0)
                {
                    Wm3::Vector3<double> pos(x,y,z);
                    double distance = (pt-pos).Length();
                    if( distance > radius )
                         continue;

                    Wm3::Vector3<double> imgPos = q.Rotate( pos - vecOrigin ) + vecOrigin;
                    float val    = m_img->linear_atXY(imgPos.X(), imgPos.Y() );
                    float newVal = computeCutOff(distance, wradius) * m_data->getValue(x,y,z);

                    if( newVal > val )
                        m_img->set_linear_atXY(newVal, imgPos.X(), imgPos.Y());
                }
            }
        }
    }

    float lowerBound =  m_wcenter - m_wwidth/2.0;
    float upperBound =  m_wcenter + m_wwidth/2.0;
    QImage img( width, height, QImage::Format_RGB32 );
    for(int y=0; y<height; y++)
    {
        for(int x=0; x<width; x++)
        {
            float fval = m_img->atXY(x,y);
            if(fval <= lowerBound ) {
                img.setPixel(x, y, qRgb(0,0,0));
            } else if ( fval >= upperBound ) {
                img.setPixel(x, y, qRgb(255,255,255));
            } else {
                int val = ( int ) ( ( fval - lowerBound ) / ( upperBound - lowerBound ) * 255.0 );
                img.setPixel(x, y, qRgb(val,val,val));
            }
        }
    }
    m_imageLabel->setPixmap( QPixmap::fromImage(img) );
}

