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

#include "imagelabel.h"

#include "Data.h"
#include "os.h"

#include <QImage>


ImageLabel::ImageLabel(QWidget* parent) :
    QLabel(parent)
{
    m_data          = 0;
    m_baseImage     = 0;
    m_overlayImage  = 0;
    m_showGrid      = false;

    m_wwidth        = 0;
    m_wcenter       = 0;

    m_mode          = ImageLabel::Resizing;
    m_offset        = QPoint(0,0);
    m_scale         = 1.0;

    this->adjustSize();

    updateCursor();

}

ImageLabel::~ImageLabel()
{
    if(m_baseImage)
		delete m_baseImage;
    if(m_overlayImage)
		delete m_overlayImage;
}

void ImageLabel::reset()
{
	if(m_baseImage)
		delete m_baseImage;
    if(m_overlayImage)
		delete m_overlayImage;

	m_baseImage = NULL;
	m_overlayImage = NULL;
	m_data = NULL;

	clear();
}

void ImageLabel::setInputData(Data* img)
{
    m_data   = img;

	if (m_overlayImage)
		delete m_overlayImage;
    if( m_baseImage )
		delete m_baseImage;

    int col = m_data->getLen(Data::COL);
    int lin = m_data->getLen(Data::LIN);
    m_baseImage     = new QImage(col, lin, QImage::Format_ARGB32_Premultiplied);
    m_overlayImage  = new QImage(col, lin, QImage::Format_ARGB32_Premultiplied);

    clearOverlayImage();

    updateBaseImage();
    updateCursor();
}

void ImageLabel::updateWindowing(float wwidth, float wcenter)
{
    m_wwidth  = wwidth;
    m_wcenter = wcenter;
    updateBaseImage();
}

void ImageLabel::updateGridPos(int horizontal, int vertical)
{
    m_horizontal = horizontal;
    m_vertical   = vertical;

    if( m_showGrid )
        drawImage();
}

void ImageLabel::showGrid(bool val)
{
    m_showGrid = val;
    drawImage();
}

void ImageLabel::updateBaseImage()
{
    if( m_baseImage )
        m_baseImage->fill(0);

    if(!m_data)
        return;

    double xOffset = double( m_offset.x() );
    double yOffset = double( m_offset.y() );
    int col = m_data->getLen(Data::COL);
    int lin = m_data->getLen(Data::LIN);
    double z = m_data->getCurrentPar();
    for(double y=0; y<lin; y++)
    {
        for(double x=0; x<col; x++)
        {
            int val = calculateGreyValue( m_data->getValue(xOffset + x*m_scale, yOffset+ y*m_scale, z ) );
            m_baseImage->setPixel(x, y, qRgb(val, val, val));
        }
    }
    drawImage();
}

void ImageLabel::updateOverlay()
{
    drawImage();
    updateCursor();
}

void ImageLabel::setOverlayImage(QImage* image)
{
	if (m_overlayImage)
		delete m_overlayImage;

    m_overlayImage = image;
}

void ImageLabel::clearOverlayImage()
{
    if( !m_overlayImage )
        return;
    m_overlayImage->fill(0);
}

void ImageLabel::drawPoint(Wm3::Vector3<double> pt, size_t size, QRgb color)
{
    if( !m_overlayImage )
        return;

    double z1 = pt.Z()-1.0;
    double z2 = pt.Z()+1.0;
    double z = m_data->getCurrentPar();

    if( !(z < z2 && z > z1) )
        return;

    pt = (pt - Wm3::Vector3<double>( m_offset.x(), m_offset.y(), 0) ) / m_scale;
    double hwidth = floor( double(size) / 2.0 );

    if(pt.X()-hwidth < 0 || pt.Y()-hwidth < 0)
        return;
    if(pt.X()+hwidth > m_overlayImage->width() || pt.Y()+hwidth > m_overlayImage->height())
        return;

    int x1 = round(pt.X()-hwidth);
    int x2 = round(pt.X()+hwidth);
    int y1 = round(pt.Y()-hwidth);
    int y2 = round(pt.Y()+hwidth);
    for(int y=y1; y<y2; y++)
    {
        for(int x=x1; x<x2; x++)
        {
            m_overlayImage->setPixel(x, y, color );
        }
    }
}

void ImageLabel::drawCross(Wm3::Vector3<double> pt, size_t size, QRgb color)
{

}

void ImageLabel::setMode(Mode m)
{
    m_mode = m;
    updateCursor();
}

void ImageLabel::setScale(double val)
{
    m_scale = val;
    updateBaseImage();
}

void ImageLabel::setOffset(const Wm3::Vector3<double>& pt)
{
    m_offset = QPoint(pt.X(), pt.Y());
    updateBaseImage();
}

int ImageLabel::calculateGreyValue(float data) const
{
    float lowerBound =  m_wcenter - m_wwidth/2.0;
    float upperBound =  m_wcenter + m_wwidth/2.0;

    if ( data <= lowerBound )
       return 0;

    if ( data >= upperBound )
       return 255;

    return static_cast<int>( ( data - lowerBound ) / ( upperBound - lowerBound ) * 255.0 );
}

void ImageLabel::drawGrid(QImage& image) const
{
    int col = m_data->getLen(Data::COL);
    int lin = m_data->getLen(Data::LIN);

    double xOffset = double( m_offset.x() );
    double yOffset = double( m_offset.y() );

    int vertical   = double(m_vertical-yOffset) / m_scale;
    int horizontal = double(m_horizontal-xOffset) / m_scale;

    if(vertical > 0 && vertical < lin)
    {
        QRgb value = qRgb(0, 255, 0);
        for(int x=0; x<col; x++)
            image.setPixel(x, vertical, value);
    }

    if(horizontal > 0 && horizontal < col)
    {
        QRgb value = qRgb(255, 0, 0);
        for(int y=0; y<lin; y++)
            image.setPixel(horizontal, y, value);
    }
}

void ImageLabel::drawOverlay(QImage& image) const
{
    if(!m_overlayImage)
        return;

    if(!m_data)
        return;

    int col = m_data->getLen(Data::COL);
    int lin = m_data->getLen(Data::LIN);

    for(int y=0; y<lin; y++)
    {
        for(int x=0; x<col; x++)
        {
            QRgb value = m_overlayImage->pixel(x,y);
            if( QColor(value) != QColor(0,0,0) )
                image.setPixel(x, y, value);
        }
    }
}

void ImageLabel::drawImage()
{
    if(!m_data)
        return;

    int col = m_data->getLen(Data::COL);
    int lin = m_data->getLen(Data::LIN);

    QImage image = m_baseImage->copy(0, 0, col, lin);

    if(m_showGrid)
        drawGrid(image);

    drawOverlay(image);

    this->setPixmap(QPixmap::fromImage(image));
}

void ImageLabel::mousePressEvent (QMouseEvent* event)
{
    m_prevMousePos = event->pos();
    m_button       = event->button();
}

void ImageLabel::mouseMoveEvent(QMouseEvent* ev)
{
    if(!m_data)
        return;

    if( m_mode == ImageLabel::Resizing )
    {
        QPoint delta = m_prevMousePos-ev->pos();
        m_prevMousePos = ev->pos();

        if( m_button == Qt::LeftButton )
        {
            m_offset += delta;

            emit offsetChanged( Wm3::Vector3<double>(m_offset.x(), m_offset.y(), 0) );
        }
        else if( m_button == Qt::RightButton )
        {
            if( abs(delta.y()) > 0)
                m_scale += double(delta.y()) / 1000.0;
            else if( abs(delta.x()) > 0)
                m_scale += double(delta.x()) / 1000.0;

            emit scaleChanged(m_scale);
        }
        updateBaseImage();
    }
}

void ImageLabel::mouseReleaseEvent(QMouseEvent* event)
{
    if(!m_data)
        return;

    if( m_mode == ImageLabel::Segmentation )
    {
        QPoint pos = event->pos();
        double x = (double)pos.x() * m_scale + m_offset.x();
        double y = (double)pos.y() * m_scale + m_offset.y();
        double z = m_data->getCurrentPar();

        Wm3::Vector3<double> wm3Vec(x, y, z);
        emit clickedAt(wm3Vec);
    }
}

void ImageLabel::updateCursor()
{
    if( m_mode == ImageLabel::Segmentation )
    {
        // Change Mouse Cursor to a red cross:
        QImage image(5, 5, QImage::Format_ARGB32_Premultiplied);
        image.fill(qRgba(0, 0, 0, 0));

        QRgb value = qRgb(255, 0, 0);
        for(int x=0; x<5; x++)
            image.setPixel(x, 2, value);
        for(int y=0; y<5; y++)
            image.setPixel(2, y, value);

        setCursor(QPixmap::fromImage(image));
    }
	if ((m_mode == ImageLabel::Resizing ) || (m_mode == ImageLabel::DisplayOnly))
    {
        setCursor( Qt::ArrowCursor );
    }
}

void ImageLabel::wheelEvent (QWheelEvent* event)
{
    if(!m_data)
        return;
    int val = event->delta()/120;
    emit sliceChanged(m_data->getCurrentPar() + double(val)*m_scale );
}

void ImageLabel::resizeEvent(QResizeEvent* event)
{
    QWidget::resizeEvent(event);

    if( !m_data )
        return;

    int col = m_data->getLen(Data::COL);
    int lin = m_data->getLen(Data::LIN);

    QRect rec = this->geometry();
    double ratioHeight = (double)rec.height() / lin;
    double ratioWidth  = (double)rec.width()  / col;

    double scale = ratioWidth;
    if(scale > ratioHeight)
        scale = ratioHeight;

    int newHeight = (int) (scale * lin );
    int newWidth  = (int) (scale * col );

    this->setMinimumSize(newWidth, newHeight);
    this->setMaximumSize(newWidth, newHeight);
}

