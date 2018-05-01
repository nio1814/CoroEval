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

#ifndef IMAGELABEL_H
#define IMAGELABEL_H

#include <QPainter>
#include <QLabel>
#include <QMouseEvent>

#include "Wm3Vector3.h"

#include <iostream>


// Forward declarations
class Data;
class QImage;


class ImageLabel : public QLabel
{
    Q_OBJECT
public:
    ImageLabel(QWidget* parent = 0);
    ~ImageLabel();

    enum Mode{
        DisplayOnly,
        Segmentation,
        Resizing
    };

	void    reset();

    void    setInputData(Data* data);
    void    updateWindowing(float wwidth, float wcenter);
    void    updateGridPos(int horizontal, int vertical);
    void    showGrid(bool val);

    void    updateOverlay();
    void    updateBaseImage();
    void    setOverlayImage(QImage* image);

    void    clearOverlayImage();
    void    drawPoint(Wm3::Vector3<double> pt, size_t size, QRgb color);
    void    drawCross(Wm3::Vector3<double> pt, size_t size, QRgb color);

    void    setMode(Mode m);
    void    setScale(double val);
    void    setOffset(const Wm3::Vector3<double>& pt);

private:
    Data*   m_data;
    QImage* m_baseImage;
    QImage* m_overlayImage;

    bool    m_showGrid;
    float   m_wwidth;
    float   m_wcenter;

    int     m_vertical;
    int     m_horizontal;

    Mode    m_mode;
    Qt::MouseButton m_button;
    QPoint  m_prevMousePos;
    float   m_scale;
    QPoint  m_offset;

    int     calculateGreyValue (float data) const;
    void    drawGrid(QImage& image) const;
    void    drawOverlay(QImage& image) const;
    void    drawImage();
    void    updateCursor();

signals:
    void    clickedAt(const Wm3::Vector3<double>&) const;
    void    sliceChanged(double) const;
    void    scaleChanged(double) const;
    void    offsetChanged(const Wm3::Vector3<double>&) const;

public slots:

protected:
    void mousePressEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);
    void resizeEvent(QResizeEvent* event);
};

#endif // IMAGELABEL_H

