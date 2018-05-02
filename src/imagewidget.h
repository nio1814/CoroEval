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

#ifndef IMAGEWIDGET_H
#define IMAGEWIDGET_H

#include <QWidget>

#include "Wm3Vector3.h"
#include "Wm3Matrix3.h"

#include "imagelabel.h"


// Forward declarations
class QSlider;
class Data;

namespace Ui {
    class ImageWidget;
}


class ImageWidget : public QWidget
{
    Q_OBJECT

public:
    explicit ImageWidget(QWidget* parent = 0);
    ~ImageWidget();

    void setInputData(Data* img, const Wm3::Matrix3<double>& rotMat);
	void reset();

    QSlider* getHorizontalSlider() const;
    QSlider* getVerticalSlider() const;

    void setMode(ImageLabel::Mode m);

public slots:
    void updateOverlay() const;
    void updateGridPos() const;
    void updateWindowing(float wwidth, float wcenter) const;
    void showGrid(bool val) const;

    void clickedAt(const Wm3::Vector3<double>& pt) const;

private slots:
    void changeSlice(int val);
    void changeSlice(double val);
    void changeSlice(const Wm3::Vector3<double>& pt);
    void updateControlPoints(const std::vector< Wm3::Vector3<double> >& pts);
    void updateSplinePoints(const std::vector< Wm3::Vector3<double> >& pts);

    void sendScale(double val) const;
    void changeScale(double val) const;

    void sendOffset(const Wm3::Vector3<double>& offset) const;
    void changeOffset(const Wm3::Vector3<double>& offset) const;
signals:
    void clicked(const Wm3::Vector3<double>&) const;
    void sliceChanged(int val) const;

    void scaleChanged(double) const;
    void offsetChanged(const Wm3::Vector3<double>&) const;

private:
    Ui::ImageWidget *ui;
    Data* m_data;

    QSlider* m_horizontalSlider;
    QSlider* m_verticalSlider;

    Wm3::Matrix3<double> m_rotMat;
    std::vector< Wm3::Vector3<double> > m_cntPts;
    std::vector< Wm3::Vector3<double> > m_splPts;
};

#endif // IMAGEWIDGET_H

