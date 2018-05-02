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

#ifndef ORTHOVIEWER_H
#define ORTHOVIEWER_H

#include <QWidget>
#include <QVBoxLayout>

#include "Wm3Vector3.h"

#define cimg_display 0
#include "CImg.h"

#include <algorithm>
#include <vector>


// Forward declarations
class Data;
class ImageWidget;
class PointWidget;
class BSplineInterpolation;
class MipDialog;

namespace Ui {
    class OrthoWidget;
}


class OrthoViewer : public QWidget
{
    Q_OBJECT
public:
    explicit OrthoViewer(QWidget* parent = 0);
    ~OrthoViewer();

	void reset();
    void setData(Data* data);
    void setInterpolator(BSplineInterpolation* bSpline);

signals:
    void setWindowing(float, float) const;
    void sendControlPoints(const std::vector< Wm3::Vector3<double> >&) const;
    void sendSplinePoints(const std::vector< Wm3::Vector3<double> >&) const;
    void changeSlice(const Wm3::Vector3<double>&) const;
    void updateProfile() const;

public slots:
    void wCenterChanged(int val);
    void wWidthChanged(int val);
    void autoWindowingChanged(bool val);
    void addPoint(const Wm3::Vector3<double>& pt);
    void showSlice(const Wm3::Vector3<double>& pt);
    void update();

    void seg_clicked() const;
    void move_clicked() const;

	void showPointsDialog();

private:
    Ui::OrthoWidget* ui;
    PointWidget*     m_pointWidget;

    Data* m_data;
    BSplineInterpolation* m_bSpline;

    std::vector< Wm3::Vector3<double> > m_sampleList;
    std::vector< Wm3::Vector3<double> > m_splineList;

    float m_maximum;
    float m_minimum;
    float m_wcenter;
    float m_wwidth;

    ImageWidget*    createDialog(Data* data, QString comment);

    MipDialog* m_mipDialog;
};

#endif // ORTHOVIEWER_H

