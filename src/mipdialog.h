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

#ifndef MIPDIALOG_H
#define MIPDIALOG_H

#include <QDialog>
#include <QWidget>

#include <vector>

#include "Wm3Vector3.h"
#include "Wm3Quaternion.h"

#define cimg_display 0
#include "CImg.h"


// Forward declarations
class QSlider;
class QLabel;
class BSplineInterpolation;
class Data;


class MipDialog : public QDialog
{
    Q_OBJECT
public:
    explicit MipDialog(QWidget* parent = 0);
    void setData(BSplineInterpolation* bSpline, Data* data);
public slots:
    void updateAlpha(int val);
    void updateWindowing(float wwidth, float wcenter);
private:
    void renderMip();
    double   m_alpha;
    QSlider* m_slider;
    QLabel* m_imageLabel;
    std::vector< Wm3::Vector3<double> > m_splinePts;
    BSplineInterpolation* m_bSpline;

    float  m_wwidth;
    float  m_wcenter;

    Data*  m_data;
    cimg_library::CImg<float>* m_img;

    inline float computeCutOff(float x, float xc, float beta = 10.0) const
    {
        return 0.5 + 1.0 * Wm3::Mathf::PI * atan(beta * (xc-x)/xc );
    }
};

#endif // MIPDIALOG_H

