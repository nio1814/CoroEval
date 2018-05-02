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

#ifndef BSPLINEINTERPOLATION_H
#define BSPLINEINTERPOLATION_H

#include "Wm3Matrix3.h"
#include "Wm3BSplineFit.h"
#include "Wm3Vector3.h"
#include "Wm3BSplineCurve3.h"
#include "Wm3Quaternion.h"

#include <vector>


class BSplineInterpolation
{
public:
    BSplineInterpolation();
    ~BSplineInterpolation();

    void reset();
    void setSmoothSamples(bool val);
    std::vector< Wm3::Vector3<double> > getInterpolatedPoints() const;
    std::vector< Wm3::Vector3<double> > getInterpolatedPoints(double sampleRate) const;
    double getLength() const;

    Wm3::Vector3<double> pointAt(double splPos) const;
    std::vector< Wm3::Vector3<double> > getNormal(double splPos) const;
    std::vector< Wm3::Vector3<double> > getProfile(double splPos) const;

    void setSamples(const std::vector< Wm3::Vector3<double> >& ptList);
    std::vector< Wm3::Vector3<double> > getSamples() const;

	static double* gaussWin(size_t size, double sigma);

private:

    // Helper structure for parsing the points file
    struct Sample {
        Sample()
            :	x(0), y(0), z(0) {}
        Sample(double px, double py, double pz)
            :	x(px), y(py), z(pz) {}

        double x;
        double y;
        double z;
    };
    typedef std::vector<Sample>         SampleVectorType;
    typedef Wm3::BSplineCurve3<double>  BSpline3Type;
    typedef Wm3::Vector3<double>        Vector3Type;
    typedef Wm3::BSplineFit<double>     BSplineFitType;
    typedef std::vector<Vector3Type>    VectorVectorType;
    typedef Wm3::Quaternion<double>     QuaternionType;


    bool             m_bSmoothSamples;
    SampleVectorType m_samples;
    BSpline3Type*    m_bSplineCurve;

    void smoothSamples(size_t windowSize = 3, double sigma = 1.0);
};

#endif // BSPLINEINTERPOLATION_H

