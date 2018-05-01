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

#include "bsplineinterpolation.h"


BSplineInterpolation::BSplineInterpolation()
{
    m_bSplineCurve   = 0;
    m_bSmoothSamples = true;
}

BSplineInterpolation::~BSplineInterpolation()
{
    if(m_bSplineCurve)
        delete m_bSplineCurve;
}

void BSplineInterpolation::reset()
{
    m_bSplineCurve = 0;
    m_samples.clear();
}

void BSplineInterpolation::setSmoothSamples(bool val)
{
    m_bSmoothSamples = val;
}

std::vector< Wm3::Vector3<double> > BSplineInterpolation::getInterpolatedPoints() const
{
    std::vector< Wm3::Vector3<double> > res;
    if( m_bSplineCurve )
        res = getInterpolatedPoints(1.0 / (getLength()*2) );
    return res;
}

std::vector< Wm3::Vector3<double> > BSplineInterpolation::getInterpolatedPoints(double sampleRate) const
{
    std::vector< Wm3::Vector3<double> > res;
    if(!m_bSplineCurve)
        return res;

    for (double i = 0.0, k = 0; i <= 1.0; i += sampleRate, ++k) {
         res.push_back( m_bSplineCurve->GetPosition(i) );
    }
    return res;
}

double BSplineInterpolation::getLength() const
{
    if(m_bSplineCurve)
        return m_bSplineCurve->GetLength(0, 1);
    return -1;
}

std::vector< Wm3::Vector3<double> > BSplineInterpolation::getSamples() const
{
    std::vector< Wm3::Vector3<double> > res;
    for(size_t i=0; i<m_samples.size(); i++)
    {
        Wm3::Vector3<double> s(m_samples.at(i).x,
                 m_samples.at(i).y,
                 m_samples.at(i).z);
        res.push_back(s);
    }
    return res;
}

std::vector< Wm3::Vector3<double> > BSplineInterpolation::getProfile(double splPos) const
{
    std::vector< Wm3::Vector3<double> > res;

    if(!m_bSplineCurve)
        return res;

    Vector3Type p, T, N, B;
    m_bSplineCurve->GetFrame(splPos, p, T, N, B);

    QuaternionType q1(T, 0);
    Vector3Type v1 = q1.Rotate(N);
    QuaternionType q2(T, Wm3::Mathd::HALF_PI);
    Vector3Type v2 = q2.Rotate(N);

    res.push_back(p);
    res.push_back(v1);
    res.push_back(v2);

    return res;
}

Wm3::Vector3<double> BSplineInterpolation::pointAt(double splPos) const
{
    if(!m_bSplineCurve)
        return Wm3::Vector3<double>::ZERO;

    return m_bSplineCurve->GetPosition(splPos);
}

std::vector< Wm3::Vector3<double> > BSplineInterpolation::getNormal(double splPos) const
{
    std::vector< Wm3::Vector3<double> > res;

    if(!m_bSplineCurve)
        return res;

    Vector3Type p, T, N, B;
    m_bSplineCurve->GetFrame(splPos, p, T, N, B);

    res.push_back(p);
    res.push_back(N);
    res.push_back(T);

    return res;
}

void BSplineInterpolation::setSamples(const std::vector< Wm3::Vector3<double> >& ptList)
{

    if(ptList.size() < 10)
        return;

    m_samples.clear();

    // Fill samples vector
    for(size_t i=0; i<ptList.size(); i++)
    {
        Sample s(ptList.at(i).X(),
                 ptList.at(i).Y(),
                 ptList.at(i).Z());
        m_samples.push_back(s);
    }

    if(m_bSmoothSamples)
        smoothSamples();

    size_t numSamples = m_samples.size();

    // Prepare data for spline fit
    double* sampleArray = new double[3 * numSamples];

    for (size_t i = 0, j = 0, k = 0; k < numSamples; ++k, ++i, j += 3) {
        sampleArray[j]     = static_cast<double>(m_samples[i].x);
        sampleArray[j + 1] = static_cast<double>(m_samples[i].y);
        sampleArray[j + 2] = static_cast<double>(m_samples[i].z);
    }

    // Heuristic by Martin Spiegel
    size_t numCtrlPoints;
    if (numSamples < 25)
        numCtrlPoints = static_cast<size_t>(std::ceil(numSamples * 0.5f));
    else if (numSamples < 45)
        numCtrlPoints = static_cast<size_t>(std::ceil(numSamples * 0.49f));
    else if (numSamples < 75)
        numCtrlPoints = static_cast<size_t>(std::ceil(numSamples * 0.46f));
    else if (numSamples < 200)
        numCtrlPoints = static_cast<size_t>(std::ceil(numSamples * 0.42f));
    else
        numCtrlPoints = static_cast<size_t>(std::ceil(numSamples * 0.1f));

    // Fit the spline
    BSplineFitType* bSplineLeastSq = new BSplineFitType(3, numSamples, sampleArray, 3, numCtrlPoints);
    Vector3Type* ctrlPointData     = new Vector3Type[numCtrlPoints];

    if (m_bSplineCurve)
		delete m_bSplineCurve;
    m_bSplineCurve = new BSpline3Type(numCtrlPoints, ctrlPointData, 3, false, true);

    const double* ctrlPoints = bSplineLeastSq->GetControlData();

    for (size_t i = 0; i < numCtrlPoints; ++i, ctrlPoints += 3) {
        Vector3Type ctrlPointData(ctrlPoints);
        m_bSplineCurve->SetControlPoint(i, ctrlPointData);
    }

    delete[] ctrlPointData;
    delete bSplineLeastSq;
    delete[] sampleArray;
}

double* BSplineInterpolation::gaussWin(size_t size, double sigma)
{
    double* gaussWin = new double[size];
    double sum = 0;
    size_t halfWindow = (size - 1) / 2;

    for (int i = -static_cast<int>(halfWindow); i <= static_cast<int>(halfWindow); ++i) {
        gaussWin[i + halfWindow] = std::exp(-0.5 * (static_cast<double>(i) / sigma) * (static_cast<double>(i) / sigma));
        sum += gaussWin[i + halfWindow];
    }

    for (size_t i = 0; i < size; ++i)
        gaussWin[i] /= sum;

    return gaussWin;
}

void BSplineInterpolation::smoothSamples(size_t windowSize, double sigma)
{
    if (!(windowSize % 2))
        ++windowSize;

    int halfWindow = (windowSize - 1) / 2;
    double* gauss = gaussWin(windowSize, sigma);

    int length = m_samples.size();

	SampleVectorType smooth(m_samples);

    for (int i = 0; i < length; ++i)
    {
        double x = 0, y = 0, z = 0;

        for (int j = -halfWindow; j <= halfWindow; ++j)
        {
            if (i + j < 0) {
                x += m_samples[i].x * gauss[j + halfWindow];
                y += m_samples[i].y * gauss[j + halfWindow];
                z += m_samples[i].z * gauss[j + halfWindow];
            } else if (i + j >= length) {
                x += m_samples[length - 1].x * gauss[j + halfWindow];
                y += m_samples[length - 1].y * gauss[j + halfWindow];
                z += m_samples[length - 1].z * gauss[j + halfWindow];
            } else {
                x += m_samples[i + j].x * gauss[j + halfWindow];
                y += m_samples[i + j].y * gauss[j + halfWindow];
                z += m_samples[i + j].z * gauss[j + halfWindow];
            }
        }

        smooth[i] = Sample(x, y, z);
    }

    delete[] gauss;

	m_samples = smooth;
}

