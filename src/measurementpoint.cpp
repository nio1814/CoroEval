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

#include "measurementpoint.h"

#include "Data.h"
#include "settings.h"
#include "bsplineinterpolation.h"
#include "vesselsharpness.h"

#include "Wm3GMatrix.h"
#include "Wm3GVector.h"
#include "Wm3Matrix3.h"
#include "Wm3Vector3.h"
#include "Wm3NoniterativeEigen3x3.h"

#include <QDebug>

#include <cmath>
#include <algorithm>


MeasurementPoint::MeasurementPoint( MeasurementPoint::Normal n, double pos, Data* data)
	:	m_settings(Settings::getSettings()){
	m_data = data;
	m_pos  = pos;
	m_cachedDiameter = 0;
	m_cachedSharpness = 0;

	m_normals.push_back( n );
}

MeasurementPoint::MeasurementPoint(double pos, BSplineInterpolation* bSpline, Data* data)
	:	m_settings(Settings::getSettings()) {
	m_data = data;

	double length = bSpline->getLength() * m_data->getVoxelSize();
	m_pos  = pos * length;

	std::vector<Wm3::Vector3<double> > vec = bSpline->getNormal(pos);
	if (vec.size() != 3)
		return;

	m_p = vec[0];
	m_x = vec[1];
	m_n = vec[2];
	
	Wm3::Quaternion<double> q(m_n, Wm3::Mathd::HALF_PI);
	m_y = q.Rotate(m_x);

	m_cachedDiameter = 0;
	m_cachedSharpness = 0;

	setPositionInVolume(m_p);
}

double MeasurementPoint::getAverageSharpness() const {
	if (m_cachedSharpness > 1e-7)
		return m_cachedSharpness;
	
	if (m_settings.getDoNormalCheck()) {
		std::vector<size_t> dummy;
		return getAverageSharpnessChecked(dummy);
	}

	std::vector<size_t> list;
	for (size_t i = 0; i < m_normals.size(); i++) {
		const Normal& n = m_normals.at(i);

		if (!n.valid)
			continue;

		list.push_back(i);
	}

	if (list.empty())
		return 0;

	return avgSharpnessHelper(list);
}

double MeasurementPoint::getAverageSharpnessChecked(std::vector<size_t>& checkedList) const {
	if (m_cachedSharpness > 1e-7)
		return m_cachedSharpness;
	
	std::vector<size_t> list;
	for (size_t i = 0; i < m_normals.size(); i++) {
		const Normal& n = m_normals.at(i);

		if (!n.valid)
			continue;

		list.push_back(i);
	}

	if (list.empty())
		return 0;

	double allMean = avgSharpnessHelper(list);
	double thr = stdSharpnessHelper(list, allMean) * 2.0;

	double mean = 0;
	size_t count = 0;
	for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
		const Normal& n = m_normals.at(*i);

		double curMean = 0;
		size_t curCount = 0;

		VesselSharpness eval;
		eval.setNormal(n.data, m_data->getVoxelSize());
		if (isPointOk(*i, Interval::Left) && (eval.getSharpness(Interval::Left) > 0)) {
			curMean += eval.getSharpness(Interval::Left);
			++curCount;
		}

        if (isPointOk(*i, Interval::Right) && (eval.getSharpness(Interval::Right) > 0)) {
			curMean += eval.getSharpness(Interval::Right);
			++curCount;
		}

		if (curCount == 0)
			continue;

		curMean /= static_cast<double>(curCount);

		if (std::abs(allMean - curMean) <= thr) {
			mean += curMean;
			++count;
			checkedList.push_back(*i);
		}
	}

	if (count == 0)
		return 0;

	return mean /= static_cast<double>(count);
}

double MeasurementPoint::getStdSharpness() const {
	if (m_settings.getDoNormalCheck())
		return getStdSharpnessChecked();

	std::vector<size_t> list;
	for (size_t i = 0; i < m_normals.size(); i++) {
		const Normal& n = m_normals.at(i);

		if (!n.valid)
			continue;

		list.push_back(i);
	}

	if (list.empty())
		return 0;

	double mean, std;
	stdSharpnessHelper(list, mean, std);

	return std;
}

double MeasurementPoint::getStdSharpnessChecked() const {
	std::vector<size_t> list;
	double mean = getAverageSharpnessChecked(list);

	return stdSharpnessHelper(list, mean);
}

void MeasurementPoint::updatePositionInVolume()
{
    m_p = m_betterPvolume;

    double cnt = 0.0;
    m_betterPvolume = Wm3::Vector3<double>(0,0,0);
    m_betterPplane  = Wm3::Vector3<double>(0,0,0);
    for(size_t i=0; i<m_normals.size(); i++)
    {
		if (!m_normals[i].valid)
			continue;

        if( isPointOk(i, Interval::Left) )
        {
            m_betterPvolume += get50PointInVolume(i, Interval::Left);
            m_betterPplane  += get50PointInPlane(i,  Interval::Left);
            cnt += 1.0;
        }
        if( isPointOk(i, Interval::Right) )
        {
            m_betterPvolume += get50PointInVolume(i, Interval::Right);
            m_betterPplane  += get50PointInPlane(i,  Interval::Right);
            cnt += 1.0;
        }
    }

    if(cnt == 0)
    {
        for(size_t i=0; i<m_normals.size(); i++)
        {
            m_betterPvolume += 0.5*(get50PointInVolume(i, Interval::Left) + get50PointInVolume(i, Interval::Right));
            m_betterPplane  += 0.5*(get50PointInPlane(i, Interval::Left)  + get50PointInPlane(i, Interval::Right));
        }
        m_betterPvolume /= double(m_normals.size());
        m_betterPplane  /= double(m_normals.size());
    } else {
        m_betterPvolume /= cnt;
        m_betterPplane  /= cnt;
    }
}

void MeasurementPoint::recalculate()
{
	sampleNormals();
    evaluateSharpness();

    m_betterPvolume = Wm3::Vector3<double>(0,0,0);
    m_betterPplane  = Wm3::Vector3<double>(0,0,0);

	size_t count = 0;
    for(size_t i=0; i<m_normals.size(); i++)
    {
		if (!m_normals[i].valid)
			continue;

		m_betterPvolume += 0.5*(get50PointInVolume(i, Interval::Left) + get50PointInVolume(i, Interval::Right));
        m_betterPplane  += 0.5*(get50PointInPlane(i, Interval::Left)  + get50PointInPlane(i, Interval::Right));
		++count;
    }
	
	m_betterPvolume /= static_cast<double>(count);
	m_betterPplane /= static_cast<double>(count);

	if (m_cachedDiameter > 1e-3)
		m_cachedDiameter = 0;

	if (m_cachedSharpness > 1e-7)
		m_cachedSharpness = 0;
}

Wm3::Vector3<double> MeasurementPoint::getTangent() const
{
	return m_n;
}


bool MeasurementPoint::isPointOk(size_t num, Interval::Position p) const
{
	if (!m_settings.getDoPointCheck())
		return true;

    double radius = getDiameter() / m_data->getVoxelSize() / 2.0;
    double thresRadius = 1.0;
    if( qAbs((get50PointInVolume(num, p)-m_betterPvolume).Length()-radius) < thresRadius )
        return true;
    return false;
}

double MeasurementPoint::getPositionInSpline() const
{
    return m_pos;
}

void MeasurementPoint::setPositionInVolume(const Wm3::Vector3<double>& p)
{
    m_p     = p;
    m_pOrig = m_p;

    recalculate();
}

void MeasurementPoint::improvePoint(double thres, size_t iter, double radius)
{
    for(size_t i=0; i<iter; i++)
    {
        if( (m_betterPvolume-m_pOrig).Length() > radius)
        {
            //qDebug() << "Distance from original position too large: " << (m_betterPvolume-m_pOrig).Length();
            break;
        }
        if( (m_p-m_betterPvolume).Length() < thres )
        {
            //qDebug() << "After " << i << " the update from previous position is too small!";
            break;
        }
        updatePositionInVolume();
        sampleNormals();
        evaluateSharpness();
    }
}

void MeasurementPoint::improvePoint2(double thres, size_t maxIter, double radius)
{
	for (size_t iter = 0; iter < maxIter; ++iter) {
		m_betterPvolume = Wm3::Vector3<double>(0,0,0);
		m_betterPplane  = Wm3::Vector3<double>(0,0,0);
		
		double valmin = getCenterInterval(0).getMaximumValue();
		double valmax = valmin;

		for (size_t i = 1; i < m_normals.size(); ++i) {
			if (!m_normals[i].valid)
				continue;

			double val = getCenterInterval(i).getMaximumValue();
			
			if (val < valmin)
				valmin = val;
			if (val > valmax)
				valmax = val;
		}

		valmax -= valmin;
		bool perfect = (valmax == 0);

		double weightsum = 0;
		
		for (size_t i = 0; i < m_normals.size(); ++i) {
			if (!m_normals[i].valid)
				continue;

			double weight = (getCenterInterval(i).getMaximumValue() - valmin) / valmax;

			if (perfect)
				weight = 1;

			if (weight <= 0)
				continue;

			assert(weight > 0);
			
			m_betterPvolume += getMaxPointInVolume(i) * weight;
			m_betterPplane += getMaxPointInPlane(i) * weight;
			
			weightsum += weight;
		}

		assert(weightsum > 0);
		m_betterPvolume /= weightsum;
		m_betterPplane /= weightsum;

		if( (m_betterPvolume-m_pOrig).Length() > radius)
        {
            qDebug() << "Distance from original position too large: " << (m_betterPvolume-m_pOrig).Length();
            break;
        }
        if( (m_p-m_betterPvolume).Length() < thres )
        {
            qDebug() << "After " << iter << " the update from previous position is too small!";
            break;
        }

		m_p = m_betterPvolume;

		sampleNormals();
		evaluateSharpness();
	}
}

void MeasurementPoint::improvePoint2d(double thres, size_t maxIter)
{
	for (size_t iter = 0; iter < maxIter; ++iter) {
		if ((m_p - m_betterPvolume).Length() < thres)
			break;

		m_p = m_betterPplane.X() * m_x + m_betterPplane.Y() * m_y + m_p;
		
		recalculate();
	}
}

const Wm3::Vector3<double>& MeasurementPoint::getPositionInVolume() const
{
    return m_p;
}

void MeasurementPoint::getPlane(Wm3::Vector3<double> &u, Wm3::Vector3<double> &v) const
{
    u = m_x;
    v = m_y;
}

const std::vector<float>& MeasurementPoint::getNormal(size_t num) const
{
    return m_normals.at(num).data;
}

double MeasurementPoint::getSharpness(size_t num) const
{
	if (!m_normals.at(num).valid)
		return 0;
	else
		return m_normals.at(num).sharpness;
}

const Interval& MeasurementPoint::getCenterInterval(size_t num) const
{
	return m_normals.at(num).centerInterval;
}

const Wm3::Vector3<double>& MeasurementPoint::getSamplingDirection(size_t num) const
{
    if (!m_normals[num].valid)
		return Wm3::Vector3<double>::ZERO;
	
	return m_normals.at(num).v;
}

double MeasurementPoint::getAngle(size_t num) const
{
    if (!m_normals[num].valid)
		return 0;
	
	return m_normals.at(num).alpha;
}

Wm3::Vector3<double> MeasurementPoint::getMaxPointInVolume(size_t num) const
{
	if (!m_normals[num].valid)
		return Wm3::Vector3<double>::ZERO;
	
	double scalar = (static_cast<double>(getCenterInterval(num).getMaximumIndex()) - static_cast<double>(getNormal(num).size()) / 2.0) * 0.5;
	return m_p + scalar * getSamplingDirection(num);
}

Wm3::Vector3<double> MeasurementPoint::getMaxPointInPlane(size_t num) const
{
	if (!m_normals[num].valid)
		return Wm3::Vector3<double>::ZERO;
	
	double scalar = (static_cast<double>(getCenterInterval(num).getMaximumIndex()) - static_cast<double>(getNormal(num).size()) / 2.0) * 0.5;
	Wm3::Vector3<double> sd = getSamplingDirection(num);
    Wm3::Vector3<double> d =  Wm3::Vector3<double>( sd.Dot(m_x) / m_x.Dot(m_x), sd.Dot(m_y) / m_y.Dot(m_y), 0);
	return scalar * d;
}

Wm3::Vector3<double> MeasurementPoint::getMinPointInVolume(size_t num, Interval::Position p) const
{
	if (!m_normals[num].valid)
		return Wm3::Vector3<double>::ZERO;
	
	double scalar = (static_cast<double>(getCenterInterval(num).getMinimumIndex(p)) - static_cast<double>(getNormal(num).size()) / 2.0) * 0.5;
	return m_p + scalar * getSamplingDirection(num);
}

Wm3::Vector3<double> MeasurementPoint::get50PointInVolume(size_t num, Interval::Position p) const
{
	if (!m_normals[num].valid)
		return Wm3::Vector3<double>::ZERO;
	
	double scalar = (static_cast<double>(getCenterInterval(num).get50Index(p)) - static_cast<double>(getNormal(num).size()) / 2.0) * 0.5;
	return m_p + scalar * getSamplingDirection(num);
}

Wm3::Vector3<double> MeasurementPoint::get50PointInPlane(size_t num, Interval::Position p) const
{
	if (!m_normals[num].valid)
		return Wm3::Vector3<double>::ZERO;
	
	double scalar = (static_cast<double>(getCenterInterval(num).get50Index(p)) - static_cast<double>(getNormal(num).size()) / 2.0) * 0.5;
	Wm3::Vector3<double> sd = getSamplingDirection(num);
    Wm3::Vector3<double> d =  Wm3::Vector3<double>( sd.Dot(m_x) / m_x.Dot(m_x), sd.Dot(m_y) / m_y.Dot(m_y), 0);
	return scalar * d;
}

Wm3::Vector3<double> MeasurementPoint::get20PointInVolume(size_t num, Interval::Position p) const
{
    if (!m_normals[num].valid)
		return Wm3::Vector3<double>::ZERO;
	
	double scalar;
    if( p==Interval::Left )
        scalar = -1.0 * (double(getNormal(num).size())/2.0 - double(getCenterInterval(num).get20Index(p))) * 0.5;
    if( p==Interval::Right )
        scalar = (double(getCenterInterval(num).get20Index(p)) - double(getNormal(num).size())/2.0) * 0.5;
    return m_p + scalar * getSamplingDirection(num);
}

Wm3::Vector3<double> MeasurementPoint::get20PointInPlane(size_t num, Interval::Position p) const
{
	if (!m_normals[num].valid)
		return Wm3::Vector3<double>::ZERO;
	
	double scalar;
    if( p==Interval::Left )
        scalar = -1.0 * (double(getNormal(num).size())/2.0 - double(getCenterInterval(num).get20Index(p))) * 0.5;
    if( p==Interval::Right )
        scalar = (double(getCenterInterval(num).get20Index(p)) - double(getNormal(num).size())/2.0) * 0.5;
    Wm3::Vector3<double> sd = getSamplingDirection(num);
    Wm3::Vector3<double> d =  Wm3::Vector3<double>( sd.Dot(m_x) / m_x.Dot(m_x), sd.Dot(m_y) / m_y.Dot(m_y), 0);
    return scalar * d;
}

const Wm3::Vector3<double>& MeasurementPoint::getImprovedPositionInVolume() const
{
    return m_betterPvolume;
}

const Wm3::Vector3<double>& MeasurementPoint::getImprovedPositionInPlane() const
{
    return m_betterPplane;
}

double MeasurementPoint::getDiameter() const {
	if (m_cachedDiameter > 1e-3)
		return m_cachedDiameter;
	
	double d = 0;

	if (m_settings.getUseEllipse())
		d = fitEllipse();

	if (d > 1e-3)
		return d;
	
	if (m_settings.getDoNormalCheck())
		return getDiameterChecked();

	std::vector<size_t> list;
	for (size_t i = 0; i < m_normals.size(); ++i) {
		if (!m_normals[i].valid)
			continue;
		
		list.push_back(i);
	}

	if (list.empty())
		return 0;

	return diameterHelper(list);
}

double MeasurementPoint::getDiameterVoxels() const
{
    return getDiameter() / m_data->getVoxelSize();
}

double MeasurementPoint::getDiameterChecked() const {
	std::vector<size_t> list;
	for (size_t i = 0; i < m_normals.size(); i++) {
		const Normal& n = m_normals.at(i);

		if (!n.valid)
			continue;

		list.push_back(i);
	}

	if (list.empty())
		return 0;

	double allMean = medDiameterHelper(list);
	double thr = stdDiameterHelper(list, allMean) * 2.0;

	Wm3::Vector3d medMax = medMaxHelper(list);
	Wm3::Vector3d maxThr = stdMaxHelper(list, medMax) * 2.0;

	double mean = 0;
	size_t count = 0;
	for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
		Wm3::Vector3<double> maxP = getMaxPointInVolume(*i);

		double len = (maxP - get50PointInVolume(*i, Interval::Left)).Length();
		len += (maxP - get50PointInVolume(*i, Interval::Right)).Length();

		if (len <= 0)
			continue;

		len *= m_data->getVoxelSize();

		if (std::abs(allMean - len) > thr)
			continue;

		Wm3::Vector3d len2 = (medMax - maxP);
		if ((std::abs(len2[0]) > maxThr[0]) || (std::abs(len2[1]) > maxThr[1]) || (std::abs(len2[2]) > maxThr[2]))
			continue;

		mean += len;
		++count;
	}

	if (count == 0)
		return 0;

	return mean /= static_cast<double>(count);
}

std::vector<Wm3::Vector3<double> > MeasurementPoint::get50Points() const {
	std::vector<Wm3::Vector3<double> > list;

	for (size_t i = 0; i < m_normals.size(); ++i) {
		if (!m_normals[i].valid)
			continue;
		
		list.push_back(getMaxPointInPlane(i));

		list.push_back(get50PointInPlane(i, Interval::Left));
		list.push_back(get50PointInPlane(i, Interval::Right));
	}

	return list;
}

static double sgn(double n) {
	if (n > 0) 
		return 1.0;
	if (n < 0)
		return -1.0;
	return 0;
}

double MeasurementPoint::fitEllipse() const {
	// Find all valid normals
	std::vector<size_t> list;
	for (size_t i = 0; i < m_normals.size(); i++) {
		const Normal& n = m_normals.at(i);

		if (!n.valid)
			continue;

		list.push_back(i);
	}

	if (list.empty())
		return 0;

	// Remove gross outliers and calculate normalisation factors
	double allMean = medDiameterHelper(list);
	double thr = stdDiameterHelper(list, allMean) * 2.0;

	Wm3::Vector3d medMax = medMaxHelper(list);
	Wm3::Vector3d maxThr = stdMaxHelper(list, medMax) * 2.0;

	std::vector<Wm3::Vector3<double> > x;
	double meanX = 0;
	double meanY = 0;
	double maxX = -1000;
	double maxY = -1000;
	double minX = 1000;
	double minY = 1000;

	for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
		Wm3::Vector3<double> maxP = getMaxPointInVolume(*i);

		double len = (maxP - get50PointInVolume(*i, Interval::Left)).Length();
		len += (maxP - get50PointInVolume(*i, Interval::Right)).Length();

		if (len <= 0)
			continue;

		len *= m_data->getVoxelSize();

		if (std::abs(allMean - len) > thr)
			continue;

		Wm3::Vector3d len2 = (medMax - maxP);
		if ((std::abs(len2[0]) > maxThr[0]) || (std::abs(len2[1]) > maxThr[1]) || (std::abs(len2[2]) > maxThr[2]))
			continue;

		Wm3::Vector3d pL = get50PointInPlane(*i, Interval::Left);
		Wm3::Vector3d pR = get50PointInPlane(*i, Interval::Right);

		x.push_back(pL);
		x.push_back(pR);

		meanX += pL.X();
		meanX += pR.X();
			
		meanY += pL.Y();
		meanY += pR.Y();

		if (pL.X() < minX)
			minX = pL.X();
		if (pL.Y() < minY)
			minY = pL.Y();
		if (pR.X() < minX)
			minX = pR.X();
		if (pR.Y() < minY)
			minY = pR.Y();

		if (pL.X() > maxX)
			maxX = pL.X();
		if (pL.Y() > maxY)
			maxY = pL.Y();
		if (pR.X() > maxX)
			maxX = pR.X();
		if (pR.Y() > maxY)
			maxY = pR.Y();
	}

	if (x.empty())
		return 0;

	// Normalise input points
	size_t n = x.size();

	meanX /= static_cast<double>(n);
	meanY /= static_cast<double>(n);

	double sx = (maxX - minX) / 2.0;
	double sy = (maxY - minY) / 2.0;

	for (size_t i = 0; i < x.size(); ++i) {
		x[i][0] = (x[i][0] - meanX) / sx;
		x[i][1] = (x[i][1] - meanY) / sy;
	}

	// Fill design matrix
	Wm3::GMatrixd D1(n, 3);
	Wm3::GMatrixd D2(n, 3);
	Wm3::GVectord tmp(3);

	for (size_t i = 0; i < n; ++i) {
		double vX = x[i][0];
		double vY = x[i][1];
		
		tmp[0] = vX * vX;
		tmp[1] = vX * vY;
		tmp[2] = vY * vY;

		D1.SetRow(i, tmp);

		tmp[0] = vX;
		tmp[1] = vY;
		tmp[2] = 1.0;

		D2.SetRow(i, tmp);
	}

	// Build scatter matrices
	Wm3::GMatrixd tmp1 = D1.TransposeTimes(D1);
	Wm3::GMatrixd tmp2 = D1.TransposeTimes(D2);
	Wm3::GMatrixd tmp3 = D2.TransposeTimes(D2);

	Wm3::Matrix3d S1;
	Wm3::Matrix3d S2;
	Wm3::Matrix3d S3;

	for (size_t r = 0; r < 3; ++r)
		for (size_t c = 0; c < 3; ++c) {
			S1(r,c) = tmp1(r,c);
			S2(r,c) = tmp2(r,c);
			S3(r,c) = tmp3(r,c);
		}
	
	// Inverse constraint matrix
	Wm3::Matrix3d Cinv(Wm3::Matrix3d::ZERO);
	Cinv(0, 2) = 0.5;
	Cinv(1, 1) = -1.0;
	Cinv(2, 0) = 0.5;

	// Solve eigensystem
	Wm3::Matrix3d E = -S3.Inverse() * S2.Transpose();
	Wm3::Matrix3d M = Cinv * (S1 + S2 * E);
	
	Wm3::NoniterativeEigen3x3d eig(M);
	const Wm3::Vector3d* ev = eig.GetEigenvectors();

	size_t idx = 4;
	for (size_t i = 0; i < 3; ++i) {
		double cond = 4 * ev[i][0] * ev[i][2] - ev[i][1] * ev[i][1];
		if (cond > 0) {
			if (idx != 4)
				return 0;
			idx = i;
		}
	}

	if (idx == 4)
		return 0;

	Wm3::Vector3d a1(ev[idx]);
	Wm3::Vector3d a2 = E * a1;
	
	// De-normalise
	Wm3::GVectord par(6);

	par[0] = a1[0] * sy * sy;
	par[1] = a1[1] * sx * sy;
	par[2] = a1[2] * sx * sx;
	par[3] = -2.0 * a1[0] * sy * sy * meanX - a1[1] * sx * sy * meanY + a2[0] * sx * sy * sy;
	par[4] = -a1[1] * sx * sy * meanX - 2.0 * a1[2] * sx * sx * meanY + a2[1] * sx * sx * sy;
	par[5] = a1[0] * sy * sy * meanX * meanX + a1[1] * sx * sy * meanX * meanY
		   + a1[2] * sx * sx * meanY * meanY - a2[0] * sx * sy * sy * meanX
		   - a2[1] * sx * sx * sy * meanY + a2[2] * sx * sx * sy * sy;

	double phi = 0.5 * std::atan2(par[1], par[0] - par[2]);
	double cost = std::cos(phi);
	double sint = std::sin(phi);
	double sin_squared = sint * sint;
	double cos_squared = cost * cost;
	double cos_sin = sint * cost;

	double Au = par[3] * cost + par[4] * sint;
	double Av = -par[3] * sint + par[4] * cost;
	double Auu = par[0] * cos_squared + par[2] * sin_squared + par[1] * cos_sin;
	double Avv = par[0] * sin_squared + par[2] * cos_squared - par[1] * cos_sin;

	double tuCentre = -Au / (2.0 * Auu);
	double tvCentre = -Av / (2.0 * Avv);

	double wCentre = par[5] - Auu * tuCentre * tuCentre - Avv * tvCentre * tvCentre;

	double Ru = -wCentre / Auu;
	double Rv = -wCentre / Avv;
	Ru = std::sqrt(std::abs(Ru)) * sgn(Ru);
	Rv = std::sqrt(std::abs(Rv)) * sgn(Rv);

	double d = std::sqrt(Ru * Rv) * 2.0 * m_data->getVoxelSize();
	
	return d;
}

void MeasurementPoint::setCachedDiameter(double d) {
	if (d > 0)
		m_cachedDiameter = d;
}

void MeasurementPoint::setCachedSharpness(double s) {
	if (s > 0)
		m_cachedSharpness = s;
}

void MeasurementPoint::sampleNormals()
{
    m_normals.clear();

	double hwidth  = static_cast<double>(m_settings.getSizeNormal()) / m_data->getVoxelSize() / 2.0;
    size_t numNormals = m_settings.getNumNormalsAtPoint();

    double angle = Wm3::Mathd::PI / static_cast<double>(numNormals);
    for(size_t i = 0; i < numNormals; ++i)
    {
        Normal n;

        n.alpha = i * angle;
        Wm3::Quaternion<double> q(m_n, n.alpha);
        n.v = q.Rotate(m_x);

        // Sample the normal (2x oversampled)
        for (double j = -hwidth; j <= hwidth; j += 0.5)
        {
            Wm3::Vector3<double> currPt = m_p + n.v * j;
            float val = m_data->getValue( currPt.X(), currPt.Y(), currPt.Z() );
            n.data.push_back(val);
        }

        m_normals.push_back(n);
    }
}

void MeasurementPoint::evaluateSharpness()
{
    for(size_t i=0; i<m_normals.size(); i++)
    {
        Normal& n = m_normals.at(i);

        VesselSharpness eval;
		eval.setNormal( n.data, m_data->getVoxelSize() );

		n.valid = eval.isValid();
		if (!n.valid)
			continue;

        n.sharpness = eval.getSharpness(Interval::Both);
        n.diffMaxToCenter = eval.getCenterDistance()*0.5;
        n.centerInterval = eval.getCenterInterval();
    }
}

double MeasurementPoint::avgSharpnessHelper(const std::vector<size_t>& list) const {
	double meanSharpness = 0;
	size_t cnt = 0;
	for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
		const Normal& n = m_normals.at(*i);

		VesselSharpness eval;
		eval.setNormal(n.data, m_data->getVoxelSize());
		if (isPointOk(*i, Interval::Left) && (eval.getSharpness(Interval::Left) > 0)) {
			meanSharpness += eval.getSharpness(Interval::Left);
			++cnt;
		}

        if (isPointOk(*i, Interval::Right) && (eval.getSharpness(Interval::Right) > 0)) {
			meanSharpness += eval.getSharpness(Interval::Right);
			++cnt;
		}
	}

	if (cnt == 0)
		return 0;

	return meanSharpness /= static_cast<double>(cnt);
}

void MeasurementPoint::stdSharpnessHelper(const std::vector<size_t>& list, double& mean, double& std) const {
	if (m_cachedSharpness > 1e-7)
		mean = m_cachedSharpness;
	else
		mean = avgSharpnessHelper(list);

	std = stdSharpnessHelper(list, mean);
}

double MeasurementPoint::stdSharpnessHelper(const std::vector<size_t>& list, double mean) const {
	double stdSharpness = 0;
	size_t cnt = 0;
    for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
		const Normal& n = m_normals.at(*i);

		VesselSharpness eval;
		eval.setNormal(n.data, m_data->getVoxelSize());
		if (isPointOk(*i, Interval::Left) && (eval.getSharpness(Interval::Left) > 0)) {
			double val = (eval.getSharpness(Interval::Left) - mean);
			stdSharpness += val * val;
			++cnt;
		}

		if (isPointOk(*i, Interval::Right) && (eval.getSharpness(Interval::Right) > 0)) {
			double val = (eval.getSharpness(Interval::Right) - mean);
			stdSharpness += val * val;
			++cnt;
		}
	}

	if (cnt == 0)
		return 0;

	if (cnt > 1)
		return sqrt(stdSharpness / static_cast<double>((cnt - 1)));
	else
		return sqrt(stdSharpness / static_cast<double>(cnt));
}

double MeasurementPoint::diameterHelper(const std::vector<size_t>& list) const {
	// New code: full width at half maximum and actual min/max points for each normal
	double radius = 0;
	size_t count = 0;

	for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
		Wm3::Vector3<double> maxP = getMaxPointInVolume(*i);

		double len = (maxP - get50PointInVolume(*i, Interval::Left)).Length();
		len += (maxP - get50PointInVolume(*i, Interval::Right)).Length();

		if (len <= 0)
			continue;

		radius += len;
		++count;
	}

	if (count == 0)
		return 0;

	double diameter = radius / static_cast<double>(count);
    
	return diameter * m_data->getVoxelSize();
}

double MeasurementPoint::medDiameterHelper(const std::vector<size_t>& list) const {
	// New code: full width at half maximum and actual min/max points for each normal
	std::vector<double> dia;
	
	for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
		Wm3::Vector3<double> maxP = getMaxPointInVolume(*i);

		double len = (maxP - get50PointInVolume(*i, Interval::Left)).Length();
		len += (maxP - get50PointInVolume(*i, Interval::Right)).Length();

		if (len <= 0)
			continue;

		dia.push_back(len);
	}

	if (dia.size() == 0)
		return 0;

	std::sort(dia.begin(), dia.end());
	double diameter = 0;
	if (dia.size() % 2 == 0) {
		diameter = dia[dia.size() / 2 - 1];
		diameter += dia[dia.size() / 2];
		diameter /= 2.0;
	} else
		diameter = dia[(dia.size() + 1) / 2 - 1];
    
	return diameter * m_data->getVoxelSize();
}

double MeasurementPoint::stdDiameterHelper(const std::vector<size_t>& list, double mean) const {
	double stdDiameter = 0;
	size_t cnt = 0;
    for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
		Wm3::Vector3<double> maxP = getMaxPointInVolume(*i);

		double len = (maxP - get50PointInVolume(*i, Interval::Left)).Length();
		len += (maxP - get50PointInVolume(*i, Interval::Right)).Length();

		if (len <= 0)
			continue;

        len *= m_data->getVoxelSize();

		double val = (len - mean);
		stdDiameter += val * val;
		++cnt;
	}

	if (cnt == 0)
		return 0;

	if (cnt > 1)
		return sqrt(stdDiameter / static_cast<double>((cnt - 1)));
	else
		return sqrt(stdDiameter / static_cast<double>(cnt));
}

Wm3::Vector3<double> MeasurementPoint::medMaxHelper(const std::vector<size_t>& list) const {
	std::vector<double> maxX;
	std::vector<double> maxY;
	std::vector<double> maxZ;
	
	for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
		Wm3::Vector3<double> maxP = getMaxPointInVolume(*i);

		maxX.push_back(maxP.X());
		maxY.push_back(maxP.Y());
		maxZ.push_back(maxP.Z());
	}

	std::sort(maxX.begin(), maxX.end());
	std::sort(maxY.begin(), maxY.end());
	std::sort(maxZ.begin(), maxZ.end());
		
	Wm3::Vector3<double> medMax;
	if (maxX.size() % 2 == 0) {
		medMax[0] = maxX[maxX.size() / 2 - 1];
		medMax[1] = maxY[maxX.size() / 2 - 1];
		medMax[2] = maxZ[maxX.size() / 2 - 1];
		
		medMax[0] += maxX[maxX.size() / 2];
		medMax[1] += maxY[maxX.size() / 2];
		medMax[2] += maxZ[maxX.size() / 2];
		
		medMax /= 2.0;
	} else {
		medMax[0] = maxX[(maxX.size() + 1) / 2 - 1];
		medMax[1] = maxY[(maxX.size() + 1) / 2 - 1];
		medMax[2] = maxZ[(maxX.size() + 1) / 2 - 1];
	}
    
	return medMax;
}

Wm3::Vector3<double> MeasurementPoint::stdMaxHelper(const std::vector<size_t>& list, const Wm3::Vector3<double>& mean) const {
	Wm3::Vector3<double> stdMax(Wm3::Vector3<double>::ZERO);
	size_t cnt = 0;
    for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
		Wm3::Vector3<double> maxP = getMaxPointInVolume(*i);

		Wm3::Vector3<double> val = (maxP - mean);
		stdMax[0] += val[0] * val[0];
		stdMax[1] += val[1] * val[1];
		stdMax[2] += val[2] * val[2];

		++cnt;
	}

	if (cnt == 0)
		return 0;

	if (cnt > 1)
		stdMax /= static_cast<double>(cnt - 1);
	else
		stdMax /= static_cast<double>(cnt);

	stdMax[0] = std::sqrt(stdMax[0]);
	stdMax[1] = std::sqrt(stdMax[1]);
	stdMax[2] = std::sqrt(stdMax[2]);

	return stdMax;
}

void MeasurementPoint::updateDiameterOutlierStatus()
{
    std::vector<size_t> list;
    for (size_t i = 0; i < m_normals.size(); i++) {
        Normal& n = m_normals.at(i);

        if (!n.valid)
        {
            n.diameterOutlier[0] = n.diameterOutlier[1] = true;
            continue;
        }

        list.push_back(i);
    }

    if (list.empty())
        return;

    // Eliminate normals where the normal max point is too far from the average normal max point

    std::vector<size_t> meanList;

    {
        Wm3::Vector3d meanMaxPosition = Wm3::Vector3d::ZERO;
        for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
            meanMaxPosition += getMaxPointInVolume(*i);
        }
        meanMaxPosition /= static_cast<double>(list.size());

        double meanDistMaxPosition = 0.0;
        for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
            meanDistMaxPosition += (getMaxPointInVolume(*i) - meanMaxPosition).Length();
        }
        meanDistMaxPosition /= static_cast<double>(list.size());

        double stdDistMaxPosition = 0.0;
        for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
            double dist = (getMaxPointInVolume(*i) - meanMaxPosition).Length() - meanDistMaxPosition;
            stdDistMaxPosition += dist * dist;
        }
        stdDistMaxPosition /= static_cast<double>(list.size() > 1 ? list.size() - 1 : list.size());
        stdDistMaxPosition = sqrt(stdDistMaxPosition);

        double thr = 2.0 * stdDistMaxPosition;

        for (std::vector<size_t>::const_iterator i = list.begin(); i != list.end(); ++i) {
            if (std::abs((getMaxPointInVolume(*i) - meanMaxPosition).Length() - meanDistMaxPosition) <= thr)
            {
                meanList.push_back(*i);
            }
            else
            {
                Normal& n = m_normals.at(*i);
                n.diameterOutlier[0] = n.diameterOutlier[1] = true;
            }
        }
    }

    // Eliminate surface points which deviate too much from the median radius

    if (meanList.empty())
        meanList = list;

    double medianRadius;

    {
        std::vector<double> rad;

        for (std::vector<size_t>::const_iterator i = meanList.begin(); i != meanList.end(); ++i) {
            Wm3::Vector3<double> maxP = getMaxPointInVolume(*i);

            double len = (maxP - get50PointInVolume(*i, Interval::Left)).Length();

            if (len <= 0)
                continue;

            rad.push_back(len);

            len = (maxP - get50PointInVolume(*i, Interval::Right)).Length();

            if (len <= 0)
                continue;

            rad.push_back(len);
        }

        if (rad.size() == 0)
            medianRadius = 0;
        else {
            std::sort(rad.begin(), rad.end());
            double medianRadiusVoxels = 0;
            if (rad.size() % 2 == 0) {
                medianRadiusVoxels = rad[rad.size() / 2 - 1];
                medianRadiusVoxels += rad[rad.size() / 2];
                medianRadiusVoxels /= 2.0;
            } else
                medianRadiusVoxels = rad[(rad.size() + 1) / 2 - 1];

            medianRadius = medianRadiusVoxels * m_data->getVoxelSize();
        }
    }

    double thr;

    // Use the global center point instead of the profile line max
    {
        double stdRadius = 0;
        size_t cnt = 0;
        for (std::vector<size_t>::const_iterator i = meanList.begin(); i != meanList.end(); ++i) {
            Wm3::Vector3<double> maxP = getMaxPointInVolume(*i);

            double len = (maxP - get50PointInVolume(*i, Interval::Left)).Length();

            if (len <= 0)
                continue;

            len *= m_data->getVoxelSize();

            double val = (len - medianRadius);
            stdRadius += val * val;
            ++cnt;

            len = (maxP - get50PointInVolume(*i, Interval::Right)).Length();

            if (len <= 0)
                continue;

            len *= m_data->getVoxelSize();

            val = (len - medianRadius);
            stdRadius += val * val;
            ++cnt;
        }

        if (cnt == 0)
            thr = 0;
        else if (cnt > 1)
            thr = 2.0 * sqrt(stdRadius / static_cast<double>((cnt - 1)));
        else
            thr = 2.0 * sqrt(stdRadius / static_cast<double>(cnt));
    }

    for (std::vector<size_t>::const_iterator i = meanList.begin(); i != meanList.end(); ++i) {
        Normal& n = m_normals.at(*i);

        double leftRadius  = (getMaxPointInVolume(*i) - get50PointInVolume(*i, Interval::Left )).Length() * m_data->getVoxelSize();
        double rightRadius = (getMaxPointInVolume(*i) - get50PointInVolume(*i, Interval::Right)).Length() * m_data->getVoxelSize();

        n.diameterOutlier[0] = !(isPointOk(*i, Interval::Left)  && std::abs(medianRadius -  leftRadius) <= thr);
        n.diameterOutlier[1] = !(isPointOk(*i, Interval::Right) && std::abs(medianRadius - rightRadius) <= thr);
    }
}

bool MeasurementPoint::isPointDiameterOutlier(size_t num, Interval::Position p) const
{
    return m_normals.at(num).diameterOutlier[p == Interval::Right];
}
