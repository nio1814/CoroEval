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

#ifndef MEASUREMENTPOINT_H
#define MEASUREMENTPOINT_H

#include "Wm3Vector3.h"
#include "Wm3Quaternion.h"

#include "Interval.h"

#include <vector>


// Forward declarations
class Data;
class Settings;
class BSplineInterpolation;


class MeasurementPoint
{
public:
	 struct Normal {
        bool                 valid;
        bool                 diameterOutlier[2];
        std::vector<float>   data;
        Wm3::Vector3<double> v;
        Interval             centerInterval;
        double               alpha;
        double               sharpness;
        double               diffMaxToCenter;
    };

	MeasurementPoint(MeasurementPoint::Normal n, double pos, Data* data);
    MeasurementPoint(double pos, BSplineInterpolation* bSpline, Data* data);

    void improvePoint(double thres = 0.01, size_t iter = 50, double radius = 3.0);
	void improvePoint2(double thres = 0.01, size_t maxIter = 50, double radius = 10.0);
	void improvePoint2d(double thres = 0.1, size_t maxIter = 50);

    void updatePositionInVolume();
    void updateDiameterOutlierStatus();

	void recalculate();

    double getAverageSharpness() const;
    double getStdSharpness() const;
    double getPositionInSpline() const;
    double getSharpness(size_t num) const;
    const Interval& getCenterInterval(size_t num) const;
    const std::vector<float>& getNormal(size_t num) const;

    void setPositionInVolume(const Wm3::Vector3<double>& p);
    const Wm3::Vector3<double>& getPositionInVolume() const;
    const Wm3::Vector3<double>& getSamplingDirection(size_t num) const;
    double getAngle(size_t num) const;
    void getPlane(Wm3::Vector3<double> &u, Wm3::Vector3<double> &v) const;
	Wm3::Vector3<double> getTangent() const;
    bool isPointOk(size_t num, Interval::Position p) const;
    bool isPointDiameterOutlier(size_t num, Interval::Position p) const;

	Wm3::Vector3<double> getMaxPointInVolume(size_t num) const;
	Wm3::Vector3<double> getMaxPointInPlane(size_t num) const;

	Wm3::Vector3<double> getMinPointInVolume(size_t num, Interval::Position p) const;
	
	Wm3::Vector3<double> get50PointInVolume(size_t num, Interval::Position p) const;
	Wm3::Vector3<double> get50PointInPlane(size_t num, Interval::Position p) const;

    Wm3::Vector3<double> get20PointInVolume(size_t num, Interval::Position p) const;
    Wm3::Vector3<double> get20PointInPlane(size_t num, Interval::Position p) const;

    const Wm3::Vector3<double>& getImprovedPositionInVolume() const;
    const Wm3::Vector3<double>& getImprovedPositionInPlane() const;

    double getDiameter() const;
    double getDiameterVoxels() const;

	std::vector<Wm3::Vector3<double> > get50Points() const;
	double fitEllipse() const;
	void setCachedDiameter(double d);
	void setCachedSharpness(double s);
private:
	const Settings&      m_settings;
    double               m_pos;
    Data*                m_data;

    Wm3::Vector3<double> m_p;
    Wm3::Vector3<double> m_pOrig;
    Wm3::Vector3<double> m_betterPvolume;
    Wm3::Vector3<double> m_betterPplane;

    Wm3::Vector3<double> m_x;
    Wm3::Vector3<double> m_y;
    Wm3::Vector3<double> m_n;

    std::vector<Normal>  m_normals;

	double               m_cachedDiameter;
	double               m_cachedSharpness;

    void sampleNormals();
    void evaluateSharpness();

	double avgSharpnessHelper(const std::vector<size_t>& list) const;
	void stdSharpnessHelper(const std::vector<size_t>& list, double& mean, double& std) const;
	double stdSharpnessHelper(const std::vector<size_t>& list, double mean) const;
	double diameterHelper(const std::vector<size_t>& list) const;
	double medDiameterHelper(const std::vector<size_t>& list) const;
	double stdDiameterHelper(const std::vector<size_t>& list, double mean) const;
	Wm3::Vector3<double> medMaxHelper(const std::vector<size_t>& list) const;
	Wm3::Vector3<double> stdMaxHelper(const std::vector<size_t>& list, const Wm3::Vector3<double>& mean) const;
	
	double getAverageSharpnessChecked(std::vector<size_t>& checkedList) const;
    double getStdSharpnessChecked() const;
	double getDiameterChecked() const;
};

#endif // MEASUREMENTPOINT_H

