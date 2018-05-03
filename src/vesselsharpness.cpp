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

#include "vesselsharpness.h"

#include "settings.h"
#include "os.h"

#define cimg_display 0
#include "CImg.h"

#include <limits>
#include <cmath>
#include <deque>
#include <stdexcept>


#undef max


VesselSharpness::VesselSharpness()
	:	m_settings(Settings::getSettings()),
		m_isValid(false) {}

bool VesselSharpness::isValid() const {
	return m_isValid;
}

void VesselSharpness::setNormal(const std::vector<float>& normal, double voxelSize)
{
    m_normal = normal;
	m_voxelSize = voxelSize;
    findExtrema();
    evaluateSharpness();
}

double VesselSharpness::getSharpness(Interval::Position p) const
{
	if (!m_isValid)
		return 0;

    switch(p)
    {
    case Interval::Left:
        return m_lSharpness;
    case Interval::Right:
        return m_rSharpness;
    case Interval::Both:
        int cnt = 0;
        double sharpness = 0;

        if(m_lSharpness != -1)
        {
            sharpness += m_lSharpness;
            cnt++;
        }

        if(m_rSharpness != -1)
        {
            sharpness += m_rSharpness;
            cnt++;
        }

        if( cnt > 1)
            return sharpness / (double)cnt;

        return 0;
    }

	return 0;
}

const Interval& VesselSharpness::getCenterInterval() const
{
    return m_center;
}

double VesselSharpness::getCenterDistance() const
{
	if (!m_isValid)
		return 0;

    double center = double(m_normal.size()) / 2.0;
    return m_center.getMaximumIndex() - center;
}

static const size_t kernelSize = 5;
static const size_t kernelHalfSize = 2;
static const double gaussKernel[5] = {
	// sigma = 1.0
	0.0544886845496429,
	0.244201342003233,
	0.402619946894247,
	0.244201342003233,
	0.0544886845496429
};

static int sgn(float n) {
	if (n > 0) 
		return 1;
	if (n < 0)
		return -1;
	return 0;
}

void VesselSharpness::findExtrema()
{
	// Prepare a filtered version of the original normal
	std::vector<float> normalFilt;
	normalFilt.resize(m_normal.size());
	
	for (size_t i = 0; i < kernelHalfSize; ++i)
		normalFilt[i] = m_normal[i];
	for (size_t i = m_normal.size() - kernelHalfSize; i < m_normal.size(); ++i)
		normalFilt[i] = m_normal[i];

	for (size_t i = kernelHalfSize; i < m_normal.size() - kernelHalfSize; ++i) {
		for (int j = -static_cast<int>(kernelHalfSize); j <= static_cast<int>(kernelHalfSize); ++j)
			normalFilt[i] += gaussKernel[kernelHalfSize + j] * m_normal[i + j];
	}

	// Find extrema
    std::deque<size_t> maxPoints, minPoints;
	size_t plateauP = 0;

    float prevDiff = normalFilt.at(0) - normalFilt.at(1);
    int prevDir = sgn(prevDiff);
	if (std::abs(prevDiff) < m_settings.getPlateauThreshold())
		prevDir = 0;
	
	int plateauDir = 1;
	bool plateau = false;
	if (prevDir == 0)
		plateau = true;

	float globalMax = 0;
	size_t globalMaxPos = 0;

	for (size_t i = 0; i < 2; i++) {
		if (m_normal[i] > globalMax) {
			globalMax = m_normal[i];
			globalMaxPos = i;
		}
	}

	for(size_t i=2; i<m_normal.size(); i++)
    {
        if (m_normal[i] > globalMax) {
			globalMax = m_normal[i];
			globalMaxPos = i;
		}
		
		float curDiff = normalFilt.at(i-1) - normalFilt.at(i);
		int curDir = sgn(curDiff);

		if (std::abs(curDiff) < m_settings.getPlateauThreshold()) {
			if (!plateau) {
				plateauDir = prevDir;
				prevDir = 0;
				plateauP = i - 1;
				plateau = true;
			}
			continue;
		} else if (plateau) {
			if (plateauDir > 0)
				minPoints.push_back(plateauP);
			else if (plateauDir < 0)
				maxPoints.push_back(plateauP);

			if (curDir < 0)
				minPoints.push_back(i - 1);
			else if (curDir > 0)
				maxPoints.push_back(i - 1);

			plateau = false;
			continue;
		}


        if ((prevDir > 0) && (curDir < 0)) {
            minPoints.push_back(i - 1);
        } else if ((prevDir < 0) && (curDir > 0)) {
            maxPoints.push_back(i - 1);
        }

        prevDir = curDir;
    }

	if (minPoints.empty() || maxPoints.empty()) {
		m_isValid = false;
		return;
	}

	// Correct position of extrema
	for (size_t i = 0; i < minPoints.size(); ++i) {
		size_t start = static_cast<size_t>(std::max<int>(0, static_cast<int>(minPoints[i]) - static_cast<int>(kernelHalfSize)));
		size_t end   = std::min<size_t>(m_normal.size() - 1, minPoints[i] + kernelHalfSize);
		float val = m_normal[start];
		size_t pos = start;

		for (size_t j = start; j <= end; ++j) {
			if (m_normal[j] < val) {
				val = m_normal[j];
				pos = j;
			}
		}

		minPoints[i] = pos;
	}

	for (size_t i = 0; i < maxPoints.size(); ++i) {
		size_t start = static_cast<size_t>(std::max<int>(0, static_cast<int>(maxPoints[i]) - static_cast<int>(kernelHalfSize) - 1));
		size_t end   = std::min<size_t>(m_normal.size() - 1, maxPoints[i] + kernelHalfSize + 1);
		float val = m_normal[start];
		size_t pos = start;

		for (size_t j = start; j <= end; ++j) {
			if (m_normal[j] > val) {
				val = m_normal[j];
				pos = j;
			}
		}

		maxPoints[i] = pos;
	}

	// If the global maximum is not in the list, add it
	bool globalMaxFound = false;
	for (size_t i = 0; i < maxPoints.size(); i++) {
    if (std::abs(maxPoints[i] - globalMaxPos) < kernelSize) {
			globalMaxFound = true;
			break;
		}
		if (std::abs(m_normal[maxPoints[i]] - globalMax) <= m_settings.getPlateauThreshold()) {
			globalMaxFound = true;
			break;
		}
	}

	if (!globalMaxFound) {
		bool maxIn = false;

		for (size_t i = 0; i < maxPoints.size(); i++) {
			if (maxPoints[i] > globalMaxPos) {
				maxPoints.insert(maxPoints.begin() + i, globalMaxPos);
				maxIn = true;
				break;
			}
		}

		if (!maxIn)
			maxPoints.push_back(globalMaxPos);
	}

	// Normal has to start with a minimum
    if( minPoints.front() > maxPoints.front() )
    {
        while( minPoints.front() > maxPoints.front() ) {
            maxPoints.pop_front();
			if (maxPoints.empty()) {
				m_isValid = false;
				return;
			}
		}
    }

	// We need a minimum right of the last maximum
	if (maxPoints.back() > minPoints.back()) {
		size_t minPos = maxPoints.back() + 1;

		if (minPos < m_normal.size()) {
			float minVal = m_normal[minPos];

			for (size_t pos = maxPoints.back() + 1; pos < m_normal.size(); pos++) {
				if (m_normal[pos] < minVal) {
					minVal = m_normal[pos];
					minPos = pos;
				}
			}

			minPoints.push_back(minPos);
		}
	}

	// Between two minima, there needs to be a maximum
	for (size_t i = 1; i < minPoints.size(); ++i) {
		size_t minAPos = minPoints[i - 1];
		size_t minBPos = minPoints[i];

		bool found = false;
		for (size_t j = 0; j < maxPoints.size(); ++j) {
			if ((maxPoints[j] > minAPos) && (maxPoints[j] < minBPos)) {
				found = true;
				break;
			}
		}

		if (found)
			continue;

		size_t maxPos = minAPos + 1;
		float maxVal = m_normal[maxPos];
		for (size_t j = maxPos; j < minBPos; ++j) {
			if (m_normal[j] > maxVal) {
				maxVal = m_normal[j];
				maxPos = j;
			}
		}

		bool maxIn = false;

		for (size_t j = 0; j < maxPoints.size(); j++) {
			if (maxPoints[j] > maxPos) {
				maxPoints.insert(maxPoints.begin() + j, maxPos);
				maxIn = true;
				break;
			}
		}

		if (!maxIn)
			maxPoints.push_back(maxPos);
	}

	size_t minDist = static_cast<size_t>(round(m_settings.getMinVesselDiameter() / m_voxelSize)); // / 2 for radius, * 2 for subsampling = * 1
	
	// One exception to the rule above: If right of the last minimum the profile continuously decreases,
	// we need a minimum at the end of the profile
	if ((m_normal.size() - minPoints.back() > minDist + 1) &&
		(m_normal[minPoints.back()] > m_normal.back())) {
		minPoints.push_back(m_normal.size() - 1);
	}

	// Attach adjacent minima to maxima
	std::vector<Interval> allIntervals;
    for(size_t i=0; i<maxPoints.size(); i++)
    {
        if (maxPoints[i] < minDist)
			continue;

		Interval curInt;
        curInt.setMaximum(maxPoints.at(i),       m_normal.at(maxPoints.at(i)));
		bool found = false;
        for(size_t j=0; j<minPoints.size(); j++)
        {
            if (minPoints.at(j) > maxPoints.at(i) + minDist)
            {
                curInt.setMinimum(minPoints.at(j), m_normal.at(minPoints.at(j)), Interval::Right );
                
				bool noLeft = false;
        for (int k = j - 1; k >= 0; --k) {
					if (minPoints[k] < maxPoints[i] - minDist) {
						curInt.setMinimum(minPoints[k], m_normal.at(minPoints[k]), Interval::Left );
						break;
					}
					if (k == 0) {
						noLeft = true;
						break;
					}
				}
				
				if (noLeft)
					break;

				found = true;
                break;
            }
        }
        if(found)
        {
            allIntervals.push_back(curInt);
        }
    }

	if (allIntervals.empty()) {
		m_isValid = false;
		return;
	}

	// Check and correct for local minima
	for (size_t i = 0; i < allIntervals.size(); i++) {
		Interval& curInt = allIntervals[i];
		
		double magnL = curInt.getMagnitude(Interval::Left);
        double magnR = curInt.getMagnitude(Interval::Right);
		double magnThr = std::max<double>(magnL, magnR) * m_settings.getMinRelMagniude();

		if ((magnL < magnThr) &&
			(curInt.getMinimumIndex(Interval::Left) > minPoints[0])) {
			size_t minIdx = 0;
			
			while ((minIdx < minPoints.size()) && (curInt.getMinimumIndex(Interval::Left) != minPoints[minIdx]))
				++minIdx;
			
			if (minIdx >= minPoints.size())
				throw std::runtime_error("Something is very wrong");

			while ((minIdx > 0) && (curInt.getMagnitude(Interval::Left) < magnThr)) {
				--minIdx;
				curInt.setMinimum(minPoints[minIdx], m_normal[minPoints[minIdx]], Interval::Left);
			}
		} else if ((magnR < magnThr) &&
		           (curInt.getMinimumIndex(Interval::Right) < minPoints[minPoints.size() - 1])) {
			size_t minIdx = 0;
			
			while ((minIdx < minPoints.size()) && (curInt.getMinimumIndex(Interval::Right) != minPoints[minIdx]))
				++minIdx;
			
			if (minIdx >= minPoints.size())
				throw std::runtime_error("Something is very wrong");

			while ((minIdx < minPoints.size() - 1) && (curInt.getMagnitude(Interval::Right) < magnThr)) {
				++minIdx;
				curInt.setMinimum(minPoints[minIdx], m_normal[minPoints[minIdx]], Interval::Right);
			}
		}
	}

	// Find maximum closest to centre
	double center = double(m_normal.size()) / 2.0;
    size_t centerMax = 0;
    double centerDist = center;
    for(size_t i=0; i<allIntervals.size(); i++)
    {
        double dist = qAbs(double(allIntervals.at(i).getMaximumIndex()) - center);
        if(dist < centerDist)
        {
            centerDist = dist;
            centerMax  = i;
        }
    }

    m_center = allIntervals.at(centerMax);

    if (allIntervals.size() == 1) {
		m_isValid = true;
		return;
	}

	// Compute the magnitude for all possible intervals
    double meanMagn = 0;
	double minMagn  = std::numeric_limits<double>::max();
    double maxMagn  = 0;
	    
    for(size_t i=0; i<allIntervals.size(); i++)
    {
        double magnL = allIntervals.at(i).getMagnitude(Interval::Left);
        double magnR = allIntervals.at(i).getMagnitude(Interval::Right);
        meanMagn += magnL;
        meanMagn += magnR;

        minMagn = std::min<double>( minMagn, magnL );
        minMagn = std::min<double>( minMagn, magnR );
        maxMagn = std::max<double>( maxMagn, magnL );
        maxMagn = std::max<double>( maxMagn, magnR );
    }
    meanMagn /= (double)(allIntervals.size()*2);

    double stdMagn = 0;
    for(size_t i=0; i<allIntervals.size(); i++)
    {
        double magnL = allIntervals.at(i).getMagnitude(Interval::Left);
        double val = magnL - meanMagn;
        stdMagn += (val*val);

        double magnR = allIntervals.at(i).getMagnitude(Interval::Right);
        val = magnR - meanMagn;
        stdMagn += (val*val);
    }
    stdMagn = sqrt(stdMagn / (double)(2*allIntervals.size() - 1));

    // Correct - if neccessary
    double thresMagn = minMagn + 1.0 * stdMagn;
    for(size_t j=0; centerMax-j > 0 && m_center.getMagnitude(Interval::Left) < thresMagn; j++ )
    {
        const Interval& prevInterval = allIntervals.at( centerMax-j );
        m_center.setMinimum( prevInterval.getMinimumIndex(Interval::Left),
                             prevInterval.getMinimumValue(Interval::Left),
                             Interval::Left );

        if( m_center.getMaximumValue() < prevInterval.getMaximumValue() )
        {
            m_center.setMaximum(prevInterval.getMaximumIndex(), prevInterval.getMaximumValue());
        }
    }

    for(size_t j=0; centerMax+j < allIntervals.size() && m_center.getMagnitude(Interval::Right) < thresMagn; j++ )
    {
        const Interval& nextInterval = allIntervals.at( centerMax+j );
        m_center.setMinimum( nextInterval.getMinimumIndex(Interval::Right),
                             nextInterval.getMinimumValue(Interval::Right),
                             Interval::Right );
        if( m_center.getMaximumValue() < nextInterval.getMaximumValue() )
        {
            m_center.setMaximum(nextInterval.getMaximumIndex(), nextInterval.getMaximumValue());
        }
    }

	m_isValid = true;

    return;
    std::cout << "Found vessel at " << (m_center.getMaximumIndex() - center)/2 << " between "
              << (m_center.getMinimumIndex(Interval::Left) - center)/2 << " and " << (m_center.getMinimumIndex(Interval::Right) - center)/2
              << " with a magnitude of " << m_center.getMagnitude(Interval::Left) << " / " << m_center.getMagnitude(Interval::Right) << std::endl;
}

void VesselSharpness::evaluateSharpness()
{
	if (!m_isValid)
		return;

    cimg_library::CImg<float> line( &(m_normal[0]), m_normal.size(), 1);

    m_center.print();

    double maxInd  = double(m_center.getMaximumIndex());
    double lMinInd = double(m_center.getMinimumIndex(Interval::Left));
    double rMinInd = double(m_center.getMinimumIndex(Interval::Right));

    double ind80l, ind50l, ind20l, thr80l, thr50l, thr20l;
    thr80l = m_center.get80Value(Interval::Left);
	thr50l = m_center.get50Value(Interval::Left);
    thr20l = m_center.get20Value(Interval::Left);
    for(double x=lMinInd; x<maxInd; x += 0.01 )
    {
        float val = line.linear_atX(x);
        if( val < thr20l )
        {
            ind20l = x;
        }
		if( val < thr50l )
        {
            ind50l = x;
        }
        if( val < thr80l )
        {
            ind80l = x;
        }
    }

    if( (ind80l - ind20l) > 0 )
    {
        m_lSharpness = 1.0 / (ind80l - ind20l);
		m_lSharpness /= 2 * m_voxelSize; // 2x oversampling of normal!
        m_center.set20Index(ind20l, Interval::Left);
		m_center.set50Index(ind50l, Interval::Left);
        m_center.set80Index(ind80l, Interval::Left);
    } else {
        m_lSharpness = -1;
    }

    double  ind80r, ind50r, ind20r, thr80r, thr50r, thr20r;
    thr80r = m_center.get80Value(Interval::Right);
	thr50r = m_center.get50Value(Interval::Right);
    thr20r = m_center.get20Value(Interval::Right);
    for(double x=rMinInd; x>maxInd; x -= 0.01 )
    {
        float val = line.linear_atX(x);
        if( val < thr20r )
        {
            ind20r = x;
        }
		if( val < thr50r )
        {
            ind50r = x;
        }
        if( val < thr80r )
        {
            ind80r = x;
        }
    }
    if( (ind20r - ind80r) > 0 )
    {
        m_rSharpness = 1.0 / (ind20r - ind80r);
		m_rSharpness /= 2 * m_voxelSize; // 2x oversampling of normal!
        m_center.set20Index(ind20r, Interval::Right);
		m_center.set50Index(ind50r, Interval::Right);
        m_center.set80Index(ind80r, Interval::Right);
    } else {
        m_rSharpness = -1;
    }
}

