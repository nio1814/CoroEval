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

#ifndef VESSELSHARPNESS_H
#define VESSELSHARPNESS_H

#include "Interval.h"

#include <vector>


// Forward declarations
class Settings;


class VesselSharpness
{
public:
    VesselSharpness();

    void setNormal(const std::vector<float>& normal, double voxelSize);
    double getSharpness(Interval::Position p) const;
    const Interval& getCenterInterval() const;
    double getCenterDistance() const;
	bool isValid() const;

private:
	const Settings&    m_settings;
    std::vector<float> m_normal;
    Interval           m_center;

	double			   m_voxelSize;

    double             m_lSharpness;
    double             m_rSharpness;

	bool               m_isValid;

    void findExtrema();
    void evaluateSharpness();
};

#endif // VESSELSHARPNESS_H

