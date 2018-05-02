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

#ifndef MPFILTER_H
#define MPFILTER_H

#include <vector>
#include <cstdlib>

// Forward declarations
class MeasurementPoint;


class MPFilter
{
private:
    MPFilter() {}

public:
	static void median(std::vector<MeasurementPoint*>& mp, size_t len = 3);
	static void gauss(std::vector<MeasurementPoint*>& mp, size_t len = 3, double sigma = 0.5);
	static void purge0(std::vector<MeasurementPoint*>& mp);
};

#endif // MPFILTER_H

