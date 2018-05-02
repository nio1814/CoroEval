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

#include "MPFilter.h"

#include "measurementpoint.h"
#include "bsplineinterpolation.h"


void MPFilter::median(std::vector<MeasurementPoint*>& mp, size_t len) {
	if (len % 2 == 0)
		++len;

	size_t lenHalf = (len - 1) / 2;

	if (mp.size() < len)
		return;

	std::vector<double> newD;
	std::vector<double> newS;
	
	for (size_t i = lenHalf; i < mp.size() - lenHalf; ++i) {
		std::vector<double> d;
		std::vector<double> s;

		for (size_t j = lenHalf; j >= 1; --j) {
			d.push_back(mp[i - j]->getDiameter());
			s.push_back(mp[i - j]->getAverageSharpness());
		}

		d.push_back(mp[i]->getDiameter());
		s.push_back(mp[i]->getAverageSharpness());

		for (size_t j = 1; j <= lenHalf; ++j) {
			d.push_back(mp[i + j]->getDiameter());
			s.push_back(mp[i + j]->getAverageSharpness());
		}

		std::sort(d.begin(), d.end());
		std::sort(s.begin(), s.end());

		newD.push_back(d[lenHalf]);
		newS.push_back(s[lenHalf]);
	}

	size_t j = 0;
	for (size_t i = lenHalf; i < mp.size() - lenHalf; ++i, ++j) {
		mp[i]->setCachedDiameter(newD[j]);
		mp[i]->setCachedSharpness(newS[j]);
	}
}

void MPFilter::gauss(std::vector<MeasurementPoint*>& mp, size_t len, double sigma) {
	if (len % 2 == 0)
		++len;

	size_t lenHalf = (len - 1) / 2;

	if (mp.size() < len)
		return;

	std::vector<double> newD;
	std::vector<double> newS;
	
	double* gausswin = BSplineInterpolation::gaussWin(len, sigma);

	for (size_t i = lenHalf; i < mp.size() - lenHalf; ++i) {
		size_t gI = 0;
		double d = 0;
		double s = 0;

		for (size_t j = lenHalf; j >= 1; --j, ++gI) {
			d += mp[i - j]->getDiameter() * gausswin[gI];
			s += mp[i - j]->getAverageSharpness() * gausswin[gI];
		}

		d += mp[i]->getDiameter() * gausswin[gI];
		s += mp[i]->getAverageSharpness() * gausswin[gI++];

		for (size_t j = 1; j <= lenHalf; ++j, ++gI) {
			d += mp[i + j]->getDiameter() * gausswin[gI];
			s += mp[i + j]->getAverageSharpness() * gausswin[gI];
		}

		newD.push_back(d);
		newS.push_back(s);
	}

	delete[] gausswin;

	size_t j = 0;
	for (size_t i = lenHalf; i < mp.size() - lenHalf; ++i, ++j) {
		mp[i]->setCachedDiameter(newD[j]);
		mp[i]->setCachedSharpness(newS[j]);
	}
}

void MPFilter::purge0(std::vector<MeasurementPoint*>& mp) {
	bool in = false;
	size_t start = 0;
	size_t end = 0;

	for (size_t i = 0; i < mp.size(); ++i) {
		if (mp[i]->getDiameter() < 1e-3) {
			if (!in) {
				start = i;
				in = true;
			}
		} else if (in) {
			in = false;
			end = i - 1;

			double startVal;
			if (start == 0)
				startVal = mp[i]->getDiameter();
			else
				startVal = mp[start - 1]->getDiameter();

			double endVal = mp[i]->getDiameter();

			double inc = (endVal - startVal) / (end - start + 2);
			size_t pos = 1;
			for (size_t j = start; j <= end; ++j)
				mp[j]->setCachedDiameter(startVal + inc * static_cast<double>(pos++));

			std::cerr << "Purged " << (end - start + 1) << " 0s" << std::endl << std::flush;
		}
	}

	if (in) {
		if (start == 0)
			return;

		double val = mp[start - 1]->getDiameter();

		for (size_t j = start; j < mp.size(); ++j)
			mp[j]->setCachedDiameter(val);

		std::cerr << "Purged " << (mp.size() - start) << " 0s" << std::endl << std::flush;
	}
}

