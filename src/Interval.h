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

#ifndef INTERVAL_H
#define INTERVAL_H

#include <algorithm>
#include <iostream>

class Interval {

public:
    enum Position{
        Left,
        Right,
        Both
    };

    Interval();

    void    setMaximum(size_t ind, double val);
    void    setMinimum(size_t ind, double val, Interval::Position p);
    void    set20Index(double ind, Interval::Position p);
	void    set50Index(double ind, Interval::Position p);
    void    set80Index(double ind, Interval::Position p);

    double  getMagnitude(Interval::Position p) const;

    size_t  getMaximumIndex() const;
    double  getMaximumValue() const;

    size_t  getMinimumIndex(Interval::Position p) const;
    double  getMinimumValue(Interval::Position p) const;

    double  get80Index(Interval::Position p) const;
    double  get80Value(Interval::Position p) const;

	double  get50Index(Interval::Position p) const;
    double  get50Value(Interval::Position p) const;

    double  get20Index(Interval::Position p) const;
    double  get20Value(Interval::Position p) const;

    void    print() const;

private:
    double m_lMinVal;
    double m_maxVal;
    double m_rMinVal;

    size_t m_lMinInd;
    size_t m_maxInd;
    size_t m_rMinInd;

    double m_l20Ind;
	double m_l50Ind;
    double m_l80Ind;
    double m_r20Ind;
	double m_r50Ind;
    double m_r80Ind;
};
#endif // INTERVAL_H

