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

#include "Interval.h"


Interval::Interval()
{
    m_lMinVal = 0.0;
    m_maxVal = 0.0;
    m_rMinVal = 0.0;

    m_lMinInd = 0;
    m_maxInd = 0;
    m_rMinInd = 0;

    m_l20Ind = 0.0;
	m_l50Ind = 0.0;
    m_l80Ind = 0.0;
    m_r20Ind = 0.0;
	m_r50Ind = 0.0;
    m_r80Ind = 0.0;
}

void Interval::setMaximum(size_t ind, double val)
{
    m_maxInd = ind;
    m_maxVal = val;
}

void Interval::setMinimum(size_t ind, double val, Interval::Position p)
{
    switch(p)
    {
    case Interval::Left:
        m_lMinInd = ind;
        m_lMinVal = val;
        break;
    case Interval::Right:
        m_rMinInd = ind;
        m_rMinVal = val;
        break;
      case Both:
        break;
    }
}

void Interval::set20Index(double ind, Interval::Position p)
{
    switch(p)
    {
    case Interval::Left:
        m_l20Ind = ind;
        break;
    case Interval::Right:
        m_r20Ind = ind;
        break;
      case Both:
        break;
    }
}

void Interval::set50Index(double ind, Interval::Position p)
{
    switch(p)
    {
    case Interval::Left:
        m_l50Ind = ind;
        break;
    case Interval::Right:
        m_r50Ind = ind;
        break;
      case Both:
        break;
    }
}

void Interval::set80Index(double ind, Interval::Position p)
{
    switch(p)
    {
    case Interval::Left:
        m_l80Ind = ind;
        break;
    case Interval::Right:
        m_r80Ind = ind;
        break;
      case Both:
        break;
    }
}

double Interval::getMagnitude(Interval::Position p) const
{
    double res = 0;
    switch(p)
    {
    case Interval::Left:
        res = m_maxVal - m_lMinVal;
        break;
    case Interval::Right:
        res = m_maxVal - m_rMinVal;
        break;
    case Interval::Both:
        res = m_maxVal - std::min<double>(m_lMinVal, m_rMinVal);
        break;
    }
    return res;
}

size_t Interval::getMaximumIndex() const
{
    return m_maxInd;
}

double Interval::getMaximumValue() const
{
    return m_maxVal;
}

size_t Interval::getMinimumIndex(Interval::Position p) const
{
    size_t res = 0;
    switch(p)
    {
    case Interval::Left:
        res = m_lMinInd;
        break;
    case Interval::Right:
        res = m_rMinInd;
        break;
      case Both:
        res = std::min(m_lMinInd, m_rMinInd);
        break;
    }
    return res;
}

double Interval::getMinimumValue(Interval::Position p) const
{
    double res = 0;
    switch(p)
    {
    case Interval::Left:
        res = m_lMinVal;
        break;
    case Interval::Right:
        res = m_rMinVal;
        break;
      case Both:
        res = std::min(m_lMinVal, m_rMinVal);
        break;
    }
    return res;
 }

double Interval::get80Index(Interval::Position p) const
{
    double res = 0;
    switch(p)
    {
    case Interval::Left:
        res = m_l80Ind;
        break;
    case Interval::Right:
        res = m_r80Ind;
        break;
    case Both:
        break;
    }
    return res;
 }

double Interval::get50Index(Interval::Position p) const
{
    double res = 0;
    switch(p)
    {
    case Interval::Left:
        res = m_l50Ind;
        break;
    case Interval::Right:
        res = m_r50Ind;
        break;
    case Both:
        break;
    }
    return res;
 }

double Interval::get20Index(Interval::Position p) const
{
    double res = 0;
    switch(p)
    {
    case Interval::Left:
        res = m_l20Ind;
        break;
    case Interval::Right:
        res = m_r20Ind;
        break;
    case Both:
        break;
    }
    return res;
 }

double Interval::get80Value(Interval::Position p) const
{
    return m_maxVal - 0.2*getMagnitude(p);
}

double Interval::get50Value(Interval::Position p) const
{
    return m_maxVal - 0.5*getMagnitude(p);
}

double Interval::get20Value(Interval::Position p) const
{
    return m_maxVal - 0.8*getMagnitude(p);
 }

void Interval::print() const
{
    if( m_lMinInd > m_maxInd )
      std::cout << "Interval: [" << m_lMinInd << ", " << m_maxInd << ", " << m_rMinInd << "]" << std::endl;
}

