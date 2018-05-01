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

#include "Data.h"


Data::Data()
{
    m_data=0;
    m_volume=0;
    m_idx_PAR=0;
    m_rotMat = Wm3::Matrix3<double>::IDENTITY;
}

Data::Data(float* data)
{
    m_volume=0;
    m_idx_PAR=0;
    m_rotMat = Wm3::Matrix3<double>::IDENTITY;

    m_data = data;
}

void Data::clone(const Data* other)
{
    create(*other);
    m_data = other->data0();
    initVolume();
}

void Data::create(const Data& other)
{
    if (m_dim2.size() > 0) m_dim2.clear();
    m_dim2.insert( std::pair<Dim, long>(COL, other.getLen(COL)) );
    m_dim2.insert( std::pair<Dim, long>(LIN, other.getLen(LIN)) );
    m_dim2.insert( std::pair<Dim, long>(PAR, other.getLen(PAR)) );

    m_idx_PAR=0;
    m_rotMat = Wm3::Matrix3<double>::IDENTITY;
}

void Data::create(std::map<Dim, long> dims)
{
    if (m_dim2.size() > 0) m_dim2.clear();
    m_dim2 = dims;
    m_rotMat = Wm3::Matrix3<double>::IDENTITY;

    initVolume();
}

void Data::create(Dim dim1, long size1, Dim dim2, long size2, Dim dim3, long size3)
{
    m_dim2.clear();
    m_dim2.insert( std::pair<Dim, long>(dim1, size1) );
    m_dim2.insert( std::pair<Dim, long>(dim2, size2) );
    m_dim2.insert( std::pair<Dim, long>(dim3, size3) );
    m_rotMat = Wm3::Matrix3<double>::IDENTITY;

    initVolume();
}

float* Data::data0() const
{
    return m_data;
}

long Data::getLen(Dim d) const
{
    if( const_cast<Data&>(*this).m_dim2[d] == 0)
        return 1;

    if(m_rotMat == Wm3::Matrix3<double>::IDENTITY)
        return const_cast<Data&>(*this).m_dim2[d];

    Wm3::Vector3<double> dim = Wm3::Vector3<double>(const_cast<Data&>(*this).m_dim2[COL],const_cast<Data&>(*this).m_dim2[LIN],const_cast<Data&>(*this).m_dim2[PAR]) * m_rotMat;

    long ret;
    switch(d)
    {
    case COL:
        ret = static_cast<long>(dim.X());
        break;
    case LIN:
        ret = static_cast<long>(dim.Y());
        break;
    case PAR:
        ret = static_cast<long>(dim.Z());
        break;
    }
    return ret;
}

void Data::setLen(Dim d, long newLen)
{
    m_dim2[d] = newLen;

    initVolume();
}

void Data::setOrientation(Wm3::Matrix3<double> rotMat)
{
    m_rotMat = rotMat;
}

Wm3::Matrix3<double> Data::getOrientation()
{
    return m_rotMat;
}

void Data::setVoxelSize(double val)
{
    if( val < 0.0 )
        return;
    m_pxSize = val;
}

double Data::getVoxelSize() const
{
    return m_pxSize;
}

long Data::size() const
{
    return getLen(LIN)*getLen(COL)*getLen(PAR);
}

long Data::sliceSize() const
{
    return getLen(LIN)*getLen(COL);
}

void Data::preset(float val)
{
    long N = size();
    for (long i=0; i<N; i++)
    {
        m_data[i] = val;
    }
}

double Data::getCurrentPar() const
{
    return m_idx_PAR;
}

void Data::setPar(double newParIdx)
{
    m_idx_PAR = MIN(getLen(Data::PAR)-1, newParIdx);
}

float Data::getValue(const double x, const double y, const double z) const
{
    return getValue( Wm3::Vector3<double>(x,y,z) );
}

float Data::getValue(const Wm3::Vector3<double> p) const
{
    if( !m_volume )
        return 0;

    Wm3::Vector3<double> realP = p * m_rotMat;
    if( realP.X() < 0 || realP.X() > const_cast<Data&>(*this).m_dim2[COL] )
        return 0.f;
    if( realP.Y() < 0 || realP.Y() > const_cast<Data&>(*this).m_dim2[LIN] )
        return 0.f;
    if( realP.Z() < 0 || realP.Z() > const_cast<Data&>(*this).m_dim2[PAR] )
        return 0.f;
    return m_volume->linear_atXYZ(static_cast<float>(realP.X()), static_cast<float>(realP.Y()), static_cast<float>(realP.Z()));
}

void Data::initVolume()
{
    if( !data0() )
        return;

    if( m_dim2.empty() )
        return;

    if( m_volume )
        delete m_volume;

    m_volume = new cimg_library::CImg<float>( data0(), m_dim2[COL], m_dim2[LIN], m_dim2[PAR], 1, true);
}

