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

#ifndef Data_h
#define Data_h

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <sstream>
#include "wm3/Wm3Vector3.h"
#include "wm3/Wm3Matrix3.h"

#define cimg_display 0
#include "CImg.h"

#ifndef MIN
#define MIN(a,b) (a < b ? a : b)
#endif

class Data 
{
public:
    enum Dim {
        COL,
        LIN,
        PAR
    };

    Data();

    Data(float* data);
    ~Data() {}

    void clone(const Data* other);
    void create(const Data& other);
    void create(Dim dim1, long size1, Dim dim2, long size2, Dim dim3=PAR, long size3=1);
    void create(std::map<Dim, long> dims);

    long getLen(Dim d) const;
    void setLen(Dim d, long newLen);

    long size() const;
    long sliceSize() const;

    void setOrientation(Wm3::Matrix3<double> rotMat);
    Wm3::Matrix3<double> getOrientation();

    void   setVoxelSize(double val);
    double getVoxelSize() const;

    void preset(float val = 0);

    double getCurrentPar() const;
    void setPar(double newParIdx);

    void setData(float* newData) { m_data = newData; initVolume(); }

    float* data0() const;
    float getValue(const Wm3::Vector3<double> p) const;
    float getValue(const double x, const double y, const double z) const;

protected:

    bool hasSameSize(const Data &other) const;
    bool hasSameSliceSize(const Data &other) const;

    Wm3::Matrix3<double>        m_rotMat;
    cimg_library::CImg<float>*  m_volume;
    float*                      m_data;
    double                      m_pxSize;
    double                        m_idx_PAR;
    std::map<Dim, long>         m_dim2;

private:
    void initVolume();
};

inline bool Data::hasSameSize(const Data &other) const
{
    if (this->size()==other.size()) return true;
    else return false;
}

inline bool Data::hasSameSliceSize(const Data &other) const
{
    if (this->sliceSize()==other.sliceSize()) return true;
    else return false;
}

#endif

