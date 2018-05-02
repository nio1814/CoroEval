// Geometric Tools, Inc.
// http://www.geometrictools.com
// Copyright (c) 1998-2006.  All Rights Reserved
//
// The Wild Magic Library (WM3) source code is supplied under the terms of
// the license agreement
//     http://www.geometrictools.com/License/WildMagic3License.pdf
// and may not be copied or disclosed except in accordance with the terms
// of that agreement.

#ifndef WM3BSPLINEFITBASIS_H
#define WM3BSPLINEFITBASIS_H

#include "Wm3FoundationLIB.h"
#include "Wm3System.h"

namespace Wm3
{

template <class Real>
class WM3_ITEM BSplineFitBasis
{
public:
    // Construction and destruction.  This class is only for open uniform
    // B-spline basis functions.  The input is the number of control points
    // for a B-spline curve using this basis and the degree of that curve.
    BSplineFitBasis (int iQuantity, int iDegree);
    ~BSplineFitBasis ();

    // Data member access.
    int GetQuantity () const;
    int GetDegree () const;

    // Evaluate the basis functions.  This function fills in the values
    // returned by GetValue(i) for 0 <= i <= degree.  The return indices iMin
    // and iMax are relative to the array of control points.  The GetValue(i)
    // are the coefficients for the control points ctrl[iMin] throught
    // ctrl[iMax] in the curve evaluation (i.e. the curve has local control).
    void Compute (Real fT, int& iMin, int& iMax) const;
    Real GetValue (int i) const;

private:
    // The number of control points and degree for the curve.
    int m_iQuantity, m_iDegree;

    // The storage for knots and basis evaluation.
    mutable Real* m_afValue;  // m_afValue[0..degree]
    mutable Real* m_afKnot;   // m_afKnot[2*degree]
};

typedef BSplineFitBasis<float> BSplineFitBasisf;
typedef BSplineFitBasis<double> BSplineFitBasisd;

}

#endif

