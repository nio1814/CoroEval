// Geometric Tools, Inc.
// http://www.geometrictools.com
// Copyright (c) 1998-2006.  All Rights Reserved
//
// The Wild Magic Library (WM3) source code is supplied under the terms of
// the license agreement
//     http://www.geometrictools.com/License/WildMagic3License.pdf
// and may not be copied or disclosed except in accordance with the terms
// of that agreement.

#ifndef WM3BSPLINEBASIS_H
#define WM3BSPLINEBASIS_H

#include "Wm3FoundationLIB.h"
#include "Wm3Math.h"

namespace Wm3
{

template <class Real>
class WM3_ITEM BSplineBasis
{
public:
    BSplineBasis ();

    // Open uniform or periodic uniform.  The knot array is internally
    // generated with equally spaced elements.
    BSplineBasis (int iNumCtrlPoints, int iDegree, bool bOpen);
    void Create (int iNumCtrlPoints, int iDegree, bool bOpen);

    // Open nonuniform.  The knot array must have n-d elements.  The elements
    // must be nondecreasing.  Each element must be in [0,1].  The caller is
    // responsible for deleting afKnot.  An internal copy is made, so to
    // dynamically change knots you must use the SetKnot function.
    BSplineBasis (int iNumCtrlPoints, int iDegree, const Real* afKnot);
    void Create (int iNumCtrlPoints, int iDegree, const Real* afKnot);

    virtual ~BSplineBasis ();

    int GetNumCtrlPoints () const;
    int GetDegree () const;
    bool IsOpen () const;
    bool IsUniform () const;

    // The knot values can be changed only if the basis function is nonuniform
    // and the input index is valid (0 <= i <= n-d-1).  If these conditions
    // are not satisfied, GetKnot returns MAX_REAL.
    void SetKnot (int i, Real fKnot);
    Real GetKnot (int i) const;

    // access basis functions and their derivatives
    Real GetD0 (int i) const;
    Real GetD1 (int i) const;
    Real GetD2 (int i) const;
    Real GetD3 (int i) const;

    // evaluate basis functions and their derivatives
    void Compute (Real fTime, unsigned int uiOrder, int& riMinIndex,
        int& riMaxIndex) const;

protected:
    int Initialize (int iNumCtrlPoints, int iDegree, bool bOpen);
    Real** Allocate () const;
    void Deallocate (Real** aafArray);

    // Determine knot index i for which knot[i] <= rfTime < knot[i+1].
    int GetKey (Real& rfTime) const;

    int m_iNumCtrlPoints;    // n+1
    int m_iDegree;           // d
    Real* m_afKnot;          // knot[n+d+2]
    bool m_bOpen, m_bUniform;

    // Storage for the basis functions and their derivatives first three
    // derivatives.  The basis array is always allocated by the constructor
    // calls.  A derivative basis array is allocated on the first call to a
    // derivative member function.
    Real** m_aafBD0;             // bd0[d+1][n+d+1]
    mutable Real** m_aafBD1;     // bd1[d+1][n+d+1]
    mutable Real** m_aafBD2;     // bd2[d+1][n+d+1]
    mutable Real** m_aafBD3;     // bd3[d+1][n+d+1]
};

typedef BSplineBasis<float> BSplineBasisf;
typedef BSplineBasis<double> BSplineBasisd;

}

#endif

