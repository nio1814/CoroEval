// Geometric Tools, Inc.
// http://www.geometrictools.com
// Copyright (c) 1998-2006.  All Rights Reserved
//
// The Wild Magic Library (WM3) source code is supplied under the terms of
// the license agreement
//     http://www.geometrictools.com/License/WildMagic3License.pdf
// and may not be copied or disclosed except in accordance with the terms
// of that agreement.

#ifndef WM3NONITERATIVEEIGEN3X3_H
#define WM3NONITERATIVEEIGEN3X3_H

#include "Wm3FoundationLIB.h"
#include "Wm3Matrix3.h"

namespace Wm3
{

template <class Real>
class WM3_ITEM NoniterativeEigen3x3
{
public:
    // The input matrix A must be symmetric.
    NoniterativeEigen3x3 (const Matrix3<Real>& rkA);
    ~NoniterativeEigen3x3 ();

    // Get the eigenvalues and eigenvectors.  The eigenvalues are stored in
    // increasing order.
    const Real GetEigenvalue (int i) const;
    const Real* GetEigenvalues () const;
    const Vector3<Real>& GetEigenvector (int i) const;
    const Vector3<Real>* GetEigenvectors () const;

private:
    // Compute the roots of the cubic polynomial.  Double precision arithmetic
    // is used to minimize the effects due to subtractive cancellation.  The
    // roots are returned in increasing order.
    void ComputeRoots (const Matrix3<Real>& rkA, double adRoot[3]);

    // Determine if M has positive rank.  The maximum-magnitude entry of M is
    // returned.  The row in which it is contained is also returned.
    bool PositiveRank (Matrix3<Real>& rkM, Real& rfMax,
        Vector3<Real>& rkMaxRow) const;

    // Compute the eigenvectors.
    void ComputeVectors (const Matrix3<Real>& rkA, Vector3<Real>& rkU2,
        int i0, int i1, int i2);

    Real m_afEigenvalue[3];
    Vector3<Real> m_akEigenvector[3];

    // For use by ComputeRoots.
    static const double ms_dInv3;
    static const double ms_dRoot3;
};

typedef NoniterativeEigen3x3<float> NoniterativeEigen3x3f;
typedef NoniterativeEigen3x3<double> NoniterativeEigen3x3d;

}

#endif
