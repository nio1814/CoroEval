// Geometric Tools, Inc.
// http://www.geometrictools.com
// Copyright (c) 1998-2006.  All Rights Reserved
//
// The Wild Magic Library (WM3) source code is supplied under the terms of
// the license agreement
//     http://www.geometrictools.com/License/WildMagic3License.pdf
// and may not be copied or disclosed except in accordance with the terms
// of that agreement.

#include "Wm3BSplineBasis.h"

namespace Wm3
{
//----------------------------------------------------------------------------
template <class Real>
BSplineBasis<Real>::BSplineBasis ()
{
}
//----------------------------------------------------------------------------
template <class Real>
BSplineBasis<Real>::BSplineBasis (int iNumCtrlPoints, int iDegree, bool bOpen)
{
    Create(iNumCtrlPoints,iDegree,bOpen);
}
//----------------------------------------------------------------------------
template <class Real>
void BSplineBasis<Real>::Create (int iNumCtrlPoints, int iDegree, bool bOpen)
{
    m_bUniform = true;

    int i, iNumKnots = Initialize(iNumCtrlPoints,iDegree,bOpen);
    Real fFactor = ((Real)1.0)/(m_iNumCtrlPoints-m_iDegree);
    if (m_bOpen)
    {
        for (i = 0; i <= m_iDegree; i++)
        {
            m_afKnot[i] = (Real)0.0;
        }

        for (/**/; i < m_iNumCtrlPoints; i++)
        {
            m_afKnot[i] = (i-m_iDegree)*fFactor;
        }

        for (/**/; i < iNumKnots; i++)
        {
            m_afKnot[i] = (Real)1.0;
        }
    }
    else
    {
        for (i = 0; i < iNumKnots; i++)
        {
            m_afKnot[i] = (i-m_iDegree)*fFactor;
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
BSplineBasis<Real>::BSplineBasis (int iNumCtrlPoints, int iDegree,
    const Real* afKnot)
{
    Create(iNumCtrlPoints,iDegree,afKnot);
}
//----------------------------------------------------------------------------
template <class Real>
void BSplineBasis<Real>::Create (int iNumCtrlPoints, int iDegree,
    const Real* afKnot)
{
    m_bUniform = false;

    int i, iNumKnots = Initialize(iNumCtrlPoints,iDegree,true);
    for (i = 0; i <= m_iDegree; i++)
    {
        m_afKnot[i] = (Real)0.0;
    }

    for (int j = 0; i < m_iNumCtrlPoints; i++, j++)
    {
        m_afKnot[i] = afKnot[j];
    }

    for (/**/; i < iNumKnots; i++)
    {
        m_afKnot[i] = (Real)1.0;
    }
}
//----------------------------------------------------------------------------
template <class Real>
BSplineBasis<Real>::~BSplineBasis ()
{
    WM3_DELETE[] m_afKnot;
    Deallocate(m_aafBD0);
    Deallocate(m_aafBD1);
    Deallocate(m_aafBD2);
    Deallocate(m_aafBD3);
}
//----------------------------------------------------------------------------
template <class Real>
int BSplineBasis<Real>::GetNumCtrlPoints () const
{
    return m_iNumCtrlPoints;
}
//----------------------------------------------------------------------------
template <class Real>
int BSplineBasis<Real>::GetDegree () const
{
    return m_iDegree;
}
//----------------------------------------------------------------------------
template <class Real>
bool BSplineBasis<Real>::IsOpen () const
{
    return m_bOpen;
}
//----------------------------------------------------------------------------
template <class Real>
bool BSplineBasis<Real>::IsUniform () const
{
    return m_bUniform;
}
//----------------------------------------------------------------------------
template <class Real>
Real BSplineBasis<Real>::GetD0 (int i) const
{
    return m_aafBD0[m_iDegree][i];
}
//----------------------------------------------------------------------------
template <class Real>
Real BSplineBasis<Real>::GetD1 (int i) const
{
    return m_aafBD1[m_iDegree][i];
}
//----------------------------------------------------------------------------
template <class Real>
Real BSplineBasis<Real>::GetD2 (int i) const
{
    return m_aafBD2[m_iDegree][i];
}
//----------------------------------------------------------------------------
template <class Real>
Real BSplineBasis<Real>::GetD3 (int i) const
{
    return m_aafBD3[m_iDegree][i];
}
//----------------------------------------------------------------------------
template <class Real>
Real** BSplineBasis<Real>::Allocate () const
{
    int iRows = m_iDegree + 1;
    int iCols = m_iNumCtrlPoints + m_iDegree;
    Real** aafArray;
    Wm3::Allocate<Real>(iCols,iRows,aafArray);
    memset(aafArray[0],0,iRows*iCols*sizeof(Real));
    return aafArray;
}
//----------------------------------------------------------------------------
template <class Real>
void BSplineBasis<Real>::Deallocate (Real** aafArray)
{
    if (aafArray)
    {
        Wm3::Deallocate<Real>(aafArray);
    }
}
//----------------------------------------------------------------------------
template <class Real>
int BSplineBasis<Real>::Initialize (int iNumCtrlPoints, int iDegree,
    bool bOpen)
{
    assert(iNumCtrlPoints >= 2);
    assert(1 <= iDegree && iDegree <= iNumCtrlPoints-1);

    m_iNumCtrlPoints = iNumCtrlPoints;
    m_iDegree = iDegree;
    m_bOpen = bOpen;

    int iNumKnots = m_iNumCtrlPoints+m_iDegree+1;
    m_afKnot = WM3_NEW Real[iNumKnots];

    m_aafBD0 = Allocate();
    m_aafBD1 = 0;
    m_aafBD2 = 0;
    m_aafBD3 = 0;

    return iNumKnots;
}
//----------------------------------------------------------------------------
template <class Real>
void BSplineBasis<Real>::SetKnot (int i, Real fKnot)
{
    if (!m_bUniform)
    {
        // access only allowed to elements d+1 <= j <= n
        int j = i + m_iDegree + 1;
        if (m_iDegree+1 <= j && j <= m_iNumCtrlPoints - 1)
        {
            m_afKnot[j] = fKnot;
        }
    }
}
//----------------------------------------------------------------------------
template <class Real>
Real BSplineBasis<Real>::GetKnot (int i) const
{
    if (!m_bUniform)
    {
        // access only allowed to elements d+1 <= j <= n
        int j = i + m_iDegree + 1;
        if (m_iDegree+1 <= j && j <= m_iNumCtrlPoints - 1)
        {
            return m_afKnot[j];
        }
    }

    return Math<Real>::MAX_REAL;
}
//----------------------------------------------------------------------------
template <class Real>
int BSplineBasis<Real>::GetKey (Real& rfTime) const
{
    if (m_bOpen)
    {
        // open splines clamp to [0,1]
        if (rfTime <= (Real)0.0)
        {
            rfTime = (Real)0.0;
            return m_iDegree;
        }
        else if (rfTime >= (Real)1.0)
        {
            rfTime = (Real)1.0;
            return m_iNumCtrlPoints-1;
        }
    }
    else
    {
        // periodic splines wrap to [0,1]
        if (rfTime < (Real)0.0 || rfTime > (Real)1.0)
        {
            rfTime -= Math<Real>::Floor(rfTime);
        }
    }


    int i;

    if (m_bUniform)
    {
        i = m_iDegree + (int)((m_iNumCtrlPoints-m_iDegree)*rfTime);
    }
    else
    {
        for (i = m_iDegree+1; i <= m_iNumCtrlPoints; i++)
        {
            if (rfTime < m_afKnot[i])
            {
                break;
            }
        }
        i--;
    }

    return i;
}
//----------------------------------------------------------------------------
template <class Real>
void BSplineBasis<Real>::Compute (Real fTime, unsigned int uiOrder,
    int& riMinIndex, int& riMaxIndex) const
{
    // only derivatives through third order currently supported
    assert(uiOrder <= 3);

    if (uiOrder >= 1)
    {
        if (!m_aafBD1)
        {
            m_aafBD1 = Allocate();
        }

        if (uiOrder >= 2)
        {
            if (!m_aafBD2)
            {
                m_aafBD2 = Allocate();
            }

            if (uiOrder >= 3)
            {
                if (!m_aafBD3)
                {
                    m_aafBD3 = Allocate();
                }
            }
        }
    }

    int i = GetKey(fTime);
    m_aafBD0[0][i] = (Real)1.0;

    if (uiOrder >= 1)
    {
        m_aafBD1[0][i] = (Real)0.0;
        if (uiOrder >= 2)
        {
            m_aafBD2[0][i] = (Real)0.0;
            if (uiOrder >= 3)
            {
                m_aafBD3[0][i] = (Real)0.0;
            }
        }
    }

    Real fN0 = fTime-m_afKnot[i], fN1 = m_afKnot[i+1]-fTime;
    Real fInvD0, fInvD1;
    int j;
    for (j = 1; j <= m_iDegree; j++)
    {
        fInvD0 = ((Real)1.0)/(m_afKnot[i+j]-m_afKnot[i]);
        fInvD1 = ((Real)1.0)/(m_afKnot[i+1]-m_afKnot[i-j+1]);

        m_aafBD0[j][i] = fN0*m_aafBD0[j-1][i]*fInvD0;
        m_aafBD0[j][i-j] = fN1*m_aafBD0[j-1][i-j+1]*fInvD1;

        if (uiOrder >= 1)
        {
            m_aafBD1[j][i] = (fN0*m_aafBD1[j-1][i]+m_aafBD0[j-1][i])*fInvD0;
            m_aafBD1[j][i-j] = (fN1*m_aafBD1[j-1][i-j+1]-m_aafBD0[j-1][i-j+1])
                *fInvD1;

            if (uiOrder >= 2)
            {
                m_aafBD2[j][i] = (fN0*m_aafBD2[j-1][i] +
                    ((Real)2.0)*m_aafBD1[j-1][i])*fInvD0;
                m_aafBD2[j][i-j] = (fN1*m_aafBD2[j-1][i-j+1] -
                    ((Real)2.0)*m_aafBD1[j-1][i-j+1])*fInvD1;

                if (uiOrder >= 3)
                {
                    m_aafBD3[j][i] = (fN0*m_aafBD3[j-1][i] +
                        ((Real)3.0)*m_aafBD2[j-1][i])*fInvD0;
                    m_aafBD3[j][i-j] = (fN1*m_aafBD3[j-1][i-j+1] -
                        ((Real)3.0)*m_aafBD2[j-1][i-j+1])*fInvD1;
                }
            }
        }
    }

    for (j = 2; j <= m_iDegree; j++)
    {
        for (int k = i-j+1; k < i; k++)
        {
            fN0 = fTime-m_afKnot[k];
            fN1 = m_afKnot[k+j+1]-fTime;
            fInvD0 = ((Real)1.0)/(m_afKnot[k+j]-m_afKnot[k]);
            fInvD1 = ((Real)1.0)/(m_afKnot[k+j+1]-m_afKnot[k+1]);

            m_aafBD0[j][k] = fN0*m_aafBD0[j-1][k]*fInvD0 + fN1*
                m_aafBD0[j-1][k+1]*fInvD1;

            if (uiOrder >= 1)
            {
                m_aafBD1[j][k] = (fN0*m_aafBD1[j-1][k]+m_aafBD0[j-1][k])*
                    fInvD0 + (fN1*m_aafBD1[j-1][k+1]-m_aafBD0[j-1][k+1])*
                    fInvD1;

                if (uiOrder >= 2)
                {
                    m_aafBD2[j][k] = (fN0*m_aafBD2[j-1][k] +
                        ((Real)2.0)*m_aafBD1[j-1][k])*fInvD0 +
                        (fN1*m_aafBD2[j-1][k+1]- ((Real)2.0)*
                        m_aafBD1[j-1][k+1])*fInvD1;

                    if (uiOrder >= 3)
                    {
                        m_aafBD3[j][k] = (fN0*m_aafBD3[j-1][k] +
                            ((Real)3.0)*m_aafBD2[j-1][k])*fInvD0 +
                            (fN1*m_aafBD3[j-1][k+1] - ((Real)3.0)*
                            m_aafBD2[j-1][k+1])*fInvD1;
                    }
                }
            }
        }
    }

    riMinIndex = i - m_iDegree;
    riMaxIndex = i;
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// explicit instantiation
//----------------------------------------------------------------------------
template WM3_ITEM
class BSplineBasis<float>;

template WM3_ITEM
class BSplineBasis<double>;
//----------------------------------------------------------------------------
}
