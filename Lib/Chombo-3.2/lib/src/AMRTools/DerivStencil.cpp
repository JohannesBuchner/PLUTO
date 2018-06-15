#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DerivStencil.H"
#include "NamespaceHeader.H"

DerivStencil::DerivStencil()
{
    clear();
}
void
DerivStencil::define()
{
    clear();
}
DerivStencil::DerivStencil(const DerivStencil& a_dsin)
{
    clear();
    *this = a_dsin;
}

///make derivstencil empty
void
DerivStencil::clear()
{
    m_vectIV.resize(0);
    m_vectWgt.resize(0);
    isdefined = true;
}
DerivStencil::~DerivStencil()
{
    clear();
}

///return length of vectors
int
DerivStencil::size() const
{
    CH_assert(isdefined);
    CH_assert(m_vectWgt.size() == m_vectIV.size());
    return(m_vectIV.size());
}

///get iv at ivec
const IntVect&
DerivStencil::getIndex(int a_ivec) const
{
    CH_assert(isdefined);
    CH_assert(a_ivec >= 0);
    CH_assert(a_ivec < m_vectIV.size());
    CH_assert(m_vectWgt.size() == m_vectIV.size());
    return(m_vectIV[a_ivec]);
}

///get weight at ivec
const Real&
DerivStencil::getWeight(int a_ivec) const
{
    CH_assert(isdefined);
    CH_assert(a_ivec >= 0);
    CH_assert(a_ivec < m_vectIV.size());
    CH_assert(m_vectIV.size()==m_vectWgt.size());
    return(m_vectWgt[a_ivec]);
}

/**
    add another set if the BAI is not in the
    stencil already.  \\
    **Add the weight to the existing weight otherwise**
    */
void
DerivStencil::accumulate(
    const IntVect& a_iv,
    Real a_weight
    )
{
    CH_assert(isdefined);
    CH_assert(m_vectIV.size()==m_vectWgt.size());
    bool isalreadythere = false;
    int ithere= -1;
    for (int ivect = 0; ivect < m_vectIV.size(); ivect++)
    {
        if (a_iv == m_vectIV[ivect])
        {
            isalreadythere = true;
            ithere = ivect;
            break;
        }
    }
    if (!isalreadythere)
    {
        m_vectWgt.push_back(a_weight);
        m_vectIV.push_back(a_iv);
    }
    else
    {
        m_vectWgt[ithere] += a_weight;
    }
}
const DerivStencil&
DerivStencil::operator=(const DerivStencil& a_dsin)
{
    CH_assert(isdefined);
    if (&a_dsin != this)
    {
        m_vectIV = a_dsin.m_vectIV;
        m_vectWgt = a_dsin.m_vectWgt;
    }
    return *this;
}

const DerivStencil&
DerivStencil::operator+=(Real a_facin)
{
    CH_assert(isdefined);
    CH_assert(m_vectIV.size()==m_vectWgt.size());
    for (int ivec = 0; ivec < m_vectWgt.size(); ivec++)
    {
        m_vectWgt[ivec] += a_facin;
    }
    return *this;
}

const DerivStencil&
DerivStencil::operator-=(Real a_facin)
{
    CH_assert(isdefined);
    CH_assert(m_vectIV.size()==m_vectWgt.size());
    for (int ivec = 0; ivec < m_vectWgt.size(); ivec++)
    {
        m_vectWgt[ivec] -= a_facin;
    }
    return *this;
}

const DerivStencil&
DerivStencil::operator*=(Real a_facin)
{
    CH_assert(isdefined);
    CH_assert(m_vectIV.size()==m_vectWgt.size());
    for (int ivec = 0; ivec < m_vectWgt.size(); ivec++)
    {
        m_vectWgt[ivec] *= a_facin;
    }
    return *this;
}

const DerivStencil&
DerivStencil::operator/=(Real a_facin)
{
    CH_assert(isdefined);
    CH_assert(m_vectIV.size()==m_vectWgt.size());
    for (int ivec = 0; ivec < m_vectWgt.size(); ivec++)
    {
        m_vectWgt[ivec] /= a_facin;
    }
    return *this;
}
#include "NamespaceFooter.H"
