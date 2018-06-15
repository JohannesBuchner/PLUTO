#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Correct1D2D.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "NamespaceHeader.H"


class EBAddOpCorr : public LDOperator<EBCellFAB>
{
public:
  EBAddOpCorr()
  {
  }

  virtual void linearIn(EBCellFAB& arg,  void* buf, const Box& R,
                        const Interval& comps) const
  {
    EBCellFAB incr;
    incr.clone(arg);
    incr.setVal(0.);
    incr.linearIn(buf, R, comps);

    EBCellFAB&       dst = arg;
    const EBCellFAB& src = incr;

    int isrc = comps.begin();
    int idst = comps.begin();
    int inco = comps.size();
    dst.plus(src, R, isrc, idst, inco);
  }

  void op(EBCellFAB& dst,
          const Box& RegionFrom,
          const Interval& Cdst,
          const Box& RegionTo,
          const EBCellFAB& src,
          const Interval& Csrc) const
  {
    CH_assert(Csrc.size() == Cdst.size());
    int isrc = Csrc.begin();
    int idst = Cdst.begin();
    int inco = Csrc.size();
    dst.plus(src, RegionTo, isrc, idst, inco);
  }
};
/****/
void
Correct1D2D::
makeIntMap(LevelData< BaseFab<int> >& a_intmap,
           const LayoutData<bool>   & a_is1D,
           const DisjointBoxLayout  & a_dbl)
{
  a_intmap.define(a_dbl, 1 , IntVect::Unit);
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      //set intmap to zero including ghost cells
      a_intmap[dit()].setVal(0);
      //set values ONLY in interior cells according to true false
      const Box& intbox = a_dbl[dit()];
      int comp = 0;
      if (a_is1D[dit()])
        {
          a_intmap[dit()].setVal(1, intbox, comp);
        }
      else
        {
          a_intmap[dit()].setVal(2, intbox, comp);
        }
    }
  //do exchange so box boundaries will have the right values
  //so now we can do the calculus to know  where stuff needs to be changed.
  a_intmap.exchange();
}
/****/
Correct1D2D::
Correct1D2D(const EBLevelGrid&      a_eblg,
            const LayoutData<bool>& a_is1D,
            int                     a_nvar)
{
  //make a map of where stuff is 1d vs 2d vs coarse fine interface.
  //this will end up being
  // 1 where a_is1d = true
  // 2 where a_is1d = false
  // 0 where undefined (coarse-fine interface or domain boundary)
  m_eblg = a_eblg;
  m_nvar = a_nvar;
  LevelData< BaseFab<int> > intmap;
  makeIntMap(intmap, a_is1D, m_eblg.getDBL());

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          //the 1d is in the interior of the source box
          //the 2d is in the exterior of the source box
          int iloc = index(idir, sit());
          if (sit() == Side::Lo)
            {
              adjCellLo(m_dbl1d[iloc], m_eblg.getDBL(), idir, -1);
              adjCellLo(m_dbl2d[iloc], m_eblg.getDBL(), idir,  1);
            }
          else
            {
              adjCellHi(m_dbl1d[iloc], m_eblg.getDBL(), idir, -1);
              adjCellHi(m_dbl2d[iloc], m_eblg.getDBL(), idir,  1);
            }

          //now define the buffers
          EBCellFactory fact(m_eblg.getEBISL());
          m_deltaU1d[iloc].define(m_dbl1d[iloc], m_nvar, IntVect::Zero, fact);
          m_deltaU2d[iloc].define(m_dbl2d[iloc], m_nvar, IntVect::Zero, fact);

          //now get the actual sets that will be operated upon
          m_sets1d[iloc].define(m_dbl1d[iloc]);
          m_sets2d[iloc].define(m_dbl2d[iloc]);

          for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
            {
              IntVectSet ivs1d(m_dbl1d[iloc][dit]);
              IntVectSet ivs2d(m_dbl2d[iloc][dit]);
              ivs1d &= m_eblg.getDomain().domainBox();
              ivs2d &= m_eblg.getDomain().domainBox();
              for (VoFIterator vofit(ivs1d, m_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
                {
                  const IntVect& thisIV = vofit().gridIndex();
                  IntVect thatIV = thisIV + sign(sit())*BASISV(idir);
                  bool found = ((intmap[dit()](thisIV, 0) == 1) && (intmap[dit()](thatIV, 0) == 2));
                  if (found)
                    {
                      m_sets1d[iloc][dit()].push_back(vofit());
                    }
                }
              for (VoFIterator vofit(ivs2d, m_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
                {
                  const IntVect& thisIV = vofit().gridIndex();
                  // negative because we are looking BACK at the box
                  // we started from
                  IntVect thatIV = thisIV - sign(sit())*BASISV(idir);

                  bool found = ((intmap[dit()](thisIV, 0) == 1) && (intmap[dit()](thatIV, 0) == 2));
                  if (found)
                    {
                      m_sets2d[iloc][dit()].push_back(vofit());
                    }
                }
            }
        }
    }
  setToZero();
}

///
/**
   sets buffers to zero
*/
void
Correct1D2D::
setToZero()
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int iloc = index(idir, sit());
          for (DataIterator dit = m_dbl1d[iloc].dataIterator(); dit.ok(); ++dit)
            {
              m_deltaU1d[iloc][dit()].setVal(0.);
            }
          for (DataIterator dit = m_dbl2d[iloc].dataIterator(); dit.ok(); ++dit)
            {
              m_deltaU2d[iloc][dit()].setVal(0.);
            }
        }
    }
}

///
/**
   increments the 1D (losing) buffer by -flux*scale*sign(side)
   (side is which side of the changed cell we are talking about)
   typically scale = 1/dx[idir]
*/
void
Correct1D2D::
increment1D(const EBFaceFAB& a_1DFlux,
            const Real&      a_scale,
            const DataIndex& a_dit)
{
  const EBISBox& ebisBox= m_eblg.getEBISL()[a_dit];
  for (SideIterator sit; sit.ok(); ++sit)
    {
      int idir = a_1DFlux.direction();
      int iloc = index(idir, sit());
      const Vector<VolIndex>& vofs = m_sets1d[iloc][a_dit];
      for (int ivof = 0; ivof < vofs.size(); ivof++)
        {
          Vector<FaceIndex> faces = ebisBox.getFaces(vofs[ivof], idir, sit());
          for (int iface = 0; iface < faces.size(); iface++)
            {
              Real area = ebisBox.areaFrac(faces[iface]);
              for (int ivar = 0; ivar < m_nvar; ivar++)
                {
                  m_deltaU1d[iloc][a_dit](vofs[ivof], ivar) -= a_1DFlux(faces[iface], ivar)*area*a_scale;
                }
            }
        }
    }
}


///
/**
   increments the 2D (winning) buffer by flux*scale*sign(side)
   (side is which side of the changed cell we are talking about)
   typically scale = 1/dx[idir]
*/
void
Correct1D2D::
increment2D(const EBFaceFAB& a_2DFlux,
                   const Real&      a_scale,
                   const DataIndex& a_dit)
{
  const EBISBox& ebisBox= m_eblg.getEBISL()[a_dit];
  for (SideIterator sit; sit.ok(); ++sit)
    {
      int idir = a_2DFlux.direction();
      int iloc = index(idir, sit());
      const Vector<VolIndex>& vofs = m_sets2d[iloc][a_dit];
      for (int ivof = 0; ivof < vofs.size(); ivof++)
        {
          //flip because of the way the set was constructed put the
          //box OUTside of the source box
          Vector<FaceIndex> faces = ebisBox.getFaces(vofs[ivof], idir, flip(sit()));
          for (int iface = 0; iface < faces.size(); iface++)
            {
              Real area = ebisBox.areaFrac(faces[iface]);
              for (int ivar = 0; ivar < m_nvar; ivar++)
                {
                  m_deltaU2d[iloc][a_dit](vofs[ivof], ivar) += a_2DFlux(faces[iface], ivar)*area*a_scale;
                }
            }
        }
    }
}


///
/**
   subtracts off change in solution due to losing flux
   and adds in change in solution due to winning flux.
*/
void
Correct1D2D::
correctSolution(LevelData<EBCellFAB>& a_U)
{
  int isrc = 0;  int idst = 0; int inc = m_nvar;
  //coarse stuff is already on the same layout as a_U
  //because it was shifted into its own box
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int iloc = index(idir, sit());
          for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
            {
              a_U[dit()].plus(m_deltaU1d[iloc][dit()], m_dbl1d[iloc][dit()], isrc, idst, inc);
            }
          //the two d change might be coming from different proc---hence the copyto
          Copier copier(m_dbl2d[iloc], a_U.disjointBoxLayout());
          EBAddOpCorr op;
          m_deltaU2d[iloc].copyTo(a_U, copier, op);
        }
    }
}

#include "NamespaceFooter.H"
