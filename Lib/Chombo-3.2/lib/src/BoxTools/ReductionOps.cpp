#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ReductionOps.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

SumOp::SumOp():scale(1.0)
{
}

SumOp::SumOp(int a_summingDir ):scale(1.0)
{
  m_summingDir.resize(1);
  m_summingDir[0] = a_summingDir;
}

SumOp::SumOp(const Vector<int>& a_summingDir ):scale(1.0)
{
  m_summingDir = a_summingDir;
}

void
SumOp::linearOut(const FArrayBox& arg, void* buf, const Box& R,
                 const Interval& comps) const
{
  Real* buffer = (Real*)buf;

  Box flattenedBox(R);
  for (int n=0; n<m_summingDir.size(); n++)
    {
      flattenedBox.setBig(m_summingDir[n],R.smallEnd(m_summingDir[n]));
    }

  // single summing direction
  if (m_summingDir.size() == 1)
    {
      int transverseLo = R.smallEnd(m_summingDir[0]);
      int transverseHi = R.bigEnd(m_summingDir[0]);
      // don't apply scale here, since we'll do it in the linearIn side of
      // things.
      ForAllXCBNN(Real, arg, flattenedBox, comps.begin(), comps.size())
        {
          // this gets rid of the "unused variable" warning
          // while we also initialize the buffer
          //*buffer = 0;
          *buffer = 0*argR;

          IntVect iv(D_DECL6(iR, jR, kR,
                             _iv[3], _iv[4], _iv[5]));
          for (int itrans=transverseLo; itrans<=transverseHi; itrans++)
            {
              iv[m_summingDir[0]] = itrans;
              *buffer += arg(iv,_n);
            }
          buffer++;
        } EndFor
    }
  else
    {
      // multiple summing directions -- need to do this a bit differently
      Box transverseBox(IntVect::Zero, IntVect::Zero);
      for (int n=0; n<m_summingDir.size(); n++)
        {
          transverseBox.setBig(m_summingDir[n],
                               (R.bigEnd(m_summingDir[n])-R.smallEnd(m_summingDir[n])));
          //transverseBox.setSmall(m_summingDir[n],R.smallEnd(m_summingDir[n]));
        }
      // don't apply scale here, since we'll do it in the linearIn side of
      // things.
      ForAllXCBNN(Real, arg, flattenedBox, comps.begin(), comps.size())
        {
          // this gets rid of the "unused variable" warning
          // while we also initialize the buffer
          //*buffer = 0;
          *buffer = 0*argR;

          IntVect iv(D_DECL6(iR, jR, kR,
                             _iv[3], _iv[4], _iv[5]));
          BoxIterator transverseIt(transverseBox);
          for (transverseIt.begin(); transverseIt.ok(); ++transverseIt)
                {
                  IntVect shift = transverseIt();
                  IntVect srcLoc = iv + shift;
                  *buffer += arg(srcLoc,_n);
                } // end loop over transverse directions
              buffer++;
        } EndFor
    }
}

void
SumOp::linearIn(FArrayBox& arg,  void* buf, const Box& R,
                const Interval& comps) const
{
  Real* buffer = (Real*)buf;
  if (scale != 1.0)
    {
      ForAllXBNNnoindx(Real, arg, R, comps.begin(), comps.size())
        {
          argR+=*buffer * scale;
          buffer++;
        } EndFor
            }
  else
    {
      ForAllXBNNnoindx(Real, arg, R, comps.begin(), comps.size())
        {
          argR+=*buffer;
          buffer++;
        } EndFor
            }
}

void
SumOp::op(FArrayBox& dest,
          const Box& RegionFrom,
          const Interval& Cdest,
          const Box& RegionTo,
          const FArrayBox& src,
          const Interval& Csrc) const
{
//int numComp = Cdest.size();
  CH_assert(Cdest.size() == Csrc.size());
  CH_assert(dest.nComp() > Cdest.end());
  CH_assert(src.nComp() > Csrc.end());

  int offset = Cdest.begin() - Csrc.begin();

  // start by doing this with a boxIterator -- can make this more efficient
  // later if necessary
  if (m_summingDir.size() == 1)
    {
      int toBoxLo = RegionTo.smallEnd(m_summingDir[0]);
      int toBoxHi = RegionTo.bigEnd(m_summingDir[0]);

      BoxIterator fromBit(RegionFrom);
      for (fromBit.begin(); fromBit.ok(); ++fromBit)
        {
          IntVect fromIV = fromBit();
          //int summingIndex = fromIV[m_summingDir];

          if (toBoxLo == toBoxHi)
            {
              IntVect toIV = fromIV;
              toIV[m_summingDir[0]] = toBoxLo;
              for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                {
                  dest(toIV, comp+offset) += scale*src(fromIV, comp);
                } // end loop over components
            }
          else
            {
              // this is in case the toBox is more than one cell wide in
              // the summing direction.  This is most likely unnecessary,
              // to be honest.
              Box toBox(fromIV, fromIV);
              toBox.setSmall(m_summingDir[0], toBoxLo);
              toBox.setBig(m_summingDir[0], toBoxHi);

              BoxIterator toBit(toBox);
              for (toBit.begin(); toBit.ok(); ++toBit)
                {
                  IntVect toIV = toBit();

                  for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                    {
                      dest(toIV, comp+offset) += scale*src(fromIV, comp);
                    } // end loop over components
                } // end loop over toBox
            }
        } // end loop over from cells
    } // end if single summing direction
  else
    {
      int width = 1;
      for (int n=0; n<m_summingDir.size(); n++)
        {
          int toBoxLo = RegionTo.smallEnd(m_summingDir[n]);
          int toBoxHi = RegionTo.bigEnd(m_summingDir[n]);
          width *= (toBoxHi - toBoxLo+1);
        }

      BoxIterator fromBit(RegionFrom);
      for (fromBit.begin(); fromBit.ok(); ++fromBit)
        {
          IntVect fromIV = fromBit();
          //int summingIndex = fromIV[m_summingDir];

          if (width == 1)
            {
              IntVect toIV = fromIV;
              for (int n=0; n<m_summingDir.size(); n++)
                {
                  toIV[m_summingDir[n]] = RegionTo.smallEnd(m_summingDir[n]);
                }
              for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                {
                  dest(toIV, comp+offset) += scale*src(fromIV, comp);
                } // end loop over components
            }
          else
            {
              // this is in case the toBox is more than one cell wide in
              // the summing direction.  This is most likely unnecessary,
              // to be honest.
              Box toBox(fromIV, fromIV);
              for (int n=0; n<m_summingDir.size(); n++)
                {
                  toBox.setSmall(m_summingDir[n], RegionTo.smallEnd(m_summingDir[n]));
                  toBox.setBig(m_summingDir[n], RegionTo.bigEnd(m_summingDir[n]));
                }

              BoxIterator toBit(toBox);
              for (toBit.begin(); toBit.ok(); ++toBit)
                {
                  IntVect toIV = toBit();

                  for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                    {
                      dest(toIV, comp+offset) += scale*src(fromIV, comp);
                    } // end loop over components
                } // end loop over toBox
            } // end if toRegion is more than one cell "wide"
        } // end loop over from cells
    } // end if more than one summing direction
}

// (DFM 11/13/08) as currently implemented, FaceSumOp doesn't
// do the right thing for multiple adjoining grids -- it will
// double-count overlying faces where boxes adjoin. Since we don't
// actually need the face-centered summing operator at the moment,
// take the cowardly path of just commenting it out to revisit if it
// becomes a needed member of the Chombo family

#if 0

// --------------------------------------------
// face-centered summing operator
// --------------------------------------------

FaceSumOp::FaceSumOp():scale(1.0), m_summingDir(-1)
{
}

FaceSumOp::FaceSumOp(int a_summingDir ):scale(1.0)
{
  m_summingDir = a_summingDir;
}

void
FaceSumOp::linearIn(FluxBox& arg,  void* buf, const Box& R,
                    const Interval& comps) const
{
  Real* buffer = (Real*)buf;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& argDir = arg[dir];
      if (scale != 1.0)
        {
          ForAllXBNNnoindx(Real, argDir, R, comps.begin(), comps.size())
            {
              argDirR+=*buffer * scale;
              buffer++;
            } EndFor
                }
      else
        {
          ForAllXBNNnoindx(Real, argDir, R, comps.begin(), comps.size())
            {
              argDirR+=*buffer;
              buffer++;
        } EndFor
            }
    } // end loop over dir
}

void
FaceSumOp::op(FluxBox& dest,
              const Box& RegionFrom,
              const Interval& Cdest,
              const Box& RegionTo,
              const FluxBox& src,
              const Interval& Csrc) const
{
//int numComp = Cdest.size();
  CH_assert(Cdest.size() == Csrc.size());
  CH_assert(dest.nComp() > Cdest.end());
  CH_assert(src.nComp() > Csrc.end());

  int offset = Cdest.begin() - Csrc.begin();

  // start by doing this with a boxIterator -- can make this more efficient
  // later if necessary

  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& destDir = dest[dir];
      const FArrayBox& srcDir = src[dir];
      Box RegionToDir(RegionTo);
      RegionToDir.surroundingNodes(dir);
      Box RegionFromDir(RegionFrom);
      RegionFromDir.surroundingNodes(dir);

      int toBoxLo = RegionToDir.smallEnd(m_summingDir);
      int toBoxHi = RegionToDir.bigEnd(m_summingDir);
      // account for the fact that we get an extra face
      // in the dir direction -- for the purposes of this
      // op, sum to the low-side face only.
      if (dir == m_summingDir) toBoxHi-= 1;

      BoxIterator fromBit(RegionFromDir);
      for (fromBit.begin(); fromBit.ok(); ++fromBit)
        {
          IntVect fromIV = fromBit();
          //int summingIndex = fromIV[m_summingDir];

          if (toBoxLo == toBoxHi)
            {
              IntVect toIV = fromIV;
              toIV[m_summingDir] = toBoxLo;
              for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                {
                  destDir(toIV, comp+offset) += scale*srcDir(fromIV, comp);
                } // end loop over components
            }
          else
            {
              // this is in case the toBox is more than one cell wide in
              // the summing direction.  This is most likely unnecessary,
              // to be honest.
              Box toBox(fromIV, fromIV);
              toBox.setSmall(m_summingDir, toBoxLo);
              toBox.setBig(m_summingDir, toBoxHi);

              BoxIterator toBit(toBox);
              for (toBit.begin(); toBit.ok(); ++toBit)
                {
                  IntVect toIV = toBit();

                  for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                    {
                      destDir(toIV, comp+offset) += scale*srcDir(fromIV, comp);
                    } // end loop over components
                } // end loop over toBox
            }
        } // end loop over from cells
    } // end loop over face dirs
}

#endif // end of commenting out FaceSumOp

// --------------------------------------------
// spreading operator
// --------------------------------------------

SpreadingOp::SpreadingOp():scale(1.0)
{
}

SpreadingOp::SpreadingOp(int a_spreadingDir ):scale(1.0)
{
  m_spreadingDir.resize(1);
  m_spreadingDir[0] = a_spreadingDir;
}

SpreadingOp::SpreadingOp(const Vector<int>& a_spreadingDir ):scale(1.0)
{
  m_spreadingDir = a_spreadingDir;
}

void
SpreadingOp::linearIn(FArrayBox& arg,  void* buf, const Box& R,
                      const Interval& comps) const
{
  Real* buffer = (Real*)buf;
  Box flattenedBox(R);
  for (int n=0; n<m_spreadingDir.size(); n++)
    {
      flattenedBox.setBig(m_spreadingDir[n],R.smallEnd(m_spreadingDir[n]));
    }

  if (scale != 1.0)
    {
      // this copies from buf to the low-end of arg in the spreadingDir
      ForAllXBNNnoindx(Real, arg, flattenedBox, comps.begin(), comps.size())
        {
          argR=*buffer * scale;
          buffer++;
        } EndFor
            }
  else
    {
      ForAllXBNNnoindx(Real, arg, flattenedBox, comps.begin(), comps.size())
        {
          argR=*buffer;
          buffer++;
        } EndFor
            }
  // now need to copy from low-end of R to fill the rest of the box
  // do this by calling SpreadingOp::applyOp with scale = 1 (since
  // we already applied the scale above)
  this->applyOp(arg, flattenedBox, comps, R, arg, comps, 1.0);

}

void
SpreadingOp::op(FArrayBox& dest,
                const Box& RegionFrom,
                const Interval& Cdest,
                const Box& RegionTo,
                const FArrayBox& src,
                const Interval& Csrc) const
{
  applyOp(dest, RegionFrom, Cdest, RegionTo, src, Csrc, scale);
}

void
SpreadingOp::applyOp(FArrayBox& dest,
                     const Box& RegionFrom,
                     const Interval& Cdest,
                     const Box& RegionTo,
                     const FArrayBox& src,
                     const Interval& Csrc,
                     Real a_scale) const
{
//int numComp = Cdest.size();
  CH_assert(Cdest.size() == Csrc.size());
  CH_assert(dest.nComp() > Cdest.end());
  CH_assert(src.nComp() > Csrc.end());

  int offset = Cdest.begin() - Csrc.begin();

  // start by doing this with a boxIterator -- can make this more efficient
  // later if necessary
  if (m_spreadingDir.size() == 1)
    {
      int fromBoxLo = RegionFrom.smallEnd(m_spreadingDir[0]);
      int fromBoxHi = RegionFrom.bigEnd(m_spreadingDir[0]);

      // this doesn't make any sense if regionFrom is more than one wide
      CH_assert(fromBoxLo == fromBoxHi);

      BoxIterator toBit(RegionTo);
      for (toBit.begin(); toBit.ok(); ++toBit)
        {
          IntVect toIV = toBit();

          IntVect fromIV = toIV;
          fromIV[m_spreadingDir[0]] = fromBoxLo;
          for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
            {
              dest(toIV, comp+offset) = a_scale*src(fromIV, comp);
            } // end loop over components

        } // end loop over from cells
    } // end if only one spreading direction
  else
    {
      // multiple spreading directions
      // this doesn't make any sense if fromBox is more than
      // one cell wide in any of the transverse directions
      for (int n=0; n<m_spreadingDir.size(); n++)
        {
          CH_assert(RegionFrom.size(m_spreadingDir[n]) == 1);
        }

      BoxIterator toBit(RegionTo);
      for (toBit.begin(); toBit.ok(); ++toBit)
        {
          IntVect toIV = toBit();
          IntVect fromIV = toIV;
          for (int n=0; n<m_spreadingDir.size(); n++)
            {
              fromIV[m_spreadingDir[n]] =RegionFrom.smallEnd(m_spreadingDir[n]);
            }
          for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
            {
              dest(toIV, comp+offset) = a_scale*src(fromIV, comp);
            } // end loop over components
        } // end loop over from cells
    } // end if more than one spreading direction
}

// --------------------------------------------
// face-centered spreading operator
// --------------------------------------------

FaceSpreadingOp::FaceSpreadingOp():scale(1.0)
{
}

FaceSpreadingOp::FaceSpreadingOp(int a_spreadingDir ):scale(1.0)
{
  m_spreadingDir.resize(1);
  m_spreadingDir[0] = a_spreadingDir;
}

FaceSpreadingOp::FaceSpreadingOp(const Vector<int>& a_spreadingDir ):scale(1.0)
{
  m_spreadingDir = a_spreadingDir;
}

void
FaceSpreadingOp::linearIn(FluxBox& arg,  void* buf, const Box& R,
                          const Interval& comps) const
{
  Real* buffer = (Real*)buf;

  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& argDir = arg[dir];
      Box faceR(R);
      Box flattenedBox(faceR);
      for (int n=0; n<m_spreadingDir.size(); n++)
        {
          flattenedBox.setBig(m_spreadingDir[n],faceR.smallEnd(m_spreadingDir[n]));
        }

      flattenedBox.surroundingNodes(dir);
      if (scale != 1.0)
        {
          ForAllXBNNnoindx(Real, argDir, flattenedBox, comps.begin(), comps.size())
            {
              argDirR=*buffer * scale;
              buffer++;
            } EndFor
                }
      else
        {
          ForAllXBNNnoindx(Real, argDir, flattenedBox, comps.begin(), comps.size())
            {
              argDirR=*buffer;
              buffer++;
            } EndFor
                }
    }

  // now need to copy from low-end of R to fill the rest of the box
  // do this by calling SpreadingOp::applyOp with scale = 1 (since
  // we already applied the scale above)
  Box flattenedBox(R);
  for (int n=0; n<m_spreadingDir.size(); n++)
    {
      flattenedBox.setBig(m_spreadingDir[n],R.smallEnd(m_spreadingDir[n]));
    }
  this->applyOp(arg, flattenedBox, comps, R, arg, comps, 1.0);

}

void
FaceSpreadingOp::op(FluxBox& dest,
                    const Box& RegionFrom,
                    const Interval& Cdest,
                    const Box& RegionTo,
                    const FluxBox& src,
                    const Interval& Csrc) const
{
  applyOp(dest, RegionFrom, Cdest, RegionTo, src, Csrc, scale);
}

void
FaceSpreadingOp::applyOp(FluxBox& dest,
                         const Box& RegionFrom,
                         const Interval& Cdest,
                         const Box& RegionTo,
                         const FluxBox& src,
                         const Interval& Csrc,
                         Real a_scale) const
{
//int numComp = Cdest.size();
  CH_assert(Cdest.size() == Csrc.size());
  CH_assert(dest.nComp() > Cdest.end());
  CH_assert(src.nComp() > Csrc.end());

  int offset = Cdest.begin() - Csrc.begin();

  if (m_spreadingDir.size() == 1)
    {
      // start by doing this with a boxIterator -- can make this more efficient
      // later if necessary
      for (int dir=0; dir<SpaceDim; dir++)
        {

          FArrayBox& destDir = dest[dir];
          const FArrayBox& srcDir = src[dir];
          Box RegionFromDir(RegionFrom);
          RegionFromDir.surroundingNodes(dir);

          // since a 1-wide cell-centered box will wind up with 2 faces
          // in the normal direction, make a design decision to only spread
          // from the lower face in the spreadingDir direction
          if (dir == m_spreadingDir[0]) RegionFromDir.growHi(dir,-1);

          Box RegionToDir(RegionTo);
          RegionToDir.surroundingNodes(dir);

          int fromBoxLo = RegionFromDir.smallEnd(m_spreadingDir[0]);
          int fromBoxHi = RegionFromDir.bigEnd(m_spreadingDir[0]);

          // this doesn't make any sense if regionFrom is more than one wide
          CH_assert(fromBoxLo == fromBoxHi);

          BoxIterator toBit(RegionToDir);
          for (toBit.begin(); toBit.ok(); ++toBit)
            {
              IntVect toIV = toBit();

              IntVect fromIV = toIV;
              fromIV[m_spreadingDir[0]] = fromBoxLo;
              for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                {
                  destDir(toIV, comp+offset) = a_scale*srcDir(fromIV, comp);
                } // end loop over components

            } // end loop over from cells
        } // end loop over face dirs
    } // end if only one spreading direction
  else
    {
      // multiple directions
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& destDir = dest[dir];
          const FArrayBox& srcDir = src[dir];
          Box RegionToDir(RegionTo);
          RegionToDir.surroundingNodes(dir);

          // this doesn't make any sense if the FromRegion is more than
          // one cell wide in the spreading dir
          for (int n=0; n<m_spreadingDir.size(); n++)
            {
              CH_assert(RegionFrom.size(m_spreadingDir[n]) == 1);
            }

          BoxIterator toBit(RegionToDir);
          for (toBit.begin(); toBit.ok(); ++toBit)
            {
              IntVect toIV = toBit();
              IntVect fromIV = toIV;

              for (int n=0; n<m_spreadingDir.size(); n++)
                {
                  fromIV[m_spreadingDir[n]] = RegionFrom.smallEnd(m_spreadingDir[n]);
                }
              for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                {
                  destDir(toIV, comp+offset) = a_scale*srcDir(fromIV, comp);
                } // end loop over components
            } // end loop over "to" cells
        } // end loop over face directions
    } // end if more than one spreading direction
}

#include "NamespaceFooter.H"
