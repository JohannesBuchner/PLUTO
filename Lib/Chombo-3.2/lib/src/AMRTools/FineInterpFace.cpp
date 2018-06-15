#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "REAL.H"
#include "DataIterator.H"
#include "Tuple.H"
#include "InterpFace_F.H"
#include "InterpF_F.H"
#include "AverageFaceF_F.H"

#include "FineInterpFace.H"

#include "NamespaceHeader.H"

FineInterpFace::FineInterpFace()
  :
  is_defined(false)
{
}


FineInterpFace::~FineInterpFace()
{
}


FineInterpFace::FineInterpFace(const DisjointBoxLayout& a_fine_domain,
                               const int&  a_numcomps,
                               const int& a_ref_ratio,
                               const Box& a_coarse_problem_domain)
  :
  is_defined(false)
{
  ProblemDomain physdomain(a_coarse_problem_domain);
  define(a_fine_domain, a_numcomps, a_ref_ratio, physdomain);
}


FineInterpFace::FineInterpFace(const DisjointBoxLayout& a_fine_domain,
                               const int&  a_numcomps,
                               const int& a_ref_ratio,
                               const ProblemDomain& a_fine_problem_domain)
  :
  is_defined(false)
{
  define(a_fine_domain, a_numcomps, a_ref_ratio, a_fine_problem_domain);
}


void
FineInterpFace::define(const DisjointBoxLayout& a_fine_domain,
                       const int& a_numcomps,
                       const int& a_ref_ratio,
                       const Box& a_problem_domain)
{
  ProblemDomain physdomain(a_problem_domain);
  define(a_fine_domain, a_numcomps, a_ref_ratio, physdomain);
}


void
FineInterpFace::define(const DisjointBoxLayout& a_fine_domain,
                       const int& a_numcomps,
                       const int& a_ref_ratio,
                       const ProblemDomain& a_problem_domain)
{
  CH_assert (a_fine_domain.checkPeriodic(a_problem_domain));
  m_ref_ratio = a_ref_ratio;
  m_coarse_problem_domain = coarsen(a_problem_domain,m_ref_ratio);
//
// create the work array
  DisjointBoxLayout coarsened_fine_domain;
  coarsen ( coarsened_fine_domain,
            a_fine_domain,
            m_ref_ratio );
  m_coarsened_fine_data.define ( coarsened_fine_domain,
                                 a_numcomps,
                                 IntVect::Unit );
  is_defined = true;
}



bool
FineInterpFace::isDefined() const
{
  return ( is_defined );
}


// interpolate from coarse level to fine level
void
FineInterpFace::interpToFine(LevelData<FluxBox>& a_fine_data,
                             const LevelData<FluxBox>& a_coarse_data,
                             bool a_averageFromDest)
{
  CH_assert(is_defined);

  if (a_averageFromDest)
    {
      // average down fine data -- this is a local operation
      DataIterator dit = a_fine_data.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          for (int dir=0; dir<SpaceDim; dir++)
            {
              FArrayBox& fineFab = a_fine_data[dit][dir];
              FArrayBox& crseFab = m_coarsened_fine_data[dit][dir];
              const Box& crseBox = crseFab.box();

              // set up refinement box
              int boxHi = m_ref_ratio-1;
              IntVect hiVect(D_DECL6(boxHi,boxHi,boxHi,
                                     boxHi,boxHi,boxHi));
              // don't want to index at all in dir direction --
              // instead, want to just march along face.
              hiVect.setVal(dir,0);
              IntVect loVect(D_DECL6(0,0,0,0,0,0));
              Box refBox(loVect, hiVect);              
              int refFactor = m_ref_ratio;

              FORT_AVERAGEFACE( CHF_FRA(crseFab),
                                CHF_CONST_FRA(fineFab),
                                CHF_BOX(crseBox),
                                CHF_CONST_INT(dir),
                                CHF_CONST_INT(m_ref_ratio),
                                CHF_CONST_INT(refFactor),
                                CHF_BOX(refBox));
            } // end loop over face directions
        } // end loop over grid boxes
    } // endif average from dest

  a_coarse_data.copyTo(a_coarse_data.interval(),
                       m_coarsened_fine_data,
                       m_coarsened_fine_data.interval() );

  const BoxLayout fine_domain = a_fine_data.boxLayout();
  DataIterator dit = fine_domain.dataIterator();

  //int nbox = 0;
  for (dit.begin(); dit.ok(); ++dit)
  {
    // now loop over directions
    for (int dir=0; dir<SpaceDim; dir++)
    {
      const BaseFab<Real>& coarsened_fine = m_coarsened_fine_data[dit()][dir];
      //      const Box& coarsened_fine_box =
      // m_coarsened_fine_data[dit()][dir].box();
      Box coarsened_fine_box(a_fine_data.getBoxes()[dit()]);
      coarsened_fine_box.coarsen(m_ref_ratio);
      coarsened_fine_box.surroundingNodes(dir);
      FluxBox& fineFluxbox = a_fine_data[dit()];
      BaseFab<Real>& fine = fineFluxbox[dir];
      // interpGridData interpolates from an entire coarse grid onto an
      // entire fine grid.
      interpGridData(fine,
                     coarsened_fine,
                     coarsened_fine_box,
                     m_ref_ratio, dir);
    } // end loop over directions
    //pout () << "end of box loop" << nbox << endl;
    //++nbox;
  }
}

// interpolate from fine grid to coarse grid.  prerequisite:
// coarsened.box contains coarsen(fine.box).
//
// uses piecewise bilinear interpolation with multidimensional-limited
// slopes.  see design document for details.
void
FineInterpFace::interpGridData(BaseFab<Real>& a_fine,
                               const BaseFab<Real>& a_coarse,
                               const Box& a_coarsened_fine_box,
                               const int& a_ref_ratio,
                               const int& a_faceDir)
  const
{
  // fill fine data with piecewise constant coarse data
  const Box& b = a_coarsened_fine_box;
  const int num_comp = a_fine.nComp ();

  // define refbox for edges
  Box facerefbox(IntVect::Zero,
             (a_ref_ratio-1)*IntVect::Unit);
  facerefbox.surroundingNodes(a_faceDir);
  facerefbox.setBig(a_faceDir,0);

  FORT_INTERPFACECONSTANT ( CHF_FRA(a_fine),
                            CHF_CONST_FRA(a_coarse),
                            CHF_BOX(b),
                            CHF_CONST_INT(a_ref_ratio),
                            CHF_BOX(facerefbox),
                            CHF_CONST_INT(a_faceDir)
                            );
  //  Tuple<BaseFab<Real>, SpaceDim> slopes;
  //  for (int dir = 0; dir < SpaceDim; ++dir)
  // hardwired to 6 due to lack of variable number of arguments in chfpp
  // this is designed to accomodate the Chombo max-dimension of 6
#define MAXDIM 6
  BaseFab<Real> slopes[MAXDIM];


  // this is a trick to make domain face-centered while still
  // respecting periodic BC's; take domain, then grow
  // by one (to allow for centered differences), then intersect
  // with domain.  this will result in a box that is larger than
  // the computational domain in periodic directions, but is
  // the same as the periodic domain in non-periodic directions.
  // so when we later intersect boxes with this domain, it will
  // cut off cells in non-periodic directions while leaving
  // periodic directions intact.
  Box domainFaceBox(m_coarse_problem_domain.domainBox());
  // using a grow radius of two here, but probably could have used 1
  domainFaceBox.grow(2);
  domainFaceBox &= m_coarse_problem_domain;
  domainFaceBox.surroundingNodes(a_faceDir);

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      BaseFab<Real>& dir_slope = slopes[dir];
      dir_slope.resize(b, num_comp);
    }
  // define the extras over a unit box in order to avoid issues with null
  // pointers and undefined FAB's
  for (int dir=SpaceDim; dir<MAXDIM; dir++)
    {
      Box unitBox(IntVect::Zero, IntVect::Zero);
      BaseFab<Real>& dir_slope = slopes[dir];
      dir_slope.resize(unitBox, num_comp);
    }

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      // only do slopes if not doing normal direction
      if (dir != a_faceDir)
      {
        BaseFab<Real>& dir_slope = slopes[dir];

        const Box bcenter = b & grow(domainFaceBox,-BASISV(dir));
        if (!bcenter.isEmpty())
          {
            FORT_INTERPCENTRALSLOPE ( CHF_FRA ( dir_slope ),
                                      CHF_CONST_FRA ( a_coarse ),
                                      CHF_BOX ( bcenter ),
                                      CHF_CONST_INT ( dir )
                                      );
          }
        const Box blo = b & adjCellLo(grow(domainFaceBox,-BASISV(dir)),dir);
        if (!blo.isEmpty())
          {
            FORT_INTERPHISIDESLOPE ( CHF_FRA ( dir_slope ),
                                     CHF_CONST_FRA ( a_coarse ),
                                     CHF_BOX ( blo ),
                                     CHF_CONST_INT ( dir )
                                     );
          }
        const Box bhi = b & adjCellHi(grow(domainFaceBox,
                                           -BASISV(dir)),dir);
        if (!bhi.isEmpty())
          {
            FORT_INTERPLOSIDESLOPE ( CHF_FRA ( dir_slope ),
                                     CHF_CONST_FRA ( a_coarse ),
                                     CHF_BOX ( bhi ),
                                     CHF_CONST_INT ( dir )
                                     );
          }
      }
      else
      {
        // if this is the normal direction, just fill with 0's
        slopes[dir].setVal(0.0);
      }
    } // end loop over directions

  // need a box which covers all neighbors in the plane of the
  // coarse face. this should result in a box that's (-1,1)
  // in the plane of the face, and (0,0) normal to the face)
  Box neighborBox(-1*IntVect::Unit,
                  IntVect::Unit);
  neighborBox.grow(a_faceDir,-1);

  // this will be a box which defines which cells in a_coarse
  // are valid (or, at least, which ones exist!)

  // first construct box of valid edges in physical domain
  Box validBox(m_coarse_problem_domain.domainBox());
  validBox.grow(2);
  validBox &= m_coarse_problem_domain;
  validBox.surroundingNodes(a_faceDir);
  // now intersect with existing cells in a_coarse
  validBox &= a_coarse.box();


  FORT_INTERPLIMITFACE ( CHF_FRA ( slopes[0] ),
                         CHF_FRA ( slopes[1] ),
                         CHF_FRA ( slopes[2] ),
                         CHF_FRA ( slopes[3] ),
                         CHF_FRA ( slopes[4] ),
                         CHF_FRA ( slopes[5] ),
                         CHF_CONST_FRA ( a_coarse ),
                         CHF_BOX ( b ),
                         CHF_BOX ( neighborBox),
                         CHF_BOX (validBox),
                         CHF_CONST_INT ( a_faceDir)
                         );

  // do (bi)linear interpolation on fine faces which overlie coarse
  // faces

  // need a cell-centered coarse box
  Box coarseCellBox = b;
  coarseCellBox.enclosedCells(a_faceDir);

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (dir != a_faceDir)
        {
          BaseFab<Real>& dir_slope = slopes[dir];
          FORT_INTERPLINEARFACE ( CHF_FRA ( a_fine ),
                                    CHF_CONST_FRA ( dir_slope ),
                                    CHF_BOX ( b ),
                                    CHF_CONST_INT ( dir ),
                                    CHF_CONST_INT ( a_ref_ratio ),
                                    CHF_BOX ( facerefbox )
                                    );
        }
    }

  // finally, do interiors

  // this box will contain the interior faces which do
  // not overlie a coarse face
  Box interiorRefBox(IntVect::Zero,
                     (a_ref_ratio-1)*IntVect::Unit);
  interiorRefBox.surroundingNodes(a_faceDir);
  // remove outer faces from interior box...
  interiorRefBox.grow(a_faceDir, -1);

  FORT_INTERPLINEARINTERIORFACE(CHF_FRA(a_fine),
                                CHF_BOX(coarseCellBox),
                                CHF_CONST_INT(a_ref_ratio),
                                CHF_CONST_INT(a_faceDir),
                                CHF_BOX(interiorRefBox));

}


#include "NamespaceFooter.H"




