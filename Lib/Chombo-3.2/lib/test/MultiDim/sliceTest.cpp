#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_MPI
#include <mpi.h>
#endif

#include "Slicing.H.transdim"
#include "Injection.H.transdim"
#include "FArrayBox.H.multidim"
#include "FluxBox.H.multidim"
#include "DebugOut.H.multidim"
#include "BRMeshRefine.H.multidim"
#include "LoadBalance.H.multidim"
#include "FABView.H.multidim"

#include "ReductionCopier.H.multidim"
#include "SpreadingCopier.H.multidim"
#include "ReductionOps.H.multidim"

//#include "unidim.H"
#include <iostream>

using namespace Chombo;
using namespace CH_MultiDim;




template<typename T> std::ostream&
operator<<( std::ostream& a_out, const D2::LevelData<T>& a_ld)
{
  D2::DataIterator dit = a_ld.dataIterator();
  for ( ; dit.ok(); ++dit )
    {
      a_out << a_ld[dit()] << '\n';
    }

  return a_out;
}
template<typename T> std::ostream&
operator<<( std::ostream& a_out, const D3::LevelData<T>& a_ld)
{
  D3::DataIterator dit = a_ld.dataIterator();
  for ( ; dit.ok(); ++dit )
    {
      a_out << a_ld[dit()] << '\n';
    }

  return a_out;
}

/** FAB_T is LevelData's template argument.  This will only work if it's
 *  some kind of FArrayBox -- D1::FArrayBox, D2::FArrayBox, etc.
*/
template<typename FAB_T> struct LevelDataApplyHelper
  : public ObjTraits<DimTraits<FAB_T>::dim>::LevelData::ApplyFunctor
{
    LevelDataApplyHelper( int i ) : m_i(i)
    {
    }

    enum
    {
      DIM=DimTraits<FAB_T>::dim
    };
    virtual void operator()( const typename ObjTraits<DIM>::Box& box,
                             int nComps,
                             FAB_T& fab ) const
    {
        for ( int i=0;i<fab.box().numPts();++i )
        {
            fab.dataPtr()[i] = i;
        }
    }
  private:
    int m_i;
};


/** FAB_T is LevelData's template argument.  This will only work if it's
 *  some kind of FArrayBox -- D1::FArrayBox, D2::FArrayBox, etc.
*/
template<typename FAB_T> struct LevelFluxApplyHelper
  : public ObjTraits<DimTraits<FAB_T>::dim>::LevelFlux::ApplyFunctor
{
  LevelFluxApplyHelper( int i )
    :
    m_i(i)
  {
  }

  enum
  {
    DIM=DimTraits<FAB_T>::dim
  };

  virtual void operator()( const typename ObjTraits<DIM>::Box& box,
                           int nComps,
                           FAB_T& fab ) const
  {
    for (int dir=0; dir<DIM; dir++)
      {
        for ( int i=0;i<fab[dir].box().numPts();++i )
          {
            fab[dir].dataPtr()[i] = i;
          }
      }
  }
  private:
  int m_i;
};

// these functions used to be in Chombo/lib/src/MultiDim/SlicingI.H.transdim,
// but I removed them at some point in a bout of possibly misguided cleanup

#ifdef CH_USE1D
std::ostream&
operator<<( std::ostream& a_out, const D1::FluxBox& a_fab)
{
    a_out << a_fab.box() << "; ";
    a_out << "dim = 0: " << a_fab[0].box() << "   : ";
    int I = 1 + a_fab[0].box().hiVect()[0] - a_fab[0].box().loVect()[0];
    for ( int ii=0; ii<I; ++ii ) a_out << a_fab[0].dataPtr()[ii] << " ";
    return a_out;
}
#endif

#ifdef CH_USE2D
std::ostream&
operator<<( std::ostream& a_out, const D2::FluxBox& a_fab)
{
    a_out << a_fab.box() << "; ";
    for (int dir=0; dir<2; dir++)
      {
        a_out << "    dir = " << dir << ": ";
        a_out << a_fab[dir].box() << "; ";
        int I = 1 + a_fab[dir].box().hiVect()[0] - a_fab[dir].box().loVect()[0];
        int J = 1 + a_fab[dir].box().hiVect()[1] - a_fab[dir].box().loVect()[1];
        for ( int ii=0; ii<I*J; ++ii ) a_out << a_fab[dir].dataPtr()[ii] << " ";
      }
    return a_out;
}
#endif

#ifdef CH_USE3D
std::ostream&
operator<<( std::ostream& a_out, const D3::FluxBox& a_fab)
{
    a_out << a_fab.box() << "; ";
    for (int dir=0; dir<3; dir++)
      {
        a_out << "    dir = " << dir << ": ";
        a_out << a_fab[dir].box() << "; ";
        int I = 1 + a_fab[dir].box().hiVect()[0] - a_fab[dir].box().loVect()[0];
        int J = 1 + a_fab[dir].box().hiVect()[1] - a_fab[dir].box().loVect()[1];
        int K = 1 + a_fab[dir].box().hiVect()[2] - a_fab[dir].box().loVect()[2];
        for ( int ii=0; ii<I*J*K; ++ii ) a_out << a_fab[dir].dataPtr()[ii] << " ";
      }
    return a_out;
}
#endif


int testLevelDataSlicing()
{
  int status = 0;
  D3::DisjointBoxLayout dbl3;
  D3::Vector<D3::Box> boxes;
  const int nBoxes = 2;
  for ( int b=0;b<nBoxes;++b )
    {
      boxes.push_back( D3::Box( D3::IntVect(1+b*10,2,3),
                                D3::IntVect(2+b*10,3,4) ) );
    }
  dbl3.defineAndLoadBalance( boxes, 0 );
  D3::LevelData<D3::FArrayBox> ld3( dbl3, 1 );
  D3::LevelData<D3::FluxBox> lf3( dbl3, 1 );
  LevelDataApplyHelper< D3::FArrayBox > helper(0);
  LevelFluxApplyHelper< D3::FluxBox > fluxhelper(0);
  ld3.apply( helper );
  lf3.apply( fluxhelper);
  D2::LevelData<D2::FArrayBox> ld2;
  D2::LevelData<D2::FluxBox> lf2;
  D3::SliceSpec slicespec(1,2);
  sliceLevelData( ld2, ld3, slicespec );
  sliceLevelFlux( lf2, lf3, slicespec );
  pout()  << "Sliced \n" << ld3 << "down to...\n";
  pout()  << ld2 << '\n';

  pout()  << "Also Sliced \n" << lf3 << "down to...\n";
  pout()  << lf2 << '\n';

  pout()  << "Calling dumpLDF...\n";
  CH_2D_dumpLDFPar( &ld2 );
  pout()  << "...called dumpLDF.\n";

  D3::LevelData<D3::FArrayBox> ld3_injectee;
  D3::LevelData<D3::FluxBox> lf3_injectee;
  injectLevelData( ld3_injectee, ld2, slicespec );
  injectLevelFlux( lf3_injectee, lf2, slicespec );
  pout()  << "and injected back to...\n" << ld3_injectee << '\n';

  // Test the version where we prescribe a DBL.
  D2::DisjointBoxLayout predefinedDBL2;
  sliceDisjointBoxLayout( predefinedDBL2, ld3.disjointBoxLayout(), slicespec );

  D2::IntVect toGhost;
  sliceIntVect( toGhost, ld3.ghostVect(), slicespec );
  D2::LevelData<D2::FArrayBox> predefinedLD2( predefinedDBL2, ld3.nComp(),
                                              toGhost );
  D2::LevelData<D2::FluxBox> predefinedLDF2( predefinedDBL2, ld3.nComp(),
                                             toGhost );
  sliceLevelData( predefinedLD2, ld3, slicespec );
  sliceLevelFlux( predefinedLDF2, lf3, slicespec );
  pout()  << "predefinedLD2:\n";
  pout()  << predefinedLD2 << '\n';

  return status;
}


int testDisjointBoxLayoutSlicing()
{
  int status = 0;
    D3::DisjointBoxLayout dbl3;
    D3::Vector<D3::Box> boxes;
    const int nBoxes = 2;
    for ( int b=0;b<nBoxes;++b )
    {
        boxes.push_back( D3::Box( D3::IntVect(1+b*10,2,3),
                                  D3::IntVect(2+b*10,3,4) ) );
    }
    dbl3.defineAndLoadBalance( boxes, 0 );
    pout()  << "Starting from D3::DisjointBoxLayout:" << dbl3 << '\n';

    D2::DisjointBoxLayout dbl2;
    sliceDisjointBoxLayout( dbl2, dbl3, D3::SliceSpec( 1, 2 ) );
    pout()  << "  sliced down to:" << dbl2 << '\n';

    D3::DisjointBoxLayout dbl3_injectee;
    injectDisjointBoxLayout( dbl3_injectee, dbl2, D3::SliceSpec( 1, 2 ) );
    pout()  << "  injected back to:" << dbl3_injectee << '\n';


    D1::DisjointBoxLayout dbl1;
    sliceDisjointBoxLayout( dbl1, dbl2, D2::SliceSpec( 1,3 ) );
    pout()  << "  And further sliced:" << dbl1 << '\n';

    D2::DisjointBoxLayout dbl2_injectee;
    injectDisjointBoxLayout( dbl2_injectee, dbl1, D2::SliceSpec( 1,3 ) );
    pout()  << "  And injected back to:" << dbl2_injectee << '\n';

    return status;
}


template<typename T> int testBaseFabSlicing()
{
  int status =0;
    #define Isize 2
    #define Jsize 3
    #define Ksize 4
    D3::Box box3( D3::IntVect(1,2,3), D3::IntVect(Isize,1+Jsize,2+Ksize) );
    T rawdata[Isize*Jsize*Ksize];
    for ( int p=1;p<=Isize*Jsize*Ksize;++p )
    {
        rawdata[p-1] = T( p*(1.01) );
    }

    D3::BaseFab<T> fab3( box3, 1, rawdata );
    D3::BaseFab<T> fab3_injectee;

    for ( int dir=0;dir<3;++dir )
    {
        int slicePos = box3.loVect()[dir];
        pout()  << "Slicing " << fab3 << "\n" << "  along dir=" << dir << " "
                  << "at position " << slicePos << " and obtaining: \n";
        D2::BaseFab<T> fab2;
        sliceBaseFab( fab2, fab3, D3::SliceSpec(dir,slicePos) );
        pout()  << fab2 << '\n';

        injectBaseFab( fab3_injectee, fab2, D3::SliceSpec(dir,slicePos) );
        pout()  << "and injecting back to get: " << fab3_injectee << '\n';

        if ( dir<2 )
        {
            D2::Box box2( fab2.box() );
            slicePos = box2.loVect()[dir];
            pout()  << "...and further sliced, at position " << slicePos;
            D1::BaseFab<T> fab1;
            sliceBaseFab( fab1, fab2, D2::SliceSpec(dir,slicePos) );
            pout()  << ", down to " << fab1 << '\n';

            D2::BaseFab<T> fab2_injectee;
            injectBaseFab( fab2_injectee, fab1, D2::SliceSpec(dir,slicePos) );
            pout()  << "and injecting back to get: " << fab2_injectee << '\n';
        }
        pout()  << "-------------------------------------------\n";
    }

    return status;
}


int testMultiDimInjection()
{
  int status = 0;

  // create a 1D LevelData<FArrayBox>
  int domainSize = 8;
  D1::Box DomainBox1(D1::IntVect::Zero, (domainSize-1)*D1::IntVect::Unit);
  D1::ProblemDomain Domain1(DomainBox1);
  //Domain1.setPeriodic(0,true);

  //  int maxBoxSize = domainSize/2;
  int maxBoxSize = domainSize;

  Vector<D1::Box> Boxes1d;

  domainSplit(Domain1, Boxes1d, maxBoxSize);

  Vector<int> procAssign1d(Boxes1d.size());

  LoadBalance(procAssign1d, Boxes1d);

  D1::DisjointBoxLayout grids1d(Boxes1d, procAssign1d, Domain1);

  D1::LevelData<D1::FArrayBox> ldf1d(grids1d, 1, D1::IntVect::Unit);

  D1::DataIterator dit1d = ldf1d.dataIterator();
  for (dit1d.begin(); dit1d.ok(); ++dit1d)
    {
      D1::FArrayBox& fab = ldf1d[dit1d];
      D1::BoxIterator bit(fab.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          D1::IntVect iv = bit();
          // set to x-index
          if (Domain1.contains(iv))
            {
              // fab(iv,0) = 1.0;
              fab(iv,0) = iv[0];
            }
          else
            {
              //fab(iv,0) = 2.0;
              fab(iv,0) = iv[0];
            }
        }
    } // end loop over 1d data boxes

  // now inject this to 3d
  D2::LevelData<D2::FArrayBox> ldf2d;
  D2::SliceSpec slicespec2d(1,0);
  injectLevelData(ldf2d, ldf1d, slicespec2d);

  D3::LevelData<D3::FArrayBox> ldf3dinjectee;
  D3::SliceSpec slicespec3d(2,0);
  injectLevelData(ldf3dinjectee, ldf2d, slicespec3d);


  // now spread to entire 3d domain
  D3::Box domainBox3d(D3::IntVect::Zero, (domainSize-1)*D3::IntVect::Unit);
  D3::ProblemDomain domain3d(domainBox3d);

  //domain3d.setPeriodic(0,true);

  maxBoxSize = domainSize/2;
  Vector<D3::Box> boxes3d;
  domainSplit(domain3d, boxes3d, maxBoxSize);

  Vector<int> procAssign3d;
  LoadBalance(procAssign3d, boxes3d);

  D3::DisjointBoxLayout grids3d(boxes3d, procAssign3d, domain3d);

  D3::IntVect ghostVect(D3::IntVect::Unit);
  D3::LevelData<D3::FArrayBox> full3dData(grids3d, 1, ghostVect);

  Vector<int> transverseDir(2);
  transverseDir[0] = 1;
  transverseDir[1] = 2;

  D3::SpreadingCopier spreadCopier(ldf3dinjectee.getBoxes(), grids3d,
                                   domain3d, ghostVect, transverseDir);

  Real spreadingScale = 1.0;
  D3::SpreadingOp spreadOp(transverseDir);
  spreadOp.scale = spreadingScale;

  ldf3dinjectee.copyTo(ldf3dinjectee.interval(), full3dData,
                       full3dData.interval(), spreadCopier, spreadOp);

  full3dData.exchange();

  // now check results
  D3::DataIterator dit3d = full3dData.dataIterator();
  for (dit3d.begin(); dit3d.ok(); ++dit3d)
    {
      D3::FArrayBox& fab = full3dData[dit3d];
      const D3::Box& gridBox = grids3d[dit3d];

      D3::BoxIterator bit(gridBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          D3::IntVect iv = bit();
          //if (fab(iv,0) != 1.0)
          if (fab(iv,0) != iv[0])
            {
              status += 1;
              pout() << "   -- at location " << iv
                     << ", incorrect spreading value of "
                     << fab(iv,0) << " was found!" << endl;
            }
        }

      // now check for ghost cells which lay outside domain in the x-direction
      D3::Box loDomainBox = adjCellLo(domainBox3d, 0, 1);
      D3::Box loBox = adjCellLo(gridBox, 0, 1);
      loBox &= loDomainBox;
      loBox &= fab.box();

      if (!loBox.isEmpty())
        {
          D3::BoxIterator bitlo(loBox);
          for (bitlo.begin(); bitlo.ok(); ++bitlo)
            {
              D3::IntVect iv = bitlo();
              CH_assert(!domainBox3d.contains(iv));

              //if (fab(iv,0) != 2.0)
              if (fab(iv,0) != iv[0])
                {
                  status += 100;
                  pout() << "   -- at location " << iv
                         << ", incorrect ghost-cell spreading value of "
                         << fab(iv,0) << " was found!" << endl;
                }
            } // end loop over domain ghost cells
        } // end if there are domain ghost cells in the low direction

      D3::Box hiDomainBox = adjCellHi(domainBox3d, 0, 1);
      D3::Box hiBox = adjCellHi(gridBox, 0, 1);
      hiBox &= hiDomainBox;
      hiBox &= fab.box();

      if (!hiBox.isEmpty())
        {
          D3::BoxIterator bithi(hiBox);
          for (bithi.begin(); bithi.ok(); ++bithi)
            {
              D3::IntVect iv = bithi();
              CH_assert(!domainBox3d.contains(iv));

              //if (fab(iv,0) != 2.0)
              if (fab(iv,0) != iv[0])
                {
                  status += 100;
                  pout() << "   -- at location " << iv
                         << ", incorrect ghost-cell spreading value of "
                         << fab(iv,0) << " was found!" << endl;
                }
            } // end loop over domain ghost cells
        } // end if there are domain ghost cells in the high direction
    } // end loop over full 3d grids


  return status;
}

int main(int argc, char* argv[])
{
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
#endif

    pout() .flags(std::ios::unitbuf);

    int status = 0;

#if 0
    // haven't figured out how to compile unidim.cpp yet
    Unidimmer u;
    u.Func();
    pout()  << '\n';
#endif
    //
    // IntVect slicing
    //
    D3::IntVect iv3(1,2,3);
    D3::IntVect iv3_injectee;
    D2::IntVect iv2;
    pout()  << "******* IntVect slicing *******\n";
    for ( int dir=0;dir<3;++dir )
    {
        sliceIntVect(
            iv2, iv3, D3::SliceSpec(dir,0) );
        pout()  << "Slicing " << iv3 << ", and obtaining " << iv2 << '\n';

        injectIntVect(
            iv3_injectee, iv2, D3::SliceSpec(dir,173) );
        pout()  << "Injecting back " << iv2 << ", and obtaining "
                  << iv3_injectee << '\n';

        if ( dir<2 )
        {
            D1::IntVect iv1;
            D2::IntVect iv2_injectee;
            sliceIntVect( iv1, iv2, D2::SliceSpec(dir,0) );
            pout()  << "...and further sliced down to " << iv1 << '\n';

            injectIntVect( iv2_injectee, iv1, D2::SliceSpec(dir,371) );
            pout()  << "...and injected back to " << iv2_injectee << '\n';
        }
    }
    pout()  << "********************************\n";

    //
    // Box slicing
    //
    D3::Box box3( D3::IntVect(1,2,3), D3::IntVect(2,3,4) );
    D3::Box box3_injectee;
    D3::ProblemDomain domain3( box3 );

    pout()  << "\n******* Box slicing *******\n";
    for ( int dir=0;dir<3;++dir )
    {
        D2::Box box2;
        sliceBox( box2, box3, D3::SliceSpec(dir,2) );
        pout()  << "Slicing " << box3 << ", and obtaining " << box2 << '\n';

        injectBox( box3_injectee, box2, D3::SliceSpec(dir,2) );
        pout()  << "Injecting back " << box2 << ", and obtaining "
                  << box3_injectee << '\n';

        bool * dummy(0);
        D2::ProblemDomain domain2( sliceDomain( domain3, D3::SliceSpec(dir,2), dummy ) );
        pout()  << "Slicing ProblemDomain " << domain3.domainBox() << ", and obtaining "
                  << domain2.domainBox() << '\n';

        if ( dir<2 )
        {
            D1::Box box1;
            D2::Box box2_injectee;
            sliceBox( box1, box2, D2::SliceSpec(dir,2) );
            pout()  << "...and further sliced down to " << box1 << '\n';
            injectBox( box2_injectee, box1, D2::SliceSpec(dir,2) );
            pout()  << "...and injected back to " << box2_injectee << '\n';
        }
    }
    pout()  << "********************************\n\n";

    //
    // BaseFab slicing
    //
    int localStatus = 0;
    pout()  << "\n******* BaseFab slicing<Real> *******\n";
    localStatus = testBaseFabSlicing<Real>();
    if (localStatus != 0)
      {
        pout() << "\n BaseFab<Real> slicing FAILED!! \n";
      }
    status += localStatus;

    pout()  << "\n******* BaseFab slicing<int> *******\n";
    localStatus = testBaseFabSlicing<int>();
    if (localStatus != 0)
      {
        pout() << "\n BaseFab<int> slicing FAILED!! \n";
      }
    status += localStatus;
    pout()  << "********************************\n\n";


    //
    // DisjointBoxLayout slicing
    //
    pout()  << "\n******* DisjointBoxLayout slicing *******\n";
    localStatus = testDisjointBoxLayoutSlicing();
    if (localStatus != 0)
      {
        pout() << "\n DisjointBoxLayout slicing FAILED!! \n";
      }
    status += localStatus;

    //
    // LevelData slicing
    //
    pout()  << "\n******* LevelData slicing *******\n";
    localStatus = testLevelDataSlicing();
    if (localStatus != 0)
      {
        pout() << "\n LevelData slicing FAILED!! \n";
      }
    status += localStatus;

    //
    // true multi-dimensional injection and spreading
    pout()  << "\n******* MultiDim injection *******\n";
    localStatus = testMultiDimInjection();
    if (localStatus != 0)
      {
        pout() << "\n MultiDim injection and spreading FAILED!! \n";
      }
    status += localStatus;

#ifdef CH_MPI
    MPI_Finalize();
#endif

    return status;
}
