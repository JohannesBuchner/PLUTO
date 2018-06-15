#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
#include "SPMD.H"

#include "CH_Timer.H"

// Constants for Berger-Rigoutsos algorithm
#if ! defined(_BR_MIN_INFLECTION_MAG_)
#define _BR_MIN_INFLECTION_MAG_ ( 3 )
#endif

#define MAXBOXES     130000
#define P_BUFFERSIZE (130000 * 8 * CH_SPACEDIM)

// Berger-Rigoutsos Mesh refinement class
// ---------
//  class BRMeshRefine
//
///////////////////////////////////////////////////////////////////////////////

// Include files:

#include "BRMeshRefine.H"
#include "MayDay.H"
#include "parstream.H"
#include "NamespaceHeader.H"
//recursive function to enforce max size of boxes in a given direction
void
breakBoxes(Vector<Box>& a_vboxin,  const int& a_maxBoxSize, const int& a_idir)
{
  int nboxes = a_vboxin.size();
  //need to use STL vector for bools.
  using std::vector;
  vector<bool> splitDec(nboxes);
  bool anyBoxesToSplit = false;
  //find out which boxes need to be chopped and in what direction
  for (int ibox = 0; ibox < nboxes; ibox++)
    {
      if ( a_vboxin[ibox].size(a_idir ) > a_maxBoxSize )
        {
          splitDec[ibox] = true;
          anyBoxesToSplit = true;
        }
      else
        {
          splitDec[ibox] = false;
        }
    }
  //if there are no boxes to split, just return
  //otherwise, split all the boxes that need to be
  //split ONCE and then call function recursively
  // and set the return vector to the temporary after
  // the recursion
  if (anyBoxesToSplit)
    {
      Vector<Box> vboxtemp;
      for (int ibox = 0; ibox < nboxes; ibox++)
        {
          Box boxone = a_vboxin[ibox];
          if (splitDec[ibox])
            {
              //int len = (boxone.smallEnd(a_idir) + boxone.bigEnd(a_idir))/2;
              int mid = boxone.smallEnd(a_idir)+a_maxBoxSize ;
              Box boxtwo = boxone.chop(a_idir, mid);
              vboxtemp.push_back(boxone);
              vboxtemp.push_back(boxtwo);
            }
          else
            {
              vboxtemp.push_back(boxone);
            }
        }
      breakBoxes(vboxtemp, a_maxBoxSize, a_idir);
      a_vboxin = vboxtemp;
    }
  return;
}

BRMeshRefine::BRMeshRefine()
{
  m_fillRatio = 1.0;
}

BRMeshRefine::~BRMeshRefine()
{
}

BRMeshRefine::BRMeshRefine(const Box& a_baseDomain,
                           const Vector<int>& a_refRatios,
                           const Real a_fillRatio,
                           const int a_blockFactor,
                           const int a_bufferSize,
                           const int a_maxBoxSize)
{
  ProblemDomain crseDom(a_baseDomain);
  define(crseDom, a_refRatios, a_fillRatio, a_blockFactor,
         a_bufferSize, a_maxBoxSize);
}

BRMeshRefine::BRMeshRefine(const ProblemDomain& a_baseDomain,
                           const Vector<int>& a_refRatios,
                           const Real a_fillRatio,
                           const int a_blockFactor,
                           const int a_bufferSize,
                           const int a_maxBoxSize)
{
  define(a_baseDomain, a_refRatios, a_fillRatio, a_blockFactor,
         a_bufferSize, a_maxBoxSize);
}

void
BRMeshRefine::define(const Box& a_baseDomain,
                     const Vector<int>& a_refRatios,
                     const Real a_fillRatio,
                     const int a_blockFactor,
                     const int a_bufferSize,
                     const int a_maxBoxSize)
{
  ProblemDomain crseDom(a_baseDomain);
  define(crseDom, a_refRatios, a_fillRatio, a_blockFactor,
         a_bufferSize, a_maxBoxSize);
}

void
BRMeshRefine::define(const ProblemDomain& a_baseDomain,
                     const Vector<int>& a_refRatios,
                     const Real a_fillRatio,
                     const int a_blockFactor,
                     const int a_bufferSize,
                     const int a_maxBoxSize)
{
  // nothing new here
  MeshRefine::define(a_baseDomain, a_refRatios, a_fillRatio,
                     a_blockFactor, a_bufferSize, a_maxBoxSize);

}

///////////////////////////////////////////////////////////////////////////////
//
// Function: makeBoxes
//
// Short Description:
//  Construct a set of boxes that cover a set of tagged cells.
//
// Purpose:
//  Given a set of tagged cells defined on a single level of an AMR grid,
//  construct a BoxArray that covers all these cells that minimizes the
//  number of boxes and maximizes the ratio of tagged cells in each box.
//  This is part of the process of refining a mesh hierarchy based on error
//  estimates.
//
// Caveat:
//   The mesh that is created must lie within the proper nesting domain (Pnd)
//
// Usage Notes:
//
// Numerical Algorithm:
//  <Berger-Rigoutsos algorithm>
//
// Implementation Method:
//
// Implementation Notes:
//  This gets tricky when the minbox around a set of tags has a non-zero offset.
//  Have to be careful to remember this when looking at box sizes.
//
// References:
//  M. Berger, I. Rigoutsos, "An Algorithm for Point Clustering and Grid Generation",
//  _IEEE_Transactions_on_Systems,_Man,_and_Cybernetics_, Vol. 21, No.5, pp. 1278-1286,
//  Sept./Oct. 1991.
//
///////////////////////////////////////////////////////////////////////////////

void
BRMeshRefine::makeBoxes(Vector<Box>&         a_mesh,
                        const IntVectSet&    a_tags,
                        const IntVectSet&    a_pnd,
                        const ProblemDomain& a_domain,
                        const int            a_maxBoxSize,
                        const int            a_totalBufferSize
                        ) const
{
  CH_TIME("BRMeshRefine::makeBoxes");
  std::list<Box> boxes;

#ifdef CH_MPI
  int size;
  MPI_Comm_size (Chombo_MPI::comm, &size );
  Interval interval(0, size-1);
  makeBoxesParallel(boxes, (IntVectSet&)a_tags, a_pnd, a_domain, a_maxBoxSize, 0, a_totalBufferSize,
                    100, interval);
#else
  makeBoxes(boxes, (IntVectSet&)a_tags, a_pnd, a_domain, a_maxBoxSize, 0, a_totalBufferSize);
#endif

  //boxes.sort();
  a_mesh.resize(boxes.size());
  std::list<Box>::iterator it = boxes.begin();
  for (int i=0; i<a_mesh.size(); ++i, ++it) a_mesh[i]=*it;
}

void
BRMeshRefine::makeBoxes(std::list<Box>&      a_mesh,
                        IntVectSet&    a_tags,
                        const IntVectSet&    a_pnd,
                        const ProblemDomain& a_domain,
                        const int            a_maxBoxSize,
                        const int            a_depth,
                        const int            a_totalBufferSize
                        ) const
{
  long long int Ntags  = a_tags.numPts();
  Box minbx ;                                  //min box around Tags
  std::list<Box> mesh_hi ;                     //boxes from recursion
  IntVectSet& tags_lo = a_tags;
  IntVectSet  tags_hi ;                //tags for recursion
  a_mesh.clear();

  //
  // Handle special cases of no tags
  //
  if ( a_tags.isEmpty() )
  {
    // return null box
    return;
  }

  //
  // Validate inputs
  //

  // The number of tagged cells in a box cannot exceed the number
  // of cells in the box so enforce an upper bound on \var{FillRatio}.
  //[NOTE: 0 or negative values are allowed -- they mean any box
  //       is acceptable.  This will probably be a mistake by the
  //       caller, but there is no obvious lower valid value for
  //       this variable (1e-6 is just as likely a mistake as 0).]
  CH_assert ( m_fillRatio <= 1.0 );

  //
  // Handle the case of all tags fitting in one box.
  // This always happens at the bottom of the recursion.
  //
  Box minbox = a_tags.minBox() ;
  bool nested = properlyNested(minbox, a_domain, a_pnd, a_totalBufferSize);

  if (!nested && Ntags == 1)
  {
    CH_assert(minbox.numPts() == 1);
    return; // IntVect wasn't in PND after all (bvs)
  }

  // If minbox has enough tagged cells and is properly nested, want
  // to add this minbox to the mesh.
  // If not, continue with splitting it and recursing on the pieces.
  if ( (Ntags >= (minbox.numPts() * m_fillRatio)) && nested )
  {
    // no IntVects in the box are outside the PND so this box can
    // be accepted.  If it is larger than the maximum box size, split
    // it into two boxes
    //[NOTE: one of the boxes may violate the FillRatio requirement, but
    //       we will ignore this.]
    //  pout()<< "split Box "<<minbox<<"  maxsize "<<a_maxBoxSize<<std::endl;
    a_mesh.push_front(minbox) ;
    if (a_maxBoxSize > 0)
      {
        for (std::list<Box>::iterator it=a_mesh.begin(); it != a_mesh.end(); ++it)
        {
          splitBox( a_mesh, it, a_maxBoxSize ) ;
        }
      } // end if we are enforcing a max box size
  }
  else
  {
    // if efficiency criterion not met or box not properly nested...

    // Note tags_lo contains a_tags going in
    splitTagsInBestDimension(tags_lo, tags_hi, a_maxBoxSize);

    // Recurse on the two halves of the Tags
    if (!(tags_lo.isEmpty()))
    {
      // low interval
      makeBoxes( a_mesh, tags_lo, a_pnd, a_domain, a_maxBoxSize, a_depth+1, a_totalBufferSize);
    }
    if (!(tags_hi.isEmpty()))
    {
      // high interval
      makeBoxes( mesh_hi, tags_hi, a_pnd, a_domain, a_maxBoxSize, a_depth+1, a_totalBufferSize);
    }

    // combine the results into a single mesh
    a_mesh.splice(a_mesh.begin(), mesh_hi);
  } // done if we need to split the box

  // Done
} //end of makeBoxes

void
BRMeshRefine::makeBoxesParallel(std::list<Box>&      a_mesh,
                                IntVectSet&    a_tags,
                                const IntVectSet&    a_pnd,
                                const ProblemDomain& a_domain,
                                const int            a_maxBoxSize,
                                const int            a_depth,
                                const int            a_totalBufferSize,
                                const int            a_minSize,
                                const Interval&      a_procInterval
                                ) const
{
  if (a_procInterval.size() == 1)
    {
      makeBoxes(a_mesh, a_tags, a_pnd, a_domain, a_maxBoxSize, a_depth, a_totalBufferSize);
      return;
    }

  long long int Ntags  = a_tags.numPts();
  std::list<Box> mesh_hi ;             //boxes from recursion
  IntVectSet& tags_lo = a_tags;
  IntVectSet  tags_hi ;                //tags for recursion
  a_mesh.clear() ;

  //
  // Handle special cases of no tags
  //
  if ( a_tags.isEmpty() )
  {
    //return null box
    return;
  }

  //
  // Validate inputs
  //

  // The number of tagged cells in a box cannot exceed the number
  // of cells in the box so enforce an upper bound on \var{FillRatio}.
  //[NOTE: 0 or negative values are allowed -- they mean any box
  //       is acceptable.  This will probably be a mistake by the
  //       caller, but there is no obvious lower valid value for
  //       this variable (1e-6 is just as likely a mistake as 0).]
  CH_assert ( m_fillRatio <= 1.0 );

  //
  // Handle the case of all tags fitting in one box.
  // This always happens at the bottom of the recursion.
  //
  Box minbox = a_tags.minBox() ;
  bool nested = properlyNested(minbox, a_domain, a_pnd, a_totalBufferSize);

  if (!nested && Ntags == 1)
  {
    CH_assert(minbox.numPts() == 1);
    return; // IntVect was not in PND after all
  }

  // If minbox has enough tagged cells and is properly nested, want
  // to add this minbox to the mesh.
  // If not, continue with splitting it and recursing on the pieces.
  if ( (Ntags >= minbox.numPts() * m_fillRatio) && nested)
  {
    // no IntVects in the box are outside the PND so this box can
    // be accepted.  If it is larger than the maximum box size, split
    // it into two boxes
    //[NOTE: one of the boxes may violate the FillRatio requirement, but
    //       we will ignore this.]
    //  pout()<< "split Box "<<minbox<<"  maxsize "<<a_maxBoxSize<<std::endl;
    a_mesh.push_front(minbox) ;
    if (a_maxBoxSize > 0)
      {
        for (std::list<Box>::iterator it=a_mesh.begin(); it != a_mesh.end(); ++it)
          {
            splitBox( a_mesh,  it, a_maxBoxSize ) ;
          }
      } // end if we are enforcing a max box size
  }
  else
  {
    // if efficiency criterion not met or box not properly nested...
    //

    // Note tags_lo = a_tags going in
    splitTagsInBestDimension(tags_lo, tags_hi, a_maxBoxSize);
    //a_tags.makeEmpty();

    if (a_procInterval.size() == 1 ||
       Ntags <= a_minSize  ) // do the regular algorithm for BRMeshRefine
      {
        // Recurse on the two halves of the Tags
        if ( !tags_lo.isEmpty() )
        {
          makeBoxes( a_mesh, tags_lo, a_pnd, a_domain, a_maxBoxSize, a_depth+1, a_totalBufferSize);
        }
        if ( !tags_hi.isEmpty() )
        {
          makeBoxes( mesh_hi, tags_hi, a_pnd, a_domain, a_maxBoxSize, a_depth+1, a_totalBufferSize);
        }
      }
    else  // do makeBoxes in Parallel
      {
        //pout()<<"depth, interval "<<a_depth<<" "<<a_procInterval.begin()
        //      <<a_procInterval.end()<<std::endl;

        // first, split interval in two
        Interval lo_interval(a_procInterval.begin(),
                             (a_procInterval.end()+a_procInterval.begin()-1)/2);
        Interval hi_interval(lo_interval.end()+1, a_procInterval.end());
        // pout()<<"lo "<<lo_interval.begin()<<lo_interval.end()
        //       <<"\nhi "<<hi_interval.begin()<<hi_interval.end()<<std::endl;
        if (lo_interval.contains(procID()) && !tags_lo.isEmpty())
          {

            makeBoxesParallel( a_mesh, tags_lo, a_pnd, a_domain, a_maxBoxSize,
                               a_depth+1, a_totalBufferSize, a_minSize, lo_interval);
            sendBoxesParallel( a_mesh, a_depth);

          }
        if (hi_interval.contains(procID()) &&  !tags_hi.isEmpty()    )
          {

            makeBoxesParallel( mesh_hi, tags_hi, a_pnd, a_domain, a_maxBoxSize,
                               a_depth+1, a_totalBufferSize, a_minSize, hi_interval);
            sendBoxesParallel( mesh_hi, a_depth);

          }
        if (hi_interval.contains(procID()) &&  !tags_lo.isEmpty()    )
          {
            receiveBoxesParallel(hi_interval, lo_interval, a_mesh, a_depth);

          }
        if (lo_interval.contains(procID()) &&  !tags_hi.isEmpty()    )
          {
            receiveBoxesParallel(lo_interval, hi_interval, mesh_hi, a_depth);

          }
      }
    // combine the results into a single mesh
    a_mesh.splice(a_mesh.begin(), mesh_hi);
  }

  // Done
} //end of makeBoxesParallel

// Compute the traces (signatures) of the minbox in each direction, and
// find a hole in the trace (zero value) and an inflection point (zero
// Laplacian) in each direction; keep the best of each.
//[NOTE: both \func{find*} functions return -1 if nothing was found.]
//[NOTE: \var{infl_val} is modified by \func{findMaxInflectionPoint}.]
//[NOTE: \var{trace} indexes from 0, so indices into \var{trace} must
//       be shifted to be used as indices into the \var{a_tags}.]
// a_tags_inout_lo holds the input tags coming in, and the lo tags going out.
// maxBoxSize is the maxboxsize and is used in determining where
// to begin searching for a split index.
void BRMeshRefine::splitTagsInBestDimension(IntVectSet& a_tags_inout_lo,
                                            IntVectSet& a_tags_hi,
                                            const int a_maxBoxSize) const
{

  int hole_indx[SpaceDim], best_hole_dim;     //holes in traces
  int infl_indx[SpaceDim], best_infl_dim;     //inflection points in traces
  int infl_val [SpaceDim] ;                   //magnitudes of infl.points
  Vector<int> traces[SpaceDim];

  //  makeTraces( a_tags_inout_lo, traces) ;
  Box minbox = a_tags_inout_lo.minBox() ;
  IntVect offset = minbox.smallEnd() ;
  const IntVect& size = minbox.size();
  D_TERM6(traces[0].resize(size[0],0);,
          traces[1].resize(size[1],0);,
          traces[2].resize(size[2],0);,
          traces[3].resize(size[3],0);,
          traces[4].resize(size[4],0);,
          traces[5].resize(size[5],0););
  IntVect iv;
  IVSIterator i(a_tags_inout_lo);
  for (i.begin() ; i.ok() ; ++i )
  {
    iv = i() - offset;
    D_TERM6(traces[0][iv[0]]++,; traces[1][iv[1]]++,;traces[2][iv[2]]++,;
            traces[3][iv[3]]++,; traces[4][iv[4]]++,;traces[5][iv[5]]++);
  }

  for ( int idim = 0 ; idim < SpaceDim ; idim++ )
  {
    //hole_indx[idim] = findSplit( traces[idim] ) ;
    //infl_indx[idim] = findMaxInflectionPoint(traces[idim], infl_val[idim] ) ;
    // The following two functions, with the a_maxBoxSize argument,
    //  help balance the tag splitting by changing where to begin
    //  searching for a split or inflection index.
    hole_indx[idim] = findSplit( traces[idim], a_maxBoxSize);
    infl_indx[idim] = findMaxInflectionPoint(traces[idim], infl_val[idim], a_maxBoxSize) ;
  }
  // Take the highest index as the best one because we want to take as large
  // a box as possible  (fewer large boxes are better than many small ones)
  best_hole_dim = maxloc( hole_indx, SpaceDim ) ;
  best_infl_dim = maxloc( infl_indx, SpaceDim ) ;

  //
  // Split the Tag set at a hole in one of the traces, if there is one, or an
  // inflection point in the Laplacian of the traces.  Failing that, split
  // at the middle of the longest dimension of the enclosing box.
  int split_dim, split_index;
  if ( hole_indx[best_hole_dim] >= 0 )
  {
    // split at a hole in the trace, adjusting the trace index for the
    // offset into \var{a_tags_inout_lo} indices
    split_dim = best_hole_dim;
    split_index = hole_indx[best_hole_dim] + minbox.smallEnd(best_hole_dim);
  }
  else if ( infl_indx[best_infl_dim] >= 0 )
  {
    // split at an inflection point in the trace, adjusting the trace
    // index for the offset into \var{a_tags_inout_lo} indices
    split_dim = best_infl_dim;
    split_index = infl_indx[best_infl_dim] + minbox.smallEnd(best_infl_dim);
  }
  else
  {
    // split on the midpoint of the longest side of \var{minbox}, rounding up,
    // allowing for \var{minbox} to have a non-zero offset
    minbox.longside(split_dim); //[NOTE: split_dim is set by \func(longside)]
    split_index = (minbox.smallEnd(split_dim) + minbox.bigEnd(split_dim) +1)/2;
  }

  splitTagsInPlace( split_dim, split_index, a_tags_inout_lo, a_tags_hi );

}

///////////////////////////////////////////////////////////////////////////////
//
// Function: splitBox()
//
// Short description:
//  If the box in the given Vector pointed to by the given index is larger
//  than the given maximum size in the given dimension, split it into two
//  boxes and recurse on each.
//
///////////////////////////////////////////////////////////////////////////////

void
BRMeshRefine::splitBox( std::list<Box>&                 a_boxes,
                        const std::list<Box>::iterator& a_box,
                        const int a_dim,
                        const int a_maxBoxSize) const
{
  // See if the box needs to be split in this dimension.  If not, do nothing.
  Box& b = *a_box;
  if ( b.size( a_dim ) > a_maxBoxSize )
    {
     // pout() <<"splitting box "<<a_boxes[a_boxIndex]<<std::endl;
      // chop() returns the upper half as a new box (which is appended to
      // the vector) and modifies the original box to be the lower half
      int len =
        ( b.smallEnd( a_dim )+
          b.bigEnd(   a_dim )+1 )/2;

      std::list<Box>::iterator c = a_boxes.insert(a_box, b.chop(a_dim, len));

      // recurse on the two boxes

      splitBox( a_boxes, a_box , a_dim, a_maxBoxSize ) ;
      splitBox( a_boxes, c     , a_dim, a_maxBoxSize ) ;
    }

  // Done
}
void
BRMeshRefine::splitBox( std::list<Box>&                 a_boxes,
                        const std::list<Box>::iterator& a_box,
                        const int a_maxBoxSize) const
{
  // See if the box needs to be split in this dimension.  If not, do nothing.
  Box& b = *a_box;
  int dir;
  int longside = b.longside(dir);
  if ( longside > a_maxBoxSize )
    {
     // pout() <<"splitting box "<<a_boxes[a_boxIndex]<<std::endl;
      // chop() returns the upper half as a new box (which is appended to
      // the vector) and modifies the original box to be the lower half
      int len =
        ( b.smallEnd( dir )+
          b.bigEnd(   dir )+1 )/2;

      std::list<Box>::iterator c = a_boxes.insert(a_box, b.chop(dir, len));

      // recurse on the two boxes

      splitBox( a_boxes, a_box , a_maxBoxSize ) ;
      splitBox( a_boxes, c     , a_maxBoxSize ) ;
    }

  // Done
}

///////////////////////////////////////////////////////////////////////////////
//
// Function: makeTrace()
//
// Short description:
//  Count the number of cells in the IVS at every index in the specified
//  coordinate direction.  In other words, make a histogram using the
//  \var{Dir}'th coordinate as data.
//
// Note: Vectors index from 0 so the trace indices have to be shifted from
//       the box indices by the value of the lower bound.
//
// References:
//  The trace computed here is the \Sigma value computed in Berger-Rigoutsos.
//
///////////////////////////////////////////////////////////////////////////////

Vector<int>
BRMeshRefine::makeTrace( const IntVectSet& a_Ivs, int a_dir ) const
{
  //
  // Validate inputs
  //
  CH_assert( a_dir >= 0 && a_dir < SpaceDim ) ;

  //
  // Histogram the coordinates of \var{a_Ivs} in direction \var{a_dir} into
  // \var{trace}.  The index range of \var{trace} must start at 0
  // regardless of the bounds of the minbox of \var{a_Ivs}
  //
  Box minbox = a_Ivs.minBox() ;
  int offset = minbox.smallEnd(a_dir) ;
  Vector<int> trace( minbox.size(a_dir), 0 ) ;  //initialize to 0
  IVSIterator i(a_Ivs) ;

  for ( i.begin() ; i.ok() ; ++i )
  {
    trace[i()[a_dir]-offset] += 1 ;
  }

  // Done
  return( trace ) ;
}

void BRMeshRefine::makeTraces ( const IntVectSet& a_Ivs, Vector<int>* traces) const
{
  Box minbox = a_Ivs.minBox() ;
  IntVect offset = minbox.smallEnd() ;
  const IntVect& size = minbox.size();
  D_TERM6(traces[0].resize(size[0],0);,
          traces[1].resize(size[1],0);,
          traces[2].resize(size[2],0);,
          traces[3].resize(size[3],0);,
          traces[4].resize(size[4],0);,
          traces[5].resize(size[5],0););
  IVSIterator i(a_Ivs) ;
  for ( i.begin() ; i.ok() ; ++i )
  {
    IntVect iv = i();
    iv-=offset;
    D_TERM6(traces[0][iv[0]]++,; traces[1][iv[1]]++,;traces[2][iv[2]]++,;
            traces[3][iv[3]]++,; traces[4][iv[4]]++,;traces[5][iv[5]]++);
  }

}

///////////////////////////////////////////////////////////////////////////////
//
// Function: findSplit()
//
// Short Description:
//  Given a trace (ie. signature) of a Tag set, find a place in the trace
//  with a zero count.  Assumes the first and last elements of the trace
//  are non-zero.
//
// Usage Note:
//  If there are no acceptable zeros in the trace, this returns -1.
//  \var{Vector}s index from 0, so the result may have to be shifted to be
//  used as an index into the \var{IntVectSect} which produced the trace.
//
///////////////////////////////////////////////////////////////////////////////

int
BRMeshRefine::findSplit( const Vector<int>& a_trace ) const
{
  // look for a place to split at -- begin at the index after the first nonzero
  for ( int i = 1 ; i < a_trace.size() -1 ; i++ )
  {
    if ( a_trace[i] == 0 ) return( i ) ;
  }

  // nothing acceptable
  return( -1 ) ;
}

// If the traces are large enough, begin looking for a split index in the center
int
BRMeshRefine::findSplit( const Vector<int>& a_trace, const int a_maxBoxSize ) const
{
  // If length of trace is larger than 2 times the maxboxsize (arbitrary),
  // do somthing different
  if (a_trace.size() > 2*a_maxBoxSize)
  {
    // look for a place to split at -- begin in the middle, count down to 1,
    // look in both directions
    int ibegin = a_trace.size()/2;
    int iforward = ibegin+1;
    for ( int i = ibegin ; (i >= 1 && iforward < a_trace.size()-1 ); i-- )
    {
      if ( a_trace[i] == 0 ) return( i );
      if ( a_trace[iforward] == 0 ) return( iforward );
      iforward++;
    }
  }
  else
  {
    int firstHole = findSplit(a_trace); // perform original algorithm
    return(firstHole);
  }

  // nothing acceptable
  return( -1 ) ;
}

///////////////////////////////////////////////////////////////////////////////
//
// Function: findMaxInflectionPoint
//
// Short Description:
//  Find the largest inflection point in the given trace (ie. signature)
//  and return the index and set the value of the inflection in the last
//  argument.  The index will be the edge-centered location where a
//  chop should be made in the BR algorithm.
//
//  If there are no inflection points or if none of them are larger than the
//  cutoff tolerance, -1 is returned.
//
// Implementation Note:
//  The smallest acceptable inflection magnitude is the global constant
//  _MIN_INFLECTION_MAG_
//
///////////////////////////////////////////////////////////////////////////////

int
BRMeshRefine::findMaxInflectionPoint(const Vector<int>& a_trace,
                                     int& a_maxVal) const
{
  // first compute the discrete Laplacian of the trace
  Vector<int> d2Trace( a_trace.size(), 0 ) ;
  for ( int i = 1 ; i < a_trace.size()-1 ; i++ )
  {
    d2Trace[i] = a_trace[i-1] - 2*a_trace[i] + a_trace[i+1] ;
  }

  // find inflection points and save one with the largest magnitude
  int absval, imax=0;
  a_maxVal = -1 ;
  for ( int i = 2 ; i < a_trace.size()-1 ; i++ )
  {
    absval = abs( d2Trace[i-1]  - d2Trace[i] ) ;
    if ( d2Trace[i-1] * d2Trace[i] < 0  && absval > a_maxVal )
    {
      imax = i ; a_maxVal = absval ;
    }
  }

  // Find the largest inflection point, if one exists and
  // has magnitude large enough and return its index.
  // return edge-centered location of chop point.
  if ( a_maxVal == -1 )                          return( -1 ) ;
  else if ( a_maxVal < _BR_MIN_INFLECTION_MAG_ ) return( -1 ) ;
  else                                        return( imax);
}

// If the traces are large enough, begin looking for a split index in the center
int
BRMeshRefine::findMaxInflectionPoint(const Vector<int>& a_trace,
                                     int& a_maxVal,
                                     const int a_maxBoxSize) const
{
  int imax=0;
  // If length of trace is larger than 2 times the maxboxsize (arbitrary),
  // do somthing different
  if (a_trace.size() > 2*a_maxBoxSize)
    {
      // first compute the discrete Laplacian of the trace
      Vector<int> d2Trace( a_trace.size(), 0 ) ;
      for ( int i = 1 ; i < a_trace.size()-1 ; i++ )
        {
          d2Trace[i] = a_trace[i-1] - 2*a_trace[i] + a_trace[i+1] ;
        }

      // find inflection points and save one with the largest magnitude
      int absval;
      a_maxVal = -1;

      // look for inflection -- begin in the middle, count down to 1,
      // look in both directions
      int ibegin = a_trace.size()/2;
      int iforward = ibegin+1;
      for ( int i = ibegin ; (i >= 2 && iforward < a_trace.size()-1); i-- )
        {
          absval = abs( d2Trace[i-1] - d2Trace[i] ) ;
          if ( d2Trace[i-1] * d2Trace[i] < 0  && absval > a_maxVal )
          {
            imax = i ;
            a_maxVal = absval ;
          }

          // same thing in the forward direction
          absval = abs( d2Trace[iforward-1] - d2Trace[iforward] ) ;
          if ( d2Trace[iforward-1] * d2Trace[iforward] < 0  && absval > a_maxVal )
          {
            imax = iforward ;
            a_maxVal = absval ;
          }
          iforward++;
        }
    }
  else
    { // call original algorithm
      int inflectionPoint = findMaxInflectionPoint(a_trace, a_maxVal);
      return(inflectionPoint);
    }

  // Find the largest inflection point, if one exists and
  // has magnitude large enough and return its index.
  // return edge-centered location of chop point.
  if ( a_maxVal == -1 )                          return( -1 ) ;
  else if ( a_maxVal < _BR_MIN_INFLECTION_MAG_ ) return( -1 ) ;
  else                                        return( imax);
}

///////////////////////////////////////////////////////////////////////////////
//
// Function: splitTags
//
// Purpose:
//  Given a direction and an index in that direction, split the IntVectSet of
//  tags into two sets (lo,hi) on either side of the index.  Tags that match
//  the index exactly go into the hi set.
//
// Implementation Method:
//  For each tag in the set, look at the coordinate in the direction specified,
//  and if it is <= the specified index, copy the tag into the lo set.  Else
//  copy it into the hi set.
//
// Usage Note:
//  It is acceptable for the split index to be outside of the box containing Tags.
//  This forces all the Tags to end up in one of the two output sets.
//
///////////////////////////////////////////////////////////////////////////////

void
BRMeshRefine::splitTags(const IntVectSet& a_tags,
                        const int a_split_dir, const int a_split_indx,
                        IntVectSet& a_tags_lo, IntVectSet& a_tags_hi ) const
{
  // Validate inputs
  //[NOTE: it is ok if the a_split_indx is outside the box containing a_tags.]
  CH_assert( a_split_dir >= 0 && a_split_dir < SpaceDim );
  CH_assert( !a_tags.isEmpty() );

  // Copy the whole set into a_tags_lo and extract the cells that aren't
  // below the split index into a_tags_hi, keeping the leftover in a_tags_lo
  // (i.e., the a_split_indx'th cells go into the high set).
  a_tags_lo = a_tags ;
  a_tags_hi = a_tags_lo.chop( a_split_dir, a_split_indx );

}

void BRMeshRefine::splitTagsInPlace(const int a_split_dir, const int a_split_indx,
                                    IntVectSet& a_tags_inout_lo,
                                    IntVectSet& a_tags_hi) const
{
  a_tags_inout_lo.chop(a_split_dir, a_split_indx, a_tags_hi);
}

//recursive function to enforce max size of boxes in a given direction
void
BRMeshRefine::breakBoxes(Vector<Box>& a_vboxin,
                         const int& a_maxBoxSize, const int& a_idir) const
{
  int nboxes = a_vboxin.size();
  //need to use STL vector for bools.
  using std::vector;
  vector<bool> splitDec(nboxes);
  bool anyBoxesToSplit = false;
  //find out which boxes need to be chopped and in what direction
  for (int ibox = 0; ibox < nboxes; ibox++)
    {
      if ( a_vboxin[ibox].size(a_idir ) > a_maxBoxSize )
        {
          splitDec[ibox] = true;
          anyBoxesToSplit = true;
        }
      else
        {
          splitDec[ibox] = false;
        }
    }
  //if there are no boxes to split, just return
  //otherwise, split all the boxes that need to be
  //split ONCE and then call function recursively
  // and set the return vector to the temporary after
  // the recursion
  if (anyBoxesToSplit)
    {
      Vector<Box> vboxtemp;
      for (int ibox = 0; ibox < nboxes; ibox++)
        {
          Box boxone = a_vboxin[ibox];
          if (splitDec[ibox])
            {
//               int len = (boxone.smallEnd(a_idir) +
//                          boxone.bigEnd(a_idir))/2;
              int len = boxone.smallEnd(a_idir)+a_maxBoxSize-1;
              Box boxtwo = boxone.chop(a_idir, len);
              vboxtemp.push_back(boxone);
              vboxtemp.push_back(boxtwo);
            }
          else
            {
              vboxtemp.push_back(boxone);
            }
        } // end loop over boxes
      breakBoxes(vboxtemp, a_maxBoxSize, a_idir);
      a_vboxin = vboxtemp;
    }
  return;
}

///////////////////////////////////////////////////////////////////////////////
//
// Function: maxloc
//
// Purpose: find the maximum value in a Vector<int> or an int[] and return its index.
//
// Usage Notes:
//  Returns -1 if the Vector has no entries.
//
///////////////////////////////////////////////////////////////////////////////

int
BRMeshRefine::maxloc( const int* a_V, const int a_Size ) const
{
  int imax = 0 ;
  for ( int i=1 ; i<a_Size ; i++ ) if ( a_V[i] > a_V[imax] ) imax = i ;
  return( imax ) ;
}

// at the moment, this is a pretty hokey approach
// (first build a BRMeshRefine, then "generate" grids
// hopefully will have a better version soon...
void
domainSplit(const ProblemDomain& a_domain, Vector<Box>& a_vbox,
            int a_maxBoxSize, int a_blockFactor)
{
  const Box& domBox = a_domain.domainBox();
  domainSplit(domBox, a_vbox, a_maxBoxSize, a_blockFactor);
}
void
domainSplit(const Box& a_domain, Vector<Box>& a_vbox, int a_maxBoxSize, int a_blockFactor)
{
  a_vbox.resize(0);
  if (a_maxBoxSize == 0)
    {
      a_vbox.push_back(a_domain);
      return;
    }
  int ratio = a_maxBoxSize/a_blockFactor;

  Box d(a_domain);
  d.coarsen(a_blockFactor);
  if (refine(d, a_blockFactor) != a_domain)
    {
      MayDay::Error("domainSplit: a_domain not coarsenable by blockingFactor");
    }
  a_vbox.push_back(d);
  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      breakBoxes(a_vbox, ratio, i);
    }
  for (int i=0; i<a_vbox.size(); ++i)
    {
      a_vbox[i].refine(a_blockFactor);
    }
}

// void
// domainSplit(const Box& a_domain, Vector<Box>& a_vbox, int a_maxBoxSize, int a_blockFactor)
// {
//   CH_assert(!a_domain.isEmpty());
//   CH_assert(a_maxBoxSize > 0);
//   int iref = 1;
//   Box coardom = coarsen(a_domain, iref);
//   Vector<Box> domains(2);
//   domains[0] = coardom;
//   domains[1] = a_domain;
//   Vector<Box> oldGrid0(1, coardom);
//   Vector<Vector<Box> > oldMeshes(1);
//   oldMeshes[0] = oldGrid0;
//   if (a_maxBoxSize/a_blockFactor/iref < 1)
//     MayDay::Error("DomainSplit: maxsize and blocking factor incompatible");
//   Vector<Vector<Box> > newMeshes;
//   Vector<int> refratio(1,iref);
//   IntVectSet tags(coardom);
//   // make irreg_boxes using berger-rigoutsos
//   Real fillrat = 0.7;
//   int buffersize = 1;
//   int block_factor= a_blockFactor;
//   Box testdom = coarsen(a_domain, iref*block_factor);
//   testdom.refine(iref*block_factor);
//   if (a_domain != testdom)
//     MayDay::Error("DomainSplit:domain and Blocking Factor incompatible");
//   BRMeshRefine mrobject(domains[0], refratio, fillrat,
//                         block_factor, buffersize, a_maxBoxSize);

//   mrobject.regrid(newMeshes, tags, 0, 0,
//                                   oldMeshes);

//   CH_assert(newMeshes.size() == 2);
//   a_vbox = newMeshes[1];
// }

std::list<int*> sendBuffers;
std::list<int>  ch_count;

void
BRMeshRefine::sendBoxesParallel( const std::list<Box>& a_mesh,
                                 int tag) const
{
#ifdef CH_MPI
  const int boxSize = 2*CH_SPACEDIM;
  int numBox = a_mesh.size();
  int* next = new int[numBox*boxSize];
  sendBuffers.push_back(next);
  ch_count.push_back(numBox*boxSize);
  std::list<Box>::const_iterator it = a_mesh.begin();
  for (; it!=a_mesh.end(); ++it, next+=boxSize)
    {
      Box b = *it;
      D_TERM6(next[0]=b.smallEnd(0);next[1]=b.bigEnd(0);,
              next[2]=b.smallEnd(1);next[3]=b.bigEnd(1);,
              next[4]=b.smallEnd(2);next[5]=b.bigEnd(2);,
              next[6]=b.smallEnd(3);next[7]=b.bigEnd(3);,
              next[8]=b.smallEnd(4);next[9]=b.bigEnd(4);,
              next[10]=b.smallEnd(5);next[11]=b.bigEnd(5););
    }
#endif
}

void
BRMeshRefine::receiveBoxesParallel(const Interval& a_from,
                                   const Interval& a_to,
                                   std::list<Box>& a_mesh,
                                   int tag) const
{
#ifdef CH_MPI

  // pout()<<"from "<<a_from.begin()<<a_from.end()<<"\n"
  //      <<"to   "<<a_to.begin()<<a_to.end()<<"\n";
  // NOTE: This static malloc is currently not being freed anywhere (ndk)
  static int* recBuffer = (int*)mallocMT(P_BUFFERSIZE);
  int* next = recBuffer;
  MPI_Status status;
  const int boxSize = 2*CH_SPACEDIM;

  int source, dest;

  source = procID()-a_from.begin() + a_to.begin();
  dest = source;

  bool hang=false;
  if (a_from.size() > a_to.size() && procID()== a_from.end()) hang = true;
  if (a_to.size() > a_from.size() && procID()== a_to.end()) hang = true;

  if (a_from.size() == a_to.size() || !hang)
    {
      //pout()<<"expecting boxes from "<<source<<"\n";
      //pout()<<"sending "<<ch_count.back()/boxSize<<" boxes to "<<dest<<std::endl;
      MPI_Sendrecv( sendBuffers.back(), ch_count.back(), MPI_INT,
                    dest, tag,
                    recBuffer, MAXBOXES*boxSize, MPI_INT,
                    source, tag, Chombo_MPI::comm, &status );
    }
  // OK, need to pick up the oddball stragler here.
  if (a_from.size() < a_to.size() && procID()==a_from.end())
    {
      dest = a_to.end();
      //pout()<<"SEnding "<<ch_count.back()/boxSize<<" boxes to "<<dest<<std::endl;
      MPI_Send(sendBuffers.back(), ch_count.back(), MPI_INT,
               a_to.end(), tag, Chombo_MPI::comm);
    }
  if (a_from.size() > a_to.size() && procID()==a_from.end())
    {
      source = a_to.end();
      //pout()<<"EXpecting boxes from "<<source<<"\n";
      MPI_Recv(recBuffer, MAXBOXES*boxSize, MPI_INT,
               source, tag, Chombo_MPI::comm, &status );
    }

  delete [] sendBuffers.back();
  sendBuffers.pop_back();
  ch_count.pop_back();

  int numInt;
  MPI_Get_count(&status, MPI_INT, &numInt);

  int* end = next+numInt;
  Box b;
  IntVect& lo = (IntVect&)(b.smallEnd());
  IntVect& hi = (IntVect&)(b.bigEnd());
  for (; next<end; next+=(2*CH_SPACEDIM))
    {
      D_TERM6(lo[0]=next[0];hi[0]=next[1];,
              lo[1]=next[2];hi[1]=next[3];,
              lo[2]=next[4];hi[2]=next[5];,
              lo[3]=next[6];hi[3]=next[7];,
              lo[4]=next[8];hi[4]=next[9];,
              lo[5]=next[10];hi[5]=next[11];);

      b.computeBoxLenNotEmpty();
      a_mesh.push_back(b);
    }
  //pout()<<"received "<<a_mesh.size()<<" boxes from "<<source<<std::endl;
#endif
}
#include "NamespaceFooter.H"
