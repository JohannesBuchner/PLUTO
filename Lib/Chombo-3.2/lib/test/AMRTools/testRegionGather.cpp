#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif



#include "DebugDump.H"
#include "RegionGather.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "MayDay.H"
#include "UsingNamespace.H"

/*****************/
/*****************/

/// Global variables for handling output:
static const char *pgmname = "testRegionGather" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = false ;
/*****************/
/*****************/

int testRegionGather();

void buildLayout(DisjointBoxLayout& layout, const Box& domain, int boxSize)
{

  Vector<Box> boxes;
  domainSplit(domain, boxes, boxSize, boxSize);
  Vector<int> procs;
  LoadBalance(procs, boxes);
  layout.define(boxes, procs);
}



int
main(int argc, char *argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int eekflag = testRegionGather();
  if (eekflag != 0)
    {
      pout() << indent << pgmname << ": failed with error code: "
           << eekflag << endl;
    }
  else
    {
      pout() << indent << pgmname << ": passed. " << endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}

// trying to write this test to work in 1D as well...we'll see.

int testRegionGather()
{

  IntVect low = 32*IntVect::Unit;
  IntVect hi  = 79*IntVect::Unit;
  int boxSize = 8;

  Box domain(low, hi);
  DisjointBoxLayout layout;
  buildLayout(layout, domain, boxSize);

  ProblemDomain pDomain(domain);

  RegionGather gather;
  gather.define(pDomain, layout, 24);


  LayoutData<Real> monopoles(layout);

  DataIterator dit = monopoles.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& b = layout.get(dit);
      monopoles[dit] = b.smallEnd()[0] + (b.bigEnd()[0] - b.smallEnd()[0])/2;

    }

  //LayoutData<Vector<GatherObject<Real> > > monopoleMoments(layout);
  LayoutData<Vector<GatherObject<Real> > > monopoleMoments;

  // gather.dump();

  regionGather(monopoles, gather, monopoleMoments);

  if (verbose)
  {
    for (dit.begin(); dit.ok(); ++dit)
      {
        const Vector<GatherObject<Real> >& obj = monopoleMoments[dit];
        pout()<<layout.get(dit);
        for (int i=0; i<obj.size(); ++i)
          {
            pout()<<" "<<obj[i].m_value<<" "<<obj[i].m_offset;
          }
        pout()<<"\n\n";
      }
  }

  return 0;
}
