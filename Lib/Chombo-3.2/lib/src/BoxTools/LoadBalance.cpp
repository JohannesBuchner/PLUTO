#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <list>
#include <set>
using std::cout;

#include "parstream.H"

#include "DataIterator.H"
#include "Misc.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "LayoutIterator.H"
#include "CH_Timer.H"

// Write a text file per call to LoadBalance()
//#define PRINT_EXTRA_LB_FILE
#ifdef PRINT_EXTRA_LB_FILE
#include <fstream>
using std::fstream;
#endif

#include "NamespaceHeader.H"

// local prototypes
int
min_element( const Vector<long>& Vect );
void
min_max_elements( int& imin ,int& imax ,const Vector<long>& Vect );

void
min_max_elements( int& imin ,int& imax ,const Vector<long long>& Vect );

// Local class definition (needed to sort loads)
class Load
{
public:
  Load():load(0) ,grid_index(0)
  {
  }

  bool operator < (const Load& rhs) const
  {
    return load < rhs.load;
  }

  long load;      //actual load on this box
  int grid_index; //link to Grids[]
};

// Code:

///
// This version takes a Vector<BoxLayout> and builds a matching Vector of
// Vector<Box> and uses it to call the full version.
///
int
LoadBalance( Vector<BoxLayout>&    Grids            //in-out: input grids to balance
             ,Real&                 effRatio         //output: ratio of min load
             ,const Vector<Vector<long> >&  ComputeLoads     //input: computational cost
             ,const Vector<int>&           RefRatios        //input: refinement ratio
             ,int nProc
             )
{
  CH_TIME("LoadBalance:VectorBoxLayout");
  Vector<Vector<Box> > boxes( Grids.size() );

  for ( int i = 0; i < Grids.size(); ++i )
    {
      for ( LayoutIterator j = Grids[i].layoutIterator() ; j.ok() ; ++j )
        {
          boxes[i].push_back( Grids[i][j()] );
        }
    }

  Vector<Vector<int> > procIDs( Grids.size() );
  int status = LoadBalance( procIDs ,effRatio
                            ,boxes ,ComputeLoads ,RefRatios, nProc );

  if ( status < 0 ) return status;  //LoadBalance() failed

  for ( int i = 0; i < Grids.size(); ++i )
    {
      LayoutIterator j = Grids[i].layoutIterator();
      int p;
      for ( j.reset() ,p = 0; j.ok(); ++j ,++p )
        {
          Grids[i].setProcID( j() ,procIDs[i][p] );
        }
    }

  return status;
}

///
// This version takes a single BoxLayout (i.e. one level)
// and uses the box volumes to construct the compute loads
// and calls the full version.
// a_LBnumProc specifies number of procs to use in LoadBalancing. defaults to numProc().
///
int LoadBalance(Vector<int>& a_procAssignments, const Vector<Box>& a_boxes,
                const int a_LBnumProc)
{
  CH_TIME("LoadBalance:VectorBoxEntry");
  Vector<long long> computeLoads(a_boxes.size());

  for (int i = 0; i < a_boxes.size(); ++i)
    {
      computeLoads[i] = a_boxes[i].numPts();
    }

  int status;
  status = LoadBalance(a_procAssignments, computeLoads, a_boxes, a_LBnumProc);

  // save this older commented code
  //   Vector<Vector<Box> > layouts(1, boxes);
  //   Vector<Vector<long> > computeLoads(1, Vector<long>(boxes.size()));
  //   Vector<int> refRatios(1,1);
  //   Vector<Vector<int> > assignments(1,Vector<int>(boxes.size(), -1));
  //   Real effRatio;

  //   for (int index = 0; index < boxes.size(); ++index)
  //     {
  //       computeLoads[0][index] = layouts[0][index].numPts();
  //     }

  //   int ret = LoadBalance(assignments, effRatio, layouts, computeLoads, refRatios);

  //   if (ret == 0)
  //     procs = assignments[0];

  return status;
}

///
// Accepts "long int" computeLoads and then just calls the "long long" version.
///
int LoadBalance(Vector<int>&             a_procAssignments,
                const Vector<long int>&  a_computeLoads,
                const Vector<Box>&       a_boxes,
                const int                a_LBnumProc)
{
  Vector<long long> cloads(a_computeLoads.size());
  for (int i=0; i<a_computeLoads.size(); i++)
    {
      cloads[i] = (long long)a_computeLoads[i];
    }
  int status = LoadBalance(a_procAssignments, cloads, a_boxes, a_LBnumProc);
  return status;
}

///
// This version takes compute loads directly and does the "simple"
// load balance algorithm (modified knapsack algorithm).
// Boxes are also included, but only used if SWAP_SAME_SIZES_BOXES
// is enabled
///
int LoadBalance(Vector<int>&             a_procAssignments,
                const Vector<long long>& a_computeLoads,
                const Vector<Box>&       a_boxes,
                const int                a_LBnumProc)
{
  CH_TIME("LoadBalance:VectorBoxSimple");
  // Phase 1  modified knapsack algorithm.  this one doesn't use
  // load sorting followed by round robin.  This one does bin packing
  // first by finding vector divisors, then does regular knapsack
  // optimization.

  int status = 0;
  a_procAssignments.resize(0);
  a_procAssignments.resize(a_computeLoads.size(), 0);

  const int Nboxes = a_procAssignments.size();

  if (a_LBnumProc == 1 || Nboxes == 1)
  {
    return 0;
  }

  // first, compute 'total load'
  double totalLoad = 0;
  for (int ibox = 0; ibox < Nboxes; ++ibox)
    {
      totalLoad += a_computeLoads[ibox];
      a_procAssignments[ibox] = -ibox;
    }

  // figure out roughly single processor load goal
  double goal = totalLoad/a_LBnumProc;

  // debugging prints
  //pout()<<"--- LoadBalance  --- Number of Boxes: " << Nboxes << std::endl;
  //pout()<<"       a_computeLoads.size()= "<< a_computeLoads.size() << std::endl;
  //pout()<<"       a_boxes.size()=        "<< a_boxes.size() << std::endl;
  //pout() << " Box sizes: " << std::endl;
  //for (int i=0; i < Nboxes; i++) pout() <<" "<<a_computeLoads[i];
  //pout() << std::endl;

  // now, first bin assignments

  // change goal to dynamicgoal (ndk)
  // We desire equal amount of cells per processor, which
  // is the total number of cells divided by the number of processors: ie the goal.
  // But because the cells are grouped into boxes of sometimes very different sizes,
  // and because we assume that the boxes are in a locality-preferable vector,
  // using a constant goal can allow for several consecutive processors to be
  // assigned a number of cells _less_ than (or _more_ than) the goal.
  // This can lead to the very last processor receiving an abnormally large
  // (or small) workload.  Changing the algorithm to compute a new dynamic
  // goal after each processor-box assignment eliminates this problem.
  // However, there are still improvements that can be made when the
  // number of boxes to be farmed out is not much larger than the total number
  // of processors doing the work.  (ndk)

  Vector<long long> loads(a_LBnumProc, 0);
  int bin = 0;
  double remainingLoad = totalLoad;
  double dynamicGoal = goal;
  long long localLoad = 0;
  //pout() << " Begin loop over boxes: remainingLoad =" << remainingLoad
  //     << " dynamicGoal = " << dynamicGoal << " a_LBnumProc=" << a_LBnumProc
  //     << " Nboxes=" << Nboxes << std::endl;
  bool NboxesLeftEqualsNprocsLeft=false;
  //char ctmp[205];
  for (int ibox=0; ibox < Nboxes; ibox++)
    {
      // If the number of boxes left to sort out is <= to the number of proc bins left,
      // then force each remaining proc to have one box -- by setting NboxesLeftEqualsNprocsLeft=true
      if ( Nboxes-ibox <= a_LBnumProc-bin)
        {
          // Once set true, remains true for the rest of the boxes.
          NboxesLeftEqualsNprocsLeft = true;
        }

      if (NboxesLeftEqualsNprocsLeft)
        {
          // a_procAssignments initialized to -ibox above
          // Already something in this bin, so go to next bin.
          //  Note: Only need to set loads[bin] here for any diagnostics after loop
          if (a_procAssignments[ibox] > 0) bin++;

          if (bin > a_LBnumProc)
          {
            // Not sure this could even happen...
            MayDay::Abort("Problem in LoadBalance");
          }
          a_procAssignments[ibox] = bin;
          loads[bin] += a_computeLoads[ibox];
          bin++;
        }
      else
        {

          if (bin != a_LBnumProc-1)  // NOT last rank
            {
              localLoad += a_computeLoads[ibox];
              if (localLoad >= dynamicGoal)
                {
                  if ( ( localLoad - dynamicGoal ) > a_computeLoads[ibox]/2
                      && loads[bin]!=0)//added this last condition so no-load procs don't get skipped. mfb
                    {
                      a_procAssignments[ibox] = bin+1;
                      loads[bin+1] += a_computeLoads[ibox];
                      localLoad = a_computeLoads[ibox];
                    }
                  else
                    {
                      a_procAssignments[ibox] = bin;
                      loads[bin] += a_computeLoads[ibox];
                      localLoad = 0;
                    }
                  remainingLoad -= (double)loads[bin];
                  dynamicGoal = remainingLoad/(double)(a_LBnumProc-(bin+1));
                  bin++;
                }
              else
                {
                  loads[bin] += a_computeLoads[ibox];
                  a_procAssignments[ibox] = bin;
                }
            }
          else
            { // last rank
              a_procAssignments[ibox] = bin;
              loads[bin] += a_computeLoads[ibox];
            }
        }

      //sprintf(ctmp, " box#%04d localLoad=%10lld computeLoads[ibox]=%10lld remainingLoad=%12.1f dynamicGoal=%12.1f a_procAssignments[i]=%4d\n",
      //  ibox, localLoad, a_computeLoads[ibox], remainingLoad, dynamicGoal, a_procAssignments[ibox]);
      //pout() << ctmp;
    }

#ifdef PRINT_EXTRA_LB_FILE
  // Write a file per call to LoadBalance
  static int NglobalLBcall=0;
  if (procID() == 0)
  {
    char temp[1024];
    sprintf(temp, "LB%03d.txt", NglobalLBcall);
    NglobalLBcall++;
    fstream FILE(temp, ios::out);

    FILE << "  LoadBalance  sizes: a_computeLoads=" << a_computeLoads.size()
         << " a_boxes=" << a_boxes.size()
         << " a_LBnumProc=" << a_LBnumProc  << "\n";

    //int NminNumberBoxesToDiagPrint = a_LBnumProc + a_LBnumProc + a_LBnumProc/2;
    //int NminNumberBoxesToDiagPrint = a_LBnumProc + a_LBnumProc/2;
    int NminNumberBoxesToDiagPrint = 0;

    // first , scan for min/max processor load
    long long maxload = loads[0];
    long long minload = loads[0];
    int rankmin=0, rankmax=0;
    for (int i=1; i<loads.size(); ++i)
    {
      if (loads[i] > maxload)
      {
        maxload=loads[i];
        rankmax=i;
      }
      if (loads[i] < minload)
      {
        minload=loads[i];
        rankmin=i;
      }
    }
    sprintf(temp, " minload=%12lld on rank%4d   maxload=%12lld on rank%4d\n",
            minload, rankmin, maxload, rankmax);
    FILE << temp;

    //parallel inefficiency: sum of idle time of processors waiting for max load processor
    // dividedd by total load
    double eff = 1 - ((double)(maxload*a_LBnumProc - totalLoad))/totalLoad;
    FILE << "loadbalance efficiency = " << eff << "  Nboxes=" << Nboxes << std::endl;

    // only print if there's something interesting to see...
    if (Nboxes > NminNumberBoxesToDiagPrint)
    {
      for (int i=0; i<a_LBnumProc; i++)
      {
        int NboxCount=0;
        for (int jbox=0; jbox < Nboxes; jbox++)
        {
          // count up number of boxes on proc i
          if (a_procAssignments[jbox] == i)
          {
            NboxCount++;
          }
        }

        int bigbox=0, totalpoints=0;
        for (int jbox=0; jbox < Nboxes; jbox++)
        {
          if (a_procAssignments[jbox] == i)
          {
            totalpoints += a_computeLoads[jbox];
            if (a_computeLoads[jbox] > bigbox) bigbox=a_computeLoads[jbox];
          }
        }
        bool allSameSize=false;
        // need to fix this -- cuz this is only per proc and
        // what i want is no printintg if every box on all procs
        // are the same size...
        //if (totalpoints/NboxCount == Nboxes) allSameSize=true;

        bool printEveryBox = false;
        if (allSameSize)
        {
          // don't print
          FILE << NboxCount << "  boxes of same size" << Nboxes << std::endl;
        }
        else
        {
          char tempstring2[1024];
          sprintf(tempstring2, "rank%04d load=%-11lld boxes(%4d) [%6d] ",
                  i, loads[i], NboxCount, bigbox);
          FILE << tempstring2;
          //pout() << "rank" << i << "  load=" << (long)loads[i] << " boxes("
          //     << NboxCount << ")=";
          if (printEveryBox)
          {
            for (int jbox=0; jbox < Nboxes; jbox++)
            {
              if (a_procAssignments[jbox] == i)
              {
                FILE << a_computeLoads[jbox] << " ";
              }
            }
          }
          FILE << std::endl;
        }
      }
    }
    FILE.close();

  }
#endif

#ifdef SWAP_SAME_SIZE_BOXES
  // Phase 3 of loadbalance:  exhanges of identical sized boxes in order
  // to improve locality.  Need to handle periodicity somehow also, but
  // that will have to wait for loadbalance to be redesigned.

  bool sorted = true;
  for (int i=1; i<Nboxes; i++)
    {
      if (a_boxes[i] < a_boxes[i-1]) sorted = false;
    }
  Vector<std::list<int> > neighbours(Nboxes);

  for (int i=0; i<Nboxes; i++)
    {
      Box b = a_boxes[i];
      b.grow(1);
      std::list<int>& list = neighbours[i];
      for (int ii=0; ii<Nboxes; ii++)
        {
          if (ii == i); //don't bother with self intersect
          else
            {
              const Box& bb = a_boxes[ii];

              if (b.intersects(bb))
                {
                  list.push_front(ii);

                }
              if (sorted && bb.smallEnd(0) > b.bigEnd(0)) ii=Nboxes;
            }
        }
    }

  int swapcount = 6;
  int diffBoxes;
  int sameBoxes;
  int iter = 0;
  while (iter < 20 && swapcount > 0)
    {
      swapcount = 0;
      diffBoxes = 0;
      sameBoxes = 0;
      iter++;

      for (int i=0; i<neighbours.size(); i++)
        {

          std::list<int>& list = neighbours[i];
          for (int ii=i+1; ii<neighbours.size(); ii++)
            {
              int p = a_procAssignments[i];
              int pii = a_procAssignments[ii];
              if (a_boxes[i].numPts() != a_boxes[ii].numPts()) diffBoxes++; // don't swap
              else if (p == pii); //do nothing, no swaps on same proc
              else
                {
                  sameBoxes++;
                  std::list<int>& listii = neighbours[ii];

                  std::list<int>::iterator it;
                  // count : number of "local" exchanges before swap
                  // countafter: number of "local" exchanges after swap
                  int count = 0; int countafter = 0;
                  for (it = list.begin();it != list.end(); it++)
                    {
                      if (a_procAssignments[*it] == p)
                        count++;
                      if (a_procAssignments[*it] == pii)
                        countafter++;
                    }
                  for (it = listii.begin();it != listii.end(); it++)
                    {
                      if (a_procAssignments[*it] == pii)
                        count++;
                      if (a_procAssignments[*it] == p)
                        countafter++;
                    }
                  if (countafter > count)
                    {
                      swapcount++;
                      a_procAssignments[i] = pii;
                      a_procAssignments[ii] = p;
                      //pout()<<"count : "<<count<<" countafter :"<<countafter
                      //    <<std::endl;
                    }
                }
            }
        }
      //pout()<<"Boxes: "<<Nboxes<<" same: "<<sameBoxes
      //    <<" swaps: "<<swapcount<<std::endl;
    }
#endif
  return status;
}

///
// This version takes compute loads directly and does the "simple"
// load balance algorithm (modified knapsack algorithm).
// Boxes are also included, but only used if SWAP_SAME_SIZES_BOXES
// is enabled
// I think we can remove this as it simply calls the above LoadBalance() (ndk 8.4.2008)
///
int UnLongLongLoadBalance(Vector<int>&                      a_procAssignments,
                          const Vector<unsigned long long>& a_computeLoads,
                          const Vector<Box>&                a_boxes,
                          const int                         a_LBnumProc)
{
  CH_TIME("UnLongLongLoadBalance:VectorBoxSimple");

  Vector<long long> cloads(a_computeLoads.size());
  for (int i=0; i<a_computeLoads.size(); i++)
    {
      cloads[i] = (long long)a_computeLoads[i];
    }
  int status = LoadBalance(a_procAssignments, cloads, a_boxes, a_LBnumProc);
  return status;
}

///
// This version does the real work.
///
int
LoadBalance(Vector<int>&          a_procAssignments
            ,Real&                a_effRatio
            ,const Vector<Box>&   a_grids
            ,const Vector<long>&  a_computeLoads)
{
  CH_TIME("LoadBalance:VectorBoxWork");
  Vector<Vector<Box> >   grids(1, a_grids);
  Vector<Vector<long> >  computeLoads(1, a_computeLoads);
  Vector<int>            refRatios(1,2);
  Vector<Vector<int> > procs;
  int nproc = numProc();
  int retval = LoadBalance(procs, a_effRatio, grids, computeLoads,
                           refRatios, nproc);

  a_procAssignments = procs[0];
  return retval;
}

int
LoadBalance(Vector<Vector<int> >& procAssignments  //output: processor number
            ,Real&                 effRatio         //output: ratio of min load
            ,const Vector<Vector<Box> >&  Grids            //input: meshes to balance
            ,const Vector<Vector<long> >& ComputeLoads     //input: computational cost
            ,const Vector<int>&           RefRatios        //input: refinement ratio
            ,int nProc                            // number of procs to assugn to
            )
{
  CH_TIME("LoadBalance:VectorBoxRealWork");
  // local variables
  Real eff_ratio; // efficiency ratio on a level
  int status = 0; // return code

  // Validate inputs
  if ( Grids.size() != ComputeLoads.size() )
    { return -1011; }
  if ( Grids.size() != RefRatios.size() )
    { return -1013; }
  for ( int lvl=0; lvl<Grids.size(); ++lvl )
    {
      if ( Grids[lvl].size() != ComputeLoads[lvl].size() )
        { return -1012; }
    }

  // set the number of elements in the output vector to the number
  // of levels and the number of elements in each element to the
  // number of boxes on each level and set the value of each
  // element to zero
  procAssignments.resize( Grids.size() );
  for ( int lvl=0; lvl<Grids.size(); ++lvl )
    {
      procAssignments[lvl].resize( Grids[lvl].size(),0 );
    }

  // check for special case of all loads on 1 processor
  if ( nProc == 1 )
    {
      for ( int lvl=0; lvl<Grids.size(); ++lvl )
        {
          for ( int i=0; i<Grids[lvl].size(); ++i )
            {
              procAssignments[lvl][i] = 0;
            }
        }
      effRatio = 1.0;
      status = 0;
    }
  else
    {
      // general case: loads on more than one processor
      effRatio = 1.0;

      // balance each level separately
      for ( int lvl=0; lvl<Grids.size(); ++lvl )
        {
          // first, build the load structure and sort by compute_load
          Vector<Load> loads( Grids[lvl].size() );
          for ( int i=0; i<Grids[lvl].size(); ++i )
            {
              loads[i].load = ComputeLoads[lvl][i];
              loads[i].grid_index = i;
            }
          //          std::sort( loads.begin() ,loads.end() );
          loads.sort();
          // do the initial assignments by sequentially
          // `handing out' the loads from largest to smallest
          Vector<long> total_loads( nProc,0 ); //total load per processor
          Vector<Vector<Load> > proc_loads( nProc ); //loads per processor
          int iproc_minload = 0; //processor with lowest load
          // loads are sorted in increasing order, so work backwards through the vector
          for ( int i=loads.size()-1; i>=0; --i )
            {
              // put the next load on the processor with the lowest total load
              proc_loads[iproc_minload].push_back( loads[i] );
              total_loads[iproc_minload] += loads[i].load;

              // recompute which processor has the lowest load
              //[NOTE: this would be faster if the loads were sorted]
              iproc_minload = min_element( total_loads );
            }
          // compute average load per processor, truncated to int
          long avg_load = 0;
          for ( int i=0; i<total_loads.size(); ++i ) avg_load += total_loads[i];
          avg_load /= nProc;

          // optimize the assignments by swapping a load off the
          // processor with the max load onto another processor
          // such that the load balance is improved
          int iter_count = 0, swap_count = 0;
          int iproc_maxload;
          long max_change; //largest change in load balance
          int ibmax,jbmax,ipmax,jpmax;  //box and processor indices corresponding to max_change

          while ( 1 )
            {
              max_change = 0;

              // find the processor that has the largest deviation from perfect load balance
              min_max_elements( iproc_minload ,iproc_maxload ,total_loads );
              if ( iproc_minload == iproc_maxload )
                {
                  // load balance is perfect
                  // (this won't happen except in test cases)
                  break;
                }
              CH_assert( total_loads[iproc_minload] <= avg_load &&
                      avg_load <= total_loads[iproc_maxload] );
              if ( avg_load - total_loads[iproc_minload] > total_loads[iproc_maxload] - avg_load )
                ipmax = iproc_minload;
              else
                ipmax = iproc_maxload;

              //[NOTE: dont need this here, but it may be useful for debugging.]
              eff_ratio = (Real)total_loads[iproc_minload] / (Real)total_loads[iproc_maxload];

              // deviation from perfect load balance for this proc
              long devib = total_loads[ipmax] - avg_load;

              // search all the other processors for the swap that has the maximum
              // reduction in the total deviation from the perfect load balance
              for ( int j=0; j<proc_loads.size(); ++j )
                {
                  if ( j != ipmax )
                    {
                      long devjb = total_loads[j] - avg_load;

                      // loop over all boxes on both processors
                      for ( int ibox=0; ibox<proc_loads[ipmax].size(); ++ibox )
                        {
                          for ( int jbox=0; jbox<proc_loads[j].size(); ++jbox )
                            {
                              iter_count++;
                              // how much bigger is the ibox load than the jbox load?
                              long diff = proc_loads[ipmax][ibox].load
                                - proc_loads[  j  ][jbox].load;
                              // change in total deviation from swapping boxes
                              long change = Abs( devib ) + Abs( devjb )
                                - Abs( devib - diff ) - Abs( devjb + diff );
                              // remember this pair of boxes if the change is better
                              //[NOTE: max_change starts at 0, so this is always an improvement]
                              if ( change > max_change )
                                {
                                  max_change = change;
                                  ibmax = ibox; jbmax = jbox; jpmax = j;
                                }
                            }
                        }
                    }
                }
              // if there is a swap that improves load balance, take it; else stop
              if ( max_change > 0 )
                {
                  // adjust the total loads on each processor
                  long load_diff = proc_loads[ipmax][ibmax].load
                    - proc_loads[jpmax][jbmax].load;
                  CH_assert( load_diff != 0 );
                  total_loads[ipmax] -= load_diff;
                  total_loads[jpmax] += load_diff;
                  // swap the loads
                  Load tmp = proc_loads[ipmax][ibmax];
                  proc_loads[ipmax][ibmax] = proc_loads[jpmax][jbmax];
                  proc_loads[jpmax][jbmax] = tmp;

                  swap_count++;
                }
              else
                {
                  break;
                }
            }

          // Done with this level.

          // Compute the final efficiency ratio and save it if appropriate.
          min_max_elements( iproc_minload ,iproc_maxload ,total_loads );
          eff_ratio = (Real)total_loads[iproc_minload] / (Real)total_loads[iproc_maxload];
          if ( eff_ratio < effRatio ) effRatio = eff_ratio;

          // Assign boxes to processors for this level.
          for ( int ip=0; ip<proc_loads.size(); ++ip )
            {
              for ( int jb=0; jb<proc_loads[ip].size(); ++jb )
                {
                  procAssignments[lvl][proc_loads[ip][jb].grid_index] = ip;
                }
            }

#ifndef NDEBUG
          //           if ( iter_count > 0 )
          //             {
          //               cout << "    debug: LoadBalance: level " << lvl << " used "
          //                    << iter_count << " iterations and "
          //                    << swap_count << " swaps to get efficiency ratio "
          //                    << eff_ratio << std::endl;
          //             }
#endif
        }

      // Done with all levels.
      // We could try to permute processors assignments between levels to
      // reduce communication, but it is probably not worth the effort
      // since it probably would have O(N^4) cost (N==#boxes).
    }

  return status;
}

/// convenience function to gather a distributed set of Boxes with their corresponding processor assignment
/** Assumption is that each processor has at most one valid box. This is useful when interacting with other distributed codes which might not have the entire set of distributed boxes on all processors.
 */
int LoadBalance(Vector<int>& a_procAssignments, 
                Vector<Box>& a_boxes,
                const Box&   a_localGridBox,
                const int    a_numProc)
{
  int status = 0;

  const int numBox = numProc();

  Vector<Box> tempBoxes;
  Vector<int> tempProcAssign;
  tempBoxes.resize(numBox);

  if (numBox == 1) tempBoxes[0] = a_localGridBox;

#ifdef CH_MPI
  Box* nonConstBox = const_cast<Box*>(&a_localGridBox);  
  int boxSize = sizeof(Box);
  status = MPI_Allgather(nonConstBox,  boxSize,  MPI_BYTE, &(tempBoxes[0]), 
                         boxSize , MPI_BYTE , Chombo_MPI::comm);
#endif

  tempProcAssign.resize(tempBoxes.size());
  for (int i=0; i< tempBoxes.size(); i++)
    {
      tempProcAssign[i] = i;
    }

  // filter out any empty boxes (it's OK if there are fewer 
  // boxes than processors, but DisjointBoxLayout will 
  // choke if given an empty box, so remove them now...
  
  // this isn't the most efficient way to do this -- assumption
  // is that we're not dealing with all that many boxes, and that 
  // this is only done once anyway...
  for (int i=0; i<tempBoxes.size(); i++)
    {
      if (!tempBoxes[i].isEmpty() )
        {
          a_boxes.push_back(tempBoxes[i]);
          a_procAssignments.push_back(tempProcAssign[i]);
        }
    }
  
  return status;
}

////////////////////////////////////////////////////////////////
//                utility functions                           //
////////////////////////////////////////////////////////////////

//
// Find the index of the small value in a non-empty (long) Vector
//
int
min_element( const Vector<long>& Vect )
{
  CH_assert( Vect.size() > 0 );
  int imin = 0;
  for ( int i=1; i<Vect.size(); ++i )
    {
      if ( Vect[i] < Vect[imin] ) imin = i;
    }
  return imin;
}

//
// Find the indices of the smallest and largest values in a non-empty (long) Vector
//
void
min_max_elements( int& imin ,int& imax ,const Vector<long>& Vect )
{
  CH_assert( Vect.size() > 0 );
  imin = 0; imax = 0;
  for ( int i=1; i<Vect.size(); ++i )
    {
      if ( Vect[i] < Vect[imin] ) imin = i;
      if ( Vect[i] > Vect[imax] ) imax = i;
    }
  return;
}
void
min_max_elements( int& imin ,int& imax ,const Vector<long long>& Vect )
{
  CH_assert( Vect.size() > 0 );
  imin = 0; imax = 0;
  for ( int i=1; i<Vect.size(); ++i )
    {
      if ( Vect[i] < Vect[imin] ) imin = i;
      if ( Vect[i] > Vect[imax] ) imax = i;
    }
  return;
}

#include "NamespaceFooter.H"
