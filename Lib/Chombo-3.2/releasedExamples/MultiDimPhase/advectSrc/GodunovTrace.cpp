#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*****************/
/*****************/
#include "GodunovTrace.H"
#include "GodunovTraceF_F.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "UsingNamespace.H"

//#define SIMPLEUPWIND
#define SET_BOGUS_VALUES
#define BOGUS_VALUE 1.0e5

void TraceState(/// state at time t+dt/2 on edges in direction dir
                FArrayBox& a_stateHalf,
                /// cell-centered state at time t
                const FArrayBox& a_state,
                /// cell-centered velocity at time t
                const FArrayBox& a_cellVel,
                /// edge-centered advection velocity at time t+dt/2
                const FluxBox& a_advectionVel,
                /// cell-centered source
                const FArrayBox& a_source,
                /// Physical domain
                const ProblemDomain& a_dProblem,
                /// interior of grid patch
                const Box& a_gridBox,
                /// timeStep
                const Real a_dt,
                /// cell-spacing
                const Real a_dx,
                /// direction in which to perform tracing
                const int a_dir,
                /// which components to trace
                const Interval& a_srcComps,
                ///  where to put traced components in a_stateHalf,
                const Interval& a_destComps)

{

  int ncomp = a_srcComps.size();
  CH_assert (ncomp == a_destComps.size());
  int offset = a_destComps.begin() - a_srcComps.begin();
  Box edgeBox = a_gridBox;
  edgeBox.surroundingNodes(a_dir);

#ifdef SIMPLEUPWIND
  // simple cell-to-edge averaging may be causing us problems --
  // instead do simple upwinding

  for (int comp=a_srcComps.begin(); comp <= a_srcComps.end(); comp++)
  {
    int destcomp = comp+offset;
    FORT_UPWINDCELLTOEDGE(CHF_FRA1(a_stateHalf, destComp),
                          CHF_CONST_FRA1(a_state,comp),
                          CHF_CONST_FRA(a_advectionVel[a_dir]),
                          CHF_BOX(edgeBox),
                          CHF_CONST_INT(a_dir));
  }

#else


  // first compute slopes
  Box slopesBox = grow(a_gridBox,1);

  // these are debugging changes to sync with old code
#ifdef MATCH_OLDCODE
  FluxBox tempEdgeVel(slopesBox,1);
  tempEdgeVel.setVal(0.0);

  FArrayBox tempCellVel(slopesBox,SpaceDim);
  tempCellVel.setVal(0.0);


  CellToEdge(a_cellVel, tempEdgeVel);
  EdgeToCell(tempEdgeVel, tempCellVel);
  tempEdgeVel.clear(); // done with this, so reclaim memory
  EdgeToCell(a_advectionVel, tempCellVel);
#endif

  FArrayBox delS(slopesBox,1);
  FArrayBox sHat(slopesBox,1);
  FArrayBox sTilde(a_state.box(),1);

  for (int comp=a_srcComps.begin(); comp<=a_srcComps.end(); comp++)
  {
    int destComp = comp+offset;

#ifdef SET_BOGUS_VALUES
    delS.setVal(BOGUS_VALUE);
    sHat.setVal(BOGUS_VALUE);
    sTilde.setVal(BOGUS_VALUE);
#endif

    // compute van leer limited slopes in the normal direction
    // for now, always limit slopes, although might want to make
    // this a parameter
    int limitSlopes = 1;
    FORT_SLOPES(CHF_FRA1(delS,0),
                CHF_CONST_FRA1(a_state,comp),
                CHF_BOX(slopesBox),
                CHF_CONST_INT(a_dir),
                CHF_CONST_INT(limitSlopes));

    // this is a way to incorporate the minion correction --
    // stateTilde = state + (dt/2)*source
    sTilde.copy(a_source,comp,0,1);
    Real sourceFactor = a_dt/2.0;
    sTilde *= sourceFactor;
    sTilde.plus(a_state, comp,0,1);
    sHat.copy(sTilde,slopesBox);

    // now loop over directions, adding transverse components
    for (int localDir=0; localDir < SpaceDim; localDir++)
    {
      if (localDir != a_dir)
      {
        // add simple transverse components
        FORT_TRANSVERSE(CHF_FRA1(sHat,0),
                        CHF_CONST_FRA1(sTilde,0),
#ifdef MATCH_OLDCODE
                        CHF_CONST_FRA1(tempCellVel, localDir),
#else
                        CHF_CONST_FRA1(a_cellVel,localDir),
#endif
                        CHF_BOX(slopesBox),
                        CHF_CONST_REAL(a_dt),
                        CHF_CONST_REAL(a_dx),
                        CHF_CONST_INT(localDir));
                        //CHF_FRA1(temp,0));  // not used in function

      }
      else if (SpaceDim == 3)
      {
        // only add cross derivative transverse piece if we're in 3D
        // note that we need both components of velocity for this one
        FORT_TRANSVERSECROSS(CHF_FRA1(sHat,0),
                             CHF_CONST_FRA1(sTilde,0),
                             CHF_CONST_FRA(a_cellVel),
                             CHF_BOX(slopesBox),
                             CHF_CONST_REAL(a_dt),
                             CHF_CONST_REAL(a_dx),
                             CHF_CONST_INT(localDir));
      }
    } // end loop over directions

    // now compute left and right states and resolve the Reimann
    // problem to get single set of edge-centered values
    FORT_PREDICT(CHF_FRA1(a_stateHalf,destComp),
                 CHF_CONST_FRA1(sHat,0),
                 CHF_CONST_FRA1(delS,0),
                 CHF_CONST_FRA1(a_cellVel,a_dir),
                 CHF_CONST_FRA1(a_advectionVel[a_dir],0),
                 CHF_BOX(edgeBox),
                 CHF_CONST_REAL(a_dt),
                 CHF_CONST_REAL(a_dx),
                 CHF_CONST_INT(a_dir));
  } // end loop over components

#endif


}

// deprecated interface
void TraceState(/// state at time t+dt/2 on edges in direction dir
                FArrayBox& a_stateHalf,
                /// cell-centered state at time t
                const FArrayBox& a_state,
                /// cell-centered velocity at time t
                const FArrayBox& a_cellVel,
                /// edge-centered advection velocity at time t+dt/2
                const FluxBox& a_advectionVel,
                /// cell-centered source
                const FArrayBox& a_source,
                /// Physical domain
                const Box& a_dProblem,
                /// interior of grid patch
                const Box& a_gridBox,
                /// timeStep
                const Real a_dt,
                /// cell-spacing
                const Real a_dx,
                /// direction in which to perform tracing
                const int a_dir,
                /// which components to trace
                const Interval& a_srcComps,
                ///  where to put traced components in a_stateHalf,
                const Interval& a_destComps)
{
  ProblemDomain physdomain(a_dProblem);

  TraceState(a_stateHalf, a_state, a_cellVel, a_advectionVel,
             a_source, physdomain, a_gridBox, a_dt, a_dx, a_dir,
             a_srcComps, a_destComps);
}

void TraceState(/// state at time t+dt/2 on edges in direction dir
                FArrayBox& a_stateHalf,
                /// cell-centered state at time t
                const FArrayBox& a_state,
                /// cell-centered velocity at time t
                const FArrayBox& a_cellVel,
                /// cell-centered source
                const FArrayBox& a_source,
                /// Physical domain
                const ProblemDomain& a_dProblem,
                /// interior of grid patch
                const Box& a_gridBox,
                /// timeStep
                const Real a_dt,
                /// cell-spacing
                const Real a_dx,
                /// direction in which to perform tracing
                const int a_dir,
                /// which components to trace
                const Interval& a_srcComps,
                ///  where to put traced components in a_stateHalf,
                const Interval& a_destComps)
{

  int ncomp = a_srcComps.size();
  CH_assert (ncomp == a_destComps.size());
  int offset = a_destComps.begin() - a_srcComps.begin();

  // think we need to do better than simple averaging.
  // first average cell-centered vel to edges
  Box edgeBox(a_gridBox);
  edgeBox.surroundingNodes(a_dir);
  FArrayBox edgeVel(edgeBox,1);
#ifdef SET_BOGUS_VALUES
  edgeVel.setVal(BOGUS_VALUE);
#endif

  CellToEdge(a_cellVel, a_dir, edgeVel, 0,a_dir);


#ifdef SIMPLEUPWIND
  // instead do simple upwinding

  for (int comp=a_srcComps.begin(); comp <= a_srcComps.end(); comp++)
  {
    int destComp = comp+offset;
    FORT_UPWINDCELLTOEDGE(CHF_FRA1(a_stateHalf, destcomp),
                          CHF_CONST_FRA1(a_state,comp),
                          CHF_CONST_FRA(edgeVel),
                          CHF_BOX(edgeBox),
                          CHF_CONST_INT(a_dir));
  }

#else

  // first compute slopes
  Box slopesBox = grow(a_gridBox,1);
  FArrayBox delS(slopesBox,1);
  FArrayBox sHat(slopesBox,1);
  FArrayBox sTilde(a_state.box(),1);

#ifdef MATCH_OLDCODE
  FArrayBox tempCellVel(slopesBox,SpaceDim);
  FluxBox tempEdgeVel(slopesBox,1);

  // these are debugging tools to sync with old code
  CellToEdge(a_cellVel, tempEdgeVel);
  EdgeToCell(tempEdgeVel, tempCellVel);
  tempEdgeVel.clear(); // reclaim memory
#endif

  for (int comp=a_srcComps.begin(); comp<=a_srcComps.end(); comp++)
  {
    int destComp = comp+offset;

#ifdef SET_BOGUS_VALUES
    delS.setVal(BOGUS_VALUE);
    sHat.setVal(BOGUS_VALUE);
    sTilde.setVal(BOGUS_VALUE);
#endif

    // compute van leer limited slopes in the normal direction
    // for now, always limit slopes, although might want to make
    // this a parameter
    int limitSlopes = 1;
    FORT_SLOPES(CHF_FRA1(delS,0),
                CHF_CONST_FRA1(a_state,comp),
                CHF_BOX(slopesBox),
                CHF_CONST_INT(a_dir),
                CHF_CONST_INT(limitSlopes));

    // this is a way to incorporate the minion correction --
    // stateTilde = state + (dt/2)*source
    sTilde.copy(a_source,comp,0,1);
    Real sourceFactor = a_dt/2.0;
    sTilde *= sourceFactor;
    sTilde.plus(a_state, comp,0,1);
    sHat.copy(sTilde,slopesBox);


    // now loop over directions, adding transverse components
    for (int localDir=0; localDir < SpaceDim; localDir++)
    {
      if (localDir != a_dir)
      {
        // add simple transverse components
        FORT_TRANSVERSE(CHF_FRA1(sHat,0),
                        CHF_CONST_FRA1(sTilde,0),
#ifdef MATCH_OLDCODE
                        CHF_CONST_FRA1(tempCellVel,localDir),
#else
                        CHF_CONST_FRA1(a_cellVel,localDir),
#endif
                        CHF_BOX(slopesBox),
                        CHF_CONST_REAL(a_dt),
                        CHF_CONST_REAL(a_dx),
                        CHF_CONST_INT(localDir));
                        //CHF_FRA1(temp,0)); // not used
      }
      else if (SpaceDim == 3)
      {
        // only add cross derivative transverse piece if we're in 3D
        // note that we need both components of velocity for this one
        FORT_TRANSVERSECROSS(CHF_FRA1(sHat,0),
                             CHF_CONST_FRA1(sTilde,0),
                             CHF_CONST_FRA(a_cellVel),
                             CHF_BOX(slopesBox),
                             CHF_CONST_REAL(a_dt),
                             CHF_CONST_REAL(a_dx),
                             CHF_CONST_INT(localDir));
      }
    } // end loop over directions

    // now compute left and right states and resolve the Reimann
    // problem to get single set of edge-centered values
    FORT_PREDICT(CHF_FRA1(a_stateHalf,destComp),
                 CHF_CONST_FRA1(sHat,0),
                 CHF_CONST_FRA1(delS,0),
                 CHF_CONST_FRA1(a_cellVel,a_dir),
                 CHF_CONST_FRA1(edgeVel,0),
                 CHF_BOX(edgeBox),
                 CHF_CONST_REAL(a_dt),
                 CHF_CONST_REAL(a_dx),
                 CHF_CONST_INT(a_dir));
  } // end loop over components

#endif


}

// deprecated interface
void TraceState(/// state at time t+dt/2 on edges in direction dir
                FArrayBox& a_stateHalf,
                /// cell-centered state at time t
                const FArrayBox& a_state,
                /// cell-centered velocity at time t
                const FArrayBox& a_cellVel,
                /// cell-centered source
                const FArrayBox& a_source,
                /// Physical domain
                const Box& a_dProblem,
                /// interior of grid patch
                const Box& a_gridBox,
                /// timeStep
                const Real a_dt,
                /// cell-spacing
                const Real a_dx,
                /// direction in which to perform tracing
                const int a_dir,
                /// which components to trace
                const Interval& a_srcComps,
                ///  where to put traced components in a_stateHalf,
                const Interval& a_destComps)
{
  ProblemDomain physdomain(a_dProblem);

  TraceState(a_stateHalf, a_state, a_cellVel, a_source, physdomain,
             a_gridBox, a_dt, a_dx, a_dir, a_srcComps, a_destComps);
}










