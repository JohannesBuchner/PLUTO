#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELRK4_H__
#define _LEVELRK4_H__

/// Basic templated implementation of RK4 advance for a single AMRLevel
/**
    This implements fourth-order Runge-Kutta time advance of an ODE.
    ODE is d(Soln)/dt = RHS

    Template types:
    \li TSoln is the datatype for the solution
    \li TFR is a flux-register datatype
    \li TOp is an object encapsulating the actual ODE operations

    TOp requirements: TOp must have the following functions:
\code
        // evaluate d(soln)/dt at current time based on soln
        void evalRHS(TSoln& rhs, // d(soln)/dt based on soln
                     TSoln& soln, // soln at current time
                     TFR& fineFR,  // flux register w/ finer level
                     TFR& crseFR,  // flux register w/ crse level
                     const TSoln& oldCrseSoln, // old-time crse solution
                     Real oldCrseTime,    // old crse time
                     const TSoln& newCrseSoln,  // new-time crse solution
                     Real newCrseTime,   // new crse time
                     Real time,   // current time centering of soln
                     Real fluxWeight // weight to apply to fluxRegister updates
                    )

        // implements soln += dt*rhs
        void updateODE(TSoln& soln,
                       const TSoln& rhs,
                       Real dt)

        // define data holder newSoln based on existingSoln,
        // including ghost cell specification
        void defineSolnData(TSoln& newSoln,
                            const TSoln& existingSoln)

        // define data holder for RHS based on existingSoln
        // including ghost cell specification
        // (which in most cases is no ghost cells)
        void defineRHSData(TSoln& newRHS, const TSoln& existingSoln)

        /// copy data in src into dest
        void copySolnData(TSoln& dest, const TSoln& src)
\endcode
*/

template <class TSoln, class TFR, class TOp>
void RK4LevelAdvance(/// the cell-centered solution at the new time (a_time + a_dt)
                     TSoln& a_newSoln,
                     /// the cell-centered solution at the old time (a_time)
                     TSoln& a_oldSoln,
                     /// old-time coarse-level solution (if it exists)
                     const TSoln& a_oldCrseSoln,
                     /// time-centering of a_oldCrseSoln
                     Real a_oldCrseTime,
                     /// new-time coarse-level solution (if it exists)
                     const TSoln& a_newCrseSoln,
                     /// time-centering of a_newCrseSoln
                     Real a_newCrseTime,
                     /// the pointer to the flux register between this level and the next coarser level
                     TFR& a_crseFRPtr,
                     /// the pointer to the flux register between this level and the next finer level
                     TFR& a_fineFRPtr,
                     /// time centering of a_oldSoln
                     Real a_time,
                     /// time step
                     Real a_dt,
                     /// object which encapsulates the ODE functionality (must meet the requirements detailed above)
                     TOp& a_op,
                     /// whether a_newSoln is initialized with a_oldSoln (default true)
                     const bool a_initializeNewSoln = true)
{
  // Advance by one time step using 4th-order Runge-Kutta.

  // allocate and define temp storage
  TSoln RHSTmp;
  TSoln solnTmp;
  a_op.defineRHSData(RHSTmp,a_oldSoln);
  a_op.defineSolnData(solnTmp,a_oldSoln);
  a_op.copySolnData(solnTmp,a_oldSoln);

  // start by copying old solution into new solution
  if (a_initializeNewSoln)
    {
      a_op.copySolnData(a_newSoln, a_oldSoln);
    }

  // RK Stage 1: compute RHS.at old time
  // RHSTmp gets RHS(tOld)
  // flux registers are also incremented if appropriate
  a_op.evalRHS(RHSTmp,a_oldSoln,
               a_fineFRPtr,a_crseFRPtr,
               a_oldCrseSoln, a_oldCrseTime,
               a_newCrseSoln, a_newCrseTime,
               a_time,a_dt/6);

  // RK Stage 1: compute update {phi}_1, partial update to new vals.

  // first part of update...
  a_op.updateODE(a_newSoln,RHSTmp,a_dt/6);


  // predict half-time value using LofPhiTmp
  a_op.updateODE(solnTmp, RHSTmp,a_dt/2);

  // RK Stage 2: compute RHS based on first predicted value
  // also increment flux registers
  a_op.evalRHS(RHSTmp,solnTmp,
               a_fineFRPtr,a_crseFRPtr,
               a_oldCrseSoln, a_oldCrseTime,
               a_newCrseSoln, a_newCrseTime,
               a_time+a_dt/2,a_dt/3);


  // RK Stage 2: compute update {phi}_2, partial update to new vals.

  a_op.updateODE(a_newSoln,RHSTmp,a_dt/3);

  // now predict half-time value using latest estimate of RHS
  a_op.copySolnData(solnTmp,a_oldSoln);
  a_op.updateODE(solnTmp,RHSTmp,a_dt/2);

  // RK Stage 3: compute new RHS.

  a_op.evalRHS(RHSTmp,solnTmp,
               a_fineFRPtr,a_crseFRPtr,
               a_oldCrseSoln, a_oldCrseTime,
               a_newCrseSoln, a_newCrseTime,
               a_time+a_dt/2,a_dt/3);


  // RK Stage 3: compute update {phi}_3, partial update to new vals.

  a_op.updateODE(a_newSoln,RHSTmp,a_dt/3);

  // predict value at new time
  a_op.copySolnData(solnTmp,a_oldSoln);
  a_op.updateODE(solnTmp,RHSTmp,a_dt);

  // RK Stage 4: compute RHS.

  a_op.evalRHS(RHSTmp,solnTmp,
               a_fineFRPtr,a_crseFRPtr,
               a_oldCrseSoln, a_oldCrseTime,
               a_newCrseSoln, a_newCrseTime,
               a_time+a_dt,a_dt/6);

  // RK Stage 4: compute final update of solution.

  a_op.updateODE(a_newSoln,RHSTmp,a_dt/6);
}


/// Templated implementation of RK4 advance for a single AMRLevel, allowing interpolation in time
/**
    This implements fourth-order Runge-Kutta time advance of an ODE,
    and saves intermediate results for use in interpolation in time.
    ODE is d(Soln)/dt = RHS

    Template types:
    \li TSoln is the datatype for the solution
    \li TInterp is a storage class for intermediate values
    \li TFR is a flux-register datatype
    \li TOp is an object encapsulating the actual ODE operations

    TInterp requirements: TInterp must have the following functions:
\code
        // set time step
        void setDt(Real dt)

        // save initial solution
        void saveInitialSoln(TSoln& soln)

        // save this value of d(soln)/dt
        void saveRHS(TSoln& rhs)
\endcode
    TOp requirements: TOp must have the following functions:
\code
        // evaluate d(soln)/dt at current time based on soln
        void evalRHS(TSoln& rhs, // d(soln)/dt based on soln
                     TSoln& soln, // soln at current time
                     TFR& fineFR,  // flux register w/ finer level
                     TFR& crseFR,  // flux register w/ crse level
                     const TSoln& oldCrseSoln, // old-time crse solution
                     Real oldCrseTime,    // old crse time
                     const TSoln& newCrseSoln,  // new-time crse solution
                     Real newCrseTime,   // new crse time
                     Real time,   // current time centering of soln
                     Real fluxWeight // weight to apply to fluxRegister updates
                    )

        // implements soln += dt*rhs
        void updateODE(TSoln& soln,
                       const TSoln& rhs,
                       Real dt)

        // define data holder newSoln based on existingSoln,
        // including ghost cell specification
        void defineSolnData(TSoln& newSoln,
                            const TSoln& existingSoln)

        // define data holder for RHS based on existingSoln
        // including ghost cell specification
        // (which in most cases is no ghost cells)
        void defineRHSData(TSoln& newRHS, const TSoln& existingSoln)

        /// copy data in src into dest
        void copySolnData(TSoln& dest, const TSoln& src)
\endcode
*/

template <class TSoln, class TInterp, class TFR, class TOp>
void RK4LevelAdvance(/// the cell-centered solution at the new time (a_time + a_dt)
                     TSoln& a_newSoln,
                     /// the cell-centered solution at the old time (a_time)
                     TSoln& a_oldSoln,
                     /// object that encapsulates time interpolation
                     TInterp& a_interpPtr,
                     /// old-time coarse-level solution (if it exists)
                     const TSoln& a_oldCrseSoln,
                     /// time-centering of a_oldCrseSoln
                     Real a_oldCrseTime,
                     /// new-time coarse-level solution (if it exists)
                     const TSoln& a_newCrseSoln,
                     /// time-centering of a_newCrseSoln
                     Real a_newCrseTime,
                     /// the pointer to the flux register between this level and the next coarser level
                     TFR& a_crseFRPtr,
                     /// the pointer to the flux register between this level and the next finer level
                     TFR& a_fineFRPtr,
                     /// time centering of a_oldSoln
                     Real a_time,
                     /// time step
                     Real a_dt,
                     /// object which encapsulates the ODE functionality (must meet the requirements detailed above)
                     TOp& a_op,
                     /// whether a_newSoln is initialized with a_oldSoln (default true)
                     const bool a_initializeNewSoln = true)


{

  // Advance by one time step using 4th-order Runge-Kutta.

  // allocate and define temp storage
  TSoln RHSTmp;
  TSoln solnTmp;
  a_op.defineRHSData(RHSTmp,a_oldSoln);
  a_op.defineSolnData(solnTmp,a_oldSoln);
  a_op.copySolnData(solnTmp,a_oldSoln);

  // start by copying old solution into new solution
  if (a_initializeNewSoln)
    {
      a_op.copySolnData(a_newSoln, a_oldSoln);
    }
  a_interpPtr.setDt(a_dt);
  a_interpPtr.saveInitialSoln(a_oldSoln);

  // RK Stage 1: compute RHS.at old time
  // RHSTmp gets RHS(tOld)
  // flux registers are also incremented if appropriate
  a_op.evalRHS(RHSTmp,a_oldSoln,
               a_fineFRPtr,a_crseFRPtr,
               a_oldCrseSoln, a_oldCrseTime,
               a_newCrseSoln, a_newCrseTime,
               a_time,a_dt/6);
  a_interpPtr.saveRHS(RHSTmp);

  // RK Stage 1: compute update {phi}_1, partial update to new vals.

  // first part of update...
  a_op.updateODE(a_newSoln,RHSTmp,a_dt/6);


  // predict half-time value using LofPhiTmp
  a_op.updateODE(solnTmp, RHSTmp,a_dt/2);

  // RK Stage 2: compute RHS based on first predicted value
  // also increment flux registers
  a_op.evalRHS(RHSTmp,solnTmp,
               a_fineFRPtr,a_crseFRPtr,
               a_oldCrseSoln, a_oldCrseTime,
               a_newCrseSoln, a_newCrseTime,
               a_time+a_dt/2,a_dt/3);
  a_interpPtr.saveRHS(RHSTmp);


  // RK Stage 2: compute update {phi}_2, partial update to new vals.
  a_op.updateODE(a_newSoln,RHSTmp,a_dt/3);

  // now predict half-time value using latest estimate of RHS
  a_op.copySolnData(solnTmp,a_oldSoln);
  a_op.updateODE(solnTmp,RHSTmp,a_dt/2);

  // RK Stage 3: compute new RHS.

  a_op.evalRHS(RHSTmp,solnTmp,
               a_fineFRPtr,a_crseFRPtr,
               a_oldCrseSoln, a_oldCrseTime,
               a_newCrseSoln, a_newCrseTime,
               a_time+a_dt/2,a_dt/3);
  a_interpPtr.saveRHS(RHSTmp);


  // RK Stage 3: compute update {phi}_3, partial update to new vals.

  a_op.updateODE(a_newSoln,RHSTmp,a_dt/3);

  // predict value at new time
  a_op.copySolnData(solnTmp,a_oldSoln);
  a_op.updateODE(solnTmp,RHSTmp,a_dt);

  // RK Stage 4: compute RHS.

  a_op.evalRHS(RHSTmp,solnTmp,
               a_fineFRPtr,a_crseFRPtr,
               a_oldCrseSoln, a_oldCrseTime,
               a_newCrseSoln, a_newCrseTime,
               a_time+a_dt,a_dt/6);
  a_interpPtr.saveRHS(RHSTmp);

  // RK Stage 4: compute final update of solution.

  a_op.updateODE(a_newSoln,RHSTmp,a_dt/6);
}

#endif // end multiple-include preventer
