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

#include "AMRNavierStokes.H"
#include "AMRIO.H"
#include "Divergence.H"
#include "EdgeToCell.H"
#include "SetValLevel.H"

#include "AMRNSF_F.H"
#include "probF_F.H"

#ifdef CH_USE_HDF5
#include "CH_HDF5.H"
#endif

#ifdef CH_USE_HDF5
// ---------------------------------------------------------------
void
AMRNavierStokes::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::writeCheckpointHeader" << endl;
    }

  HDF5HeaderData header;
  // since number of velocity components is obviously SpaceDim,
  // only need to do number of scalar components here
  CH_assert (s_num_vel_comps == SpaceDim);
  header.m_int["num_components"] = s_num_scal_comps;
  char comp_str[30];
  for (int comp=0; comp<s_num_scal_comps; ++comp)
    {
      sprintf (comp_str, "component_%d", comp);
      header.m_string[comp_str] = s_scal_names[comp];
    }

  // now write lambda name
  header.m_string["lambda_component"] = "lambda";

  // now write velocity names..
  for (int comp=0; comp<s_num_vel_comps; ++comp)
    {
      sprintf (comp_str, "vel_component_%d", comp);
      header.m_string[comp_str] = s_vel_names[comp];
    }

  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout () << header << endl;
    }
}

// ---------------------------------------------------------------
void
AMRNavierStokes::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::writeCheckpointLevel" << endl;
    }

  char level_str[20];
  sprintf (level_str, "%d", m_level);
  const std::string label = std::string("level_") + level_str;

  a_handle.setGroup(label);

  HDF5HeaderData header;

  header.m_int ["ref_ratio"] = m_ref_ratio;
  header.m_real["dx"]        = m_dx;
  header.m_real["dt"]        = m_dt;
  header.m_real["time"]      = m_time;
  header.m_real["cfl"]       = m_cfl;

  header.m_int["finest_level"] = m_finest_level;
  header.m_int["is_empty"] = m_is_empty;

  // don't write out static variables -- get them from inputs files?
  // or should we checkpoint those?

  header.m_box["prob_domain"] = m_problem_domain.domainBox();

  // write out periodicity info
  D_TERM(
         if (m_problem_domain.isPeriodic(0))
         {
           header.m_int ["is_periodic_0"] = 1;
         }
         else
         {
           header.m_int ["is_periodic_0"] = 0;
         }
  ,
         if (m_problem_domain.isPeriodic(1))
         {
           header.m_int ["is_periodic_1"] = 1;
         }
         else
         {
           header.m_int ["is_periodic_1"] = 0;
         }
  ,
         if (m_problem_domain.isPeriodic(2))
         {
           header.m_int ["is_periodic_2"] = 1;
         }
         else
         {
           header.m_int ["is_periodic_2"] = 0;
         }
  );

  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout () << header << endl;
    }

  if (!isEmpty())
    {
      // now dump out this level's data
      write (a_handle, m_vel_new_ptr->boxLayout());
      write (a_handle, *m_vel_new_ptr, "new_velocity");
      write (a_handle, *m_lambda_new_ptr, "new_lambda");

      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          char scal_str[20];
          sprintf(scal_str,"scalar_component_%d",comp);
          write (a_handle, *m_scal_new[comp], scal_str);
        }

      // now do level checkpoint for projection
      m_projection.writeCheckpointLevel(a_handle);
    }
}

// ---------------------------------------------------------------
void
AMRNavierStokes::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::readCheckpointHeader" << endl;
    }

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << "hdf5 header data:" << endl;
      pout () << header << endl;
    }

  // read number of components
  if (header.m_int.find("num_components") == header.m_int.end())
    {
      MayDay::Error("AMRNavierStokes::readCheckpointHeader: checkfile does not have num_components");
    }
  int num_comps = header.m_int["num_components"];
  if (num_comps != s_num_scal_comps)
    {
      MayDay::Error("AMRNavierStokes::readCheckpointHeader: num_components in input file does not match solver");
    }

  // read component names
  std::string state_name;
  char comp_str[60];
  for (int comp=0; comp < s_num_scal_comps; ++comp)
    {
      sprintf (comp_str, "component_%d", comp);
      if (header.m_string.find(comp_str) == header.m_string.end())
        {
          MayDay::Error("AMRNavierStokes::readCheckpointHeader: checkfile does not have enough component names");
        }
      state_name = header.m_string[comp_str];
      if (state_name != s_scal_names[comp])
        {
          MayDay::Error("AMRNavierStokes::readCheckpointHeader: state_name in checkpoint does not match solver");
        }
    }

  sprintf(comp_str, "lambda_component");
  if (header.m_string.find(comp_str) == header.m_string.end())
    {
      MayDay::Error("AMRNavierStokes::readCheckpointHeader: checkfile does not have lambda name");
    }
  state_name = header.m_string[comp_str];
  if (state_name != "lambda")
    {
      MayDay::Error("AMRNavierStokes::readCheckpointHeader: lambda_name in checkfile does not match that in solver");
    }

  for (int comp=0; comp < SpaceDim; ++comp)
    {
      sprintf(comp_str, "vel_component_%d", comp);
      if (header.m_string.find(comp_str) == header.m_string.end())
        {
          MayDay::Error("AMRNavierStokes::readCheckpointHeader: checkfile does not have enough velocity names");
        }
      state_name = header.m_string[comp_str];
      if (state_name != s_vel_names[comp])
        {
          MayDay::Error("AMRNavierStokes::readCheckpointHEader: vel_name in checkfile does not match solver");
        }
    }
}

// ---------------------------------------------------------------
// read level data from checkpoint file
void
AMRNavierStokes::readCheckpointLevel (HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::readCheckpointLevel " << m_level << endl;
    }

  char level_str[20];
  sprintf(level_str, "%d", m_level);
  const std::string label = std::string("level_") + level_str;

  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout () << "hdf5 header data:" << endl;
      pout () << header << endl;
    }

  // read refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
      MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];

  // read dx
  if (header.m_real.find("dx") == header.m_real.end())
    {
      MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain dx");
    }
  m_dx = header.m_real["dx"];

  // read dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  // read time
  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain time");
    }
  m_time = header.m_real["time"];

  // read cfl
  if (header.m_real.find("cfl") == header.m_real.end())
    {
      MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain cfl");
    }
  m_cfl = header.m_real["cfl"];

  // read finest_level
  if (header.m_int.find("finest_level") == header.m_int.end())
    {
      MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain finest_level");
    }
  m_finest_level = header.m_int["finest_level"];

  // read is_empty
  if (header.m_int.find("is_empty") == header.m_int.end())
    {
      MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain is_empty");
    }
  m_is_empty = header.m_int["is_empty"];

  // read problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
    {
      MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain prob_domain");
    }
  Box domainBox = header.m_box["prob_domain"];

  // read in periodicity info -- this is more complicated than
  // it really needs to be in order to preserve backward compatibility
  bool is_periodic[SpaceDim];

  D_TERM(
         if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
         {
           is_periodic[0] =  (header.m_int["is_periodic_0"] == 1);
         }
         else
         {
           is_periodic[0] = false;
         }
  ,
         if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
         {
           is_periodic[1] =  (header.m_int["is_periodic_1"] == 1);
         }
         else
         {
           is_periodic[1] = false;
         }
  ,
         if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
         {
           is_periodic[2] =  (header.m_int["is_periodic_2"] == 1);
         }
         else
         {
           is_periodic[2] = false;
         }
  );

  m_problem_domain = ProblemDomain(domainBox, is_periodic);

  // read grids
  Vector<Box> grids;
  const int grid_status = read(a_handle, grids);
  if (grid_status != 0)
    {
      MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain a Vector<Box>");
    }
  // create level domain
  const DisjointBoxLayout level_domain = loadBalance(grids);

  LayoutIterator lit = level_domain.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = level_domain[lit];
      m_level_grids.push_back(b);
    }

  if (s_verbosity >= 4)
    {
      pout () << "read level domain: " << endl;
      LayoutIterator lit = level_domain.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box& b = level_domain[lit];
          pout() << lit().intCode() << ": " << b << endl;
        }
      pout () << endl;
    }

  if (!isEmpty())
    {
      // no need to do this if no data present on this level

      // reshape data holders with new grids
      // reshape state with new grids
      IntVect ghostVect(D_DECL(1,1,1));
      m_vel_new_ptr = new LevelData<FArrayBox>(level_domain,
                                               s_num_vel_comps, ghostVect);
      m_vel_old_ptr = new LevelData<FArrayBox>(level_domain,
                                               s_num_vel_comps, ghostVect);

      m_lambda_new_ptr = new LevelData<FArrayBox>(level_domain,
                                                  1, ghostVect);
      m_lambda_old_ptr = new LevelData<FArrayBox>(level_domain,
                                                  1, ghostVect);

#ifdef DEBUG
      m_vel_save_ptr = new LevelData<FArrayBox>(level_domain,
                                                s_num_vel_comps, ghostVect);

#endif
      m_scal_new.resize(s_num_scal_comps,NULL);
      m_scal_old.resize(s_num_scal_comps,NULL);
      // also resize flux registers here
      m_scal_fluxreg_ptrs.resize(s_num_scal_comps,NULL);
#if DEBUG
      m_scal_save.resize(s_num_scal_comps,NULL);
#endif

      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          m_scal_new[comp] = new LevelData<FArrayBox>(level_domain,
                                                      1, ghostVect);

          m_scal_old[comp] = new LevelData<FArrayBox>(level_domain,
                                                      1, ghostVect);

#if DEBUG
          m_scal_save[comp] = new LevelData<FArrayBox>(level_domain,
                                                       1, ghostVect);
#endif

          // this should get defined later
          m_scal_fluxreg_ptrs[comp] = new LevelFluxRegister;
        }

      setAllBogus();

      // read data
      LevelData<FArrayBox>& new_vel = *m_vel_new_ptr;

      const int velData_status = read<FArrayBox> (a_handle,
                                                  new_vel,
                                                  "new_velocity",
                                                  level_domain);

      if (velData_status != 0)
        {
          MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain new_velocity data");
        }

      LevelData<FArrayBox>& new_lambda = *m_lambda_new_ptr;

      const int lambdaData_status = read<FArrayBox> (a_handle,
                                                     new_lambda,
                                                     "new_lambda",
                                                     level_domain);

      if (lambdaData_status != 0)
        {
          MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain new_lambda data");
        }

      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          LevelData<FArrayBox>& new_scal = *m_scal_new[comp];

          char scal_str[20];
          sprintf(scal_str,"scalar_component_%d",comp);
          const int scalData_status = read<FArrayBox> (a_handle,
                                                       new_scal,
                                                       scal_str,
                                                       level_domain);

          if (scalData_status != 0)
            {
              MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain new_scalar data");
            }
        } // end loop over scalar data

    } // end if level is not empty
  else
    {
      m_vel_new_ptr = NULL;
      m_vel_old_ptr = NULL;
      m_lambda_new_ptr = NULL;
      m_lambda_old_ptr = NULL;
    } // end if level is empty

  levelSetup(level_domain);

  m_projection.readCheckpointLevel(a_handle);

  // finally, set physical and C/F boundary conditions on data
  for (int comp=0; comp<s_num_scal_comps; comp++)
    {
      LevelData<FArrayBox>& thisNewScal = newScal(comp);
      LevelData<FArrayBox>& thisOldScal = oldScal(comp);

      if (m_level > 0)
        {
          AMRNavierStokes& crseLevel = *crseNSPtr();
          LevelData<FArrayBox>& oldCrseScal = crseLevel.oldScal(comp);
          LevelData<FArrayBox>& newCrseScal = crseLevel.newScal(comp);
          const DisjointBoxLayout& crseGrids = oldCrseScal.getBoxes();
          const ProblemDomain& crseDomain = crseLevel.problemDomain();
          int nRefCrse = crseLevel.refRatio();

          // interpolate coarse-level data into ghost cells
          int scalGrow = thisNewScal.ghostVect()[0];
          PiecewiseLinearFillPatch filpatcher(level_domain, crseGrids,
                                              thisNewScal.nComp(),
                                              crseDomain, nRefCrse,
                                              scalGrow);
          // do old time, then new time
          Real timeInterpCoeff = 0;
          filpatcher.fillInterp(thisOldScal, oldCrseScal, newCrseScal,
                                timeInterpCoeff, 0,0,thisNewScal.nComp());

          timeInterpCoeff = 1.0;
          filpatcher.fillInterp(thisNewScal, oldCrseScal, newCrseScal,
                                timeInterpCoeff, 0,0,thisNewScal.nComp());
        }

      // do exchange
      Interval scalComps = thisNewScal.interval();
      thisNewScal.exchange(scalComps);

      BCHolder thisBC = m_physBCPtr->scalarTraceFuncBC(comp);
      // finally, set physical BC's
      DataIterator dit = thisNewScal.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          const Box& bx = level_domain[dit];
          thisBC(thisNewScal[dit], bx,
                 m_problem_domain, m_dx,
                 false); // inhomogeneous
          thisBC(thisOldScal[dit], bx,
                 m_problem_domain, m_dx,
                 false); // inhomogeneous
        }
    }
}

// ---------------------------------------------------------------
void
AMRNavierStokes::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::writePlotHeader " << m_level << endl;
    }

  HDF5HeaderData header;

  int numcomp = numPlotComps();

  header.m_int ["num_components"] = numcomp;
  char comp_str[30];
  int comp = 0;

  sprintf(comp_str, "component_%d", comp);
  header.m_string[comp_str] = "xVel";
  comp++;

  if (SpaceDim >1)
    {
      sprintf(comp_str, "component_%d", comp);
      header.m_string[comp_str] = "yVel";
      comp++;

      if (SpaceDim >2)
        {
          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "zVel";
          comp++;
        }
    }

  if (s_write_divergence)
    {
      sprintf(comp_str, "component_%d", comp);
      header.m_string[comp_str] = "divergence";
      comp++;
    }

  if (s_write_lambda)
    {
      sprintf(comp_str, "component_%d", comp);
      header.m_string[comp_str] = "lambda";
      comp++;
    }

  if (s_write_time_derivatives)
    {
      sprintf(comp_str, "component_%d", comp);
      header.m_string[comp_str] = "duDt";
      comp++;

      if (SpaceDim >1)
        {
          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "dvDt";
          comp++;

        if (SpaceDim >2)
          {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "dwDt";
            comp++;
          }
        }

      // dLambda/dt
      sprintf(comp_str, "component_%d", comp);
      header.m_string[comp_str] = "dLambda/dt";
      comp++;
    } // end if writing time derivatives

  if (s_write_vorticity)
    {
      if (SpaceDim == 2)
        {
          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "vorticity";
          comp++;
        }
      else if (SpaceDim == 3)
        {
          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "x_vort";
          comp++;

          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "y_vort";
          comp++;

          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "z_vort";
          comp++;

          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "mag_vort";
          comp++;
        } // end if 3d
    } // end if writing vorticity

  if (s_write_scalars)
    {
      for (int scal=0; scal<s_num_scal_comps; scal++)
        {
          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = s_scal_names[scal];
          comp++;
        }
    }

  // dS/dt
  if (s_write_dScalar_dt)
    {
      char scal_str[30];
      for (int scal=0; scal<s_num_scal_comps; scal++)
        {
          sprintf(comp_str, "component_%d", comp);
          sprintf(scal_str, "d(%s)/dt",s_scal_names[scal]);
          header.m_string[comp_str] = scal_str;
          comp++;
        }
    }

  // strain rates
  if (s_write_strains)
    {
      if (SpaceDim == 2)
        {
          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "du/dy + dv/dx";
          comp++;

          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "du/dx";
          comp++;

          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "dv/dy";
          comp++;
        }
    } // end if writing strains

  if (s_write_grad_eLambda)
    {
      sprintf(comp_str, "component_%d", comp);
      header.m_string[comp_str] = "grad_eLambdaX";
      comp++;

      if (SpaceDim >1)
        {
          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "grad_eLambdaY";
          comp++;

          if (SpaceDim >2)
            {
              sprintf(comp_str, "component_%d", comp);
              header.m_string[comp_str] = "grad_eLambdaZ";
              comp++;
            }
        }
    } // end if writing grad_eLambda

  if (s_write_error)
    {
      sprintf(comp_str, "component_%d", comp);
      header.m_string[comp_str] = "xVelErr";
      comp++;

      if (SpaceDim >1)
        {
          sprintf(comp_str, "component_%d", comp);
          header.m_string[comp_str] = "yVelErr";
          comp++;

          if (SpaceDim >2)
            {
              sprintf(comp_str, "component_%d", comp);
              header.m_string[comp_str] = "zVelErr";
              comp++;
            }
        }
    }

  // write procID's
  if (s_write_proc_ids)
    {
      sprintf(comp_str, "component_%d", comp);
      header.m_string[comp_str] = "procIDs";
      comp++;
    }

  header.writeToFile(a_handle);

    if (s_verbosity >= 3)
      {
        pout () << header << endl;
      }
}

// ---------------------------------------------------------------
void
AMRNavierStokes::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::writePlotLevel " << m_level << endl;
    }
  char level_str[20];
  sprintf (level_str, "%d", m_level);
  const std::string label = std::string("level_") + level_str;

  a_handle.setGroup(label);

  HDF5HeaderData header;

  header.m_int ["ref_ratio"]    = m_ref_ratio;
  header.m_real ["dx"]          = m_dx;
  header.m_real ["dt"]          = m_dt;
  header.m_real ["time"]        = m_time;
  header.m_box  ["prob_domain"] = m_problem_domain.domainBox();
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout () << header << endl;
  }

  const DisjointBoxLayout& levelGrids = m_vel_new_ptr->getBoxes();
  int numcomp = numPlotComps();

  LevelData<FArrayBox> plotData(levelGrids, numcomp, IntVect::Unit);
  getPlotData(plotData);

  plotData.exchange(plotData.interval());

  write (a_handle, levelGrids);
  write (a_handle, plotData, "data", IntVect::Unit);
}
#endif

// ---------------------------------------------------------------
int
AMRNavierStokes::numPlotComps() const
{
  // for now, hardwire number of components:
  // velocity
  int numcomp = SpaceDim;
  // divergence
  if (s_write_divergence)
  {
    ++numcomp;
  }
  // lambda
  if (s_write_lambda)
  {
    ++numcomp;
  }
  // du/dt
  if (s_write_time_derivatives)
  {
    numcomp += SpaceDim;
    // dLambda/dt
    ++numcomp;
  }

  // vorticity
  if (s_write_vorticity)
    {
      if (SpaceDim == 2)
        {
          ++numcomp;
        }
      else if (SpaceDim == 3)
        {
          numcomp += SpaceDim;
          // also include mag(vort)
          ++numcomp;
        }
    }

  // scalars
  if (s_write_scalars)
  {
    numcomp += s_num_scal_comps;
  }

  // dS/dt
  if (s_write_dScalar_dt)
  {
    numcomp += s_num_scal_comps;
  }

  // add in strain-rate stuff
  if (s_write_strains)
  {
    if (SpaceDim == 2)
    {
      numcomp += 3;
    }
  }

  // add in grad_eLambda
  if (s_write_grad_eLambda)
  {
    numcomp += SpaceDim;
  }

  // error
  if (s_write_error)
  {
    numcomp += SpaceDim;
  }

  // procID's
  if (s_write_proc_ids)
  {
    ++numcomp;
  }

  return numcomp;
}

// ---------------------------------------------------------------
void
AMRNavierStokes::getPlotData(LevelData<FArrayBox>& a_plot_data) const
{
  CH_assert (a_plot_data.nComp() == numPlotComps());
  const DisjointBoxLayout& levelGrids = a_plot_data.getBoxes();

  int plot_data_counter = 0;

  LevelData<FArrayBox> temp(levelGrids, 1);
  Interval srcComp(0,0);

  // first do velocities
  Interval velComps(0, SpaceDim-1);
  Interval destVelComps(plot_data_counter, plot_data_counter +SpaceDim-1);

  m_vel_new_ptr->copyTo(velComps, a_plot_data, destVelComps);
  plot_data_counter += SpaceDim;

  LevelData<FArrayBox>* fineVelPtr = NULL;
  LevelData<FArrayBox>* crseVelPtr = NULL;
  LevelData<FArrayBox> velTemp(levelGrids, SpaceDim);

  int nRefCrse = -1;

  if (m_level > 0)
    {
      crseVelPtr = crseNSPtr()->m_vel_new_ptr;
      nRefCrse = m_coarser_level_ptr->refRatio();
    }

  if (!finestLevel())
    {
      fineVelPtr = finerNSPtr()->m_vel_new_ptr;
    }

  // now compute divergence
  if (s_write_divergence)
    {
      m_vel_new_ptr->exchange(velComps);
      Divergence::compDivergenceCC(temp, *m_vel_new_ptr,
                                   crseVelPtr, fineVelPtr,
                                   m_dx, nRefCrse, m_ref_ratio,
                                   m_problem_domain,
                                   true);

      Interval divDestComp(plot_data_counter, plot_data_counter);
      temp.copyTo(srcComp, a_plot_data, divDestComp);
      plot_data_counter += 1;
    } // end divergence

  if (s_write_lambda)
    {
      // now copy lambda
      Interval lambdaDestComp(plot_data_counter, plot_data_counter);
      m_lambda_new_ptr->copyTo(srcComp, a_plot_data, lambdaDestComp);
      plot_data_counter += 1;
    }

  if (s_write_time_derivatives)
    {
      // now do du/dt and dLambda/dt
      m_vel_new_ptr->copyTo(velComps, velTemp, velComps);
      LevelData<FArrayBox> lambdaTemp(levelGrids,1);
      m_lambda_new_ptr->copyTo(m_lambda_new_ptr->interval(),
                               lambdaTemp, lambdaTemp.interval());

#ifdef DEBUG
      const LevelData<FArrayBox>& old_vel = *m_vel_save_ptr;
      const LevelData<FArrayBox>& old_lambda = *m_lambda_save_ptr;
#else
      const LevelData<FArrayBox>& old_vel = *m_vel_old_ptr;
      const LevelData<FArrayBox>& old_lambda = *m_lambda_old_ptr;
#endif

      if (m_dt_save != 0.0)
        {
          DataIterator dit = velTemp.dataIterator();
          for (dit.reset(); dit.ok(); ++dit)
            {
              velTemp[dit].minus(old_vel[dit], levelGrids[dit],
                                 0, 0, SpaceDim);

              lambdaTemp[dit].minus(old_lambda[dit], levelGrids[dit],
                                    0, 0, 1);

#ifdef DEBUG
              Real dt = m_time - m_saved_time;
              velTemp[dit].divide(dt);
              lambdaTemp[dit].divide(dt);
#else
              velTemp[dit].divide(m_dt_save);
              lambdaTemp[dit].divide(m_dt_save);
#endif
            }
        }
      else
        {
          setValLevel(velTemp, 0.0);
        }

      // set covered regions to 0
      if (!finestLevel())
        {
          fineVelPtr = finerNSPtr()->m_vel_new_ptr;
          const DisjointBoxLayout& fineGrids = fineVelPtr->getBoxes();

          DataIterator ditFine = fineGrids.dataIterator();
          DataIterator dit = velTemp.dataIterator();
          for (ditFine.begin(); ditFine.ok(); ++ditFine)
            {
              Box coarsenedFineBox(fineGrids[ditFine]);
              coarsenedFineBox.coarsen(m_ref_ratio);
              for (dit.reset(); dit.ok(); ++dit)
                {
                  Box intersectBox(levelGrids[dit]);
                  intersectBox &= coarsenedFineBox;
                  if (!intersectBox.isEmpty())
                    {
                      velTemp[dit].setVal(0.0, intersectBox, 0, SpaceDim);
                      lambdaTemp[dit].setVal(0.0, intersectBox, 0, 1);
                    }
                }
            }
        } // end setting covered regions to 0

      Interval duDtDestComps(plot_data_counter, plot_data_counter +SpaceDim-1);
      velTemp.copyTo(velComps, a_plot_data, duDtDestComps);
      plot_data_counter += SpaceDim;
      Interval dLambdaDtDestComps(plot_data_counter, plot_data_counter);
      lambdaTemp.copyTo(lambdaTemp.interval(), a_plot_data, dLambdaDtDestComps);
      plot_data_counter++;
    } // end d()/dt computations

  // finally do vorticity
  if (s_write_vorticity)
    {
      if (SpaceDim == 2)
        {
          computeVorticity(temp);
          Interval vortComp(plot_data_counter, plot_data_counter);
          temp.copyTo(srcComp, a_plot_data, vortComp);
          ++plot_data_counter;
        }
      else if (SpaceDim == 3)
        {
          computeVorticity(velTemp);
          Interval vortComps(plot_data_counter, plot_data_counter +SpaceDim-1);
          velTemp.copyTo(velComps, a_plot_data, vortComps);
          plot_data_counter += SpaceDim;
          // now compute norm of vorticity
          DataIterator dit = velTemp.dataIterator();
          for (dit.reset(); dit.ok(); ++dit)
            {
              // computes mag(velTemp) and places it in temp
              FORT_MAGVECT(CHF_FRA1(temp[dit],0),
                           CHF_CONST_FRA(velTemp[dit]),
                           CHF_BOX(levelGrids[dit]));
            }
          Interval magVortComps(plot_data_counter, plot_data_counter);
          temp.copyTo(srcComp,a_plot_data, magVortComps);
          ++plot_data_counter;
        }
    } // end vorticity

  // and last, but not least, do scalars
  if (s_write_scalars)
    {
      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          Interval scalComps = newScal(comp).interval();
          Interval scalDestComps(plot_data_counter,
                                 plot_data_counter+scalComps.size()-1);
          newScal(comp).copyTo(scalComps, a_plot_data, scalDestComps);
          plot_data_counter += scalComps.size();
        }
    }

  // and now more -- do ds/dt
  if (s_num_scal_comps > 0 && s_write_dScalar_dt)
    {
      LevelData<FArrayBox> scalTemp(levelGrids, s_num_scal_comps);

      // loop over components
      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          Interval scalSrcComps = newScal(comp).interval();
          Interval scalDestComps(comp,comp+scalSrcComps.size()-1);
          newScal(comp).copyTo(scalSrcComps, scalTemp, scalDestComps);
#ifdef DEBUG
          const LevelData<FArrayBox>& old_scal = *m_scal_save[comp];
#else
          const LevelData<FArrayBox>& old_scal = *m_scal_old[comp];
#endif

          if (m_dt_save != 0.0)
            {
              DataIterator dit = scalTemp.dataIterator();
              for (dit.reset(); dit.ok(); ++dit)
                {
                  scalTemp[dit].minus(old_scal[dit], levelGrids[dit],
                                      0, comp, scalSrcComps.size());
#ifdef DEBUG
                  Real dt = m_time - m_saved_time;
                  scalTemp[dit].divide(dt,comp,scalSrcComps.size());
#else
                  scalTemp[dit].divide(m_dt_save,comp,scalSrcComps.size());
#endif
                }
            }
          else
            {
              setValLevel(scalTemp, 0.0);
            }
        } // end loop over scalar components

      Interval dsDtDestComps(plot_data_counter,
                             plot_data_counter+s_num_scal_comps-1);
      scalTemp.copyTo(scalTemp.interval(), a_plot_data, dsDtDestComps);
      plot_data_counter += s_num_scal_comps;
    } // end if there are even scalars here and we want to write them

  // now put in strain-rate stuff.
  if (SpaceDim == 2 && s_write_strains)
    {
      // note that quadratic C/F BC's and extrap physical BCs
      // were set in computeVorticity

      // first add 1/2(du/dy + dv/dx) component
      DataIterator dit = temp.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = levelGrids[dit];
          // du/dy
          int derivDir = 1;
          FORT_CENTERED_DERIV(CHF_FRA1(temp[dit],0),
                              CHF_CONST_FRA1(newVel()[dit],0),
                              CHF_BOX(gridBox),
                              CHF_CONST_REAL(m_dx),
                              CHF_INT(derivDir));
          temp[dit] *= 0.5;
        }

      Interval currentComp(plot_data_counter, plot_data_counter);
      temp.copyTo(temp.interval(), a_plot_data, currentComp);

      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = levelGrids[dit];
          // dv/dx
          int derivDir = 0;
          FORT_CENTERED_DERIV(CHF_FRA1(temp[dit],0),
                              CHF_CONST_FRA1(newVel()[dit],1),
                              CHF_BOX(gridBox),
                              CHF_CONST_REAL(m_dx),
                              CHF_INT(derivDir));

          temp[dit] *= 0.5;
          a_plot_data[dit].plus(temp[dit],0,plot_data_counter,1);
        }
      ++plot_data_counter;

      // now add du/dx
      Interval dudxcomp(plot_data_counter, plot_data_counter);
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = levelGrids[dit];
          int derivDir = 0;
          FORT_CENTERED_DERIV(CHF_FRA1(temp[dit],0),
                              CHF_CONST_FRA1(newVel()[dit],0),
                              CHF_BOX(gridBox),
                              CHF_CONST_REAL(m_dx),
                              CHF_INT(derivDir));
        }
      temp.copyTo(temp.interval(), a_plot_data, dudxcomp);
      ++plot_data_counter;

      // now dv/dy
      Interval dvdycomp(plot_data_counter, plot_data_counter);
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = levelGrids[dit];
          int derivDir = 1;
          FORT_CENTERED_DERIV(CHF_FRA1(temp[dit],0),
                              CHF_CONST_FRA1(newVel()[dit],1),
                              CHF_BOX(gridBox),
                              CHF_CONST_REAL(m_dx),
                              CHF_INT(derivDir));
        }
      temp.copyTo(temp.interval(), a_plot_data, dvdycomp);
      ++plot_data_counter;

    } // end strain-rate stuff
  else
    if (s_write_strains)
      {
        pout() << "Strain-rate components for plotfiles only in 2D!!!"
               << endl;
      }

  // finally, add in freestream preservation correction
  if (s_write_grad_eLambda)
    {
      const LevelData<FluxBox>& grad_eLambda = m_projection.grad_eLambda();

      // now do edge-to-cell averaging
      DataIterator dit = a_plot_data.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& thisPlotData = a_plot_data[dit];
          const FluxBox& thisGradE = grad_eLambda[dit];
          for (int dir=0; dir<SpaceDim; ++dir)
            {
              int edgeComp = 0;
              EdgeToCell(thisGradE,
                         edgeComp,
                         thisPlotData,
                         plot_data_counter+dir,
                         dir);
            }
        }

      plot_data_counter += SpaceDim;
    } // end writing freestream preservation correction

  if (s_write_error)
    {
      // first compute exact solution -- use velTemp storage for exact soln
      const LevelData<FArrayBox>& vel = newVel();
      DataIterator dit = a_plot_data.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& thisExact = velTemp[dit];
          FORT_COMPUTEEXACT(CHF_FRA(thisExact),
                            CHF_CONST_REAL(s_nu),
                            CHF_CONST_REAL(m_time),
                            CHF_CONST_REAL(m_dx));

          // now subtract velocity
          velTemp[dit].minus(vel[dit],0,0,SpaceDim);
        }

      // should we set covered regions to 0?

      // now copy error into big fab
      Interval errComps(plot_data_counter, plot_data_counter+SpaceDim-1);
      velTemp.copyTo(velTemp.interval(), a_plot_data, errComps);

      plot_data_counter += SpaceDim;
    }

  // procID's
  if (s_write_proc_ids)
    {
      Real procVal = (Real) procID();

      DataIterator dit = a_plot_data.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& thisPlotData = a_plot_data[dit];
          thisPlotData.setVal(procVal, plot_data_counter);
        }

      ++plot_data_counter;
    }
}
