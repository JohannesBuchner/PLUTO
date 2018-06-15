#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Scheduler.H"
#include <float.h>
#include "AMR.H"
#include <fstream>
#include <dirent.h> // For reading plot files from the cwd. -JNJ
#include "NamespaceHeader.H"
using namespace std;

//-----------------------------------------------------------------------
Scheduler::PeriodicFunction::
PeriodicFunction()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Scheduler::PeriodicFunction::
~PeriodicFunction()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
Scheduler::PeriodicFunction::
setUp(AMR& a_AMR, int a_interval)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
Scheduler::PeriodicFunction::
setUp(AMR& a_AMR, Real a_interval)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
Scheduler::PeriodicFunction::
conclude(int a_step, Real a_time)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Scheduler::
Scheduler():
  m_stepTriggeredFunctions(),
  m_timeTriggeredFunctions()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Scheduler::
~Scheduler()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
Scheduler::
schedule(RefCountedPtr<PeriodicFunction> a_function,
         int a_interval)
{
  m_stepTriggeredFunctions[a_function] = a_interval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
Scheduler::
schedule(RefCountedPtr<PeriodicFunction> a_function,
         Real a_interval)
{
  m_timeTriggeredFunctions[a_function] = pair<Real, Real>(a_interval, -FLT_MAX);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
Scheduler::
execute(int a_step, Real a_time) const
{
  // Step through our step-triggered functions and call those that match.
  for (map<RefCountedPtr<PeriodicFunction>, int, PeriodicFunctionLessThan>::iterator
       iter = m_stepTriggeredFunctions.begin(); iter != m_stepTriggeredFunctions.end(); ++iter)
  {
    if ((a_step % iter->second) == 0)
    {
      PeriodicFunction& function = const_cast<PeriodicFunction&>((*iter->first));
      function(a_step, a_time);
    }
  }

  // Step through our time-triggered functions and call those that match.
  for (map<RefCountedPtr<PeriodicFunction>, pair<Real, Real>, PeriodicFunctionLessThan>::iterator
       iter = m_timeTriggeredFunctions.begin(); iter != m_timeTriggeredFunctions.end(); ++iter)
  {
    Real interval = iter->second.first;
    Real tLast = iter->second.second;
    if ((a_time - tLast) >= interval)
    {
      // Call the function.
      PeriodicFunction& function = const_cast<PeriodicFunction&>((*iter->first));
      function(a_step, a_time);

      // Mark this as the last time the function was called.
      iter->second.second = a_time;
    }
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
Scheduler::
setUp(AMR& a_AMR)
{
  // Step through our step-triggered functions and call setUp().
  for (map<RefCountedPtr<PeriodicFunction>, int, PeriodicFunctionLessThan>::iterator
       iter = m_stepTriggeredFunctions.begin(); iter != m_stepTriggeredFunctions.end(); ++iter)
  {
    PeriodicFunction& function = const_cast<PeriodicFunction&>((*iter->first));
    int interval = iter->second;
    function.setUp(a_AMR, interval);
  }

  // Step through our time-triggered functions and call setUp().
  for (map<RefCountedPtr<PeriodicFunction>, pair<Real, Real>, PeriodicFunctionLessThan>::iterator
       iter = m_timeTriggeredFunctions.begin(); iter != m_timeTriggeredFunctions.end(); ++iter)
  {
    PeriodicFunction& function = const_cast<PeriodicFunction&>((*iter->first));
    Real interval = iter->second.first;
    function.setUp(a_AMR, interval);
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
Scheduler::
conclude(int a_step, Real a_time) const
{
  // Step through our step-triggered functions and call conclude().
  for (map<RefCountedPtr<PeriodicFunction>, int, PeriodicFunctionLessThan>::iterator
       iter = m_stepTriggeredFunctions.begin(); iter != m_stepTriggeredFunctions.end(); ++iter)
  {
    PeriodicFunction& function = const_cast<PeriodicFunction&>((*iter->first));
    function.conclude(a_step, a_time);
  }

  // Step through our time-triggered functions and call conclude().
  for (map<RefCountedPtr<PeriodicFunction>, pair<Real, Real>, PeriodicFunctionLessThan>::iterator
       iter = m_timeTriggeredFunctions.begin(); iter != m_timeTriggeredFunctions.end(); ++iter)
  {
    PeriodicFunction& function = const_cast<PeriodicFunction&>((*iter->first));
    function.conclude(a_step, a_time);
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PlotterPeriodicFunction::
PlotterPeriodicFunction(const std::string& a_prefix):
  Scheduler::PeriodicFunction(),
  m_prefix(a_prefix)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PlotterPeriodicFunction::
operator()(int a_step, Real a_time)
{
  // This doesn't actually write the plot file, since that is done by
  // the AMR object. However, it does aggregate all files into one file
  // understandable by Visit.

  // Find all the plot files with this prefix in the current working
  // directory.
  vector<string> plotFiles;
  DIR* dir = opendir(".");
  struct dirent* entry;
  CH_assert(dir != NULL);
  entry = readdir(dir);
  while (entry != NULL)
  {
    // Get the filename and add it to our list of plot files if
    // it begins with the prefix and ends in .hdf5.
    // NOTE: We exclude files containing '.map.hdf5' at the end.
    string name(entry->d_name);
    if ((name.find(m_prefix) == 0) &&
        (name.find(".hdf5") == name.length() - 5) &&
        (name.find(".map.hdf5") == string::npos))
      plotFiles.push_back(name);

    // Get the next one.
    entry = readdir(dir);
  }
  closedir(dir);

#ifdef CH_MPI
  int rank = 0;
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  if (rank == 0)
  {
#endif
  // Write all the plot file names into an aggregate file.
  if (!plotFiles.empty())
  {
    // Sort the list of files by their names.
    sort(plotFiles.begin(), plotFiles.end());

    // The aggregate file ends in .visit. Avoid double dots.
    string filename;
    if (m_prefix[m_prefix.length()-1] == '.')
      filename = m_prefix + string("visit");
    else
      filename = m_prefix + string(".visit");

    ofstream aggFile(filename.c_str());
    for (int i = 0; i < plotFiles.size(); ++i)
      aggFile << plotFiles[i] << endl;
    aggFile.close();
  }
#ifdef CH_MPI
  }
#endif
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PlotterPeriodicFunction::
setUp(AMR& a_AMR, int a_interval)
{
  a_AMR.plotPrefix(m_prefix);
  a_AMR.plotInterval(a_interval);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PlotterPeriodicFunction::
setUp(AMR& a_AMR, Real a_interval)
{
  a_AMR.plotPrefix(m_prefix);
  a_AMR.plotPeriod(a_interval);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
PlotterPeriodicFunction::
conclude(int a_step, Real a_time)
{
  // One last plot!
  (*this)(a_step, a_time);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
CheckpointPeriodicFunction::
CheckpointPeriodicFunction(const std::string& a_prefix):
  Scheduler::PeriodicFunction(),
  m_prefix(a_prefix)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
CheckpointPeriodicFunction::
setUp(AMR& a_AMR, int a_interval)
{
  // Set the checkpoint prefix and interval.
  a_AMR.checkpointPrefix(m_prefix);
  a_AMR.checkpointInterval(a_interval);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
CheckpointPeriodicFunction::
operator()(int a_step, Real a_time)
{
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
