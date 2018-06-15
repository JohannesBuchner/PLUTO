#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "VisItChomboDriver.H"
#include "MayDay.H"
#include <cstdio>

#include "NamespaceHeader.H"

VisItChomboDriver::VisItChomboDriver()
{
    initialized = false;
}

VisItChomboDriver::~VisItChomboDriver()
{
}

void
VisItChomboDriver::VisualizeFile(const char *fname)
{
    VisualizeFileHelper(fname, true, Pseudocolor);
}

void
VisItChomboDriver::BrowseFile(const char *fname)
{
    VisualizeFileHelper(fname, true, Spreadsheet);
}

void
VisItChomboDriver::VisualizeFileHelper(const char *fname,
                                       bool allowRetry,
                                       VisualizationType vistype)
{
    bool success = true;
    if (!visit.IsOpen())
    {
        initialized = false;
        success = visit.Open();
        if (!success)
        {
            visit.Close();
            initialized = false;
            MayDay::Warning("VisualizeFile: could not open connection to VisIt");
            return;
        }
    }

    char opencmd[2000];
    if (!initialized)
    {
        success = success && visit.SendCommand("OpenGUI()");

        // If you wanted to set some default plot options, like
        // turning off the tracer plane in spreadsheet plots, here
        // is the perfect place to do that:
        //success = success && visit.SendCommand("spreadsheetatts = SpreadsheetAttributes()");
        //success = success && visit.SendCommand("spreadsheetatts.showTracerPlane=0");
        //success = success && visit.SendCommand("SetDefaultPlotOptions(spreadsheetatts)");
        initialized = true;
    }
    else
    {
        success = success && visit.SendCommand("AddWindow()");
        success = success && visit.SendCommand("DeleteAllPlots()");
    }
    sprintf(opencmd, "OpenDatabase('%s', 0, 'Chombo_1.0')", fname);
    success = success && visit.SendCommand(opencmd);
    switch (vistype)
    {
      case Pseudocolor:
        success = success && visit.SendCommand("AddPlot('Pseudocolor','component_0')");
        break;
      case Spreadsheet:
        success = success && visit.SendCommand("AddPlot('Spreadsheet','component_0')");
        break;
    }
    success = success && visit.SendCommand("DrawPlots()");
    if (!success)
    {
        visit.Close();
        initialized = false;
        if (allowRetry)
        {
            VisualizeFileHelper(fname, false, vistype);
        }
        else
        {
            MayDay::Warning("VisualizeFile: could not communicate commands to VisIt");
        }
    }
}

void
VisItChomboDriver::Reset()
{
    visit.Close();
    initialized = false;
}

#include "NamespaceFooter.H"
