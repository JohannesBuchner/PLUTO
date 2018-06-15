#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VISITCHOMBODRIVER_H_
#define _VISITCHOMBODRIVER_H_

#include "VisItPythonConnection.H"

#include "NamespaceHeader.H"

class VisItChomboDriver
{
  public:
          VisItChomboDriver();

         ~VisItChomboDriver();

    void  VisualizeFile(const char *fname);

    void  BrowseFile(const char *fname);

    void  Reset();

  private:
    enum VisualizationType
    {
      Pseudocolor,
      Spreadsheet
    };

    void VisualizeFileHelper(const char *fname,
                             bool allowRetry,
                             VisualizationType vistype);

    VisItPythonConnection visit;

    bool                  initialized;
};

#include "NamespaceFooter.H"

#endif
