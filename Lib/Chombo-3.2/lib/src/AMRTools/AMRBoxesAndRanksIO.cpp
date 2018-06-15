#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Vector.H"
#include "Box.H"
#include "SPMD.H"
#include "MayDay.H"

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include "parstream.H"

#ifdef MPI
#include <mpi.h>
#endif

#include "NamespaceHeader.H"

using std::cout;
using std::cerr;
using std::ifstream;
using namespace std;

bool parseLine(string s, const string key, float& value, bool& flag);
bool parseLine(string s, const string key, int& value, bool& flag);
bool parseLine(string s, const string key, Box& value, bool& flag);
bool parseMultiValueLine(string s, const string key, vector<int>& value, const int nvalues, bool& flag);
bool isCommentOrBlankLine(const string s);
void checkForDuplicateKey(bool flag, const string key);
bool getNextValidLine(std::ifstream& ifs, string& sline);
int determineNumberBoxesOnLine(const string sline);
int countChar(const string s, const char* c);

void writeABRfile(const Vector<Vector<Box> >&      a_amrBoxes,
                  const Vector<Vector<int> >&      a_amrRanks,
                  const Vector<int>&               a_refinementRatio,
                  const int                        a_numLevels,
                  const int                        a_targetNumProcs,
                  const Box&                       a_baseLevelBoundingBox,
                  const string a_filename)
{
  // assumes this will be run in serial....

  pout() << "writeABRfile  filename=" << a_filename << "\n";

  float versionABRfileFormat = 0.1;
  const int nboxesPerLine = 3;

  std::fstream ABR(a_filename.c_str(), ios::out);
  ABR << "! Written by gridGenerator.cpp via writeABRfile()" << "\n";
  ABR << "ABR file format version " << versionABRfileFormat << "\n";
  ABR << "\n";
  ABR << "Number of Levels = " << a_numLevels << "\n";
  ABR << "\n";
  ABR << "! nlevels-1 values" << "\n";
  ABR << "Refinement Ratio = " << a_refinementRatio << "\n";
  ABR << "\n";
  ABR << "Cell Centering = 0" << "\n";
  ABR << "\n";
  ABR << "Number of Time Steps = 1" << "\n";
  ABR << "\n";
  ABR << "Number of Processors = " << a_targetNumProcs << "\n";
  ABR << "\n";
  ABR << "Level 0 Domain = " << a_baseLevelBoundingBox.smallEnd()
      << a_baseLevelBoundingBox.bigEnd() << "\n";
  ABR << "\n";
  ABR << "TimeStep: 0" << "\n";
  ABR << "\n";

  int rank;
  for (int ilev=0; ilev < a_numLevels; ilev++)
    {
      const int Nboxes = a_amrBoxes[ilev].size();
      //if (Nboxes <= 0) break;
      pout() << "Number of Boxes on Level " << ilev << ": " << Nboxes << "\n";
      ABR << "Number of Boxes on Level " << ilev << ": " << Nboxes << "\n";
      CH_assert(Nboxes > 0); // ???
      for (int ibox=0; ibox<Nboxes; ibox++)
      {
        ABR << a_amrBoxes[ilev][ibox].smallEnd() << a_amrBoxes[ilev][ibox].bigEnd();
        rank = a_amrRanks[ilev][ibox];
        ABR << "[" << rank << "]";
        if ((ibox+1)==Nboxes)
        {
          // if last one
          ABR << "\n";
          break;
        }
        if ((ibox+1)%nboxesPerLine == 0)
        {
          ABR << "\n";
        }
        else
        {
          ABR << " ";
        }
      }
      ABR << "\n";
    }

  //(0,0,0)(31,31,31)[0]

  ABR.close();
  pout() << "writeABRfile done" << "\n";
}

// Default is one ABR file for all procs -- so proc0 reads in and broadcasts info.
// Set dataPerRank to true to read in entire ABR file data to be used for your proc,
//   ie read one ABR file per processor
//
void readABRfile(Vector<Vector<Box> >&      a_amrBoxes,
                 Vector<Vector<int> >&      a_amrRanks,
                 Vector<int>&               a_refinementRatio,
                 Box&                       a_baseLevelBoundingBox,
                 const string               a_filename,
                 const bool                 a_dataPerRank=false)
{
  int verbose = 1;

  // answers:
  float version;
  int nlevels;
  int cellcentering;
  vector<int> refratio;
  int nprocs;
  Box level0BoundingBox;
  int ntimesteps;
  //Vector<Vector<Box> > amrBoxes;
  //Vector<Vector<int> > amrRanks;

  //barrier();
  bool troubleOpeningFile=false;

  if (a_dataPerRank || procID() == uniqueProc(SerialTask::compute))
    {
      // read boxes (and proc assignments) from file
      std::ifstream ifs(a_filename.c_str(), std::ios_base::in);
      if (ifs.fail())
      {
        troubleOpeningFile=true;
      }
      ifs.close();
    }

  if (troubleOpeningFile)
    {
      char buffer[300];
      sprintf(buffer,"readABRfile:  Cannot open grids file: %s",a_filename.c_str());
      MayDay::Error(buffer);
    }

  if (a_dataPerRank || procID() == uniqueProc(SerialTask::compute))
    {
      std::ifstream ifs(a_filename.c_str(), std::ios_base::in);
      if (verbose) pout() << "procID: " << procID() << "  opening gridfile" << "\n";

      int maxRankFromFile = -1;
      bool readyForBoxes = false;
      string sline;

      // ABR file format version 0.1
      const string sversion0 = "DBL file format version";
      const string sversion1 = "ABR file format version";
      bool bversion=false;

      // Number of Levels = 3
      const string snlevels = "Number of Levels =";
      nlevels=-1;
      bool bnlevels=false;

      // Refinement Ratio = 2 2
      const string srefratio = "Refinement Ratio =";
      bool brefratio=false;

      // Cell Centering = 0
      const string scellcentering = "Cell Centering =";
      cellcentering=0;
      bool bcellcentering=false;

      // Number of Processors = 1
      const string snprocs = "Number of Processors =";
      nprocs=1;
      bool bnprocs=false;

      //Problem is Periodic in Direction = 0 0 0

      //Level 0 Domain = (0,0,0)(31,31,31)
      const string slevel0BoundingBox = "Level 0 Domain =";
      bool blevel0BoundingBox=false;

      // Number of Time Steps = 1
      const string sntimesteps = "Number of Time Steps =";
      ntimesteps=1;
      bool bntimesteps=false;

      // TimeStep: 0
      const string stimestep = "TimeStep:";
      int timestep;
      bool btimestep=false;

      string saverefratio="";
      string saveL0domain;

      while (!readyForBoxes)
      {
        getNextValidLine(ifs, sline);
        //cout << sline << "\n";

        if (parseLine(sline, sversion0, version, bversion))
        {
          //cout << " version = " << version << "\n";
        } else if (parseLine(sline, sversion1, version, bversion))
        {
          //cout << " version = " << version << "\n";
        } else if (parseLine(sline, snlevels, nlevels, bnlevels))
        {
          //cout << " nlevels = " << nlevels << "\n";
        } else if (parseLine(sline, scellcentering, cellcentering, bcellcentering))
        {
          //cout << " cellcentering = " << cellcentering << "\n";
        } else if (parseLine(sline, sntimesteps, ntimesteps, bntimesteps))
        {
          //cout << " ntimesteps = " << ntimesteps << "\n";
        } else if (parseLine(sline, snprocs, nprocs, bnprocs))
        {
          //cout << " nprocs = " << nprocs << "\n";
        } else if (parseLine(sline, slevel0BoundingBox, level0BoundingBox, blevel0BoundingBox))
        {
          //cout << " level0BoundingBox = " << level0BoundingBox << "\n";
          stringstream ist1;
          ist1 << level0BoundingBox.smallEnd() << " " << level0BoundingBox.bigEnd();
          saveL0domain = ist1.str();
        }
        else if (parseMultiValueLine(sline, srefratio, refratio, nlevels-1, brefratio))
        {

          stringstream ist1;
          for (int ilev=0; ilev<nlevels-1; ilev++)
          {
            ist1 << refratio[ilev] << " ";
          }
          saverefratio = ist1.str();
          // gotta be a cleaner way to do this...
          //char cc[10];
          //for (int ilev=0; ilev<nlevels-1; ilev++)
          //{
          //  sprintf(cc, " %d", refratio[ilev]);
          //  saverefratio += cc;
          //  cout << refratio[ilev] << " (" << saverefratio << ")" <<"\n";
          //}
          //cout << " refratio = " << saverefratio << "\n";
        }
        else
        {
          pout() << sline << "\n";
          MayDay::Error("Found other line. what is this?");
        }

        readyForBoxes = (bversion && bnlevels && brefratio && blevel0BoundingBox
                         && bcellcentering && bnprocs && bntimesteps);
        //cout << " readyForBoxes=" << readyForBoxes << "\n";
      }

      stringstream ist1;
      ist1 << level0BoundingBox.smallEnd() << " " << level0BoundingBox.bigEnd();
      string sL0domain = ist1.str();

      // print nifty one-liner
      char summary[160];
      sprintf(summary,
              "%s v%-5.1f Nlevels=%d Nprocs=%d centering=%d ts=%d L0domain=%s rr=%s\n",
              a_filename.c_str(), version, nlevels, nprocs, cellcentering, ntimesteps,
              sL0domain.c_str(), saverefratio.c_str());
      pout() << summary << "\n";

      // Now read all box info.  For each timestep and each level
      for (int ts=0; ts < ntimesteps; ts++)
      {
        // Look for TimeStep:
        while (getNextValidLine(ifs, sline))
        {
          if (parseLine(sline, stimestep, timestep, btimestep))
          {
            if (timestep != ts)
            {
              MayDay::Error("TimeStep values are currently required to be consecutive and start at 0");
            }
            //pout() << "TimeStep: " << ts << "\n";
            break;
          }
        }
        if (btimestep == false)
        {
          MayDay::Error("TimeStep line not found. btimestep==false");
        }

        a_amrBoxes.resize(nlevels);
        a_amrRanks.resize(nlevels);

        for (int ilev=0; ilev < nlevels; ilev++)
        {

          int nboxes=0;
          bool bnboxes=false;

          // Number of Boxes on Level 0: 1
          char cc[100];
          sprintf(cc, "Number of Boxes on Level %d:", ilev);
          const string snboxes(cc);

          // Look for Number of Boxes on Level ilev:
          while (getNextValidLine(ifs, sline))
          {

            //cout << " nob line sline =" << sline << "\n";
            if (parseLine(sline, snboxes, nboxes, bnboxes))
            {
              pout() << "Nboxes on level " << ilev << " = " << nboxes << "\n";
              break;
            }
          }
          //CH_assert(nboxes > 0);  //???

          //cout << "ilev=" << ilev << " nboxes=" << nboxes << "\n";
          a_amrBoxes[ilev].resize(nboxes);
          a_amrRanks[ilev].resize(nboxes);

          //(16,56,72)(31,71,79)[0] (16,56,80)(31,71,95)[0]

          int ibox=0;
          while (ibox < nboxes)
          {
            //cout << "ibox=" << ibox << "\n";

            getNextValidLine(ifs, sline);
            //cout << "box line sline=" << sline << "\n";

            int npl = determineNumberBoxesOnLine(sline);
            //cout << " npl=" << npl << "\n";

            // first weak attempt at ensuring data is in correct format for stream parsing
            CH_assert(countChar(sline, "(") == countChar(sline, ")"));
            CH_assert(countChar(sline, "[") == countChar(sline, "]"));

            istringstream ist(sline);
            for (int i=0; i<npl; i++)
            {
              IntVect smallend, bigend;
              ist >> smallend;
              ist >> bigend;

              Box tb(smallend, bigend);
              //cout << "box tb=" << tb << "\n";

              int rank;
              while (char ch = ist.get()) if (ch=='[') break;
              ist >> rank;
              while (char ch = ist.get()) if (ch==']') break;
              //cout << "rank=" << rank << "\n";

              if (rank > maxRankFromFile) maxRankFromFile=rank;

              a_amrBoxes[ilev][ibox] = tb;
              a_amrRanks[ilev][ibox] = rank;
              ibox++;
            }
          } // ibox

        } // ilev
      } // ts


      a_baseLevelBoundingBox = level0BoundingBox;

      a_refinementRatio.resize(nlevels-1);
      for (int ilev=0; ilev < nlevels-1; ilev++)
        {
          a_refinementRatio[ilev] = refratio[ilev];
        }


      pout() << " Max rank in grids file = " << maxRankFromFile
             << "  num processors = " << nprocs << "\n";

    } // rank0


  // answers:
  //float version;
  //int nlevels;
  //int cellcentering;
  //vector<int> refratio;
  //int nprocs;
  //Box level0BoundingBox;
  //int ntimesteps;
  //Vector<Vector<Box> > amrBoxes;
  //Vector<Vector<int> > amrRanks;

  // Not fully tested code, but it appears to be working
#ifdef CH_MPI
  if (!a_dataPerRank)
  {
    // Otherwise, don't need to broadcast
    //pout() << " broadcast nlevels & nprocs, etc" << "\n";
    broadcast(nlevels, 0);
    broadcast(nprocs, 0);

    broadcast(a_refinementRatio, 0);
    broadcast(cellcentering, 0);
    broadcast(a_baseLevelBoundingBox, 0);
    broadcast(ntimesteps, 0);

    a_amrBoxes.resize(nlevels);
    a_amrRanks.resize(nlevels);

    for (int ilev=0; ilev < nlevels; ilev++)
      {
        int nboxes;
        if (procID() == 0)
        {
          nboxes = a_amrBoxes[ilev].size();
        }
        broadcast(nboxes, 0);
        a_amrBoxes[ilev].resize(nboxes);
        a_amrRanks[ilev].resize(nboxes);

        broadcast(a_amrBoxes[ilev], 0);
        broadcast(a_amrRanks[ilev], 0);

        //for (int b=0; b<a_amrBoxes[ilev].size(); b++)
        //{
        //  broadcast(a_amrBoxes[ilev][b], 0);
        //  broadcast(a_amrRanks[ilev][b], 0);
        //}
      }
  }
#endif

  pout() << "Done reading in ABR file (AMRBoxesAndRanksIO.cpp)\n";

}


int determineNumberBoxesOnLine(const string sline)
{
  // count pairs of open ('s for now...
  int npl = countChar(sline, "(");
  if (npl%2 != 0)
  {
    MayDay::Error("Expected pairs of open parens");
  }
  npl /= 2;
  if (npl <= 0)
  {
    MayDay::Error("Expected at least one box on line");
  }
  //CH_assert(nc%2 == 0);
  //CH_assert(npl > 0);
  //cout << " npl=" << npl << "\n";
  return npl;
}

bool getNextValidLine(std::ifstream& ifs, string& sline)
{
  bool ok;
  while ((ok = (bool)getline(ifs, sline, '\n')))
  {
    if (isCommentOrBlankLine(sline))
    {
      continue;
    }
    else
    {
      //cout << "getNextValidLine=(" << sline << ")" << "\n";
      break;
    }
    //if EOF FAIL?
  }
  return ok;
}

// fix this to account for tabs as well.
bool isCommentOrBlankLine(const string s)
{
  int i=0;
  while (s[i++] == ' '); i--;
  if (s[i] == '!' || i == s.size())
  {
    return true;
  }
  else
  {
    return false;
  }
}

int countChar(const string s, const char* c)
{
  int nc=0;
  string::size_type st=0;
  while ( (st = s.find(c[0], st)) != string::npos )
  {
    //cout << " st=" << st << "\n";
    st++;  // so we don't find same one again?
    nc++;
    //if (nc > 10) break;
  }
  //CH_assert(nc > 0);
  //nc /= 2;
  //cout << "nc = " << nc << "\n";
  return nc;
}

// float
bool parseLine(string s, const string key, float& value, bool& flag)
{
  int i1;
  if ( (i1 = s.find(key)) >= 0)
  {
    //cout << "Found key:" << key << " @ position" << i1 << "\n";
    s.replace(i1, key.size(), "");
    value = atof(s.c_str());
    //cout << " left with " << " (" << s << ") value = " << value << "\n";
    checkForDuplicateKey(flag, key);
    flag = true;
    return true;
  }
  else
  {
    return false;
  }
}

// int
bool parseLine(string s, const string key, int& value, bool& flag)
{
  int i1;
  if ( (i1 = s.find(key)) >= 0)
  {
    //cout << "Found key:" << key << " @ position" << i1 << "\n";
    s.replace(i1, key.size(), "");
    value = atoi(s.c_str());
    //cout << " left with " << " (" << s << ") value = " << value << "\n";
    checkForDuplicateKey(flag, key);
    flag = true;
    return true;
  }
  else
  {
    return false;
  }
}

// Box
bool parseLine(string s, const string key, Box& value, bool& flag)
{
  int i1;
  if ( (i1 = s.find(key)) >= 0)
  {
    //cout << "Found key:" << key << " @ position" << i1 << "\n";
    s.replace(i1, key.size(), "");
    istringstream ist(s);
    IntVect smallend, bigend;
    ist >> smallend;
    ist >> bigend;
    Box tb(smallend, bigend);
    value = tb;
    //cout << " left with " << " (" << s << ") value = " << value << "\n";
    checkForDuplicateKey(flag, key);
    flag = true;
    return true;
  }
  else
  {
    return false;
  }
}


bool parseMultiValueLine(string s, const string key, vector<int>& value, const int nvalues, bool& flag)
{

  if (nvalues < 0)
  {
    MayDay::Error("nvalues <= 0  (Must set Number of Levels before Refine Ratio)");
  }

  int i1;
  if ( (i1 = s.find(key)) >= 0)
  {
    //cout << "Found key:" << key << " @ position" << i1 << "\n";
    s.replace(i1, key.size(), "");

    // needs much more work!
    int j=0;
    vector<string> vst(nvalues);
    for (int ilev=0; ilev<nvalues; ilev++)
    {
      while (s[j++] == ' ');
      j--;
      while (s[j] != ' ')
      {
        vst[ilev] += s[j];
        j++;
      }
      j--;
    }

    value.resize(nvalues);
    for (int ilev=0; ilev<nvalues; ilev++)
    {
      value[ilev] = atoi(vst[ilev].c_str());
      //cout << "ilev= " << ilev << " value=" << value[ilev] << "\n";
    }
    //cout << " left with " << " (" << s << ") value = " << value << "\n";
    checkForDuplicateKey(flag, key);
    flag = true;
    return true;
  }
  else
  {
    return false;
  }
}

void checkForDuplicateKey(bool flag, const string key)
{
  if (flag != false)
  {
    pout() << " Found duplicate key = " << key << "\n";
    MayDay::Error("Found input parameter duplicated");
  }
}



// //  from IntVect.cpp
// #define CH_IGNORE_MAX 100000

// istream&
// operator>> (istream& is,
//             IntVect& p)
// {
//     is >> ws;
//     char c;
//     is >> c;
//     is.putback(c);
//     if (c == '(')
//     {
//         D_EXPR(is.ignore(CH_IGNORE_MAX, '(') >> p[0],
//                is.ignore(CH_IGNORE_MAX, ',') >> p[1],
//                is.ignore(CH_IGNORE_MAX, ',') >> p[2]);
//         is.ignore(CH_IGNORE_MAX, ')');
//     }
//     else if (c == '<')
//     {
//         D_EXPR(is.ignore(CH_IGNORE_MAX, '<') >> p[0],
//                is.ignore(CH_IGNORE_MAX, ',') >> p[1],
//                is.ignore(CH_IGNORE_MAX, ',') >> p[2]);
//         is.ignore(CH_IGNORE_MAX, '>');
//     }
//     else
//         MayDay::Error("operator>>(istream&,IntVect&): expected \'(\' or \'<\'");

//     if (is.fail())
//         MayDay::Error("operator>>(istream&,IntVect&) failed");

//     return is;
// }

#include "NamespaceFooter.H"
