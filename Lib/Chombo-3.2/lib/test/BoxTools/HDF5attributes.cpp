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
#include <cstring>

#include "CH_HDF5.H"
#include "SPMD.H"
#include "parstream.H"

#include "UsingNamespace.H"

using std::endl;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
test();

/// Global variables for handling output:
static const char *pgmname = "HDF5attributes" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

#ifdef CH_USE_HDF5
int compare(const HDF5HeaderData& lhs, const HDF5HeaderData& rhs);
#endif

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions( argc ,argv ) ;

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int icode = test();
  if (icode != 0)
    {
      pout() << indent << pgmname <<" failed"<<endl;
    }
  else
    {
      pout() << indent << pgmname <<" passed"<<endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}

// returns 0 on all tests passed.

int test()
{
#ifdef CH_USE_HDF5

  int error;
  HDF5Handle testFile;

 CH_assert(!testFile.isOpen());

  error = testFile.open("attributeTest.h5", HDF5Handle::CREATE);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "File creation failed "<<error<<endl;
      return error;
    }

 CH_assert(testFile.isOpen());

  HDF5HeaderData set1;

  Real dx=0.004;
  Box b1(IntVect(D_DECL6(1,2,1,1,2,1)), IntVect(D_DECL6(4,4,4,4,4,4)));
  Box b2(IntVect(D_DECL6(5,2,1,5,2,1)), IntVect(D_DECL6(12,4,4,12,4,4)));
  int currentStep = 2332;

  set1.m_string["name"] = "set1";
  set1.m_real["dx"] = dx;
  set1.m_int["currentStep"] = currentStep;
  set1.m_intvect["some intvect or other"] = b1.smallEnd();
  set1.m_box["b1"] = b1;
  set1.m_box["b2"] = b2;

  error = set1.writeToFile(testFile); //write information to file root.
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "attribute write failed "<<error<<endl;
      return error;
    }
  testFile.setGroupToLevel(2);

  error = set1.writeToFile(testFile);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "attribute write failed "<<error<<endl;
      return error;
    }

  error = testFile.setGroup("/myspecialgroup");
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "set group failed "<<error<<endl;
      return error;
    }

  testFile.pushGroup("inner_group1");
  int i = 31;
  Real x = 1.23456;
  // herr_t status;
  {
    HDF5HeaderData header;
    header.m_int["an_integer"] = i;
    header.writeToLocation( testFile.groupID() );
  }
  testFile.pushGroup("inner_group2");
  {
    HDF5HeaderData header;
    header.m_real["a_real"] = x;
    header.writeToLocation( testFile.groupID() );
  }
  testFile.popGroup();
  testFile.popGroup();
  {
    HDF5HeaderData header;
    header.m_string["a_CXXstring"] = std::string("howdydoody");
    header.writeToLocation( testFile.groupID() );
  }
  testFile.close();
  CH_assert(!testFile.isOpen());
  //=================================================================

  HDF5HeaderData readData;

  error = testFile.open("attributeTest.h5", HDF5Handle::OPEN_RDONLY);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "File open failed "<<error<<endl;
      return error;
    }
 CH_assert(testFile.isOpen());

  error = readData.readFromFile(testFile);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "attribute read failed "<<error<<endl;
      return error;
    }

  error = compare(set1, readData);
  if (error != 0)
    {
      if ( verbose )
        {
          pout() << indent2 << "attribute read does not match write in root"<<error<<"\n";
          pout() << set1;
          pout() <<"readData:\n";
          pout() <<readData<<endl;
        }
      return error;
    }

  testFile.setGroupToLevel(2);

  readData.clear();

  if (readData.m_int.find("currentStep") != readData.m_int.end())
    {
      if ( verbose )
        pout() << indent2 << "HDF5HeaderData.clear() not functioning correct"<<endl;
      return 1;
    }

  testFile.close();

#endif  // CH_USE_HDF5
  return 0;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
  {
    if ( argv[i][0] == '-' ) //if it is an option
    {
      // compare 3 chars to differentiate -x from -xx
      if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
      {
        verbose = true ;
        // argv[i] = "" ;
      }
      else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
        // argv[i] = "" ;
      }
      else
      {
        break ;
      }
    }
  }
  return ;
}

#ifdef CH_USE_HDF5
/** returns 0 if they are the same
            1 if an attribute exists in lhs that does not exist in rhs
            2 if an attribute exists in rhs that does not exist in lhs
            3 if a common attribute has a different value
*/
int compare(const HDF5HeaderData& lhs, const HDF5HeaderData& rhs)
{

#define LEFTVSRIGHT(left, right,  mapName, Ttype, code)                          \
  for (map<std::string, Ttype>::const_iterator p = left.mapName.begin();         \
      p!= left.mapName.end(); ++p)                                               \
    {                                                                            \
      map<std::string, Ttype>::const_iterator  q = right.mapName.find(p->first); \
      if (q == right.mapName.end())                                              \
      {                                                                          \
         return code;                                                            \
      }                                                                          \
      if (q->second != p->second) return 3;                                      \
    }

  LEFTVSRIGHT(lhs, rhs, m_real, Real, 1);
  LEFTVSRIGHT(lhs, rhs, m_int, int, 1);
  LEFTVSRIGHT(lhs, rhs, m_string, std::string, 1);
  LEFTVSRIGHT(lhs, rhs, m_intvect, IntVect, 1);
  LEFTVSRIGHT(lhs, rhs, m_box, Box, 1);

  LEFTVSRIGHT(rhs, lhs, m_real, Real, 2);
  LEFTVSRIGHT(rhs, lhs, m_int, int, 2);
  LEFTVSRIGHT(rhs, lhs, m_string, std::string, 2);
  LEFTVSRIGHT(rhs, lhs, m_intvect, IntVect, 2);
  LEFTVSRIGHT(rhs, lhs, m_box, Box, 2);

  return 0;
}
#endif // CH_USE_HDF5
