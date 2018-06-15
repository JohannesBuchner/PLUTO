#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include <cstring>

#include "REAL.H"
#include "RefCountedPtr.H"
#include "parstream.H"
#include "Vector.H"
#include "CH_Timer.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingBaseNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testRCP();

/// Global variables for handling output:
static const char *pgmname = "testRefCountedPtr" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
  CH_TIMERS("mymain");
  CH_TIMER("wholething", t1);
  CH_START(t1);

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;


  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testRCP() ;
#ifdef CH_MPI
  MPI_Finalize();
#endif

  double timeForCode=0;
  CH_STOPV(t1, timeForCode);

  return ret;
}

class Bar
{
public:
  virtual ~Bar()
  {
  }

  int cherry;
};

class Foo : public Bar
{
public:
  virtual ~Foo()
  {
  }

  int tangerine;
};


RefCountedPtr<Bar> returnByValueFunc(RefCountedPtr<Bar> input)
{
  return input;
}

void breakpointHook()
{
}

int
testRCP()
{
  CH_TIME("testRPC");
  int return_code = 0;

  RefCountedPtr<Vector<Foo> > rv(new Vector<Foo>(5));
  RefCountedPtr<Vector<Foo> > rv2;
  if (rv2.isNonUnique() || rv.isNonUnique())
  {
    if (verbose)
    {
      pout() << indent2 << "RefCountedPtr failed isNonUnique" << endl;
    }
    return_code = -1;
  }

  rv2=rv;
  if (!rv.isNonUnique() && !rv2.isNonUnique())
  {
    if (verbose)
    {
      pout() << indent2 << "RefCountedPtr failed isNonUnique" << endl;
    }
    return_code = -1;
  }

  {
    RefCountedPtr<Vector<Foo> > rv3;
    rv3=rv2;
    if (rv3.refCount() != 3)
      {
        if (verbose)
          {
            pout() << indent2 << "RefCountedPtr::refCount is wrong" << endl;
          }
        return_code = -1;
      }
  } //rv3 now out of scope
  if (rv.refCount() != 2)
    {
      if (verbose)
        {
          pout() << indent2 << "RefCountedPtr::refCount is wrong" << endl;
        }
      return_code = -1;
    }

  RefCountedPtr<Foo> myFoo(new Foo);
  myFoo->cherry = 5;
  RefCountedPtr<Bar> myBar(new Bar), myBar2, myBar3;
  myBar->cherry = 3;
  myBar2 = myBar;
  myBar = myFoo;
  myBar3 = myBar;
  if (myBar->cherry != 5 || myBar2.refCount() != 1 || myFoo.refCount() != 3 || myBar.refCount() != 3)
    {
      if (verbose)
        {
          pout() << indent2 << "RefCountedPtr template conversion is wrong" << endl;
        }
      return_code = -1;
    }

  myBar =  returnByValueFunc(myFoo);
  if (myBar->cherry != 5 || myBar2.refCount() != 1 || myFoo.refCount() != 3 || myBar.refCount() != 3)
    {
      if (verbose)
        {
          pout() << indent2 << "RefCountedPtr template conversion is wrong" << endl;
        }
      return_code = -1;
    }

//
  if (return_code == 0)
  {
    pout() << indent << "testRefCountedPtr test passed." << endl;
  }
  else
  {
    pout() << indent << "testRefCountedPtr test FAILED!!!" << endl;
  }

  breakpointHook();

  return return_code;
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
