#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Test program for AMR/LoadBalance
// The program simulates different processor counts so
// it should not be run under MPI.

// Test 1: single level, small number of loads.
// Test 2: single level, large number of loads.
// Test 3: single level, loads from Brian's HDF test code

#include <limits.h>

#include <iostream>

#include "SPMD.H"  // to get num_procs global variable
#include "LoadBalance.H"
#include "Misc.H"
#include "parstream.H"
#include "UsingNamespace.H"

/// Prototypes:
int
testLB1(void);
int
testLB2(void);
int
testLB3(void);
int
testLB4(void);

using std::endl;

void
parseTestOptions(int argc ,char* argv[]) ;

/// Global variables for handling output:
static const char* pgmname = "testLoadBalance" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;

int
main(int argc ,char* argv[])
{
  int stat_all = 0 ;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << endl ;

  int status = testLB1();

  if ( status == 0 )
  {
    if ( verbose ) pout() << indent << pgmname << " passed test 1." << endl ;
  }
  else
  {
    pout() << indent << pgmname << " failed test 1 with return code " << status << endl ;
    stat_all = status ;
  }

  status = testLB2();

  if ( status == 0 )
  {
    if ( verbose ) pout() << indent << pgmname << " passed test 2." << endl ;
  }
  else
  {
    pout() << indent << pgmname << " failed test 2 with return code " << status << endl ;
    stat_all = status ;
  }

  status = testLB3();

  if ( status == 0 )
  {
    if ( verbose ) pout() << indent << pgmname << " passed test 3." << endl ;
  }
  else
  {
    pout() << indent << pgmname << " failed test 3 with return code " << status << endl ;
    stat_all = status ;
  }

  status = testLB4();

  if ( status == 0 )
  {
    if ( verbose ) pout() << indent << pgmname << " passed test 4." << endl ;
  }
  else
  {
    pout() << indent << pgmname << " failed test 4 with return code " << status << endl ;
    stat_all = status ;
  }

  if ( stat_all == 0 )
    pout() << indent << pgmname << ": passed all tests." << endl ;
  else
    pout() << indent << pgmname << ": failed one or more tests." << endl ;

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return stat_all ;
}

int
testLB1()
{
  const int numLoads = 9 ;
  Vector<Vector<long> > loads( 1 ) ;
  Real total_load = 0.0 ;
  for ( int i=0 ; i<numLoads ; ++i )
    {
      loads[0].push_back( (i+1) ) ;
      total_load += (i+1) ;
    }
  Vector<Vector<Box> > grids( 1 ) ;
  grids[0].resize( numLoads ) ;
  Vector<int> refratios( 1,2 ) ;
  Vector<Vector<int> > assignments( 1 ) ;
  Real eff_ratio ;
  int status ;

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios,1 ) ;
  if ( status != 0 ) return status ;
  if ( verbose )
  {
    pout() << indent2 << "test1: 1 processors, total load " << total_load
         << ", efficiency " << (int)eff_ratio*100.0 << "%, assignments are" << endl ;
    pout() << indent2 ;
    for ( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  if ( eff_ratio != 1.0 )
    {
      if ( verbose )
        pout() << indent2 << "test1: 1 processor, expected efficiency 100%, got "
             << eff_ratio*100 << endl ;
      return -11 ;
    }
  for ( int i=0 ; i<numLoads ; ++i )
    {
      if ( assignments[0][i] != 0 )
        {
          if ( verbose )
            pout() << indent2 << "test1: 1 processor, expected assignment 0, got "
                 << assignments[0][i] << ", for load " << i << endl ;
          return -12 ;
        }
    }

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios, 2 ) ;
  if ( status != 0 ) return status ;
  if ( verbose )
  {
    pout() << indent2 << "test1: 2 processors, total load " << total_load
         << ", efficiency " << (int)(eff_ratio*100.0) << "%, assignments are" << endl ;
    pout() << indent2 ;
    for ( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  if ( eff_ratio < 0.95 )
    {
      if ( verbose )
        pout() << indent2 << "test1: 2 processor, expected efficiency >95%, got "
             << eff_ratio*100 << endl ;
      return -21 ;
    }
  for ( int i=0 ; i<numLoads ; ++i )
    {
      if ( assignments[0][i] < 0 || assignments[0][i] > 1 )
        {
          if ( verbose )
            pout() << indent2 << "test1: 2 processor, got assignment " << assignments[0][i] << ", for load " << i << endl ;
          return -22 ;
        }
    }

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios, 4 ) ;
  if ( verbose )
  {
    pout() << indent2 << "test1: 4 processors, total load " << total_load
         << ", efficiency " << (int)(eff_ratio*100.0) << "%, assignments are" << endl ;
    pout() << indent2 ;
    for ( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  if ( eff_ratio < 0.80 )
    {
      if ( verbose )
        pout() << indent2 << "test1: 4 processor, expected efficiency >80%, got "
             << eff_ratio*100 << endl ;
      return -41 ;
    }
  for ( int i=0 ; i<numLoads ; ++i )
    {
      if ( assignments[0][i] < 0 || assignments[0][i] > 3 )
        {
          if ( verbose )
            pout() << indent2 << "test1: 4 processor: got assignment " << assignments[0][i] << ", for load " << i << endl ;
          return -42 ;
        }
    }

  return status ;
}

int
testLB2()
{
  const int numLoads = 8 ;
  Vector<Vector<long> > loads( 1 ) ;
  Real total_load = 0.0 ;
  loads[0].push_back( 7 ) ; total_load += 7 ;
  loads[0].push_back( 6 ) ; total_load += 6 ;
  loads[0].push_back( 5 ) ; total_load += 5 ;
  loads[0].push_back( 4 ) ; total_load += 4 ;
  loads[0].push_back( 4 ) ; total_load += 4 ;
  loads[0].push_back( 3 ) ; total_load += 3 ;
  loads[0].push_back( 3 ) ; total_load += 3 ;
  loads[0].push_back( 2 ) ; total_load += 2 ;

  Vector<Vector<Box> > grids( 1 ) ;
  grids[0].resize( numLoads ) ;
  Vector<int> refratios( 1,2 ) ;
  Vector<Vector<int> > assignments( 1 ) ;
  Real eff_ratio ;
  int status ;

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios, 3 ) ;
  if ( status != 0 ) return status ;
  if ( verbose )
  {
    pout() << indent2 << "test2: 3 processors, total load " << total_load
         << ", efficiency " << (int)(eff_ratio*100.0) << "%, assignments are" << endl ;
    pout() << indent2 ;
    for ( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  Real correct_eff = (Real)11 / (Real)12 ;
  if ( Abs(eff_ratio - correct_eff) > (Real)1e-5 )
    {
      if ( verbose )
        pout() << indent2 << "test2: 3 processor, expected efficiency " << correct_eff
             << ", got " << eff_ratio*100 << endl ;
      return -221 ;
    }
  for ( int i=0 ; i<numLoads ; ++i )
    {
      if ( assignments[0][i] < 0 || assignments[0][i] > 2 )
        {
          if ( verbose )
            pout() << indent2 << "test2: 3 processor, assignment " << assignments[0][i] << ", for load " << i
                 << " is out of range" << endl ;
          return -222 ;
        }
    }

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios, 4 ) ;
  if ( verbose )
  {
    pout() << indent2 << "test2: 4 processors, total load " << total_load
         << ", efficiency " << (int)(eff_ratio*100.0) << "%, assignments are" << endl ;
    pout() << indent2 ;
    for ( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  if ( eff_ratio < 0.80 )
    {
      if ( verbose )
        pout() << indent2 << "test2: 4 processor, expected efficiency >80%, got "
             << eff_ratio*100 << endl ;
      return -241 ;
    }
  for ( int i=0 ; i<numLoads ; ++i )
    {
      if ( assignments[0][i] < 0 || assignments[0][i] > 3 )
        {
          if ( verbose )
            pout() << indent2 << "test2: 4 processor: got assignment " << assignments[0][i] << ", for load " << i << endl ;
          return -242 ;
        }
    }

  return status ;
}

int
testLB3()
{
  const int numLoads = 16 ;
  Vector<Vector<long> > loads( 1 ) ;
  Real total_load = 0.0 ;
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((4,10) (4,10) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((5,7) (5,7) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((5,13) (5,13) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((6,6) (6,6) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((6,14) (6,14) (0,0))
  loads[0].push_back( 3 ) ; total_load += 3 ; //box ((7,5) (9,5) (0,0))
  loads[0].push_back( 3 ) ; total_load += 3 ; //box ((7,15) (9,15) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((10,4) (10,4) (0,0))
  loads[0].push_back( 8 ) ; total_load += 8 ; //box ((10,15) (13,16) (0,0))
  loads[0].push_back( 2 ) ; total_load += 2 ; //box ((11,5) (12,5) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((13,5) (13,5) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((14,6) (14,6) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((14,14) (14,14) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((15,7) (15,7) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((15,13) (15,13) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((16,10) (16,10) (0,0))

  Vector<Vector<Box> > grids( 1 ) ;
  grids[0].resize( numLoads ) ;
  Vector<int> refratios( 1,2 ) ;
  Vector<Vector<int> > assignments( 1 ) ;
  Real eff_ratio ;
  int status ;

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios, 3 ) ;
  if ( status != 0 ) return status ;
  if ( verbose )
  {
    pout() << indent2 << "test3: 3 processors, total load " << total_load
         << ", efficiency " << (int)(eff_ratio*100.0) << "%, assignments are" << endl ;
    pout() << indent2 ;
    for ( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  Real correct_eff = (Real)9 / (Real)10 ;
  if ( Abs(eff_ratio - correct_eff) > (Real)1e-5 )
    {
      if ( verbose )
        pout() << indent2 << "test3: 3 processor, expected efficiency " << correct_eff
             << ", got " << eff_ratio*100 << endl ;
      return -321 ;
    }
  for ( int i=0 ; i<numLoads ; ++i )
    {
      if ( assignments[0][i] < 0 || assignments[0][i] > 3 )
        {
          if ( verbose )
            pout() << indent2 << "test3: 3 processor, assignment " << assignments[0][i] << ", for load " << i
                 << " is out of range" << endl ;
          return -322 ;
        }
    }

  return 0 ;
}
int
testLB4()
{
  const int numBoxes = 483 ;

  Vector<long> loads(numBoxes);

  //the following data is from a real run
  loads[  0]=     13824 ;
  loads[  1]=     24576 ;
  loads[  2]=       512 ;
  loads[  3]=      4096 ;
  loads[  4]=     18432 ;
  loads[  5]=     32768 ;
  loads[  6]=     32768 ;
  loads[  7]=     32768 ;
  loads[  8]=     12288 ;
  loads[  9]=     18432 ;
  loads[ 10]=     16384 ;
  loads[ 11]=     32768 ;
  loads[ 12]=     32768 ;
  loads[ 13]=      8192 ;
  loads[ 14]=      8192 ;
  loads[ 15]=     12288 ;
  loads[ 16]=      8192 ;
  loads[ 17]=      9216 ;
  loads[ 18]=     13824 ;
  loads[ 19]=      9216 ;
  loads[ 20]=      8192 ;
  loads[ 21]=     12288 ;
  loads[ 22]=      8192 ;
  loads[ 23]=     12288 ;
  loads[ 24]=      8192 ;
  loads[ 25]=     12288 ;
  loads[ 26]=     13824 ;
  loads[ 27]=      9216 ;
  loads[ 28]=     13824 ;
  loads[ 29]=     12288 ;
  loads[ 30]=      8192 ;
  loads[ 31]=     12288 ;
  loads[ 32]=      8192 ;
  loads[ 33]=     12288 ;
  loads[ 34]=      8192 ;
  loads[ 35]=     12288 ;
  loads[ 36]=      9216 ;
  loads[ 37]=     13824 ;
  loads[ 38]=      9216 ;
  loads[ 39]=     13824 ;
  loads[ 40]=      8192 ;
  loads[ 41]=     12288 ;
  loads[ 42]=      8192 ;
  loads[ 43]=     12288 ;
  loads[ 44]=      8192 ;
  loads[ 45]=     12288 ;
  loads[ 46]=      8192 ;
  loads[ 47]=      9216 ;
  loads[ 48]=     13824 ;
  loads[ 49]=      9216 ;
  loads[ 50]=      8192 ;
  loads[ 51]=     12288 ;
  loads[ 52]=      8192 ;
  loads[ 53]=     12288 ;
  loads[ 54]=     18432 ;
  loads[ 55]=     12288 ;
  loads[ 56]=     12288 ;
  loads[ 57]=     18432 ;
  loads[ 58]=     12288 ;
  loads[ 59]=     12288 ;
  loads[ 60]=      8192 ;
  loads[ 61]=     12288 ;
  loads[ 62]=     13824 ;
  loads[ 63]=      9216 ;
  loads[ 64]=     13824 ;
  loads[ 65]=     12288 ;
  loads[ 66]=     24576 ;
  loads[ 67]=      8192 ;
  loads[ 68]=     12288 ;
  loads[ 69]=     12288 ;
  loads[ 70]=     18432 ;
  loads[ 71]=     12288 ;
  loads[ 72]=      9216 ;
  loads[ 73]=     13824 ;
  loads[ 74]=      9216 ;
  loads[ 75]=     12288 ;
  loads[ 76]=     18432 ;
  loads[ 77]=     12288 ;
  loads[ 78]=     18432 ;
  loads[ 79]=     12288 ;
  loads[ 80]=     18432 ;
  loads[ 81]=     18432 ;
  loads[ 82]=     12288 ;
  loads[ 83]=     18432 ;
  loads[ 84]=     12288 ;
  loads[ 85]=     18432 ;
  loads[ 86]=     12288 ;
  loads[ 87]=     18432 ;
  loads[ 88]=     12288 ;
  loads[ 89]=     18432 ;
  loads[ 90]=     12288 ;
  loads[ 91]=     18432 ;
  loads[ 92]=     18432 ;
  loads[ 93]=     12288 ;
  loads[ 94]=     18432 ;
  loads[ 95]=     13824 ;
  loads[ 96]=      9216 ;
  loads[ 97]=     13824 ;
  loads[ 98]=     18432 ;
  loads[ 99]=     12288 ;
  loads[100]=     18432 ;
  loads[101]=     12288 ;
  loads[102]=     18432 ;
  loads[103]=     12288 ;
  loads[104]=     18432 ;
  loads[105]=      9216 ;
  loads[106]=     13824 ;
  loads[107]=      9216 ;
  loads[108]=     13824 ;
  loads[109]=     12288 ;
  loads[110]=     18432 ;
  loads[111]=     12288 ;
  loads[112]=     18432 ;
  loads[113]=     12288 ;
  loads[114]=     18432 ;
  loads[115]=     12288 ;
  loads[116]=     12288 ;
  loads[117]=     18432 ;
  loads[118]=     12288 ;
  loads[119]=     18432 ;
  loads[120]=     12288 ;
  loads[121]=     18432 ;
  loads[122]=     18432 ;
  loads[123]=     12288 ;
  loads[124]=     18432 ;
  loads[125]=     12288 ;
  loads[126]=     18432 ;
  loads[127]=     12288 ;
  loads[128]=      9216 ;
  loads[129]=     13824 ;
  loads[130]=      9216 ;
  loads[131]=     12288 ;
  loads[132]=     18432 ;
  loads[133]=     12288 ;
  loads[134]=     18432 ;
  loads[135]=     12288 ;
  loads[136]=     18432 ;
  loads[137]=     13824 ;
  loads[138]=      9216 ;
  loads[139]=     13824 ;
  loads[140]=     18432 ;
  loads[141]=     12288 ;
  loads[142]=     18432 ;
  loads[143]=     12288 ;
  loads[144]=     18432 ;
  loads[145]=     12288 ;
  loads[146]=     18432 ;
  loads[147]=     13824 ;
  loads[148]=     12288 ;
  loads[149]=     18432 ;
  loads[150]=     12288 ;
  loads[151]=     18432 ;
  loads[152]=     12288 ;
  loads[153]=     18432 ;
  loads[154]=     18432 ;
  loads[155]=     12288 ;
  loads[156]=     12288 ;
  loads[157]=     12288 ;
  loads[158]=     12288 ;
  loads[159]=      2048 ;
  loads[160]=      2048 ;
  loads[161]=     12288 ;
  loads[162]=      8192 ;
  loads[163]=      1024 ;
  loads[164]=     12288 ;
  loads[165]=      9216 ;
  loads[166]=     13824 ;
  loads[167]=     12288 ;
  loads[168]=     12288 ;
  loads[169]=      9216 ;
  loads[170]=     12288 ;
  loads[171]=     12288 ;
  loads[172]=     12288 ;
  loads[173]=     12288 ;
  loads[174]=     13824 ;
  loads[175]=      9216 ;
  loads[176]=     12288 ;
  loads[177]=     13824 ;
  loads[178]=     16384 ;
  loads[179]=     12288 ;
  loads[180]=     16384 ;
  loads[181]=      9216 ;
  loads[182]=     13824 ;
  loads[183]=     12288 ;
  loads[184]=     12288 ;
  loads[185]=      9216 ;
  loads[186]=     13824 ;
  loads[187]=     12288 ;
  loads[188]=     12288 ;
  loads[189]=     12288 ;
  loads[190]=     12288 ;
  loads[191]=      9216 ;
  loads[192]=     12288 ;
  loads[193]=     13824 ;
  loads[194]=      9216 ;
  loads[195]=     12288 ;
  loads[196]=     12288 ;
  loads[197]=     12288 ;
  loads[198]=     12288 ;
  loads[199]=     12288 ;
  loads[200]=     13824 ;
  loads[201]=     12288 ;
  loads[202]=      9216 ;
  loads[203]=     13824 ;
  loads[204]=     16384 ;
  loads[205]=     12288 ;
  loads[206]=     16384 ;
  loads[207]=     24576 ;
  loads[208]=      9216 ;
  loads[209]=     13824 ;
  loads[210]=     12288 ;
  loads[211]=     18432 ;
  loads[212]=      9216 ;
  loads[213]=     12288 ;
  loads[214]=     12288 ;
  loads[215]=     18432 ;
  loads[216]=     12288 ;
  loads[217]=     12288 ;
  loads[218]=     18432 ;
  loads[219]=     12288 ;
  loads[220]=     12288 ;
  loads[221]=     18432 ;
  loads[222]=     12288 ;
  loads[223]=     13824 ;
  loads[224]=      9216 ;
  loads[225]=     18432 ;
  loads[226]=     12288 ;
  loads[227]=     13824 ;
  loads[228]=     18432 ;
  loads[229]=     18432 ;
  loads[230]=     12288 ;
  loads[231]=     18432 ;
  loads[232]=     18432 ;
  loads[233]=     12288 ;
  loads[234]=     18432 ;
  loads[235]=     18432 ;
  loads[236]=     12288 ;
  loads[237]=     18432 ;
  loads[238]=      9216 ;
  loads[239]=     13824 ;
  loads[240]=     12288 ;
  loads[241]=     18432 ;
  loads[242]=      9216 ;
  loads[243]=     13824 ;
  loads[244]=     12288 ;
  loads[245]=     18432 ;
  loads[246]=     12288 ;
  loads[247]=     18432 ;
  loads[248]=     12288 ;
  loads[249]=     18432 ;
  loads[250]=     12288 ;
  loads[251]=     18432 ;
  loads[252]=     12288 ;
  loads[253]=     18432 ;
  loads[254]=     12288 ;
  loads[255]=     18432 ;
  loads[256]=     12288 ;
  loads[257]=     18432 ;
  loads[258]=      9216 ;
  loads[259]=     12288 ;
  loads[260]=     13824 ;
  loads[261]=      9216 ;
  loads[262]=     18432 ;
  loads[263]=     12288 ;
  loads[264]=     12288 ;
  loads[265]=     18432 ;
  loads[266]=     12288 ;
  loads[267]=     12288 ;
  loads[268]=     18432 ;
  loads[269]=     12288 ;
  loads[270]=     12288 ;
  loads[271]=     18432 ;
  loads[272]=     12288 ;
  loads[273]=     13824 ;
  loads[274]=     18432 ;
  loads[275]=      9216 ;
  loads[276]=     13824 ;
  loads[277]=     12288 ;
  loads[278]=     18432 ;
  loads[279]=     18432 ;
  loads[280]=     12288 ;
  loads[281]=     18432 ;
  loads[282]=     18432 ;
  loads[283]=     12288 ;
  loads[284]=     18432 ;
  loads[285]=     18432 ;
  loads[286]=     12288 ;
  loads[287]=     18432 ;
  loads[288]=     12288 ;
  loads[289]=     18432 ;
  loads[290]=     12288 ;
  loads[291]=     18432 ;
  loads[292]=     12288 ;
  loads[293]=     18432 ;
  loads[294]=      6144 ;
  loads[295]=      9216 ;
  loads[296]=     13824 ;
  loads[297]=      6144 ;
  loads[298]=      9216 ;
  loads[299]=     18432 ;
  loads[300]=      6144 ;
  loads[301]=      2048 ;
  loads[302]=      2048 ;
  loads[303]=       512 ;
  loads[304]=      9216 ;
  loads[305]=      6144 ;
  loads[306]=      6144 ;
  loads[307]=     12288 ;
  loads[308]=      9216 ;
  loads[309]=      3072 ;
  loads[310]=      9216 ;
  loads[311]=      6144 ;
  loads[312]=      6144 ;
  loads[313]=      6144 ;
  loads[314]=      1024 ;
  loads[315]=      9216 ;
  loads[316]=      1536 ;
  loads[317]=      9216 ;
  loads[318]=      1024 ;
  loads[319]=      9216 ;
  loads[320]=      6144 ;
  loads[321]=      8192 ;
  loads[322]=      1536 ;
  loads[323]=      9216 ;
  loads[324]=      1024 ;
  loads[325]=      1536 ;
  loads[326]=     12288 ;
  loads[327]=      6144 ;
  loads[328]=      6144 ;
  loads[329]=      6144 ;
  loads[330]=      1024 ;
  loads[331]=      9216 ;
  loads[332]=      1536 ;
  loads[333]=      9216 ;
  loads[334]=      1024 ;
  loads[335]=      1536 ;
  loads[336]=      9216 ;
  loads[337]=      6144 ;
  loads[338]=      6144 ;
  loads[339]=      6144 ;
  loads[340]=      9216 ;
  loads[341]=      1024 ;
  loads[342]=      1536 ;
  loads[343]=      9216 ;
  loads[344]=      1024 ;
  loads[345]=      9216 ;
  loads[346]=      6144 ;
  loads[347]=      8192 ;
  loads[348]=      1536 ;
  loads[349]=      9216 ;
  loads[350]=      1024 ;
  loads[351]=     12288 ;
  loads[352]=      1536 ;
  loads[353]=      4096 ;
  loads[354]=      6144 ;
  loads[355]=      4096 ;
  loads[356]=      9216 ;
  loads[357]=     13824 ;
  loads[358]=      9216 ;
  loads[359]=      6144 ;
  loads[360]=      4096 ;
  loads[361]=      6144 ;
  loads[362]=     13824 ;
  loads[363]=      9216 ;
  loads[364]=     13824 ;
  loads[365]=      4096 ;
  loads[366]=      6144 ;
  loads[367]=      4096 ;
  loads[368]=      6144 ;
  loads[369]=      9216 ;
  loads[370]=     13824 ;
  loads[371]=      9216 ;
  loads[372]=     13824 ;
  loads[373]=      4096 ;
  loads[374]=      6144 ;
  loads[375]=      4096 ;
  loads[376]=      9216 ;
  loads[377]=     13824 ;
  loads[378]=      9216 ;
  loads[379]=      6144 ;
  loads[380]=      4096 ;
  loads[381]=      6144 ;
  loads[382]=     13824 ;
  loads[383]=      9216 ;
  loads[384]=     13824 ;
  loads[385]=     13824 ;
  loads[386]=      3072 ;
  loads[387]=      4608 ;
  loads[388]=     16384 ;
  loads[389]=     24576 ;
  loads[390]=     18432 ;
  loads[391]=     16384 ;
  loads[392]=     24576 ;
  loads[393]=      6144 ;
  loads[394]=      9216 ;
  loads[395]=      9216 ;
  loads[396]=      6144 ;
  loads[397]=      9216 ;
  loads[398]=      8192 ;
  loads[399]=      6144 ;
  loads[400]=      9216 ;
  loads[401]=      6144 ;
  loads[402]=      6144 ;
  loads[403]=      1024 ;
  loads[404]=      9216 ;
  loads[405]=      1536 ;
  loads[406]=      9216 ;
  loads[407]=      6144 ;
  loads[408]=      1024 ;
  loads[409]=      9216 ;
  loads[410]=      6144 ;
  loads[411]=      1536 ;
  loads[412]=      9216 ;
  loads[413]=      1024 ;
  loads[414]=      8192 ;
  loads[415]=      1536 ;
  loads[416]=     12288 ;
  loads[417]=      6144 ;
  loads[418]=      6144 ;
  loads[419]=      1024 ;
  loads[420]=      9216 ;
  loads[421]=      1536 ;
  loads[422]=      9216 ;
  loads[423]=      6144 ;
  loads[424]=      1024 ;
  loads[425]=      1536 ;
  loads[426]=      9216 ;
  loads[427]=      6144 ;
  loads[428]=      9216 ;
  loads[429]=      1024 ;
  loads[430]=      6144 ;
  loads[431]=      6144 ;
  loads[432]=      1536 ;
  loads[433]=      9216 ;
  loads[434]=      1024 ;
  loads[435]=      9216 ;
  loads[436]=      6144 ;
  loads[437]=      1536 ;
  loads[438]=      9216 ;
  loads[439]=      8192 ;
  loads[440]=      1024 ;
  loads[441]=     12288 ;
  loads[442]=      1536 ;
  loads[443]=      9216 ;
  loads[444]=     13824 ;
  loads[445]=      4096 ;
  loads[446]=      6144 ;
  loads[447]=      9216 ;
  loads[448]=      4096 ;
  loads[449]=     13824 ;
  loads[450]=      9216 ;
  loads[451]=      6144 ;
  loads[452]=      4096 ;
  loads[453]=     13824 ;
  loads[454]=      6144 ;
  loads[455]=      9216 ;
  loads[456]=     13824 ;
  loads[457]=      4096 ;
  loads[458]=      6144 ;
  loads[459]=      9216 ;
  loads[460]=     13824 ;
  loads[461]=      4096 ;
  loads[462]=      6144 ;
  loads[463]=      9216 ;
  loads[464]=      4096 ;
  loads[465]=     13824 ;
  loads[466]=      9216 ;
  loads[467]=      6144 ;
  loads[468]=      4096 ;
  loads[469]=     13824 ;
  loads[470]=      6144 ;
  loads[471]=      9216 ;
  loads[472]=     13824 ;
  loads[473]=      4096 ;
  loads[474]=      6144 ;
  loads[475]=      9216 ;
  loads[476]=     13824 ;
  loads[477]=      6144 ;
  loads[478]=      9216 ;
  loads[479]=     13824 ;
  loads[480]=      9216 ;
  loads[481]=     13824 ;
  loads[482]=     18432 ;

  int total_load = 0;
  int minLoadN = INT_MAX;
  int maxLoadN = -1;
  int minLoadBox = -1;
  int maxLoadBox = -1;
  Vector<Box> grids;
  for (int ibox = 0;ibox<numBoxes;ibox++)
    {
      grids.push_back(Box(IntVect::Zero,IntVect::Unit));
      const int thisLoad = loads[ibox];
      total_load += thisLoad;
      if (thisLoad>maxLoadN)
        {
          maxLoadN = thisLoad;
          maxLoadBox = ibox;
        }
      if (thisLoad<minLoadN)
        {
          minLoadN = thisLoad;
          minLoadBox = ibox;
        }
    }
  if (verbose)
    {
      pout() << indent2 << "test4: load dataset has: " << std::endl;
      pout() << indent2 << "numBoxes = " << numBoxes << "; total load = " << total_load << std::endl;
      pout() << indent2 << "Min box load = " << minLoadN << "; ibox = " << minLoadBox << std::endl;
      pout() << indent2 << "Max box load = " << maxLoadN << "; ibox = " << maxLoadBox << std::endl;
      pout() << std::endl;
    }

  Vector<int> assignments;
  // Real eff_ratio;
  int status;

  //  int numProcs = 384;//this was breaking LoadBalance on Jan 5, 2011
  for ( int numProcs = 12;numProcs<2*numBoxes;numProcs *= 2)
    {

      status = LoadBalance( assignments, loads, grids, numProcs ) ;
      if ( status != 0 ) return status ;
      if ( verbose )
        {
          pout() << indent2 << "test4: run with " << numProcs << " processors" << std::endl;
        }

      int minNum = INT_MAX;
      int maxNum = 0;
      int procMin = -1;
      int procMax = -1;
      int minLoad = INT_MAX;
      int maxLoad = -1;
      int minLoadProc = -1;
      int maxLoadProc = -1;
      for ( int iproc = 0;iproc<numProcs;iproc++)
        {
          //count boxes for this proc
          int thisNum = 0;
          long int thisLoad = 0;
          for ( int ibox=0 ; ibox<numBoxes ; ++ibox )
            {
              if (assignments[ibox] == iproc)
                {
                  thisLoad += loads[ibox];
                  thisNum++;
                }
            }
          if (thisNum<minNum)
            {
              minNum = thisNum;
              procMin = iproc;
            }
          if (thisNum>maxNum)
            {
              maxNum = thisNum;
              procMax = iproc;
            }
          if (thisLoad>maxLoad)
            {
              maxLoad = thisLoad;
              maxLoadProc = iproc;
            }
          if (thisLoad<minLoad)
            {
              minLoad = thisLoad;
              minLoadProc = iproc;
            }
        }
      if ( verbose )
        {
          pout() << indent2 << indent2 << "Min proc load      = " << minLoad  << "; proc = " << minLoadProc        << std::endl;
          pout() << indent2 << indent2 << "Max proc load      = " << maxLoad  << "; proc = " << maxLoadProc        << std::endl;
          pout() << indent2 << indent2 << "Min boxes per proc = " << minNum   << "; proc = " << procMin            << std::endl;
          pout() << indent2 << indent2 << "Max boxes per proc = " << maxNum   << "; proc = " << procMax            << std::endl;
          pout() << indent2 << indent2 << "Goal = total_load/numProcs    = " << (Real)total_load / (Real)numProcs  << std::endl;
          pout() << indent2 << indent2 << "Efficiency = minLoad/maxLoad  = " << (Real)minLoad    / (Real)(maxLoad) << std::endl;
          pout() << std::endl;
        }

      if (numBoxes>=numProcs && minNum==0)
        {
          if ( verbose )
            {
              pout() << indent2 << indent2 << "ERROR:: found processor without load when numBoxes>=numProcs" << std::endl;
            }
          return -421;

        }

      for ( int i=0 ; i<numBoxes ; ++i )
        {
          if ( assignments[i] < 0 || assignments[i] > numProcs )
            {
              if ( verbose )
                {
                  pout() << indent2 << "test4: " << numProcs << " processor, assignment " << assignments[i] << ", for load " << i << " is out of range" << std::endl;
                }
              return -422;
            }
        }
    }

  return 0 ;
}

///
// Parse the standard test options (-v -q) out of the command line.
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
        }
    }
  return ;
}
