/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Miscellaneous random number generator functions.

  Some useful links:

  - General discussion on pseudo-random-number generator:
    http://stackoverflow.com/questions/1167253/implementation-of-rand

  - A C Library for PRNG
    http://www.cs.wm.edu/~va/software/park/park.html

  - Implementation of the C rand function:
    http://www.mscs.dal.ca/~selinger/random/

  - GSL random numbers:
    https://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html

  - Mersenne-Twister algorithms in C and C++:
    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/c-lang.html
    http://www.mcs.anl.gov/~kazutomo/hugepage-old/twister.c
  

  \author A. Mignone (mignone@ph.unito.it)

  \date   April 21, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static long int rnd_seq_size=-1;
static long int rnd_seq_seed;

/* ********************************************************************* */
double RandomNumber (double rmin, double rmax)
/*!
 * Generate and return a random number between [rmin, rmax] using
 * the PLUTO random number generator. 
 *
 * Sequence must already have been seeded using 
 * RandomSeed() function (see below).
 *
 * Generate a random number between a and b:
 *
 * \code
 *    static int first_call = 1;
 *
 *    if (first_call){
 *      RandomSeed(seed,offset);
 *      first_call = 0;
 *    }
 *   
 *    for (k = 0; k < N; k++) {
 *      print ("rand = %f\n", RandomNumber(a,b));
 *    }
 * \endcode
 *
 * In parallel, you must ensure different sequences are generated on
 * different processors.
 * - A first possibility is to use the same seed and then to offset
 *   the sequence:
 *   \code
 *     Random_Seed(time(NULL), 32768);
 *   \endcode
 *   Each processor then generates the the same sequence but
 *   the starting point is offsetted by 32768*prank prns ahead.
 *
 * - A second possibility is to seed a different sequence on each
 *   processor and use the same offset:
 *   \code
 *     Random_Seed(SeedGenerator(), 0);
 *   \endcode
 *
 *********************************************************************** */
{
  static long int count = 0;
  double rnd;

/* -------------------------------------------------
   0. Check we're not exceeding prng period
      when sequence size has been specified.
   ------------------------------------------------- */

  count++;
  if (count >= rnd_seq_size){
    print ("! RandomNumber(): maximum number of prns [%ld] exceeded\n", rnd_seq_size);
    print ("                  count = %ld, rnd_seq_size = %ld\n",
                              count, rnd_seq_size);
    QUIT_PLUTO(1); 
  }

#if PRNG == PRNG_DEFAULT
  rnd = drand48();
#elif PRNG == PRNG_ECUYER
  rnd = NR_ran2(&rnd_seq_seed);
#elif PRNG == PRNG_MT
  rnd = genrand64_real1();
#endif

  return rmin + rnd*(rmax - rmin);
}

/* ********************************************************************* */
void RandomSeed (long int seed, long int offset)
/*!
 * Seed the random number generator.
 * Important: the seed must be called only once
 *            per sequence.
 *
 * \param [in] seed    long integer used to seed the random
 *                     number generator
 * \param [in] offset  long integer specifying the offset where 
 *                     to start the sequence (only when offset > 0).
 *                     Otherwise (offset <= 0) no offset is used.
 * 
 *
 * In parallel, each processor uses a different sequence obtained by
 * offsetting the initial one by offset*prank elements.
 * For example, if offset = 8 and 2 processors are used, a total of
 * 16 numbers (8 per processor) are generated.
 * Each processor starts at a different point in the sequence:
 *
 * \verbatim
 
   r1  r2  r3  r4  r5  r6  r7  r8  r9  r10  r11  r12  r13  r14  r15  r16
 
   <---------- rank 0 ---------->  <-------------- rank 1 -------------->
 
   \endverbatim
 *
 * Example:
 *
 * \code
 *
 *   static int first_time = 1;
 *
 *   //  Generate a sequence of 16384 prns (for each thread) 
 *   if (first_time) {
 *     RandomSeed (-34, 16384); 
 *     first_time = 0;                 
 *   }
 * 
 *   r1 = RandomNumber();
 *   .
 *   .
 *   .
 *   rn = RandomNumber();
 *
 * \endcode
 * 
 *********************************************************************** */
{
  int  nprocs = 1;
  long int i;
  double scrh;

#ifdef PARALLEL  
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif

/* ----------------------------------------------------
   1. Seed sequence
   ---------------------------------------------------- */

#if PRNG == PRNG_DEFAULT
  print ("> RandomSeed(): seeding sequence");
  print ("  [PRNG_DEFAULT; seed = %ld, offset = %ld]\n", seed, offset);
  rnd_seq_seed = seed;
  rnd_seq_size = (long)pow(2.0,48);  /* Period of 48-bit prng */
  srand48(rnd_seq_seed);
#elif PRNG == PRNG_ECUYER
  print ("> RandomSeed(): seeding sequence");
  print ("  [PRNG_ECUYER; seed = %ld, offset = %ld]\n", seed, offset);
  rnd_seq_seed = -fabs(seed);
  rnd_seq_size = (long)2.e18;  /* Approximate period of L'Ecuyer prng */
  NR_ran2(&rnd_seq_seed);
#elif PRNG == PRNG_MT
  print ("> RandomSeed(): seeding sequence");
  print ("  [PRNG_MT; seed = %ld, offset = %ld]\n", seed, offset);
  rnd_seq_seed = seed;
  rnd_seq_size = (long)pow(2,19937);  /* Mersenne-Twister prng period */ 
  init_genrand64(rnd_seq_seed);
#else
  print ("! RandomSeed(): unknown PRNG\n");
  QUIT_PLUTO(1);
#endif

/* -------------------------------------------------------
   2. Skip ahead in the sequence by offset*prank elements
      only if offset is > 0.
      Use direct call rather than RandomNumber() since
      counter must not be incremented yet.
   ------------------------------------------------------- */

  if (offset > 0){
    rnd_seq_size = offset;
    for (i = 0L; i < offset*prank; i++) {
      #if PRNG == PRNG_DEFAULT
      scrh = drand48();
      #elif PRNG == PRNG_ECUYER
      scrh = NR_ran2(&rnd_seq_seed);
      #elif PRNG == PRNG_MT
      scrh = genrand64_real1();
      #endif
    }
  }

}

/* ********************************************************************* */
double GaussianRandomNumber(double mu, double sigma)
/*!
 * Generate random deviates with a normal Gaussian distribution with
 * mean mu and standard deviation sigma.
 * (Adapted from Numerical Recipe, see "Normal (Gaussian) Deviates in
 *  Chap. 7)
 *********************************************************************** */
{
  static int iset = 0;
  double u1, u2, rsq, f;
  static double gset;

  if (iset == 0){
    do {
      u1 = 2.0*RandomNumber(0.0,1.0) - 1.0;
      u2 = 2.0*RandomNumber(0.0,1.0) - 1.0;
      rsq = u1*u1 + u2*u2;
    } while (rsq >= 1.0 || rsq == 0.0);
    f = sigma*sqrt(-2.0*log(rsq)/rsq);
    gset = u1*f;
    iset = 1;
    return u2*f + mu;
  } else{
    iset = 0;
    return gset + mu;
  }

}

/* ********************************************************************* */
double PowerLawRandomNumber(double xmin, double xmax, double n)
/*!
 * Generate random deviates with a power-law distribution.
 * Here xmin and xmax are the distribution range and n is the power-law
 * index (n != -1)
 *
 *********************************************************************** */
{
  double r, x;

  if (fabs(n+1.0) < 1.e-12) {
    print ("! PowerLawRandomNumber(): n = -1 not allowed\n");
    QUIT_PLUTO(1);
  }
  r = RandomNumber(0,1);
  x = (pow(xmax,n+1) - pow(xmin,n+1))*r + pow(xmin,n+1);
  x = pow(x,1.0/(n+1));
  return x;
}

/* ********************************************************************* */
unsigned int SeedGenerator(void )
/*
 *  Seed random sequence on multi-core system.
 *
 * \b Reference
 *    - "Random Numbers in Scientific Computing: An Introduction" \n
 *       Katzgrabber  (http://arxiv.org/abs/1005.4117), Sec. 7.1 
 *
 *
 *********************************************************************** */
{
  long s, seed,pid;

  pid  = prank;
  s    = time (NULL); /* get CPU seconds since 01/01/1970 */
  seed = abs(((s*181)*((pid-83)*359))%104729);
  return seed;
}



#if PRNG == PRNG_ECUYER

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#ifndef DBL_EPSILON
  #define DBL_EPSILON 2.2204460492503131e-16
#endif
//#define EPS 1.2e-7
#define RNMX (1.0-DBL_EPSILON)
/* ********************************************************************* */
double NR_ran2(long int *idum)
/*!
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer with
 * Bays-Durham shuffle and added safeguards.
 * Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
 * the endpoint values). Call with idum a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a sequence.
 * RNMX should approximate the largest floating value that is less than 1.
 *
 * Example:
 * \code
 *   
 *   if (first_time){
 *     long dummy=-34
 *     NR_ran2(&dummy);
 *   }
 *
 *   for (k = 0; k < N; k++){
 *     print ("rand = %10.3e\n",NR_ran2(&dummy));
 *   }  
 * \endcode
 * 
 * \b Reference:
 *    - "Numerical Recipes in C", Chap 7
 *********************************************************************** */
{
  int  j;
  long int k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0) {             /* Initialize.  */
    if (-(*idum) < 1) *idum=1;  /* Be sure to prevent idum = 0. */
    else              *idum = -(*idum);

    idum2=(*idum);
    for (j = NTAB+7; j >= 0; j--){ /* Load the shuffle table (after 8 warm-ups). */

      k     = (*idum)/IQ1;
      *idum = IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k     = (*idum)/IQ1;              /* Start here when not initializing. */
  *idum = IA1*(*idum-k*IQ1)-k*IR1;  /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;      /* overflows by SchrageÕs method.        */
  k     = idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2;  /* Compute idum2=(IA2*idum) % IM2 likewise.*/

  if (idum2 < 0) idum2 += IM2;
  j     = iy/NDIV;        /* Will be in the range 0..NTAB-1. */
  iy    = iv[j]-idum2;    /* Here idum is shuffled, idum and idum2  */
  iv[j] = *idum;          /* are combined to generate output.       */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* Because users donÕt expect
                                           endpoint values.*/
  else return temp;
}
#endif /* PRNG == PRNG_ECUYER */


#if PRNG == PRNG_MT
/* ===================================================================== 
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)  
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

 
   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and 
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and 
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)

   =====================================================================

o  Taken from
  
  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
  ==> MT fo 64-bit machines ==>  mt19937-64.c 
  
o For some explanation see
  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/readme-mt.txt
  
o  For parallel/dynamic creation:
  
  ==> http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/DC/dc.html
  
o Alternatively, the PCG:
  http://www.pcg-random.org/download.html

*/
#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */

static unsigned long long mt[NN]; /* The array for the state vector */
static int mti=NN+1;              /* mti==NN+1 means mt[NN] is not initialized */

/* initializes mt[NN] with a seed */
void init_genrand64(unsigned long long seed)
{
    mt[0] = seed;
    for (mti=1; mti<NN; mti++) 
        mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long genrand64_int64(void)
{
  int i;
  unsigned long long x;
  static unsigned long long mag01[2]={0ULL, MATRIX_A};

  if (mti >= NN) { /* generate NN words at one time */

    /* if init_genrand64() has not been called, */
    /* a default initial seed is used     */
    if (mti == NN+1) init_genrand64(5489ULL); 

    for (i=0;i<NN-MM;i++) {
      x = (mt[i]&UM)|(mt[i+1]&LM);
      mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
    }
    for (;i<NN-1;i++) {
      x = (mt[i]&UM)|(mt[i+1]&LM);
      mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
    }
    x = (mt[NN-1]&UM)|(mt[0]&LM);
    mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

    mti = 0;
  }
  
  x = mt[mti++];

  x ^= (x >> 29) & 0x5555555555555555ULL;
  x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
  x ^= (x << 37) & 0xFFF7EEE000000000ULL;
  x ^= (x >> 43);

  return x;
}

/* generates a random number on [0,1]-real-interval */
double genrand64_real1(void)
{
  return (genrand64_int64() >> 11) * (1.0/9007199254740991.0);
}

/*
int main(void)
{
    int i;
    init_genrand64(3175671);

    printf("\n20 outputs of genrand64_real2()\n");
    for (i=0; i<2000; i++) {
      printf("%10.8f\n", genrand64_real1());
    }
    return 0;
}
*/
#endif /* PRNG == PRNG_MT */
