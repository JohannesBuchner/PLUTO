#ifndef _AVERAGEF_F_H_
#define _AVERAGEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAVERAGE 
#define GUARDAVERAGE 
// Prototype for Fortran procedure average ...
//
void FORTRAN_NAME( AVERAGE ,average )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_BOX(bref) );

#define FORT_AVERAGE FORTRAN_NAME( inlineAVERAGE, inlineAVERAGE)
#define FORTNT_AVERAGE FORTRAN_NAME( AVERAGE, average)

inline void FORTRAN_NAME(inlineAVERAGE, inlineAVERAGE)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_AVERAGE");
 FORTRAN_NAME( AVERAGE ,average )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(refRatio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDAVERAGE 

#ifndef GUARDAVERAGEHARMONIC 
#define GUARDAVERAGEHARMONIC 
// Prototype for Fortran procedure averageharmonic ...
//
void FORTRAN_NAME( AVERAGEHARMONIC ,averageharmonic )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_BOX(bref) );

#define FORT_AVERAGEHARMONIC FORTRAN_NAME( inlineAVERAGEHARMONIC, inlineAVERAGEHARMONIC)
#define FORTNT_AVERAGEHARMONIC FORTRAN_NAME( AVERAGEHARMONIC, averageharmonic)

inline void FORTRAN_NAME(inlineAVERAGEHARMONIC, inlineAVERAGEHARMONIC)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_AVERAGEHARMONIC");
 FORTRAN_NAME( AVERAGEHARMONIC ,averageharmonic )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(refRatio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDAVERAGEHARMONIC 

#ifndef GUARDAVERAGEINTVECTREF 
#define GUARDAVERAGEINTVECTREF 
// Prototype for Fortran procedure averageintvectref ...
//
void FORTRAN_NAME( AVERAGEINTVECTREF ,averageintvectref )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INTVECT(ref)
      ,CHFp_BOX(bref) );

#define FORT_AVERAGEINTVECTREF FORTRAN_NAME( inlineAVERAGEINTVECTREF, inlineAVERAGEINTVECTREF)
#define FORTNT_AVERAGEINTVECTREF FORTRAN_NAME( AVERAGEINTVECTREF, averageintvectref)

inline void FORTRAN_NAME(inlineAVERAGEINTVECTREF, inlineAVERAGEINTVECTREF)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INTVECT(ref)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_AVERAGEINTVECTREF");
 FORTRAN_NAME( AVERAGEINTVECTREF ,averageintvectref )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INTVECT(ref)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDAVERAGEINTVECTREF 

}

#endif
