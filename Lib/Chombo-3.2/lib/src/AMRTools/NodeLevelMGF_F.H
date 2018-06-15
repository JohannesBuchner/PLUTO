#ifndef _NODELEVELMGF_F_H_
#define _NODELEVELMGF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDNODEINTERPMG_GETWEIGHTS 
#define GUARDNODEINTERPMG_GETWEIGHTS 
// Prototype for Fortran procedure NODEINTERPMG_GETWEIGHTS ...
//
void FORTRAN_NAME( NODEINTERPMG_GETWEIGHTS ,nodeinterpmg_getweights )(
      CHFp_CONST_INT(nref)
      ,CHFp_BOX(bref)
      ,CHFp_FRA(wtcrnr) );

#define FORT_NODEINTERPMG_GETWEIGHTS FORTRAN_NAME( inlineNODEINTERPMG_GETWEIGHTS, inlineNODEINTERPMG_GETWEIGHTS)
#define FORTNT_NODEINTERPMG_GETWEIGHTS FORTRAN_NAME( NODEINTERPMG_GETWEIGHTS, nodeinterpmg_getweights)

inline void FORTRAN_NAME(inlineNODEINTERPMG_GETWEIGHTS, inlineNODEINTERPMG_GETWEIGHTS)(
      CHFp_CONST_INT(nref)
      ,CHFp_BOX(bref)
      ,CHFp_FRA(wtcrnr) )
{
 CH_TIMELEAF("FORT_NODEINTERPMG_GETWEIGHTS");
 FORTRAN_NAME( NODEINTERPMG_GETWEIGHTS ,nodeinterpmg_getweights )(
      CHFt_CONST_INT(nref)
      ,CHFt_BOX(bref)
      ,CHFt_FRA(wtcrnr) );
}
#endif  // GUARDNODEINTERPMG_GETWEIGHTS 

#ifndef GUARDNODEINTERPMG 
#define GUARDNODEINTERPMG 
// Prototype for Fortran procedure NODEINTERPMG ...
//
void FORTRAN_NAME( NODEINTERPMG ,nodeinterpmg )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(crse)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(nref)
      ,CHFp_BOX(bref)
      ,CHFp_FRA(wtcrnr) );

#define FORT_NODEINTERPMG FORTRAN_NAME( inlineNODEINTERPMG, inlineNODEINTERPMG)
#define FORTNT_NODEINTERPMG FORTRAN_NAME( NODEINTERPMG, nodeinterpmg)

inline void FORTRAN_NAME(inlineNODEINTERPMG, inlineNODEINTERPMG)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(crse)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(nref)
      ,CHFp_BOX(bref)
      ,CHFp_FRA(wtcrnr) )
{
 CH_TIMELEAF("FORT_NODEINTERPMG");
 FORTRAN_NAME( NODEINTERPMG ,nodeinterpmg )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(crse)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(nref)
      ,CHFt_BOX(bref)
      ,CHFt_FRA(wtcrnr) );
}
#endif  // GUARDNODEINTERPMG 

#ifndef GUARDWTSIDE 
#define GUARDWTSIDE 
// Prototype for Fortran procedure WTSIDE ...
//
void FORTRAN_NAME( WTSIDE ,wtside )(
      CHFp_CONST_INT(i)
      ,CHFp_CONST_REAL(frac)
      ,CHFp_REAL(wt) );

#define FORT_WTSIDE FORTRAN_NAME( inlineWTSIDE, inlineWTSIDE)
#define FORTNT_WTSIDE FORTRAN_NAME( WTSIDE, wtside)

inline void FORTRAN_NAME(inlineWTSIDE, inlineWTSIDE)(
      CHFp_CONST_INT(i)
      ,CHFp_CONST_REAL(frac)
      ,CHFp_REAL(wt) )
{
 CH_TIMELEAF("FORT_WTSIDE");
 FORTRAN_NAME( WTSIDE ,wtside )(
      CHFt_CONST_INT(i)
      ,CHFt_CONST_REAL(frac)
      ,CHFt_REAL(wt) );
}
#endif  // GUARDWTSIDE 

#ifndef GUARDMAXB 
#define GUARDMAXB 
// Prototype for Fortran procedure MAXB ...
//
void FORTRAN_NAME( MAXB ,maxb )(
      CHFp_CONST_INT(iref)
      ,CHFp_INT(ibmax) );

#define FORT_MAXB FORTRAN_NAME( inlineMAXB, inlineMAXB)
#define FORTNT_MAXB FORTRAN_NAME( MAXB, maxb)

inline void FORTRAN_NAME(inlineMAXB, inlineMAXB)(
      CHFp_CONST_INT(iref)
      ,CHFp_INT(ibmax) )
{
 CH_TIMELEAF("FORT_MAXB");
 FORTRAN_NAME( MAXB ,maxb )(
      CHFt_CONST_INT(iref)
      ,CHFt_INT(ibmax) );
}
#endif  // GUARDMAXB 

}

#endif
