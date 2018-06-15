#ifndef _LEVELTGAF_F_H_
#define _LEVELTGAF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRESGHOSTBC 
#define GUARDRESGHOSTBC 
// Prototype for Fortran procedure RESGHOSTBC ...
//
void FORTRAN_NAME( RESGHOSTBC ,resghostbc )(
      CHFp_FRA(state)
      ,CHFp_BOX(bcBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(ncomp) );

#define FORT_RESGHOSTBC FORTRAN_NAME( inlineRESGHOSTBC, inlineRESGHOSTBC)
#define FORTNT_RESGHOSTBC FORTRAN_NAME( RESGHOSTBC, resghostbc)

inline void FORTRAN_NAME(inlineRESGHOSTBC, inlineRESGHOSTBC)(
      CHFp_FRA(state)
      ,CHFp_BOX(bcBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(ncomp) )
{
 CH_TIMELEAF("FORT_RESGHOSTBC");
 FORTRAN_NAME( RESGHOSTBC ,resghostbc )(
      CHFt_FRA(state)
      ,CHFt_BOX(bcBox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(ncomp) );
}
#endif  // GUARDRESGHOSTBC 

#ifndef GUARDHORESGHOSTBC 
#define GUARDHORESGHOSTBC 
// Prototype for Fortran procedure HORESGHOSTBC ...
//
void FORTRAN_NAME( HORESGHOSTBC ,horesghostbc )(
      CHFp_FRA(state)
      ,CHFp_BOX(bcBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(ncomp) );

#define FORT_HORESGHOSTBC FORTRAN_NAME( inlineHORESGHOSTBC, inlineHORESGHOSTBC)
#define FORTNT_HORESGHOSTBC FORTRAN_NAME( HORESGHOSTBC, horesghostbc)

inline void FORTRAN_NAME(inlineHORESGHOSTBC, inlineHORESGHOSTBC)(
      CHFp_FRA(state)
      ,CHFp_BOX(bcBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(ncomp) )
{
 CH_TIMELEAF("FORT_HORESGHOSTBC");
 FORTRAN_NAME( HORESGHOSTBC ,horesghostbc )(
      CHFt_FRA(state)
      ,CHFt_BOX(bcBox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(ncomp) );
}
#endif  // GUARDHORESGHOSTBC 

}

#endif
