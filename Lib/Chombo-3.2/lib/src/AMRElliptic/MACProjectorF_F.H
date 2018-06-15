#ifndef _MACPROJECTORF_F_H_
#define _MACPROJECTORF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDREGCORRECTTANVEL 
#define GUARDREGCORRECTTANVEL 
// Prototype for Fortran procedure regcorrecttanvel ...
//
void FORTRAN_NAME( REGCORRECTTANVEL ,regcorrecttanvel )(
      CHFp_FRA1(vel)
      ,CHFp_CONST_FRA1(grad)
      ,CHFp_BOX(interiorbox)
      ,CHFp_INT(veldir)
      ,CHFp_INT(graddir) );

#define FORT_REGCORRECTTANVEL FORTRAN_NAME( inlineREGCORRECTTANVEL, inlineREGCORRECTTANVEL)
#define FORTNT_REGCORRECTTANVEL FORTRAN_NAME( REGCORRECTTANVEL, regcorrecttanvel)

inline void FORTRAN_NAME(inlineREGCORRECTTANVEL, inlineREGCORRECTTANVEL)(
      CHFp_FRA1(vel)
      ,CHFp_CONST_FRA1(grad)
      ,CHFp_BOX(interiorbox)
      ,CHFp_INT(veldir)
      ,CHFp_INT(graddir) )
{
 CH_TIMELEAF("FORT_REGCORRECTTANVEL");
 FORTRAN_NAME( REGCORRECTTANVEL ,regcorrecttanvel )(
      CHFt_FRA1(vel)
      ,CHFt_CONST_FRA1(grad)
      ,CHFt_BOX(interiorbox)
      ,CHFt_INT(veldir)
      ,CHFt_INT(graddir) );
}
#endif  // GUARDREGCORRECTTANVEL 

#ifndef GUARDMACDIVERGEF 
#define GUARDMACDIVERGEF 
// Prototype for Fortran procedure macdivergef ...
//
void FORTRAN_NAME( MACDIVERGEF ,macdivergef )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(divf)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx) );

#define FORT_MACDIVERGEF FORTRAN_NAME( inlineMACDIVERGEF, inlineMACDIVERGEF)
#define FORTNT_MACDIVERGEF FORTRAN_NAME( MACDIVERGEF, macdivergef)

inline void FORTRAN_NAME(inlineMACDIVERGEF, inlineMACDIVERGEF)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(divf)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_MACDIVERGEF");
 FORTRAN_NAME( MACDIVERGEF ,macdivergef )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(divf)
      ,CHFt_CONST_FRA(flux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nconserved)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDMACDIVERGEF 

#ifndef GUARDMACGRADPHI 
#define GUARDMACGRADPHI 
// Prototype for Fortran procedure macgradphi ...
//
void FORTRAN_NAME( MACGRADPHI ,macgradphi )(
      CHFp_FRA1(gradphi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox) );

#define FORT_MACGRADPHI FORTRAN_NAME( inlineMACGRADPHI, inlineMACGRADPHI)
#define FORTNT_MACGRADPHI FORTRAN_NAME( MACGRADPHI, macgradphi)

inline void FORTRAN_NAME(inlineMACGRADPHI, inlineMACGRADPHI)(
      CHFp_FRA1(gradphi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_MACGRADPHI");
 FORTRAN_NAME( MACGRADPHI ,macgradphi )(
      CHFt_FRA1(gradphi)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDMACGRADPHI 

}

#endif
