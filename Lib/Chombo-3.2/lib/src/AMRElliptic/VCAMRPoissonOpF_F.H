#ifndef _VCAMRPOISSONOPF_F_H_
#define _VCAMRPOISSONOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGSRBHELMHOLTZVC1D 
#define GUARDGSRBHELMHOLTZVC1D 
// Prototype for Fortran procedure GSRBHELMHOLTZVC1D ...
//
void FORTRAN_NAME( GSRBHELMHOLTZVC1D ,gsrbhelmholtzvc1d )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBHELMHOLTZVC1D FORTRAN_NAME( inlineGSRBHELMHOLTZVC1D, inlineGSRBHELMHOLTZVC1D)
#define FORTNT_GSRBHELMHOLTZVC1D FORTRAN_NAME( GSRBHELMHOLTZVC1D, gsrbhelmholtzvc1d)

inline void FORTRAN_NAME(inlineGSRBHELMHOLTZVC1D, inlineGSRBHELMHOLTZVC1D)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBHELMHOLTZVC1D");
 FORTRAN_NAME( GSRBHELMHOLTZVC1D ,gsrbhelmholtzvc1d )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(lambda)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBHELMHOLTZVC1D 

#ifndef GUARDGSRBHELMHOLTZVC2D 
#define GUARDGSRBHELMHOLTZVC2D 
// Prototype for Fortran procedure GSRBHELMHOLTZVC2D ...
//
void FORTRAN_NAME( GSRBHELMHOLTZVC2D ,gsrbhelmholtzvc2d )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBHELMHOLTZVC2D FORTRAN_NAME( inlineGSRBHELMHOLTZVC2D, inlineGSRBHELMHOLTZVC2D)
#define FORTNT_GSRBHELMHOLTZVC2D FORTRAN_NAME( GSRBHELMHOLTZVC2D, gsrbhelmholtzvc2d)

inline void FORTRAN_NAME(inlineGSRBHELMHOLTZVC2D, inlineGSRBHELMHOLTZVC2D)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBHELMHOLTZVC2D");
 FORTRAN_NAME( GSRBHELMHOLTZVC2D ,gsrbhelmholtzvc2d )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_CONST_FRA(lambda)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBHELMHOLTZVC2D 

#ifndef GUARDGSRBHELMHOLTZVC3D 
#define GUARDGSRBHELMHOLTZVC3D 
// Prototype for Fortran procedure GSRBHELMHOLTZVC3D ...
//
void FORTRAN_NAME( GSRBHELMHOLTZVC3D ,gsrbhelmholtzvc3d )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBHELMHOLTZVC3D FORTRAN_NAME( inlineGSRBHELMHOLTZVC3D, inlineGSRBHELMHOLTZVC3D)
#define FORTNT_GSRBHELMHOLTZVC3D FORTRAN_NAME( GSRBHELMHOLTZVC3D, gsrbhelmholtzvc3d)

inline void FORTRAN_NAME(inlineGSRBHELMHOLTZVC3D, inlineGSRBHELMHOLTZVC3D)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBHELMHOLTZVC3D");
 FORTRAN_NAME( GSRBHELMHOLTZVC3D ,gsrbhelmholtzvc3d )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_CONST_FRA(bCoef2)
      ,CHFt_CONST_FRA(lambda)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBHELMHOLTZVC3D 

#ifndef GUARDVCCOMPUTEOP1D 
#define GUARDVCCOMPUTEOP1D 
// Prototype for Fortran procedure VCCOMPUTEOP1D ...
//
void FORTRAN_NAME( VCCOMPUTEOP1D ,vccomputeop1d )(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_VCCOMPUTEOP1D FORTRAN_NAME( inlineVCCOMPUTEOP1D, inlineVCCOMPUTEOP1D)
#define FORTNT_VCCOMPUTEOP1D FORTRAN_NAME( VCCOMPUTEOP1D, vccomputeop1d)

inline void FORTRAN_NAME(inlineVCCOMPUTEOP1D, inlineVCCOMPUTEOP1D)(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_VCCOMPUTEOP1D");
 FORTRAN_NAME( VCCOMPUTEOP1D ,vccomputeop1d )(
      CHFt_FRA(lofphi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDVCCOMPUTEOP1D 

#ifndef GUARDVCCOMPUTEOP2D 
#define GUARDVCCOMPUTEOP2D 
// Prototype for Fortran procedure VCCOMPUTEOP2D ...
//
void FORTRAN_NAME( VCCOMPUTEOP2D ,vccomputeop2d )(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_VCCOMPUTEOP2D FORTRAN_NAME( inlineVCCOMPUTEOP2D, inlineVCCOMPUTEOP2D)
#define FORTNT_VCCOMPUTEOP2D FORTRAN_NAME( VCCOMPUTEOP2D, vccomputeop2d)

inline void FORTRAN_NAME(inlineVCCOMPUTEOP2D, inlineVCCOMPUTEOP2D)(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_VCCOMPUTEOP2D");
 FORTRAN_NAME( VCCOMPUTEOP2D ,vccomputeop2d )(
      CHFt_FRA(lofphi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDVCCOMPUTEOP2D 

#ifndef GUARDVCCOMPUTEOP3D 
#define GUARDVCCOMPUTEOP3D 
// Prototype for Fortran procedure VCCOMPUTEOP3D ...
//
void FORTRAN_NAME( VCCOMPUTEOP3D ,vccomputeop3d )(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_VCCOMPUTEOP3D FORTRAN_NAME( inlineVCCOMPUTEOP3D, inlineVCCOMPUTEOP3D)
#define FORTNT_VCCOMPUTEOP3D FORTRAN_NAME( VCCOMPUTEOP3D, vccomputeop3d)

inline void FORTRAN_NAME(inlineVCCOMPUTEOP3D, inlineVCCOMPUTEOP3D)(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_VCCOMPUTEOP3D");
 FORTRAN_NAME( VCCOMPUTEOP3D ,vccomputeop3d )(
      CHFt_FRA(lofphi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_CONST_FRA(bCoef2)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDVCCOMPUTEOP3D 

#ifndef GUARDVCCOMPUTERES1D 
#define GUARDVCCOMPUTERES1D 
// Prototype for Fortran procedure VCCOMPUTERES1D ...
//
void FORTRAN_NAME( VCCOMPUTERES1D ,vccomputeres1d )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_VCCOMPUTERES1D FORTRAN_NAME( inlineVCCOMPUTERES1D, inlineVCCOMPUTERES1D)
#define FORTNT_VCCOMPUTERES1D FORTRAN_NAME( VCCOMPUTERES1D, vccomputeres1d)

inline void FORTRAN_NAME(inlineVCCOMPUTERES1D, inlineVCCOMPUTERES1D)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_VCCOMPUTERES1D");
 FORTRAN_NAME( VCCOMPUTERES1D ,vccomputeres1d )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDVCCOMPUTERES1D 

#ifndef GUARDVCCOMPUTERES2D 
#define GUARDVCCOMPUTERES2D 
// Prototype for Fortran procedure VCCOMPUTERES2D ...
//
void FORTRAN_NAME( VCCOMPUTERES2D ,vccomputeres2d )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_VCCOMPUTERES2D FORTRAN_NAME( inlineVCCOMPUTERES2D, inlineVCCOMPUTERES2D)
#define FORTNT_VCCOMPUTERES2D FORTRAN_NAME( VCCOMPUTERES2D, vccomputeres2d)

inline void FORTRAN_NAME(inlineVCCOMPUTERES2D, inlineVCCOMPUTERES2D)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_VCCOMPUTERES2D");
 FORTRAN_NAME( VCCOMPUTERES2D ,vccomputeres2d )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDVCCOMPUTERES2D 

#ifndef GUARDVCCOMPUTERES3D 
#define GUARDVCCOMPUTERES3D 
// Prototype for Fortran procedure VCCOMPUTERES3D ...
//
void FORTRAN_NAME( VCCOMPUTERES3D ,vccomputeres3d )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_VCCOMPUTERES3D FORTRAN_NAME( inlineVCCOMPUTERES3D, inlineVCCOMPUTERES3D)
#define FORTNT_VCCOMPUTERES3D FORTRAN_NAME( VCCOMPUTERES3D, vccomputeres3d)

inline void FORTRAN_NAME(inlineVCCOMPUTERES3D, inlineVCCOMPUTERES3D)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_VCCOMPUTERES3D");
 FORTRAN_NAME( VCCOMPUTERES3D ,vccomputeres3d )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_CONST_FRA(bCoef2)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDVCCOMPUTERES3D 

#ifndef GUARDRESTRICTRESVC1D 
#define GUARDRESTRICTRESVC1D 
// Prototype for Fortran procedure RESTRICTRESVC1D ...
//
void FORTRAN_NAME( RESTRICTRESVC1D ,restrictresvc1d )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_RESTRICTRESVC1D FORTRAN_NAME( inlineRESTRICTRESVC1D, inlineRESTRICTRESVC1D)
#define FORTNT_RESTRICTRESVC1D FORTRAN_NAME( RESTRICTRESVC1D, restrictresvc1d)

inline void FORTRAN_NAME(inlineRESTRICTRESVC1D, inlineRESTRICTRESVC1D)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_RESTRICTRESVC1D");
 FORTRAN_NAME( RESTRICTRESVC1D ,restrictresvc1d )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDRESTRICTRESVC1D 

#ifndef GUARDRESTRICTRESVC2D 
#define GUARDRESTRICTRESVC2D 
// Prototype for Fortran procedure RESTRICTRESVC2D ...
//
void FORTRAN_NAME( RESTRICTRESVC2D ,restrictresvc2d )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_RESTRICTRESVC2D FORTRAN_NAME( inlineRESTRICTRESVC2D, inlineRESTRICTRESVC2D)
#define FORTNT_RESTRICTRESVC2D FORTRAN_NAME( RESTRICTRESVC2D, restrictresvc2d)

inline void FORTRAN_NAME(inlineRESTRICTRESVC2D, inlineRESTRICTRESVC2D)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_RESTRICTRESVC2D");
 FORTRAN_NAME( RESTRICTRESVC2D ,restrictresvc2d )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDRESTRICTRESVC2D 

#ifndef GUARDRESTRICTRESVC3D 
#define GUARDRESTRICTRESVC3D 
// Prototype for Fortran procedure RESTRICTRESVC3D ...
//
void FORTRAN_NAME( RESTRICTRESVC3D ,restrictresvc3d )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_RESTRICTRESVC3D FORTRAN_NAME( inlineRESTRICTRESVC3D, inlineRESTRICTRESVC3D)
#define FORTNT_RESTRICTRESVC3D FORTRAN_NAME( RESTRICTRESVC3D, restrictresvc3d)

inline void FORTRAN_NAME(inlineRESTRICTRESVC3D, inlineRESTRICTRESVC3D)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_FRA(aCoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoef0)
      ,CHFp_CONST_FRA(bCoef1)
      ,CHFp_CONST_FRA(bCoef2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_RESTRICTRESVC3D");
 FORTRAN_NAME( RESTRICTRESVC3D ,restrictresvc3d )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_FRA(aCoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoef0)
      ,CHFt_CONST_FRA(bCoef1)
      ,CHFt_CONST_FRA(bCoef2)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDRESTRICTRESVC3D 

#ifndef GUARDSUMFACES 
#define GUARDSUMFACES 
// Prototype for Fortran procedure SUMFACES ...
//
void FORTRAN_NAME( SUMFACES ,sumfaces )(
      CHFp_FRA(lhs)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoefs)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(scale) );

#define FORT_SUMFACES FORTRAN_NAME( inlineSUMFACES, inlineSUMFACES)
#define FORTNT_SUMFACES FORTRAN_NAME( SUMFACES, sumfaces)

inline void FORTRAN_NAME(inlineSUMFACES, inlineSUMFACES)(
      CHFp_FRA(lhs)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_FRA(bCoefs)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(scale) )
{
 CH_TIMELEAF("FORT_SUMFACES");
 FORTRAN_NAME( SUMFACES ,sumfaces )(
      CHFt_FRA(lhs)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_FRA(bCoefs)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(scale) );
}
#endif  // GUARDSUMFACES 

}

#endif
