#ifndef _NEWPOISSONOP4F_F_H_
#define _NEWPOISSONOP4F_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDOPERATORLAP4 
#define GUARDOPERATORLAP4 
// Prototype for Fortran procedure OPERATORLAP4 ...
//
void FORTRAN_NAME( OPERATORLAP4 ,operatorlap4 )(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_OPERATORLAP4 FORTRAN_NAME( inlineOPERATORLAP4, inlineOPERATORLAP4)
#define FORTNT_OPERATORLAP4 FORTRAN_NAME( OPERATORLAP4, operatorlap4)

inline void FORTRAN_NAME(inlineOPERATORLAP4, inlineOPERATORLAP4)(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_OPERATORLAP4");
 FORTRAN_NAME( OPERATORLAP4 ,operatorlap4 )(
      CHFt_FRA(lofphi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDOPERATORLAP4 

#ifndef GUARDOPERATORLAPRES4 
#define GUARDOPERATORLAPRES4 
// Prototype for Fortran procedure OPERATORLAPRES4 ...
//
void FORTRAN_NAME( OPERATORLAPRES4 ,operatorlapres4 )(
      CHFp_FRA(r)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_OPERATORLAPRES4 FORTRAN_NAME( inlineOPERATORLAPRES4, inlineOPERATORLAPRES4)
#define FORTNT_OPERATORLAPRES4 FORTRAN_NAME( OPERATORLAPRES4, operatorlapres4)

inline void FORTRAN_NAME(inlineOPERATORLAPRES4, inlineOPERATORLAPRES4)(
      CHFp_FRA(r)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_OPERATORLAPRES4");
 FORTRAN_NAME( OPERATORLAPRES4 ,operatorlapres4 )(
      CHFt_FRA(r)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDOPERATORLAPRES4 

#ifndef GUARDRESTRICTRES4 
#define GUARDRESTRICTRES4 
// Prototype for Fortran procedure RESTRICTRES4 ...
//
void FORTRAN_NAME( RESTRICTRES4 ,restrictres4 )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_RESTRICTRES4 FORTRAN_NAME( inlineRESTRICTRES4, inlineRESTRICTRES4)
#define FORTNT_RESTRICTRES4 FORTRAN_NAME( RESTRICTRES4, restrictres4)

inline void FORTRAN_NAME(inlineRESTRICTRES4, inlineRESTRICTRES4)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_RESTRICTRES4");
 FORTRAN_NAME( RESTRICTRES4 ,restrictres4 )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDRESTRICTRES4 

#ifndef GUARDGSRBLAPLACIAN4 
#define GUARDGSRBLAPLACIAN4 
// Prototype for Fortran procedure GSRBLAPLACIAN4 ...
//
void FORTRAN_NAME( GSRBLAPLACIAN4 ,gsrblaplacian4 )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_FRA(tmp)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBLAPLACIAN4 FORTRAN_NAME( inlineGSRBLAPLACIAN4, inlineGSRBLAPLACIAN4)
#define FORTNT_GSRBLAPLACIAN4 FORTRAN_NAME( GSRBLAPLACIAN4, gsrblaplacian4)

inline void FORTRAN_NAME(inlineGSRBLAPLACIAN4, inlineGSRBLAPLACIAN4)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_FRA(tmp)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBLAPLACIAN4");
 FORTRAN_NAME( GSRBLAPLACIAN4 ,gsrblaplacian4 )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_FRA(tmp)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBLAPLACIAN4 

#ifndef GUARDSORLAPLACIAN4 
#define GUARDSORLAPLACIAN4 
// Prototype for Fortran procedure SORLAPLACIAN4 ...
//
void FORTRAN_NAME( SORLAPLACIAN4 ,sorlaplacian4 )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_SORLAPLACIAN4 FORTRAN_NAME( inlineSORLAPLACIAN4, inlineSORLAPLACIAN4)
#define FORTNT_SORLAPLACIAN4 FORTRAN_NAME( SORLAPLACIAN4, sorlaplacian4)

inline void FORTRAN_NAME(inlineSORLAPLACIAN4, inlineSORLAPLACIAN4)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_SORLAPLACIAN4");
 FORTRAN_NAME( SORLAPLACIAN4 ,sorlaplacian4 )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDSORLAPLACIAN4 

#ifndef GUARDGSRBHELMHOLTZ4 
#define GUARDGSRBHELMHOLTZ4 
// Prototype for Fortran procedure GSRBHELMHOLTZ4 ...
//
void FORTRAN_NAME( GSRBHELMHOLTZ4 ,gsrbhelmholtz4 )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_FRA(tmp)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBHELMHOLTZ4 FORTRAN_NAME( inlineGSRBHELMHOLTZ4, inlineGSRBHELMHOLTZ4)
#define FORTNT_GSRBHELMHOLTZ4 FORTRAN_NAME( GSRBHELMHOLTZ4, gsrbhelmholtz4)

inline void FORTRAN_NAME(inlineGSRBHELMHOLTZ4, inlineGSRBHELMHOLTZ4)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_FRA(tmp)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBHELMHOLTZ4");
 FORTRAN_NAME( GSRBHELMHOLTZ4 ,gsrbhelmholtz4 )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_FRA(tmp)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBHELMHOLTZ4 

#ifndef GUARDNEWGETFLUX4 
#define GUARDNEWGETFLUX4 
// Prototype for Fortran procedure NEWGETFLUX4 ...
//
void FORTRAN_NAME( NEWGETFLUX4 ,newgetflux4 )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(beta_dx)
      ,CHFp_CONST_INT(a_idir) );

#define FORT_NEWGETFLUX4 FORTRAN_NAME( inlineNEWGETFLUX4, inlineNEWGETFLUX4)
#define FORTNT_NEWGETFLUX4 FORTRAN_NAME( NEWGETFLUX4, newgetflux4)

inline void FORTRAN_NAME(inlineNEWGETFLUX4, inlineNEWGETFLUX4)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(beta_dx)
      ,CHFp_CONST_INT(a_idir) )
{
 CH_TIMELEAF("FORT_NEWGETFLUX4");
 FORTRAN_NAME( NEWGETFLUX4 ,newgetflux4 )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(beta_dx)
      ,CHFt_CONST_INT(a_idir) );
}
#endif  // GUARDNEWGETFLUX4 

#ifndef GUARDPROLONGLINEAR 
#define GUARDPROLONGLINEAR 
// Prototype for Fortran procedure PROLONGLINEAR ...
//
void FORTRAN_NAME( PROLONGLINEAR ,prolonglinear )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(fineBox)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(r) );

#define FORT_PROLONGLINEAR FORTRAN_NAME( inlinePROLONGLINEAR, inlinePROLONGLINEAR)
#define FORTNT_PROLONGLINEAR FORTRAN_NAME( PROLONGLINEAR, prolonglinear)

inline void FORTRAN_NAME(inlinePROLONGLINEAR, inlinePROLONGLINEAR)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(fineBox)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(r) )
{
 CH_TIMELEAF("FORT_PROLONGLINEAR");
 FORTRAN_NAME( PROLONGLINEAR ,prolonglinear )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(fineBox)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_INT(r) );
}
#endif  // GUARDPROLONGLINEAR 

}

#endif
