#ifndef _AMRPOISSONOPF_F_H_
#define _AMRPOISSONOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGSMCAMRPOP 
#define GUARDGSMCAMRPOP 
// Prototype for Fortran procedure gsmcamrpop ...
//
void FORTRAN_NAME( GSMCAMRPOP ,gsmcamrpop )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_BOX(coloredbox)
      ,CHFp_REAL(dx)
      ,CHFp_REAL(alpha)
      ,CHFp_REAL(beta) );

#define FORT_GSMCAMRPOP FORTRAN_NAME( inlineGSMCAMRPOP, inlineGSMCAMRPOP)
#define FORTNT_GSMCAMRPOP FORTRAN_NAME( GSMCAMRPOP, gsmcamrpop)

inline void FORTRAN_NAME(inlineGSMCAMRPOP, inlineGSMCAMRPOP)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_BOX(coloredbox)
      ,CHFp_REAL(dx)
      ,CHFp_REAL(alpha)
      ,CHFp_REAL(beta) )
{
 CH_TIMELEAF("FORT_GSMCAMRPOP");
 FORTRAN_NAME( GSMCAMRPOP ,gsmcamrpop )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_BOX(coloredbox)
      ,CHFt_REAL(dx)
      ,CHFt_REAL(alpha)
      ,CHFt_REAL(beta) );
}
#endif  // GUARDGSMCAMRPOP 

#ifndef GUARDGSRBHELMHOLTZ 
#define GUARDGSRBHELMHOLTZ 
// Prototype for Fortran procedure GSRBHELMHOLTZ ...
//
void FORTRAN_NAME( GSRBHELMHOLTZ ,gsrbhelmholtz )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBHELMHOLTZ FORTRAN_NAME( inlineGSRBHELMHOLTZ, inlineGSRBHELMHOLTZ)
#define FORTNT_GSRBHELMHOLTZ FORTRAN_NAME( GSRBHELMHOLTZ, gsrbhelmholtz)

inline void FORTRAN_NAME(inlineGSRBHELMHOLTZ, inlineGSRBHELMHOLTZ)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBHELMHOLTZ");
 FORTRAN_NAME( GSRBHELMHOLTZ ,gsrbhelmholtz )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBHELMHOLTZ 

#ifndef GUARDGSRBLAPLACIAN 
#define GUARDGSRBLAPLACIAN 
// Prototype for Fortran procedure GSRBLAPLACIAN ...
//
void FORTRAN_NAME( GSRBLAPLACIAN ,gsrblaplacian )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBLAPLACIAN FORTRAN_NAME( inlineGSRBLAPLACIAN, inlineGSRBLAPLACIAN)
#define FORTNT_GSRBLAPLACIAN FORTRAN_NAME( GSRBLAPLACIAN, gsrblaplacian)

inline void FORTRAN_NAME(inlineGSRBLAPLACIAN, inlineGSRBLAPLACIAN)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBLAPLACIAN");
 FORTRAN_NAME( GSRBLAPLACIAN ,gsrblaplacian )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBLAPLACIAN 

#ifndef GUARDGSRBLAZY 
#define GUARDGSRBLAZY 
// Prototype for Fortran procedure gsrblazy ...
//
void FORTRAN_NAME( GSRBLAZY ,gsrblazy )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(lphi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_BOX(coloredbox)
      ,CHFp_REAL(alpha)
      ,CHFp_REAL(beta)
      ,CHFp_REAL(dx) );

#define FORT_GSRBLAZY FORTRAN_NAME( inlineGSRBLAZY, inlineGSRBLAZY)
#define FORTNT_GSRBLAZY FORTRAN_NAME( GSRBLAZY, gsrblazy)

inline void FORTRAN_NAME(inlineGSRBLAZY, inlineGSRBLAZY)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(lphi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_BOX(coloredbox)
      ,CHFp_REAL(alpha)
      ,CHFp_REAL(beta)
      ,CHFp_REAL(dx) )
{
 CH_TIMELEAF("FORT_GSRBLAZY");
 FORTRAN_NAME( GSRBLAZY ,gsrblazy )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(lphi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_BOX(coloredbox)
      ,CHFt_REAL(alpha)
      ,CHFt_REAL(beta)
      ,CHFt_REAL(dx) );
}
#endif  // GUARDGSRBLAZY 

#ifndef GUARDAMRPMULTICOLOR 
#define GUARDAMRPMULTICOLOR 
// Prototype for Fortran procedure amrpmulticolor ...
//
void FORTRAN_NAME( AMRPMULTICOLOR ,amrpmulticolor )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(coloredbox) );

#define FORT_AMRPMULTICOLOR FORTRAN_NAME( inlineAMRPMULTICOLOR, inlineAMRPMULTICOLOR)
#define FORTNT_AMRPMULTICOLOR FORTRAN_NAME( AMRPMULTICOLOR, amrpmulticolor)

inline void FORTRAN_NAME(inlineAMRPMULTICOLOR, inlineAMRPMULTICOLOR)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(coloredbox) )
{
 CH_TIMELEAF("FORT_AMRPMULTICOLOR");
 FORTRAN_NAME( AMRPMULTICOLOR ,amrpmulticolor )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(weight)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(coloredbox) );
}
#endif  // GUARDAMRPMULTICOLOR 

#ifndef GUARDOPERATORLAP 
#define GUARDOPERATORLAP 
// Prototype for Fortran procedure OPERATORLAP ...
//
void FORTRAN_NAME( OPERATORLAP ,operatorlap )(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_OPERATORLAP FORTRAN_NAME( inlineOPERATORLAP, inlineOPERATORLAP)
#define FORTNT_OPERATORLAP FORTRAN_NAME( OPERATORLAP, operatorlap)

inline void FORTRAN_NAME(inlineOPERATORLAP, inlineOPERATORLAP)(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_OPERATORLAP");
 FORTRAN_NAME( OPERATORLAP ,operatorlap )(
      CHFt_FRA(lofphi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDOPERATORLAP 

#ifndef GUARDOPERATORLAPRES 
#define GUARDOPERATORLAPRES 
// Prototype for Fortran procedure OPERATORLAPRES ...
//
void FORTRAN_NAME( OPERATORLAPRES ,operatorlapres )(
      CHFp_FRA(r)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_OPERATORLAPRES FORTRAN_NAME( inlineOPERATORLAPRES, inlineOPERATORLAPRES)
#define FORTNT_OPERATORLAPRES FORTRAN_NAME( OPERATORLAPRES, operatorlapres)

inline void FORTRAN_NAME(inlineOPERATORLAPRES, inlineOPERATORLAPRES)(
      CHFp_FRA(r)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_OPERATORLAPRES");
 FORTRAN_NAME( OPERATORLAPRES ,operatorlapres )(
      CHFt_FRA(r)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDOPERATORLAPRES 

#ifndef GUARDRESTRICTRES 
#define GUARDRESTRICTRES 
// Prototype for Fortran procedure RESTRICTRES ...
//
void FORTRAN_NAME( RESTRICTRES ,restrictres )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_RESTRICTRES FORTRAN_NAME( inlineRESTRICTRES, inlineRESTRICTRES)
#define FORTNT_RESTRICTRES FORTRAN_NAME( RESTRICTRES, restrictres)

inline void FORTRAN_NAME(inlineRESTRICTRES, inlineRESTRICTRES)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_RESTRICTRES");
 FORTRAN_NAME( RESTRICTRES ,restrictres )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDRESTRICTRES 

#ifndef GUARDPROLONG 
#define GUARDPROLONG 
// Prototype for Fortran procedure PROLONG ...
//
void FORTRAN_NAME( PROLONG ,prolong )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(m) );

#define FORT_PROLONG FORTRAN_NAME( inlinePROLONG, inlinePROLONG)
#define FORTNT_PROLONG FORTRAN_NAME( PROLONG, prolong)

inline void FORTRAN_NAME(inlinePROLONG, inlinePROLONG)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(m) )
{
 CH_TIMELEAF("FORT_PROLONG");
 FORTRAN_NAME( PROLONG ,prolong )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(m) );
}
#endif  // GUARDPROLONG 

#ifndef GUARDPROLONG_2 
#define GUARDPROLONG_2 
// Prototype for Fortran procedure PROLONG_2 ...
//
void FORTRAN_NAME( PROLONG_2 ,prolong_2 )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(m) );

#define FORT_PROLONG_2 FORTRAN_NAME( inlinePROLONG_2, inlinePROLONG_2)
#define FORTNT_PROLONG_2 FORTRAN_NAME( PROLONG_2, prolong_2)

inline void FORTRAN_NAME(inlinePROLONG_2, inlinePROLONG_2)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(m) )
{
 CH_TIMELEAF("FORT_PROLONG_2");
 FORTRAN_NAME( PROLONG_2 ,prolong_2 )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(m) );
}
#endif  // GUARDPROLONG_2 

#ifndef GUARDNEWGETFLUX 
#define GUARDNEWGETFLUX 
// Prototype for Fortran procedure NEWGETFLUX ...
//
void FORTRAN_NAME( NEWGETFLUX ,newgetflux )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(beta_dx)
      ,CHFp_CONST_INT(a_idir) );

#define FORT_NEWGETFLUX FORTRAN_NAME( inlineNEWGETFLUX, inlineNEWGETFLUX)
#define FORTNT_NEWGETFLUX FORTRAN_NAME( NEWGETFLUX, newgetflux)

inline void FORTRAN_NAME(inlineNEWGETFLUX, inlineNEWGETFLUX)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(beta_dx)
      ,CHFp_CONST_INT(a_idir) )
{
 CH_TIMELEAF("FORT_NEWGETFLUX");
 FORTRAN_NAME( NEWGETFLUX ,newgetflux )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(beta_dx)
      ,CHFt_CONST_INT(a_idir) );
}
#endif  // GUARDNEWGETFLUX 

}

#endif
