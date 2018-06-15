#ifndef _INTERPF_F_H_
#define _INTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDINTERPCONSTANT 
#define GUARDINTERPCONSTANT 
// Prototype for Fortran procedure interpconstant ...
//
void FORTRAN_NAME( INTERPCONSTANT ,interpconstant )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(bref) );

#define FORT_INTERPCONSTANT FORTRAN_NAME( inlineINTERPCONSTANT, inlineINTERPCONSTANT)
#define FORTNT_INTERPCONSTANT FORTRAN_NAME( INTERPCONSTANT, interpconstant)

inline void FORTRAN_NAME(inlineINTERPCONSTANT, inlineINTERPCONSTANT)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_INTERPCONSTANT");
 FORTRAN_NAME( INTERPCONSTANT ,interpconstant )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(ref_ratio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDINTERPCONSTANT 

#ifndef GUARDINTERPCENTRALSLOPE 
#define GUARDINTERPCENTRALSLOPE 
// Prototype for Fortran procedure interpcentralslope ...
//
void FORTRAN_NAME( INTERPCENTRALSLOPE ,interpcentralslope )(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) );

#define FORT_INTERPCENTRALSLOPE FORTRAN_NAME( inlineINTERPCENTRALSLOPE, inlineINTERPCENTRALSLOPE)
#define FORTNT_INTERPCENTRALSLOPE FORTRAN_NAME( INTERPCENTRALSLOPE, interpcentralslope)

inline void FORTRAN_NAME(inlineINTERPCENTRALSLOPE, inlineINTERPCENTRALSLOPE)(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_INTERPCENTRALSLOPE");
 FORTRAN_NAME( INTERPCENTRALSLOPE ,interpcentralslope )(
      CHFt_FRA(slope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDINTERPCENTRALSLOPE 

#ifndef GUARDINTERPHISIDESLOPE 
#define GUARDINTERPHISIDESLOPE 
// Prototype for Fortran procedure interphisideslope ...
//
void FORTRAN_NAME( INTERPHISIDESLOPE ,interphisideslope )(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) );

#define FORT_INTERPHISIDESLOPE FORTRAN_NAME( inlineINTERPHISIDESLOPE, inlineINTERPHISIDESLOPE)
#define FORTNT_INTERPHISIDESLOPE FORTRAN_NAME( INTERPHISIDESLOPE, interphisideslope)

inline void FORTRAN_NAME(inlineINTERPHISIDESLOPE, inlineINTERPHISIDESLOPE)(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_INTERPHISIDESLOPE");
 FORTRAN_NAME( INTERPHISIDESLOPE ,interphisideslope )(
      CHFt_FRA(slope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDINTERPHISIDESLOPE 

#ifndef GUARDINTERPLOSIDESLOPE 
#define GUARDINTERPLOSIDESLOPE 
// Prototype for Fortran procedure interplosideslope ...
//
void FORTRAN_NAME( INTERPLOSIDESLOPE ,interplosideslope )(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) );

#define FORT_INTERPLOSIDESLOPE FORTRAN_NAME( inlineINTERPLOSIDESLOPE, inlineINTERPLOSIDESLOPE)
#define FORTNT_INTERPLOSIDESLOPE FORTRAN_NAME( INTERPLOSIDESLOPE, interplosideslope)

inline void FORTRAN_NAME(inlineINTERPLOSIDESLOPE, inlineINTERPLOSIDESLOPE)(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_INTERPLOSIDESLOPE");
 FORTRAN_NAME( INTERPLOSIDESLOPE ,interplosideslope )(
      CHFt_FRA(slope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDINTERPLOSIDESLOPE 

#ifndef GUARDINTERPLIMIT 
#define GUARDINTERPLIMIT 
// Prototype for Fortran procedure interplimit ...
//
void FORTRAN_NAME( INTERPLIMIT ,interplimit )(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(bcoarse)
      ,CHFp_BOX(bn)
      ,CHFp_BOX(physdomain) );

#define FORT_INTERPLIMIT FORTRAN_NAME( inlineINTERPLIMIT, inlineINTERPLIMIT)
#define FORTNT_INTERPLIMIT FORTRAN_NAME( INTERPLIMIT, interplimit)

inline void FORTRAN_NAME(inlineINTERPLIMIT, inlineINTERPLIMIT)(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(bcoarse)
      ,CHFp_BOX(bn)
      ,CHFp_BOX(physdomain) )
{
 CH_TIMELEAF("FORT_INTERPLIMIT");
 FORTRAN_NAME( INTERPLIMIT ,interplimit )(
      CHFt_FRA(islope)
      ,CHFt_FRA(jslope)
      ,CHFt_FRA(kslope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(bcoarse)
      ,CHFt_BOX(bn)
      ,CHFt_BOX(physdomain) );
}
#endif  // GUARDINTERPLIMIT 

#ifndef GUARDINTERPLINEAR 
#define GUARDINTERPLINEAR 
// Prototype for Fortran procedure interplinear ...
//
void FORTRAN_NAME( INTERPLINEAR ,interplinear )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(bref) );

#define FORT_INTERPLINEAR FORTRAN_NAME( inlineINTERPLINEAR, inlineINTERPLINEAR)
#define FORTNT_INTERPLINEAR FORTRAN_NAME( INTERPLINEAR, interplinear)

inline void FORTRAN_NAME(inlineINTERPLINEAR, inlineINTERPLINEAR)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_INTERPLINEAR");
 FORTRAN_NAME( INTERPLINEAR ,interplinear )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(slope)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(ref_ratio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDINTERPLINEAR 

#ifndef GUARDINTERPHOMO_OLD 
#define GUARDINTERPHOMO_OLD 
// Prototype for Fortran procedure INTERPHOMO_OLD ...
//
void FORTRAN_NAME( INTERPHOMO_OLD ,interphomo_old )(
      CHFp_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(x1)
      ,CHFp_CONST_REAL(dxCrse)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ihilo) );

#define FORT_INTERPHOMO_OLD FORTRAN_NAME( inlineINTERPHOMO_OLD, inlineINTERPHOMO_OLD)
#define FORTNT_INTERPHOMO_OLD FORTRAN_NAME( INTERPHOMO_OLD, interphomo_old)

inline void FORTRAN_NAME(inlineINTERPHOMO_OLD, inlineINTERPHOMO_OLD)(
      CHFp_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(x1)
      ,CHFp_CONST_REAL(dxCrse)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ihilo) )
{
 CH_TIMELEAF("FORT_INTERPHOMO_OLD");
 FORTRAN_NAME( INTERPHOMO_OLD ,interphomo_old )(
      CHFt_FRA(phi)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(x1)
      ,CHFt_CONST_REAL(dxCrse)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(ihilo) );
}
#endif  // GUARDINTERPHOMO_OLD 

#ifndef GUARDINTERPHOMO 
#define GUARDINTERPHOMO 
// Prototype for Fortran procedure INTERPHOMO ...
//
void FORTRAN_NAME( INTERPHOMO ,interphomo )(
      CHFp_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(x1)
      ,CHFp_CONST_REAL(dxCrse)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ihilo) );

#define FORT_INTERPHOMO FORTRAN_NAME( inlineINTERPHOMO, inlineINTERPHOMO)
#define FORTNT_INTERPHOMO FORTRAN_NAME( INTERPHOMO, interphomo)

inline void FORTRAN_NAME(inlineINTERPHOMO, inlineINTERPHOMO)(
      CHFp_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(x1)
      ,CHFp_CONST_REAL(dxCrse)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ihilo) )
{
 CH_TIMELEAF("FORT_INTERPHOMO");
 FORTRAN_NAME( INTERPHOMO ,interphomo )(
      CHFt_FRA(phi)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(x1)
      ,CHFt_CONST_REAL(dxCrse)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(ihilo) );
}
#endif  // GUARDINTERPHOMO 

#ifndef GUARDINTERPHOMOLINEAR 
#define GUARDINTERPHOMOLINEAR 
// Prototype for Fortran procedure INTERPHOMOLINEAR ...
//
void FORTRAN_NAME( INTERPHOMOLINEAR ,interphomolinear )(
      CHFp_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(x1)
      ,CHFp_CONST_REAL(dxCrse)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ihilo) );

#define FORT_INTERPHOMOLINEAR FORTRAN_NAME( inlineINTERPHOMOLINEAR, inlineINTERPHOMOLINEAR)
#define FORTNT_INTERPHOMOLINEAR FORTRAN_NAME( INTERPHOMOLINEAR, interphomolinear)

inline void FORTRAN_NAME(inlineINTERPHOMOLINEAR, inlineINTERPHOMOLINEAR)(
      CHFp_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(x1)
      ,CHFp_CONST_REAL(dxCrse)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ihilo) )
{
 CH_TIMELEAF("FORT_INTERPHOMOLINEAR");
 FORTRAN_NAME( INTERPHOMOLINEAR ,interphomolinear )(
      CHFt_FRA(phi)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(x1)
      ,CHFt_CONST_REAL(dxCrse)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(ihilo) );
}
#endif  // GUARDINTERPHOMOLINEAR 

}

#endif
