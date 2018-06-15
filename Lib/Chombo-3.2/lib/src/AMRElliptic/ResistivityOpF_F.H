#ifndef _RESISTIVITYOPF_F_H_
#define _RESISTIVITYOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDPROLONGROP 
#define GUARDPROLONGROP 
// Prototype for Fortran procedure PROLONGROP ...
//
void FORTRAN_NAME( PROLONGROP ,prolongrop )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(m) );

#define FORT_PROLONGROP FORTRAN_NAME( inlinePROLONGROP, inlinePROLONGROP)
#define FORTNT_PROLONGROP FORTRAN_NAME( PROLONGROP, prolongrop)

inline void FORTRAN_NAME(inlinePROLONGROP, inlinePROLONGROP)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(m) )
{
 CH_TIMELEAF("FORT_PROLONGROP");
 FORTRAN_NAME( PROLONGROP ,prolongrop )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(m) );
}
#endif  // GUARDPROLONGROP 

#ifndef GUARDRESTRICTRESROP 
#define GUARDRESTRICTRESROP 
// Prototype for Fortran procedure RESTRICTRESROP ...
//
void FORTRAN_NAME( RESTRICTRESROP ,restrictresrop )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(resFine)
      ,CHFp_BOX(region) );

#define FORT_RESTRICTRESROP FORTRAN_NAME( inlineRESTRICTRESROP, inlineRESTRICTRESROP)
#define FORTNT_RESTRICTRESROP FORTRAN_NAME( RESTRICTRESROP, restrictresrop)

inline void FORTRAN_NAME(inlineRESTRICTRESROP, inlineRESTRICTRESROP)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(resFine)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_RESTRICTRESROP");
 FORTRAN_NAME( RESTRICTRESROP ,restrictresrop )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(resFine)
      ,CHFt_BOX(region) );
}
#endif  // GUARDRESTRICTRESROP 

#ifndef GUARDCELLGRADROP 
#define GUARDCELLGRADROP 
// Prototype for Fortran procedure cellgradrop ...
//
void FORTRAN_NAME( CELLGRADROP ,cellgradrop )(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(divdir) );

#define FORT_CELLGRADROP FORTRAN_NAME( inlineCELLGRADROP, inlineCELLGRADROP)
#define FORTNT_CELLGRADROP FORTRAN_NAME( CELLGRADROP, cellgradrop)

inline void FORTRAN_NAME(inlineCELLGRADROP, inlineCELLGRADROP)(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(divdir) )
{
 CH_TIMELEAF("FORT_CELLGRADROP");
 FORTRAN_NAME( CELLGRADROP ,cellgradrop )(
      CHFt_FRA1(grad)
      ,CHFt_CONST_FRA1(vel)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(divdir) );
}
#endif  // GUARDCELLGRADROP 

#ifndef GUARDADDGRADTOFLUXROP 
#define GUARDADDGRADTOFLUXROP 
// Prototype for Fortran procedure addgradtofluxrop ...
//
void FORTRAN_NAME( ADDGRADTOFLUXROP ,addgradtofluxrop )(
      CHFp_FRA(flux)
      ,CHFp_INT(fluxComp)
      ,CHFp_FRA(grad)
      ,CHFp_INT(gradComp)
      ,CHFp_REAL(sign)
      ,CHFp_BOX(regionface) );

#define FORT_ADDGRADTOFLUXROP FORTRAN_NAME( inlineADDGRADTOFLUXROP, inlineADDGRADTOFLUXROP)
#define FORTNT_ADDGRADTOFLUXROP FORTRAN_NAME( ADDGRADTOFLUXROP, addgradtofluxrop)

inline void FORTRAN_NAME(inlineADDGRADTOFLUXROP, inlineADDGRADTOFLUXROP)(
      CHFp_FRA(flux)
      ,CHFp_INT(fluxComp)
      ,CHFp_FRA(grad)
      ,CHFp_INT(gradComp)
      ,CHFp_REAL(sign)
      ,CHFp_BOX(regionface) )
{
 CH_TIMELEAF("FORT_ADDGRADTOFLUXROP");
 FORTRAN_NAME( ADDGRADTOFLUXROP ,addgradtofluxrop )(
      CHFt_FRA(flux)
      ,CHFt_INT(fluxComp)
      ,CHFt_FRA(grad)
      ,CHFt_INT(gradComp)
      ,CHFt_REAL(sign)
      ,CHFt_BOX(regionface) );
}
#endif  // GUARDADDGRADTOFLUXROP 

#ifndef GUARDGETFACEGRADROP 
#define GUARDGETFACEGRADROP 
// Prototype for Fortran procedure getfacegradrop ...
//
void FORTRAN_NAME( GETFACEGRADROP ,getfacegradrop )(
      CHFp_FRA1(gradvelface)
      ,CHFp_CONST_FRA1(gradvelcell)
      ,CHFp_CONST_FRA1(velcomp)
      ,CHFp_BOX(region)
      ,CHFp_BOX(centerbox)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(divdir) );

#define FORT_GETFACEGRADROP FORTRAN_NAME( inlineGETFACEGRADROP, inlineGETFACEGRADROP)
#define FORTNT_GETFACEGRADROP FORTRAN_NAME( GETFACEGRADROP, getfacegradrop)

inline void FORTRAN_NAME(inlineGETFACEGRADROP, inlineGETFACEGRADROP)(
      CHFp_FRA1(gradvelface)
      ,CHFp_CONST_FRA1(gradvelcell)
      ,CHFp_CONST_FRA1(velcomp)
      ,CHFp_BOX(region)
      ,CHFp_BOX(centerbox)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(divdir) )
{
 CH_TIMELEAF("FORT_GETFACEGRADROP");
 FORTRAN_NAME( GETFACEGRADROP ,getfacegradrop )(
      CHFt_FRA1(gradvelface)
      ,CHFt_CONST_FRA1(gradvelcell)
      ,CHFt_CONST_FRA1(velcomp)
      ,CHFt_BOX(region)
      ,CHFt_BOX(centerbox)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(divdir) );
}
#endif  // GUARDGETFACEGRADROP 

#ifndef GUARDCELLDIVINCRROP 
#define GUARDCELLDIVINCRROP 
// Prototype for Fortran procedure celldivincrrop ...
//
void FORTRAN_NAME( CELLDIVINCRROP ,celldivincrrop )(
      CHFp_FRA1(divvel)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(divdir)
      ,CHFp_BOX(region) );

#define FORT_CELLDIVINCRROP FORTRAN_NAME( inlineCELLDIVINCRROP, inlineCELLDIVINCRROP)
#define FORTNT_CELLDIVINCRROP FORTRAN_NAME( CELLDIVINCRROP, celldivincrrop)

inline void FORTRAN_NAME(inlineCELLDIVINCRROP, inlineCELLDIVINCRROP)(
      CHFp_FRA1(divvel)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(divdir)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_CELLDIVINCRROP");
 FORTRAN_NAME( CELLDIVINCRROP ,celldivincrrop )(
      CHFt_FRA1(divvel)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(divdir)
      ,CHFt_BOX(region) );
}
#endif  // GUARDCELLDIVINCRROP 

#ifndef GUARDFACEDIVINCRROP 
#define GUARDFACEDIVINCRROP 
// Prototype for Fortran procedure facedivincrrop ...
//
void FORTRAN_NAME( FACEDIVINCRROP ,facedivincrrop )(
      CHFp_FRA1(divvel)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA(gradvel)
      ,CHFp_BOX(region)
      ,CHFp_BOX(centerbox)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(divdir)
      ,CHFp_CONST_INT(gradcomp) );

#define FORT_FACEDIVINCRROP FORTRAN_NAME( inlineFACEDIVINCRROP, inlineFACEDIVINCRROP)
#define FORTNT_FACEDIVINCRROP FORTRAN_NAME( FACEDIVINCRROP, facedivincrrop)

inline void FORTRAN_NAME(inlineFACEDIVINCRROP, inlineFACEDIVINCRROP)(
      CHFp_FRA1(divvel)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA(gradvel)
      ,CHFp_BOX(region)
      ,CHFp_BOX(centerbox)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(divdir)
      ,CHFp_CONST_INT(gradcomp) )
{
 CH_TIMELEAF("FORT_FACEDIVINCRROP");
 FORTRAN_NAME( FACEDIVINCRROP ,facedivincrrop )(
      CHFt_FRA1(divvel)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_CONST_FRA(gradvel)
      ,CHFt_BOX(region)
      ,CHFt_BOX(centerbox)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(divdir)
      ,CHFt_CONST_INT(gradcomp) );
}
#endif  // GUARDFACEDIVINCRROP 

#ifndef GUARDDECRINVLAMBDAROP 
#define GUARDDECRINVLAMBDAROP 
// Prototype for Fortran procedure decrinvlambdarop ...
//
void FORTRAN_NAME( DECRINVLAMBDAROP ,decrinvlambdarop )(
      CHFp_FRA(lambda)
      ,CHFp_FRA(eta)
      ,CHFp_BOX(box)
      ,CHFp_REAL(beta)
      ,CHFp_REAL(dx)
      ,CHFp_INT(idir) );

#define FORT_DECRINVLAMBDAROP FORTRAN_NAME( inlineDECRINVLAMBDAROP, inlineDECRINVLAMBDAROP)
#define FORTNT_DECRINVLAMBDAROP FORTRAN_NAME( DECRINVLAMBDAROP, decrinvlambdarop)

inline void FORTRAN_NAME(inlineDECRINVLAMBDAROP, inlineDECRINVLAMBDAROP)(
      CHFp_FRA(lambda)
      ,CHFp_FRA(eta)
      ,CHFp_BOX(box)
      ,CHFp_REAL(beta)
      ,CHFp_REAL(dx)
      ,CHFp_INT(idir) )
{
 CH_TIMELEAF("FORT_DECRINVLAMBDAROP");
 FORTRAN_NAME( DECRINVLAMBDAROP ,decrinvlambdarop )(
      CHFt_FRA(lambda)
      ,CHFt_FRA(eta)
      ,CHFt_BOX(box)
      ,CHFt_REAL(beta)
      ,CHFt_REAL(dx)
      ,CHFt_INT(idir) );
}
#endif  // GUARDDECRINVLAMBDAROP 

#ifndef GUARDINVERTLAMBDAROP 
#define GUARDINVERTLAMBDAROP 
// Prototype for Fortran procedure invertlambdarop ...
//
void FORTRAN_NAME( INVERTLAMBDAROP ,invertlambdarop )(
      CHFp_FRA(lambda)
      ,CHFp_REAL(safety)
      ,CHFp_BOX(box) );

#define FORT_INVERTLAMBDAROP FORTRAN_NAME( inlineINVERTLAMBDAROP, inlineINVERTLAMBDAROP)
#define FORTNT_INVERTLAMBDAROP FORTRAN_NAME( INVERTLAMBDAROP, invertlambdarop)

inline void FORTRAN_NAME(inlineINVERTLAMBDAROP, inlineINVERTLAMBDAROP)(
      CHFp_FRA(lambda)
      ,CHFp_REAL(safety)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_INVERTLAMBDAROP");
 FORTRAN_NAME( INVERTLAMBDAROP ,invertlambdarop )(
      CHFt_FRA(lambda)
      ,CHFt_REAL(safety)
      ,CHFt_BOX(box) );
}
#endif  // GUARDINVERTLAMBDAROP 

#ifndef GUARDGSRBROP 
#define GUARDGSRBROP 
// Prototype for Fortran procedure gsrbrop ...
//
void FORTRAN_NAME( GSRBROP ,gsrbrop )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(lphi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_BOX(coloredbox) );

#define FORT_GSRBROP FORTRAN_NAME( inlineGSRBROP, inlineGSRBROP)
#define FORTNT_GSRBROP FORTRAN_NAME( GSRBROP, gsrbrop)

inline void FORTRAN_NAME(inlineGSRBROP, inlineGSRBROP)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(lphi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_BOX(coloredbox) )
{
 CH_TIMELEAF("FORT_GSRBROP");
 FORTRAN_NAME( GSRBROP ,gsrbrop )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(lphi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA(lambda)
      ,CHFt_BOX(coloredbox) );
}
#endif  // GUARDGSRBROP 

#ifndef GUARDADDDIVFLUXDIRROP 
#define GUARDADDDIVFLUXDIRROP 
// Prototype for Fortran procedure ADDDIVFLUXDIRROP ...
//
void FORTRAN_NAME( ADDDIVFLUXDIRROP ,adddivfluxdirrop )(
      CHFp_FRA(lhs)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir) );

#define FORT_ADDDIVFLUXDIRROP FORTRAN_NAME( inlineADDDIVFLUXDIRROP, inlineADDDIVFLUXDIRROP)
#define FORTNT_ADDDIVFLUXDIRROP FORTRAN_NAME( ADDDIVFLUXDIRROP, adddivfluxdirrop)

inline void FORTRAN_NAME(inlineADDDIVFLUXDIRROP, inlineADDDIVFLUXDIRROP)(
      CHFp_FRA(lhs)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_ADDDIVFLUXDIRROP");
 FORTRAN_NAME( ADDDIVFLUXDIRROP ,adddivfluxdirrop )(
      CHFt_FRA(lhs)
      ,CHFt_CONST_FRA(flux)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDADDDIVFLUXDIRROP 

}

#endif
