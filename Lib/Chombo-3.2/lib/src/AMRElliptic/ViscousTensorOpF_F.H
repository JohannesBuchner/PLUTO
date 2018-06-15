#ifndef _VISCOUSTENSOROPF_F_H_
#define _VISCOUSTENSOROPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDPROLONGVTOP 
#define GUARDPROLONGVTOP 
// Prototype for Fortran procedure prolongvtop ...
//
void FORTRAN_NAME( PROLONGVTOP ,prolongvtop )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(m) );

#define FORT_PROLONGVTOP FORTRAN_NAME( inlinePROLONGVTOP, inlinePROLONGVTOP)
#define FORTNT_PROLONGVTOP FORTRAN_NAME( PROLONGVTOP, prolongvtop)

inline void FORTRAN_NAME(inlinePROLONGVTOP, inlinePROLONGVTOP)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(m) )
{
 CH_TIMELEAF("FORT_PROLONGVTOP");
 FORTRAN_NAME( PROLONGVTOP ,prolongvtop )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(m) );
}
#endif  // GUARDPROLONGVTOP 

#ifndef GUARDRESTRICTRESVTOP 
#define GUARDRESTRICTRESVTOP 
// Prototype for Fortran procedure restrictresvtop ...
//
void FORTRAN_NAME( RESTRICTRESVTOP ,restrictresvtop )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(resfine)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(ncomp) );

#define FORT_RESTRICTRESVTOP FORTRAN_NAME( inlineRESTRICTRESVTOP, inlineRESTRICTRESVTOP)
#define FORTNT_RESTRICTRESVTOP FORTRAN_NAME( RESTRICTRESVTOP, restrictresvtop)

inline void FORTRAN_NAME(inlineRESTRICTRESVTOP, inlineRESTRICTRESVTOP)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(resfine)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(ncomp) )
{
 CH_TIMELEAF("FORT_RESTRICTRESVTOP");
 FORTRAN_NAME( RESTRICTRESVTOP ,restrictresvtop )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(resfine)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(ncomp) );
}
#endif  // GUARDRESTRICTRESVTOP 

#ifndef GUARDCELLGRADVTOP 
#define GUARDCELLGRADVTOP 
// Prototype for Fortran procedure cellgradvtop ...
//
void FORTRAN_NAME( CELLGRADVTOP ,cellgradvtop )(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(divdir) );

#define FORT_CELLGRADVTOP FORTRAN_NAME( inlineCELLGRADVTOP, inlineCELLGRADVTOP)
#define FORTNT_CELLGRADVTOP FORTRAN_NAME( CELLGRADVTOP, cellgradvtop)

inline void FORTRAN_NAME(inlineCELLGRADVTOP, inlineCELLGRADVTOP)(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(divdir) )
{
 CH_TIMELEAF("FORT_CELLGRADVTOP");
 FORTRAN_NAME( CELLGRADVTOP ,cellgradvtop )(
      CHFt_FRA1(grad)
      ,CHFt_CONST_FRA1(vel)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(divdir) );
}
#endif  // GUARDCELLGRADVTOP 

#ifndef GUARDADDGRADTOFLUXVTOP 
#define GUARDADDGRADTOFLUXVTOP 
// Prototype for Fortran procedure addgradtofluxvtop ...
//
void FORTRAN_NAME( ADDGRADTOFLUXVTOP ,addgradtofluxvtop )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA1(eta)
      ,CHFp_CONST_INT(fluxcomp)
      ,CHFp_CONST_FRA(grad)
      ,CHFp_CONST_INT(gradcomp)
      ,CHFp_CONST_INT(gradtran)
      ,CHFp_BOX(regionface) );

#define FORT_ADDGRADTOFLUXVTOP FORTRAN_NAME( inlineADDGRADTOFLUXVTOP, inlineADDGRADTOFLUXVTOP)
#define FORTNT_ADDGRADTOFLUXVTOP FORTRAN_NAME( ADDGRADTOFLUXVTOP, addgradtofluxvtop)

inline void FORTRAN_NAME(inlineADDGRADTOFLUXVTOP, inlineADDGRADTOFLUXVTOP)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA1(eta)
      ,CHFp_CONST_INT(fluxcomp)
      ,CHFp_CONST_FRA(grad)
      ,CHFp_CONST_INT(gradcomp)
      ,CHFp_CONST_INT(gradtran)
      ,CHFp_BOX(regionface) )
{
 CH_TIMELEAF("FORT_ADDGRADTOFLUXVTOP");
 FORTRAN_NAME( ADDGRADTOFLUXVTOP ,addgradtofluxvtop )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA1(eta)
      ,CHFt_CONST_INT(fluxcomp)
      ,CHFt_CONST_FRA(grad)
      ,CHFt_CONST_INT(gradcomp)
      ,CHFt_CONST_INT(gradtran)
      ,CHFt_BOX(regionface) );
}
#endif  // GUARDADDGRADTOFLUXVTOP 

#ifndef GUARDGETFACEGRADVTOP 
#define GUARDGETFACEGRADVTOP 
// Prototype for Fortran procedure getfacegradvtop ...
//
void FORTRAN_NAME( GETFACEGRADVTOP ,getfacegradvtop )(
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

#define FORT_GETFACEGRADVTOP FORTRAN_NAME( inlineGETFACEGRADVTOP, inlineGETFACEGRADVTOP)
#define FORTNT_GETFACEGRADVTOP FORTRAN_NAME( GETFACEGRADVTOP, getfacegradvtop)

inline void FORTRAN_NAME(inlineGETFACEGRADVTOP, inlineGETFACEGRADVTOP)(
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
 CH_TIMELEAF("FORT_GETFACEGRADVTOP");
 FORTRAN_NAME( GETFACEGRADVTOP ,getfacegradvtop )(
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
#endif  // GUARDGETFACEGRADVTOP 

#ifndef GUARDCELLDIVINCRVTOP 
#define GUARDCELLDIVINCRVTOP 
// Prototype for Fortran procedure celldivincrvtop ...
//
void FORTRAN_NAME( CELLDIVINCRVTOP ,celldivincrvtop )(
      CHFp_FRA1(divvel)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(divdir)
      ,CHFp_BOX(region) );

#define FORT_CELLDIVINCRVTOP FORTRAN_NAME( inlineCELLDIVINCRVTOP, inlineCELLDIVINCRVTOP)
#define FORTNT_CELLDIVINCRVTOP FORTRAN_NAME( CELLDIVINCRVTOP, celldivincrvtop)

inline void FORTRAN_NAME(inlineCELLDIVINCRVTOP, inlineCELLDIVINCRVTOP)(
      CHFp_FRA1(divvel)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(divdir)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_CELLDIVINCRVTOP");
 FORTRAN_NAME( CELLDIVINCRVTOP ,celldivincrvtop )(
      CHFt_FRA1(divvel)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(divdir)
      ,CHFt_BOX(region) );
}
#endif  // GUARDCELLDIVINCRVTOP 

#ifndef GUARDFACEDIVINCRVTOP 
#define GUARDFACEDIVINCRVTOP 
// Prototype for Fortran procedure facedivincrvtop ...
//
void FORTRAN_NAME( FACEDIVINCRVTOP ,facedivincrvtop )(
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

#define FORT_FACEDIVINCRVTOP FORTRAN_NAME( inlineFACEDIVINCRVTOP, inlineFACEDIVINCRVTOP)
#define FORTNT_FACEDIVINCRVTOP FORTRAN_NAME( FACEDIVINCRVTOP, facedivincrvtop)

inline void FORTRAN_NAME(inlineFACEDIVINCRVTOP, inlineFACEDIVINCRVTOP)(
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
 CH_TIMELEAF("FORT_FACEDIVINCRVTOP");
 FORTRAN_NAME( FACEDIVINCRVTOP ,facedivincrvtop )(
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
#endif  // GUARDFACEDIVINCRVTOP 

#ifndef GUARDDECRINVRELCOEFVTOP 
#define GUARDDECRINVRELCOEFVTOP 
// Prototype for Fortran procedure decrinvrelcoefvtop ...
//
void FORTRAN_NAME( DECRINVRELCOEFVTOP ,decrinvrelcoefvtop )(
      CHFp_FRA(relcoef)
      ,CHFp_FRA(eta)
      ,CHFp_FRA(lambda)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_REAL(dx)
      ,CHFp_INT(idir)
      ,CHFp_INT(ncomp) );

#define FORT_DECRINVRELCOEFVTOP FORTRAN_NAME( inlineDECRINVRELCOEFVTOP, inlineDECRINVRELCOEFVTOP)
#define FORTNT_DECRINVRELCOEFVTOP FORTRAN_NAME( DECRINVRELCOEFVTOP, decrinvrelcoefvtop)

inline void FORTRAN_NAME(inlineDECRINVRELCOEFVTOP, inlineDECRINVRELCOEFVTOP)(
      CHFp_FRA(relcoef)
      ,CHFp_FRA(eta)
      ,CHFp_FRA(lambda)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_REAL(dx)
      ,CHFp_INT(idir)
      ,CHFp_INT(ncomp) )
{
 CH_TIMELEAF("FORT_DECRINVRELCOEFVTOP");
 FORTRAN_NAME( DECRINVRELCOEFVTOP ,decrinvrelcoefvtop )(
      CHFt_FRA(relcoef)
      ,CHFt_FRA(eta)
      ,CHFt_FRA(lambda)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(box)
      ,CHFt_REAL(dx)
      ,CHFt_INT(idir)
      ,CHFt_INT(ncomp) );
}
#endif  // GUARDDECRINVRELCOEFVTOP 

#ifndef GUARDINITIALIZERELAXCOEF 
#define GUARDINITIALIZERELAXCOEF 
// Prototype for Fortran procedure initializerelaxcoef ...
//
void FORTRAN_NAME( INITIALIZERELAXCOEF ,initializerelaxcoef )(
      CHFp_FRA(relcoef)
      ,CHFp_FRA1(acoef)
      ,CHFp_REAL(alpha)
      ,CHFp_BOX(box)
      ,CHFp_INT(ncomp) );

#define FORT_INITIALIZERELAXCOEF FORTRAN_NAME( inlineINITIALIZERELAXCOEF, inlineINITIALIZERELAXCOEF)
#define FORTNT_INITIALIZERELAXCOEF FORTRAN_NAME( INITIALIZERELAXCOEF, initializerelaxcoef)

inline void FORTRAN_NAME(inlineINITIALIZERELAXCOEF, inlineINITIALIZERELAXCOEF)(
      CHFp_FRA(relcoef)
      ,CHFp_FRA1(acoef)
      ,CHFp_REAL(alpha)
      ,CHFp_BOX(box)
      ,CHFp_INT(ncomp) )
{
 CH_TIMELEAF("FORT_INITIALIZERELAXCOEF");
 FORTRAN_NAME( INITIALIZERELAXCOEF ,initializerelaxcoef )(
      CHFt_FRA(relcoef)
      ,CHFt_FRA1(acoef)
      ,CHFt_REAL(alpha)
      ,CHFt_BOX(box)
      ,CHFt_INT(ncomp) );
}
#endif  // GUARDINITIALIZERELAXCOEF 

#ifndef GUARDINVERTLAMBDAVTOP 
#define GUARDINVERTLAMBDAVTOP 
// Prototype for Fortran procedure invertlambdavtop ...
//
void FORTRAN_NAME( INVERTLAMBDAVTOP ,invertlambdavtop )(
      CHFp_FRA(lambda)
      ,CHFp_REAL(safety)
      ,CHFp_BOX(box)
      ,CHFp_INT(ncomp) );

#define FORT_INVERTLAMBDAVTOP FORTRAN_NAME( inlineINVERTLAMBDAVTOP, inlineINVERTLAMBDAVTOP)
#define FORTNT_INVERTLAMBDAVTOP FORTRAN_NAME( INVERTLAMBDAVTOP, invertlambdavtop)

inline void FORTRAN_NAME(inlineINVERTLAMBDAVTOP, inlineINVERTLAMBDAVTOP)(
      CHFp_FRA(lambda)
      ,CHFp_REAL(safety)
      ,CHFp_BOX(box)
      ,CHFp_INT(ncomp) )
{
 CH_TIMELEAF("FORT_INVERTLAMBDAVTOP");
 FORTRAN_NAME( INVERTLAMBDAVTOP ,invertlambdavtop )(
      CHFt_FRA(lambda)
      ,CHFt_REAL(safety)
      ,CHFt_BOX(box)
      ,CHFt_INT(ncomp) );
}
#endif  // GUARDINVERTLAMBDAVTOP 

#ifndef GUARDGSRBVTOP 
#define GUARDGSRBVTOP 
// Prototype for Fortran procedure gsrbvtop ...
//
void FORTRAN_NAME( GSRBVTOP ,gsrbvtop )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(lphi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_BOX(coloredbox)
      ,CHFp_INT(ncomp) );

#define FORT_GSRBVTOP FORTRAN_NAME( inlineGSRBVTOP, inlineGSRBVTOP)
#define FORTNT_GSRBVTOP FORTRAN_NAME( GSRBVTOP, gsrbvtop)

inline void FORTRAN_NAME(inlineGSRBVTOP, inlineGSRBVTOP)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(lphi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_BOX(coloredbox)
      ,CHFp_INT(ncomp) )
{
 CH_TIMELEAF("FORT_GSRBVTOP");
 FORTRAN_NAME( GSRBVTOP ,gsrbvtop )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(lphi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA(lambda)
      ,CHFt_BOX(coloredbox)
      ,CHFt_INT(ncomp) );
}
#endif  // GUARDGSRBVTOP 

#ifndef GUARDADDDIVFLUXDIRVTOP 
#define GUARDADDDIVFLUXDIRVTOP 
// Prototype for Fortran procedure adddivfluxdirvtop ...
//
void FORTRAN_NAME( ADDDIVFLUXDIRVTOP ,adddivfluxdirvtop )(
      CHFp_FRA(lhs)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_CONST_INT(facedir) );

#define FORT_ADDDIVFLUXDIRVTOP FORTRAN_NAME( inlineADDDIVFLUXDIRVTOP, inlineADDDIVFLUXDIRVTOP)
#define FORTNT_ADDDIVFLUXDIRVTOP FORTRAN_NAME( ADDDIVFLUXDIRVTOP, adddivfluxdirvtop)

inline void FORTRAN_NAME(inlineADDDIVFLUXDIRVTOP, inlineADDDIVFLUXDIRVTOP)(
      CHFp_FRA(lhs)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_ADDDIVFLUXDIRVTOP");
 FORTRAN_NAME( ADDDIVFLUXDIRVTOP ,adddivfluxdirvtop )(
      CHFt_FRA(lhs)
      ,CHFt_CONST_FRA(flux)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(ncomp)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDADDDIVFLUXDIRVTOP 

#ifndef GUARDLINEAREXTRAPVTOP 
#define GUARDLINEAREXTRAPVTOP 
// Prototype for Fortran procedure linearextrapvtop ...
//
void FORTRAN_NAME( LINEAREXTRAPVTOP ,linearextrapvtop )(
      CHFp_FRA(phi)
      ,CHFp_BOX(ghostBox)
      ,CHFp_INT(dir)
      ,CHFp_INT(hiLo) );

#define FORT_LINEAREXTRAPVTOP FORTRAN_NAME( inlineLINEAREXTRAPVTOP, inlineLINEAREXTRAPVTOP)
#define FORTNT_LINEAREXTRAPVTOP FORTRAN_NAME( LINEAREXTRAPVTOP, linearextrapvtop)

inline void FORTRAN_NAME(inlineLINEAREXTRAPVTOP, inlineLINEAREXTRAPVTOP)(
      CHFp_FRA(phi)
      ,CHFp_BOX(ghostBox)
      ,CHFp_INT(dir)
      ,CHFp_INT(hiLo) )
{
 CH_TIMELEAF("FORT_LINEAREXTRAPVTOP");
 FORTRAN_NAME( LINEAREXTRAPVTOP ,linearextrapvtop )(
      CHFt_FRA(phi)
      ,CHFt_BOX(ghostBox)
      ,CHFt_INT(dir)
      ,CHFt_INT(hiLo) );
}
#endif  // GUARDLINEAREXTRAPVTOP 

#ifndef GUARDSLOPESVTOP 
#define GUARDSLOPESVTOP 
// Prototype for Fortran procedure SLOPESVTOP ...
//
void FORTRAN_NAME( SLOPESVTOP ,slopesvtop )(
      CHFp_FRA(slopes)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(cBox)
      ,CHFp_INT(dir) );

#define FORT_SLOPESVTOP FORTRAN_NAME( inlineSLOPESVTOP, inlineSLOPESVTOP)
#define FORTNT_SLOPESVTOP FORTRAN_NAME( SLOPESVTOP, slopesvtop)

inline void FORTRAN_NAME(inlineSLOPESVTOP, inlineSLOPESVTOP)(
      CHFp_FRA(slopes)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(cBox)
      ,CHFp_INT(dir) )
{
 CH_TIMELEAF("FORT_SLOPESVTOP");
 FORTRAN_NAME( SLOPESVTOP ,slopesvtop )(
      CHFt_FRA(slopes)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(cBox)
      ,CHFt_INT(dir) );
}
#endif  // GUARDSLOPESVTOP 

}

#endif
