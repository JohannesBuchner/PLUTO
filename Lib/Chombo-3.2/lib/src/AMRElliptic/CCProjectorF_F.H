#ifndef _CCPROJECTORF_F_H_
#define _CCPROJECTORF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCCPEXTRAPTODOMFACE 
#define GUARDCCPEXTRAPTODOMFACE 
// Prototype for Fortran procedure ccpextraptodomface ...
//
void FORTRAN_NAME( CCPEXTRAPTODOMFACE ,ccpextraptodomface )(
      CHFp_FRA1(facevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(ishift)
      ,CHFp_BOX(facebox) );

#define FORT_CCPEXTRAPTODOMFACE FORTRAN_NAME( inlineCCPEXTRAPTODOMFACE, inlineCCPEXTRAPTODOMFACE)
#define FORTNT_CCPEXTRAPTODOMFACE FORTRAN_NAME( CCPEXTRAPTODOMFACE, ccpextraptodomface)

inline void FORTRAN_NAME(inlineCCPEXTRAPTODOMFACE, inlineCCPEXTRAPTODOMFACE)(
      CHFp_FRA1(facevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(ishift)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_CCPEXTRAPTODOMFACE");
 FORTRAN_NAME( CCPEXTRAPTODOMFACE ,ccpextraptodomface )(
      CHFt_FRA1(facevel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(ishift)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDCCPEXTRAPTODOMFACE 

#ifndef GUARDCCPAVECELLTOFACE 
#define GUARDCCPAVECELLTOFACE 
// Prototype for Fortran procedure ccpavecelltoface ...
//
void FORTRAN_NAME( CCPAVECELLTOFACE ,ccpavecelltoface )(
      CHFp_FRA1(facevel)
      ,CHFp_CONST_FRA1(cellvel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox) );

#define FORT_CCPAVECELLTOFACE FORTRAN_NAME( inlineCCPAVECELLTOFACE, inlineCCPAVECELLTOFACE)
#define FORTNT_CCPAVECELLTOFACE FORTRAN_NAME( CCPAVECELLTOFACE, ccpavecelltoface)

inline void FORTRAN_NAME(inlineCCPAVECELLTOFACE, inlineCCPAVECELLTOFACE)(
      CHFp_FRA1(facevel)
      ,CHFp_CONST_FRA1(cellvel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_CCPAVECELLTOFACE");
 FORTRAN_NAME( CCPAVECELLTOFACE ,ccpavecelltoface )(
      CHFt_FRA1(facevel)
      ,CHFt_CONST_FRA1(cellvel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDCCPAVECELLTOFACE 

#ifndef GUARDCCPAVEFACETOCELL 
#define GUARDCCPAVEFACETOCELL 
// Prototype for Fortran procedure ccpavefacetocell ...
//
void FORTRAN_NAME( CCPAVEFACETOCELL ,ccpavefacetocell )(
      CHFp_FRA1(cellgrad)
      ,CHFp_CONST_FRA1(facegrad)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(cellbox) );

#define FORT_CCPAVEFACETOCELL FORTRAN_NAME( inlineCCPAVEFACETOCELL, inlineCCPAVEFACETOCELL)
#define FORTNT_CCPAVEFACETOCELL FORTRAN_NAME( CCPAVEFACETOCELL, ccpavefacetocell)

inline void FORTRAN_NAME(inlineCCPAVEFACETOCELL, inlineCCPAVEFACETOCELL)(
      CHFp_FRA1(cellgrad)
      ,CHFp_CONST_FRA1(facegrad)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(cellbox) )
{
 CH_TIMELEAF("FORT_CCPAVEFACETOCELL");
 FORTRAN_NAME( CCPAVEFACETOCELL ,ccpavefacetocell )(
      CHFt_FRA1(cellgrad)
      ,CHFt_CONST_FRA1(facegrad)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(cellbox) );
}
#endif  // GUARDCCPAVEFACETOCELL 

#ifndef GUARDCCPCELLGRADFROMFACEDATA 
#define GUARDCCPCELLGRADFROMFACEDATA 
// Prototype for Fortran procedure ccpcellgradfromfacedata ...
//
void FORTRAN_NAME( CCPCELLGRADFROMFACEDATA ,ccpcellgradfromfacedata )(
      CHFp_FRA1(cellgrad)
      ,CHFp_CONST_FRA1(facedata)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(cellbox) );

#define FORT_CCPCELLGRADFROMFACEDATA FORTRAN_NAME( inlineCCPCELLGRADFROMFACEDATA, inlineCCPCELLGRADFROMFACEDATA)
#define FORTNT_CCPCELLGRADFROMFACEDATA FORTRAN_NAME( CCPCELLGRADFROMFACEDATA, ccpcellgradfromfacedata)

inline void FORTRAN_NAME(inlineCCPCELLGRADFROMFACEDATA, inlineCCPCELLGRADFROMFACEDATA)(
      CHFp_FRA1(cellgrad)
      ,CHFp_CONST_FRA1(facedata)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(cellbox) )
{
 CH_TIMELEAF("FORT_CCPCELLGRADFROMFACEDATA");
 FORTRAN_NAME( CCPCELLGRADFROMFACEDATA ,ccpcellgradfromfacedata )(
      CHFt_FRA1(cellgrad)
      ,CHFt_CONST_FRA1(facedata)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(cellbox) );
}
#endif  // GUARDCCPCELLGRADFROMFACEDATA 

}

#endif
