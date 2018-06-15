#ifndef _AMRNODEOPF_F_H_
#define _AMRNODEOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDNODEOPLAP 
#define GUARDNODEOPLAP 
// Prototype for Fortran procedure NODEOPLAP ...
//
void FORTRAN_NAME( NODEOPLAP ,nodeoplap )(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_NODEOPLAP FORTRAN_NAME( inlineNODEOPLAP, inlineNODEOPLAP)
#define FORTNT_NODEOPLAP FORTRAN_NAME( NODEOPLAP, nodeoplap)

inline void FORTRAN_NAME(inlineNODEOPLAP, inlineNODEOPLAP)(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_NODEOPLAP");
 FORTRAN_NAME( NODEOPLAP ,nodeoplap )(
      CHFt_FRA(lofphi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDNODEOPLAP 

#ifndef GUARDNODEOPLAPPOINT 
#define GUARDNODEOPLAPPOINT 
// Prototype for Fortran procedure NODEOPLAPPOINT ...
//
void FORTRAN_NAME( NODEOPLAPPOINT ,nodeoplappoint )(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_INTVECT(pt)
      ,CHFp_CONST_REAL(dx) );

#define FORT_NODEOPLAPPOINT FORTRAN_NAME( inlineNODEOPLAPPOINT, inlineNODEOPLAPPOINT)
#define FORTNT_NODEOPLAPPOINT FORTRAN_NAME( NODEOPLAPPOINT, nodeoplappoint)

inline void FORTRAN_NAME(inlineNODEOPLAPPOINT, inlineNODEOPLAPPOINT)(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_INTVECT(pt)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_NODEOPLAPPOINT");
 FORTRAN_NAME( NODEOPLAPPOINT ,nodeoplappoint )(
      CHFt_FRA(lofphi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_INTVECT(pt)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDNODEOPLAPPOINT 

#ifndef GUARDNODEGRAD 
#define GUARDNODEGRAD 
// Prototype for Fortran procedure NODEGRAD ...
//
void FORTRAN_NAME( NODEGRAD ,nodegrad )(
      CHFp_FRA(grdphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_NODEGRAD FORTRAN_NAME( inlineNODEGRAD, inlineNODEGRAD)
#define FORTNT_NODEGRAD FORTRAN_NAME( NODEGRAD, nodegrad)

inline void FORTRAN_NAME(inlineNODEGRAD, inlineNODEGRAD)(
      CHFp_FRA(grdphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_NODEGRAD");
 FORTRAN_NAME( NODEGRAD ,nodegrad )(
      CHFt_FRA(grdphi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDNODEGRAD 

#ifndef GUARDNODEGRADPOINT 
#define GUARDNODEGRADPOINT 
// Prototype for Fortran procedure NODEGRADPOINT ...
//
void FORTRAN_NAME( NODEGRADPOINT ,nodegradpoint )(
      CHFp_FRA(grdphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_INTVECT(pt)
      ,CHFp_CONST_REAL(dx) );

#define FORT_NODEGRADPOINT FORTRAN_NAME( inlineNODEGRADPOINT, inlineNODEGRADPOINT)
#define FORTNT_NODEGRADPOINT FORTRAN_NAME( NODEGRADPOINT, nodegradpoint)

inline void FORTRAN_NAME(inlineNODEGRADPOINT, inlineNODEGRADPOINT)(
      CHFp_FRA(grdphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_INTVECT(pt)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_NODEGRADPOINT");
 FORTRAN_NAME( NODEGRADPOINT ,nodegradpoint )(
      CHFt_FRA(grdphi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_INTVECT(pt)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDNODEGRADPOINT 

#ifndef GUARDNODEGSRBLEVELLAP 
#define GUARDNODEGSRBLEVELLAP 
// Prototype for Fortran procedure NODEGSRBLEVELLAP ...
//
void FORTRAN_NAME( NODEGSRBLEVELLAP ,nodegsrblevellap )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_NODEGSRBLEVELLAP FORTRAN_NAME( inlineNODEGSRBLEVELLAP, inlineNODEGSRBLEVELLAP)
#define FORTNT_NODEGSRBLEVELLAP FORTRAN_NAME( NODEGSRBLEVELLAP, nodegsrblevellap)

inline void FORTRAN_NAME(inlineNODEGSRBLEVELLAP, inlineNODEGSRBLEVELLAP)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_NODEGSRBLEVELLAP");
 FORTRAN_NAME( NODEGSRBLEVELLAP ,nodegsrblevellap )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDNODEGSRBLEVELLAP 

}

#endif
