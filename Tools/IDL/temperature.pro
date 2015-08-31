;+
;
; NAME:        TEMPERATURE
;
; AUTHOR:      A. Mignone (mignone@ph.unito.it)
;
; SYNTAX:      Tmp = TEMPERATURE(Unit_Vel=Unit_Vel)
;
; PURPOSE:     Compute and return the temperature (in K) from PLUTO Data files.
;              Temperature is calculated as T = p/rho * mu(X) * V0*mu/kB,
;              where mu(x) is the mean molecular weight (re-implemented here),
;              mu is the mean molecular weight, kB is the Boltzmann constant and
;              V0 is the unit velocity (in c.g.s units).
;              It employs pressure, density and chemical composition (when 
;              SNEq or MINEq) are used.
;
; ARGUMENTS:
;
;   Unit_Vel   The unit velocity (in cm/s) used by PLUTO
;
; KEYWORDS:
;
; LAST MODIFIED:   July 12, 2014 by A. Mignone (mignone@ph.unito.it)
;
;-
FUNCTION TEMPERATURE,Unit_Vel = Unit_vel

 COMMON PLUTO_VAR
 COMMON PLUTO_GRID

 K = 1.66053886e-24/1.3806506e-16*Unit_Vel^2

 IF (size(fHI,/N_DIMENSIONS) GT 0) THEN BEGIN ; MINEq is being used
   
   IF (size(fCII,/N_DIMENSIONS)  LE 0) THEN fCII = 0.0
   IF (size(fCIII,/N_DIMENSIONS) LE 0) THEN fCIII = 0.0
   if (size(fCIV,/N_DIMENSIONS)  LE 0) THEN fCIV = 0.0
   IF (size(fCV,/N_DIMENSIONS)   LE 0) THEN fCV = 0.0
   IF (size(fNII,/N_DIMENSIONS)  LE 0) THEN fNII = 0.0
   IF (size(fNIII,/N_DIMENSIONS) LE 0) THEN fNIII = 0.0
   IF (size(fNIV,/N_DIMENSIONS)  LE 0) THEN fNIV = 0.0
   IF (size(fNV,/N_DIMENSIONS)   LE 0) THEN fNV = 0.0
   IF (size(fOII,/N_DIMENSIONS)  LE 0) THEN fOII = 0.0
   IF (size(fOIII,/N_DIMENSIONS) LE 0) THEN fOIII = 0.0
   IF (size(fOIV,/N_DIMENSIONS)  LE 0) THEN fOIV = 0.0
   IF (size(fOV,/N_DIMENSIONS)   LE 0) THEN fOV = 0.0
   IF (size(fNeII,/N_DIMENSIONS)  LE 0) THEN fNeII = 0.0
   if (size(fNeIII,/N_DIMENSIONS) LE 0) THEN fNeIII = 0.0
   IF (size(fNeIV,/N_DIMENSIONS)  LE 0) THEN fNeIV = 0.0
   IF (size(fNeV,/N_DIMENSIONS)   LE 0) THEN fNeV = 0.0
   IF (size(fSII,/N_DIMENSIONS)   LE 0) THEN fSII = 0.0
   IF (size(fSIII,/N_DIMENSIONS)  LE 0) THEN fSIII = 0.0
   IF (size(fSIV,/N_DIMENSIONS)   LE 0) THEN fSIV = 0.0
   IF (size(fSV,/N_DIMENSIONS)    LE 0) THEN fSV = 0.0

   elem_ab   = [ 0.93 , 0.074, 3.e-4, 5.e-5, 4.0e-4, 7.0e-5, 1.5e-5, 2.69e-5 ]
   elem_mass = [ 1.007, 4.002, 12.01, 14.01, 15.99 , 20.18 , 32.07, 55.845 ]
   elem_part = [ 0, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7 ]
   rad_rec_z = [ 1., 1., 2., 1., 2., 3., 4., 5.,$
                       1., 2., 3., 4., 5., 1., 2., 3., 4., 5.,$
                       1., 2., 3., 4., 5., 1., 2., 3., 4., 5., 1., 2., 3.];

   ; -- H

   num = 1.007*0.93*fHI
   den = 0.93*fHI

   ; -- He

   num += 4.002*0.074*(fHeI+fHeII)
   den += 0.074*(fHeI + 2.0*fHeII)

   ; -- C

   num += 12.01*3.e-4*(fCI+fCII+fCIII+fCIV+fCV)
   den += 3.e-4*(fCI + 2.0*fCII + 3.0*fCIII + 4.0*fCIV + 5.0*fCV)

   ; -- N

   num += 14.01*5.e-5*(fNI+fNII+fNIII+fNIV+fNV)
   den += 5.e-5*(fNI + 2.0*fNII + 3.0*fNIII + 4.0*fNIV + 5.0*fNV)

   ; -- O

   num += 15.99*4.e-4*(fOI+fOII+fOIII+fOIV+fOV)
   den += 4.e-4*(fOI + 2.0*fOII + 3.0*fOIII + 4.0*fOIV + 5.0*fOV)

   ; -- Ne

   num += 20.18*7.e-5*(fNeI+fNeII+fNeIII+fNeIV+fNeV)
   den += 7.e-5*(fNeI + 2.0*fNeII + 3.0*fNeIII + 4.0*fNeIV + 5.0*fNeV)

   ; -- S

   num += 32.07*1.5e-5*(fSI+fSII+fSIII+fSIV+fSV)
   den += 7.e-5*(fSI + 2.0*fSII + 3.0*fSIII + 4.0*fSIV + 5.0*fSV)

  ; --------------------------------------------------
  ;      add now contribution from ionized H  
  ; --------------------------------------------------  

   num += 1.007*0.93*(1.0 - fHI);
   den += 0.93*(1.0 - fHI)*2.

   mu = num/den;
   T = K*prs/rho*mu
   RETURN,T
 ENDIF

 IF (size(fneut,/n_el) NE 0) THEN BEGIN  ; SNEq is being used
   fZ   =  1.e-3
   fHe  = 0.082  
   A_Z  = 30.0
   A_He = 4.004
   A_H  = 1.008 

   mu = (A_H + fHe*A_He + fZ*A_Z)/(2.0 + fHe + 2.0*fZ - fneut)
   T = K*prs/rho*mu
   RETURN,T
 ENDIF ELSE BEGIN
   mu = 1.237190
   RETURN, K*prs/rho*mu
 ENDELSE
 
 PRINT,"! Wrong network"
 STOP
END
