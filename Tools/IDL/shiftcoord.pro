;+
;
; NAME:      ShiftCoord
;
; AUTHOR:    Andrea Mignone (mignone@ph.unito.it)
;
; PURPOSE:   Shift the coordinate origin of AMR Data by an arbitrary amount.
;
; SYNTAX:    ShiftAMRData, offset
;
;
; ARGUMENTS:  offset is 1, 2 or 3D array specifying the coordinate shift(s).
;        
; KEYWORDS:
;
;      NONE
;
;  EXAMPLES:
;
;    * Example #1: Shift the origin by -0.5 in the x-direction:
;
;      IDL> ShiftAMRData, [-0.5,0]
;
;    * Example #2: Shift the origin in the point 1,3,10:
;
;      IDL> ShiftAMRData, [1, 3, 10]
;
;- 
PRO ShiftCoord, r0

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

 sr0 = size(r0, /dimensions)

 x0 = 0.0 & y0 = 0.0 & z0 = 0.0

 IF (sr0 EQ 1) THEN x0 = r0[0]
 IF (sr0 EQ 2) THEN BEGIN
   x0 = r0[0]
   y0 = r0[1]
 ENDIF 
 IF (sr0 EQ 3) THEN BEGIN
   x0 = r0[0]
   y0 = r0[1]
   z0 = r0[2]
 ENDIF 

 ; ---------------------------------------
 ;  get the number of refinement levels
 ; ---------------------------------------

 nlev = (size(AMRLevel))[1]
 FOR nl = 0,nlev-1 DO BEGIN
   level = AMRLevel[nl]
   FOR nb = 0, level.nbox-1 DO BEGIN
     box = level.box[nb]
     AMRLevel[nl].box[nb].x0 += x0
     AMRLevel[nl].box[nb].x1 += x0
     AMRLevel[nl].box[nb].y0 += y0
     AMRLevel[nl].box[nb].y1 += y0
     AMRLevel[nl].box[nb].z0 += z0
     AMRLevel[nl].box[nb].z1 += z0
   ENDFOR
 ENDFOR

; -------------------------------------
;   now shift grid coordinates
; -------------------------------------

  x1 += x0
  x2 += y0
  x3 += z0

END


