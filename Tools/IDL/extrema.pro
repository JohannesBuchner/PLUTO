;+
;
; NAME:      EXTREMA
;
; AUTHOR:    Andrea Mignone
;
; SYNTAX:    indx = EXTREMA(f)
;
; PURPOSE:   Find maxima in a 1D array f[].
;            On output, indx is an integer array whose elements are the
;            indices of the extrema of f.
;            If no maximum is found, return -1.
;
; ARGUMENTS:
;
;   f        a 1D array
;
; KEYWORD:   none
;
; LAST MODIFIED:   Dec 21, 2016 by A. Mignone 
;
;-
FUNCTION EXTREMA, f

  f   = REFORM(f)
  sf  = SIZE(f)
  n   = sf(1)
  eps = 0.e-8

  imax = INTARR(n)
  fmax = FLTARR(n)

  nmax = 0; number of maxima counter
  FOR i = 1, n - 2 DO BEGIN
    IF (f(i) GT f(i-1)) THEN BEGIN
      FOR k = i + 1, n - 1 do begin
        IF (f[k] GT f[i]-eps) THEN BREAK
        IF (f[k] LE f[i]) THEN BEGIN
          i0   = (k - 1 + i)/2
          imax[nmax] = i0
          fmax[nmax] = f[i0]
;          print, "max found at", imax[nmax], fmax[nmax]
          nmax = nmax + 1
          BREAK
        ENDIF
      ENDFOR
    ENDIF
  ENDFOR
  
  IF (nmax EQ 0) THEN RETURN,-1
  RETURN, imax[0:nmax-1]
;
;  Sort maxima
;

 ;print, "N elements", n
 ;print, "N max", n_max
 IF (n_max GE 0) THEN BEGIN

   is = imax(reverse(sort(fmax(0:n_max))))
;   fs   = fmax(reverse(sort(fmax(0:n_max))))
 ENDIF

 
 IF (n_max EQ -1) THEN return,-1
 RETURN,is
END
