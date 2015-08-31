;+
;
; NAME:      extrema
;
; AUTHOR:    Andrea Mignone
;
; DATE:      April 21, 2004
;
; PURPOSE:   Find maxima/minima in a function f. 
;            f is a 1-D vector.
;
; SYNOPSIS:  Result = extrema(f). On output
;            result is an integer array whose
;            elements are the indexes k of the
;            extrema of f.  
;
; KEYWORD:   none
;
;
;-
FUNCTION extrema, f

  sf  = SIZE(f)
  n   = sf(1)
  eps = 1.e-6

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
          nmax = nmax + 1
          ;print, "max found at", imax[nmax], fmax[nmax]
          BREAK
        ENDIF
      ENDFOR
    ENDIF
  ENDFOR

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
