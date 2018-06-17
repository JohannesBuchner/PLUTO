;+
;
; NAME:         FIND_SHOCK
;
; AUTHOR:       Andrea Mignone (mignone@to.astro.it)
;
; PURPOSE:      Find shocks in a multidimensional plane
;
; SYNTAX:       array = FIND_SHOCK(p, vx, vy, vz
;                                  [,eps_min=eps_min][,eps_max=eps_max])
;
; DESCRIPTION:  shockfind returns a multi-D array (with dimensionality
;               equal to the input arrays) whose elements are
;               numbers between 0 (no shock) and 1 (strong shock).
;
;               A shock is detected when div.V < 0
;               and |grad p|*dl_min/p > eps_min. 
;
;               The strength of the shock is controlled by
;
;               eps_min <  |grad p|/p < eps_max
;
;               This function requires PTOOLS.
;
; ARGUMENTS and KEYWORDS:
;
;              p:    an array giving the pressure distribution
;              vx:   an array giving the x-component of velocity
;              vy:   an array giving the y-component of velocity
;              vz:   an array giving the z-component of velocity
;              eps_min:  the minimum threshold for a shock to exist
;              eps_max:  the maximum strength of a shock (in units
;                        of |grad p|/p).
;
;
; LAST MODIFIED:   May 30, 2016 by A.Mignone (mignone@to.astro.it)
;
;-

FUNCTION FIND_SHOCK, prs, vx1, vx2, vx3, $
                     eps_min = eps_min,eps_max = eps_max

  COMMON PLUTO_GRID

  IF (NOT KEYWORD_SET(eps_min)) THEN eps_min = 0.33
  IF (NOT KEYWORD_SET(eps_max)) THEN eps_max = 5.0

  sa = size(prs)
  ndim = SIZE(prs,/N_DIMENSIONS)

; -----------------------------------------------------
;  Compute shock strenght in 2D
; -----------------------------------------------------

  IF (ndim EQ 2) THEN BEGIN
    sh    = FLTARR(nx1,nx2)

    dl_min = FLTARR(nx1,nx2)
    FOR j = 0,nx2-1 DO BEGIN
    FOR i = 0,nx1-1 DO BEGIN
      dl_min[i,j] = MIN([dx1[i], dx2[j]])
    ENDFOR
    ENDFOR

    divV  = PDIV (vx1,vx2)
    gradP = PGRAD(prs)
    dp    = sqrt(gradP[*,*,0]^2 + gradP[*,*,1]^2)*dl_min/prs
    sh    = (dp - eps_min)/(eps_max - eps_min)

    sh    = sh > 0.0
    sh    = sh < 1.0
    sh    = sh*(divV LT 0.0)
  ENDIF

; -----------------------------------------------------
;  Compute shock strenght in 3D
; -----------------------------------------------------

  IF (ndim EQ 3) THEN BEGIN
    sh    = FLTARR(nx1,nx2,nx3)

    dl_min = FLTARR(nx1,nx2,nx3)
    FOR k = 0,nx3-1 DO BEGIN
    FOR j = 0,nx2-1 DO BEGIN
    FOR i = 0,nx1-1 DO BEGIN
      dl_min[i,j,k] = MIN([dx1[i], dx2[j], dx3[k]])
    ENDFOR
    ENDFOR
    ENDFOR

    divV  = PDIV (vx1,vx2,vx3)
    gradP = PGRAD(prs)
    dp    = sqrt(gradP[*,*,*,0]^2 + gradP[*,*,*,1]^2 + gradP[*,*,*,2]^2)*dl_min/prs
    sh    = (dp - eps_min)/(eps_max - eps_min)

    sh    = sh > 0.0
    sh    = sh < 1.0
    sh    = sh*(divV LT 0.0)
  ENDIF

  RETURN,sh
END
