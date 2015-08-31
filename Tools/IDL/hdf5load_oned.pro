
;+
;
; NAME:      HDF5LOAD_ONED
;
; AUTHOR:    A. Mignone, O. Tesileanu, C. Zanni
;
; REQUIRES:  HDF5 support.
;
; PURPOSE:   Read a one-dimensional HDF5 file and store its content on a
;            non-uniform grid that glues together different blocks at the
;            local larger level of refinement available.
;            Data will be stored sequentially using array concatenation, starting
;            from the finest level to the base one.
;            It is normally by PLOAD
;
; SYNTAX:    HDF5_1DLOAD, nout, dir, level=level
;
; KEYWORDS:  silent
;
; LAST MODIFIED
;
;   March 24, 2015 by A. Mignone (mignone@ph.unito.it)
;
;-
PRO HDF5LOAD_ONED, nout, dir, level=level,silent=silent, x1range=x1range

  COMMON PLUTO_GRID
  COMMON PLUTO_VAR

; ----------------------------------------------
;  Check keywrods
; ----------------------------------------------

  IF (NOT KEYWORD_SET(level))  THEN level = 0
  IF (NOT KEYWORD_SET(silent)) THEN silent = 0

 
  qrho  = [0.0] & qpr    = [0.0]
  qtr1  = [0.0] & qtr2   = [0.0] & qtr3  = [0.0]
  qvx   = [0.0] & qvy    = [0.0] & qvz   = [0.0]
  qx_HI = [0.0] & qx_HII = [0.0] & qx_H2 = [0.0]
  qbx   = [0.0] & qby    = [0.0] & qbz   = [0.0]
  xx    = [0.0] 

  HDF5LOAD, nout, dir, level = level, silent=silent, x1range=x1range
  indx = WHERE(AMRBoxes EQ level)

  IF (n_elements(rho) GT 0) THEN qrho = [rho(indx)]
  IF (n_elements(vx1) GT 0) THEN qvx1 = [vx1(indx)]
  IF (n_elements(vx2) GT 0) THEN qvx2 = [vx2(indx)]
  IF (n_elements(vx3) GT 0) THEN qvx3 = [vx3(indx)]

  IF (n_elements(prs) GT 0) THEN qprs = [prs(indx)]

  IF (n_elements(tr1) GT 0) THEN qtr1 = [tr1(indx)]
  IF (n_elements(tr2) GT 0) THEN qtr2 = [tr2(indx)]
  IF (n_elements(tr3) GT 0) THEN qtr3 = [tr3(indx)]

  IF (n_elements(x_HI)  GT 0) THEN qx_HI  = [x_HI(indx)]
  IF (n_elements(x_HII) GT 0) THEN qx_HII = [x_HII(indx)]
  IF (n_elements(x_H2)  GT 0) THEN qx_H2  = [x_H2(indx)]

  IF (n_elements(bx1) GT 0) THEN qbx1 = [bx1(indx)]
  IF (n_elements(bx2) GT 0) THEN qbx2 = [bx2(indx)]
  IF (n_elements(bx3) GT 0) THEN qbx3 = [bx3(indx)]
 
  xx   = [x1(indx)]
  dx   = [dx1(indx)]

  AMRBoxesMax = AMRBoxes
  xmax        = x1
  FOR n = level-1,0,-1 DO BEGIN
    PRINT,"> HDF5LOAD_ONED: Building lev ",n
    HDF5LOAD, nout, dir, level = n, silent=silent, x1range=x1range
    indx = WHERE((AMRBoxes EQ n) AND REBIN(AMRBoxesmax,nx1) EQ n)
    IF (n_elements(rho) GT 0) THEN qrho = [qrho, rho(indx)]

    IF (n_elements(vx1) GT 0) THEN qvx1 = [qvx1, vx1(indx)]
    IF (n_elements(vx2) GT 0) THEN qvx2 = [qvx2, vx2(indx)]
    IF (n_elements(vx3) GT 0) THEN qvx3 = [qvx3, vx3(indx)]

    IF (n_elements(prs) GT 0) THEN qprs = [qprs, prs(indx)]

    IF (n_elements(tr1) GT 0) THEN qtr1 = [qtr1, tr1(indx)]
    IF (n_elements(tr2) GT 0) THEN qtr2 = [qtr2, tr2(indx)]
    IF (n_elements(tr3) GT 0) THEN qtr3 = [qtr3, tr3(indx)]

    IF (n_elements(x_HI)  GT 0) THEN qx_HI  = [qx_HI,  x_HI(indx)]
    IF (n_elements(x_HII) GT 0) THEN qx_HII = [qx_HII, x_HII(indx)]
    IF (n_elements(x_H2)  GT 0) THEN qx_H2  = [qx_H2,  x_H2(indx)]

    IF (n_elements(bx1) GT 0) THEN qbx1 = [qbx1, bx1(indx)]
    IF (n_elements(bx2) GT 0) THEN qbx2 = [qbx2, bx2(indx)]
    IF (n_elements(bx3) GT 0) THEN qbx3 = [qbx3, bx3(indx)]

    xx = [xx,x1(indx)]
    dx = [dx,dx1(indx)]
  ENDFOR

  isort = SORT(xx)
  x1  = xx[isort]
  dx1 = dx[isort]

  IF(n_elements(rho) GT 0) THEN rho = qrho[isort]

  IF(n_elements(vx1) GT 0) THEN vx1 = qvx1[isort]
  IF(n_elements(vx2) GT 0) THEN vx2 = qvx2[isort]
  IF(n_elements(vx3) GT 0) THEN vx3 = qvx3[isort]

  IF(n_elements(prs) GT 0) THEN prs = qprs[isort]

  IF(n_elements(tr1) GT 0) THEN tr1 = qtr1[isort]
  IF(n_elements(tr2) GT 0) THEN tr2 = qtr2[isort]
  IF(n_elements(tr3) GT 0) THEN tr3 = qtr3[isort]

  IF(n_elements(x_HI)  GT 0) THEN x_HI  = qx_HI[isort]
  IF(n_elements(x_HII) GT 0) THEN x_HII = qx_HII[isort]
  IF(n_elements(x_H2)  GT 0) THEN x_H2  = qx_H2[isort]

  IF(n_elements(bx1) GT 0) THEN bx1 = qbx1[isort]
  IF(n_elements(bx2) GT 0) THEN bx2 = qbx2[isort]
  IF(n_elements(bx3) GT 0) THEN bx3 = qbx3[isort]
  
END

