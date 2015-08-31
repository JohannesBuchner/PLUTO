;+
; NAME:      PLOT_ASCII
;
; AUTHOR:    Andrea Mignone (mignone@ph.unito.it)
;
; PURPOSE:   plot selected columns from an ascii data files.
;
; SYNTAX:  PLOT_ASCII,filename,nx,ny[,/OPLOT][,_EXTRA=extra]
;
; ARGUMENTS:
;
;   filename       ascii data file
; 
;   nx,ny          the abscissa and ordinata column index to be plotted
;
; KEYWORDS:
;
;   /OPLOT       overplot without erasing
;
; LAST MODIFIED: March 25, 2014  
;
; -
PRO PLOT_ASCII,filename,nx,ny, _EXTRA = extra, oplot=oplot

  IF (FILE_TEST(filename)) THEN BEGIN
    Q = READ_ASCII(filename,comm="#")
  ENDIF ELSE BEGIN
    PRINT,"! file "+filename+" does not exist"
    RETURN
  ENDELSE

  Q = Q.(0)
  IF (KEYWORD_SET(oplot)) THEN BEGIN
    OPLOT,Q[nx,*], Q[ny,*],_EXTRA=extra
  ENDIF ELSE BEGIN
    PLOT,Q[nx,*], Q[ny,*],_EXTRA=extra
  ENDELSE

END

