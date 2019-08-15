; --------------------------------
;  make the userdef symbol filled
;  circle
; --------------------------------

PRO SYMCIRCLE,fill=fill
  npt  = 17
  xsym = fltarr(npt+1)
  ysym = fltarr(npt+1)

  FOR i = 0,npt DO BEGIN
    xsym(i) = sin(2.0*!PI*i/(npt*1.0))
    ysym(i) = cos(2.0*!PI*i/(npt*1.0))
  ENDFOR
  IF (KEYWORD_SET(FILL)) THEN BEGIN
    usersym,xsym,ysym,/fill 
  ENDIF ELSE BEGIN
    usersym,xsym,ysym
  ENDELSE

END
