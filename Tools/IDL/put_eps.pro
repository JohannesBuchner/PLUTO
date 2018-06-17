;+
;
;  NAME:      put_eps
;
;  AUTHOR:    Andrea Mignone
;
;  PURPOSE:   Put a single encapsulated postscript figure in
;             a multi-plot file.
;
;  SYNTAX:    put_eps, data, x, y, pos
;             [,xtitle=xtitle][,ytitle=ytitle][,title=title][,charsize=charsize]
;             [,dmax = dmax][,dmin = dmin][,/vbar][,/hbar][,color=color]
;             [,cbwidth=cbwidth][,cbcharsize=cbcharsize][,cbformat=cbformat]
;             [,/noxticks][,/noyticks][,xrange=xrange][,yrange=yrange]
;
;             where data is a two-dimensional array, x and y
;             are the coordinates, pos=[x0,y0,x1,y1] is the position
;             of the figure in the multiplot file,
;             calculated with set_multi_plot_pos.
;
;  KEYWORDS:
;
;
;     xtitle = the string name for the x-axis
;
;     ytitle = the string name for the y-axis
;
;     title  = the string name for the title
;
;     dmax   = a floating point number giving the maximum value
;              data value to which the image should be scaled to.
;
;     dmin   = a floating point number giving the minimum value
;              data value to which the image should be scaled to.
;
;     vbar   = set this keyword to put a vertical colorbar to the right
;              of the image.
;
;     hbar   = set this keyword to put a horizontal colorbar to
;              the top of the image.
;
;     color   = an integer giving the color (in the current color table)
;               of the axis.
;
;     charsize = size of characters used for colorbar
;
;     noxticks = suppress x-tick marks
;
;     noyticks = suppress y-tick marks
;
;     sample   = increase resolution by combining sub-sampling and
;                interpolation.
;
;     xrange = a two element vector giving the range in the x direction
;
;     yrange = a two element vector giving the range in the y direction
;
;  COLORBAR CONTROL KEYWORDS:
;
;     cbwidth    = the horizontal extent of the colorbar in units
;                  of the image (only if /vbar or /hbar is given).
;
;     cbcharsize = size of characters used for colorbar
;
;     cbdiv      = the number of divisions of the colorbar.
;
;     cbformat   = the string format used in the colorbar
;
;     cbposition = lower left position of the colorbar
;
;  LAST MODIFIED:  Oct 31, 2017 by A. Mignone (mignone@ph.unito.it)
;
;-

PRO put_eps, a0, x0, y0, pos, $
             xtitle=xtitle,ytitle=ytitle,title=title,$
             dmax = dmax,dmin = dmin, vbar=vbar, hbar=hbar, color=color,$
             cbdiv=cbdiv, cbwidth = cbwidth, $
             cbcharsize=cbcharsize, cbformat=cbformat, cbposition=cbposition,$
             charsize = charsize, $
             noxticks=noxticks, noyticks=noyticks,sample=sample,$
             xrange=xrange, yrange=yrange,_EXTRA=extra

a = a0
x = x0
y = y0
sa = size(a)
nx = sa(1)
ny = sa(2)

; ------------------------------------------
;      max and min data values
; ------------------------------------------

IF (NOT KEYWORD_SET(color)) THEN color = 0

IF (NOT KEYWORD_SET(title)) THEN title = " "

IF (NOT KEYWORD_SET(dmax)) THEN dmax = max(a)

IF (NOT KEYWORD_SET(dmin)) THEN dmin = min(a)

IF (NOT KEYWORD_SET(charsize)) THEN charsize = !P.CHARSIZE

xtickformat = !X.TICKFORMAT
ytickformat = !Y.TICKFORMAT
IF (KEYWORD_SET (noxticks)) THEN BEGIN 
  xtickname   = REPLICATE(' ',10)
  xtickformat = "(A1)"
ENDIF

IF (KEYWORD_SET (noyticks)) THEN BEGIN
  ytickname   = REPLICATE(' ',10)
  ytickformat = "(A1)"
ENDIF

IF (KEYWORD_SET(xrange)) THEN BEGIN

  xs = xrange[0]
  xe = xrange[1]

  is = min(WHERE(x GT xs))
  ie = max(WHERE(x LT xe))

  is = max([0,is])
  ie = min([nx - 1,ie])

  a  = extrac(a0, is, 0, ie - is, ny)
  x  = extrac(x0, is, ie-is)
  sa = size(a)
  nx = sa(1)
ENDIF ELSE BEGIN
  xrange = [min(x),max(x)]
ENDELSE

IF (KEYWORD_SET(yrange)) THEN BEGIN

  ys = yrange[0]
  ye = yrange[1]

  is = min(WHERE(y GT ys))
  ie = max(WHERE(y LT ye))

  is = max([0,is])
  ie = min([ny - 1,ie])

  a  = extrac(a, 0, is, nx, ie-is)
  y  = extrac(y0, is, ie-is)
  sa = size(a)
  ny = sa(2)
ENDIF ELSE BEGIN
  yrange = [min(y),max(y)]
ENDELSE

; ------------------------------------------------
;   set offsets for the figures.
;   All units are in normalized units, i.e. [0,1]
; ------------------------------------------------

imsize_x = pos(2) - pos(0)
imsize_y = pos(3) - pos(1)

; ------------------------------------------
;   make a tv byte-scaled image
; ------------------------------------------

q = reform(a)

xx = x
yy = y

IF (KEYWORD_SET (sample)) THEN BEGIN
  q  = CONGRID(q, sample*nx, sample*ny, cubic=-.5)
  xx = CONGRID(x, sample*nx,/interp)
  yy = CONGRID(x, sample*ny,/interp)
ENDIF

q = BYTSCL(q,max=dmax, min=dmin)
TV, q, pos(0), pos(1), xsize = imsize_x, ysize = imsize_y, /normal
TVLCT,r,g,b,/GET; -- get current colortable
loadct,0  ; -- load black & white
CONTOUR, q, xx, yy, color = color, $
         xstyle = 1, ystyle = 1,$
         xtitle = xtitle, ytitle = ytitle, title=title,$
         /noerase,/nodata, pos=pos, charsize=charsize,$
         xtickname = xtickname, ytickname = ytickname,$
         xtickformat = xtickformat, ytickformat = ytickformat,$
         xrange=xrange,yrange=yrange,_EXTRA=extra
TVLCT,r,g,b

; -------------------------------------------
;   Put colorbar
; -------------------------------------------

IF (NOT KEYWORD_SET(cbcharsize)) THEN cbcharsize = 1.2
IF (NOT KEYWORD_SET(cbwidth))    THEN cbwidth = 0.04
IF (NOT KEYWORD_SET(cbdiv))      THEN cbdiv   = 4
IF (NOT KEYWORD_SET(cbformat))   THEN cbformat = '(F8.2)'

cbdiv = cbdiv - 1

IF (KEYWORD_SET (vbar)) THEN BEGIN

  IF (NOT KEYWORD_SET(cbposition)) THEN cbposition = [pos[2] + 0.01, pos[1]]


  poscb = [cbposition[0], cbposition[1],$
           cbposition[0] + cbwidth*imsize_x,pos[3]]

  COLORBAR2,/vertical, /right,range=[dmin,dmax],division=cbdiv,$
            position = poscb, chars=cbcharsize,color=color,$
            format = cbformat, _EXTRA=extra
ENDIF

IF (KEYWORD_SET (hbar)) THEN BEGIN

  IF (NOT KEYWORD_SET(cbposition)) THEN cbposition = [pos[0], pos[3] + 0.01]

  poscb = [cbposition[0], cbposition[1],$
           pos(2), cbposition[1] + cbwidth*imsize_y]

  COLORBAR2,/horizontal, /bot, range=[dmin,dmax],division=cbdiv,$
            position = poscb, format = cbformat,chars=cbcharsize,color=color,$
            _EXTRA=extra
ENDIF

END
