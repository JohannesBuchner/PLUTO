;+
;
;  NAME:      REGRID
;
;  AUTHOR:    Andrea Mignone (mignone@ph.unito.it)
;
;  PURPOSE:   Linearly interpolate irregularly-gridded data to a
;             a regular grid. Interpolation is done
;             separately row by row and column by column.
;             On output the original arrays are replaced with the
;             the new ones.
;
;  SYNTAX:    REGRID, img, x1, x2[,n1=integer][,n2=integer]
;                     [, {X | Y}RANGE = [min,max]]
;
;
;
;    img = a 2-D array.
;
;    x1  = a 1-D array with the abscissa of img
;
;    x2  = a 1-D array with the ordinata of img
;
;
;  NOTICE: On output, img, x1 and x2 are replaced by the new
;          regridded arrays. x1 and x2 will have uniform spacing.
;
;
;
;  KEYWORDS:
;
;    n1  = the number of points in the regular grid (x direction)
;
;    n2  = the number of points in the regular grid (y direction)
;
;    XRANGE = a two dimensional vector [xbeg, xend] defining
;             the lower and upper points of the new grid in
;             the x direction
;
;    YRANGE = a two dimensional vector [ybeg, yend] defining
;             the lower and upper points of the new grid in
;             the y direction
;
;
;
;  LAST MODIFIED:  Feb 4, 2006
;
;-

PRO REGRID, img, x, y, n1 = n1, n2 = n2, $
            XRANGE = XRANGE, YRANGE = YRANGE


 simg  = size(img)

 nx = simg(1)
 ny = simg(2)

; ---------------------------------------------
;           set defaults
; ---------------------------------------------

 nx_reg = nx
 ny_reg = ny

 xmax = x(nx - 1)
 xmin = x(0)

 ymax = y(ny - 1)
 ymin = y(0)

; ---------------------------------------------
;          check keywords
; ---------------------------------------------

 IF (KEYWORD_SET (n1)) THEN nx_reg = n1
 IF (KEYWORD_SET (n2)) THEN ny_reg = n2

 IF (KEYWORD_SET (XRANGE)) THEN BEGIN
   xmin = XRANGE[0]
   xmax = XRANGE[1]
 ENDIF

 IF (KEYWORD_SET (YRANGE)) THEN BEGIN
   ymin = YRANGE[0]
   ymax = YRANGE[1]
 ENDIF

; ---------------------------------------------
;     define leftmost and rightmost
;     points for each grid.
; ---------------------------------------------

 dx_reg = (xmax - xmin)/(nx_reg - 1.0)
 dy_reg = (ymax - ymin)/(ny_reg - 1.0)

 x_reg = fltarr(nx_reg)
 y_reg = fltarr(ny_reg)

 FOR i = 0, nx_reg - 1 DO x_reg(i) = xmin + dx_reg*i
 FOR j = 0, ny_reg - 1 DO y_reg(j) = ymin + dy_reg*j

; print," > Regridding image in the x-direction; ",$
;       "n1 = " ,strcompress(string(nx),/remove_all),$
;       " --> ", strcompress(string(nx_reg),/remove_all)


 scrh = fltarr(nx_reg, ny)
 FOR j = 0, ny - 1 DO BEGIN
   scrh(*, j) = INTERPOL(REFORM(img(*,j)), x, x_reg)
 ENDFOR
 img = scrh
 x  = x_reg

; --------------------------------------------
;   interpolate along y
; --------------------------------------------

; print," > Regridding image in the y-direction; ",$
;       "n2 = " ,ARG2STR(ny)," --> ",ARG2STR(ny_reg)

 scrh = fltarr(nx_reg, ny_reg)
 FOR i = 0, nx_reg - 1 DO BEGIN
   scrh(i, *) = INTERPOL(REFORM(img(i,*)), y, y_reg)
 ENDFOR

; img2 = transpose(img)
; FOR i = 0, nx_reg - 1 DO BEGIN
;   scrh(i, *) = INTERPOL(REFORM(img2(*,i)), y, y_reg)
; ENDFOR

 img = scrh
 y   = y_reg
END


