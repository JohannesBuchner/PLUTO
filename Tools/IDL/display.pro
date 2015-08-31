;+
;
; NAME:      DISPLAY
;
; AUTHOR:    Andrea Mignone (mignone@ph.unito.it)
;
; PURPOSE:   Make a tvscl image with axis, colorbar, and all the
;            nice stuffs.
;
; SYNTAX:  display, img[,x1=x1][,x2=x2][,BACKGROUND=value]
;                  [,CBDIV=integer][,CHARSIZE=value][,COLOR=value][,/EPS]
;                  [,FILENAME=string][,IMAX=value][,IMIN=value]
;                  [,IMSIZE=value][,LABEL1=string][,LABEL2=string][,/LOG_GRID]
;                  [,/SMOOTH][,TITLE=string][,/HBAR][,/VBAR][,/XSYM][,/YSYM]
;                  [,XRANGE=[min,max]][,YRANGE=[min,max]][,NWIN=integer]
;                  [,/POLAR][,{LFT,RGT,TOP,BOT}=value]
;
;
; KEYWORDS:
;
;
;      x1 =    a one 1-D array with the horizontal coordinates
;              dim(x1) = first dimension of img;
;
;      x2 =    a one 1-D array with the vertical coordinates
;              dim(x2) = second dimension of img; both
;              x1 and x2 must be present one wish to
;              scale axis to coordinates
;
;      BACKGROUND = in integer value in the range [1,255] giving the 
;                   background color.
;
;      CBDIV      = an integer number specifying the number of 
;                   divisions to divide the colorbar into.
;
;      CHARSIZE   = size of characters.
;
;      COLOR      = an integer in the range [1,255] giving contour and 
;                   colorbar should be > 0 and <= 255
;
;      /EPS       = when this keyword is set, it produces eps output.
;                   the image keeps the aspect ratio with a default
;                   horizontal extent of 10cm. Use imsize to shrink/enlarge it.
;
;      FILENAME   = a string giving the name of the file where output is saved;
;                   only works when /eps is used.
;
;      IMAX       = maximum floating value to which img should be scaled to;
;
;      IMIN       = minimum floating value to which img should be scaled to;
;
;      IMSIZE     = a magnification/shrinking scale factor;
;                   for example, imsize=2.0 produces an image twice
;                   as big as the original.
;                   When graphics is NOT eps, it can also be a two
;                   element arrays giving the new image sizes
;
;      LABEL1     = a string label for the x1 - axis
;
;      LABEL2     = a string label for the x2 - axis
;
;      NWIN       = an integer selecting the window.
;
;      /POLAR     = maps a polar (r,phi) domain into cartesian
;                   coordinates. The user must also supply the
;                   x1 (= r) and x2 (= phi) coordinates. It also
;                   recomended that imsize = [nx,ny] be given, where
;                   nx and ny are the final number of points in the
;                   cartesian map.
;
;      /SMOOTH    = smooth the image using cubic interpolation when resizing
;
;      TITLE      = string title of the plot
;
;      /HBAR      = add horizontal color bar at the bottom of the plot
;
;      /VBAR      = add a vertical color bar to the right of the plot
;
;      /XSYM      = symmetrize the image with respect to the x-axis
;
;      /YSYM      = symmetrize the image with respect to the y-axis
;
;      XRANGE     = a two-element vector specifying the abscissa of the
;                    lower and upper boundary for the displayed image.
;
;      YRANGE     = a two-element vector specifying the ordinata of the
;                   lower and upper boundary for the displayed image.
;
;      /XLOG      = set a logarithmic scaling for the x-axis. 
;                   CAUTION: the input grid must be regularly spaced in log scale.
;
;      /YLOG      = set a logarithmic scaling for the y-axis.
;                   CAUTION: the input grid must be regularly spaced in log scale.
;
;     /UNIFORM    = do not attempt to regrid the image and display as it is,
;                   without interpolating
;     LFT,RGT     = these parameters can be used to specify the margin (in pixels)
;     BOT,TOP       between the plot borders and the left (lft), right (rgt),
;                   top (top) and bottom (bot) side of the window.
;
;
;  EXAMPLES:
;
;    * Example #1: Display the logarithm of the intial density and put axis, 
;                  horizontal colorbar and title:
;
;      IDL> display, alog(rho(0)), x1 = x1, x2 = x2, title = "My title", /hbar
;
;    * Example #2: Display another image in a different window with a 
;                  vertical colorbar using the red - blue color table:
;
;      IDL> loadct,33
;      IDL> display, img2, x1 = x1, x2 = x2, /vbar, nwin = 1
;
;    * Example #3: Display a subdomain of the original density map and leave 
;                  extra space to the right of the plot for extra annotations:
;
;      IDL> display, alog(rho(0)), x1 = x1, x2 = x2, $
;                    xrange=[0.3,0.7], yrange=[0.4,0.8],rgt = 120
;
;    * Example #4: Mapping an image from polar to Cartesian on a 400x400 grid:
;
;      IDL> display, img, x1 = x1, x2 = x2, /polar, ims = [400,400]
;
;    * Example #5: Mapping an image from spherical to Cartesian:
;
;      IDL> display, img, x1 = x1, x2 = 0.5*!PI - x2, /polar, ims = [400,400]
;
;
;  LAST MODIFIED:   June 6, 2015 by A.Mignone (mignone@ph.unito.it)
;
;-
PRO DISPLAY, a_in, x1 = x1, x2 = x2, $
           xlog=xlog, ylog=ylog,uniform=uniform,$
           title =title, label1 = label1, label2 = label2,$
           imax = imax, imin = imin, hbar = hbar,$
           charsize=charsize, vbar = vbar, color = color,$
           background = background, eps  = eps,$
           imsize=imsize,filename = filename,xsym=xsym,ysym=ysym,$
           xrange=xrange, yrange=yrange, nwin=nwin, polar=polar,$
           cbdiv = cbdiv, smooth=smooth, noclose=noclose,$
           lft=lft,rgt=rgt,top=top,bot=bot


 old_device = !D.NAME

 ; check input argument

 IF ( (size(a_in))[0] EQ 0) THEN BEGIN
   print,"! DISPLAY: image is not valid"
   return
 ENDIF
 
 ; +++++++++++++++++++++++++++++++++++++++++++++++
 ;        Store original position vector
 ; +++++++++++++++++++++++++++++++++++++++++++++++

 pos0 = !P.POSITION
 
 ;
 ; make copies of argument in case
 ; it's modified
 ;

 a  = reform(a_in)  ; remove dimensions of size 1 (useful for YZ or XZ 3-D slices)
 sa = size(a)

 IF (sa(0) LE 1) THEN BEGIN
   print," ! DISPLAY: image is not 2-D !"
   return
 ENDIF

 IF (NOT KEYWORD_SET(nwin))   THEN nwin = 0

; -------------------------------------
; Check if array dimensions make sense
; Also, if the grid is not uniform,
; set the regrid_flag = 1.
; -------------------------------------

 regrid_flag = 0;
 IF (KEYWORD_SET (x1)) THEN BEGIN
   sx1 = size(x1)
   IF (NOT (sx1(1) EQ sa(1))) THEN BEGIN
     print,' > ERROR: X1 and image 1st dimension are incompatible !'
     return
   ENDIF
   x1_saved  = x1
   dx1 = x1 - shift(x1,1)
   maxdx1 = max(dx1(1:*))
   mindx1 = min(dx1(1:*))
   IF (ABS(maxdx1/mindx1-1.) GT 1.e-3 AND NOT KEYWORD_SET(polar)) THEN BEGIN
     regrid_flag = 1
   ENDIF
 ENDIF ELSE BEGIN
   x1 = 1 + findgen(sa(1))
 ENDELSE

 IF (KEYWORD_SET (x2)) THEN BEGIN
   sx2 = size(x2)
   IF (NOT (sx2(1) EQ sa(2))) THEN BEGIN
     PRINT,' > ERROR: X2 and image 2nd dimension are incompatible !'
     RETURN
   ENDIF
   x2_saved  = x2
   dx2 = x2 - shift(x2,1)
   maxdx2 = max(dx2(1:*))
   mindx2 = min(dx2(1:*))
   IF (ABS(maxdx2/mindx2-1.) GT 1.e-3 AND NOT KEYWORD_SET(polar)) THEN BEGIN
     regrid_flag = 1
   ENDIF
 ENDIF ELSE BEGIN
   x2 = 1 + findgen(sa(2))
 ENDELSE

; ---------------------
;  KEYWORD: xlog, ylog
; ---------------------

 IF (NOT KEYWORD_SET(xlog)) THEN xlog = 0 
 IF (NOT KEYWORD_SET(ylog)) THEN ylog = 0 

 IF (xlog OR ylog) THEN regrid_flag = 0

; --------------------
;  KEYWORD: imax,imin
; --------------------

 IF (NOT KEYWORD_SET(imax))    THEN imax = max(a)
 IF (NOT KEYWORD_SET(imin))    THEN imin = min(a)

 IF (KEYWORD_SET(polar)) THEN BEGIN
   POLAR, a, x1, x2, missing=imin
   sa = size(a)
   sx1 = size(x1)
   sx2 = size(x2)
 ENDIF

; -------------------
;  KEYWORD: xrange
; -------------------

 IF (KEYWORD_SET(xrange)) THEN BEGIN

   xs = xrange[0]
   xe = xrange[1]

   is = min(WHERE(x1 GT xs))
   ie = max(WHERE(x1 LT xe))

 ; limit bounds

   is = max([0,is])
   ie = min([sa(1) - 1,ie])

   a  = extrac(a, is, 0, ie - is, sa(2))
   x1 = extrac(x1, is, ie-is)
   sa = size(a)

 ENDIF

; -------------------
;  KEYWORD: yrange
; -------------------

 IF (KEYWORD_SET(yrange)) THEN BEGIN

   ys = yrange[0]
   ye = yrange[1]

   js = min(WHERE(x2 GT ys))
   je = max(WHERE(x2 LT ye))

 ; limit bounds

   js = max([0,js])
   je = min([sa(2) - 1,je])

   a  = extrac(a, 0, js, sa(1), je - js)
   x2 = extrac(x2, js, je-js)
   sa = size(a)

 ENDIF

; ----------------------
;  Need to regrid ?
; ----------------------

 IF (KEYWORD_SET(uniform)) THEN regrid_flag = 0

 IF (KEYWORD_SET(regrid_flag)) THEN BEGIN
   print,"Regridding Image..."
   REGRID, a, x1, x2
   sx = size(xc)
   sy = size(yc)
   sa = size(a)
 ENDIF

; ----------------------
;  KEYWORDS: xsym, ysym
; ----------------------

 IF (KEYWORD_SET(xsym)) THEN BEGIN
   scrh = mirror(a,bottom=1)
   a    = scrh
   IF (KEYWORD_SET(x2)) THEN BEGIN
     scrh = mirror(x2,left=-1)
     x2   = scrh
   ENDIF
   sa   = size(a)
 ENDIF
 IF (KEYWORD_SET(ysym)) THEN BEGIN
   scrh = mirror(a,left=1)
   a    = scrh
   scrh = mirror(x1,left=-1)
   x1   = scrh
   sa   = size(a)
 ENDIF

; -------------------
;  KEYWORD: charsize
; -------------------

 IF (NOT KEYWORD_SET(charsize)) THEN charsize = MAX([!P.CHARSIZE,1])

; -------------------
;  KEYWORD: imsize
; -------------------

 magn = 1.0
 IF (KEYWORD_SET (imsize)) THEN BEGIN
   scrh = size(imsize)
   IF (scrh[0] EQ 0) THEN BEGIN
     magn   = imsize
     imsize = intarr(2)
     imsize[0] = magn*sa(1)
     imsize[1] = magn*sa(2)
   ENDIF

   IF (NOT KEYWORD_SET (eps)) THEN BEGIN
     IF (KEYWORD_SET(smooth)) THEN a = congrid(a,imsize(0),imsize(1),cubic=-0.5) $
     ELSE a = congrid(a,imsize(0),imsize(1))
     x1 = congrid(x1, imsize(0))
     x2 = congrid(x2, imsize(1))
     sa   = size(a)
   ENDIF ELSE BEGIN

   ENDELSE
 ENDIF
 imsize = [sa(1), sa(2)]

; -------------------
;  KEYWORD: eps
; -------------------

 IF (KEYWORD_SET (eps)) THEN BEGIN        ; ----  Postscript ----

   IF (NOT KEYWORD_SET(filename)) THEN filename = 'idl.eps'

 ENDIF

; -------------------------------------------------------
;  KEYWORDS: color, background, title, label1, label2
; -------------------------------------------------------

 IF (NOT KEYWORD_SET(background))  THEN BEGIN
    background = (!D.NAME EQ 'X') ? 255:0 ; -- use black background on windows --
 ;  IF (KEYWORD_SET(eps)) THEN background = 0
 ENDIF

 IF (NOT KEYWORD_SET(color)) THEN BEGIN
   color = 255-background
 ;  IF (KEYWORD_SET(eps)) THEN color = 255
 ENDIF

 IF (NOT KEYWORD_SET(TITLE))   THEN  TITLE = ' '
 IF (NOT KEYWORD_SET(label1))  THEN label1 = ' '
 IF (NOT KEYWORD_SET(label2))  THEN label2 = ' '

 IF (KEYWORD_SET (cbdiv)) THEN cbdiv = cbdiv-1 ELSE cbdiv = 6

; -----------------------------------------------------------------
; Make coordinate positions for
;
;   * plot area
;   * color bar
;
;  This is an example when /vbar is used:
;
;   +------------------------------------+
;   |              |                     |
;   |              4                     |
;   |              |                     |
;   |   +--------------------+   +---+   |
;   |   |                    |   |   |   |
;   |   |                    |   |   |   |
;   |   |                    |   |   |   |
;   |   |                    |   |   |   |
;   |   |                    |   |   |   |
;   |-1-|                    |-5-|-6-|-2-|
;   |   |                    |   |   |   |
;   |   |                    |   |   |   |
;   |   |                    |   |   |   |
;   |   |                    |   |   |   |
;   |   +--------------------+   +---+   |
;   |              |                     |
;   |              7                     |
;   |              |                     |
;   |   +--------------------+           |
;   |   |          |         |           |
;   |   |          8         |           |
;   |   |          |         |           |
;   |   +--------------------+           |
;   |              |                     |
;   |              3                     |
;   |              |                     |
;   +------------------------------------+
;
;
;
;
; 1 = lft:      horizontal offset between left margin and image
; 2 = rgt:      horizontal offset between image and right margin
; 3 = bot:      vertical offset between bottom margin and image
; 4 = top:      vertical offset between top margin and image
; 5 = bar_rgt_margin
; 6 = bar_rgt_width
; 7 = bar_bot_margin
; 8 = bar_bot_width
;
;
;  When label2 is included xl[0] is augmented by some extra space
;  When label1 is included xl[1] is augmented by some extra space
;
; -----------------------------------------------------------------

 IF (NOT KEYWORD_SET (lft)) THEN lft = 40. + (charsize)*20
 IF (NOT KEYWORD_SET (rgt)) THEN rgt = 24.
 IF (NOT KEYWORD_SET (bot)) THEN bot = 30. + (charsize)*10
 IF (NOT KEYWORD_SET (top)) THEN top = 25.

 bar_rgt_width  = 0.
 bar_bot_width  = 0.
 bar_rgt_margin = 0.
 bar_bot_margin = 0.

 IF (KEYWORD_SET(hbar)) THEN BEGIN
   bar_bot_width  = 20.0
   bar_bot_margin = 25.0 + (charsize)^2*10
   bot = bot - 10
 ENDIF

 IF (KEYWORD_SET(vbar)) THEN BEGIN
   bar_rgt_width  = 20.0
   bar_rgt_margin = 16.0
   rgt = rgt + 18 + (charsize)^2*10  ; allow space for colorbar ticknames
 ENDIF

 IF ( label1 NE ' ') THEN BEGIN
   bar_bot_margin = bar_bot_margin + (charsize)^2*10
 ENDIF

 IF ( label2 NE ' ') THEN BEGIN
   lft = lft + (charsize)^2*15
 ENDIF

 IF ( title NE ' ' ) THEN BEGIN
   top = top + (charsize)^2*10
 ENDIF

 impos = [lft            , bot + bar_bot_width + bar_bot_margin, $
          lft + imsize[0], bot + bar_bot_width + bar_bot_margin + imsize[1]]

 Lx = lft + imsize[0] + bar_rgt_margin + bar_rgt_width + rgt
 Ly = bot + bar_bot_width + bar_bot_margin + imsize[1] + top

; -----------------------------------------------------------------
;
;                            NOW PLOT
;
; -----------------------------------------------------------------

 IF (KEYWORD_SET(eps)) THEN BEGIN

;-------------------------------
;  DEVICE: POSTSCRIPT
;-------------------------------

   SET_PLOT,'ps'

   DEVICE,xsize=10*magn,ysize=10*magn*Ly/Lx, /color, bits_per_pixel=8, $
         /portrait,/encapsulated,filename=filename

   TV, bytscl(a, max=imax, min= imin), impos[0]/Lx, impos[1]/Ly,$
       xsize = imsize[0]/Lx, ysize = imsize[1]/Ly,/normal

 ENDIF ELSE BEGIN

;-------------------------------
;  DEVICE: SCREEN
;-------------------------------

   ; ------------------------
   ;  get current colortable
   ; ------------------------

   TVLCT,r,g,b,/get; -- store current color table in the r,g,b variables --
   LOADCT,0,/SILENT; -- background in black and white --

   old_background = !P.BACKGROUND
   !P.BACKGROUND  = background

   ; -----------------------------------------
   ;  if the current window has the same size
   ;  erase it rather than creating a new one
   ; -----------------------------------------

   IF (!D.WINDOW EQ nwin AND !D.X_SIZE EQ Lx AND !D.Y_SIZE EQ Ly) THEN BEGIN
     ERASE
   ENDIF ELSE BEGIN
     WINDOW, nwin, xsize=Lx, ysize=Ly
   ENDELSE

   TVLCT,r,g,b  ; -- reload color table --
   TV, bytscl(a, max=imax, min=imin), impos(0)+1,impos(1)+1

 ENDELSE

 position = [impos(0)/Lx, impos(1)/Ly, impos(2)/Lx, impos(3)/Ly]
 !P.POSITION=position

; --------------------------------------------------------
 TVLCT,r,g,b,/GET; -- store current color table in the r,g,b variables --
 LOADCT,0,/SILENT; -- background in black and white --

 CONTOUR, a, x1, x2, /nodata,/noerase,xstyle=1,ystyle=1,$
         position = position, xlog=xlog, ylog=ylog,$
         xtitle = label1, ytitle = label2, title = TITLE,$
         color = color,charsize=charsize
 TVLCT,r,g,b

; ------------------------------------------------------
;
;                 PUT   COLORBAR
;
; ------------------------------------------------------

 IF (imin EQ imax) THEN BEGIN
  imin = imin - 0.00001
  imax = imax + 0.00001
 ENDIF

 ; ----------
 ; Horizontal
 ; ----------

 IF (KEYWORD_SET (hbar)) THEN BEGIN

; ---------------------------------------------
; get colorbar position as [xb0, yb0, xb1, yb1]
; y position is set to yoff from the bottom
;
; align x-positions with the plot,
; when bar is horizontal
; ---------------------------------------------

   xb0 = lft/Lx
   yb0 = bot/Ly

   xb1 = xb0 + imsize[0]/Lx
   yb1 = yb0 + bar_bot_width/Ly

   COLORBAR2,position = [xb0, yb0, xb1, yb1],$
            range=[imin,imax], format = '(f6.2)',$
            color = color, divisions = cbdiv, charsize=charsize
 ENDIF


 ; --------
 ; Vertical
 ; --------

 IF (KEYWORD_SET (vbar)) THEN BEGIN

; scale bar with image size

; ---------------------------------------------
; get colorbar position as [xb0, yb0, xb1, yb1]
;
; x position is set to xoff to the right
; of the image
;
; align y-positions with the plot,
; when bar is vertical
; ---------------------------------------------

   xb0 = (lft + imsize[0] + bar_rgt_margin)/Lx
   yb0 =  impos(1)/Ly

   xb1 = xb0 + bar_rgt_width/Lx
   yb1 = yb0 + imsize[1]/Ly
   COLORBAR2,/vertical,position=[xb0,yb0,xb1,yb1],$
            range=[imin,imax], format = '(f6.2)',$
            /right, color = color, divisions = cbdiv,charsize=charsize
 ENDIF

;
; Restore previous device if eps graphics was produced
;

 IF (KEYWORD_SET(eps)) THEN BEGIN
   IF (NOT KEYWORD_SET(noclose)) THEN BEGIN
     DEVICE,/close
     SET_PLOT,old_device
   ENDIF
 ENDIF

 IF (KEYWORD_SET(x1_saved)) THEN x1 = x1_saved
 IF (KEYWORD_SET(x2_saved)) THEN x2 = x2_saved

 IF (NOT KEYWORD_SET (eps)) THEN !P.BACKGROUND = old_background

; --------------------------------------------
;   restore original position vector
; --------------------------------------------

; !P.POSITION = pos0

END



