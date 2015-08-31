;+
;
;
;  NAME:      get_frame
;
;  AUTHOR:    Andrea Mignone
;
;  PURPOSE:   grab window's content and write
;             an image file. The output file name is required.
;
;  SYNOPSIS:  get_frame,name=name[,/jpg][,/tiff][,/png][,/ppm]
;                                [,/gif][,quality=quality]
;
;
;  KEYWORDS:  describe the image format.
;
;
;  EXAMPLES
;  
;    #1 Read the current window content and dump it to rho.jpg:
;    IDL> get_frame, name = "rho", /jpg
;
;
;  LAST MODIFIED: Jan 13 2006 by A. Mignone (mignone@to.astro.it)
;
;
;-

PRO GET_FRAME, name=fname, jpg=jpg,$
               jpeg=jpeg, gif=gif, tiff=tiff,$
               png=png,ppm=ppm,$
               quality = quality, _EXTRA=extra

 IF (NOT KEYWORD_SET(name))    THEN name = 'tmp'
 IF (NOT KEYWORD_SET(quality)) THEN quality = 100

 DONE = 0

; -------------------
;    gif output
; -------------------

 IF (KEYWORD_SET(gif)) THEN BEGIN
   name = fname+'.gif'
   img  = tvrd(true=1)
   img2 = Color_Quan(img,1,r,g,b)
   write_gif,name,img2,r,g,b,_EXTRA=extra
   DONE = 1
 ENDIF

; ----------------
;  jpeg output
; ----------------

 IF (KEYWORD_SET(jpg) OR KEYWORD_SET(jpeg)) THEN BEGIN
   name = fname+'.jpg'
   img = tvrd(true=3)
   write_jpeg,name,img,quality=quality,true=3
   DONE = 1
 ENDIF

; -------------------
;   tiff output
; -------------------

 IF (KEYWORD_SET(tiff)) THEN BEGIN
   name = fname+'.tiff'
   img  = tvrd(true=1,/order)
   img2 = Color_Quan(img,1,r,g,b)
   write_tiff,name,img2,0,red=r,green=g,blue=b
   DONE = 1

;
;  alternative formulation:
;

;tvlct,r,g,b,/get
;dd = tvrd(true=1,/order)
;write_image,name+'tiff','tiff',dd,r,g,b

 ENDIF

; -------------------
;     png output
; -------------------

 IF (KEYWORD_SET(png)) THEN BEGIN
   name = fname+'.png'
   img  = tvrd(true=1)
   img2 = Color_Quan(img,1,r,g,b)
   write_png,name,img2,r,g,b
   DONE = 1
 ENDIF

; -------------------
;     ppm output
; -------------------

 IF (KEYWORD_SET(ppm)) THEN BEGIN
   name = fname+'.ppm'
   img  = tvrd(true=1)
   img2 = Color_Quan(img,1,r,g,b)
   write_image,name,'ppm',img2,r,g,b
   DONE = 1
 ENDIF


if (NOT DONE) then print,"Image Format Unkown."
END
