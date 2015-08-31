pro make_frame, img, name=name, number=number,$
               jpg=jpg, jpeg=jpeg,bmp=bmp,$
               max_val = max_val, min_val = min_val, $
               quality = quality,_EXTRA=extra


 IF (NOT KEYWORD_SET(name))    THEN name = 'tmp'
 IF (NOT KEYWORD_SET(number))  THEN number = 0
 IF (NOT KEYWORD_SET(quality)) THEN quality = 100
 IF (NOT KEYWORD_SET(max_val)) THEN max_val = max(img)
 IF (NOT KEYWORD_SET(min_val)) THEN min_val = min(img)


 scrh = string(format='(I4.4)',number)
 scrh = strcompress(scrh,/remove_all)
 name = name+'.'+scrh+'.'

 DONE = 0

;;;;;;;;;;;;;;;;;;;;;;;
;
;   JPEG   OUTPUT
;
;;;;;;;;;;;;;;;;;;;;;;;

 IF (KEYWORD_SET(jpg) OR KEYWORD_SET(jpeg)) THEN BEGIN

   name = name+'jpg'
   ;
   ; Scale the 2D image into the number of colors in the
   ; IDL session.
   ;

   scrh = Bytscl(img,Top=!D.Table_Size-1,max=max_val,min=min_val)

   ; Get the color table vectors

   sz    = Size(scrh,/Dimensions)
   img24 = bytarr(3,sz[0],sz[1])

   ; get the color table vectors

   tvlct,r,g,b,/get

   img24[0,*,*] = r[scrh]
   img24[1,*,*] = g[scrh]
   img24[2,*,*] = b[scrh]

   write_jpeg,name,img24,true=1,quality=quality
   DONE = 1
 endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  BMP   OUTPUT
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF (KEYWORD_SET(bmp)) THEN BEGIN

   name = name+'bmp'
   ;
   ; Scale the 2D image into the number of colors in the
   ; IDL session.
   ;

   scrh = Bytscl(img,Top=!D.Table_Size-1)

   ; Get the color table vectors

   sz    = Size(scrh,/Dimensions)
   img24 = bytarr(3,sz[0],sz[1])

   ; get the color table vectors

   tvlct,r,g,b,/get

   img24[0,*,*] = b[scrh]
   img24[1,*,*] = g[scrh]
   img24[2,*,*] = r[scrh]

   write_bmp,name,img24,r,g,b
   DONE = 1
 endif

;
;  ext NOT AVAILABLE
;

 if (NOT DONE) then print,"Don't know "+ext+$
              " extension"

end
