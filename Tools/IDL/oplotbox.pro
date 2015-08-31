;+
;
; NAME:      OPLOTBOX
;
; AUTHOR:    A. Mignone
;
; PURPOSE:   overplot the AMR box layout on the current window.
;            The axis of the plot should be scaled to physical units.
;            (e.g. display should have been called with x1=x1,x2=x2 keywords).
;
;
; SYNTAX:    OPLOTBOX, [, CTAB=integer][, CVAL=vector][, LRANGE=[min,max]][, /POLAR]
;                      [,_EXTRA=extra]
;                      
;
; KEYWORDS
;
;   CTAB = an integer giving the color table to be used when plotting
;          the boxes. Default is the current color table.
;          When CTAB=12 or CTAB=23, color indices are pre-defined to
;          default values so that level 0 is always red, level 1 is green,
;          and so forth.
;
;   CVAL = on output, it returns an array of integers with the 
;             corresponding color indices used to draw the boxes.
;             Useful for legends.
;   LRANGE = a two element integer array of the form [minlev, maxlev] giving
;            the minimum and maximum levels to be plotted.
;            
;   /POLAR  = plot the AMR box layout over a /POLAR display (spherical and polar coordinates) 
;
;   OPLOTBOX also accepts keyword inheritance via the _EXTRA keyword to be passed to plots
;
; EXAMPLES:
;
;  * Example #1: plot boxes on the current graphic device: 
;
;                IDL> OPLOTBOX
;
;  * Example #2: display density maps leaving extra space for a
;                legend plot (rgt=140), overplot the boxes and then 
;                add a legend:
;
;              IDL> display,q,x1=x1,x2=x2,rgt = 140
;              IDL> oplot,cval=vals
;              IDL> legend, ["Level "+strcompress(string(indgen(4)),/remove_all)],$
;                   line=0,colors=vals,pos=[0.83,0.9],/normal,spac=0.1
;
;
; LAST MODIFIED
;
;   May 17th, 2011 by A. Mignone (mignone@ph.unito.it)
;
;-
PRO OPLOTBOX, cval=cval, ctab=ctab, lrange=lrange, polar=polar, $
              islice = islice, jslice=jslice, kslice=kslice, _EXTRA=extra

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

 nlev = (size(AMRLevel))[1]
 IF (nlev EQ 0) THEN BEGIN
   print,"! AMRLevel structure not defined"
   RETURN
 ENDIF

; -----------------------------------
;         check keywords 
; -----------------------------------

 IF (KEYWORD_SET(lrange)) THEN BEGIN
   dim_lrange = size(lrange,/dimensions)
   IF (dim_lrange NE 2) THEN BEGIN
     print,"! lrange must be a 2 element array, lrange=[minlev, maxlev]"
     RETURN
   ENDIF
 ENDIF ELSE BEGIN
   lrange = [0,nlev-1]
 ENDELSE

; ----------------------------------------------------
;  check slice keywords. 
;  Add an offset since the block numbering
;  always starts from beginning of whole domain
;  whereas X|Y|ZRANGE keywords previously given 
;  to HDF5LOAD start counting from the given range
; ----------------------------------------------------

 XYCUT = 1 & XZCUT = 0 & YZCUT = 0
 IF (N_ELEMENTS(islice) GT 0) THEN BEGIN
   islice += AMRLevel[0].ibeg
   XYCUT = 0 & XZCUT = 0 & YZCUT = 1
 ENDIF  
 IF (N_ELEMENTS(jslice) GT 0) THEN BEGIN
   jslice += AMRLevel[0].jbeg
   XYCUT = 0 & XZCUT = 1 & YZCUT = 0
 ENDIF
 IF (N_ELEMENTS(kslice) GT 0) THEN BEGIN
   kslice += AMRLevel[0].kbeg
   XYCUT = 1 & XZCUT = 0 & YZCUT = 0
 ENDIF ELSE BEGIN
   kslice = 0
 ENDELSE

; ------------------------------------
;  save current color table and set
;  color indices
; ------------------------------------

 IF (NOT KEYWORD_SET(CVAL)) THEN CVAL = 255*(1+FINDGEN(nlev))/nlev

 IF (KEYWORD_SET(ctab)) THEN BEGIN
   TVLCT,r,g,b,/GET
   LOADCT,ctab
   CASE ctab OF
              ;red,grn,blu,pr1,cyn,wht
     12:cval = [200, 40,100,140, 90,250, 10,1]
               ;red,grn,ylw,prp,cyn,orn,blu
     23:cval = [255,170,214,  6,130,235,100]
   ELSE:
   ENDCASE
 ENDIF

; ------------------------------------------
;   start plotting 
; ------------------------------------------

 lgnd_str = strarr(nlev)
 FOR nl = lrange[0], lrange[1] DO BEGIN

   lgnd_str[nl] = "Level "+strcompress(string(nl),/remove_all)
   level        = AMRLevel[nl]

  IF (KEYWORD_SET(POLAR) AND ((geometry EQ 'spherical') OR (geometry EQ 'polar'))) THEN BEGIN
  
   IF (geometry EQ 'spherical') THEN BEGIN

    IF (XYCUT) THEN BEGIN  ; default
      FOR nb = 0, level.nbox-1 DO BEGIN
        box = level.box[nb]
        IF ((kslice-box.kb)*(box.ke-kslice) GE 0 ) THEN BEGIN
           xx0 = box.x0*sin(box.y0)
           xx1 = box.x1*sin(box.y0)
           xx2 = box.x1*sin(box.y1)
           xx3 = box.x0*sin(box.y1)

           yy0 = box.x0*cos(box.y0)
           yy1 = box.x1*cos(box.y0)
           yy2 = box.x1*cos(box.y1)
           yy3 = box.x0*cos(box.y1)
        
           plots,[xx0,xx1],[yy0,yy1],color=cval[nl], noclip=0, _EXTRA=extra 
           plots,[xx2,xx3],[yy2,yy3],color=cval[nl], noclip=0, _EXTRA=extra

           theta = findgen(101)/100.*(box.y1-box.y0)+box.y0
           oplot,box.x0*sin(theta),box.x0*cos(theta),color=cval[nl], noclip=0, _EXTRA=extra
           oplot,box.x1*sin(theta),box.x1*cos(theta),color=cval[nl], noclip=0, _EXTRA=extra
 
        ENDIF
      ENDFOR
    ENDIF ELSE BEGIN
     print,'POLAR keyword in spherical curvilinear coordinates only available for R-Theta cuts'
    ENDELSE

   ENDIF
  
   IF (geometry EQ 'polar') THEN BEGIN
    IF (XYCUT) THEN BEGIN  ; default
      FOR nb = 0, level.nbox-1 DO BEGIN
        box = level.box[nb]
        IF ((kslice-box.kb)*(box.ke-kslice) GE 0 ) THEN BEGIN
           xx0 = box.x0*cos(box.y0)
           xx1 = box.x1*cos(box.y0)
           xx2 = box.x1*cos(box.y1)
           xx3 = box.x0*cos(box.y1)

           yy0 = box.x0*sin(box.y0)
           yy1 = box.x1*sin(box.y0)
           yy2 = box.x1*sin(box.y1)
           yy3 = box.x0*sin(box.y1)

           plots,[xx0,xx1],[yy0,yy1],color=cval[nl], noclip=0, _EXTRA=extra
           plots,[xx2,xx3],[yy2,yy3],color=cval[nl], noclip=0, _EXTRA=extra

           theta = findgen(101)/100.*(box.y1-box.y0)+box.y0
           oplot,box.x0*cos(theta),box.x0*sin(theta),color=cval[nl], noclip=0, _EXTRA=extra
           oplot,box.x1*cos(theta),box.x1*sin(theta),color=cval[nl], noclip=0, _EXTRA=extra

        ENDIF
      ENDFOR
    ENDIF ELSE BEGIN
     print,'POLAR keyword in polar curvilinear coordinates only available for R-Phi cuts'
    ENDELSE
   ENDIF
 
  ENDIF ELSE BEGIN

   IF (XYCUT) THEN BEGIN  ; default
     FOR nb = 0, level.nbox-1 DO BEGIN
       box = level.box[nb]
       IF ((kslice-box.kb)*(box.ke-kslice) GE 0 ) THEN BEGIN
         plots, [box.x0, box.x1, box.x1, box.x0, box.x0], $
                [box.y0, box.y0, box.y1, box.y1, box.y0], $
                color = cval[nl], noclip=0,_EXTRA=extra
       ENDIF
     ENDFOR
   ENDIF

   IF (XZCUT) THEN BEGIN  
     FOR nb = 0, level.nbox-1 DO BEGIN
       box = level.box[nb]
       IF ((jslice-box.jb)*(box.je-jslice) GE 0 ) THEN BEGIN
         plots, [box.x0, box.x1, box.x1, box.x0, box.x0], $
                [box.z0, box.z0, box.z1, box.z1, box.z0], $
                color = cval[nl], noclip=0,_EXTRA=extra
       ENDIF
     ENDFOR
   ENDIF

   IF (YZCUT) THEN BEGIN  
     FOR nb = 0, level.nbox-1 DO BEGIN
       box = level.box[nb]
       IF ((islice-box.ib)*(box.ie-islice) GE 0 ) THEN BEGIN
         plots, [box.y0, box.y1, box.y1, box.y0, box.y0], $
                [box.z0, box.z0, box.z1, box.z1, box.z0], $
                color = cval[nl], noclip=0,_EXTRA=extra
       ENDIF
     ENDFOR
   ENDIF
  
  ENDELSE ; POLAR 

 ENDFOR

 ; *******************************************
 ;  restore original color table if necessary
 ; *******************************************

 IF (KEYWORD_SET (ctab)) THEN TVLCT,r,g,b

END
