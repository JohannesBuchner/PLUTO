;+
;
; NAME:      polar
;
; AUTHOR:    Andrea Mignone
;
; PURPOSE:   Interpolates the surface array a from polar coordinates (r, theta) to
;            rectangular coordinates (x, y).
;            On output, a is replaced by the array in cartesian coordinates,
;            while r becomes x, theta becomes y. The array a can be
;            irregularly gridded. This function uses the polar_surface
;            routine of IDL.
;
;            Note: if spherical polar coordinates
;                  are use, then simply let
;                  x2 = !pi/2. - x2
;
; SYNTAX:    polar, a, x1, x2[,sample=sample][,missing=missing]
;
; KEYWORD & PARAMETERS:
;
;        see IDL manual page for polar_surface.
;
;
; DATE:     Aug 20, 2004
;
;-

PRO polar, a, x1, x2, bound=bound, sample=sample,$
           missing = missing,inside = inside

 print," > interpolating from polar to cartesian coordinates..."
 sa    = size(a)
 n1    = sa(1)
 n2    = sa(2)
 r_min = min(x1)

 amin = min(a)

 IF (KEYWORD_SET (missing)) THEN amin = missing
 IF (NOT KEYWORD_SET (inside)) THEN inside = amin

 ;r2 = x1
 ;t2 = x2

 r2  = fltarr(n1,n2)
 t2  = fltarr(n1,n2)
 for i = 0, n1 - 1 do t2(i,*) = x2
 for j = 0, n2 - 1 do r2(*,j) = x1

 xmin = min(x1#cos(x2))
 xmax = max(x1#cos(x2))
 ymin = min(x1#sin(x2))
 ymax = max(x1#sin(x2))

 IF (KEYWORD_SET(sample)) THEN n1 = n1*sample

 dx = (xmax - xmin)/(n1)
 dy = (ymax - ymin)/(n1)

 a = polar_surface(a, r2, t2, spacing=[dx,dy], $
                   missing = amin, bounds=[xmin,ymin, xmax,ymax])

 sa = size(a)

 x1 = fltarr(sa(1))
 x2 = fltarr(sa(2))

 x1(0) = xmin
 x2(0) = ymin

 for i = 1, sa(1)-1 do x1(i) = x1(i - 1) + dx
 for j = 1, sa(2)-1 do x2(j) = x2(j - 1) + dy

 ; -----------------------------------------------
 ;  get rid of that annoying straight
 ;  line that IDL produces when 'polar_surface'
 ;  is called. As a simple correction we
 ;  super-impose a black circle
 ; -----------------------------------------------


 ie = max(WHERE(x1 LT r_min))
 je = max(WHERE(x2 LT r_min))

 FOR i = 0, ie DO BEGIN
 FOR j = 0, je DO BEGIN
   scrh = sqrt(x1(i)*x1(i) + x2(j)*x2(j))
   IF (scrh LT r_min) THEN a(i,j) = inside
 ENDFOR
 ENDFOR



; a  = rotate(reverse(a),3)
; scrh = x2
; x2   = x1
; x1   = scrh
; sa = size(a)


end
