;+
; NAME:   MIRROR
; 
; AUTHOR:   A. Mignone
;
; PURPOSE:  Build a symmetric image out of a bidimensional array a.
;
; SYNTAX:  MIRROR, a[, LEFT=left][, RIGHT=right][, TOP=top][, BOTTOM=bottom]
;
; ARGUMENTS:  
;
;   a   = 2D array
;
; KEYWORDS
;
;  One keyword between LEFT, ..., BOTTOM should be given with values equal to
;
;  1 = reflect the image w/ respect to the axis
;     [useful w/ reflective boundary conditions]
;
;  2 = reflect the image w/ respect to the axis midpoint
;      [useful w/ point-reflective boundary conditions]
;
;  3 = do not reflect, simply copy the array
;      [useful w/ periodic boundary conditions]
;
; LAST MODIFIED:   Aug 21, 2012 by A. Mignone (mignone@ph.unito.it)
;
;-

FUNCTION mirror, a, left=left, right=right, top=top, bottom=bottom

; -- determine array size --

ss   = size(a)
ndim = ss[0]
nx   = ss[1]
ny   = ss[2]

IF (ndim EQ 1) THEN ny = 1
tmp  = a

; -- determine array type (4 = float, 5 = double, etc...) --

type = SIZE(a, /type)

; -----------------------------
;        bottom 
; -----------------------------

IF (KEYWORD_SET(bottom)) THEN BEGIN
  ng = bottom/abs(bottom)
  bottom = abs(bottom)

  scrh = MAKE_ARRAY(nx, 2*ny, type=type)

  IF (bottom eq 1) THEN BEGIN

    scrh(*,ny:*)   =  tmp(*,*)
    scrh(*,0:ny-1) = ng*reverse(tmp(*,*),2)
    
  endif
  if (bottom eq 2) then begin

    scrh(*,ny:*)   =  tmp(*,*)
    scrh(*,0:ny-1) = ng*rotate(tmp(*,*),2)
  endif
  if (bottom eq 3) then begin
     scrh(*,ny:*)   =  tmp(*,*)
     scrh(*,0:ny-1) =  ng*tmp(*,*)
  endif


  nx  = nx
  ny  = 2*ny
  tmp = scrh
ENDIF

; ----------------
;  Symmetrize top
; ----------------

if (KEYWORD_SET(top)) then begin
  ng = top/abs(top)
  top = abs(top)
  scrh = MAKE_ARRAY(nx, 2*ny, type=type)

  if (top eq 1) then begin
    scrh(*,ny:*)   = ng*reverse(tmp(*,*),2)
    scrh(*,0:ny-1) = tmp(*,*)
  endif
  if (top eq 2) then begin
    scrh(*,ny:*)   = ng*rotate(tmp(*,*),2)
    scrh(*,0:ny-1) = tmp(*,*)
  endif
  if (top eq 3) then begin
    scrh(*,ny:*)   = ng*tmp(*,*)
    scrh(*,0:ny-1) = tmp(*,*)
  endif

  nx  = nx
  ny  = 2*ny
  tmp = scrh
endif

;
; Symmetrize left
;

if (KEYWORD_SET(left)) then begin
  ng = left/abs(left)
  left = abs(left)
;  scrh = fltarr(2*nx,ny)
  scrh = MAKE_ARRAY(2*nx, ny, type=type)
  if (left eq 1) then begin
    scrh(nx:*,*)   = tmp(*,*)
    scrh(0:nx-1,*) = ng*reverse(tmp(*,*),1)
  endif
  if (left eq 2) then begin
    scrh(nx:*,*)   = tmp(*,*)
    scrh(0:nx-1,*) = ng*rotate(tmp(*,*),2)
  endif
  if (left eq 3) then begin
    scrh(nx:*,*)   = tmp(*,*)
    scrh(0:nx-1,*) = ng*tmp(*,*)
  endif

  nx = 2*nx
  ny = ny
  tmp = scrh
endif

;
; Symmetrize right
;

if (KEYWORD_SET(right)) then begin
  ng = right/abs(right)
  right = abs(right)
;  scrh = fltarr(2*nx,ny)
  scrh = MAKE_ARRAY(2*nx,ny, type=type)
  if (right eq 1) then begin
    scrh(nx:*,*)   = ng*reverse(tmp(*,*),1)
    scrh(0:nx-1,*) = tmp(*,*)
  endif
  if (right eq 2) then begin
    scrh(nx:*,*)   = ng*rotate(tmp(*,*),2)
    scrh(0:nx-1,*) = tmp(*,*)
  endif
  if (right eq 3) then begin
    scrh(nx:*,*)   = ng*tmp(*,*)
    scrh(0:nx-1,*) = tmp(*,*)
  endif

endif

return,scrh
end
