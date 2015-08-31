;+
;
; NAME:       shockfind
;
; AUTHOR:    Andrea Mignone (mignone@to.astro.it)
;
; PURPOSE:    find shocks in a multidimensional plane
;
; SYNTAX:      array = shockfind(p, vx,vy[,eps_min=eps_min][,eps_max=eps_max])
;
; DESCRIPTION:  shockfind returns a multi-D array (with dimensionality
;               equal to the input arrays) whose elements are
;               numbers between 0 (no shock) and 1 (strong shock).
;
;               A shock is detected when div.V < 0
;               and |grad p|/p > eps. Divergence and gradients
;               are undivided differences.
;
;               The strength of the shock is controlled by
;
;               eps_min <  |grad p|/p < eps_max
;
;
; ARGUMENTS and KEYWORDS:
;
;              p:    an array giving the pressure distribution
;              vx:   an array giving the x-component of velocity
;              vy:   an array giving the y-component of velocity
;              eps_min:  the minimum threshold for a shock to exist
;              eps_max:  the maximum strength of a shock (in units
;                        of |grad p|/p).
;
;
; LAST MODIFIED:   Aug 6 2006 by A.Mignone (mignone@to.astro.it)
;
;-

function shockfind, pr, vx, vy, $
                    eps_min = eps_min,eps_max = eps_max


IF (NOT KEYWORD_SET(eps_min)) THEN eps_min = 0.33
IF (NOT KEYWORD_SET(eps_max)) THEN eps_max = 5.0

sa = size(pr)
nx = sa(1)
ny = sa(2)

sh = fltarr(nx,ny)

; -----------------------------------------------------
;   Compute undived differences in the whole domain
; -----------------------------------------------------

dvx = shift (vx,-1,0) - shift(vx,1,0)
dvy = shift (vy,0,-1) - shift(vy,0,-1)
dpx = shift (pr,-1,0) - shift(pr,1,0)
dpy = shift (pr,0,-1) - shift(pr,0,-1)

; -----------------------------------------------------
;  do boundaries reverting to one-sided first
;  order derivative
; -----------------------------------------------------

for i = 0,nx - 1 do begin
  dvy(i,0) = vy(i,1) - vy(i,0)
  dpy(i,0) = pr(i,1) - pr(i,0)
  dvy(i,ny-1) = vy(i,ny-1) - vy(i,ny-2)
  dpy(i,ny-1) = pr(i,ny-1) - pr(i,ny-2)
endfor

for j = 0,ny - 1 do begin
  dvx(0,j) = vx(1,j) - vy(0,j)
  dpx(0,j) = pr(1,j) - pr(0,j)
  dvx(nx-1,j) = vx(nx-1,j) - vx(nx-2,j)
  dpx(nx-1,j) = pr(nx-1,j) - pr(nx-2,j)
endfor

; ----------------------------------------------
;   Compute |grad p| /p
; ----------------------------------------------

dp = sqrt(dpx*dpx + dpy*dpy)/pr
x  = (dp - eps_min)/(eps_max - eps_min)
x  = x < 1.


sh = ( ((dvx + dvy) LT 0.0) AND (dp GT eps_min))
sh = sh*x

IF (max(sh) eq 0.0) THEN print,"No shock detected" ELSE print,"Shocks detected"
return,sh
end
