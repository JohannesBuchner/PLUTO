;+
; NAME:      PARTICLES_OPLOT
;
; AUTHOR:    Andrea Mignone
;
; PURPOSE:  
;
;   Overlay point particles on a pre-existing 2D plot.
;
; SYNTAX:
;
;   PARTICLES_OPLOT[, cvalue_max = cvalue_max][, cvalue_min = cvalue_min]
;                  [, field=field][, symsize=symsize]
;
;
; ARGUMENTS:
;
;   none
;
; KEYWORDS:
;
;   symsize      = the symbol size. By default particles are plotted with
;                  circles.
;
;
;   color        = defines the color of the particle. To color all particles
;		   with same value, put one number. To color each particle 
;		   differently based on some attribute e.g. particle kinetic energy
;                  or part.color, give an array of size number of particles. The 
;                  color is normalised to lie in range 0-255 of the session's color
;                  table. Default 255.
;
;
;   cmax         = the maximum value of the color to be plotted. Values higher are
;		   set to 255.
;
;
;   cmin         = the minimum value of the color to be plotted. Values lower are
;		   set to 0.
;
;
;   x1p          = An array of dimension: nparticles, giving the x1 coordinates
;		   to be plotted
;
;
;   x2p          = An array of dimension: nparticles, giving the x2 coordinates 
;                  to be plotted. 
;		   For 2D plots, default is x1=x1,x2=x2. Need not be specified.
;		   For 3D runs, x1,x2 should be specified if x3 coordinate
;		   need to be plotted.
;		   
;
; EXAMPLES:
;
;  * Example #1: Load fluid and particles data. Overplot them.
;
;    IDL> PLOAD, 10                    ; Load fluid data file (in .dbl format)
;    IDL> DISPLAY,rho, x1=x1, x2=x2    ; Display density
;    IDL> PARTICLES_LOAD, 10           ; Load particle data file (int .dbl format)
;    IDL> PARTICLES_OPLOT, particles   ; Overlay particles on the plot generated
;                                      ; by previous DISPLAY call.
;
;  * Example #2: Overplot just a subsamples of particles
;
;    IDL> Ek    = 0.5*(particle.vx1^2 + particle.vx2^2 + particle.vx3^2)
;    IDL> indx  = WHERE (Ek GT 0.9*MAX(Ek)); Select most energetic particles
;    IDL> part1 = particle(indx) 
;    IDL> PARTICLES_OPLOT, part1        
;
;
;  * Example #3: Overplot for 3D PLUTO data array to make x-z slice
;
;    IDL> PLOAD, 1, var='vx3'                                        
;    IDL> PARTICLES_LOAD,1                                         
;    IDL> DISPLAY, x1=x1, x2=x3, (vx3[*,nx2/2,*]), /vbar   
;
;    ;---Select particles in the y=0 plane. See particles_slice.pro for more.
;    IDL> part1 = PARTICLES_SLICE(particles, xar=x2, dx=dx2, /xz) 
;
;    ;---Plot all particles with color=255 
;    IDL> PARTICLES_OPLOT, part, x1p=part1.x1, x2p=part.x3, color=255
;
;    ;---Plot particles with color based on value of field color
;    IDL> PARTICLES_OPLOT,part, x1p=part.x1, x2p=part.x3, $
;        color=part.color, cvalue_min=0, cvalue_max=4
;
;    ;---Plot particles with color based on kinetic energy
;    IDL>  ek = 0.5*(part.vx1^2.+part.vx2^2.+part.vx3^2.)         
;    IDL>  PARTICLES_OPLOT, part, x1p=part.x1, x2p=part.x3, color=ek 
;
;
; LAST MODIFIED:
;
; June 15, D. Mukherjee (dipanjan.mukherjee@unito.it) 
;
;-
PRO PARTICLES_OPLOT, parts, symsize=symsize, color=color, cmax = cmax, cmin = cmin, $
                     x1p = x1p, x2p = x2p

  nparts = N_ELEMENTS(parts)

; --------------------------------------
;   Set keywords
; --------------------------------------

  IF (NOT KEYWORD_SET (symsize)) THEN symsize = 1.0
  IF (NOT KEYWORD_SET (x1p))     THEN x1p     = parts.x1
  IF (NOT KEYWORD_SET (x2p))     THEN x2p     = parts.x2
  IF (NOT KEYWORD_SET (color))   THEN color   = 255
  
  ncolor = N_ELEMENTS(color)
  ;----ncolor should be either 1, or equal to no. of particles--
  IF ((ncolor ne 1) AND (ncolor ne nparts)) THEN BEGIN
     PRINT,'color must be of size 1 or equal to no. of particles.'
     STOP
  ENDIF

  IF (N_ELEMENTS(x1p) NE N_ELEMENTS(x2p)) THEN BEGIN
     print,"! x1 and x2 arrays must have same dimensions."
     help,x1p
     help,x2p
     STOP
  ENDIF

  IF (N_ELEMENTS (x1p) NE nparts) THEN BEGIN
     print,"! x1 and x2 must have same size as particles"
     STOP
  ENDIF 

; ----------------------------------------------------
;  Set the color value depending on the chosen field
; ----------------------------------------------------

  IF (NOT KEYWORD_SET(cmax)) THEN cmax = MAX(color)
  IF (NOT KEYWORD_SET(cmin)) THEN cmin = MIN(color)

  IF (ncolor gt 1) THEN BEGIN
     cvalue = FLTARR(nparts)
     PRINT,"> Color range: ",ARG2STR(cmin),", ",ARG2STR(cmax)
     FOR i = 0, nparts-1 DO cvalue[i] = color[i]
     cvalue = (cvalue - cmin)/(cmax - cmin)*255.0; Range [0,255]
     cvalue = cvalue < 255;  Saturate large cvalues
     cvalue = cvalue > 0;    Saturate small cvalues
  ENDIF
; ---------------------------------------------------
;   Start plotting
; ---------------------------------------------------
 
;--USERSYM CIRCLE----
  ang = 2*!PI*findgen(49)/48.
  xarr = symsize*cos(ang)  
  yarr = symsize*sin(ang)

  FOR i = 0ULL, nparts-1 DO BEGIN
      IF (ncolor eq 1) THEN BEGIN
         usersym, xarr, yarr, /FILL, color = color 
      ENDIF ELSE  BEGIN
         usersym, xarr, yarr, /FILL, thick = thick,$
                  color = cvalue[i]
      ENDELSE
    xp = x1p[i]
    yp = x2p[i]
    PLOTS, xp, yp, psym = 8
  ENDFOR

END

