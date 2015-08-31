;+
;
;  NAME:      set_multi_plot_pos
;
;  AUTHOR:    Andrea Mignone (mignone@to.astro.it)
;
;  PURPOSE:   Compute position coordinates for multi plot figures
;
;  SYNTAX:    pos = set_multi_plot_pos(ncol,nrow),
;                   [,top=top][,bot=bot][,lft=lft][,rgt=rgt]
;                   [,xsep=xsep][,ysep=ysep]
;
;
;             where nrow and ncol are the number of rows and columns
;             in the plot. Note: this function only computes
;             the coordinates for the position vector; it does not
;             produce any plot!
;             On output, pos(i,j,*) will be a 1-D array giving
;             the position of the plot as (x,y,x+xsep,y+ysep).
;             For example, in a 2 x 2 multiplot figure you have
;
;             pos(0,0)    pos(1,0)
;
;             pos(0,1)    pos(1,1)
;
;  KEYWORDS:
;
;
;       top = a scalar specifying the space between the top margin and
;             the first row
;
;       bot = a scalar specifying the space between the bottom margin and
;             the last row
;
;       lft = a scalar specifying the space between the left margin and
;             the first column
;
;       rgt = a scalar specifying the space between the right margin and
;             the last column
;
;       xsep = horizontal separation between plots
;
;       ysep = vertical separation between plots
;
;
;  LAST MODIFIED:  March 12, 2007
;-

function set_multi_plot_pos, ncol, nrow,$
         bot = bot, top=top, lft = lft, rgt = rgt, $
         xsep = xsep, ysep = ysep

IF (NOT KEYWORD_SET (top)) THEN top = 0.1
IF (NOT KEYWORD_SET (bot)) THEN bot = 0.1
IF (NOT KEYWORD_SET (lft)) THEN lft = 0.1
IF (NOT KEYWORD_SET (rgt)) THEN rgt = 0.1
IF (NOT KEYWORD_SET (xsep)) THEN xsep = 0.1
IF (NOT KEYWORD_SET (ysep)) THEN ysep = 0.1

xl   = lft
xr   = rgt
yl   = bot
yr   = top

xsz = (1.0 - xl - xr - (ncol-1.0)*xsep)/(ncol*1.0)
ysz = (1.0 - yl - yr - (nrow-1.0)*ysep)/(nrow*1.0)

pos = fltarr(ncol, nrow, 4)

for j = 0, nrow - 1 do begin
for i = 0, ncol - 1 do begin

  x = xl + i*(xsz + xsep)
  y = yl + j*(ysz + ysep)

  ; ----------------------------------------
  ;   pos(i,j) starts from top left corner
  ; ----------------------------------------

  pos(i, nrow - j - 1,*) = [x, y, x+xsz, y+ysz]

endfor
endfor

return,pos
end
