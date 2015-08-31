pro plotbox,filename,x1 = x1, y1 = y1, x2 = x2, y2 = y2, resf = resf, rot_mat = rot_mat

;  Routine PLOTBOX
;  called by HDF5PLOT to draw the patch binding boxes
;  
;  for parameters explanation see HDF5PLOT
;
;
;  Last updated: 31 May 2008 by Ovidiu Tesileanu



 checkfile,filename

 ifil = H5F_OPEN(filename)
 igen = H5G_OPEN(ifil,'/')

 inlev = H5A_OPEN_NAME(igen,'num_levels')
 nlev  = H5A_READ(inlev)
 H5A_CLOSE,inlev

 lev = strarr(nlev)
 for i=0,nlev-1 do lev[i] = 'level_'+strcompress(string(i),/REMOVE_ALL)

 freb = intarr(nlev)

for nl=nlev-1,0,-1 do begin

 ilev = H5G_OPEN(ifil,lev[nl])

if (nl eq nlev-1) then begin
   freb[nl] = 1
   iprobdom = H5A_OPEN_NAME(ilev,'prob_domain')
   pdom = H5A_READ(iprobdom)
   H5A_CLOSE,iprobdom
   ilo = pdom.lo_i
   jlo = pdom.lo_j
   ihi = pdom.hi_i
   jhi = pdom.hi_j
   nx = ihi-ilo+1
   ny = jhi-jlo+1
   idx = H5A_OPEN_NAME(ilev,'dx')
   dx = H5A_READ(idx)
   H5A_CLOSE,idx
 endif else begin
  irat = H5A_OPEN_NAME(ilev,'ref_ratio')
  rat = H5A_READ(irat)
  H5A_CLOSE,irat
  freb[nl] = freb[nl+1]*rat
 endelse

 H5G_CLOSE,ilev

endfor

IF (NOT KEYWORD_SET (x1)) THEN x1 = 0
IF (NOT KEYWORD_SET (x2)) THEN x2 = nx-1
IF (NOT KEYWORD_SET (y1)) THEN y1 = 0
IF (NOT KEYWORD_SET (y2)) THEN y2 = ny-1
IF (NOT KEYWORD_SET (resf)) THEN resf = 1

tagn = intarr(nx,ny)
tagn[*,*] = 0
color=0
for ll=nlev-1,1,-1 do begin

  color = 255.*ll/(nlev-1.)
; if (ll eq 0) then color = 0
; if (ll eq 1) then color = 90
; if (ll eq 2) then color = 255
 tag = tagn

 ilev = H5G_OPEN(ifil,lev[ll])

 iboxes = H5D_OPEN(ilev,'boxes')
 boxes = H5D_READ(iboxes)
 H5D_CLOSE,iboxes


 nbox = n_elements(boxes.lo_i)

 iref = freb[ll]
 jref = freb[ll]

for ib = 0,nbox-1 do begin

 illoc = boxes[ib].lo_i
 ihloc = boxes[ib].hi_i
 jlloc = boxes[ib].lo_j
 jhloc = boxes[ib].hi_j

 IF ( (illoc*iref GT x2) OR (ihloc*iref LT x1) OR (jlloc*iref GT y2) OR (jhloc*iref LT y1) ) THEN continue


 tagn[illoc*iref:(ihloc+1)*iref-1,jlloc*jref:(jhloc+1)*jref-1] = 1

;; Plot left side of the box...

    IF ( (illoc*iref) LT x1) THEN illoc = FIX(x1/iref)
    IF ( (ihloc*iref) GT x2) THEN ihloc  = FIX(x2/iref)
    IF ( (jlloc*iref) LT y1) THEN jlloc = FIX(y1/iref)
    IF ( (jhloc*iref) GT y2) THEN jhloc  = FIX(y2/iref)


 xloc = (illoc*iref)*dx
 first = 1
 for j = jlloc,jhloc do begin
  for jsub = 0,jref-1 do begin
  if (tag[illoc*iref,j*jref+jsub] eq 0) then begin
   if (first eq 1) then begin
    IF (rot_mat EQ 0) THEN plots,[xloc,xloc],[(j*jref+jsub)*dx,(j*jref+jsub+1)*dx],color=color
    IF (rot_mat EQ 1) THEN plots,[-(j*jref+jsub)*dx,-(j*jref+jsub+1)*dx],[xloc,xloc],color=color
    IF (rot_mat EQ 2) THEN plots,[-xloc,-xloc],[-(j*jref+jsub)*dx,-(j*jref+jsub+1)*dx],color=color
    IF (rot_mat EQ 3) THEN plots,[(j*jref+jsub)*dx,(j*jref+jsub+1)*dx],[-xloc,-xloc],color=color
    IF (rot_mat EQ 4) THEN plots,[(j*jref+jsub)*dx,(j*jref+jsub+1)*dx],[xloc,xloc],color=color
    first = 0
   endif else begin
    IF (rot_mat EQ 0) THEN plots,xloc,(j*jref+jsub+1)*dx,/continue,color=color
    IF (rot_mat EQ 1) THEN plots,-(j*jref+jsub+1)*dx,xloc,/continue,color=color
    IF (rot_mat EQ 2) THEN plots,-xloc,-(j*jref+jsub+1)*dx,/continue,color=color
    IF (rot_mat EQ 3) THEN plots,(j*jref+jsub+1)*dx,-xloc,/continue,color=color
    IF (rot_mat EQ 4) THEN plots,(j*jref+jsub+1)*dx,xloc,/continue,color=color
   endelse
  endif else begin
   first = 1
  endelse
  endfor
 endfor
;; End plot left

;; Plot upper side of the box...
 yloc = ((jhloc+1)*jref)*dx
 first = 1
 for i = illoc,ihloc do begin
  for isub = 0,iref-1 do begin
  if (tag[i*iref+isub,(jhloc+1)*jref-1] eq 0) then begin
   if (first eq 1) then begin
    IF (rot_mat EQ 0) THEN plots,[(i*iref+isub)*dx,(i*iref+isub+1)*dx],[yloc,yloc],color=color
    IF (rot_mat EQ 1) THEN plots,[-yloc,-yloc],[(i*iref+isub)*dx,(i*iref+isub+1)*dx],color=color
    IF (rot_mat EQ 2) THEN plots,[-(i*iref+isub)*dx,-(i*iref+isub+1)*dx],[-yloc,-yloc],color=color
    IF (rot_mat EQ 3) THEN plots,[yloc,yloc],[-(i*iref+isub)*dx,-(i*iref+isub+1)*dx],color=color
    IF (rot_mat EQ 4) THEN plots,[yloc,yloc],[(i*iref+isub)*dx,(i*iref+isub+1)*dx],color=color
    first = 0
   endif else begin
    IF (rot_mat EQ 0) THEN plots,(i*iref+isub+1)*dx,yloc,/continue,color=color
    IF (rot_mat EQ 1) THEN plots,-yloc,(i*iref+isub+1)*dx,/continue,color=color
    IF (rot_mat EQ 2) THEN plots,-(i*iref+isub+1)*dx,-yloc,/continue,color=color
    IF (rot_mat EQ 3) THEN plots,yloc,-(i*iref+isub+1)*dx,/continue,color=color
    IF (rot_mat EQ 4) THEN plots,yloc,(i*iref+isub+1)*dx,/continue,color=color
   endelse
  endif else begin
   first = 1
  endelse
  endfor
 endfor
;; End plot upper

;; Plot right side of the box...
 xloc = ((ihloc+1)*iref)*dx
 first = 1
 for j = jlloc,jhloc do begin
  for jsub = 0,jref-1 do begin
  if (tag[(ihloc+1)*iref-1,j*jref+jsub] eq 0) then begin
   if (first eq 1) then begin
    IF (rot_mat EQ 0) THEN plots,[xloc,xloc],[(j*jref+jsub)*dx,(j*jref+jsub+1)*dx],color=color
    IF (rot_mat EQ 1) THEN plots,[-(j*jref+jsub)*dx,-(j*jref+jsub+1)*dx],[xloc,xloc],color=color
    IF (rot_mat EQ 2) THEN plots,[-xloc,-xloc],[-(j*jref+jsub)*dx,-(j*jref+jsub+1)*dx],color=color
    IF (rot_mat EQ 3) THEN plots,[(j*jref+jsub)*dx,(j*jref+jsub+1)*dx],[-xloc,-xloc],color=color
    IF (rot_mat EQ 4) THEN plots,[(j*jref+jsub)*dx,(j*jref+jsub+1)*dx],[xloc,xloc],color=color
    first = 0
   endif else begin
    IF (rot_mat EQ 0) THEN plots,xloc,(j*jref+jsub+1)*dx,/continue,color=color
    IF (rot_mat EQ 1) THEN plots,-(j*jref+jsub+1)*dx,xloc,/continue,color=color
    IF (rot_mat EQ 2) THEN plots,-xloc,-(j*jref+jsub+1)*dx,/continue,color=color
    IF (rot_mat EQ 3) THEN plots,(j*jref+jsub+1)*dx,-xloc,/continue,color=color
    IF (rot_mat EQ 4) THEN plots,(j*jref+jsub+1)*dx,xloc,/continue,color=color
   endelse
  endif else begin
   first = 1
  endelse
  endfor
 endfor
;; End plot right

;; Plot bottom side of the box...
 yloc = (jlloc*jref)*dx
 first = 1
 for i = illoc,ihloc do begin
  for isub = 0,iref-1 do begin
  if (tag[i*iref+isub,jlloc*jref] eq 0) then begin
   if (first eq 1) then begin
    IF (rot_mat EQ 0) THEN plots,[(i*iref+isub)*dx,(i*iref+isub+1)*dx],[yloc,yloc],color=color
    IF (rot_mat EQ 1) THEN plots,[-yloc,-yloc],[(i*iref+isub)*dx,(i*iref+isub+1)*dx],color=color
    IF (rot_mat EQ 2) THEN plots,[-(i*iref+isub)*dx,-(i*iref+isub+1)*dx],[-yloc,-yloc],color=color
    IF (rot_mat EQ 3) THEN plots,[yloc,yloc],[-(i*iref+isub)*dx,-(i*iref+isub+1)*dx],color=color
    IF (rot_mat EQ 4) THEN plots,[yloc,yloc],[(i*iref+isub)*dx,(i*iref+isub+1)*dx],color=color
    first = 0
   endif else begin
    IF (rot_mat EQ 0) THEN plots,(i*iref+isub+1)*dx,yloc,/continue,color=color
    IF (rot_mat EQ 1) THEN plots,-yloc,(i*iref+isub+1)*dx,/continue,color=color
    IF (rot_mat EQ 2) THEN plots,-(i*iref+isub+1)*dx,-yloc,/continue,color=color
    IF (rot_mat EQ 3) THEN plots,yloc,-(i*iref+isub+1)*dx,/continue,color=color
    IF (rot_mat EQ 4) THEN plots,yloc,(i*iref+isub+1)*dx,/continue,color=color
   endelse
  endif else begin
   first = 1
  endelse
  endfor
 endfor
;; End plot bottom

endfor

H5G_CLOSE,ilev
endfor

 H5G_CLOSE,igen
 H5F_CLOSE,ifil

return
end
