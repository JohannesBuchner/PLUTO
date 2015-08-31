PRO HDF5CUT, nout, _EXTRA=extra, xcut=xcut, ycut=ycut, lrange=lrange, var=var

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

 maxlev=max(lrange)

 start = 1;  to initialize arrays
 IF (KEYWORD_SET(xcut)) THEN BEGIN
   FOR l = maxlev,0,-1 DO BEGIN
     PLOAD, nout, _EXTRA=extra,/hdf5,lev=l,/silent,$
            xrange=[xcut*0.9,xcut*1.1], /NO_VAR_MATCH
     i = min(WHERE(x1 GE xcut))
     IF (l EQ maxlev) THEN AMRBoxesMax = REFORM(AMRBoxes(i,*))
     indx = WHERE( (AMRBoxes(i,*) EQ l) AND REBIN(AMRBoxesMax,n2) EQ l)
     IF (SIZE(indx, /N_ELEM) EQ 1) THEN BEGIN
       print,"Skipping level ",l
       CONTINUE
     ENDIF
     q1 = REFORM(vars(i,  indx, 0, 0:NVAR-1))
     IF (start) THEN BEGIN
       qq    = q1
       xx    = x2(indx)
       dx    = dx2(indx)
       start = 0
     ENDIF ELSE BEGIN
       qq = [qq, q1]
       xx = [xx,x2(indx)]
       dx = [dx,dx2(indx)]
     ENDELSE
   ENDFOR
   x2 = xx
   dx2 = dx
 ENDIF 

 IF (KEYWORD_SET(ycut)) THEN BEGIN
   FOR l = maxlev,0,-1 DO BEGIN
     PLOAD, nout, _EXTRA=extra,/hdf5,lev=l,/silent,$
            yrange=[ycut*0.9, ycut*1.1], /NO_VAR_MATCH
     j = min(WHERE(x2 GE ycut))
     IF (l EQ maxlev) THEN AMRBoxesMax = REFORM(AMRBoxes(*,j))
     indx = WHERE( (AMRBoxes(*,j) EQ l) AND REBIN(AMRBoxesMax,n1) EQ l)
     IF (SIZE(indx, /N_ELEM) EQ 1) THEN BEGIN
       print,"Skipping level ",l
       CONTINUE
     ENDIF
     q1 = REFORM(vars(indx, j, 0, 0:NVAR-1))
     IF (start) THEN BEGIN
       qq    = q1
       xx    = x1(indx)
       dx    = dx1(indx)
       start = 0
     ENDIF ELSE BEGIN
       qq = [qq, q1]
       xx = [xx,x1(indx)]
       dx = [dx,dx1(indx)]
     ENDELSE
   ENDFOR
   x1 = xx
   dx1 = dx
 ENDIF 
 FOR n=0,NVAR-1 DO MATCH_VARNAME, qq(*,n), vname[n]
END
