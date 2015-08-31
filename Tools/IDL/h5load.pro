;+
;
; NAME:      H5LOAD
;
; AUTHOR:    C. Zanni
;
; REQUIRES:  HDF5 support.
;
; PURPOSE:   Read PLUTO variables stored in a (pixie) HDF5 data file 
;            (static grid).  It is automatically called by PLOAD when the 
;            /H5 keyword is supplied.
;
; SYNTAX:    H5LOAD, filename
;
; ARGUMENTS
;
;   filename   = the data file name
;
; KEYWORDS 
;
;   /SILENT    Set this keyword to suppress output
;
; LAST MODIFIED
;
;   Sept 25, 2012 by C. Zanni (zanni@oato.inaf.it)
;
;-
PRO H5LOAD, filename, silent=silent

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

 checkfile,filename
 ifile = H5F_OPEN(filename)
 
 foundtstep = 0
 foundstvar = 0

 numgrp = H5G_GET_NMEMBERS(ifile,"/")
 FOR ng=0,numgrp-1 DO BEGIN
  strg = H5G_GET_MEMBER_NAME(ifile,"/",ng)
  IF (strcmp(strg,"Timestep",8)) THEN BEGIN
   foundtstep = 1
   strstep = strg 
  ENDIF
 ENDFOR
 
 IF (foundtstep eq 1) THEN BEGIN
  istep = H5G_OPEN(ifile,strstep)
  numgrp = H5G_GET_NMEMBERS(ifile,strstep)
  IF (numgrp gt 1) THEN foundstvar = 1
 ENDIF ELSE BEGIN
  istep = ifile
  IF (numgrp gt 3) THEN foundstvar = 1
 ENDELSE

  igrp  = H5G_OPEN(istep,"vars")
  numvars = H5G_GET_NMEMBERS(istep,"vars")

  FOR nv=0,numvars-1 DO BEGIN
   vname = H5G_GET_MEMBER_NAME(istep,"vars",nv)
   idata = H5D_OPEN(igrp,vname)
   vpt = H5D_READ(idata)
   MATCH_VARNAME, vpt, vname, silent=silent
   H5D_CLOSE,idata
  ENDFOR

  H5G_CLOSE,igrp

  IF (foundstvar eq 1) THEN BEGIN
   igrp  = H5G_OPEN(istep,"stag_vars")
   numvars = H5G_GET_NMEMBERS(istep,"stag_vars")

   FOR nv=0,numvars-1 DO BEGIN
     vname = H5G_GET_MEMBER_NAME(istep,"stag_vars",nv)
     idata = H5D_OPEN(igrp,vname)
     vpt = H5D_READ(idata)
     MATCH_VARNAME, vpt, vname, silent=silent
     H5D_CLOSE,idata
   ENDFOR
   H5G_CLOSE,igrp
  ENDIF

  IF (foundtstep eq 1) THEN H5G_CLOSE,istep
  H5F_CLOSE,ifile

 RETURN
END
