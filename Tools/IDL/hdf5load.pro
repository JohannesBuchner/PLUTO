;+
;
; NAME:      HDF5LOAD
;
; AUTHOR:    A. Mignone, O. Tesileanu, C. Zanni
;
; REQUIRES:  HDF5 support.
;
; PURPOSE:   Read a HDF5 Chombo data-file and store its content on the 
;            variables vars. It is automatically called by PLOAD when the 
;            /HDF5 keyword is supplied.
;
; SYNTAX:    HDF5LOAD, nout[, VAR=string][, FILENAME=string][, /INTERP]
;                          [,LEVEL=value][,{X | Y | Z}RANGE=[min,max]]
;                          [, /NO_VAR_MATCH][,/NODATA]
;
; KEYWORDS
;
;   VAR      = the variable name, e.g. "rho" or "bx2"
;   FILENAME = the file name. If omitted, "data.nnnn.hdf5" is assumed.
;   LEVEL    = the output level resolution (default = 0)
;   INTERP   = set this keyword if you wish to interpolate a coarse level
;              down to the resolution of the selected level using bilinear
;              interpolation provided by REBIN. Otherwise nearest REBIN uses
;              nearest neighbor sampling for both magnification and minification.
;   X1RANGE = a 2 element vector giving the x-range of the portion
;             to be extracted.
;   X2RANGE = a 2 element vector giving the y-range of the portion
;             to be extracted.
;   X3RANGE = a 2 element vector giving the z-range of the portion
;             to be extracted.
;   NO_VAR_MATCH = do not try to match variable names with those provides 
;                  by PLUTO and copy the corresponding storage areas
;                  (e.g. rho, vx1, vx2, etc...). In other words, all variables 
;                  read from the HDF file are stores inside vars(i,j,k,nv). 
;                  Useful when HDF5LOAD is called by another function (e.g.
;                  HDF5CUT).
;
; EXAMPLES:
;
;  * Example #1: read "data.0002.hdf5" and set the output file resolution
;                equivalent to the 4th level of refinement:
;
;                IDL> HDF5LOAD, 2, lev=4
;
;
;  * Example #2: read "data.0005.hdf5" and set the output file resolution
;                equivalent to the 2th level of refinement:
;
;                IDL> HDF5LOAD, 5, lev=2, var="v1"
;
;
; LAST MODIFIED
;
;   Oct 19, 2012 by A. Mignone (mignone@ph.unito.it) & C. Zanni (zanni@oato.inaf.it)
;-
PRO HDF5LOAD, nout, dir, FILENAME=FILENAME, var = var, INTERP=INTERP,$
              NO_VAR_MATCH=NO_VAR_MATCH, NODATA=NODATA,$
              level = level, silent=silent,$
              x1range=x1range, x2range=x2range, x3range=x3range

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

; *********************************************
;      check if we have a valid file
; *********************************************

 IF (NOT KEYWORD_SET(filename)) THEN BEGIN
   filename = DIR+"data."+string(nout,format='(I4.4)')+".hdf5"
 ENDIF

 checkfile,filename

 ifil = H5F_OPEN(filename)

; *********************************************
;  the number of spatial dimensions can
;  directly be read from the file
; *********************************************

 igen = H5G_OPEN(ifil,'/Chombo_global')
 dim  = H5A_READ(H5A_OPEN_NAME(igen,'SpaceDim'))

; *********************************************
;  get number of levels (NLEV) and set the
;  resolution level at which performing output
; *********************************************

 igen = H5G_OPEN(ifil,'/')

 inlev = H5A_OPEN_NAME(igen,'num_levels')
 NLEV = H5A_READ(inlev)
 H5A_CLOSE,inlev

 IF (NOT KEYWORD_SET(level)) THEN level = 0
 IF (level GE NLEV) THEN BEGIN
   print,"! Maximum level is ",NLEV-1
   RETURN
 ENDIF

; *********************************************
;    obtain current simulation time (t)
; *********************************************

 intime = H5A_OPEN_NAME(igen,'time')
 t = H5A_READ(intime)
 H5A_CLOSE,intime

; ******************************************
;   get number of variables (NVAR)
;   and their names         (VNAME[nv])
; ******************************************

 icomp = H5A_OPEN_NAME(igen,'num_components')
 NVAR  = H5A_READ (icomp)
 H5A_CLOSE,icomp

 scomp = STRARR(NVAR)
 VNAME = STRARR(NVAR)
 FOR nv=0,NVAR - 1 DO BEGIN
   scomp[nv] = 'component_'+strcompress(string(nv),/REMOVE_ALL)
   icomp = H5A_OPEN_NAME(igen, scomp[nv])
   VNAME[nv] = H5A_READ(icomp)
   H5A_CLOSE,icomp
 ENDFOR

; ******************************************
;     select variables to load
; ******************************************

 varbeg = 0
 varend = NVAR-1
 IF (KEYWORD_SET(var)) THEN BEGIN
   FOR nv = 0, NVAR-1 DO BEGIN
     IF (var EQ VNAME(nv)) THEN BEGIN
       varbeg = nv
       varend = nv
     ENDIF
   ENDFOR
 ENDIF

; ******************************************************
;  find the refinement ratio between levels (freb[nl]),
;  and the output grid corresponding to the given level.
;  freb[level] = 1 always.
; ******************************************************

 lev = strarr(NLEV)
 for i=0,NLEV-1 do lev[i] = 'level_'+strcompress(string(i),/REMOVE_ALL)
 freb = intarr(NLEV)

 FOR nl = level, 0, -1 DO BEGIN

   ilev = H5G_OPEN(ifil,lev[nl])

   IF (nl EQ level) THEN BEGIN
     iprobdom = H5A_OPEN_NAME(ilev,'prob_domain')
     pdom = H5A_READ(iprobdom)
     H5A_CLOSE,iprobdom

     jbeg = 0 & jend = 0 & ny = 1
     kbeg = 0 & kend = 0 & nz = 1
     SWITCH dim OF
       3: BEGIN & kbeg = pdom.lo_k & kend = pdom.hi_k & nz = kend - kbeg + 1 & END
       2: BEGIN & jbeg = pdom.lo_j & jend = pdom.hi_j & ny = jend - jbeg + 1 & END
       1: BEGIN & ibeg = pdom.lo_i & iend = pdom.hi_i & nx = iend - ibeg + 1 & END
     ENDSWITCH

   ; **** get grid spacing ****

     idx = H5A_OPEN_NAME(ilev,'dx')
     dx  = H5A_READ(idx)
     H5A_CLOSE,idx

     x1b = 0.0
     x2b = 0.0
     x3b = 0.0

     freb[nl] = 1

     ystr = 1.
     zstr = 1.

     igeom = 1
     geometry = 'cartesian'
     str_geom = "> Some information about the geometry is missing ... "
     logr = 0

   ; **** get coordinate origin ****

     ; **** catch errors when trying to read old hdf5 files **** 

     CATCH, err_domBeg
     IF (err_domBeg NE 0) THEN BEGIN
       print, "! HDF5LOAD: coordinate origins cannot be found... skipping"
       CATCH,/cancel
       CONTINUE
     ENDIF

     ixb = H5A_OPEN_NAME(ilev,'domBeg1')
     x1b = H5A_READ(ixb)
     H5A_CLOSE,ixb

     IF (dim GE 2) THEN BEGIN
       ixb = H5A_OPEN_NAME(ilev,'domBeg2')
       x2b = H5A_READ(ixb)
       H5A_CLOSE,ixb
     ENDIF

     IF (dim EQ 3) THEN BEGIN
       ixb = H5A_OPEN_NAME(ilev,'domBeg3')
       x3b = H5A_READ(ixb)
       H5A_CLOSE,ixb
     ENDIF

   ; **** get stretch factor ****

     ; **** catch errors when trying to read old hdf5 files ****

     CATCH, err_stretch
     IF (err_stretch NE 0) THEN BEGIN
       print, "! HDF5LOAD: stretch factors cannot be found... skipping"
       CATCH,/cancel
       CONTINUE
     ENDIF

     IF (dim GE 2) THEN BEGIN
       istr = H5A_OPEN_NAME(ilev,'g_x2stretch')
       ystr = H5A_READ(istr)
       H5A_CLOSE,istr
     ENDIF
 
     IF (dim EQ 3) THEN BEGIN
       istr = H5A_OPEN_NAME(ilev,'g_x3stretch')
       zstr = H5A_READ(istr)
       H5A_CLOSE,istr
     ENDIF

   ; **** get geometry ****

     ; **** catch errors when trying to read old hdf5 files ****

     CATCH, err_geom
     IF (err_geom NE 0) THEN BEGIN
       print, "! HDF5LOAD: geometry can not be found... skipping"
       CATCH,/cancel
       CONTINUE
     ENDIF

     ig    = H5A_OPEN_NAME(ilev,'geometry')
     igeom = H5A_READ(ig)
     H5A_CLOSE,ig

     CASE igeom OF
      1: BEGIN geometry = 'cartesian'   & str_geom = "> Geometry:    Cartesian"   & END
      2: BEGIN geometry = 'cylindrical' & str_geom = "> Geometry:    Cylindrical" & END
      3: BEGIN geometry = 'polar'       & str_geom = "> Geometry:    Polar"       & END
      4: BEGIN geometry = 'spherical'   & str_geom = "> Geometry:    Spherical"   & END
     ELSE:
     ENDCASE

   ; **** get logarithmic factor ****

     ; **** catch errors when trying to read old hdf5 files ****

     CATCH, err_logr
     IF (err_logr NE 0) THEN BEGIN
       print, "! HDF5LOAD: logr can not be found... skipping"
       CATCH,/cancel
       CONTINUE
     ENDIF

     ilog = H5A_OPEN_NAME(ilev,'logr')
     logr = H5A_READ(ilog)
     H5A_CLOSE,ilog

     IF (igeom GE 3) THEN BEGIN
      CASE logr OF
      0: str_geom += " - linear radius"
      1: str_geom += " - logarithmic radius"
      ELSE:
      ENDCASE
     ENDIF
  
   ENDIF ELSE BEGIN
     irat = H5A_OPEN_NAME(ilev,'ref_ratio')
     rat  = H5A_READ(irat)
     H5A_CLOSE,irat
     freb[nl] = freb[nl+1]*rat
   ENDELSE
   H5G_CLOSE,ilev

    
 ENDFOR

 dx0 = dx*freb[0]; ** dx of the base (level 0) grid **

; *******************************************************
;  when either X1RANGE, X2RANGE or X3RANGE keyword are
;  given, ibeg and iend at the given level should align
;  with base grid to make overlapping between blocks at
;  different resolution always possible:
;
;
;                 ib0     ie0
;     |-------|-------|-------|-------|
;             |-|-|-|-|-|-|-|-|
;              b             e
;
;   b=ib0*freb[0], e = (ie0+1)*freb[0]-1
;
; *****************************************************

 IF (KEYWORD_SET (x1range)) THEN BEGIN
   IF (logr EQ 0) THEN BEGIN
    x1range = x1range-x1b
   ENDIF ELSE BEGIN
    xbeg = alog(x1range[0]/x1b)
    xend = alog(x1range[1]/x1b)
    x1range=[xbeg,xend] 
   ENDELSE
   ibeg0 = LONG(min(x1range)/dx0)     & iend0 = LONG(max(x1range)/dx0)
   ibeg  = MAX([ibeg, ibeg0*freb[0]]) & iend  = MIN([iend, (iend0+1)*freb[0]-1])
   nx    = LONG(iend-ibeg+1L)
 ENDIF

 IF (KEYWORD_SET (x2range)) THEN BEGIN
   x2range = (x2range - x2b)/ystr
   jbeg0 = LONG(min(x2range)/dx0)     & jend0 = LONG(max(x2range)/dx0)
   jbeg  = MAX([jbeg, jbeg0*freb[0]]) & jend  = MIN([jend, (jend0+1)*freb[0]-1])
   ny    = LONG(jend-jbeg+1L)
 ENDIF

 IF (KEYWORD_SET (x3range) AND dim EQ 3) THEN BEGIN
   x3range = (x3range - x3b)/zstr
   kbeg0 = LONG(min(x3range)/dx0)     & kend0 = LONG(max(x3range)/dx0)
   kbeg  = MAX([kbeg, kbeg0*freb[0]]) & kend  = MIN([kend, (kend0+1)*freb[0]-1])
   nz    = LONG(kend-kbeg+1L)
 ENDIF

; ***********************************************************************
;  allocate memory for:
;
;   - variables,
;   - AMRLevel structure, containing the collection
;     of all boxes from all levels rescaled at the
;     requested resolution.
;
;     AMRLevel is an array of structures, each containing
;     a collection of rectangular boxes with lower and
;     upper corner coordinates given by
;
;     [(box.x0, box.y0, box.z0), (box.x1, box.y1, box.z1)]
;
;     and indices
;
;     [(box.ib, box.jb, box.kb), (box.ie, box.je, box.ke)]
;
;     Each AMRLevel[nl] also contains the beginning and ending indices
;     (e.g. AMRLevel[0].ibeg ... AMRLevel[0].kend) of the portion of
;     the domain being loaded through X1|X2|X3RANGE keywords.
;
;     IMPORTANT:
;       - the structure will contain information up to
;         the requested level.
;       - the indices box.ib, etc.. are rescaled down
;         to the level equivalent resolution.
;
;     AMRBoxes is a single byte array giving the level
;     integer number as function of the coordinate indices.
;
; ***********************************************************************

 if (NOT KEYWORD_SET(NODATA)) THEN vars = fltarr(nx,ny,nz,varend-varbeg+1)
 NMAX_BOXES = 16384
 AMRbox = {x0:0.0, x1: 0.0, ib: 0L, ie: 0L,$
           y0:0.0, y1: 0.0, jb: 0L, je: 0L,$
           z0:0.0, z1: 0.0, kb: 0L, ke: 0L}
 AMRLevel = REPLICATE({nbox:0, $
                       ibeg:ibeg, iend:iend, jbeg:jbeg, jend:jend, kbeg:kbeg, kend:kend,$
                       box: REPLICATE(AMRbox,NMAX_BOXES)},level+1)
 AMRBoxes = BYTARR(nx,ny,nz)

; ******************************************
;        do some printing here...
; ******************************************

 IF (NOT KEYWORD_SET(SILENT)) THEN BEGIN
   print,"> Filename:", filename,$
         "  [dim = ", ARG2STR(dim),"]",$
         "  [nvar = ",ARG2STR(NVAR),"]",$
         "  [max ref level = ", ARG2STR(NLEV-1),"]"

   domstring = " "
   IF (logr EQ 0) THEN BEGIN
    xbeg = x1b+ibeg*dx
    xend = x1b+(iend+1)*dx
   ENDIF ELSE BEGIN 
    xbeg = x1b*exp(ibeg*dx)
    xend = x1b*exp((iend+1)*dx)
   ENDELSE

   SWITCH dim OF
     3: domstring = " x [ "+ARG2STR(x3b+kbeg*dx*zstr)+", "+ARG2STR(x3b+(kend+1)*dx*zstr)+"]"
     2: domstring = " x [ "+ARG2STR(x2b+jbeg*dx*ystr)+", "+ARG2STR(x2b+(jend+1)*dx*ystr)+"]"+domstring
     1: domstring = "   [ "+ARG2STR(xbeg)+", "+ARG2STR(xend)+"]"+domstring
   ENDSWITCH

   resstring = " "
   SWITCH dim OF
     3: resstring = " x "+ARG2STR(nz)
     2: resstring = " x "+ARG2STR(ny)+resstring
     1: resstring =       ARG2STR(nx)+resstring
   ENDSWITCH

   print,str_geom
   print,"> Domain: ",domstring
   print,"> output size: ",resstring,"  [required level: ",ARG2STR(level),"]"
  ENDIF

; ******************************************
;               read data
; ******************************************

 FOR nl = 0, level DO BEGIN;  ** Loop on Levels **
   ilev  = H5G_OPEN(ifil,lev[nl])
   IF (NOT KEYWORD_SET(NODATA)) THEN BEGIN
     idata = H5D_OPEN(ilev,'data:datatype=0')
     data  = H5D_READ(idata)
     H5D_CLOSE,idata 
   ENDIF

   iboxes = H5D_OPEN(ilev,'boxes')
   boxes  = H5D_READ(iboxes)

   H5D_CLOSE,iboxes

   nbox   = n_elements(boxes.lo_i)
   ncount = long(0)
   AMRLevel[nl].nbox = nbox

   jb = 0 & je = 0 & nby = 1
   kb = 0 & ke = 0 & nbz = 1
   IF (nbox GT NMAX_BOXES) THEN BEGIN
     print, "! number of boxes exceeds max default size",nbox
     return
   ENDIF
   FOR nb = 0L,nbox-1 DO BEGIN;   ** Loop on all boxes of a given level **

     ; ** box indexes **

     SWITCH dim OF
       3: BEGIN & kb = boxes[nb].lo_k & ke = boxes[nb].hi_k & nbz = ke - kb + 1 & END
       2: BEGIN & jb = boxes[nb].lo_j & je = boxes[nb].hi_j & nby = je - jb + 1 & END
       1: BEGIN & ib = boxes[nb].lo_i & ie = boxes[nb].hi_i & nbx = ie - ib + 1 & END
     ENDSWITCH
     szb = LONG(nbx*nby*nbz*NVAR); ** box size

     ; ** rescale box bounds to current level's **

     SWITCH dim OF
       3: BEGIN & kb = kb*freb[nl] & ke = (ke+1)*freb[nl] - 1 & END
       2: BEGIN & jb = jb*freb[nl] & je = (je+1)*freb[nl] - 1 & END
       1: BEGIN & ib = ib*freb[nl] & ie = (ie+1)*freb[nl] - 1 & END
     ENDSWITCH

     ; ** skip boxes lying outside ranges **

     IF ( (ib GT iend) OR (ie LT ibeg) OR $
          (jb GT jend) OR (je LT jbeg) OR $
          (kb GT kend) OR (ke LT kbeg)) THEN BEGIN
       ncount += szb
       CONTINUE
     ENDIF

     IF (NOT KEYWORD_SET(NODATA)) THEN Q = REFORM(data(ncount:ncount+szb-1), nbx, nby, nbz, NVAR)

     ; ** find boxes intersection with current domain ranges **

     ib0 = max([ibeg, ib]) & ie0 = min([iend, ie])
     jb0 = max([jbeg, jb]) & je0 = min([jend, je])
     kb0 = max([kbeg, kb]) & ke0 = min([kend, ke])

     ; ** store box corners in the AMRLevel structure **
  
     IF (logr EQ 0) THEN BEGIN
      AMRLevel[nl].box[nb].x0 = x1b + dx*( ib0 ) 
      AMRLevel[nl].box[nb].x1 = x1b + dx*(ie0+1)
     ENDIF ELSE BEGIN
      AMRLevel[nl].box[nb].x0 = x1b*exp(( ib0 )*dx)
      AMRLevel[nl].box[nb].x1 = x1b*exp((ie0+1)*dx)
     ENDELSE

     AMRLevel[nl].box[nb].y0 = x2b + dx*( jb0 )*ystr 
     AMRlevel[nl].box[nb].y1 = x2b + dx*(je0+1)*ystr

     AMRLevel[nl].box[nb].z0 = x3b + dx*( kb0 )*zstr 
     AMRlevel[nl].box[nb].z1 = x3b + dx*(ke0+1)*zstr

     AMRLevel[nl].box[nb].ib = ib0 & AMRLevel[nl].box[nb].ie = ie0
     AMRLevel[nl].box[nb].jb = jb0 & AMRLevel[nl].box[nb].je = je0
     AMRLevel[nl].box[nb].kb = kb0 & AMRLevel[nl].box[nb].ke = ke0

     AMRBoxes[ib0-ibeg:ie0-ibeg, jb0-jbeg:je0-jbeg, kb0-kbeg:ke0-kbeg] = nl

     ; ** extract the box intersection from original data Q **

     IF (NOT KEYWORD_SET(NODATA)) THEN BEGIN
       cib0 = (ib0-ib)/freb[nl] & cie0 = (ie0-ib)/freb[nl]
       cjb0 = (jb0-jb)/freb[nl] & cje0 = (je0-jb)/freb[nl]
       ckb0 = (kb0-kb)/freb[nl] & cke0 = (ke0-kb)/freb[nl]
       ;Q1 = Q[cib0:cie0, cjb0:cje0, ckb0:cke0,*]
       Q1 = Q[cib0:cie0, cjb0:cje0, ckb0:cke0,varbeg:varend]

       ; ** remap the extracted portion onto vars **

       IF (KEYWORD_SET(INTERP)) THEN BEGIN
         FOR nv = varbeg, varend DO BEGIN
           vars(ib0-ibeg:ie0-ibeg, jb0-jbeg:je0-jbeg, kb0-kbeg:ke0-kbeg,nv) = $
           CONGRID(Q1(*,*,*,nv), ie0-ib0+1, je0-jb0+1, ke0-kb0+1,/INTERP)
         ENDFOR
       ENDIF ELSE BEGIN
         vars(ib0-ibeg:ie0-ibeg, jb0-jbeg:je0-jbeg, kb0-kbeg:ke0-kbeg,*) = $
    ;     REBIN(Q1, ie0-ib0+1, je0-jb0+1, ke0-kb0+1,NVAR)
         REBIN(Q1, ie0-ib0+1, je0-jb0+1, ke0-kb0+1,varend-varbeg+1,/sample)
       ENDELSE
     ENDIF

     ncount += szb

   ENDFOR;  ** on boxes **
   H5G_CLOSE,ilev

 ENDFOR; ** on levels **

 H5G_CLOSE,igen
 H5F_CLOSE,ifil

; *********************************************************
;  Build all necessary information for PLUTO common blocks
; *********************************************************

 IF (logr EQ 0) THEN BEGIN
  x1 = x1b + (ibeg + findgen(nx) + 0.5)*dx
 ENDIF ELSE BEGIN
  x1 = x1b*(exp((ibeg+findgen(nx)+1)*dx)+exp((ibeg+findgen(nx))*dx))*0.5
 ENDELSE
 x2 = x2b + (jbeg + findgen(ny) + 0.5)*dx*ystr
 x3 = x3b + (kbeg + findgen(nz) + 0.5)*dx*zstr

 NX1 = nx
 NX2 = ny
 NX3 = nz

 IF (logr EQ 0) THEN BEGIN
  dx1 = replicate(dx, NX1)
 ENDIF ELSE BEGIN
  dx1 = x1b*(exp((ibeg+findgen(nx)+1)*dx)-exp((ibeg+findgen(nx))*dx))
 ENDELSE
 dx2 = replicate(dx, NX2)*ystr
 dx3 = replicate(dx, NX3)*zstr

 IF (geometry EQ 'spherical') THEN BEGIN
     xpos = fltarr(NX1,NX2)  ; define projection on cartesian coordinates
     ypos = fltarr(NX1,NX2)
     for i = 0,NX1-1 do begin
     for j = 0,NX2-1 do begin
       xpos(i,j) = x1(i)*sin(x2(j))
       ypos(i,j) = x1(i)*cos(x2(j))
     endfor
     endfor
 ENDIF

 IF (geometry EQ 'polar') THEN BEGIN
     xpos = fltarr(NX1,NX2)  ; define projection on cartesian coordinates
     ypos = fltarr(NX1,NX2)
     for i = 0,NX1-1 do begin
     for j = 0,NX2-1 do begin
       xpos(i,j) = x1(i)*cos(x2(j))
       ypos(i,j) = x1(i)*sin(x2(j))
     endfor
     endfor
 ENDIF

; ** match variables with corresponding name **

 IF (KEYWORD_SET(NO_VAR_MATCH)) THEN RETURN 
 IF (KEYWORD_SET(NODATA))      THEN RETURN

; FOR nv = varbeg, varend DO MATCH_VARNAME, vars(*,*,*,nv-varbeg), vname[nv]
;return
;FOR nv =  0, nvar-1 DO BEGIN
;  MATCH_VARNAME, vars(*,*,*,0), vname[nv]
;  IF (nv LT NVAR-1) THEN vars = vars(*,*,*,1:nvar-1-nv)
;ENDFOR

; ** progressively free up memory by removing variables from vars **

 FOR nv = varbeg, varend DO BEGIN
   MATCH_VARNAME, vars(*,*,*,0), vname[nv], SILENT=SILENT
   IF (nv LT varend) THEN vars = vars(*,*,*,1:varend-varbeg-nv)
 ENDFOR
 vars = 0

END
