;+
;
; NAME:      VTK_LOAD
;
; AUTHOR:    A. Mignone
;
; PURPOSE:   Read PLUTO variables stored in a VTK data file. 
;            It is automatically called by PLOAD when the /VTK keyword
;            is supplied.
;
; SYNTAX:    VTK_LOAD, filename[, KSLICE=integer][, /SILENT][, VAR=string]
;            [, {X1 | X2 | X3}RANGE=[left,right]]
;
; ARGUMENTS
;
;   filename   the data file name
;
; KEYWORDS 
;
;   KSLICE        Set this keyword to an integer specifying the array index
;                 of the 3rd coordinate. This will read and store a 2D array
;                 (e.g. rho[*,*] at z = z[k]) rather then the full 3D array.
;                 Useful for large datasets.
;
;   /SILENT       Set this keyword to suppress output.
;
;   VAR           the name of a particular variable to be loaded, e.g., 
;                 "prs". The remainig ones are ignored.
;
;   {X1 | X2 | X3}RANGE    A 2 elements array giving the desired data range of 
;                          the axis.
;                          
;
; LAST MODIFIED
;
;   Dec 5, 2013 by A. Mignone (mignone@ph.unito.it)
;
;-
PRO VTK_LOAD, filename, KSLICE=KSLICE, SILENT=SILENT, VAR=VAR,$
              X1RANGE=X1RANGE, X2RANGE=X2RANGE, X3RANGE=X3RANGE, domain_indx

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

 ; -------------------------------------
 ;  check keywords 
 ; -------------------------------------
 
 domain_indx = 0
 get_kslice = N_ELEMENTS(KSLICE); = 0 or on 1 depending on whether 
                                ; kslice keyword is present 
 IF (NOT KEYWORD_SET(VAR)) THEN var = 'all'; means "read all variables"

 extract_sub_arr = 0
 IF (KEYWORD_SET (X1RANGE) OR KEYWORD_SET (X2RANGE) OR KEYWORD_SET (X3RANGE)) THEN BEGIN 
   IF (NOT KEYWORD_SET (X1RANGE)) THEN x1range=[MIN(x1),MAX(x1)]
   IF (NOT KEYWORD_SET (X2RANGE)) THEN x2range=[MIN(x2),MAX(x2)]
   IF (NOT KEYWORD_SET (X3RANGE)) THEN x3range=[MIN(x3),MAX(x3)]
   extract_sub_arr = 1
 ENDIF

 ndim = 3
 IF (nx3 EQ 1) THEN ndim = 2
 
 ; -------------------------------------
 ;  Open for reading
 ; -------------------------------------

 GET_LUN, unit
 OPENR, unit, filename, /SWAP_IF_LITTLE_ENDIAN

 data_size = ULONG64(nx1)*ULONG64(nx2)*ULONG64(nx3)*4; in bytes

 ; -------------------------------------------------
 ;  Read header file to get the dimensions
 ;  (this part can be skipped since grid
 ;  information has already been read previously
 ;  from PLOAD)   
 ; -------------------------------------------------

; WHILE ~ EOF(unit) DO BEGIN
;   sline = ''
;   READF, unit, sline
;   words = STRSPLIT(sline," ", /EXTRACT)
;   IF (words[0] EQ "DIMENSIONS") THEN BEGIN
;     PRINT,sline
;     nx1 = FIX(words[1])-1
;     nx2 = FIX(words[2])-1
;     nx3 = FIX(words[3])-1
;     BREAK
;   ENDIF
; ENDWHILE
 
 ;IF (nx3 EQ 0) THEN nx3 = 1
 ; PRINT,"> nx1, nx2, nx3 = ",nx1,nx2,nx3

 ; -------------------------------------------------
 ;         Read scalar and vector data
 ; -------------------------------------------------

 WHILE ~ EOF(unit) DO BEGIN

   sline = ''
   READF, unit, sline
   words = STRSPLIT(sline," ", /EXTRACT)
   IF (SIZE(words,/N_ELEMENT) LE 1) THEN CONTINUE
   
 ; -----------------------------------------------------------------------------
 ;                    Read Scalar Data
 ; -----------------------------------------------------------------------------

   IF (words[0] EQ "SCALARS") THEN BEGIN

     READF, unit, sline; read one more line (i.e., skip "LOOKUP_TABLE")
     vname     = words[1]    ; the variable name after "SCALARS"
     POINT_LUN, -unit, offset; obtain the current position in the file
     read_var  = (var EQ 'all') OR (var EQ vname); read all or just one variable ? 
  
     CASE 1 OF 

       (~read_var): BEGIN
           ; -- skip this file section if there's no variable to read --
       END
       
       get_kslice: BEGIN;  read a 2D slice at constant k
         vpt = ASSOC(unit, FLTARR(nx1,nx2), offset)
         MATCH_VARNAME, vpt[KSLICE], vname, silent=silent
       END
       
       extract_sub_arr: BEGIN 
         vpt = ASSOC(unit, FLTARR(nx1,nx2), offset)
         ARRAY_EXTRACT, vpt, x1range, x2range, x3range,domain_indx,silent=silent
         MATCH_VARNAME, vpt, vname, silent=silent
       END
       
       ELSE: BEGIN; read the whole record
         vpt = ASSOC(unit, FLTARR(nx1,nx2,nx3), offset)


       ; ------------------------------------------------------
       ;  Version 1:    standard call
       ; ------------------------------------------------------
         
         MATCH_VARNAME, vpt[0], vname, silent=silent
         vpt = 0

       ; ------------------------------------------------------
       ;  Version 2:    avoid using MATCH_VARNAME 
       ;               (which is quite memory consuming)
       ; ------------------------------------------------------

;         CASE vname OF
;           "rho":  rho = TEMPORARY(vpt[0])
;           "prs":  prs = TEMPORARY(vpt[0])
;           "tr1":  tr1 = TEMPORARY(vpt[0])
;         ELSE: MATCH_VARNAME, vpt[0], vname, silent=silent
;         ENDCASE

       END
       
     ENDCASE
   
     POINT_LUN, unit, offset + data_size;  move file pointer to the end
                                        ;  of this record 
   ENDIF; words[0] EQ "SCALARS"
   

 ; -----------------------------------------------------------------------------
 ;                     Read Vector Data
 ; -----------------------------------------------------------------------------

   IF (words[0] EQ "VECTORS") THEN BEGIN

     vname    = words[1]
     POINT_LUN, -unit, offset; obtain the current position in the file

   ; -- check how many variables do we need to read --

     IF (STRMATCH(vname,'*Velocity*')) THEN component = ['vx1', 'vx2', 'vx3']
     IF (STRMATCH(vname,'*Magnetic*')) THEN component = ['bx1', 'bx2', 'bx3']
     
     var_indx = WHERE(component EQ var)
     read_var = (var EQ 'all') OR (var_indx GE 0)
       
     nbeg = 0  
     nend = 2
     IF (read_var AND var NE 'all') THEN BEGIN; read either all three
                                              ; components (0,1,2) or just one
       nbeg = var_indx[0]
       nend = var_indx[0]
     ENDIF
     
     CASE 1 OF

       (~read_var): BEGIN
           ; -- skip this file section if there's no variable to read --
       END
           
       get_kslice: BEGIN;  read a 2D slice at constant k
         vpt = ASSOC(unit, FLTARR(3, nx1, nx2), offset)
         vpt = vpt[KSLICE]
         FOR n = nbeg, nend DO MATCH_VARNAME, REFORM(vpt[n,*,*]), $
                                              component[n], silent=silent
       END
     
       extract_sub_arr: BEGIN
         vpt = ASSOC(unit, FLTARR(3, nx1, nx2), offset)
         ARRAY_EXTRACT, vpt, x1range, x2range, x3range,domain_indx,$
                        silent=silent,/VECTOR
         FOR n = nbeg, nend DO MATCH_VARNAME, REFORM(vpt[n,*,*,*]), $
                                              component[n], silent=silent
       END
              
       ELSE: BEGIN; read the whole file

version = 1
       ; ------------------------------------------------------
       ; Version 1: standard version
       ; ------------------------------------------------------

IF (version EQ 1) THEN BEGIN
         vpt = ASSOC(unit, FLTARR(3, nx1, nx2, nx3), offset)
         vpt = vpt[0]
         FOR n = nbeg, nend DO MATCH_VARNAME, REFORM(vpt[n,*,*,*],/OVERWRITE),$
                                              component[n], silent=silent
ENDIF

       ; ------------------------------------------------------
       ;  Version 2: avoid using MATCH_VARNAME for memory
       ;             saving
       ; ------------------------------------------------------

IF (version EQ 2) THEN BEGIN
         vpt = ASSOC(unit,FLTARR(3,nx1,nx2), offset)
         FOR n = nbeg, nend DO BEGIN
           Q = FLTARR(nx1,nx2,nx3)
           FOR k = 0,nx3-1 DO Q[*,*,k] = vpt(n,*,*,k)
           IF (component[n] EQ "vx1") THEN vx1 = TEMPORARY(Q)
           IF (component[n] EQ "vx2") THEN vx2 = TEMPORARY(Q)
           IF (component[n] EQ "vx3") THEN vx3 = TEMPORARY(Q)
           IF (component[n] EQ "bx1") THEN Bx1 = TEMPORARY(Q)
           IF (component[n] EQ "bx2") THEN Bx2 = TEMPORARY(Q)
           IF (component[n] EQ "bx3") THEN Bx3 = TEMPORARY(Q)
           Q = 0
         ENDFOR 
ENDIF

       ; ------------------------------------------------------
       ;              Version 3
       ; ------------------------------------------------------

IF (version EQ 3) THEN BEGIN
         vpt = FLTARR(3,nx1)
         FOR n = nbeg, nend DO BEGIN
           Q   = FLTARR(nx1,nx2,nx3)
           FOR k = 0,nx3-1 DO BEGIN
           FOR j = 0,nx2-1 DO BEGIN
             READU,unit,vpt
             Q(*,j,k) = vpt[n,*]
           ENDFOR
           ENDFOR
           
           MATCH_VARNAME, Q, component[n], silent=silent
           Q = 0
         ENDFOR
         vpt = 0
ENDIF
       END
         
     ENDCASE
     
     POINT_LUN, unit, offset + data_size*3;  move file pointer to the end
                                          ;  of this record 
   ENDIF

   vpt = 0
 ENDWHILE
 CLOSE, unit
 FREE_LUN,unit

END
