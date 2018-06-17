;+
;
; NAME:      PARTICLES_BINARY_LOAD
;
; AUTHOR:    Andrea Mignone
;
; PURPOSE:
;
;   Read PLUTO particle file in .flt or .dbl binary data format. 
;   Called from PARTICLES_LOAD with the /FLOAT (or no) keywords.
;   
; SYNTAX:
;
;   PARTICLES_BINARY_LOAD, fname[, /FLOAT][, /SILENT], endianity, err
;   
; ARGUMENTS:
;
;   fname      The name of a valid .flt or .dbl particle file.
;
; KEYWORDS:
;
;
; EXAMPLES:
;
;  * Example #1: Load the 10-th particles output file:
;
;    IDL> PARTICLES_LOAD, 10
;
;
; LAST MODIFIED:
;
;  May 30, 2017 by A. Mignone (mignone@ph.unito.it)
;
;-
PRO PARTICLES_BINARY_LOAD, fname, FLOAT=FLOAT, SILENT=SILENT, err

  COMMON PLUTO_GRID
  COMMON PLUTO_VAR
  COMMON PLUTO_RUN

; --------------------------------------------
; 1. Read header file section
; --------------------------------------------

  OPENR, unit, fname, /GET_LUN, ERROR = err
  IF (err NE 0) THEN BEGIN
    PRINT,"! Cannot open file "+fname
    err = 1
    RETURN
  ENDIF

  header = ''
  READF,unit,header
PRINT,HEADER
  s = STRMID(header,0,1); Get first character of the string
  WHILE (s EQ '#') DO BEGIN ;  Loop on header strings
    POINT_LUN, -unit, pos ;  Get current file position

  ; -- Extract dimensions (important for reduced format) --
    substr = STRSPLIT(header,/EXTRACT,count=count)
PRINT,count    
    FOR i=0,count-1 DO BEGIN
      IF (substr[i] EQ "dimensions") THEN BEGIN
PRINT,"i=",i      
        dimensions = FIX(substr[i+1])
        BREAK
      ENDIF
    ENDFOR

  ; -- Extract number of particles --
    substr = STRSPLIT(header,/EXTRACT,count=count)
    FOR i=0,count-1 DO BEGIN
      IF (substr[i] EQ "nparticles") THEN BEGIN
PRINT,"i=",i      
        nparticles = ULONG64(substr[i+1])
        BREAK
      ENDIF
    ENDFOR

  ; -- Extract endianity --
    substr = STRSPLIT(header,/EXTRACT,count=count)
    FOR i=0,count-1 DO BEGIN
      IF (substr[i] EQ "endianity") THEN BEGIN
PRINT,"i=",i      
        endianity = substr[i+1]
        BREAK
      ENDIF
    ENDFOR

    READF,unit,header     ;  Read next string
    s = STRMID(header,0,1);  Get first character
  ENDWHILE


STOP

  POINT_LUN,unit,pos;  Set file pointer to beginning of binary section
    
; --------------------------------------------
; 0. Print some info
; --------------------------------------------

  IF (NOT KEYWORD_SET(SILENT)) THEN BEGIN
    PRINT,"> Binary data file: ",fname
    IF (KEYWORD_SET (FLOAT)) THEN BEGIN
      PRINT,"> Precision: float"
    ENDIF ELSE BEGIN
      PRINT,"> Precision: double"
    ENDELSE
    PRINT,"> Endianity:  ",endianity
    PRINT,"> nparticles: ", ARG2STR(nparticles)
  ENDIF
  
; --------------------------------------------------------
; 2. Allocate memory for an array of particles structures
; --------------------------------------------------------

  particles = 0.0;  Free memory (?)
  particles = REPLICATE({particle_struct},nparticles)
  nfields   = 11;  Number of structure members to be read

; ---------------------------------------------------------------
; 3a. Read data (using assoc):
;     For a full particle:
;       id / x[3] / v[3] / color / tinj
; ---------------------------------------------------------------

  POINT_LUN,-unit,pos
  IF (KEYWORD_SET(FLOAT)) THEN BEGIN
    Q = ASSOC(unit,FLTARR(nfields,nparticles),pos)
  ENDIF ELSE BEGIN
    Q = ASSOC(unit,DBLARR(nfields,nparticles),pos)
  ENDELSE

  ; -- Swap endianity and close file --

  Q = Q[0]
  IF (endianity EQ "little") THEN BEGIN
    Q = SWAP_ENDIAN (Q, /SWAP_IF_BIG_ENDIAN)
  ENDIF ELSE BEGIN
    Q = SWAP_ENDIAN (Q, /SWAP_IF_LITTLE_ENDIAN)
  ENDELSE

  IF (format EQ "full") THEN BEGIN
    particles.(0)  = REFORM(Q[0,*]);   particles id
    particles.(1)  = REFORM(Q[3,*]);   particles x1-coordinate
    particles.(2)  = REFORM(Q[4,*]);   particles x2-coordinate
    particles.(3)  = REFORM(Q[5,*]);   particles x3-coordinate
    particles.(6)  = REFORM(Q[6,*]);   particles x1-velocity
    particles.(7)  = REFORM(Q[7,*]);   particles x2-velocity
    particles.(8)  = REFORM(Q[8,*]);   particles x3-velocity
    particles.(9)  = REFORM(Q[9,*]);   particles color
    particles.(10) = REFORM(Q[10,*]); particles injection time
  ENDIF

  CLOSE,unit
  FREE_LUN,unit
  err = 0
END

;+
;
; NAME:      PARTICLES_VTK_LOAD
;
; AUTHOR:    Andrea Mignone
;
; PURPOSE:
;
;   Read PLUTO particles file in .vtk data format. 
;   Called from PARTICLES_LOAD with the /VTK keyword.
;   
; SYNTAX:
;
;   PARTICLES_VTKLOAD, fname
;   
; ARGUMENTS:
;
;   fname      The name of a valid .vtk paticle file.
;
; KEYWORDS:
;
;
; EXAMPLES:
;
;  * Example #1: Load the 10-th particles output file:
;
;    IDL> PARTICLES_LOAD, /VTK,10
;
;
; LAST MODIFIED:
;
;  May 30, 2017 by A. Mignone (mignone@ph.unito.it)
;
;-
PRO PARTICLES_VTK_LOAD, fname, SILENT=SILENT, err

  COMMON PLUTO_GRID
  COMMON PLUTO_VAR
  COMMON PLUTO_RUN

  T0 = SYSTIME(1)
    
; ---------------------------------------------------------
; 0. Print some info
; ---------------------------------------------------------

  PRINT,"> Reading vtk data format"

; ---------------------------------------------------------
; 1. Open file for reading
; ---------------------------------------------------------

  OPENR, unit, fname, /GET_LUN, ERROR = err, /SWAP_IF_LITTLE_ENDIAN
  IF (err NE 0) THEN BEGIN
    PRINT,"! Cannot open file "+fname
    err = 1
    RETURN
  ENDIF
  endianity = "big"
  
; ---------------------------------------------------------
; 2. Start reading data until end of file
; ---------------------------------------------------------

  WHILE ~ EOF(unit) DO BEGIN
    sline = ''
    READF, unit, sline
    words = STRSPLIT(sline," ", /EXTRACT)
    IF (SIZE(words,/N_ELEMENT) LE 1) THEN CONTINUE
    
    ; -- Read particles positions -- 
    
    IF (words[0] EQ "POINTS") THEN BEGIN
      nparticles = ULONG64(words[1])
      particles = 0.0;
      IF (nparticles EQ 0) THEN RETURN
      particles = REPLICATE({particle_struct},nparticles)

      POINT_LUN, -unit, offset; obtain the current position in the file
      Q     = ASSOC(unit,FLTARR(3,nparticles),offset)
      Q     = Q[0]
      particles.(3) = REFORM(Q[0,*]);   particles x1-coordinate
      particles.(4) = REFORM(Q[1,*]);   particles x2-coordinate
      particles.(5) = REFORM(Q[2,*]);   particles x3-coordinate
    ENDIF

    ; -- Read particles identity --

    IF (words[0] EQ "SCALARS" AND words[1] EQ "Identity") THEN BEGIN
      READF,unit,sline;         skip next line ("LOOKUP_TABLE")
      POINT_LUN, -unit, offset; obtain the current position in the file
      Q     = ASSOC(unit,LONARR(nparticles),offset)
      Q     = Q[0]
      particles.(0) = REFORM(Q[*]);   particles identity
    ENDIF

    ; -- Read particles birth proc --

    IF (words[0] EQ "SCALARS" AND words[1] EQ "Birth_Proc") THEN BEGIN
      READF,unit,sline;         skip next line ("LOOKUP_TABLE")
      POINT_LUN, -unit, offset; obtain the current position in the file
      Q     = ASSOC(unit,LONARR(nparticles),offset)
      Q     = Q[0]
      particles.(1) = REFORM(Q[*]);   particles birth proc
    ENDIF

    ; -- Read particles restart id --

    IF (words[0] EQ "SCALARS" AND words[1] EQ "RestartID") THEN BEGIN
      READF,unit,sline;         skip next line ("LOOKUP_TABLE")
      POINT_LUN, -unit, offset; obtain the current position in the file
      Q     = ASSOC(unit,LONARR(nparticles),offset)
      Q     = Q[0]
      particles.(2) = REFORM(Q[*]);   particles restart id
    ENDIF

    ; -- Read particle velocity --

    IF (words[0] EQ "VECTORS" AND words[1] EQ "Velocity") THEN BEGIN
      POINT_LUN, -unit, offset; obtain the current position in the file
      Q     = ASSOC(unit,FLTARR(3, nparticles),offset)
      Q     = Q[0]
      particles.(6) = REFORM(Q[0,*]);   particles x1 velocity
      particles.(7) = REFORM(Q[1,*]);   particles x2 velocity
      particles.(8) = REFORM(Q[2,*]);   particles x3 velocity
    ENDIF

    ; -- Read particle Four-velocity --

    IF (words[0] EQ "VECTORS" AND words[1] EQ "Four-Velocity") THEN BEGIN
      POINT_LUN, -unit, offset; obtain the current position in the file
      Q     = ASSOC(unit,FLTARR(3, nparticles),offset)
      Q     = Q[0]
      particles.(6) = REFORM(Q[0,*]);   particles x1 four-velocity
      particles.(7) = REFORM(Q[1,*]);   particles x2 four-velocity
      particles.(8) = REFORM(Q[2,*]);   particles x3 four-velocity
    ENDIF

    ; -- Read particles color --

    IF (words[0] EQ "SCALARS" AND words[1] EQ "Color") THEN BEGIN
      READF,unit,sline;         skip next line ("LOOKUP_TABLE")
      POINT_LUN, -unit, offset; obtain the current position in the file
      Q     = ASSOC(unit,FLTARR(nparticles),offset)
      Q     = Q[0]
      particles.(9) = REFORM(Q[*]);   particles color
    ENDIF

    ; -- Read particles injection time --

    IF (words[0] EQ "SCALARS" AND words[1] EQ "t_inj") THEN BEGIN
      READF,unit,sline;         skip next line ("LOOKUP_TABLE")
      POINT_LUN, -unit, offset; obtain the current position in the file
      Q     = ASSOC(unit,FLTARR(nparticles),offset)
      Q     = Q[0]
      particles.(10) = REFORM(Q[*]);   particles injection time
    ENDIF

  ENDWHILE

; ---------------------------------------------------------
; 3. Print some info
; ---------------------------------------------------------

  IF (NOT KEYWORD_SET(SILENT)) THEN BEGIN
    PRINT,"> Binary data file: ",fname
    IF (KEYWORD_SET (FLOAT)) THEN BEGIN
      PRINT,"> Precision: float"
    ENDIF ELSE BEGIN
      PRINT,"> Precision: double"
    ENDELSE
    PRINT,"> Endianity:  ",endianity
    PRINT,"> nparticles: ", ARG2STR(nparticles)
  ENDIF

  CLOSE,unit
  FREE_LUN,unit
  err = 0
END

;+
;
; NAME:      PARTICLES_LOAD
;
; AUTHOR:    Andrea Mignone
;
; PURPOSE:  
;
;   Read PLUTO particles file in .dbl, .flt and .vtk format.
;   By default, if neither /VTK or /FLOAT are given, double precision 
;   raw-binary files are read (particles.nnnn.dbl). 
;   
;   Data is read and stored in the 'particle' array of structures which has
;   variables size depending on the actual number of particles (nparticles).
;   Both particle and nparticles can be shared among different routines using
;   the PLUTO_VAR common block.    
;
; SYNTAX:
;
;   PARTICLES_LOAD, nout[, DIR=DIR][, SILENT=SILENT][, FLOAT=FLOAT][, VTK = VTK]
;
; ARGUMENTS:
;
;   nout = the output number, e.g., 0,1,2, ...
;
; KEYWORDS:
;
;   DIR        the directory name where output files are located (e.g. 
;              DIR = 'Data1/'). Default is current directory.
;
;   /SILENT    set this keyword to suppress printing information
;
;   /FLOAT     set this keyword to read single precision binary datafile
;              (particles.nnnn.flt). 
;              
;   /VTK       set this keyword to read .vtk particle files.
;              In this case, the PARTICLES_VTKLOAD procedure is invoked.
;
;   /REDUCED   Assume particle data files have been written with
;              the PARTCILES_REDUCE_FLT switch turned to YES.
;              In this case, the particle structure will contain
;              fewer fields.
;
; EXAMPLES:
;
;  * Example #1: Load the 7-th particles output file in single precision
;                and print the x2-component of velocity of the 4-th particle:
;
;    IDL> PARTICLES_LOAD, /FLOAT, 10
;    IDL> PRINT,particles[4].vx2
;
;
; LAST MODIFIED:
;
;  May 30, 2017 by A. Mignone (mignone@ph.unito.it)
;
;-
PRO PARTICLES_LOAD, nout, dir=dir, SILENT=SILENT, FLOAT=FLOAT, VTK=VTK,$
                    REDUCED=REDUCED 

  COMMON PLUTO_GRID
  COMMON PLUTO_VAR
  COMMON PLUTO_RUN

  T0 = SYSTIME(1)
  IF (NOT KEYWORD_SET(DIR)) THEN dir = "."

; ---------------------------------------------------------
; 0. Define Particle structure here.
;    The main particle structure contains 11 fields placed
;    in the following order:
;
;    id ,  x1,  x2,  x3, vx1, vx2, vx3, color, tinj
;    (0),  (1), (2), (3), (4), (5), (6), (7),   (8)
; ---------------------------------------------------------

  struct1 = {particle_basic, id: 0.d0, x1:  0.d0,  x2: 0.d0, x3: 0.d0,    $
                                       vx1: 0.d0, vx2: 0.d0, vx3: 0.d0,  $
                                       color:0.d0, tinj:0.d0}
  struct2 = {particle_struct, INHERITS particle_basic, rho: 0.d0}

; ---------------------------------------------------------
; 1. Build file name & call corresponding
;    loader
; ---------------------------------------------------------

  str_nout = STRCOMPRESS(STRING(nout,FORMAT='(I4.4)'))
  fname  = dir+"/particles."+str_nout

  IF (KEYWORD_SET(FLOAT) OR KEYWORD_SET(FLT)) THEN BEGIN
    PARTICLES_BINARY_LOAD,fname+".flt",/FLOAT,SILENT=SILENT, err
  ENDIF ELSE IF (KEYWORD_SET(VTK)) THEN BEGIN
    PARTICLES_VTK_LOAD,fname+".vtk",SILENT=SILENT, err
  ENDIF ELSE BEGIN
    PARTICLES_BINARY_LOAD,fname+".dbl",SILENT=SILENT, err
  ENDELSE  
    
  IF (err NE 0) THEN RETURN

;  PRINT,"Function took ",SYSTIME(1)-T0, "seconds"
END
