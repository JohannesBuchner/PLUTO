;+
; NAME:      DISTRIB_FUNCTION
;
; AUTHOR:    A. Mignone (mignone@ph.unito.it)
;
; PURPOSE:   Compute the density function of an array E[i] using
;            the HISTOGRAM function.
;            On output, fE[i]=dN[i]/dEbin (the density array) and the bin
;            location array Ebin[i] (at bin center) are returned.
;            Bins are equally spaced in the linear coordinate, or in
;            log space (if the /log keyword is supplied).
;            Bin spacing is optionally returned by supplying the BINSIZE
;            keyword.
;
;            Note that the integral of fE*dEbin should equal the total number
;            of elements of the array E between Emin and Emax.
;            If Emin and Emax enclose the entire range (or if these keywords
;            are omitted), then TOTAL(fE*dEbin) = N_ELEMENTS(E).
;
; SYNTAX:    DISTRIB_FUNCTION, E, fE, Ebins, dEbins [,/LOG]
;                              [,MIN=MIN][,MAX=MAX][,NBINS=NBINS]
;
; ARGUMENTS:
;
;   E        (input)  The array for which the density function is required.
;   fE       (output) The density function of E[].
;   Ebins    (output) An array containing the center location of each bin.
;   dEbins   (output) An array containing the bin widths.
;   
;   
; KEYWORDS:
;
;   LOG      Enable this keyword to compute equally spaced bins in log space.
;            The location of the bins (left-interface) will then be
;
;              E[i] = 10^Xi[i],     where   Xi[i] = Log(Emin) + i*dXi,
;
;            and i = 0,...nbins-1, while
;     
;              dXi = (Log(Emax) - Log(Emin))/nbins
;
;            The linear spacing is determined as
;
;              dE[i] = 10^Xi[i]*(10^dXi - 1)
;                    = E[i]*[(Emax/Emin)^q - (Emin/Emax)^q] 
;
;            where q = 1/(2*nbins).
;
;   MIN      Set this keyword to the minimum value to consider.
;
;   MAX      Set this keyword to the maximum value to consider.
;            No values will be consider beyond this limit (note that IDL
;            interprets "max" as the *starting* location of the last bin).
;
;   NBINS    The number of bins used to compute the distribution function.
;
;
; EXAMPLES:
;
;   * Example #1: Compute the density array of v
;
;     IDL> DISTRIB_FUNCTION, v, df, vbins, nbins=12
;     IDL> PLOT, vbins, df
;     
;   * Example #2: Compute the density array of v, compare totals:
;
;     IDL> DISTRIB_FUNCTION, a, fa, abins, BINSIZE=da_bins, nbins=nbins, /log
;     IDL> PRINT,"Total(fa*da) = ",N_ELEMENTS(a)
;     IDL> PRINT,"Total(fa*da) = ",TOTAL(fa*da_bins)
;
;
; LAST MODIFIED:  Oct 13, 2017 by A. Mignone
;-
PRO DISTRIB_FUNCTION, E, fE, Ebins, dEbins, nbins=nbins,$
                      MIN=Emin, MAX=Emax, LOG=LOG
 
  IF (NOT KEYWORD_SET(nbins)) THEN nbins = 10
  
  IF (N_ELEMENTS(Emin) EQ 0)  THEN Emin  = MIN(E)
  IF (N_ELEMENTS(Emax) EQ 0)  THEN Emax  = MAX(E)


  IF (KEYWORD_SET(log)) THEN BEGIN

    log_Emax = ALOG10(Emax)
    log_Emin = ALOG10(Emin)

    bin_width = (log_Emax - log_Emin)/DOUBLE(nbins)

  ; ------------------------------------------------------------
  ;  Call HISTOGRAM with max = Emax-bin_width since IDL
  ;  interprets "max" as the starting location of
  ;  the last bin
  ; ------------------------------------------------------------

    dN  = HISTOGRAM (ALOG10(E), LOCATIONS=logE, $
                     MIN=log_Emin, MAX=log_Emax-bin_width,$
                     NBINS=nbins, /L64)

    dlogE  = logE[1] - logE[0]        ; Uniform log spacing
    dEbins = 10^logE*(10^dlogE - 1.0) ; Spacing in linear coordinate
    Ebins  = 10^(logE+0.5*dlogE)      ; Bin coordinate
    fE     = (1.0*dN)/dEbins          ; Trasform to float

;PRINT,"logE  = ",logE
;PRINT,"dlogE = ",dlogE
;PRINT,"Ebins = ",Ebins
;PRINT,"dE    = ",dEbins
;PRINT,"dE    = ",Ebins*( (Emax/Emin)^(1.0/(2.0*nbins))   $
;                        -(Emax/Emin)^(-1.0/(2.0*nbins)) )
;PRINT,"dN    = ",dN*1.0

  ENDIF ELSE BEGIN

    bin_width = (Emax-Emin)/DOUBLE(nbins)

  ; ------------------------------------------------------------
  ;  Call HISTOGRAM with max = Emax-bin_width since IDL
  ;  interprets "max" as the starting location of
  ;  the last bin
  ; ------------------------------------------------------------
  
    dN     = HISTOGRAM (E, LOCATIONS=Ebins, $
                        MIN=Emin, MAX=Emax-bin_width, NBINS=nbins, /L64)
    dEbins = Ebins[1] - Ebins[0]; Spacing in linear coordinate
    Ebins  = Ebins + 0.5*dEbins ; Bin center
    fE     = (1.0*dN)/dEbins    ; Trasform to float

;PRINT,"Ebins  = ",Ebins
;PRINT,"dEbins = ",dEbins
;PRINT,"dN     = ",fE*dEbins
  ENDELSE
  
END


