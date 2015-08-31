;+
;
;  NAME:      fourier 
;
;  AUTHOR:    Andrea Mignone
;
;  DATE:      Feb 10, 2005
;
;  PURPOSE:   Take the fourier transform of some signal and make
;             a log-log plot of its power in frequency space.
;
;  SYNOPSIS:  fourier, time, sp, nbeg, nend[, title = title]
;
;  ARGUMENTS / KEYWORDS:
;
;                        time:  a 1-D vector containing the
;                               times at which the signal 
;                               is sampled.
;
;                        sp  :  a 1-D vector containing the 
;                               sampled data.
;
;                        nbeg:  the initial index of the time
;                               window to be transformed.
;                               In other words, only the interval
;                               time(nbeg) < t < time(nend) 
;                               is considered.
;
;                        nend:  the final index. 
;
;                        title: an optional keyword specifying the 
;                               title of the plot
;
; VERSION:           beta
;
;-
pro FOURIER, time, sp, nbeg, nend, title = title

 IF (NOT KEYWORD_SET(title))  THEN title = 'Power Spectrum'

 tstring1 = strcompress(string(time(nbeg),format='(f5.1)'),/remove_all)
 tstring2 = strcompress(string(time(nend),format='(f5.1)'),/remove_all)
 twindow_string = tstring1+" < t < "+tstring2
 title = title + ", "+twindow_string

; ----------------------------------------------
;  compute power spectrum for the entire signal
; ----------------------------------------------

 print," > Nyquist critical frequency is:",0.5/(time(1)-time(0)) * 2.0*!PI

 ds = sp(nbeg:nend) - sp(0)  ; -- signal to be analyzed -- 

 Nds  = nend - nbeg + 1
 delt = (time(nend) - time(nbeg))/(Nds - 1)
 tg   = delt*findgen(Nds) + time(nbeg)
 U    = interpol(ds, time(nbeg:nend), tg)

 ;V = FFT(hanning(n)*U)
; V = FFT(hanning(Nds, alpha = 0.56)*U)
 V = FFT(U)

 tot_pow = abs(V(0:Nds/2))^2
 freq = FINDGEN(Nds/2+1) / (Nds*delt)

; -----------------------------------------
;    normalize power to its total 
; ----------------------------------------- 

 tot_pow = tot_pow/total(tot_pow)
 freq    = freq*2.0*!PI

;   -------------------------------------
;           Plot total power
;   -------------------------------------

 fmin = 2.0*!PI/(time(nend) - time(nbeg))
 fmax = 0.5/(time(1)-time(0)) * 2.0*!PI
 pmin = min(alog10(tot_pow))
 
 plot, freq, alog10(tot_pow),/xlog,$
       xtitle='!4x!3',xrange = [fmin, fmax], xstyle = 1, $
       yrange = [pmin, 0.0] , ystyle = 1, $
       ytitle = 'Log(Normalized Power)', $
;       ytickformat='(f4.0)',$
       chars = 1.2, title = title


end
