;+
; Check if filename is a valid HDF5 file
;-
PRO CHECKFILE,filename

 finfo = file_info(filename)

 IF (finfo.exists eq 0) THEN BEGIN
   print,'! File '+strcompress(filename,/REMOVE_ALL)+' does not exist'
   RETALL
 ENDIF ELSE BEGIN
   ishdf = H5F_IS_HDF5(filename)
   IF (ishdf eq 0) THEN BEGIN
     print,'! File '+strcompress(filename,/REMOVE_ALL)+' is not an HDF5 file'
     RETALL
   ENDIF
 ENDELSE
 RETURN
END
