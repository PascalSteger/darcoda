;============================= ARRAYDUMP ===============================
;This routine takes a number, narr, of equal length arrays:
;arrays=[[arr1], [arr2]...] and dumps them to a specified
;file (fname) in columnated data format.

PRO arraydump,fname,arrays,narr

;Make sure the file unit is closed!!
close,12

;Find length of each column ( = no. of rows):
totsize = N_ELEMENTS(arrays)
lenarr = totsize/narr

;Set string output type:
fmt='(A)'

;Open the output file and send the columnated data:
openw,12,fname
FOR i=0L,lenarr-1 DO BEGIN
    outstr = ''
    FOR j=i,totsize-1,lenarr DO BEGIN
        IF (j LT totsize-lenarr) THEN $ 
          outstr = outstr+strcompress(arrays(j),/remove_all)+' ' $
        ELSE $
          outstr = outstr+strcompress(arrays(j),/remove_all)
    ENDFOR
    printf,12,outstr,format=fmt
ENDFOR
close,12

END
