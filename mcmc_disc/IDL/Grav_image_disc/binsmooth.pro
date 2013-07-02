;============================ BINSMOOTH ============================
;This routine takes an array(r) and bins it in r bins of size bin. If
;there is no data in a particular bin, it assigns a value of
;array(r)=nanflag.

;WARNING!! THIS ROUTINE REQUIRES SORTED ACSENDING r ARRAYS.

PRO binsmooth,r,array,low,high,bin,rout,arrayout,nanflag,count_bin

COMMON constants

Pnts = ROUND((high-low)/bin)+1
rout = indgen(Pnts)*bin + low
arrayout = fltarr(Pnts)
count_bin = fltarr(Pnts)
j=0L
size = N_ELEMENTS(r)

FOR i=0L,Pnts-1 DO BEGIN
    count=0L
    WHILE (rout(i) GT r(j)) DO BEGIN
        arrayout(i)=arrayout(i)+array(j)
        IF (j LT size-1) THEN BEGIN 
            j=j+1 
            IF (array(j) NE 0) THEN count=count+1
        ENDIF ELSE GOTO, SKIP
    ENDWHILE
SKIP:
    IF (count GT 0) THEN arrayout(i)=arrayout(i)/count ELSE $
       arrayout(i)=nanflag
    count_bin(i) = count
ENDFOR

END
