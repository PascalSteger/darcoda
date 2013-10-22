;============================ BINCOUNT ============================
;This routine takes an array, r, and counts the number of elements 
;in r bins of size bin. 

;WARNING!! THIS ROUTINE REQUIRES SORTED ACSENDING r ARRAYS.

PRO bincount,r,low,high,bin,rout,arrayout,count_bin

COMMON constants

Pnts = ROUND((high-low)/bin)+1
rout = indgen(Pnts)*bin + low
arrayout = fltarr(Pnts)
count_bin = fltarr(Pnts)
error = fltarr(Pnts)
j=0L
size = N_ELEMENTS(r)

FOR i=0L,Pnts-1 DO BEGIN
   WHILE (rout(i) GT r(j)) DO BEGIN
      arrayout(i)=arrayout(i)+1.
      IF (j LT size-1) THEN BEGIN 
         j=j+1 
      ENDIF ELSE GOTO, SKIP
   ENDWHILE
   SKIP:
   count_bin(i) = arrayout(i)

ENDFOR

END
