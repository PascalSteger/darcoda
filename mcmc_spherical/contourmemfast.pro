;=========================== CONTOURMEMFAST ============================
;This routine takes a distribution of points (x,y) contained in two
;(SAME SIZE!) arrays x and y and calculates an array denxy(xnew,ynew) of
;point densities binned by radius in bins of size dx*dy. dy is
;calculated such that denxy is a square array as required by the IDL
;contour routine. It also generates length arrays xnew and ynew which
;correspond to denxy(xnew,ynew). 

PRO contourmemfast,dx,xmax,xmin,ymax,ymin,x,y,mass,xnew,ynew,denxy

COMMON constants

;Set up length arrays:
xpnts = (xmax-xmin)/dx+1
ypnts = xpnts
dy = (ymax-ymin)/(ypnts-1)
xnew = indgen(xpnts)*dx+xmin
ynew = indgen(ypnts)*dy+ymin

;Set up contour density array
denxy = fltarr(xpnts,ypnts)
 
;Do the contour calc:
pnts = N_ELEMENTS(x)

;Cut x and y to the domain size:
xc = x > xmin
xc = xc < xmax
yc = y > ymin
yc = yc < ymax
print,'Doing contour calc:'

FOR i=0L,pnts-1 DO BEGIN
    placex = (xc(i)-xmin)/dx
    placey = (yc(i)-ymin)/dy
    IF (placex NE xmin) AND (placex NE xmax) $
      AND (placey NE ymin) AND (placey NE ymax) THEN $
      denxy(placex,placey) = denxy(placex,placey) + mass(i)
ENDFOR
print,'Done'

;Fix sample bug i.e. centre on centre of cell:
xnew = xnew + dx/2.
ynew = ynew + dy/2.

END
