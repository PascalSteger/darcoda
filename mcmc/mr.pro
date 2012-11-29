;General function to calculate enclosed mass M
FUNCTION Mr, z_in, zpars, Mpars
  
;Non-parametric Kz function here. Use cumulative integral to
;ensure monotinicity in an efficient manner. Note here we
;use kz_z(0) = kzpars(0) * dz so that it can be zero, or
;otherwise just small in the case of large bins. This latter
;avoids the interpolated kz_out(0) going negative for coarse 
;binning. 

;Mirror prior:
Mparsu = abs(Mpars)
;print,Mparsu
;print,zpars
dz = zpars(2)-zpars(1)
M_r = dblarr(n_elements(zpars)+1)
M_r(1) = Mparsu(0) * dz
FOR i=2L,n_elements(M_r)-1 DO $
   M_r(i) = M_r(i-1) + Mparsu(i-1) * dz
M_r(0) = 0.
z=dblarr(n_elements(zpars)+1)
z(0)=0
z(1:n_elements(zpars))=zpars

;Then interpolate back to the input array:
M_out = interpol(M_r,z,z_in)

;Error checking. Sometimes when kz_z(0) << kz_z(1), 
;the interpolant goes negative. This just means that 
;we have resolved what should be kz_z(0) = 0. A simple
;fix should suffice: 
IF (M_out(0) LT 0) THEN M_out(0) = 0.

;This should never happen: 
;IF (min(M_out) LT 0) THEN BEGIN
;   print,'Kz sign error! stopping ... '
;   stop
;ENDIF 

RETURN, M_out

END
