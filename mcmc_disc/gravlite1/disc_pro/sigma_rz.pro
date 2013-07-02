;General function to describe the tilt profile:
FUNCTION sigma_rz, z, zpars, tpars
  
;Mirror prior
tparsu = abs(tpars)
 
;dz = zpars(2)-zpars(1)
;sig_Rz = dblarr(n_elements(zpars))
;sig_Rz(0) = tparsu(0) * dz / 2.
;FOR i=1L,n_elements(sig_Rz)-1 DO $
;   sig_Rz(i) = sig_Rz(i-1) + tparsu(i) * dz
;f = interpol(sig_Rz,zpars,z)

;Alternative here --> don't assume monotonic!  
f = interpol(tparsu,zpars,z)

RETURN, f

END
