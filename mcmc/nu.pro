;General function to describe the density profile
FUNCTION nu, z, zpars, pars
  
;Fully non-parametric monotonic *decreasing* function [c.f. Kz function].
;Note, here we explicitly avoid rnuz_z(0) = 0 since this would correspond
;to a constraint rnu_z(zmax) = 0 which is only true if zmax = infty. 
  
;Mirror prior
parsu = abs(pars)

dz = zpars(2)-zpars(1)
rnuz_z = dblarr(n_elements(zpars))
rnuz_z(0) = parsu(0) * dz
FOR i=1L,n_elements(rnuz_z)-1 DO $
   rnuz_z(i) = rnuz_z(i-1) + parsu(i) * dz
f = interpol(reverse(rnuz_z),zpars,z)

RETURN, f

END
