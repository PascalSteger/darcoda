;General function to calculate the Kz force law:
FUNCTION Kz, z_in, zpars, kzpars, blow, quadratic
  
;Non-parametric Kz function here. Use cumulative integral to
;ensure monotinicity in an efficient manner. Note here we
;use kz_z(0) = kzpars(0) * dz so that it can be zero, or
;otherwise just small in the case of large bins. This latter
;avoids the interpolated kz_out(0) going negative for coarse 
;binning. 

;Mirror prior + baryon minimum prior [blow]:
kzparsu = abs(kzpars)

;Assume here interpolation between rho "grid points", 
;[linear or quadratic]. The grid points are stored 
;in kzpars and give the *differential* increase in 
;rho(z) towards small z [monotonic rho-prior]: 
denarr = dblarr(n_elements(kzparsu))
denarr(0) = kzparsu(0)
FOR i=1L,n_elements(kzparsu)-1 DO $
   denarr(i) = denarr(i-1) + kzparsu(i)
denarr = reverse(denarr)

;Solve interpolated integral for Kz: 
IF (quadratic EQ 'no') THEN BEGIN 
   ;Linear interpolation here: 
   kz_z = dblarr(n_elements(zpars))
   FOR i=1L,n_elements(zpars)-1 DO BEGIN
      zl = zpars(i-1)
      zr = zpars(i)
      zz = zr-zl
      b = denarr(i-1)
      a = denarr(i)
      kz_z(i) = kz_z(i-1) + (a-b)/2./zz*(zr^2.-zl^2.) + $
                b*zr - a*zl
   ENDFOR
ENDIF ELSE BEGIN 
   ;Quadratic interpolation here: 
   kz_z = dblarr(n_elements(zpars))
   FOR i=1L,n_elements(zpars)-2 DO BEGIN
      z0 = zpars(i-1)
      z1 = zpars(i)
      z2 = zpars(i+1)
      f0 = denarr(i-1)
      f1 = denarr(i)
      f2 = denarr(i+1)

      h = z2-z1
      a = f0
      b = (f1 - a) / h
      c = (f2 - 2.*b*h - a) / 2. / h^2.

      z1d = z1-z0
      z2d = z1^2. - z0^2.
      z3d = z1^3. - z0^3.
   
      intbit = (a - b*z0 + c*z0*z1)*z1d + $
               (b/2. - c/2.*(z0+z1))*z2d + $
               c/3.*z3d

      kz_z(i) = kz_z(i-1) + intbit
      IF (i EQ n_elements(zpars)-2) THEN BEGIN
         ;Deal with very last bin:
         z1d = z2-z1
         z2d = z2^2. - z1^2.
         z3d = z2^3. - z1^3.

         intbit = (a - b*z0 + c*z0*z1)*z1d + $
                  (b/2. - c/2.*(z0+z1))*z2d + $
                  c/3.*z3d
         kz_z(i+1) = kz_z(i) + intbit
      ENDIF
   ENDFOR
ENDELSE 

qtest = 'no'
IF (qtest EQ 'yes') THEN BEGIN 
   !P.multi = 0
   pnts = 5000
   zmin = min(zpars)
   zmax = max(zpars)
   z = dindgen(pnts)*(zmax-zmin)/double(pnts-1) + zmin
   IF (quadratic EQ 'no') THEN $
      denarr_z = interpol(denarr,zpars,z) $
   ELSE denarr_z = interpol(denarr,zpars,z,/quadratic)
   
   kz_th = dblarr(pnts)
   FOR i=1L,pnts-1 DO $
      kz_th(i) = int_tabulated(z(0:i),denarr_z(0:i))
   
   IF (quadratic EQ 'yes') THEN $
      test = interpol(kz_z,zpars,z,/quadratic) $
   ELSE test = interpol(kz_z,zpars,z)

   testsimp = dblarr(n_elements(zpars))
   delta_z = zpars(2) - zpars(1) 
   FOR i=1L,n_elements(zpars)-1 DO $
      testsimp(i) = testsimp(i-1) + denarr(i) * delta_z
   
   plot,z,kz_th
   oplot,zpars,kz_z,color=2
   oplot,z,test,color=4
   oplot,zpars,testsimp,color=3,linestyle=1
   stop
ENDIF 
kz_z = kz_z + blow

;Then interpolate back to the input array:
IF (quadratic EQ 'no') THEN $
   kz_out = interpol(kz_z,zpars,z_in) $
ELSE kz_out = interpol(kz_z,zpars,z_in,/quadratic)

;Error checking. Sometimes when kz_z(0) << kz_z(1), 
;the interpolant can go negative. This just means that 
;we have resolved what should be kz_z(0) = 0. A simple
;fix should suffice: 
FOR jj=0L,n_elements(kz_out)-1 DO $
   IF (kz_out(jj) LT 0) THEN kz_out(jj) = 0.

RETURN, -kz_out

END
