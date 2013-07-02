;General function to describe the density profile:
FUNCTION nu, z, zpars, pars, quadratic, monotonic
  
;Fully non-parametric monotonic *decreasing* function [c.f. Kz function].
;Note, here we explicitly avoid rnuz_z(0) = 0 since this would correspond
;to a constraint rnu_z(zmax) = 0 which is only true if zmax = infty. 
  
;Mirror prior
parsu = abs(pars)

IF (monotonic EQ 'yes') THEN BEGIN 
   rnuz_z = dblarr(n_elements(zpars))
   rnuz_z(0) = parsu(0) 
   FOR i=1L,n_elements(rnuz_z)-1 DO $
      rnuz_z(i) = rnuz_z(i-1) + parsu(i)
   fun = reverse(rnuz_z)
ENDIF ELSE BEGIN 
   ;Alternative here --> don't assume monotonic!
   fun = parsu
ENDELSE

;Normalise: 
IF (quadratic EQ 'no') THEN BEGIN 
   ;Exact linear interpolation integral: 
   norm_nu = double(0)
   FOR i=0L,n_elements(zpars)-2 DO BEGIN
      zl = zpars(i)
      zr = zpars(i+1)
      zz = zr-zl
      b = fun(i)
      a = fun(i+1)
      norm_nu = norm_nu + $
                (a-b)/2./zz*(zr^2.-zl^2.) + $
                b*zr - a*zl
   ENDFOR
ENDIF ELSE BEGIN 
   qtest = 'no'
   IF (qtest EQ 'yes') THEN BEGIN 
      !P.multi = 0
      plot, zpars, pars
      test = interpol(pars,zpars,z,/quadratic)
      oplot,z,test,psym=3,color=2
   ENDIF 
   
   ;Exact quadratic interpolation:
   norm_nu = double(0)
   FOR i=1L,n_elements(zpars)-2 DO BEGIN
      z0 = zpars(i-1)
      z1 = zpars(i)
      z2 = zpars(i+1)
      f0 = fun(i-1)
      f1 = fun(i)
      f2 = fun(i+1)
      
      h = z2-z1
      a = f0
      b = (f1 - a) / h 
      c = (f2 - 2.*b*h - a) / 2. / h^2. 
      
      z1d = z1-z0
      z2d = z1^2. - z0^2.
      z3d = z1^3. - z0^3.
      
      IF (i EQ n_elements(zpars)-2) THEN BEGIN 
         ;Last bin integrate from z0 --> z2: 
         z1d = z2-z0
         z2d = z2^2. - z0^2.
         z3d = z2^3. - z0^3.
      ENDIF 
      
      intquad = (a - b*z0 + c*z0*z1)*z1d + $
                (b/2. - c/2.*(z0+z1))*z2d + $
                c/3.*z3d
      norm_nu = norm_nu + intquad
      
      IF (qtest EQ 'yes') THEN BEGIN
         print, i, h, intquad            
         testy = a + b*(z-z0) + c*(z-z0)*(z-z1)
         IF (i NE n_elements(zpars)-3) THEN BEGIN 
            zcut =  z(where(z GT z0 AND z LT z1))
            tcut = testy(where(z GT z0 AND z LT z1))
         ENDIF ELSE BEGIN 
            zcut =  z(where(z GT z0 AND z LT z2))
            tcut = testy(where(z GT z0 AND z LT z2))
         ENDELSE             
         oplot,zcut,tcut,psym=3,color=4
      ENDIF 
   ENDFOR
   IF (qtest EQ 'yes') THEN BEGIN 
      print, int_tabulated(z,test,/sort), norm_nu 
      stop
   ENDIF 
ENDELSE 
fun = fun / norm_nu

;Interpolate to input z:
IF (quadratic EQ 'no') THEN $
   f = interpol(fun,zpars,z) $
ELSE f = interpol(fun,zpars,z,/quadratic)

;Check for negative density: 
small = min(f(where(f GT 0)))
FOR jj=0L,n_elements(f)-1 DO $
   IF (f(jj) LT 0) THEN f(jj) = small

RETURN, f

END
