;General function to calculate Sigma_z:
FUNCTION sigma_z, z_in, zp, kzpars, blow, nupars, norm, tpars, tparsR, quadratic, monotonic 

;Mirror prior:
normu = abs(norm)

;Calculate density and Kz force: 
nu_z = nu(zp,zp,nupars,quadratic,monotonic)
kz_z = kz(zp,zp,kzpars,blow,quadratic)

;Add tilt correction [if required]:
IF (tparsR(0) GT 0) THEN BEGIN
   Rsun = tparsR(0)
   hr = tparsR(1)
   hsig = tparsR(2)

   tc = sigma_rz(zp,zp,tpars)
   tc = tc * (1.0/Rsun - 1.0/hr - 1.0/hsig)
   
   ;Flag when the tilt becomes significant:
   IF (abs(total(tc))/abs(total(kz_z)) GT 0.1) THEN $
      print,'Tilt > 10%!',abs(total(tc)),abs(total(kz_z))

   kz_z = kz_z - tc
ENDIF

;Do exact integral assuming linear/quad. interpolation: 
IF (quadratic EQ 'no') THEN BEGIN
   ;Linear interpolation here: 
   sigint = dblarr(n_elements(zp))
   FOR i=1L,n_elements(zp)-1 DO BEGIN 
      zl = zp(i-1)
      zr = zp(i)
      zz = zr-zl
      b = nu_z(i-1)
      a = nu_z(i)
      q = kz_z(i-1)
      p = kz_z(i)
      
      intbit = (a-b)*(p-q)/(3.0*zz^2.)*(zr^3.-zl^3.) + $
               ((a-b)/(2.0*zz)*(q - (p-q)*zl/zz)+$
                (p-q)/(2.0*zz)*(b - (a-b)*zl/zz))*(zr^2.-zl^2.)+$
               (b - (a-b)*zl/zz)*(q - (p-q)*zl/zz)*zz
      
      IF (i EQ 0) THEN $
         sigint(0) = intbit $
      ELSE sigint(i) = sigint(i-1) + intbit
   ENDFOR
ENDIF ELSE BEGIN
   sigint = dblarr(n_elements(zp))
   FOR i=1L,n_elements(zp)-2 DO BEGIN
      z0 = zp(i-1)
      z1 = zp(i)
      z2 = zp(i+1)
      f0 = nu_z(i-1)
      f1 = nu_z(i)
      f2 = nu_z(i+1)
      fd0 = kz_z(i-1)
      fd1 = kz_z(i)
      fd2 = kz_z(i+1)      

      h = z2-z1
      a = f0
      b = (f1 - a) / h
      c = (f2 - 2.*b*h - a) / 2. / h^2.
      
      ad = fd0
      bd = (fd1 - ad) / h
      cd = (fd2 - 2.*bd*h - ad) / 2. / h^2.

      AA = a*bd + ad*b
      BB = c*ad + cd*a
      CC = b*cd + bd*c

      z1d = z1-z0
      z2d = z1^2. - z0^2.
      z3d = z1^3. - z0^3.
      z4d = z1^4. - z0^4.
      z5d = z1^5. - z0^5.

      intbit = z1d * (a*ad - z0*(AA-BB*z1) + z0^2.*(b*bd-CC*z1) + $
                      c*cd*z0^2.*z1^2.) + $
               z2d * (0.5*(AA-BB*z1)-z0*(b*bd-CC*z1) - z0/2.*BB + z0^2./2.*CC - $
                      (z0*z1^2. + z0^2.*z1)*c*cd) + $
               z3d * (1.0/3.0*(b*bd-CC*z1)+1.0/3.0*BB - 2.0/3.0*z0*CC + $
                      1.0/3.0*c*cd*(z1^2. + z0^2. + 4.0*z0*z1)) + $
               z4d * (1.0/4.0*CC - c*cd/2.0 * (z1 + z0)) + $
               z5d * c*cd / 5.      

      sigint(i) = sigint(i-1) + intbit
      IF (i EQ n_elements(zp)-2) THEN BEGIN 
         ;Deal with very last bin: 
         z1d = z2-z1
         z2d = z2^2. - z1^2.
         z3d = z2^3. - z1^3.
         z4d = z2^4. - z1^4.
         z5d = z2^5. - z1^5.
         
         intbit = z1d * (a*ad - z0*(AA-BB*z1) + z0^2.*(b*bd-CC*z1) + $
                         c*cd*z0^2.*z1^2.) + $
                  z2d * (0.5*(AA-BB*z1)-z0*(b*bd-CC*z1) - z0/2.*BB + z0^2./2.*CC - $
                         (z0*z1^2. + z0^2.*z1)*c*cd) + $
                  z3d * (1.0/3.0*(b*bd-CC*z1)+1.0/3.0*BB - 2.0/3.0*z0*CC + $
                         1.0/3.0*c*cd*(z1^2. + z0^2. + 4.0*z0*z1)) + $
                  z4d * (1.0/4.0*CC - c*cd/2.0 * (z1 + z0)) + $
                  z5d * c*cd / 5.

         sigint(i+1) = sigint(i) + intbit
      ENDIF 
   ENDFOR
ENDELSE 

qtest = 'no'
IF (qtest EQ 'yes') THEN BEGIN 
   !P.multi = 0
   pnts = 5000
   zmin = min(zp) 
   zmax = max(zp)
   z = dindgen(pnts)*(zmax-zmin)/double(pnts-1) + zmin 
   nu_z = nu(z,zp,nupars,quadratic)
   kz_z = kz(z,zp,kzpars,blow,quadratic)

   sigint_th = dblarr(pnts) 
   FOR i=1L,pnts-1 DO $
      sigint_th(i) = int_tabulated(z(0:i),Kz_z(0:i)*nu_z(0:i)) 

   IF (quadratic EQ 'yes') THEN $
      test = interpol(sigint,zp,z,/quadratic) $
   ELSE test = interpol(sigint,zp,z)

   plot,z,sigint_th
   oplot,zp,sigint,color=2
   oplot,z,test,color=4
   stop
ENDIF

sig_z_t2 = 1.0/nu_z * (sigint + normu)

;Interpolate back to input array:
IF (quadratic EQ 'no') THEN $
   sig_z2 = interpol(sig_z_t2,zp,z_in) $
ELSE sig_z2 = interpol(sig_z_t2,zp,z_in,/quadratic)
f = sqrt(sig_z2)

RETURN, f

END
