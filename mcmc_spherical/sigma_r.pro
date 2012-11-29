;General function to calculate Sigma_r:
FUNCTION sigma_r, r_in, rp_nu, rp_M, rp_beta, Mpars, nupars, betapars

;print,betapars
;Mirror prior

G1 = 1.

;Set up r array for integration:
pnts = 400
pntssig = 200
rmin = 0.
rmax = max(rp_M)
dr = (rmax-rmin)/double(pnts-1)
r = dindgen(pnts) * dr + rmin

;Calculate density and M force:
nu_r = nu(r,rp_nu,nupars)
M_r = Mr(r,rp_M,Mpars)
beta_r=betar(r,rp_beta,betapars)

;print,nu_r
;print,M_r
;print,beta_r

;print,beta_tot
beta_r_t = dblarr(pnts)
FOR i=2,pnts-1 DO $
   beta_r_t(i) = int_tabulated(r(1:i),beta_r(1:i)/r(1:i))

;print,beta_r_t

;Solve cumulative intgral.
;the integrand is assumed to be zero above 2*rmax
sig_r_t = dblarr(pntssig)
FOR i=1L,pntssig-1 DO $
   sig_r_t(i) = exp(-2*beta_r_t(i)) * int_tabulated(r(i:pnts-1),G1*M_r(i:pnts-1)*nu_r(i:pnts-1)*exp(2*beta_r_t(i:pnts-1))/(r(i:pnts-1)^2)) / nu_r(i)
sig_r_t(0) = sig_r_t(1)
;print,sig_r_t
;Interpolate back to input array:
;print,sig_r_t
sig_r = interpol(sig_r_t,r(0:pntssig-1),r_in)
;print,sig_r
f = sqrt(sig_r)
;print,r(0:pntssig-1)
;print,r_in
;stop
RETURN, f

END
