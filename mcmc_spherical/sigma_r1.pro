;General function to calculate Sigma_r:
FUNCTION sigma_r1, r_in, rp_nu, rp_M, rp_beta, Mpars, nupars, betapars, mprior

;print,betapars
;Mirror prior
mprioru=abs(mprior)

G1 = 1.

;Set up r array for integration:
pnts = 100
rmin = 0.
rmax = max(r_in)
dr = (rmax-rmin)/double(pnts-1)
r = dindgen(pnts) * dr + rmin

;Calculate density and M force:
nu_r = nu(r,rp_nu,nupars)
M_r = Mr(r,rp_M,Mpars)
beta_r=betar(r,rp_beta,betapars)
;beta_1=dblarr(n_elements(rp_beta)+1)
;beta_1(0)=0
;beta_1(1:n_elements(rp_beta))=betapars
;rpbetanew=dblarr(n_elements(rp_beta)+1)
;rpbetanew(0)=0.
;rpbetanew(1:n_elements(rp_beta))=rp_beta
;beta_r = interpol(beta_1,rpbetanew,r)

;arrays of length pnts-1
;nu goes down from nu(pnts-1) at rmax to zero at 2*rmax
;M rises from M(pnts-1) at rmax up to 2*M(pnts-1) at 2*rmax)
;beta is assumed to be constant from rmax to 2*rmax
nu_outer=nu_r(pnts-1)-(dindgen(pnts-1)+1)*nu_r(pnts-1)/(rmax-rmin)*dr
M_outer=M_r(pnts-1)+(dindgen(pnts-1)+1)*mprioru*M_r(pnts-1)/(rmax-rmin)*dr
beta_outer=(dindgen(pnts-1)+1)/(dindgen(pnts-1)+1)*beta_r(pnts-1)

;merge arrays
r_tot=dblarr(2*pnts-1)
r_tot(0:pnts-1)=r
r_tot(pnts:2*pnts-2)=rmax+r(1:pnts-1)
nu_tot=dblarr(2*pnts-1)
nu_tot(0:pnts-1)=nu_r
nu_tot(pnts:2*pnts-2)=nu_outer
M_tot=dblarr(2*pnts-1)
M_tot(0:pnts-1)=M_r
M_tot(pnts:2*pnts-2)=M_outer
beta_tot=dblarr(2*pnts-1)
beta_tot(0:pnts-1)=beta_r
beta_tot(pnts:2*pnts-2)=beta_outer

;print,nu_tot
;print,M_tot
;print,beta_tot
;print,r_tot
;Solve Beta integral

;print,beta_tot
beta_r_t = dblarr(2*pnts-1)
FOR i=2,2*pnts-2 DO $
   beta_r_t(i) = int_tabulated(r_tot(1:i),beta_tot(1:i)/r_tot(1:i))

;r_int=reverse(r_tot)
;nu_int=reverse(nu_tot)
;M_int=reverse(M_tot)
;beta_int=reverse(beta_r_t)
;print,'nu'
;print,nu_tot
;print,'M'
;print,M_tot
;print,'beta'
;print,beta_tot
;print,'r'
;print,r_tot
;print,'exp'
;print,exp(2*beta_tot)
;print,'exp2'
;print,exp(-2*beta_r_t)

;Solve cumulative intgral.
;the integrand is assumed to be zero above 2*rmax
sig_r_t = dblarr(pnts)
FOR i=1L,pnts-1 DO $
   sig_r_t(i) = exp(-2*beta_r_t(i)) * int_tabulated(r_tot(i:2*pnts-2),G1*M_tot(i:2*pnts-2)*nu_tot(i:2*pnts-2)*exp(2*beta_tot(i:2*pnts-2))/(r_tot(i:2*pnts-2)^2)) / nu_r(i)
sig_r_t(0) = sig_r_t(1)
;print,sig_r_t
;Interpolate back to input array:
sig_r = interpol(sig_r_t,r,r_in)
;print,sig_r
f = sqrt(sig_r)
;stop
RETURN, f

END
