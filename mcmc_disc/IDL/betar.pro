;General function to describe the density profile
FUNCTION betar, r, rp_beta, betapars
   

beta_1=dblarr(n_elements(rp_beta)+1)
beta_1(0)=0
beta_1(1:n_elements(rp_beta))=betapars
rpbetanew=dblarr(n_elements(rp_beta)+1)
rpbetanew(0)=0.
rpbetanew(1:n_elements(rp_beta))=rp_beta
beta_r = interpol(beta_1,rpbetanew,r)

RETURN, beta_r

END
