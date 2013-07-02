;============================ SPHERICAL CLEAN ==========================
;A non-parametric method to determine the total enclosed mass of a dwarf spheroidal
;Assume form of nu,M,beta up to twice the maximal data radius

PRO spherical_clean1
  
;Set up the common block:
;setcommon
;COMMON constants

testplot='no'
populations=1
whichpop='stars'

;Set screen output for testing:
IF (testplot eq 'yes') then begin
   devtype='X'
   set_plot,devtype

   ;Set truetype fonts:
   !P.FONT=-1
   !P.CHARSIZE = 1.5
   !P.THICK = 4
   !X.STYLE = 1
   !Y.STYLE = 1

   IF devtype EQ 'X' THEN BEGIN 
      window,xsize=999,ysize=999
      device,decomposed=0
   ENDIF
   IF (populations eq 2) THEN !P.multi=[0,2,3] ELSE !P.multi=[0,2,2]
endif

;Set up my Pallette to contain the three colours I want.
mymakecolour

;Units: 
G1 = 1. ;units are [Mtot],[Rs]

;Set which to do: 
dowhich = 'sim'

;Gradient prior [-1 = no prior]:
gprior = -1
;gprior = 120.


;Central pixel prior [gprior for first point in M_r; default = 0.]: 
cprior = 0.
;cprior = -1

;Baryon minimum surfden prior:
bprior = 'no'
bfifty = 'no'

;Rising sigma_z(z) prior:
sprior = 'no'

;Mirror prior:
mirror = 'yes'

;Logarithmic prior: 
logprior = 'yes'

;Mass prior
mprior = -1

;Initial parameters [-1 = set to default]:
nupars1 = -1
if (populations eq 2) then nupars2 = -1 else nupars2 = 0
Mpars = -1
betapars1 = -1
if (populations eq 2) then betapars2 = -1 else betapars2 = 0

;Default use 10% step size: 
nuparstep1 = -1
if (populations eq 2) then nuparstep2 = -1 else nuparstep2 = 0
Mparstep = -1
betaparstep1 = -1
if (populations eq 2) then betaparstep2 = -1 else betaparstep2 = 0

;Betaprior
betap='yes'
betaprior=0.

;Default low/high-r range = min/max of data: 
rpmin = -1 
rpmax = -1 

;Min/max M for MCMC search [only affects logprior].
;If positive assume constant; if negative take fraction
;of local baryonic value for that bin: 
Mmin = 0.
Mmax = 1d10

;Lbprior (for last bin): 
lbprior = 'no'
lbtol = 0.33

dir="/scratch/psteger/software/mcmc"

;Set number of terms for enclosedmass+tracer+anisotropy models:   
IF (dowhich EQ 'sim') THEN BEGIN
   nt_M = 25
   nt_nu = 25
   nt_beta = 25
   Mno = '25'
   nuno = '25'
   betano = '25'
   whichsim = '1'
   basename = dir+'/sphericalmethod_output/dual_unit_hern_'+whichsim+'/M'+Mno+'_nu'+nuno+'_beta'+betano+'/met'
   massfile = dir+'/sphericalmethod_data/enclosedmass/dual_unit_hern_'+whichsim+'_tot_enclosedmass.txt'
   if (populations eq 2) then begin
      nufile1  = dir+'/sphericalmethod_data/densityfalloff/dual_unit_hern_'+whichsim+'_stars_falloff.txt'
      sigfile1 = dir+'/sphericalmethod_data/velocitydispersion/dual_unit_hern_'+whichsim+'_stars_veldisp.txt'
      nufile2  = dir+'/sphericalmethod_data/densityfalloff/dual_unit_hern_'+whichsim+'_dm_falloff.txt'
      sigfile2 = dir+'/sphericalmethod_data/velocitydispersion/dual_unit_hern_'+whichsim+'_dm_veldisp.txt'
   endif
   IF (populations eq 1) then begin
      nufile1  = dir+'/sphericalmethod_data/densityfalloff/dual_unit_hern_'+whichsim+'_'+whichpop+'_falloff.txt'
      sigfile1 = dir+'/sphericalmethod_data/velocitydispersion/dual_unit_hern_'+whichsim+'_'+whichpop+'_veldisp.txt'
   ENDIF
   dosim = 'yes'
ENDIF

;Iterations:
niter = LONG(900000)


;*************************************************************
;DATA

;Read in the data:
readcol,nufile1,r_dat,nu_dat1,nu_dat_err1
readcol,sigfile1,r_sigdat,sig_dat1,sig_dat_err1
IF (populations eq 2) then readcol,nufile2,r_dat,nu_dat2,nu_dat_err2
IF (populations eq 2) then readcol,sigfile2,r_sigdat,sig_dat2,sig_dat_err2
readcol,massfile,r_Mdat,M_dat,M_dat_err

IF (dosim EQ 'yes') THEN BEGIN 
   ;Fix error on central nu(z) point:
   nu_dat_err1(0) = nu_dat_err1(1)*1.5
   if (populations eq 2) then nu_dat_err2(0) = nu_dat_err2(1)*1.5
ENDIF

;populations
if (populations eq 2) then basename = basename + '_popboth'
if (populations eq 1) then basename = basename + '_pop'+whichpop+''

;Gprior: 
IF (gprior GT 0) THEN basename = basename + '_gprior'
IF (gprior LT 1) THEN gprior = 1d30

;Cprior:
IF (cprior GE 0) THEN basename = basename + '_cprior'
IF (cprior LT 0) THEN cprior = 1d30

;Mirror, log, lb:
IF (mirror EQ 'yes') THEN basename = basename + '_mirr'
IF (logprior EQ 'yes') THEN basename = basename + '_log'
IF (lbprior EQ 'yes') THEN basename = basename + '_lb'

;betaprior
IF (betap eq 'yes') THEN basename = basename + '_beta' 

if (mprior lt 0) then basename = basename + '_mslope' else basename = basename + '_mconst'

;Bprior:
IF (bprior EQ 'yes') THEN BEGIN 
   IF (bfifty EQ 'no') THEN basename = basename + '_bprior' $
   ELSE basename = basename + '_b50'

   IF (baryonmodel EQ 'sim') THEN BEGIN
      readcol,massfile,rvis,massusevis,massuseviserr
   ENDIF
   
   massvismin = massusevis - massuseviserr
ENDIF
  
;S-prior:
IF (sprior EQ 'yes') THEN basename = basename + '_sprior'
   
;Binning in z:
IF (rpmin LT 0) THEN rpmin = min(r_dat)
IF (rpmax LT 0) THEN rpmax = max(r_dat)

rp_M = dindgen(nt_M) * (rpmax - rpmin) / double(nt_M-1) + rpmin
rp_nu = dindgen(nt_nu) * (rpmax - rpmin) / double(nt_nu-1) + rpmin
rp_beta = dindgen(nt_beta) * (rpmax - rpmin) / double(nt_beta-1) + rpmin

;Convert gprior into appropriate units: 
IF (gprior GT 0) THEN $
   gpriorconv = gprior $ ;unit [Msun/kpc^3]
ELSE gpriorconv = 1d30
IF (cprior GE 0) THEN $
   cpriorconv = cprior $ ;unit [Msun/kpc^3]
ELSE cpriorconv = 1d30

;Set up bprior:
IF (bprior EQ 'yes') THEN BEGIN
   bmass = interpol(massvismin,rvis,rp_M) 	;interpolate on nodes zp_M, unit [Msun]
   dr = rp_M(2) - rp_M(1) 			;binlength
   blow = dblarr(n_elements(rp_M))
   FOR i=1L,n_elements(rp_M)-1 DO $
      blow(i) = (bmass(i) - bmass(i-1))/dr 		; slope??
ENDIF ELSE blow = dblarr(n_elements(rp_M))

;Set up kzmin/max arrays:
IF (Mmin LT 0) THEN $
   Mminarr = dblarr(nt_M) + abs(Mmin)*blow $
ELSE Mminarr = dblarr(nt_M) + Mmin
IF (Mmax LT 0) THEN $
   Mmaxarr = dblarr(nt_M) + abs(Mmax)*blow $
ELSE Mmaxarr = dblarr(nt_M) + Mmax

;Default Initial seed guess parameters [assume flat]:
IF (Mpars(0) LT 0) THEN Mpars = dblarr(nt_M) + .3  ;unit [Mtot]
IF (nupars1(0) LT 0) THEN nupars1 = dblarr(nt_nu) + .3
IF (betapars1(0) LT 0) THEN betapars1 =  dblarr(nt_beta)
IF (nupars2(0) LT 0) THEN nupars2 = dblarr(nt_nu) + .3
IF (betapars2(0) LT 0) THEN betapars2 =  dblarr(nt_beta)
IF (mprior(0) LT 0) THEN Mslopepars = 0.5 else Mslopepars = mprior
;IF (norm1 LT 0) THEN norm1 = 20.^2.
;IF (norm2 LT 0) THEN norm2 = 20.^2.

;*************************************************************
;Metropolis parameters here:
if (populations eq 2) then pars = [nupars1,nupars2,Mpars,betapars1,betapars2,Mslopepars] else pars = [nupars1,Mpars,betapars1,Mslopepars]
npars = n_elements(pars)        ;total number of parameter
metarray = dblarr(npars+1,niter)
IF (nuparstep1(0) LT 0) THEN nuparstep1 = nupars1 / 10.
IF (nuparstep2(0) LT 0) THEN nuparstep2 = nupars2 / 10.
IF (Mparstep(0) LT 0) THEN Mparstep = dblarr(nt_M) + 0.1 ;unit [Mtot]
IF (betaparstep1(0) LT 0) THEN betaparstep1 = 0.02
IF (betaparstep2(0) LT 0) THEN betaparstep2 = 0.02
IF (mprior(0) LT 0) THEN Mslopeparstep = 0.02 else Mslopeparstep = 0.
;IF (normstep1 LT 0) THEN normstep1 = norm1/10. 
;IF (normstep1 LT 0) THEN normstep2 = norm2/10. 

if (populations eq 2) then parstep = [nuparstep1,nuparstep2,Mparstep,betaparstep1,betaparstep2,Mslopeparstep] else parstep = [nuparstep1,Mparstep,betaparstep1,Mslopeparstep]

;Logarithmic prior: 
IF (logprior EQ 'yes') THEN BEGIN
   lpars=dblarr(npars)
   lpars(0:populations*nt_nu+nt_M-1) = alog10(pars(0:populations*nt_nu+nt_M-1)) ;log sampling for nu and M
   lpars(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta) = pars(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta) ;linear sampling for beta and Mslope
   lparstep = lpars 
   lparstep(0:nt_nu-1) = 0.05 ;nu1
   if (populations eq 2) then lparstep(nt_nu:2*nt_nu-1) = 0.05 ;nu2
   lparstep(populations*nt_nu:populations*nt_nu+nt_M-1) = 0.05 ;M 
   lparstep(populations*nt_nu+nt_M:populations*nt_nu+nt_M+nt_beta-1) = 0.02 ;beta1
   if (populations eq 2) then lparstep(2*nt_nu+nt_M+nt_beta:2*nt_nu+nt_M+2*nt_beta-1) = 0.02 ;beta2
   lparstep(populations*nt_nu+nt_M+populations*nt_beta) = Mslopeparstep ;Mslope
;   lparstep(2*nt_nu+nt_M+2*nt_beta) = 1.25 ;norm1
;   lparstep(2*nt_nu+nt_M+2*nt_beta+1) = 1.25 ;norm2
ENDIF

;Output data files:
outfilemet = basename+'.dat'
outfiledat = basename+'.txt'
outfileprofM = basename+'.profM'
outfileprofnu1 = basename+'.profnu1'
if (populations eq 2) then outfileprofnu2 = basename+'.profnu2'
outfileprofbeta1 = basename+'.profbeta1'
if (populations eq 2) then outfileprofbeta2 = basename+'.profbeta2'

;Write the key data parameters to a corresponding file:
openw,12,outfiledat
printf,12,'Number of terms [M,nu,beta]  & iterations:'
printf,12,nt_nu,nt_M,nt_beta,niter
printf,12,'Run parameters [gprior, cprior, bprior]:'
printf,12,gprior,cprior,bprior
close,12
arraydump,outfileprofM,[[rp_M]],1
arraydump,outfileprofnu1,[[rp_nu]],1
if (populations eq 2) then arraydump,outfileprofnu2,[[rp_nu]],1
arraydump,outfileprofbeta1,[[rp_beta]],1
if (populations eq 2) then arraydump,outfileprofbeta2,[[rp_beta]],1


;*************************************************************
;METROPOLIS

chisq = 1d300                   ;Initial chisq [large]
n = 0L                          ;counter
initphase = 'over'             ;Initialisation phase flag
rejcount = 0.                   ;Rejection count
acccount = 0.                   ;Acceptance count
accrejtollow = 0.22              ;Acceptance/rejection rate
accrejtolhigh = 0.28            ;
if (populations eq 1) then chisqtol = 10. else chisqtol = 20.     ;Parameters to end initphase 
endcount = 200                  ;
endgame = 'no'                  ;Ending flag
WHILE (n LT niter-1) DO BEGIN
   ;Store the old parameters:
   IF (initphase EQ 'over') THEN BEGIN
      metarray(0:npars-1,n) = pars
      metarray(npars,n) = chisq
   ENDIF 

AGAIN:
   ;Wiggle the parameters -> new parameters:
   IF (logprior EQ 'no') THEN BEGIN 
      ranarr = 2.0 * (0.5 - randomu(seed,/DOUBLE,npars))
      parst = pars + ranarr * parstep
   ENDIF ELSE BEGIN 
      ranarr = 2.0 * (0.5 - randomu(seed,/DOUBLE,npars))
      lparst = lpars + ranarr * lparstep 		;take one step in logspace
      parst = dblarr(npars)
      parst(0:populations*nt_nu+nt_M-1) = double(10)^lparst(0:populations*nt_nu+nt_M-1) ;log sampling for nu and M
      IF ( betap eq 'no' ) THEN parst(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta-1) = lparst(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta-1) $  ;linear sampling for beta
      ELSE parst(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta-1) = betaprior
      parst(populations*nt_nu+nt_M+populations*nt_beta) = lparst(populations*nt_nu+nt_M+populations*nt_beta) ;linear sampling for Mslope

      ;Mmin/max range prior:
      Mpars = parst(populations*nt_nu:populations*nt_nu+nt_M-1) 
      FOR jj=0L,nt_M-1 DO BEGIN
         IF (Mpars(jj) LT Mminarr(jj)) THEN GOTO,AGAIN
         IF (Mpars(jj) GT Mmaxarr(jj)) THEN GOTO,AGAIN
      ENDFOR
   ENDELSE

   ;Lbprior 
   IF (lbprior EQ 'yes') THEN BEGIN 
      totmlastb = total(Mpars(0:nt_M-2) + blow(0:nt_M-2))
      lastb = Mpars(nt_M-1) + blow(nt_M-1)
      IF (lastb / totmlastb GT lbtol) THEN GOTO,AGAIN
   ENDIF

   ;Ensure positivity --> monotinicity constraint:
   IF (mirror EQ 'no') THEN BEGIN
      FOR jj=0L,npars-1 DO $
         IF (parst(jj) LT 0.) THEN GOTO,AGAIN
   ENDIF
   ;Extract parameters:
   nupars1 = parst(0:nt_nu-1)
   if (populations eq 2) then nupars2 = parst(nt_nu:2*nt_nu-1)
   Mpars = parst(populations*nt_nu:populations*nt_nu+nt_M-1)
   betapars1 = parst(populations*nt_nu+nt_M:populations*nt_nu+nt_M+nt_beta-1)
   if (populations eq 2) then betapars2 = parst(2*nt_nu+nt_M+nt_beta:2*nt_nu+nt_M+2*nt_beta-1)
   Mslopepars = parst(populations*nt_nu+nt_M+populations*nt_beta)
;   norm1 = parst(2*nt_nu+nt_M+2*nt_beta)
;   norm2 = parst(2*nt_nu+nt_M+2*nt_beta+1)
   
   ;Apply "cprior" and "gprior" on M function:
   ;IF (cpriorconv GT 0) THEN BEGIN
   ;   IF (abs(Mpars(0)) GT cpriorconv) THEN GOTO,AGAIN
   ;ENDIF ELSE BEGIN 
   ;   Mpars(0) = Mminarr(0)
   ;   parst(2*nt_nu) = Mminarr(0)
   ;ENDELSE

   IF (abs(Mslopepars) GT 2.) THEN GOTO,AGAIN
 
   FOR jj=0L,nt_M-1 DO $
      IF (abs(Mpars(jj)) GT gpriorconv) THEN GOTO,AGAIN

   ;Calculate theoretical values: 
   nu1_r = nu(r_dat,rp_nu,nupars1)
   if (populations eq 2) then nu2_r = nu(r_dat,rp_nu,nupars2)
   sig1_r = sigma_r1(r_sigdat,rp_nu,rp_M,rp_beta,Mpars,nupars1,betapars1,Mslopepars)
   if (populations eq 2) then sig2_r = sigma_r1(r_sigdat,rp_nu,rp_M,rp_beta,Mpars,nupars2,betapars2,Mslopepars)

   ;Reject models with NaN sig_z: 
   IF (logprior EQ 'yes') THEN BEGIN 
    ;  print,'Nan condition'
      FOR jj=0L,n_elements(sig1_r)-1 DO $
         IF (sig1_r(jj) NE sig1_r(jj)) THEN GOTO, AGAIN
      if (populations eq 2) then begin
         FOR jj=0L,n_elements(sig2_r)-1 DO $
            IF (sig2_r(jj) NE sig2_r(jj)) THEN GOTO, AGAIN
      endif    
   ENDIF ELSE BEGIN 
      ;Another problem with linear sampling (best avoided). 
      ;If I do the "right" thing here and reject NaN models, 
      ;the MCMC gets stuck. So we set sig_z = 0. where it
      ;is NaN and assume that this will be penalised by the 
      ;data. This assumption appears to be very good, 
      ;but still it's not ideal ): 
      FOR jj=0L,n_elements(sig1_r)-1 DO $
         IF (sig1_r(jj) NE sig1_r(jj)) THEN sig1_r(jj) = 0.
      if (populations eq 2) then begin
         FOR jj=0L,n_elements(sig2_r)-1 DO $
            IF (sig2_r(jj) NE sig2_r(jj)) THEN sig2_r(jj) = 0.
      endif    
   ENDELSE
  ; print,'before sprior'
   ;S-prior: ensure sigma_z(z) rises: 
   IF (sprior EQ 'yes') THEN BEGIN 
      FOR jj=1L,n_elements(sig_z)-1 DO $
         IF (sig_z(jj) LT sig_z(jj-1)) THEN BEGIN 
         GOTO,AGAIN
         ENDIF
   ENDIF
   ;print,'before chisq'
   ;Calculate chi-squared
   if (populations eq 2) then chisqt = total((nu1_r - nu_dat1)^2./nu_dat_err1^2.) + $
            total((nu2_r - nu_dat2)^2./nu_dat_err2^2.) + $
            total((sig1_r - sig_dat1)^2./sig_dat_err1^2.) + $
            total((sig2_r - sig_dat2)^2./sig_dat_err2^2.)
   if (populations eq 1) then chisqt = total((nu1_r - nu_dat1)^2./nu_dat_err1^2.) + $
            total((sig1_r - sig_dat1)^2./sig_dat_err1^2.)
;   chisqt = total((nu1_r - nu_dat1)^2.) + $
;            total((nu2_r - nu_dat2)^2.) + $
;            total((sig1_r - sig_dat1)^2.) + $
;            total((sig2_r - sig_dat2)^2.)


   ;Calculate the f-function
   fnewoverf = EXP(chisq/2.0-chisqt/2.0)
    
   ;Accept the new f-function?
   ran = randomu(seed,/DOUBLE)
   
   IF (ran LT fnewoverf) THEN BEGIN
      acccount = acccount + 1.
      pars = parst
      chisq = chisqt

      IF (logprior EQ 'yes') THEN lpars = lparst

      IF (testplot EQ 'yes') THEN BEGIN
         ;Calculate profiles: 
         mprioru = abs(Mslopepars)
         rplmin = 0;
         rplmax = max(rp_M)
         rpnts = 500
         rpl = dindgen(rpnts) * (rplmax-rplmin)/double(rpnts-1) + rplmin
         dr = (rplmax-rplmin)/double(rpnts-1)
         nu1 = nu(rpl,rp_nu,nupars1)
         sig1 = sigma_r1(rpl,rp_nu,rp_M,rp_beta,Mpars,nupars1,betapars1,Mslopepars)
         M = Mr(rpl,rp_M,Mpars)
         M_outer=M(rpnts-1)+(dindgen(rpnts-1)+1)*mprioru*M(rpnts-1)/(rplmax-rplmin)*dr
         nu1_outer=nu1(rpnts-1)-(dindgen(rpnts-1)+1)*nu1(rpnts-1)/(rplmax-rplmin)*dr
         M_tot=dblarr(2*rpnts-1)
         M_tot(0:rpnts-1)=M
         M_tot(rpnts:2*rpnts-2)=M_outer
         r_tot=dblarr(2*rpnts-1)
         r_tot(0:rpnts-1)=rpl
         r_tot(rpnts:2*rpnts-2)=rplmax+rpl(1:rpnts-1)
         beta1 = betar(rpl,rp_beta,betapars1)
         nu1_tot=dblarr(2*rpnts-1)
         nu1_tot(0:rpnts-1)=nu1
         nu1_tot(rpnts:2*rpnts-2)=nu1_outer

         ;Do the plots: 
         plot,r_dat,nu_dat1,xrange=[0,rplmax],yrange=[0.000001,1.5],/ylog,$
              title='!6',xtitle='r(Rs)',ytitle='nu1 [unity]'
         oploterror,r_dat,nu_dat1,nu_dat_err1
         oplot,r_tot,nu1_tot,color=2

         plot,r_dat,sig_dat1,xrange=[0,rplmax],yrange=[0,1.],/NODATA,$
              title='!6',xtitle='r(Rs)',ytitle='sigma_r1 []'
         oploterror,r_sigdat,sig_dat1,sig_dat_err1
         oplot,rpl,sig1,color=2

         if (populations eq 2) then begin
            sig2 = sigma_r1(rpl,rp_nu,rp_M,rp_beta,Mpars,nupars2,betapars2,Mslopepars)
            nu2 = nu(rpl,rp_nu,nupars2)
            beta2 = betar(rpl,rp_beta,betapars2)

            plot,r_dat,nu_dat2,xrange=[0,rplmax],yrange=[0.000001,1.5],/ylog,$
                 title='!6',xtitle='r(Rs)',ytitle='nu2 [unity]'
            oploterror,r_dat,nu_dat2,nu_dat_err2
            oplot,rpl,nu2,color=2

            plot,r_dat,sig_dat2,xrange=[0,rplmax],yrange=[0,1.],/NODATA,$
                 title='!6',xtitle='r(Rs)',ytitle='sigma_r2 []'
            oploterror,r_sigdat,sig_dat2,sig_dat_err2
            oplot,rpl,sig2,color=2
         endif
         plot,r_Mdat,M_dat,xrange=[0,2*rplmax],yrange=[0,2.],/NODATA,$
              title='!6',xtitle='r(Rs)',ytitle='M [Mtot]'
         oploterror,r_Mdat,M_dat,M_dat_err
         oplot,r_tot,M_tot,color=2

         plot,rpl,beta1,xrange=[rplmin,rplmax],yrange=[-1.5,1.5],/NODATA,$
              title='!6',xtitle='r(Rs)',ytitle='M [Mtot]'
         if (populations eq 2) then oplot,rpl,beta2,color=3
         oplot,rpl,beta1,color=2

      ENDIF
      
      ;Output for testing: 
  ;    IF (logprior EQ 'no') THEN stepout = parstep(nt_kz-1) $
  ;    ELSE stepout = lparstep(nt_kz-1)
      print,n,chisq,lparstep(0),acccount/rejcount

      ;Decide whether to end initphase:
      IF (endgame EQ 'yes') THEN BEGIN
         endcount = endcount - 1
         IF (endcount EQ 0) THEN BEGIN
            initphase = 'over'
            print,''
            print,'********** Initialisation phase over **********'
            print,'Step size chosen:'
            IF (logprior EQ 'no' ) THEN BEGIN
               print,parstep
               print,pars
            ENDIF ELSE BEGIN
               print,lparstep
               print,lpars
               print,pars
            ENDELSE
            print,''
         ENDIF
      ENDIF
   ENDIF ELSE BEGIN 
      ;Rejected
      rejcount = rejcount + 1.
   ENDELSE 

   ;Adapt stepsize during initialisation phase: 
   IF (initphase EQ 'start') THEN BEGIN
      IF (acccount GT 0) AND (rejcount GT 0) THEN BEGIN 
         IF (acccount/rejcount LT accrejtollow) THEN BEGIN 
            parstep = parstep / 1.001
            IF (logprior EQ 'yes') THEN lparstep = lparstep / 1.01
         ENDIF 
         IF (acccount/rejcount GT accrejtolhigh) THEN BEGIN
            parstep = parstep * 1.001
            IF (logprior EQ 'yes') THEN lparstep = lparstep * 1.01
         ENDIF
      ENDIF
      IF (chisq LT chisqtol) THEN endgame = 'yes'
   ENDIF

   ;Update the counter:
   IF (initphase EQ 'over') THEN n = n + 1
ENDWHILE

;Write the data to a file:
openw,1,outfilemet
writeu,1,LONG(npars),LONG(niter)
writeu,1,metarray
close,1

print,'Finished!'

;close the output file, set the plotting back to the screen and exit.
if devtype EQ 'PS' then device,/close	
set_plot, 'X'
!P.FONT=-1
!P.multi=0

END
