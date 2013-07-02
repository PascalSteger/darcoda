;============================ SPHERICAL CLEAN ==========================
;A non-parametric method to determine the total enclosed mass of a dwarf spheroidal
;Implemented for two independent tracer populations

PRO dwarf_los
  
;Set up the common block:
;setcommon
;COMMON constants

testplot='yes'
debug='yes'
populations=2
whichpop='1'

;Set screen output for testing:
if (testplot eq 'yes') then begin
   devtype='X'
   set_plot,devtype

   ;Set truetype fonts:
   !P.FONT=-1
   !P.CHARSIZE = 1.5
   !P.THICK = 4
   !X.STYLE = 1
   !Y.STYLE = 1

   if devtype eq 'X' then begin 
      window,xsize=999,ysize=999
      device,decomposed=0
   endif
   if (populations eq 2) then !P.multi=[0,2,3] else !P.multi=[0,2,2]
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

;sigmaprior
sigmaprior1 = -1
sigmaprior2 = -1

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

;Set number of terms for enclosedmass+tracer+anisotropy models:   

if (dowhich eq 'sim') then begin
   nt_M = 12
   nt_nu = 12
   nt_beta = 12
   dir = '/home/ast/read/dark/dwarf_data/mock/'
   basename = dir+'c1_100_050_100_100_core_c2_010_050_100_100_core_003_6d.mem2'
   massfile = basename+'_totenclosedmass.txt'
   if (populations eq 2) then begin
      nufile1  = basename+'_densityfalloff_1.txt'
      sigfile1 = basename+'_velocitydispersionlos_1.txt'
      nufile2  = basename+'_densityfalloff_2.txt'
      sigfile2 = basename+'_velocitydispersionlos_2.txt'
   endif
   if (populations eq 1) then begin
      nufile1  = basename+'_densityfalloff_'+whichpop+'.txt'
      sigfile1 = basename+'_velocitydispersionlos_'+whichpop+'.txt'
   endif
   dosim = 'yes'
endif

;Iterations:
niter = LONG(900000)

;DATA
readcol,nufile1,r_dat,nu_dat1,nu_dat_err1
readcol,sigfile1,r_sigdat,sig_dat1,sig_dat_err1
if (populations eq 2) then readcol,nufile2,r_dat,nu_dat2,nu_dat_err2
if (populations eq 2) then readcol,sigfile2,r_sigdat,sig_dat2,sig_dat_err2
readcol,massfile,r_Mdat,M_dat,M_dat_err

if (dosim eq 'yes') then begin 
   ;Fix error on central nu(z) point:
   nu_dat_err1(0) = nu_dat_err1(1)*1.5
   if (populations eq 2) then nu_dat_err2(0) = nu_dat_err2(1)*1.5
endif

;populations
if (populations eq 2) then basename = basename + '_popboth'
if (populations eq 1) then basename = basename + '_pop'+whichpop+''

;Gprior: 
if (gprior gt 0) then basename = basename + '_gprior'
if (gprior lt 1) then gprior = 1d30

;Cprior:
if (cprior GE 0) then basename = basename + '_cprior'
if (cprior lt 0) then cprior = 1d30

;Mirror, log, lb:
if (mirror eq 'yes') then basename = basename + '_mirr'
if (logprior eq 'yes') then basename = basename + '_log'
if (lbprior eq 'yes') then basename = basename + '_lb'

;betaprior
if (betap eq 'yes') then basename = basename + '_beta' 
if (mprior lt 0) then basename = basename + '_mslope' else basename = basename + '_mconst'

;Bprior:
if (bprior eq 'yes') then begin 
   if (bfifty eq 'no') then basename = basename + '_bprior' $
   else basename = basename + '_b50'

   if (baryonmodel eq 'sim') then begin
      readcol,massfile,rvis,massusevis,massuseviserr
   endif
   massvismin = massusevis - massuseviserr
endif
  
;S-prior:
if (sprior eq 'yes') then basename = basename + '_sprior'
   
;Binning in z:
if (rpmin lt 0) then rpmin = min(r_dat)
if (rpmax lt 0) then rpmax = max(r_dat)

rp_M = dindgen(nt_M) * (rpmax - rpmin) / double(nt_M-1) + rpmin
rp_nu = dindgen(nt_nu) * (rpmax - rpmin) / double(nt_nu-1) + rpmin
rp_beta = dindgen(nt_beta) * (rpmax - rpmin) / double(nt_beta-1) + rpmin

;Convert gprior into appropriate units: 
if (gprior gt 0) then $
   gpriorconv = gprior $ ;unit [Msun/kpc^3]
else gpriorconv = 1d30
if (cprior GE 0) then $
   cpriorconv = cprior $ ;unit [Msun/kpc^3]
else cpriorconv = 1d30

;Set up bprior:
if (bprior eq 'yes') then begin
   bmass = interpol(massvismin,rvis,rp_M) 	;interpolate on nodes zp_M, unit [Msun]
   dr = rp_M(2) - rp_M(1) 			;binlength
   blow = dblarr(n_elements(rp_M))
   for i=1L,n_elements(rp_M)-1 do $
      blow(i) = (bmass(i) - bmass(i-1))/dr 		; slope??
endif else blow = dblarr(n_elements(rp_M))

;Set up kzmin/max arrays:
if (Mmin lt 0) then $
   Mminarr = dblarr(nt_M) + abs(Mmin)*blow $
else Mminarr = dblarr(nt_M) + Mmin
if (Mmax lt 0) then $
   Mmaxarr = dblarr(nt_M) + abs(Mmax)*blow $
else Mmaxarr = dblarr(nt_M) + Mmax

;Default Initial seed guess parameters [assume flat]:
if (Mpars(0) lt 0) then Mpars = dblarr(nt_M) + .3  ;unit [Mtot]
if (nupars1(0) lt 0) then nupars1 = dblarr(nt_nu) + .3
if (betapars1(0) lt 0) then betapars1 =  dblarr(nt_beta)
if (nupars2(0) lt 0) then nupars2 = dblarr(nt_nu) + .3
if (betapars2(0) lt 0) then betapars2 =  dblarr(nt_beta)
if (mprior(0) lt 0) then Mslopepars = 0.5 else Mslopepars = mprior
if (sigmaprior1(0) lt 0) then sigmaslopepars1 = 0. 
if (sigmaprior2(0) lt 0) then sigmaslopepars2 = 0. 
;if (norm1 lt 0) then norm1 = 20.^2.
;if (norm2 lt 0) then norm2 = 20.^2.

;*************************************************************
;Metropolis parameters here:
if (populations eq 2) then pars = [nupars1,nupars2,Mpars,betapars1,betapars2,Mslopepars,sigmaslopepars1,sigmaslopepars2] else pars = [nupars1,Mpars,betapars1,Mslopepars,sigmaslopepars1]
npars = n_elements(pars)        ;total number of parameter
metarray = dblarr(npars+1,niter)
if (nuparstep1(0) lt 0) then nuparstep1 = nupars1 / 10.
if (nuparstep2(0) lt 0) then nuparstep2 = nupars2 / 10.
if (Mparstep(0) lt 0) then Mparstep = dblarr(nt_M) + 0.1 ;unit [Mtot]
if (betaparstep1(0) lt 0) then betaparstep1 = 0.02
if (betaparstep2(0) lt 0) then betaparstep2 = 0.02
if (mprior(0) lt 0) then Mslopeparstep = 0.02 else Mslopeparstep = 0.
if (sigmaprior1(0) lt 0) then sigmaslopeparstep1 = 0.02 
if (sigmaprior1(0) lt 0) then sigmaslopeparstep2 = 0.02 
;if (normstep1 lt 0) then normstep1 = norm1/10. 
;if (normstep1 lt 0) then normstep2 = norm2/10. 

if (populations eq 2) then parstep = [nuparstep1,nuparstep2,Mparstep,betaparstep1,betaparstep2,Mslopeparstep,sigmaslopeparstep1,sigmaslopeparstep2] else parstep = [nuparstep1,Mparstep,betaparstep1,Mslopeparstep,sigmaslopeparstep1]

;Logarithmic prior: 
if (logprior eq 'yes') then begin
   lpars=dblarr(npars)
   lpars(0:populations*nt_nu+nt_M-1) = alog10(pars(0:populations*nt_nu+nt_M-1)) ;log sampling for nu and M
   lpars(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta+1) = pars(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta+1) ;linear sampling for beta, Mslope and sigmaslope
   if (populations eq 2) then lpars(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta+2) = pars(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta+2) ;linear sampling for beta, Mslope and sigmaslope
   lparstep = lpars 
   lparstep(0:nt_nu-1) = 0.04 ;nu1
   if (populations eq 2) then lparstep(nt_nu:2*nt_nu-1) = 0.04 ;nu2
   lparstep(populations*nt_nu:populations*nt_nu+nt_M-1) = 0.04 ;M 
   lparstep(populations*nt_nu+nt_M:populations*nt_nu+nt_M+nt_beta-1) = 0.02 ;beta1
   if (populations eq 2) then lparstep(2*nt_nu+nt_M+nt_beta:2*nt_nu+nt_M+2*nt_beta-1) = 0.02 ;beta2
   lparstep(populations*nt_nu+nt_M+populations*nt_beta) = Mslopeparstep ;Mslope
   lparstep(populations*nt_nu+nt_M+populations*nt_beta+1) = sigmaslopeparstep1 ;sigmaslope
   if (populations eq 2) then lparstep(populations*nt_nu+nt_M+populations*nt_beta+2) = sigmaslopeparstep2 ;sigmaslope
;   lparstep(2*nt_nu+nt_M+2*nt_beta) = 1.25 ;norm1
;   lparstep(2*nt_nu+nt_M+2*nt_beta+1) = 1.25 ;norm2
endif

;Output data files:
outfilemet = basename+'.dat'
outfiledat = basename+'.txt'
if(debug eq 'yes') then print,'outfile',outfiledat
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
if(debug eq 'yes') then print,'after arraydump'

;*************************************************************
;METROPOLIS

chisq = 1d300                   ;Initial chisq [large]
n = 0L                          ;counter
initphase = 'start'             ;Initialisation phase flag, first 'start', not: 'over'
rejcount = 0.                   ;Rejection count
acccount = 0.                   ;Acceptance count
accrejtollow = 0.22             ;Acceptance/rejection rate
accrejtolhigh = 0.28            ;
if (populations eq 1) then chisqtol = 10. else chisqtol = 20.     ;Parameters to end initphase 
endcount = 200                  ;
endgame = 'no'                  ;Ending flag
while (n lt niter-1) do begin
   ;Store the old parameters:
   if (initphase eq 'over') then begin
      metarray(0:npars-1,n) = pars
      metarray(npars,n) = chisq
   endif 
   print,'n = ',n

AGAIN:
   ;Wiggle the parameters -> new parameters:
   if (logprior eq 'no') then begin 
      ranarr = 2.0 * (0.5 - randomu(seed,/DOUBLE,npars))
      parst = pars + ranarr * parstep
   endif else begin 
      ranarr = 2.0 * (0.5 - randomu(seed,/DOUBLE,npars))
      lparst = lpars + ranarr * lparstep 		;take one step in logspace
      parst = dblarr(npars)
      parst(0:populations*nt_nu+nt_M-1) = double(10)^lparst(0:populations*nt_nu+nt_M-1) ;log sampling for nu and M
      if ( betap eq 'no' ) then parst(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta-1) = lparst(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta-1) $  ;linear sampling for beta
      else parst(populations*nt_nu+nt_M:populations*nt_nu+nt_M+populations*nt_beta-1) = betaprior
      parst(populations*nt_nu+nt_M+populations*nt_beta:populations*nt_nu+nt_M+populations*nt_beta+populations) = lparst(populations*nt_nu+nt_M+populations*nt_beta:populations*nt_nu+nt_M+populations*nt_beta+populations) ;linear sampling for Mslope

      ;Mmin/max range prior:
      Mpars = parst(populations*nt_nu:populations*nt_nu+nt_M-1) 
      for jj=0L,nt_M-1 do begin
         if (Mpars(jj) lt Mminarr(jj)) then goto,AGAIN
         if (Mpars(jj) gt Mmaxarr(jj)) then goto,AGAIN
      endfor
   endelse

   ;Lbprior 
   if (lbprior eq 'yes') then begin 
      totmlastb = total(Mpars(0:nt_M-2) + blow(0:nt_M-2))
      lastb = Mpars(nt_M-1) + blow(nt_M-1)
      if (lastb / totmlastb gt lbtol) then goto,AGAIN
   endif

   ;Ensure positivity --> monotinicity constraint:
   if (mirror eq 'no') then begin
      for jj=0L,npars-1 do $
         if (parst(jj) lt 0.) then goto,AGAIN
   endif
   if(debug eq 'yes') then print,'Extract parameters:'
   nupars1 = parst(0:nt_nu-1)
   if (populations eq 2) then nupars2 = parst(nt_nu:2*nt_nu-1)
   Mpars = parst(populations*nt_nu:populations*nt_nu+nt_M-1)
   betapars1 = parst(populations*nt_nu+nt_M:populations*nt_nu+nt_M+nt_beta-1)
   if (populations eq 2) then betapars2 = parst(2*nt_nu+nt_M+nt_beta:2*nt_nu+nt_M+2*nt_beta-1)
   Mslopepars = parst(populations*nt_nu+nt_M+populations*nt_beta)
   sigmaslopepars1 = parst(populations*nt_nu+nt_m+populations*nt_beta+1)
   if (populations eq 2) then sigmaslopepars2 = parst(populations*nt_nu+nt_m+populations*nt_beta+2)
;   norm1 = parst(2*nt_nu+nt_M+2*nt_beta)
;   norm2 = parst(2*nt_nu+nt_M+2*nt_beta+1)
   
   ;Apply "cprior" and "gprior" on M function:
   ;if (cpriorconv gt 0) then begin
   ;   if (abs(Mpars(0)) gt cpriorconv) then goto,AGAIN
   ;endif else begin 
   ;   Mpars(0) = Mminarr(0)
   ;   parst(2*nt_nu) = Mminarr(0)
   ;endelse

   if (abs(Mslopepars) gt 2.) then goto,AGAIN
   for jj=0L,nt_beta-1 do if (betapars1(jj) gt 1.) then goto, AGAIN
   if (populations eq 2) then for jj=0L,nt_beta-1 do if (betapars1(jj) gt 1.) then goto, AGAIN

   for jj=0L,nt_M-1 do $
      if (abs(Mpars(jj)) gt gpriorconv) then goto,AGAIN

   if(debug eq 'yes') then print,'Calculate theoretical values: '
   nu1_r = nu(r_dat,rp_nu,nupars1)
   if (populations eq 2) then nu2_r = nu(r_dat,rp_nu,nupars2)
   sig1_r = sigma_los1(r_sigdat,rp_nu,rp_M,rp_beta,Mpars,nupars1,betapars1,Mslopepars,sigmaslopepars1)
   if (populations eq 2) then sig2_r = sigma_los1(r_sigdat,rp_nu,rp_M,rp_beta,Mpars,nupars2,betapars2,Mslopepars,sigmaslopepars2)

   if(debug eq 'yes') then print,'Reject models with NaN sig_z: '
   if (logprior eq 'yes') then begin 
      if(debug eq 'yes') then print,'Nan condition'
      for jj=0L,n_elements(sig1_r)-1 do $
         if (sig1_r(jj) NE sig1_r(jj)) then goto, AGAIN
      if (populations eq 2) then begin
         for jj=0L,n_elements(sig2_r)-1 do $
            if (sig2_r(jj) NE sig2_r(jj)) then goto, AGAIN
      endif    
   endif else begin 
      ;Another problem with linear sampling (best avoided). 
      ;If I do the "right" thing here and reject NaN models, 
      ;the MCMC gets stuck. So we set sig_z = 0. where it
      ;is NaN and assume that this will be penalised by the 
      ;data. This assumption appears to be very good, 
      ;but still it's not ideal ): 
      for jj=0L,n_elements(sig1_r)-1 do $
         if (sig1_r(jj) NE sig1_r(jj)) then sig1_r(jj) = 0.
      if (populations eq 2) then begin
         for jj=0L,n_elements(sig2_r)-1 do $
            if (sig2_r(jj) NE sig2_r(jj)) then sig2_r(jj) = 0.
      endif    
   endelse
   if(debug eq 'yes') then print,'before sprior'
   ;S-prior: ensure sigma_z(z) rises: 
   if (sprior eq 'yes') then begin 
      for jj=1L,n_elements(sig_z)-1 do $
         if (sig_z(jj) lt sig_z(jj-1)) then begin 
         goto,AGAIN
         endif
   endif
   if(debug eq 'yes') then print,'Calculate chi-squared'
   if(debug eq 'yes') then print,'  nu_dat_err1 = ',total(1./nu_dat_err1)
   if(debug eq 'yes') then print,'  nu_dat_err2 = ',total(1./nu_dat_err2)
   if(debug eq 'yes') then print,'  sig_dat_err1 = ',total(1./sig_dat_err1)
   if(debug eq 'yes') then print,'  sig_dat_err2 = ',total(1./sig_dat_err2)

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


   if(debug eq 'yes') then print,'Calculate the f-function, chisq = ',chisq,', chisqt = ',chisqt
   fnewoverf = EXP(chisq/2.0-chisqt/2.0)
    
   if(debug eq 'yes') then print,'Accept the new f-function?'
   ran = randomu(seed,/DOUBLE)
   print,'random seed = ',ran,', fnewoverf = ',fnewoverf
   if (ran lt fnewoverf) then begin
      acccount = acccount + 1.
      pars = parst
      chisq = chisqt

      if (logprior eq 'yes') then lpars = lparst

      if (testplot eq 'yes') then begin
         ;Calculate profiles: 
         mprioru = abs(Mslopepars)
         rplmin = 0;
         rplmax = max(rp_M)
         rpnts = 200
         rpl = dindgen(rpnts) * (rplmax-rplmin)/double(rpnts-1) + rplmin
         dr = (rplmax-rplmin)/double(rpnts-1)
         nu1 = nu(rpl,rp_nu,nupars1)
         sig1 = sigma_los1(rpl,rp_nu,rp_M,rp_beta,Mpars,nupars1,betapars1,Mslopepars,sigmaslopepars1)
         M = Mr(rpl,rp_M,Mpars)
         M_outer=M(rpnts-1)+(dindgen(rpnts-1)+1)*mprioru*M(rpnts-1)/(rplmax-rplmin)*dr
         nu1_outer=nu1(rpnts-1)-(dindgen(rpnts-1)+1)*nu1(rpnts-1)/(rplmax-rplmin)*dr
         ;sigr1 = sigma_r1(rpl,rp_nu,rp_M,rp_beta,Mpars,nupars1,betapars1,Mslopepars)
         ;sigr1_outer=sigr1(rpnts-1)+(dindgen(rpnts-1)+1)*sigmaslopepars1*sigr1(rpnts-1)/(rplmax-rplmin)*dr
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
         ;sigr1_tot=dblarr(2*rpnts-1)
         ;sigr1_tot(0:rpnts-1)=sigr1
         ;sigr1_tot(rpnts:2*rpnts-2)=sigr1_outer

         ;Do the plots: 
         plot,r_dat,nu_dat1,xrange=[0,rplmax],yrange=[0.000001,1.5],/ylog,$
              title='!6',xtitle='r(Rs)',ytitle='nu1 [unity]'
         oploterror,r_dat,nu_dat1,nu_dat_err1
         oplot,r_tot,nu1_tot,color=2

         plot,r_dat,sig_dat1,xrange=[0,rplmax],yrange=[0,1.],/NODATA,$
              title='!6',xtitle='r(Rs)',ytitle='sigma_los1 []'
         oploterror,r_sigdat,sig_dat1,sig_dat_err1
         oplot,rpl,sig1,color=2

         if (populations eq 2) then begin
            sig2 = sigma_los1(rpl,rp_nu,rp_M,rp_beta,Mpars,nupars2,betapars2,Mslopepars,sigmaslopepars2)
            ;sigr2 = sigma_r1(rpl,rp_nu,rp_M,rp_beta,Mpars,nupars2,betapars2,Mslopepars)
            ;sigr2_outer=sigr2(rpnts-1)+(dindgen(rpnts-1)+1)*sigmaslopepars2*sigr2(rpnts-1)/(rplmax-rplmin)*dr
            nu2 = nu(rpl,rp_nu,nupars2)
            beta2 = betar(rpl,rp_beta,betapars2)
            ;sigr2_tot=dblarr(2*rpnts-1)
            ;sigr2_tot(0:rpnts-1)=sigr2
            nu2_outer=nu2(rpnts-1)-(dindgen(rpnts-1)+1)*nu2(rpnts-1)/(rplmax-rplmin)*dr
            ;sigr2_tot(rpnts:2*rpnts-2)=sigr2_outer
            nu2_tot=dblarr(2*rpnts-1)
            nu2_tot(0:rpnts-1)=nu2
            nu2_tot(rpnts:2*rpnts-2)=nu2_outer

            plot,r_dat,nu_dat2,xrange=[0,rplmax],yrange=[0.000001,1.5],/ylog,$
                 title='!6',xtitle='r(Rs)',ytitle='nu2 [unity]'
            oploterror,r_dat,nu_dat2,nu_dat_err2
            oplot,rpl,nu2,color=2

            plot,r_dat,sig_dat2,xrange=[0,rplmax],yrange=[0,1.],/NODATA,$
                 title='!6',xtitle='r(Rs)',ytitle='sigma_los2 []'
            oploterror,r_sigdat,sig_dat2,sig_dat_err2
            oplot,rpl,sig2,color=2
         endif
         plot,r_Mdat,M_dat,xrange=[0,2*rplmax],yrange=[0,2.],/NODATA,$
              title='!6',xtitle='r(Rs)',ytitle='M [Mtot]'
         oploterror,r_Mdat,M_dat,M_dat_err
         oplot,r_tot,M_tot,color=2

         plot,rpl,beta1,xrange=[rplmin,rplmax],yrange=[-1.5,1.5],/NODATA,$
              title='!6',xtitle='r(Rs)',ytitle='beta [unity]'
         if (populations eq 2) then oplot,rpl,beta2,color=3
         oplot,rpl,beta1,color=2

      endif
      
      if(debug eq 'yes') then print,'Output for testing: '
  ;    if (logprior eq 'no') then stepout = parstep(nt_kz-1) $
  ;    else stepout = lparstep(nt_kz-1)
      print,n,chisq,lparstep(0),acccount/rejcount

      ;Decide whether to end initphase:
      if (endgame eq 'yes') then begin
         endcount = endcount - 1
         if (endcount eq 0) then begin
            initphase = 'over'
            print,''
            print,'********** Initialisation phase over **********'
            print,'Step size chosen:'
            if (logprior eq 'no' ) then begin
               print,parstep
               print,pars
            endif else begin
               print,lparstep
               print,lpars
               print,pars
            endelse
            print,''
         endif
      endif
   endif else begin 
      print,'Rejected',rejcount
      rejcount = rejcount + 1.
   endelse

   if(debug eq 'yes') then print,'Adapt stepsize during initialisation phase: '
   if (initphase eq 'start') then begin
      if (acccount gt 0) AND (rejcount gt 0) then begin 
         if (acccount/rejcount lt accrejtollow) then begin 
            parstep = parstep / 1.001
            if (logprior eq 'yes') then lparstep = lparstep / 1.01
         endif 
         if (acccount/rejcount gt accrejtolhigh) then begin
            parstep = parstep * 1.001
            if (logprior eq 'yes') then lparstep = lparstep * 1.01
         endif
      endif
      if (chisq lt chisqtol) then endgame = 'yes'
   endif

   ;Update the counter:
   if (initphase eq 'over') then n = n + 1
endwhile

;Write the data to a file:
openw,1,outfilemet
writeu,1,LONG(npars),LONG(niter)
writeu,1,metarray
close,1

print,'Finished!'

;close the output file, set the plotting back to the screen and exit.
if devtype eq 'PS' then device,/close	
set_plot, 'X'
!P.FONT=-1
!P.multi=0

end
