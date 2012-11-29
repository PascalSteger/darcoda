;============================ PLOT CLEAN ==========================
;Program to plot the output from silvia_clean.pro

PRO plot_spherical2

setcommon
common constants

;Decide if the plot is to be output to the screen or to a
;postcript file.
devtype='X'                    ;Output device type
set_plot,devtype
;IF (devtype EQ 'PS') THEN set_thick,4.,-1,2. ELSE $
;   set_thick,1,-1,2.
IF (devtype EQ 'X') THEN filetype = '.png' ELSE filetype = '.eps'
!X.STYLE=1
!Y.STYLE=1
!P.multi=0

IF (devtype EQ 'X') THEN BEGIN
   device,decomposed=0
   fxsize = 900
   fysize = 900
   window,0,xsize=fxsize,ysize=fysize,retain=2
ENDIF

;Set colours:
mymakecolour

populations=1
whichpop='dm'

;Sim or data: 
plotsim = 'yes'

;Order chisq before burnin cut? 
chisqorder = 'no'

;Monocut? 
monocut = 'no'

;Modify sigma_s,max? [-1 = no]
modsigmax = -1
;modsigmax = 50.

;Default burnin: 
burninfac = 0.02

;Read in the metarray data file:
dir="/scratch/psteger/software/mcmc"
 
IF (plotsim EQ 'yes') THEN BEGIN
   nt_M = 25
   nt_nu = 25
   nt_beta = 25
   Mno = '25'
   nuno = '25'
   betano = '25'
   whichsim = '1'
   massfile = dir+'/sphericalmethod_data/enclosedmass/dual_unit_hern_'+whichsim+'_tot_enclosedmass_half.txt'
   if (populations eq 2) then begin
      basename = dir+'/sphericalmethod_output/dual_unit_hern_'+whichsim+'/M'+Mno+'_nu'+nuno+'_beta'+betano+'/met_popboth_cprior_mirr_log_beta'
      nufile1  = dir+'/sphericalmethod_data/densityfalloff/dual_unit_hern_'+whichsim+'_stars_falloff.txt'
      sigfile1 = dir+'/sphericalmethod_data/velocitydispersion/dual_unit_hern_'+whichsim+'_stars_veldisp_half.txt'
      nufile2  = dir+'/sphericalmethod_data/densityfalloff/dual_unit_hern_'+whichsim+'_dm_falloff.txt'
      sigfile2 = dir+'/sphericalmethod_data/velocitydispersion/dual_unit_hern_'+whichsim+'_dm_veldisp_half.txt'
   endif
   IF (populations eq 1) then begin
      basename = dir+'/sphericalmethod_output/dual_unit_hern_'+whichsim+'/M'+Mno+'_nu'+nuno+'_beta'+betano+'/met_pop'+whichpop+'_cprior_mirr_log_beta'
      nufile1  = dir+'/sphericalmethod_data/densityfalloff/dual_unit_hern_'+whichsim+'_'+whichpop+'_falloff.txt'
      sigfile1 = dir+'/sphericalmethod_data/velocitydispersion/dual_unit_hern_'+whichsim+'_'+whichpop+'_veldisp_half.txt'
   ENDIF
   dosim = 'yes'
ENDIF

infilemet = basename+'.dat'
infiledat = basename+'.txt'
infileprofM = basename+'.profM'
infileprofnu1 = basename+'.profnu1'
infileprofbeta1 = basename+'.profbeta1'
readcol,infiledat,nt_nu,nt_M,nt_beta,niter,numline=1,skipline=1
readcol,infileprofM,rp_M
readcol,infileprofnu1,rp_nu
readcol,infileprofbeta1,rp_beta
npars = LONG(populations*nt_nu+nt_M+populations*nt_beta)
npars = LONG(npars)
niter = LONG(niter)
metarray = $
  dblarr(npars+1,niter)

openr,1,infilemet
readu,1,npars,niter
readu,1,metarray
close,1

IF (plotsim EQ 'yes') THEN BEGIN  
   ;Read in simulation data here: 
   silv = [-10,0,0]
   readcol,massfile,r_Mdat,M_dat,M_dat_err

   ysigmin = 0.
   ysigmax = 0.8
   ynumin = 0.000001
   ynumax = 2.
   ybetamin = -2.
   ybetamax = 2.
   Mmax = 1.
ENDIF

readcol,nufile1,r_dat,nu_dat1,nu_dat_err1
if (populations eq 2) then readcol,nufile2,r_dat,nu_dat2,nu_dat_err2
readcol,sigfile1,r_sigdat,sig_dat1,sig_dat_err1
if (populations eq 2) then readcol,sigfile2,r_sigdat,sig_dat2,sig_dat_err2

;Renormalise density to one:
;normden = max(nu_dat)
;nu_dat = nu_dat / normden
;nu_dat_err = nu_dat_err / normden

IF (plotsim EQ 'yes') THEN BEGIN
   ;Fix central error to match silvia_clean.pro
   nu_dat_err1(0) = nu_dat_err1(1)*1.5
   if (populations eq 2) then nu_dat_err2(0) = nu_dat_err2(1)*1.5
ENDIF

;Set burnin: 
burnin = burninfac * niter

IF (chisqorder EQ 'yes') THEN BEGIN 
   ;Extract post burn-in array parmaeters:
   index = sort(metarray(npars,*))
   burn = niter - burnin 
   t = index(0:burn)
   temp = metarray(npars,*)
   chimods = temp(t)
   temp = metarray(nt_nu:nt_nu+nt_kz-1,*)
   kzpars = temp(*,t)
   temp = metarray(0:nt_nu-1,*)
   nupars = temp(*,t)
   temp = metarray(nt_nu+nt_kz:nt_nu+nt_kz,*)
   norm1 = temp(t)
ENDIF ELSE BEGIN 
   chimods = metarray(npars,burnin:niter-1)
   Mpars = metarray(populations*nt_nu:populations*nt_nu+nt_M-1,burnin:niter-1)
   nupars1 = metarray(0:nt_nu-1,burnin:niter-1)
   if (populations eq 2) then nupars2 = metarray(nt_nu:2*nt_nu-1,burnin:niter-1)
   betapars1 = metarray(populations*nt_nu+nt_M:populations*nt_nu+nt_M+nt_beta-1,burnin:niter-1)
   if (populations eq 2) then betapars2 = metarray(2*nt_nu+nt_M+nt_beta:2*nt_nu+nt_M+2*nt_beta-1,burnin:niter-1)
ENDELSE 

nn = niter-burnin
nn = nn(0)

;Set the r array:
nr=n_elements(rp_M)+1
rpl = dblarr(nr)
for i=0L,nr-1 do rpl(i)=i*max(rp_M)/double(nr-1)
rpnts = n_elements(rpl)
print,rpl

;Calcualte averaged profiles:
M_r = dblarr(rpnts,nn)
nu1_r = dblarr(rpnts,nn)
if (populations eq 2) then nu2_r = dblarr(rpnts,nn)
sig1_r = dblarr(rpnts,nn)
if (populations eq 2) then sig2_r = dblarr(rpnts,nn)
beta1_r = dblarr(rpnts,nn)
if (populations eq 2) then beta2_r = dblarr(rpnts,nn)

;Decide probability boundaries to store:
med = 0.5
int1 = 0.95
int2 = 0.68
p = [med-int1/2.,med-int2/2.,med,med+int2/2.,med+int1/2.]
pstyle = [2,1,0,1,2]
np = n_elements(p)

M_rp = dblarr(np,rpnts)
nu1_rp = dblarr(np,rpnts)
if (populations eq 2) then nu2_rp = dblarr(np,rpnts)
sig1_rp = dblarr(np,rpnts)
if (populations eq 2) then sig2_rp = dblarr(np,rpnts)
beta1_rp = dblarr(np,rpnts)
if (populations eq 2) then beta2_rp = dblarr(np,rpnts)

;Make chisq-sample plots here: 
oname = basename + '_chisqMmax' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

;Make a contour plot: 
x = alog10(metarray(npars,burnin:niter-1))
y = alog10(abs(Mpars(nt_M-1,*)))
xmax = alog10(70.)
xmin = 0.3
ymax = 6
ymin = -8
dx = 0.1
mass = dblarr(n_elements(x)) + 1.
contourmemfast,dx,xmax,xmin,ymax,ymin,x,y,mass,xnew,ynew,denxy
denxy = denxy / max(denxy)
nlevels = 10
clabels = dblarr(nlevels)
FOR jj=0L,nlevels-1,2 DO clabels(jj) = 1

contour,denxy,xnew,ynew,nlevels=nlevels,$
        title='!6',xtitle='log10[!7v!6!u2!n]',$
        ytitle='log10[M!dr!n(r!dmax!n) ; (code units)]',$
        c_labels=clabels

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

;Arrange data:
pp = 0L
FOR i=0L,nn-1 DO BEGIN
   M_r(*,pp) = mr(rpl,rp_M,Mpars(*,i))
   nu1_r(*,pp) = nu(rpl,rp_nu,nupars1(*,i))
   if (populations eq 2) then nu2_r(*,pp) = nu(rpl,rp_nu,nupars2(*,i))
   sig1_r(*,pp) = sigma_r(rpl,rp_nu,rp_M,rp_beta,Mpars(*,i),nupars1(*,i),betapars1(*,i))
   if (populations eq 2) then sig2_r(*,pp) = sigma_r(rpl,rp_nu,rp_M,rp_beta,Mpars(*,i),nupars2(*,i),betapars2(*,i))
   beta1_r(*,pp) = betar(rpl,rp_beta,betapars1(*,i))
   if (populations eq 2) then beta2_r(*,pp) = betar(rpl,rp_beta,betapars2(*,i))
 
   IF (monocut EQ 'no') THEN pp = pp + 1 ELSE BEGIN 
      flag = 0
      dzpl = zpl(3)-zpl(2)
      FOR jj=1L,n_elements(zpl)-1 DO BEGIN
         grad = (sig_z(jj,pp) - sig_z(jj-1,pp)) / dzpl
         IF (grad LT -10.) OR (grad GT 10.) THEN flag = 1
      ENDFOR
      IF (flag EQ 0) THEN pp = pp + 1
   ENDELSE
ENDFOR
IF (monocut EQ 'yes') THEN BEGIN
   print,''
   print,'*********** MONOCUT PRIOR *************'
   print,'Number of models before cut [after burnin]:', nn
   print,'Number of models after cut & burnin:', pp
   print,''
   nn = pp
ENDIF

;Calculate median 68% and 90% for each bin:
FOR i=0L,rpnts-1 DO BEGIN 
   indexM = sort(M_r(i,0:nn-1))
   indexnu1 = sort(nu1_r(i,0:nn-1))
   if (populations eq 2) then indexnu2 = sort(nu2_r(i,0:nn-1))
   indexsg1 = sort(sig1_r(i,0:nn-1))
   if (populations eq 2) then indexsg2 = sort(sig2_r(i,0:nn-1))
   indexbeta1 = sort(beta1_r(i,0:nn-1))
   if (populations eq 2) then indexbeta2 = sort(beta2_r(i,0:nn-1))

   FOR j=0L,np-1 DO BEGIN
      M_rp(j,i) = M_r(i,indexM(p(j)*(nn-1)))
      nu1_rp(j,i) = nu1_r(i,indexnu1(p(j)*(nn-1)))
      if (populations eq 2) then nu2_rp(j,i) = nu2_r(i,indexnu2(p(j)*(nn-1)))
      sig1_rp(j,i) = sig1_r(i,indexsg1(p(j)*(nn-1)))
      if (populations eq 2) then sig2_rp(j,i) = sig2_r(i,indexsg2(p(j)*(nn-1)))
      beta1_rp(j,i) = beta1_r(i,indexbeta1(p(j)*(nn-1)))
      if (populations eq 2) then beta2_rp(j,i) = beta2_r(i,indexbeta2(p(j)*(nn-1)))
   ENDFOR
ENDFOR

IF (monocut EQ 'yes') THEN basename = basename + '_mono'

; default modsigmax=-1
IF (modsigmax GT 0) THEN basename = basename + $
     strcompress(string(modsigmax),/remove_all)


;Array dump the processed data: 
rplt = transpose(rpl)
arraydump,basename+'_M_rp.txt',[[rplt],[M_rp(0,*)],[M_rp(1,*)],[M_rp(2,*)],$
                                     [M_rp(3,*)],[M_rp(4,*)]],6
arraydump,basename+'_nu1_rp.txt',[[rplt],[nu1_rp(0,*)],[nu1_rp(1,*)],[nu1_rp(2,*)],$
                                 [nu1_rp(3,*)],[nu1_rp(4,*)]],6
if (populations eq 2) then arraydump,basename+'_nu2_rp.txt',[[rplt],[nu2_rp(0,*)],[nu2_rp(1,*)],[nu2_rp(2,*)],$
                                 [nu2_rp(3,*)],[nu2_rp(4,*)]],6
arraydump,basename+'_sig1_rp.txt',[[rplt],[sig1_rp(0,*)],[sig1_rp(1,*)],[sig1_rp(2,*)],$
                                  [sig1_rp(3,*)],[sig1_rp(4,*)]],6
if (populations eq 2) then arraydump,basename+'_sig2_rp.txt',[[rplt],[sig2_rp(0,*)],[sig2_rp(1,*)],[sig2_rp(2,*)],$
                                  [sig2_rp(3,*)],[sig2_rp(4,*)]],6
arraydump,basename+'_beta1_rp.txt',[[rplt],[beta1_rp(0,*)],[beta1_rp(1,*)],[beta1_rp(2,*)],$
                                 [beta1_rp(3,*)],[beta1_rp(4,*)]],6
if (populations eq 2) then arraydump,basename+'_beta2_rp.txt',[[rplt],[beta2_rp(0,*)],[beta2_rp(1,*)],[beta2_rp(2,*)],$
                                 [beta2_rp(3,*)],[beta2_rp(4,*)]],6


;Do the plots: 
oname = basename + '_chisq' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

redchi = metarray(npars,*)

;Calculate number of parameters only over the range of the data: 
;Note: this calculation *assumes* that np_kz = np_nu and that 
;max(zp_kz) = max(z_dat).
;jl = 0L
;while (rp_M(jl) lt min(r_dat)) do jl=jl+1
;parnumber = npars(0) - 2.0 * jl
;redfac = populations*2.*n_elements(r_dat) - parnumber - 1.0

;IF (redfac LT 0) THEN BEGIN
;   print,'*** Warning *** reduced chisq not defined. Too many pars.' 
;   chilabel = '!7v!6!u2!n!dred!n'
;ENDIF ELSE BEGIN
;   redchi = redchi / redfac
;   chilabel = '!7v!6!u2!n!dred!n'
;ENDELSE
xarr = dindgen(n_elements(redchi))
plot,xarr/1.e4,redchi,yrange=[100,10000],/ylog,title='!6',$
     xtitle='iterations / 10!u4!n',ytitle=chilabel

;Overplot chosen burnin: 
oplot,[burnin,burnin]/1.e4,[1d-300,1d300],color=2

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

oname = basename + '_chisqhist' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

;Cut the burnin:
;redchicut = chimods / redfac

h = histogram(chimods,binsize=0.1,min = 0., max = 70., locations = loc)
plot,loc,h,psym=10,title='!6',xtitle=chilabel,ytitle='N',$
     xrange=[0,7]

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE


;Mr plot
oname = basename + '_Mr' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

FOR i=0L,np-1 DO BEGIN
   IF (i EQ 0) THEN BEGIN 
      plot,rpl,M_rp(i,*),title='!6',xtitle='r(Rs)',$
           ytitle='M!dr!n [Msun]',$
           yrange=[0,Mmax],xrange=[0,max(r_Mdat)],$
           linestyle = pstyle(i)
      ;Overplot max baryonic surface density:
;      oplot,[0,5],[sigs(0)-sigs(1),sigs(0)-sigs(1)],linestyle=1,color=6
;      oplot,[0,5],[sigs(0),sigs(0)],color=6
;      oplot,[0,5],[sigs(0)+sigs(1),sigs(0)+sigs(1)],linestyle=1,color=6

      IF (plotsim EQ 'yes') THEN BEGIN
         oplot,r_Mdat,M_dat,color=2
      ENDIF
   ENDIF ELSE $
      oplot,rpl,M_rp(i,*),linestyle=pstyle(i)
ENDFOR

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

;sigma1 plot
oname = basename + '_sigma1r' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

FOR i=0L,np-1 DO BEGIN
   IF (i EQ 0) THEN BEGIN
      plot,rpl,sig1_rp(i,*),title='!6',xtitle='r(Rs)',$
           ytitle='!7r!6!dr!n [unit]',$
           yrange=[ysigmin,ysigmax],xrange=[0,max(r_sigdat)],$
           linestyle = pstyle(i)

      ;Overplot data:
      oploterror,r_sigdat,sig_dat1,sig_dat_err1,color=2,errcol=2
   ENDIF ELSE $
      oplot,rpl,sig1_rp(i,*),linestyle=pstyle(i)
ENDFOR

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

;sigma2 plot
if (populations eq 2) then begin
   oname = basename + '_sigma2r' + filetype
   IF (devtype EQ 'PS') THEN $
      device,filename=oname,/color,$
             xsize=16,ysize=16,/encapsul

   FOR i=0L,np-1 DO BEGIN
      IF (i EQ 0) THEN BEGIN
         plot,rpl,sig2_rp(i,*),title='!6',xtitle='r(Rs)',$
              ytitle='!7r!6!dr!n [unit]',$
              yrange=[ysigmin,ysigmax],xrange=[0,max(r_sigdat)],$
              linestyle = pstyle(i)

         ;Overplot data:
         oploterror,r_sigdat,sig_dat2,sig_dat_err2,color=2,errcol=2
      ENDIF ELSE $
         oplot,rpl,sig2_rp(i,*),linestyle=pstyle(i)
   ENDFOR

   IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
      image=TVRD(0,0,fxsize,fysize,TRUE=1)
      write_png,oname,image
   ENDELSE
endif

;nu1 plot
oname = basename + '_nu1r' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

FOR i=0L,np-1 DO BEGIN
   IF (i EQ 0) THEN BEGIN
      plot,rpl,nu1_rp(i,*),title='!6',xtitle='r(Rs)',$
           ytitle='!7m!6 [units]',yrange=[ynumin,ynumax],/ylog,$
           xrange=[0,max(r_dat)],$
           linestyle = pstyle(i)
      
      ;Overplot data:
      oploterror,r_dat,nu_dat1,nu_dat_err1,color=2,errcol=2
   ENDIF ELSE $
      oplot,rpl,nu1_rp(i,*),linestyle=pstyle(i)
ENDFOR

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

;nu2 plot
if (populations eq 2) then begin
   oname = basename + '_nu2r' + filetype
   IF (devtype EQ 'PS') THEN $
      device,filename=oname,/color,$
             xsize=16,ysize=16,/encapsul

   FOR i=0L,np-1 DO BEGIN
      IF (i EQ 0) THEN BEGIN
         plot,rpl,nu2_rp(i,*),title='!6',xtitle='r(Rs)',$
              ytitle='!7m!6 [units]',yrange=[ynumin,ynumax],/ylog,$
              xrange=[0,max(r_dat)],$
              linestyle = pstyle(i)
         
         ;Overplot data:
         oploterror,r_dat,nu_dat2,nu_dat_err2,color=2,errcol=2
      ENDIF ELSE $
         oplot,rpl,nu2_rp(i,*),linestyle=pstyle(i)
   ENDFOR

   IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
      image=TVRD(0,0,fxsize,fysize,TRUE=1)
      write_png,oname,image
   ENDELSE
endif

;beta1 plot
oname = basename + '_beta1r' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

FOR i=0L,np-1 DO BEGIN
   IF (i EQ 0) THEN BEGIN
      plot,rpl,beta1_rp(i,*),title='!6',xtitle='r(Rs)',/ylog,$
           ytitle='beta [unity]',yrange=[ybetamin,ybetamax],$
           xrange=[0,max(r_dat)],$
           linestyle = pstyle(i)
      
   ENDIF ELSE $
      oplot,rpl,beta1_rp(i,*),linestyle=pstyle(i)
ENDFOR

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

;beta2 plot
if (populations eq 2) then begin
   oname = basename + '_beta2r' + filetype
   IF (devtype EQ 'PS') THEN $
      device,filename=oname,/color,$
             xsize=16,ysize=16,/encapsul

   FOR i=0L,np-1 DO BEGIN
      IF (i EQ 0) THEN BEGIN
         plot,rpl,beta2_rp(i,*),title='!6',xtitle='r(Rs)',/ylog,$
              ytitle='beta [unity]',yrange=[ybetamin,ybetamax],$
              xrange=[0,max(r_dat)],$
              linestyle = pstyle(i)
         
      ENDIF ELSE $
         oplot,rpl,beta2_rp(i,*),linestyle=pstyle(i)
   ENDFOR

   IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
      image=TVRD(0,0,fxsize,fysize,TRUE=1)
      write_png,oname,image
   ENDELSE
endif

;close the output file, set the plotting back to the screen and exit.
IF devtype EQ 'PS' THEN device,/close	
tidy

END
