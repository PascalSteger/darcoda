;============================ PLOT_GRAVIMAGE =========================
;Program to plot the output from grav_image_disc.pro

PRO plot_gravimage

;Set up the common block:
setcommon
COMMON constants

;Decide if the plot is to be output to the screen or to a
;postcript file.
devtype='PS'                    ;Output device type
set_plot,devtype
IF (devtype EQ 'PS') THEN set_thick,4.,-1,2. ELSE $
   set_thick,1,-1,2.
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

;Units:
G1 = 6.67d-11 * 1.989d30 / 3.086d19 / 1000.^2.

;Sim or data: 
plotwhich = 'sim'               ;What to load+plot
usetpar = 'no'                  ;Using tilt? 
bprior = 'yes'                  ;Using "bprior"?
med = 0.5                       ;Prob. boundaries 
int1 = 0.95                     ;to store 
int2 = 0.68                     ;
quadratic = 'yes'               ;Defaults 
monotonic = 'no'                ;

IF (plotwhich EQ 'sim') THEN baryonmodel = 'sim'
IF (plotwhich EQ 'kg') THEN baryonmodel = 'silvia'
IF (plotwhich EQ 'simple') THEN baryonmodel = 'simple'

;Burn-in [should not be req. if using "initphase"]:
burninfac = 0.01

;Read in the metarray data file:
IF (plotwhich EQ 'sim') THEN BEGIN
   patch = '135'
   basename = '/home/ast/user/jread/Data/Local_dm/Met_data/6_6_'+patch+$
              '/metdist_nolike8_mirr_log_bprior_constdm_rprior_kprior_quad'
   surfsim =  '/home/ast/user/jread/Data/Local_dm/Sim_data/Surfaceden/SigmaS_z_a'+patch+'.txt'
   nufile = '/home/ast/user/jread/Data/Local_dm/Sim_data/j_den_'+patch+'.txt'
   sigfile = '/home/ast/user/jread/Data/Local_dm/Sim_data/j_sigma_'+patch+'.txt'
   rhofile = '/home/ast/user/jread/Data/Local_dm/Sim_data/j_rhoSig_'+patch+'.txt'
   surfsfile = '/home/ast/user/jread/Data/Local_dm/Sim_data/Surfaceden/SigmaS_z_a'+patch+'.txt'
   surfdmfile = '/home/ast/user/jread/Data/Local_dm/Sim_data/Surfaceden/SigmaDM_z_a'+patch+'.txt'
   datafile = '/home/ast/read/dark/sgarbari/Justin/FETTAden+vz_a'+patch+'_4p_R250_r85_vphi.txt'

   readcol,datafile,mass,x_dat,y_dat,z_dat,vz_dat,vx_dat,vy_dat,pot_dat
   IF (max(mass) NE min(mass)) THEN BEGIN
      print,'**** Multimass data not yet supported ****'
      stop
   ENDIF
ENDIF
IF (plotwhich EQ 'kg') THEN BEGIN
   basename = '/home/ast/user/jread/Data/Local_dm/Met_data/6_6/met_cprior_mirr_bprior'
ENDIF
IF (plotwhich EQ 'simple') THEN BEGIN
   basename = $
      '/home/ast/user/jread/Data/Local_dm/Met_data/Simple/metdist_2compdd_1e48_mirr_log_bprior_rprior_kprior_quad'
   nufile = basename + '_nu_dat_bin.txt'
   sigfile = basename + '_sigz_dat_bin.txt'
ENDIF
infilemet = basename+'.dat'
infiledat = basename+'.txt'
infileprofkz = basename+'.profkz'
readcol,infiledat,nt_kz,nt_nu,niter,numline=1,skipline=1
readcol,infileprofkz,zp_kz
npars = LONG(nt_kz+nt_nu+1)
npars = LONG(npars)
niter = LONG(niter)
metarray = $
  dblarr(npars+1,niter)

openr,1,infilemet
readu,1,npars,niter
readu,1,metarray
close,1

;And tpar file if required: 
tpars = dblarr(4,niter)
tpars(0,*) = -1
ntpars = 1
IF (usetpar EQ 'yes') THEN BEGIN
   ntpars = LONG(0)
   dummy = LONG(0)

   openr,1,infilemet+'_tpars'
   readu,1,ntpars,dummy
   
   tpars = dblarr(ntpars+3,dummy)
   readu,1,tpars
   close,1
ENDIF

;Fix Obiwon bug [ought to fix this in grav_image_disc really]:
niter = niter - 1

;Bprior:
IF (bprior EQ 'yes') THEN BEGIN
   IF (baryonmodel EQ 'silvia') THEN BEGIN
      readcol,'/home/ast/user/jread/Data/Local_dm/Vis/Sigma_MM.txt',$
              zvis,sigexpvis,sigexpviserr,sigsecvis,sigsecviserr
      sigusevis = sigsecvis
      siguseviserr = sigsecviserr
   ENDIF
   IF (baryonmodel EQ 'sim') THEN BEGIN
      surfsfile = surfsim
      readcol,surfsfile,zvis,sigusevis,siguseviserr
   ENDIF
   IF (baryonmodel EQ 'simple') THEN BEGIN
      zmin = 0.
      zmax = 1.2
      zpnts = 5000
      zth = dindgen(zpnts) * (zmax-zmin)/double(zpnts-1) + zmin
      z0 = 0.25
      K = 1.0 * 1650.
      F = 0.1 * 1650.
      C = 17.^2.
      D = 0.25
      z02 = 0.4      

      ;Model:
      nu_zth = exp(-zth/z0)
      Kz_zth = -(K*zth/sqrt(zth^2.+D^2.) + 2.0 * F * zth)
      adddarkdisc = 'no'
      IF (adddarkdisc EQ 'yes') THEN BEGIN
         DD = 0.6
         KD = 0.15 * 1650.
         Kz_zth = Kz_zth - KD*zth/sqrt(zth^2.+DD^2.)
      ENDIF 

      zvis = zth
      sigusevis = K*zth/sqrt(zth^2.+D^2.) / (2.0*!PI*G1) / 1000^2.
      siguseviserr = sigusevis*0.01
   ENDIF
   sigvismin = sigusevis - siguseviserr
ENDIF

;Read in the actual data:
IF (plotwhich EQ 'kg') THEN BEGIN
   nufile = '/home/ast/user/jread/Data/Local_dm/density_j_mcs.txt'
   sigfile = '/home/ast/user/jread/Data/Local_dm/vdisp_j_mcs.txt'
   
   sigs = [49.4,4.6]
   rhodm = [0.01,0.005]
   silv = [0.025,0.013,0.014]

   ysigmin = 17
   ysigmax = 45
   ynumin = 0.04
   ynumax = 1.5
   Msigmax = 140
ENDIF

IF (plotwhich EQ 'sim') THEN BEGIN  
   ;Read in simulation data here: 
   silv = [-10,0,0]
   readcol,rhofile,rhodm,rhos,sigs
   readcol,surfsfile,zsurf,sigs_surf,sigs_surferr
   readcol,surfdmfile,zsurfdm,sigdm_surf,sigdm_surferr

   ysigmin = 15
   ysigmax = 27
   ynumin = 0.2
   ynumax = 1.8
   Msigmax = 80
ENDIF

IF (plotwhich EQ 'simple') THEN BEGIN
   ysigmin = 0
   ysigmax = 55
   ynumin = 1d-3
   ynumax = 1d0
   Msigmax = 90
ENDIF


readcol,nufile,z_bin,nu_bin,nu_bin_err
readcol,sigfile,z_bin,sig_bin,sig_bin_err
IF (n_elements(z_dat) EQ 0) THEN z_dat = z_bin

;Renormalise density to one:
normden = max(nu_bin)
nu_bin = nu_bin / normden
nu_bin_err = nu_bin_err / normden

IF (usetpar EQ 'yes') THEN $
   readcol,tfile,z_bin,tilt_bin,tilt_bin_err

;Set burnin: 
burnin = burninfac * niter

;Extract parameters:
chimods = metarray(npars,burnin:niter-1)
kzpars = metarray(nt_nu:nt_nu+nt_kz-1,burnin:niter-1)
nupars = metarray(0:nt_nu-1,burnin:niter-1)
norm1 = metarray(nt_nu+nt_kz,burnin:niter-1)
nn = niter-burnin
nn = nn(0)

IF (usetpar EQ 'yes') THEN $
   tpars = tpars(*,burnin:niter-1)

;Set up bprior:
IF (bprior EQ 'yes') THEN BEGIN
   bsurf = interpol(sigvismin,zvis,zp_kz)
   blow = bsurf * (2.*!PI*G1) * 1000.^2.
ENDIF ELSE blow = dblarr(n_elements(zp_kz)) + 1d-5

;Set the z array:
zpl = zp_kz
zpnts = n_elements(zpl)

;Calcualte averaged profiles:
kz_z = dblarr(zpnts,nn)
nu_z = dblarr(zpnts,nn)
sig_z = dblarr(zpnts,nn)
Msigma_z = dblarr(zpnts,nn)
dm_z = dblarr(zpnts,nn)

IF (usetpar EQ 'yes') THEN tilt_z = dblarr(zpnts,nn)

;Set probability boundaries to store:
p = [med-int1/2.,med-int2/2.,med,med+int2/2.,med+int1/2.]
pstyle = [2,1,0,1,2]
np = n_elements(p)

kz_zp = dblarr(np,zpnts)
nu_zp = dblarr(np,zpnts)
sig_zp = dblarr(np,zpnts)
Msigma_zp = dblarr(np,zpnts)
dm_zp = dblarr(np,zpnts)

IF (usetpar EQ 'yes') THEN tilt_zp = dblarr(np,zpnts)

;Arrange data:
pp = 0L
FOR i=0L,nn-1 DO BEGIN
   kz_z(*,pp) = kz(zpl,zp_kz,$
                   kzpars(*,i),blow,quadratic)
   nu_z(*,pp) = nu(zpl,zp_kz,nupars(*,i),quadratic,monotonic)
   sig_z(*,pp) = sigma_z(zpl,zp_kz,$
                         kzpars(*,i),blow,$
                         nupars(*,i),norm1(i),$
                         tpars(0:ntpars-1,i),tpars(ntpars:ntpars+2,i),$
                         quadratic,monotonic)

   IF (usetpar EQ 'yes') THEN $
      tilt_z(*,pp) = sigma_rz(zpl,zp_kz,tpars(0:ntpars-1,i))
   
   ;Convert kz_z to Sigma in Msun / kpc^2:
   Msigma_z(*,pp) = abs(kz_z(*,pp)) / (2.*!PI*G1) / 1000^2.

   kzparsu = abs(kzpars(*,i))
   denarr = dblarr(n_elements(kzparsu))
   denarr(0) = kzparsu(0)
   FOR jj=1L,n_elements(kzparsu)-1 DO $
      denarr(jj) = denarr(jj-1) + kzparsu(jj)
   dm_z(*,pp) = reverse(denarr) / (4.0*!PI*G1) / 1000^3.

   pp = pp + 1 
ENDFOR

;Calculate prob. ranges for each bin: 
FOR i=0L,zpnts-1 DO BEGIN 
   indexMs = sort(Msigma_z(i,0:nn-1))
   indexnu = sort(nu_z(i,0:nn-1))
   indexsg = sort(sig_z(i,0:nn-1))
   indexkz = sort(kz_z(i,0:nn-1))
   indexdm = sort(dm_z(i,0:nn-1))

   IF (usetpar EQ 'yes') THEN indext = sort(tilt_z(i,0:nn-1))

   FOR j=0L,np-1 DO BEGIN
      Msigma_zp(j,i) = Msigma_z(i,indexMs(p(j)*(nn-1)))
      nu_zp(j,i) = nu_z(i,indexnu(p(j)*(nn-1)))
      sig_zp(j,i) = sig_z(i,indexsg(p(j)*(nn-1)))
      kz_zp(j,i) = kz_z(i,indexkz(p(j)*(nn-1)))
      dm_zp(j,i) = dm_z(i,indexdm(p(j)*(nn-1)))

      IF (usetpar EQ 'yes') THEN tilt_zp(j,i) = tilt_z(i,indext(p(j)*(nn-1)))
   ENDFOR
ENDFOR

;Array dump the processed data: 
zplt = transpose(zpl)
arraydump,basename+'_msigma_zp.txt',[[zplt],[Msigma_zp(0,*)],[Msigma_zp(1,*)],[Msigma_zp(2,*)],$
                                     [Msigma_zp(3,*)],[Msigma_zp(4,*)]],6
arraydump,basename+'_nu_zp.txt',[[zplt],[nu_zp(0,*)],[nu_zp(1,*)],[nu_zp(2,*)],$
                                 [nu_zp(3,*)],[nu_zp(4,*)]],6
arraydump,basename+'_sig_zp.txt',[[zplt],[sig_zp(0,*)],[sig_zp(1,*)],[sig_zp(2,*)],$
                                  [sig_zp(3,*)],[sig_zp(4,*)]],6
arraydump,basename+'_kz_zp.txt',[[zplt],[kz_zp(0,*)],[kz_zp(1,*)],[kz_zp(2,*)],$
                                 [kz_zp(3,*)],[kz_zp(4,*)]],6
IF (usetpar EQ 'yes') THEN $
   arraydump,basename+'_tilt_zp.txt',[[zplt],[tilt_zp(0,*)],[tilt_zp(1,*)],[tilt_zp(2,*)],$
                                      [tilt_zp(3,*)],[tilt_zp(4,*)]],6
;Do the plots: 
oname = basename + '_likelihood' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

likeli = abs(metarray(npars,*))
likeli = likeli(where(likeli NE 0))
likeli = likeli(where(likeli EQ likeli))
xarr = dindgen(n_elements(likeli))
plot,xarr/1.e4,likeli,/ylog,yrange=[min(likeli),max(likeli)],title='!6',$
     xtitle='iterations / 10!u4!n',ytitle='abs(Log[Likelihood])'

;Overplot chosen burnin: 
oplot,[burnin,burnin]/1.e4,[1d-300,1d300],color=2

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

oname = basename + '_msigmaz' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

maxzdat = max(z_dat) 
FOR i=0L,np-1 DO BEGIN
   IF (i EQ 0) THEN BEGIN 
      plot,zpl,Msigma_zp(i,*),title='!6',xtitle='z(kpc)',$
           ytitle='!7R!6!dz!n [Msun pc!u-2!n]',$
           yrange=[0,Msigmax],xrange=[0,maxzdat],$
           linestyle = pstyle(i)

      IF (plotwhich EQ 'sim') THEN BEGIN
         oplot,zsurf,sigs_surf,color=2
         oplot,zsurf,sigdm_surf,color=3
         oplot,zsurf,sigs_surf+sigdm_surf,color=4
      ENDIF
      IF (plotwhich EQ 'simple') THEN BEGIN 
         oplot,zth,abs(Kz_zth) / (2.*!PI*G1) / 1000^2.,color=4
         oplot,zvis,sigusevis,color=2
      ENDIF
   ENDIF ELSE $
      oplot,zpl,Msigma_zp(i,*),linestyle=pstyle(i)
ENDFOR

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

oname = basename + '_sigmaz' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

FOR i=0L,np-1 DO BEGIN
   IF (i EQ 0) THEN BEGIN
      plot,zpl,sig_zp(i,*),title='!6',xtitle='z(kpc)',$
           ytitle='!7r!6!dz!n [km/s]',$
           yrange=[ysigmin,ysigmax],xrange=[0,max(z_dat)],$
           linestyle = pstyle(i)

      ;Overplot data:
      oploterror,z_bin,sig_bin,sig_bin_err,color=2,errcol=2
   ENDIF ELSE $
      oplot,zpl,sig_zp(i,*),linestyle=pstyle(i)
ENDFOR

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

oname = basename + '_nuz' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

jl=0L
while (zpl(jl) LT min(z_bin)) do jl=jl+1
FOR i=0L,np-1 DO BEGIN
   pnorm = nu_zp(i,jl)
   IF (i EQ 0) THEN BEGIN
      plot,zpl,nu_zp(i,*)/pnorm,title='!6',xtitle='z(kpc)',/ylog,$
           ytitle='!7m!6 [units]',yrange=[ynumin,ynumax],$
           xrange=[0,max(z_dat)],$
           linestyle = pstyle(i)
      
      ;Overplot data:
      oploterror,z_bin,nu_bin,nu_bin_err,color=2,errcol=2
   ENDIF ELSE $
      oplot,zpl,nu_zp(i,*)/pnorm,linestyle=pstyle(i)
ENDFOR

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

oname = basename + '_dmprof' + filetype
IF (devtype EQ 'PS') THEN $
   device,filename=oname,/color,$
          xsize=16,ysize=16,/encapsul

FOR i=0L,np-1 DO BEGIN
   IF (i EQ 0) THEN BEGIN
      plot,zp_kz,dm_zp(i,*),$
           yrange = [0,0.03],xtitle='z(kpc)',ytitle='!7q!6!ddm!n [Msun pc!u-3!n]',$
           title='!6',linestyle=pstyle(i)
   ENDIF ELSE $
      oplot,zpl,dm_zp(i,*),linestyle=pstyle(i)
ENDFOR

IF (plotwhich EQ 'simple') THEN BEGIN
   ;Overplot theory curve:
   Kz_zthd = -2.0 * F * zth
   IF (adddarkdisc EQ 'yes') THEN $
      Kz_zthd = Kz_zthd - KD*zth/sqrt(zth^2.+DD^2.)
   Sigz_zth = abs(Kz_zthd) / (2.0*!PI*G1) / 1000^2.
   denth = deriv(zth,Sigz_zth) / 1000.
   oplot,zth,denth,color=4
ENDIF 

IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
   image=TVRD(0,0,fxsize,fysize,TRUE=1)
   write_png,oname,image
ENDELSE

IF (usetpar EQ 'yes') THEN BEGIN 
   oname = basename + '_tilt' + filetype
   IF (devtype EQ 'PS') THEN $
      device,filename=oname,/color,$
             xsize=16,ysize=16,/encapsul
   
   FOR i=0L,np-1 DO BEGIN
      IF (i EQ 0) THEN BEGIN
         plot,zpl,tilt_zp(0,*),$
              title='!6',xtitle='z(kpc)',$
              ytitle='!7r!6!dRz!n [km!u2!n s!u-2!n]',$
              yrange=[0.,3250.],$
              linestyle = pstyle(i)
         
         oploterror,z_bin,tilt_bin,tilt_bin_err,color=2,errcol=2
      ENDIF ELSE oplot,zpl,tilt_zp(i,*),linestyle=pstyle(i)
   ENDFOR
   
   IF (devtype EQ 'PS') THEN device,/close ELSE BEGIN
      image=TVRD(0,0,fxsize,fysize,TRUE=1)
      write_png,oname,image
   ENDELSE
ENDIF 

stop

;close the output file, set the plotting back to the screen and exit.
IF devtype EQ 'PS' THEN device,/close	
tidy

END