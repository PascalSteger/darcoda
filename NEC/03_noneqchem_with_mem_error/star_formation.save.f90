!################################################################
!################################################################
!################################################################
!################################################################
subroutine star_formation(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH ,neq_spec,dust
  use random
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine spawns star-particle of constant mass
  ! using a Poisson probability law if some gas condition are fulfilled. 
  ! It modifies hydrodynamic variables according to mass conservation 
  ! and assumes an isothermal transformation... 
  ! On exit, the gas velocity and sound speed are unchanged.
  ! New star particles are synchronized with other collisionless particles.
  ! Array flag2 is used as temporary work space.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::t0,d0,e0
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ivar,ngrid,icpu,index_star,ndebris_tot
  integer ::igrid,ix,iy,iz,ind,i,j,n,iskip,istar,inew,nx_loc
  integer ::ntot,ntot_all,info,nstar_corrected,ideb,ndeb,ii,nn
#ifdef SOLVERhydro
  integer ::imetal=6
#endif
#ifdef SOLVERmhd
  integer ::imetal=9
#endif
  logical ::ok_free,ok_all
  real(dp)::d,x,y,z,u,v,w,e,zg,vdisp,dgas
  real(dp)::mstar,dstar,tstar,nISM,nCOM
  real(dp)::velc,uc,vc,wc
  real(dp)::vxgauss,vygauss,vzgauss,birth_epoch
  real(kind=8)::mlost,mtot,mlost_all,mtot_all
  real(kind=8)::RandNum,GaussNum,PoissMean   
  real(dp)::vsn,costheta,sintheta,phi,cosphi,sinphi,twopi
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min
  real(dp)::mdebris,vdebris,zdebris,rdebris
  real(dp)::bx1,bx2,by1,by2,bz1,bz2

  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,nstar
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  integer ,dimension(1:nvector),save::list_debris,ind_debris1,ind_debris2
  logical ,dimension(1:nvector),save::ok,ok_new=.true.,ok_true=.true.
  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all
  
  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)' Entering star_formation'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim
  twopi=2.0*ACOS(-1.0D0)
  ! Bubble radius cannot be smaller than 2 fine cells
  rdebris=MAX(2.0*dx_min*scale_l/aexp,rbubble*3.08d18)

  ! Set up supernovae parameters
  if(eta_sn>0)then
     ! Compute debris metallicity
     zdebris=yield/(1d0+f_w)
     ! Supernovae debris velocity in cgs
     vsn=sqrt(2.*1d51/(10.*2d33))/sqrt(1d0+f_w)
     ! Compute debris flight time in Myr
     t_delay=(rdebris/vsn) / (1d6*365.*24.*3600.)
     ! Compute debris velocity in code units
     vdebris=vsn/scale_v
  endif

  ! Star formation time scale from Gyr to code units
  ! SFR apply here for long lived stars only
  t0=t_star*(1d9*365.*24.*3600.)/scale_t

  ! ISM density threshold from H/cc to code units
  nISM = n_star
  nCOM = del_star*omega_b*rhoc/aexp**3*XH/mH
  nISM = MAX(nCOM,nISM)
  d0   = nISM/scale_nH

  ! ISM typical temperature (T/mu) from K to code units
  e0=T2_star/scale_T2/(gamma-1.0)

  ! Star particle mass
  mstar=MAX(del_star*omega_b*rhoc*XH/mH,n_star) &
       & /(scale_nH*aexp**3)*vol_min/(1.0+eta_sn*(1d0+f_w))
  dstar=mstar/vol_loc

  ! SN debris particle mass
  mdebris=mstar*eta_sn*(1d0+f_w)/dble(2*ndebris)

  ! Birth epoch
  birth_epoch=t

  ! Cells center position relative to grid center position
  do ind=1,twotondim  
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

#if NDIM==3
  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim  
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           e=uold(ind_cell(i),5)/d
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e-0.5d0*(u**2+v**2+w**2)
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=u
           uold(ind_cell(i),3)=v
           uold(ind_cell(i),4)=w
           uold(ind_cell(i),5)=e
        end do
        if(metal)then
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),imetal)/d
              uold(ind_cell(i),imetal)=w
           end do
        endif
     end do
  end do

  !------------------------------------------------
  ! Compute number of new stars in each cell
  !------------------------------------------------
  ntot=0
  ndebris_tot=0
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Star formation criterion ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           if(d<=d0)ok(i)=.false. ! Density criterion
        end do
        ! Calculate number of new stars in each cell using Poisson statistics
        do i=1,ngrid
           nstar(i)=0
           if(ok(i))then
              ! Compute mean number of events
              d=uold(ind_cell(i),1)
              tstar=t0*sqrt(d0/d)
              PoissMean=dtnew(ilevel)/tstar*d/dstar
              ! Compute Poisson realisation
              call poissdev(localseed,PoissMean,nstar(i))
              ! Compute depleted gas
              dgas=nstar(i)*dstar*(1.0+eta_sn*(1d0+f_w))
              ! Security to prevent more than 90% of gas depletion
              if (dgas > 0.9*d) then
                 nstar_corrected=int(0.9*d/(dstar*(1.0+eta_sn*(1d0+f_w))))
                 mstar_lost=mstar_lost+(nstar(i)-nstar_corrected)*mstar
                 nstar(i)=nstar_corrected
              endif
              ! Compute new stars local statistics
              mstar_tot=mstar_tot+nstar(i)*mstar
              if(nstar(i)>0)ntot=ntot+1
              if(eta_sn>0.0d0)ndebris_tot=ndebris_tot+2*ndebris*nstar(i)
           endif
        enddo
        ! Store nstar in array flag2
        do i=1,ngrid
           flag2(ind_cell(i))=nstar(i)
        end do
     end do
  end do

  !---------------------------------
  ! Check for free particle memory
  !---------------------------------
  ok_free=(numbp_free-ntot-ndebris_tot)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if
  
  !---------------------------------
  ! Compute global stars statistics
  !---------------------------------
#ifndef WITHOUTMPI
  mlost=mstar_lost; mtot=mstar_tot
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mtot,mtot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mlost,mlost_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
  mtot_all=mstar_tot
  mlost_all=mstar_lost
#endif
  ntot_star_cpu=0; ntot_star_all=0
  ntot_star_cpu(myid)=ntot
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot_star_cpu,ntot_star_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_star_cpu(1)=ntot_star_all(1)
#endif
  do icpu=2,ncpu
     ntot_star_cpu(icpu)=ntot_star_cpu(icpu-1)+ntot_star_all(icpu)
  end do
  nstar_tot=nstar_tot+ntot_all
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" Level = ",I6," New star = ",I6," Tot =",I8," Mass =",1PE9.3," Lost =",1Pe9.3,"%")')&
             & ilevel,ntot_all,nstar_tot,mtot_all,mlost_all/mtot_all*100.
     endif
  end if

  !------------------------------
  ! Create new star particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_star=nstar_tot-ntot_all
  else
     index_star=nstar_tot-ntot_all+ntot_star_cpu(myid-1)
  end if

  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Flag cells with at least one new star
        do i=1,ngrid
           ok(i)=flag2(ind_cell(i))>0
        end do
        
        ! Gather new star arrays
        nnew=0
        do i=1,ngrid
           if (ok(i))then
              nnew=nnew+1
              ind_grid_new(nnew)=ind_grid(i)
              ind_cell_new(nnew)=ind_cell(i)
           end if
        end do
        
        ! Update linked list
        call remove_free(ind_part,nnew)
        call add_list(ind_part,ind_grid_new,ok_new,nnew)
        
        ! Calculate new star particle and modify gas density
        do i=1,nnew
           index_star=index_star+1

           ! Get gas variables
           n=flag2(ind_cell_new(i))
           d=uold(ind_cell_new(i),1)
           u=uold(ind_cell_new(i),2)
           v=uold(ind_cell_new(i),3)
           w=uold(ind_cell_new(i),4)
           x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
           y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale
           if(metal)zg=uold(ind_cell_new(i),imetal)

           ! Set new star particle variables
           tp(ind_part(i))=birth_epoch  ! Birth epoch
           mp(ind_part(i))=n*mstar      ! Mass
           levelp(ind_part(i))=ilevel   ! Level
           idp(ind_part(i))=index_star  ! Identity
           xp(ind_part(i),1)=x
           xp(ind_part(i),2)=y
           xp(ind_part(i),3)=z
           vp(ind_part(i),1)=u
           vp(ind_part(i),2)=v
           vp(ind_part(i),3)=w
           if(metal)zp(ind_part(i))=zg  ! Initial star metallicity
           
           if(eta_sn>0.0)then

              ! Loop over debris by vector sweeps
              ndebris_tot=n*ndebris
              do ideb=1,ndebris_tot,nvector
                 ndeb=MIN(nvector,ndebris_tot-ideb+1)
                 
                 ! Get ndebris twin particles
                 do j=1,ndeb
                    list_debris(j)=ind_grid_new(i)
                 end do
                 call remove_free(ind_debris1,ndeb)
                 call add_list(ind_debris1,list_debris,ok_true,ndeb)
                 call remove_free(ind_debris2,ndeb)
                 call add_list(ind_debris2,list_debris,ok_true,ndeb)
                 
                 ! Set debris twin particle variables
                 do j=1,ndeb
                    
                    ! Monte-Carlo Sedov solution
                    call ranf(localseed,RandNum)
                    velc=vdebris*RandNum**(1./3.)*f_ek
                    call ranf(localseed,RandNum)
                    costheta=2.0*RandNum-1.0
                    sintheta=sqrt(1.0-costheta**2)
                    call ranf(localseed,RandNum)
                    cosphi=cos(twopi*RandNum)
                    sinphi=sin(twopi*RandNum)
                    uc=velc*cosphi*sintheta
                    vc=velc*sinphi*sintheta
                    wc=velc*costheta
                    
                    ! First debris
                    tp(ind_debris1(j))=birth_epoch  ! Birth epoch
                    mp(ind_debris1(j))=mdebris      ! Mass
                    levelp(ind_debris1(j))=ilevel   ! Level
                    idp(ind_debris1(j))=0           ! Identity
                    xp(ind_debris1(j),1)=x
                    xp(ind_debris1(j),2)=y
                    xp(ind_debris1(j),3)=z
                    vp(ind_debris1(j),1)=u + uc
                    vp(ind_debris1(j),2)=v + vc
                    vp(ind_debris1(j),3)=w + wc
                    if(metal)zp(ind_debris1(j))=zg !+ zdebris ! Commented out by ACB. Added again below
                    
                    ! Second debris
                    tp(ind_debris2(j))=birth_epoch  ! Birth epoch
                    mp(ind_debris2(j))=mdebris      ! Mass
                    levelp(ind_debris2(j))=ilevel   ! Level
                    idp(ind_debris2(j))=0           ! Identity
                    xp(ind_debris2(j),1)=x
                    xp(ind_debris2(j),2)=y
                    xp(ind_debris2(j),3)=z
                    vp(ind_debris2(j),1)=u - uc
                    vp(ind_debris2(j),2)=v - vc
                    vp(ind_debris2(j),3)=w - wc
                    if(metal)zp(ind_debris2(j))=zg !+ zdebris ! Commented out by ACB. Added again below                
                    
                 end do
              end do
              ! End loop over debris
           endif

        end do
        ! End loop over new star particles

        ! Modify gas density according to mass depletion
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           uold(ind_cell(i),1)=uold(ind_cell(i),1) - &
                & flag2(ind_cell(i))*dstar*(1.0+eta_sn*(1d0+f_w))
#ifdef NONEQCHEM
           !ACB make sure to adjust chemical species, too.
           d=uold(ind_cell(i),1)/d
           nn=neq_spec
           if (.not.dust)nn=nn-1
           do ii=1,nn
              uold(ind_cell(i),6+ii)=uold(ind_cell(i),6+ii)*d
           enddo
#endif
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids
  
  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim  
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)
           v=uold(ind_cell(i),3)
           w=uold(ind_cell(i),4)
           e=uold(ind_cell(i),5)
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e+0.5d0*(u**2+v**2+w**2)
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=d*u
           uold(ind_cell(i),3)=d*v
           uold(ind_cell(i),4)=d*w
           uold(ind_cell(i),5)=d*e
        end do
        if(metal)then
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),imetal)
              uold(ind_cell(i),imetal)=d*w
           end do
        endif
     end do
  end do

#endif

end subroutine star_formation 
!################################################################
!################################################################
!################################################################
!################################################################
subroutine feedback(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel,ii
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by exploding star particles, 
  ! called here debris particles, after a time delay t_delay.
  !------------------------------------------------------------------------
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0,scale,dx_min,vsn,vdebris
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

#if NDIM==3
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=(0.5D0**nlevelmax)*scale

  ! Supernovae debris velocity in cgs
  vsn=sqrt(2.0*1d51/(10.*2d33))/sqrt(1d0+f_w)
  ! Compute debris flight time in Myr
  t_delay=1.0d1  !(2.0*dx_min*scale_l/aexp/vsn) / (1d6*365.*24.*3600.)    !OSCAR 

  ! Time delay from Myr to code units
  t0=t_delay*1d6*(365.*24.*3600.)/scale_t
  ! Compute debris velocity in code units
  vdebris=vsn/scale_v
  ! Gather debris particles only.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count old enough debris particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).eq.0.and.tp(ipart).lt.(t-t0))then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather old enough debris particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ! Select only debris particles
              if(idp(ipart).eq.0.and.tp(ipart).lt.(t-t0))then
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call sn2(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,vdebris)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call sn2(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,vdebris)
  end do 
  ! End loop over cpus

#endif

111 format('   Entering feedback for level ',I2)

end subroutine feedback
!################################################################
!################################################################
!################################################################
!################################################################ 
subroutine sn2(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,vdebris)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, only:neq_spec,dust
  implicit none
  integer::ng,np,ilevel,ii
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine feedback. Each debris particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! uold.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,nn
  real(dp)::xxx,mmm,ethermal,vdebris,zdebris
  real(dp)::dx,dx_loc,scale,vol_loc,d
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok,ok_true=.true.
  real(dp),dimension(1:nvector),save::meff
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  zdebris=yield/(1d0+f_w)

#if NDIM==3
  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in sn2'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=x(j,idim)
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! If not, rescale position at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           x(j,idim)=x(j,idim)/2.0D0
        end if
     end do
  end do
  ! If not, redo NGP at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           id(j,idim)=x(j,idim)
        end if
     end do
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        else
           icd(j,idim)=id(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     else
        icell(j)=1+icd(j,1)+3*icd(j,2)+9*icd(j,3)   
     end if
  end do
        
  ! Compute parent cell adresses and particle effective mass
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
        meff(j)=mp(ind_part(j))/vol_loc 
     else
        indp(j)=nbors_father_cells(ind_grid_part(j),icell(j))
        meff(j)=mp(ind_part(j))/vol_loc/dble(twotondim)
     end if
  end do

  ! Update hydro variables due to feedback 
  do j=1,np
     ethermal=0.5d0*meff(j)*vdebris**2.0 
     d=uold(indp(j),1)
     uold(indp(j),1)=uold(indp(j),1)+meff(j)   !We must dump mass for conservation
#ifdef NONEQCHEM
           !ACB make sure to adjust chemical species, too.
           d=uold(indp(j),1)/d
           nn=neq_spec
           if (.not.dust)nn=nn-1
           do ii=1,nn
              uold(indp(j),6+ii)=uold(indp(j),6+ii)*d
           enddo
#endif
     uold(indp(j),2)=uold(indp(j),2)+meff(j)*vp(ind_part(j),1)  
     uold(indp(j),3)=uold(indp(j),3)+meff(j)*vp(ind_part(j),2)   
     uold(indp(j),4)=uold(indp(j),4)+meff(j)*vp(ind_part(j),3)
     uold(indp(j),5)=uold(indp(j),5)+0.5d0*meff(j)*(vp(ind_part(j),1)**2+&
          & vp(ind_part(j),2)**2+vp(ind_part(j),3)**2)+ethermal
  end do
  if(metal)then
     do j=1,np
        uold(indp(j),6)=uold(indp(j),6)+meff(j)*(zp(ind_part(j))+zdebris) ! always same debris
     end do
  endif

  ! Parent particle linked list
  do j=1,np
     list1(j)=ind_grid(ind_grid_part(j))
  end do 

  ! Remove debris particle
  call remove_list(ind_part,list1,ok_true,np)
  call add_free_cond(ind_part,ok_true,np)

#endif
  
end subroutine sn2
!################################################################
!################################################################
!################################################################
!################################################################
