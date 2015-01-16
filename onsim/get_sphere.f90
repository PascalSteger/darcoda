program get_sphere
  !-------------------------------------------------------------------
  ! extract particles (DM/stars only) to ASCII
  !-------------------------------------------------------------------
  implicit none
  integer::ncpu,ndim,npart,i,j,k,icpu,ipos,nstar,iii
  integer::ncpu2,npart2,ndim2,levelmin,levelmax
  integer::nx=0,ny=0,ncpu_read,n_frw
  real(kind=8)::mtot,time,time_tot,time_simu,weight,boxlen
  real(kind=8)::xc,yc,zc,r, drx, dry, drz, dd2
  integer::npart_actual,nsink
  real(kind=8)::massconv,vscale
  real(kind=8)::aexp,t,omega_m,omega_l,omega_b,omega_k,h0,unit_l,unit_t,unit_d
  real(kind=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  real(kind=8),dimension(:,:),allocatable::map
  real(kind=4),dimension(:,:),allocatable::x,v
  real(kind=4),dimension(:)  ,allocatable::m,age,star_metal,dum
  integer(kind=4),dimension(:),allocatable::id_dm
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=128)::nomfich,repository,filetype='bin'
  character::A13*13,A14*14
  logical::ok,ok_part,periodic=.false.,star=.false.,ageweight=.true.
  integer::impi
  real(kind=8)::pi
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list,id_star

  pi = acos(-1d0)
  call read_params

  !-----------------------------------------------
  ! Lecture du fichier particules au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  open(unit=10,file=nomfich,form='formatted',status='old')
  read(10,'(a13,I11)')A13,ncpu
  read(10,'(a13,I11)')A13,ndim
  read(10,'(a13,I11)')A13,levelmin
  read(10,'(a13,I11)')A13,levelmax
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,'(a13,E23.15)')A13,boxlen
  read(10,'(a13,E23.15)')A13,t
  read(10,'(a13,E23.15)')A13,aexp
  read(10,'(a13,E23.15)')A13,h0
  read(10,'(a13,E23.15)')A13,omega_m
  read(10,'(a13,E23.15)')A13,omega_l
  read(10,'(a13,E23.15)')A13,omega_k
  read(10,'(a13,E23.15)')A13,omega_b
  read(10,'(a13,E23.15)')A13,unit_l
  read(10,'(a13,E23.15)')A13,unit_d
  read(10,'(a13,E23.15)')A13,unit_t
  read(10,*)

  massconv=(h0/3.08d24*1e5)**2*3./8./pi/6.67d-8*omega_m*(boxlen*3.08d24)**3/2d33
  vscale=unit_l/unit_t

  read(10,'(a14,A80)')A14,ordering
  !write(*,'(a,a14,A20)')"#",A14,TRIM(ordering)
  read(10,*)
  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
     allocate(bound_key(0:ncpu))
     allocate(cpu_read(1:ncpu))
     cpu_read=.false.
     do impi=1,ncpu
        read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
     end do
  endif
  close(10)

  !-----------------------
  ! Cosmological model
  !-----------------------
  ! Allocate look-up tables
  n_frw=1000
  allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
  allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
  
  ! Compute Friedman model look up table
  !write(*,*)'#Computing Friedman model'
  call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
       & 1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)
  
  ! Find neighboring expansion factors
  i=1
  do while(aexp_frw(i)>aexp.and.i<n_frw)
     i=i+1
  end do
  ! Interpolate time
  time_simu=t_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
       & t_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
  !write(*,*)'#Age simu=',(time_tot+time_simu)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)

  !-----------------------
  ! Map parameters
  !-----------------------
  if(nx==0)then
     nx=2**levelmin
  endif
  if(ny==0)then
     ny=nx
  end if
  !write(*,*)'#time=',t
  !write(*,*)'#Working map =',nx,ny
  allocate(map(0:nx,0:ny))
  map=0.0d0
  ncpu_read=ncpu
  do j=1,ncpu
     cpu_list(j)=j
  end do

  npart=0
  do k=1,ncpu!ncpu_read
     icpu=k!cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     !print *, "# FILE ",nomfich
     open(unit=1,file=nomfich,status='old',form='unformatted')
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)nstar
     close(1)
     npart=npart+npart2
  end do
  !write(*,*)'#Found ',npart,' particles.'
  if(nstar>0)then
     if(star)then
        !write(*,*)'#Keeping star particles.'
     else
        !write(*,*)'#Discard star particles.'
     endif
  endif

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  npart_actual=0
  mtot=0.0d0
  do k=1,ncpu!ncpu_read
     icpu=k!cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted',ACCESS='SEQUENTIAL')
     !write(*,*)'#Processing file '//TRIM(nomfich)
     read(1)ncpu2
     !write(*,*)'#ncpu2 = ',ncpu2
     read(1)ndim2
     !write(*,*)'#ndim2 = ',ndim2
     read(1)npart2
     !write(*,*)'#npart2 = ',npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)nsink
     !write(*,*)'#nsink = ',nsink

     allocate(m(1:npart2))
     allocate(dum(1:npart2))
     allocate(id_dm(1:npart2))
     if(nstar>0)then
        !write(*,*)'#we have stars!'
        allocate(age(1:npart2))
        allocate(star_metal(1:npart2))
        allocate(id_star(1:npart2))
     endif
     allocate(x(1:npart2,1:ndim2))
     allocate(v(1:npart2,1:ndim2))
     ! Read position
     do i=1,ndim
        !write(*,*)'# read x dim',i
        read(1)m
        x(1:npart2,i)=m
     end do
     ! read velocity
     do i=1,ndim
        !write(*,*)'# read x dim',i
        read(1)m
        v(1:npart2,i)=m
     end do
     ! Read mass
     !write(*,*)'# read mass'
     read(1)m
     !write(*,*)'# read no dum'
     !read(1)dum
     if(nstar>0)then
        !write(*,*)'# read stars'
        read(1)age
        read(1)star_metal ! Skip identity
        read(1)id_star ! Skip identity
     endif
     !write(*,*)'# read id_dm'
     read(1)id_dm
     !write(*,*)'# success reading'
     close(1)

     
     do i=1,npart2
        weight=1.
        drx = x(i,1)-xc
        dry = x(i,2)-yc
        drz = x(i,3)-zc
        dd2 = drx**2 + dry**2 + drz**2
        ok_part=(dd2<r**2)
        
        if(nstar>0)then
           if(star)then
              ok_part=ok_part.and.(age(i).ne.0.0d0)
              if(ageweight)then
                 iii=1
                 do while(tau_frw(iii)>age(i).and.iii<n_frw)
                    iii=iii+1
                 end do
                 ! Interploate time
                 time=t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                      & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                 time=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                 weight=1.
                 if(time>0.01)then
                    weight=(time/0.01)**(-0.7)
                 endif
                 age(i)=time!*weight
              endif
           else
              ok_part=ok_part.and.(age(i).eq.0.0d0)
           endif
        endif
        
        if(ok_part)then
           npart_actual=npart_actual+1
           if(star)then 
              print "(6(1pe15.8,1X),I8,1X,3(1pe15.8,1X))", m(i)*unit_d*unit_l**3/2d33,x(i,1),x(i,2),x(i,3),star_metal(i),age(i),&
                   id_star(i),v(i,1)*vscale,v(i,2)*vscale,&
                   v(i,3)*vscale
           else
              !print "(I8,1X)", id_dm(i)
              print "(4(1pe15.8,1X),I8,1X,3(1pe15.8,1X))", m(i)*unit_d*unit_l**3/2d33,x(i,1),x(i,2),x(i,3),id_dm(i), &
                   v(i,1)*vscale,v(i,2)*vscale,v(i,3)*vscale

           endif
           mtot=mtot+m(i)
        end if
     end do
     deallocate(x,m,v,id_dm,dum)
     if(nstar>0) deallocate(age,star_metal,id_star)
  end do
  !write(*,*)'#Total mass=',mtot
  !if(.not. star)write(*,*)'#npart tot=',npart_actual
  
contains

  subroutine read_params

      implicit none

      integer       :: i,n
      integer       :: iargc
      character(len=4)   :: opt
      character(len=128) :: arg
      
      n = iargc()
      if (n < 4) then
         print *, 'usage: get_sphere -inp input_dir'
         print *, '                 [-xc  xcenter [Mpc/h]] '
         print *, '                 [-yc  ycenter [Mpc/h]] '
         print *, '                 [-zc  zcenter [Mpc/h]] '
         print *, '                 [-r   radius  [Mpc/h]] '
         print *, '                 [-nx  nx  ] '
         print *, '                 [-ny  ny  ] '
         print *, '                 [-per flag] '
         print *, '                 [-str flag] '
         print *, '                 [-fil filetype] '
         print *, 'ex: get_sphere -inp ../run/outpu_00050'// &
              &   '             -xc 0.5 -yc 0.5 -zc 0.5 -r 0.1 > sample.dat'
         stop
      end if

      do i = 1,n,2
         call getarg(i,opt)
         if (i == n) then
            print '("option ",a2," has no argument")', opt
            stop 2
         end if
         call getarg(i+1,arg)
         select case (opt)
         case ('-inp')
            repository = trim(arg)
         case ('-xc')
            read (arg,*) xc
         case ('-yc')
            read (arg,*) yc
         case ('-zc')
            read (arg,*) zc
         case ('-r')
            read (arg,*) r
         case ('-nx')
            read (arg,*) nx
         case ('-ny')
            read (arg,*) ny
         case ('-per')
            read (arg,*) periodic
         case ('-str')
            read (arg,*) star
         case ('-age')
            read (arg,*) ageweight
         case ('-fil')
            filetype = trim(arg)
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

      return

    end subroutine read_params

  end program get_sphere

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif

end subroutine title

!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(kind=8),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo

     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do

     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert3d
!================================================================
!================================================================
!================================================================
!================================================================
subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
     & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)

  implicit none
  integer::ntable
  real(kind=8)::O_mat_0, O_vac_0, O_k_0
  real(kind=8)::alpha,axp_min,age_tot
  real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  real(kind=8)::axp_tau, axp_t
  real(kind=8)::axp_tau_pre, axp_t_pre
  real(kind=8)::dadtau, dadt
  real(kind=8)::dtau,dt
  real(kind=8)::tau,t
  integer::nstep,nout,nskip

!  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
!     write(*,*)'Error: non-physical cosmological constants'
!     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
!     write(*,*)'The sum must be equal to 1.0, but '
!     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
!     stop
!  end if

  axp_tau = 1.0D0
  axp_t = 1.0D0
  tau = 0.0D0
  t = 0.0D0
  nstep = 0
  
  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau
     
     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
  end do

  age_tot=-t
  !write(*,666)-t
  !666 format('#Age of the Universe (in unit of 1/H0)=',1pe10.3)

  nskip=nstep/ntable
  
  axp_t = 1.d0
  t = 0.d0
  axp_tau = 1.d0
  tau = 0.d0
  nstep = 0
  nout=0
  t_out(nout)=t
  tau_out(nout)=tau
  axp_out(nout)=axp_tau
  hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau

     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
     if(mod(nstep,nskip)==0)then
        nout=nout+1
        t_out(nout)=t
        tau_out(nout)=tau
        axp_out(nout)=axp_tau
        hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
     end if

  end do
  t_out(ntable)=t
  tau_out(ntable)=tau
  axp_out(ntable)=axp_tau
  hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

end subroutine friedman

function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
  real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
  dadtau = axp_tau*axp_tau*axp_tau *  &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
       &     O_k_0   * axp_tau )
  dadtau = sqrt(dadtau)
  return
end function dadtau

function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
  real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
  dadt   = (1.0D0/axp_t)* &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_t*axp_t*axp_t + &
       &     O_k_0   * axp_t )
  dadt = sqrt(dadt)
  return
end function dadt
