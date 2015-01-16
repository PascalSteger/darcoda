program part2map
  !--------------------------------------------------------------------------
  ! dark matter density plot
  !--------------------------------------------------------------------------
  implicit none
  integer::ncpu,ndim,npart,i,j,k,icpu,ipos
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel
  integer::nx=0,ny=0,ix,iy,ixp1,iyp1,idim,jdim,ncpu_read
  real*8::mtot,ddx,ddy,dex,dey,t
  real*8::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real*8::xc=0.5,yc=0.5,zc=0.5,rc=0.5
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin,nstar
  real*8::xxmin,xxmax,yymin,yymax,dx,dy,deltax
  real*4,dimension(:,:),allocatable::toto
  real*8,dimension(:,:),allocatable::map

#ifndef NPRE
  integer,parameter::dp=kind(1.0E0) ! default: real*4
#else
#if NPRE==4
  integer,parameter::dp=kind(1.0E0) ! real*4
#else
  integer,parameter::dp=kind(1.0D0) ! real*8
#endif
#endif
  real(dp),dimension(:,:),allocatable::x
  real(dp),dimension(:)  ,allocatable::m,metal

  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering,junk
  character(LEN=128)::nomfich,repository,outfich
  logical::ok,ok_part,periodic=.false.,star=.false.
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real*8,dimension(1:8)::bounding_min,bounding_max
  real*8::dkey,order_min,dmax
  real*8,dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list

  call read_params

  !-----------------------------------------------
  ! Lecture du fichier particules au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/part/part_'//TRIM(nchar)//'.out00001'
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
  read(10,"(A13,I11)")junk,ncpu
  read(10,"(A13,I11)")junk,ndim
  read(10,"(A13,I11)")junk,levelmin
  read(10,"(A13,I11)")junk,levelmax
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,'(A13,E23.15)')junk,t
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,'(A13,A80)')junk,ordering
  write(*,'(A13,A20)')junk,TRIM(ordering)
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
  ! Map parameters
  !-----------------------
  if(nx==0)then
     nx=2**(levelmin)
  endif
  if(ny==0)then
     ny=nx
  end if
  write(*,*)'time=',t
  write(*,*)'Working map =',nx,ny
  allocate(map(0:nx,0:ny))
  map=0.0d0
  if (proj=='x')then
     idim=2
     jdim=3
     xxmin=ymin ; xxmax=ymax
     yymin=zmin ; yymax=zmax
  else if (proj=='y') then
     idim=1
     jdim=3
     xxmin=xmin ; xxmax=xmax
     yymin=zmin ; yymax=zmax
  else
     idim=1
     jdim=2
     xxmin=xmin ; xxmax=xmax
     yymin=ymin ; yymax=ymax
  end if
  dx=(xxmax-xxmin)/dble(nx)
  dy=(yymax-yymin)/dble(ny)

  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     do ilevel=1,levelmax
        deltax=0.5d0**ilevel
        if(deltax.lt.dmax)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xmin*dble(maxdom))
        imax=imin+1
        jmin=int(ymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zmin*dble(maxdom))
        kmax=kmin+1
     endif
     
     dkey=(dble(2**(levelmax+1)/dble(maxdom)))**ndim
     ndom=1
     if(bit_length>0)ndom=8
     idom(1)=imin; idom(2)=imax
     idom(3)=imin; idom(4)=imax
     idom(5)=imin; idom(6)=imax
     idom(7)=imin; idom(8)=imax
     jdom(1)=jmin; jdom(2)=jmin
     jdom(3)=jmax; jdom(4)=jmax
     jdom(5)=jmin; jdom(6)=jmin
     jdom(7)=jmax; jdom(8)=jmax
     kdom(1)=kmin; kdom(2)=kmin
     kdom(3)=kmin; kdom(4)=kmin
     kdom(5)=kmax; kdom(6)=kmax
     kdom(7)=kmax; kdom(8)=kmax
     
     do i=1,ndom
        if(bit_length>0)then
           call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+1.0D0)*dkey
     end do
     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (   bound_key(impi-1).le.bounding_min(i).and.&
                & bound_key(impi  ).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (   bound_key(impi-1).lt.bounding_max(i).and.&
                & bound_key(impi  ).ge.bounding_max(i))then
              cpu_max(i)=impi
           endif
        end do
     end do
     
     ncpu_read=0
     do i=1,ndom
        do j=cpu_min(i),cpu_max(i)
           if(.not. cpu_read(j))then
              ncpu_read=ncpu_read+1
              cpu_list(ncpu_read)=j
              cpu_read(j)=.true.
           endif
        enddo
     enddo
  else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  end  if

  npart=0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     close(1)
     npart=npart+npart2
  end do
  write(*,*)'Found ',npart,' particles'

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  mtot=0.0d0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     write(*,*)'Processing file '//TRIM(nomfich)
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)!localseed
     read(1)nstar
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     allocate(metal(1:npart2))
     allocate(x(1:npart2,1:ndim2))
     do i=1,ndim
        read(1)m
        x(1:npart2,i)=m
     end do
     do i=1,ndim
        read(1)m
     end do
     read(1)m
#ifdef OLD
     if(nstar.gt.0) then
        read(1)!age
        read(1)metal
     endif
     read(1)!id
     read(1)!level
#else
     read(1)!id
     read(1)!level
     if(nstar.gt.0)then
        read(1)!age
        read(1)metal
     endif
#endif     

     close(1)
     if(periodic)then
     do i=1,npart2
        ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
             &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
             &   x(i,3)>=zmin.and.x(i,3)<=zmax)
        if(star) then
           ok_part = ok_part.and.metal(i)>0.0
        endif
        if(ok_part)then
           ddx=(x(i,idim)-xxmin)/dx
           ddy=(x(i,jdim)-yymin)/dy
           ix=ddx
           iy=ddy
           ddx=ddx-ix
           ddy=ddy-iy
           dex=1.0-ddx
           dey=1.0-ddy
           if(ix<0)ix=ix+nx
           if(ix>=nx)ix=ix-nx
           if(iy<0)iy=iy+ny
           if(iy>=ny)iy=iy-ny
           ixp1=ix+1
           iyp1=iy+1
           if(ixp1<0)ixp1=ixp1+nx
           if(ixp1>=nx)ixp1=ixp1-nx
           if(iyp1<0)iyp1=iyp1+ny
           if(iyp1>=ny)iyp1=iyp1-ny
           map(ix  ,iy  )=map(ix  ,iy  )+m(i)*dex*dey
           map(ix  ,iyp1)=map(ix  ,iyp1)+m(i)*dex*ddy
           map(ixp1,iy  )=map(ixp1,iy  )+m(i)*ddx*dey
           map(ixp1,iyp1)=map(ixp1,iyp1)+m(i)*ddx*ddy
           mtot=mtot+m(i)
        end if
     end do
     else
     do i=1,npart2
        ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
             &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
             &   x(i,3)>=zmin.and.x(i,3)<=zmax)
        if(star)then
           ok_part = ok_part.and.metal(i)>0.0
        endif
        if(ok_part)then
           ddx=(x(i,idim)-xxmin)/dx
           ddy=(x(i,jdim)-yymin)/dy
           ix=ddx
           iy=ddy
           ddx=ddx-ix
           ddy=ddy-iy
           dex=1.0-ddx
           dey=1.0-ddy
           if(ix>=0.and.ix<nx.and.iy>=0.and.iy<ny)then
              map(ix  ,iy  )=map(ix  ,iy  )+m(i)*dex*dey
              map(ix+1,iy  )=map(ix+1,iy  )+m(i)*ddx*dey
              map(ix  ,iy+1)=map(ix  ,iy+1)+m(i)*dex*ddy
              map(ix+1,iy+1)=map(ix+1,iy+1)+m(i)*ddx*ddy
              mtot=mtot+m(i)
           endif
        end if
     end do
     endif
     deallocate(x,m,metal)
  end do
  write(*,*)'Total mass=',mtot

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)
  open(unit=10,file=nomfich,form='unformatted')
  if(periodic)then
     write(10)nx,ny
     allocate(toto(nx,ny))
     toto=map(0:nx-1,0:ny-1)
     write(10)toto
  else
     write(10)nx+1,ny+1
     allocate(toto(nx+1,ny+1))
     toto=map(0:nx,0:ny)
     write(10)toto
  endif
  close(10)

contains

  subroutine read_params
      implicit none

      integer       :: i,n
      integer       :: iargc
      character(len=4)   :: opt
      character(len=128) :: arg
      
      n = iargc()
      if (n < 4) then
         print *, 'usage: part2map  -inp  input_dir'
         print *, '                 -out  output_file'
         print *, '                 [-dir axis] '
         print *, '                 [-xc  xc  ] '
         print *, '                 [-yc  yc  ]'
         print *, '                 [-zc  zc  ] '
         print *, '                 [-rc  rad ] '
         print *, '                 [-nx  nx  ] '
         print *, '                 [-ny  ny  ] '
         print *, '                 [-per flag] '
         print *, '                 [-str flag] '
         print *, 'ex: part2map -inp output_00001 -out map.dat'// &
              &   ' -dir z -xc 0.5 -yc 0.5 -zc 0.5 -rc 0.001'
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
         case ('-out')
            outfich = trim(arg)
         case ('-dir')
            proj = trim(arg) 
         case ('-xc')
            read (arg,*) xc
         case ('-yc')
            read (arg,*) yc
         case ('-zc') 
            read (arg,*) zc
         case ('-rc') 
            read (arg,*) rc
         case ('-nx')
            read (arg,*) nx
         case ('-ny')
            read (arg,*) ny
         case ('-per')
            read (arg,*) periodic
         case ('-str')
            read (arg,*) star
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

      xmin = xc-rc
      xmax = xc+rc
      ymin = yc-rc
      ymax = yc+rc
      zmin = zc-rc
      zmax = zc+rc

      return

    end subroutine read_params

  end program part2map

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
