subroutine output_rectest
  use amr_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH, neq_spec, weight_spec
  implicit none
  character(LEN=80)::filename,filedir,fileloc,stri,tmp

  integer::i,ivar,ivarn,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,NN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_nH2,mtoth=0.0
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  real(dp),dimension(nvar)::mtot

  if(verbose)write(*,*)'Entering output_rectest'

888 format('NEC: ',10e20.10)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Chemistry constants
  scale_nH2=scale_d/mH
  NN=neq_spec-1

  do ivar=1,NN
     mtot(ivar)=0.0
  end do
  ilun=ncpu+myid+10
     
  call title(myid,nchar)
  filedir='output_'//TRIM(nchar)//'/'
  filename=TRIM(filedir)//'nec_'//TRIM(nchar)//'.out'
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=10,file=fileloc,form='unformatted')
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do ivar=1,NN
                 do i=1,ncache
                    ivarn=3+ndim+ivar ! get index of NEC species,  7,8,9,10,11,12,13,14,15
                    xdp(i)=uold(ind_grid(i)+iskip,ivarn)/uold(ind_grid(i)+iskip,1)*scale_nH2/weight_spec(ivar) ! this scales to number density
                    mtot(ivar) = mtot(ivar)+xdp(i)/(2**(3*(ilevel-1)))
                 end do
                 write(10)xdp
              end do
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(10)

    ! do i=1,NN
    !     nTot=nTot+ySingle(i)*weight_spec(i)
    ! enddo
    ! do i=1,NN
    !     y(i)=ySingle(i)*weight_spec(i)/nTot
    ! enddo

  !mtoth=mtot(1)+mtot(2)
  ! ^-- useful if wished to have abundances relative to total H budget
  write(*,888)aexp,mtot(1),mtot(2),mtot(3),mtot(4),mtot(5),mtot(6),mtot(7),mtot(8),mtot(9)
end subroutine output_rectest
