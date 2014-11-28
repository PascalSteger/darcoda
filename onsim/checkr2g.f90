  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ramses2gadget (r2g) check tool
  ! ==============================                                         
  !                                                                   
  ! This tool reads data from r2g outputs (i.e. Gadget files) and
  ! prints debugging information
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program checkr2g

implicit none
    
  ! counters

  integer        :: i, j, ilevel, idim, idm, istar, ifile, icpu, ivar, ind,col
  integer        :: nfiles, ncpus, ndim, twotondim, firstfile, lastfile
  integer*4:: npart=11132

  ! dummy variables

  integer*4,dimension(6)         :: dummy_int
  real*4,dimension(6)            :: dummy_float
  character(len=4)  :: dummy_char

  ! files
  character(len=128)  :: dir_name, filename

  i=0
  call read_args(dir_name)
  filename = trim(dir_name) // '/r2g.0'
  !filename="try"
  print *,"File read is ",filename
  dummy_char="HEAD"
!  open(unit=11,file=filename,status='old',form='unformatted', &
!           & action='readwrite', access='direct', recl=4)
  open(unit=11,file=filename,status='replace',form='unformatted',action='write',recl=4)

  print *, "sizeof(integer)", sizeof(i)
  write(11) 'POS ', 3*8*npart+8
!  write(11) gaspart_pos, (dmpart_pos(:,col),col=1,ndmpart), (starpart_pos(:,col),col=1,nstarpart)
  
  ! write particle velocities
  write(11) 'VEL '!, 3*8*npart+8
!  write(11) gaspart_vel, (dmpart_vel(:,col),col=1,ndmpart), (starpart_vel(:,col),col=1,nstarpart)

  close(11)

end program

subroutine read_args(dir_name)
  ! reads arguments from command line: filename
  implicit none
  
  integer                        :: n, iargc
  character(len=128)             :: arg1
  character(len=128),intent(out) :: dir_name
    
  ! retrieve arguments    
  n = iargc()
  if (n /= 1) then
     write (*,*) 'usage: ./ramses2gadget input_dir'
     call exit(5)
  else
     ! first argument tells us which directory to search for r2g files
     call getarg(1,arg1)
     dir_name = trim(arg1)
  end if
  return  
end subroutine read_args

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



