program testin
  !--------------------------------------------------
  ! write sample output unformatted
  !--------------------------------------------------

  !declare
  integer::nx,ny
  integer(kind=8)::dx,dy

  real(kind=4),dimension(:),allocatable::x
  real(kind=8),dimension(:),allocatable::y

  !allocate
  allocate(x(1:5))
  allocate(y(1:5))

  !read
  open(unit=1,file="test.dat",status='old',form='unformatted')
  read(1)nx,ny
  read(1)dx,dy
  read(1)x
  read(1)y
  close(1)

  print *,"nx=",nx
  print *,"ny=",ny
  print *,"dx=",dx
  print *,"dy=",dy
  do i=1,5
     print *,"x(",i,")=",x(i)
     print *,"y(",i,")=",y(i)
  end do
end program testin
