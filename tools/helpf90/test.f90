program testout
  !--------------------------------------------------
  ! write sample output unformatted
  !--------------------------------------------------

  !declare
  integer::nx,ny
  integer(kind=8)::dx,dy

  real(kind=4),dimension(:),allocatable::x
  real(kind=8),dimension(:),allocatable::y

  !set
  nx=2
  ny=4567
  dx=2
  dy=4567
  allocate(x(1:5))
  allocate(y(1:5))
  x(1)=0.3
  x(2)=0.4
  x(3)=0.5
  x(4)=0.6
  x(5)=0.7
  y(1)=0.3
  y(2)=0.4
  y(3)=0.5
  y(4)=0.6
  y(5)=0.7

  !write
  open(unit=1,file="test.dat",status='new',form='unformatted')
  write(1)nx,ny
  write(1)dx,dy
  write(1)x
  write(1)y
  close(1)
  
end program testout
