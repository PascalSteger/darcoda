program test_oob
  implicit none
  real b(10)
  integer i
  do i=1, 10
     write(*,*) "Before set_val !"
     call set_val(b(i), 30000);
  enddo
  write(*,*) "--> ", b(1)
end program test_oob

subroutine set_val(in, dim)
  implicit none
  integer dim
  real in, out(dim*dim)
  write(*,*) "Inside set_val !"
  in = out(1) + out(dim)
end subroutine set_val
