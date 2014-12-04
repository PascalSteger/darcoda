module fd

  ! Double precision real kind
  integer, parameter :: dp = selected_real_kind(15)

contains

subroutine lprsmf(th)
  implicit none
  real(dp) th
  write(*,*) 'th - fd',th
end subroutine lprsmf

end module fd
