module fmesh
  implicit none
  integer, parameter:: dp=kind(0.d0)                   ! double precision
contains

  subroutine mesh_exp(r_min, r_max, N, mesh)
    real(dp) :: r_min
    real(dp) :: r_max
    integer :: N
    real(dp) :: mesh(N)
    real(dp) :: a
    integer :: i

    a = (r_max-r_min)/N
    do i = 1, N, 1
       mesh(i) = r_min + i*a
    end do
  end subroutine mesh_exp

end module fmesh
