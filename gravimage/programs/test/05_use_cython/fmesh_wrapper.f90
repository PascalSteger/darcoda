module fmesh_wrapper
  use iso_c_binding, only: c_double, c_int
  use fmesh, only: mesh_exp

  implicit none

contains
  subroutine c_mesh_exp(r_min, r_max, N, mesh) bind(c)
    real(c_double), intent(in) :: r_min
    real(c_double), intent(in) :: r_max
    integer(c_int), intent(in) :: N
    real(c_double), intent(out) :: mesh(N)

    ! real(c_double) :: a
    ! integer(c_int) :: i
    ! a = (r_max-r_min)/N
    ! do i = 1, N, 1
    !    mesh(i) = r_min + i*a
    ! end do
    call mesh_exp(r_min, r_max, N, mesh)
  end subroutine c_mesh_exp

end module fmesh_wrapper
