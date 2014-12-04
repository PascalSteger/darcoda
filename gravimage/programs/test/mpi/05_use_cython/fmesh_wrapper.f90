module fmesh_wrapper
  use iso_c_binding, only: c_double, c_int, c_bool, c_char, c_null_char
  use fmesh, only: mesh_exp
  implicit none

contains


  subroutine c_mesh_exp(r_min, r_max, N, mesh, pp, str) bind(c, name="c_mesh_exp")
    real(c_double), intent(in) :: r_min
    real(c_double), intent(in) :: r_max
    integer(c_int), intent(in) :: N
    real(c_double), intent(out) :: mesh(N)
    logical(c_bool), intent(in) :: pp
    character(kind=c_char,len=1), intent(in) :: str(*)
    character(len=:), allocatable :: msg
    integer i, nchars
    i = 1
    do
       if (str(i) == c_null_char) exit
       i = i + 1
    end do
    nchars = i - 1  ! Exclude null character from Fortran string
    allocate(character(len=nchars) :: msg)
    msg = transfer(str(1:nchars), msg)

    print *, ' .. in fmesh_wrapper, we have ',msg
    call mesh_exp(r_min, r_max, N, mesh, logical(pp), msg, nchars)
  end subroutine c_mesh_exp

end module fmesh_wrapper
