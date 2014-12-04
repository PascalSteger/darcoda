module fmesh
  implicit none
  integer, parameter:: dp=kind(0.d0)                   ! double precision
contains

  subroutine mesh_exp(r_min, r_max, N, mesh, pp, msg, nchars, off)
    real(dp) :: r_min
    real(dp) :: r_max
    integer :: N
    real(dp) :: mesh(N)
    real(dp) :: a
    integer :: i, nchars
    logical :: pp
    character(LEN=nchars):: msg

    a = (r_max-r_min)/N
    do i = 1, N, 1
       mesh(i) = r_min + i*a + off(i)
    end do
    if ( pp ) then
       print *,' and in fmesh, we have ',msg

       do i=1, N, 1
          write(*,*)mesh(i)
       end do
    end if
  end subroutine mesh_exp

end module fmesh
