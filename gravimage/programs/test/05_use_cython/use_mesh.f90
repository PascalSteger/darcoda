program main
  use fmesh
  !integer, parameter:: dp=kind(0.d0)                   ! double precision

  real(dp) :: r_min = 0.0_dp
  real(dp) :: r_max = 1.0_dp
  integer :: N = 5
  real(dp), dimension(:),allocatable :: mesh
  integer:: i

  allocate (mesh(N))
  call mesh_exp(r_min, r_max, N, mesh)

  do i = 1, N, 1
     write(*,*)'entry ',i,' in mesh is ',mesh(i)
  end do
end program main
