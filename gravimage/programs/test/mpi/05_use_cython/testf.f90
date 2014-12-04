program test
  integer, parameter:: dp=kind(0.d0)                   ! double precision

  integer ii, jj, kk
  common/ijk/ ii, jj, kk
  real(dp)  ff
  character*32 cc

  integer :: n_dim = 8
  integer :: nPar = 4
  real(dp) :: Cube(8)
  real(dp) :: lnew = -1.0
  integer :: context_pass = 0

  real(dp) :: out = -1.0
  external :: printsample
  external :: doubleijk

  ii = 2
  jj = 3
  kk = 4
  ff = 9.0567
  cc = 'Example of a character string'

  print *, 'ii=', ii, ' ff= ',ff

  !call abc(ii, doubleijk)
  call abc(ii, printsample)

  print *, 'ii=', ii
  print *, 'lnew = ',lnew

  print *, 'after external call: ii, jj, kk = ',ii, jj, kk
  print *, 'cc=',cc
  stop
end program test

subroutine abc(jj, func)
  !call loglike(Cube, n_dim, nPar, lnew, context_pass)
  n_dim = 14
  call func(n_dim)
  return
end subroutine abc
