module nested_wrapper
  use iso_c_binding, only: c_double, c_int, c_char, c_null_char, c_bool
  use nested, only: nestRun
  implicit none

contains
    ! INTERFACE
    !    !the likelihood function
    !    subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
    !      integer n_dim,nPar,context_pass
    !      double precision lnew,Cube(nPar)
    !    end subroutine loglike
    ! end INTERFACE
    ! INTERFACE
    !    !the user dumper function
    !    subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, INSlogZ, logZerr, context_pass)
    !      integer nSamples, nlive, nPar, context_pass
    !      double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
    !      double precision maxLogLike, logZ, INSlogZ, logZerr
    !    end subroutine dumper
    ! end INTERFACE

  subroutine wrap_nestrun(nest_IS, nest_mmodal, nest_ceff, nest_nlive, nest_tol, nest_ef, nest_ndims, nest_totPar, nest_nCdims, maxClst,  nest_updInt, nest_Ztol, str, seed, nest_pWrap, nest_fb, nest_resume, nest_outfile, initMPI, nest_logZero, nest_maxIter,  context) bind(c)
    logical(c_bool), intent(in)::nest_IS, nest_mmodal, nest_ceff
    integer(c_int), intent(in):: nest_nlive, nest_ndims, nest_totPar, nest_nCdims, maxClst
    integer(c_int), intent(in):: nest_updInt, seed
    real(c_double), intent(in):: nest_tol, nest_ef, nest_Ztol
    integer(c_int), intent(in):: nest_pWrap(nest_ndims)
    logical(c_bool), intent(in)::nest_fb, nest_resume, nest_outfile, initMPI
    real(c_double), intent(in):: nest_logZero
    integer(c_int), intent(in):: nest_maxIter, context
    !CFUNCTYPE(c_void_p, c_int, c_int, c_int, POINTER(c_double),POINTER(c_double),POINTER(c_double), c_double,c_double,c_double,c_void_p) :: loglike
    !c_f_pointer, intent(in) :: loglike
    external :: loglike
    external :: dumper

    ! convert char array to string, pass number of chars
    character(kind=c_char,len=1), intent(in) :: str(*)
    character(len=:), allocatable :: nest_root
    integer i, nchars
    i = 1
    do
       if (str(i) == c_null_char) exit
       i = i + 1
    end do
    nchars = i - 1  ! Exclude null character from Fortran string
    allocate(character(len=nchars) :: nest_root)
    nest_root = transfer(str(1:nchars), nest_root)


    call nestRun(logical(nest_IS), logical(nest_mmodal), logical(nest_ceff), nest_nlive, nest_tol, nest_ef, nest_ndims, nest_totPar, nest_nCdims, maxClst, nest_updInt, nest_Ztol, nest_root, nchars, seed, nest_pWrap, logical(nest_fb), logical(nest_resume), logical(nest_outfile), logical(initMPI), nest_logZero, nest_maxIter, loglike, dumper, context)
    end subroutine wrap_nestRun


end module nested_wrapper
