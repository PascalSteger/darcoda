module itimes
  implicit none
  include 'mpif.h'

  integer mpi_status(MPI_STATUS_SIZE), errcode
  integer ( kind = 4 ) error
  integer ( kind = 4 ) id
  integer ( kind = 4 ) p_here
  real ( kind = 8 ) wtime

contains
  subroutine test(nest_mine)
    implicit none
    integer nest_mine
    nest_mine = nest_mine*2

    write(*,*)'test successfully started'

    write(*,*)'MPI test:'
    call MPI_Init ( error )
    call MPI_Comm_size ( MPI_COMM_WORLD, p_here, error )
    call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
    if ( id == 0 ) then
       wtime = MPI_Wtime ( )
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'HELLO_MPI - Master process:'
       write ( *, '(a)' ) '  FORTRAN90/MPI version'
       write ( *, '(a)' ) '  An MPI test program.'
       write ( *, '(a)' ) ' '
       write ( *, '(a,i8)' ) '  The number of processes is ', p_here
       write ( *, '(a)' ) ' '
    end if
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) '  Process ', id, ' says "Hello, world!"'
    if ( id == 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'HELLO_MPI - Master process:'
       write ( *, '(a)' ) '  Normal end of execution: "Goodbye, world!".'
       wtime = MPI_Wtime ( ) - wtime
       write ( *, '(a)' ) ' '
       write ( *, '(a,g14.6,a)' ) &
            '  Elapsed wall clock time = ', wtime, ' seconds.'
    end if
    call MPI_Finalize ( error )
    stop

  end subroutine test

  subroutine mytime(th)
    use fd
    implicit none
    real(dp) th

    write(*,*) 'th - it',th
    call lprsmf(th)
  end subroutine mytime

end module itimes
