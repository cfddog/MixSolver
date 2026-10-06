!===============================================================================
! mod_mpi.f90 -- Thin MPI wrapper for UNSSolver parallel mode
!
! Provides cached rank/nprocs and convenience helpers so that the rest of the
! code does not need to scatter 'use mpi' or pass comm/rank arguments through
! every call site.  Call mpi_bootstrap() at program start and mpi_shutdown()
! at the end.
!
! Also provides f2c_comm() to convert a Fortran MPI communicator handle to a
! C MPI_Comm pointer, needed when calling ParMETIS C functions from Fortran
! (OpenMPI Fortran/C communicator ABIs differ).
!===============================================================================
module mod_uns_mpi_core
   use mpi
   use iso_c_binding, only: c_ptr, c_int
   implicit none
   private
   public :: mpi_bootstrap, mpi_shutdown, mpi_check, myrank, nprocs, &
             mpi_comm, f2c_comm

   ! cached world communicator (kept as integer for portability)
   integer, protected :: mpi_comm = MPI_COMM_WORLD
   integer, protected :: myrank   = 0
   integer, protected :: nprocs   = 1

   ! C interface for MPI_Comm_f2c: converts Fortran comm handle to C MPI_Comm
   interface
      function c_mpi_comm_f2c( f_comm ) bind(c, name="MPI_Comm_f2c")
         import :: c_ptr, c_int
         integer(c_int), value :: f_comm
         type(c_ptr)           :: c_mpi_comm_f2c
      end function c_mpi_comm_f2c
   end interface

contains

   ! Convert a Fortran MPI communicator to a type(c_ptr) holding the C MPI_Comm.
   ! The result can be passed (by reference) to C functions expecting MPI_Comm *.
   function f2c_comm( f_comm ) result( c_comm )
      integer, intent(in)  :: f_comm
      type(c_ptr)          :: c_comm
      c_comm = c_mpi_comm_f2c( int( f_comm, c_int ) )
   end function f2c_comm

   subroutine mpi_bootstrap()
      integer :: ierr
      call MPI_Init( ierr )
      call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )
      call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )
   end subroutine mpi_bootstrap

   subroutine mpi_shutdown()
      integer :: ierr
      call MPI_Finalize( ierr )
   end subroutine mpi_shutdown

   subroutine mpi_check( ierr, msg )
      integer,          intent(in) :: ierr
      character(len=*), intent(in), optional :: msg
      integer :: errstr_len, errcls, iabort_ierr
      character(len=MPI_MAX_ERROR_STRING) :: errstr
      if ( ierr /= MPI_SUCCESS ) then
         call MPI_Error_class( ierr, errcls, errstr_len )
         call MPI_Error_string( ierr, errstr, errstr_len, errcls )
         if ( present(msg) ) then
            write(*,'(a,a,a,i0,a,a)') 'MPI ERROR [', trim(msg), &
               '] code=', ierr, ' msg=', errstr(1:errstr_len)
         else
            write(*,'(a,i0,a,a)') 'MPI ERROR code=', ierr, &
               ' msg=', errstr(1:errstr_len)
         end if
         call MPI_Abort( MPI_COMM_WORLD, ierr, iabort_ierr )
      end if
   end subroutine mpi_check

end module mod_uns_mpi_core
