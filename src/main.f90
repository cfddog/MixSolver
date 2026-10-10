!===============================================================================
! main.f90 -- mixed structured/unstructured weak-coupling driver (phase 6).
!
! Architecture:
!   - All ranks share MPI_COMM_WORLD.
!   - Ranks are split into two sub-communicators by colour:
!       STRUCT_GROUP (colour 0) -> runs the structured solver
!       UNS_GROUP    (colour 1) -> runs the unstructured solver
!   - The split is controlled by the first n_struct ranks (from mix.control or
!     command line).  Default: rank 0 = struct, everyone else = uns.
!   - Each group runs its own solver init + coupling iteration loop.  Between
!     iterations the interface state is exchanged across the groups (phase 5
!     units + exchange modules).
!
! The two group drivers were moved verbatim into mod_mix_driver so they are
! also reusable by the zero-argument self-dispatching bin/mixnsolver.  The
! legacy behaviour of this program is unchanged.
!
! Usage: mpirun -np N bin/mixsolver_mpi  [mix.control]  [struct.cas/.control]
!                                                     [uns.cas uns.control]
!===============================================================================
program mixsolver
   use mpi
   use mod_mix_driver, only: struct_group_driver, uns_group_driver
   use mod_reference_state, only: read_mix_control
   implicit none

   integer, parameter :: STRUCT_GROUP = 0, UNS_GROUP = 1

   integer :: ierr, rank, nproc, colour, key
   integer :: comm_struct, comm_uns
   integer :: n_struct_ranks
   character(len=512) :: mixfile, scas, sctl, ucas, uctl

   ! ---- MPI init -------------------------------------------------------------
   call MPI_Init( ierr )
   call MPI_Comm_rank( MPI_COMM_WORLD, rank, ierr )
   call MPI_Comm_size( MPI_COMM_WORLD, nproc, ierr )

   ! ---- arguments ------------------------------------------------------------
   mixfile = 'mix.control'
   scas    = 'grid_BC/Mesh3d.x'
   sctl    = 'grid_BC/control.ec'
   ucas    = 'grid_BC/unMesh.cas'
   uctl    = 'grid_BC/unMesh.control'
   if ( command_argument_count() >= 1 ) call get_command_argument(1, mixfile)
   if ( command_argument_count() >= 2 ) call get_command_argument(2, scas)
   if ( command_argument_count() >= 3 ) call get_command_argument(3, sctl)
   if ( command_argument_count() >= 4 ) call get_command_argument(4, ucas)
   if ( command_argument_count() >= 5 ) call get_command_argument(5, uctl)

   ! ---- reference state (all ranks read the same mix.control) ----------------
   call read_mix_control( trim(mixfile), ierr )

   ! ---- split communicator ---------------------------------------------------
   ! Default: rank 0 -> structured, rest -> unstructured.
   ! Special case: nproc==1 -> single rank runs both (struct first, then uns).
   n_struct_ranks = 1
   if ( nproc == 1 ) then
      n_struct_ranks = 0   ! all ranks -> uns (no struct side)
   end if
   if ( rank < n_struct_ranks ) then
      colour = STRUCT_GROUP
   else
      colour = UNS_GROUP
   end if
   key = rank

   call MPI_Comm_split( MPI_COMM_WORLD, colour, key, comm_struct, ierr )
   comm_uns = comm_struct

   if ( rank == 0 ) then
      write(*,'(a)') ''
      write(*,'(a)') '=========================================================='
      write(*,'(a)') '  MixNSSolver -- weak-coupling driver (phase 6)'
      write(*,'(a)') '=========================================================='
      write(*,'(a,i0)') '  total ranks   : ', nproc
      write(*,'(a,i0)') '  struct ranks  : ', n_struct_ranks
      write(*,'(a,i0)') '  uns ranks     : ', nproc - n_struct_ranks
      write(*,'(a,a)')    '  mix.control   : ', trim(mixfile)
      write(*,'(a)') '=========================================================='
   end if

   ! ---- dispatch by group ----------------------------------------------------
   select case ( colour )
   case ( STRUCT_GROUP )
      call struct_group_driver( comm_struct, scas, sctl, rank, nproc )
   case ( UNS_GROUP )
      call uns_group_driver( comm_uns, ucas, uctl, rank, nproc )
   end select

   call MPI_Finalize( ierr )

end program mixsolver