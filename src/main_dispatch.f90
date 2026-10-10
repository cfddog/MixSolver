!===============================================================================
! main_dispatch.f90 -- single self-dispatching executable (bin/mixnsolver).
!
! Takes NO command-line arguments.  The coupling / solver mode is declared in a
! REQUIRED co-located control file `mix.control` via the key `mode`:
!
!     mode = coupled | struct | uns
!
!   * COUPLED : weak-coupling; reads Mesh3d.x + control.ec (struct) and
!              unMesh.cas + unMesh.control (uns), all relative to cwd, and runs
!              the same split + group drivers as bin/mixsolver.
!   * STRUCT  : structured (compressible) solver; reads control.ec + Mesh3d.x.
!   * UNS     : unstructured (incompressible) solver; reads unMesh.cas +
!              unMesh.control.
!
! The mode is authoritative -- the dispatcher does NOT guess from which files
! happen to be present.  `mix.control` must exist and contain a valid `mode`
! key, otherwise the run aborts with a clear error.
!
! Every rank shares the cwd and parses the same mix.control, so the decision is
! deterministic on all ranks with NO MPI communication (and it completes before
! any MPI_Comm_split, because struct vs uns membership depends on the mode).
!
! Run modes:
!   - COUPLED : same split + group drivers as bin/mixsolver (mod_mix_driver).
!   - STRUCT  : mod_struct_driver::struct_solver_run on the full MPI_COMM_WORLD.
!   - UNS     : mod_uns_driver::uns_solver_init/step on the full communicator.
!               The unstructured side holds the full mesh on every rank (no
!               METIS partition), so '-np N' is a redundant-but-numerically-
!               identical run; the result equals the serial solve.
!===============================================================================
program mixnsolver
   use mpi
   use mod_precision, only: dp
   use mod_mix_driver, only: struct_group_driver, uns_group_driver
   use mod_struct_driver, only: struct_solver_run
   use mod_uns_driver, only: uns_solver_init, uns_solver_step
   use mod_reference_state, only: read_mix_control
   implicit none

   integer, parameter :: MODE_COUPLED = 1, MODE_UNS = 2, MODE_STRUCT = 3

   integer :: ierr, rank, nproc, colour, key, mode
   integer :: comm_sub, n_struct_ranks

   ! ---- MPI init -------------------------------------------------------------
   call MPI_Init( ierr )
   call MPI_Comm_rank( MPI_COMM_WORLD, rank, ierr )
   call MPI_Comm_size( MPI_COMM_WORLD, nproc, ierr )

   ! ---- mode: REQUIRED from mix.control (no file-presence guessing) ----------
   call parse_mix_mode( mode, rank )

   ! ---- dispatch -------------------------------------------------------------
   select case ( mode )
   case ( MODE_COUPLED )
      ! reference state + coupling parameters from the same mix.control
      call read_mix_control( 'mix.control', ierr )

      n_struct_ranks = 1
      if ( nproc == 1 ) n_struct_ranks = 0   ! all ranks -> uns (no struct side)
      if ( rank < n_struct_ranks ) then
         colour = 0
      else
         colour = 1
      end if
      key = rank
      call MPI_Comm_split( MPI_COMM_WORLD, colour, key, comm_sub, ierr )

      if ( rank == 0 ) then
         write(*,'(a)') '=========================================================='
         write(*,'(a)') '  MixNSSolver [mixnsolver] -- mode: COUPLED'
         write(*,'(a)') '=========================================================='
      end if
      if ( rank < n_struct_ranks ) then
         call struct_group_driver( comm_sub, 'Mesh3d.x', 'control.ec', rank, nproc )
      else
         call uns_group_driver( comm_sub, 'unMesh.cas', 'unMesh.control', rank, nproc )
      end if

   case ( MODE_STRUCT )
      if ( rank == 0 ) then
         write(*,'(a)') '=========================================================='
         write(*,'(a)') '  MixNSSolver [mixnsolver] -- mode: STRUCT (compressible)'
         write(*,'(a)') '  (reads control.ec + Mesh3d.x)'
         write(*,'(a)') '=========================================================='
      end if
      call struct_solver_run( MPI_COMM_WORLD )

   case ( MODE_UNS )
      if ( rank == 0 ) then
         write(*,'(a)') '=========================================================='
         write(*,'(a)') '  MixNSSolver [mixnsolver] -- mode: UNS (incompressible)'
         write(*,'(a)') '  Mesh file    : unMesh.cas'
         write(*,'(a)') '  Control file : unMesh.control'
         write(*,'(a)') '=========================================================='
      end if
      call run_uns( MPI_COMM_WORLD, 'unMesh.cas', 'unMesh.control' )

   case default
      ! parse_mix_mode already printed the precise failure on rank 0
      call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
   end select

   call MPI_Finalize( ierr )

contains

   !---------------------------------------------------------------------------
   ! parse_mix_mode -- read `mode = <val>` from the REQUIRED co-located
   ! mix.control and set the dispatch mode.  Deterministic on every rank (no MPI
   ! comm).  Aborts with a clear message if mix.control is missing or the mode
   ! key is absent / invalid.
   !---------------------------------------------------------------------------
   subroutine parse_mix_mode( o_mode, drank )
      implicit none
      integer,          intent(out) :: o_mode
      integer,          intent(in)  :: drank

      logical :: ex
      integer :: u, ios, ie
      integer :: i
      character(len=512) :: line, ckey, val
      character(len=16)  :: lv

      o_mode = 0
      inquire( file='mix.control', exist=ex )
      if ( .not. ex ) then
         if ( drank == 0 ) then
            write(*,'(a)') '[mixnsolver] ERROR: mix.control is REQUIRED in the'
            write(*,'(a)') '    current directory to declare the solve mode.'
            write(*,'(a)') '    Add:  mode = coupled | struct | uns'
         end if
         call MPI_Abort( MPI_COMM_WORLD, 2, ios )
         return
      end if

      open( newunit=u, file='mix.control', status='old', action='read', iostat=ios )
      if ( ios /= 0 ) then
         if ( drank == 0 ) write(*,'(a)') '[mixnsolver] ERROR: cannot open mix.control'
         call MPI_Abort( MPI_COMM_WORLD, 3, ios )
         return
      end if

      do
         read( u, '(a)', iostat=ios ) line
         if ( ios /= 0 ) exit
         ie = scan( line, '#' )
         if ( ie > 0 ) line = line(:ie-1)
         line = adjustl( line )
         if ( len_trim(line) == 0 ) cycle
         ie = scan( line, '=' )
         if ( ie == 0 ) cycle
         ckey = adjustl( line(:ie-1) )
         if ( trim(ckey) /= 'mode' ) cycle
         val = adjustl( line(ie+1:) )
         ! lowercase value
         lv = ''
         do i = 1, len_trim(val)
            if ( val(i:i) >= 'A' .and. val(i:i) <= 'Z' ) then
               lv(i:i) = achar( iachar(val(i:i)) + 32 )
            else
               lv(i:i) = val(i:i)
            end if
         end do
         select case ( trim(lv) )
         case ( 'coupled', 'coupling', 'mix', '1' )
            o_mode = MODE_COUPLED
         case ( 'struct', 'structured', 'compressible', '2' )
            o_mode = MODE_STRUCT
         case ( 'uns', 'unstruct', 'unstructured', 'incompressible', '3' )
            o_mode = MODE_UNS
         case default
            if ( drank == 0 ) write(*,'(a,a,a)') '[mixnsolver] ERROR: invalid mode "', &
                  trim(val), '" in mix.control (use coupled | struct | uns)'
            call MPI_Abort( MPI_COMM_WORLD, 4, ios )
            return
         end select
      end do
      close( u )

      if ( o_mode == 0 ) then
         if ( drank == 0 ) then
            write(*,'(a)') '[mixnsolver] ERROR: mix.control has no "mode" key.'
            write(*,'(a)') '    Add:  mode = coupled | struct | uns'
         end if
         call MPI_Abort( MPI_COMM_WORLD, 5, ios )
      end if
   end subroutine parse_mix_mode

   !---------------------------------------------------------------------------
   ! run_uns -- unstructured-only solve on the full communicator.
   ! uns_solver_init runs the FULL mesh on every rank (same as the coupling
   ! uns side, no METIS partition), so all ranks step in lockstep and the field
   ! equals the serial solve.  Only the global rank 0 writes the result file.
   !---------------------------------------------------------------------------
   subroutine run_uns( comm, casfile, ctlfile )
      use mod_uns_driver, only: uns_solver_init, uns_solver_step
      use mod_uns_mesh, only: mesh_t
      use mod_uns_connectivity, only: conn_t
      use mod_uns_geometry, only: geom_t
      use mod_uns_control, only: ctrl_t
      use mod_uns_bc, only: bc_t
      use mod_uns_fields, only: fields_t
      use mod_uns_output, only: vtk_write, tecplot_write
      implicit none
      integer, intent(in) :: comm
      character(len=*), intent(in) :: casfile, ctlfile

      type(mesh_t)   :: m
      type(conn_t)   :: c
      type(geom_t)   :: g
      type(ctrl_t)   :: ctrl
      type(bc_t)     :: bcs
      type(fields_t) :: fld
      integer :: ier, ier2, j, k
      character(len=512) :: outstem, vtufile
      character(len=8)   :: outext

      call uns_solver_init( casfile, ctlfile, m, c, g, ctrl, bcs, fld, ier )
      if ( ier /= 0 ) then
         write(*,'(a,i0)') '[mixnsolver/uns] FATAL: init failed ier=', ier
         call MPI_Abort( comm, ier, ier2 )
         return
      end if

      call uns_solver_step( m, c, g, ctrl, bcs, fld, ctrl%outer_max, ier )
      if ( ier /= 0 ) then
         write(*,'(a,i0)') '[mixnsolver/uns] FATAL: solver failed ier=', ier
         call MPI_Abort( comm, ier, ier2 )
         return
      end if

      call MPI_Barrier( comm, ier2 )

      if ( rank == 0 ) then
         ! result file name: <cas-basename>.<ext>, ext from ctrl%out_format
         outstem = 'result'
         j = scan( casfile, '/', back=.true. )
         k = scan( casfile(j+1:), '.', back=.true. )
         if ( k > 1 ) outstem = casfile(j+1:j+k-1)
         if ( ctrl%out_format == 'tecplot' ) then
            outext = '.plt'
         else
            outext = '.vtu'
         end if
         vtufile = trim(outstem)//trim(outext)
         if ( ctrl%out_format == 'tecplot' ) then
            call tecplot_write( m, c, g, ctrl, fld, trim(vtufile), ier )
         else
            call vtk_write( m, c, g, fld, trim(vtufile), ier )
         end if
         if ( ier /= 0 ) then
            write(*,'(a,i0)') '[mixnsolver/uns] FATAL: result output failed ier=', ier
            call MPI_Abort( comm, ier, ier2 )
            return
         end if
         write(*,'(a,a)') '[mixnsolver/uns] wrote ', trim(vtufile)
      end if
   end subroutine run_uns

end program mixnsolver