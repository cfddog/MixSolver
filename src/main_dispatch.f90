!===============================================================================
! main_dispatch.f90 -- single self-dispatching executable (bin/mixnsolver).
!
! Takes NO command-line arguments.  Running it inside a case directory probes
! that directory (cwd) for input files and automatically decides which solver
! to run:
!
!   * COUPLED : a 'mix.control' is present
!                  -> weak-coupling; reads Mesh3d.x/control.ec (struct) and
!                     unMesh.cas/unMesh.control (uns), all relative to cwd.
!   * UNS     : no mix.control, a unique <stem>.cas + same-stem <stem>.control
!               pair is found in cwd (canonical unMesh.cas/unMesh.control first)
!                  -> unstructured (incompressible) solver.
!   * STRUCT  : no mix.control, no uns pair, control.ec + Mesh3d.x present
!                  -> structured (compressible) solver.
!
! Every rank shares the cwd, so probing is deterministic on all ranks with NO
! MPI communication (and it must complete before any MPI_Comm_split, because
! struct vs uns membership depends on the mode).
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

   integer, parameter :: MODE_UNSET = 0, MODE_COUPLED = 1, MODE_UNS = 2, MODE_STRUCT = 3
   integer, parameter :: MAXCAND = 64

   integer :: ierr, rank, nproc, colour, key, mode
   integer :: comm_sub, n_struct_ranks, n_pair
   character(len=512) :: scas, sctl, ucas, uctl
   character(len=512) :: stems(MAXCAND)

   ! ---- MPI init -------------------------------------------------------------
   call MPI_Init( ierr )
   call MPI_Comm_rank( MPI_COMM_WORLD, rank, ierr )
   call MPI_Comm_size( MPI_COMM_WORLD, nproc, ierr )

   ! ---- mode detection (cwd, no MPI comm) ------------------------------------
   mode = MODE_UNSET
   call detect_in_cwd( mode, scas, sctl, ucas, uctl, n_pair, stems, rank )

   ! ---- dispatch -------------------------------------------------------------
   select case ( mode )
   case ( MODE_COUPLED )
      ! reference state: all ranks read the same mix.control
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
         write(*,'(a)') '  MixNSSolver [mixnsolver] -- auto-dispatch: COUPLED'
         write(*,'(a)') '=========================================================='
      end if
      if ( rank < n_struct_ranks ) then
         call struct_group_driver( comm_sub, scas, sctl, rank, nproc )
      else
         call uns_group_driver( comm_sub, ucas, uctl, rank, nproc )
      end if

   case ( MODE_STRUCT )
      if ( rank == 0 ) then
         write(*,'(a)') '=========================================================='
         write(*,'(a)') '  MixNSSolver [mixnsolver] -- auto-dispatch: STRUCT (compressible)'
         write(*,'(a)') '=========================================================='
      end if
      call struct_solver_run( MPI_COMM_WORLD )

   case ( MODE_UNS )
      if ( rank == 0 ) then
         write(*,'(a)') '=========================================================='
         write(*,'(a)') '  MixNSSolver [mixnsolver] -- auto-dispatch: UNS (incompressible)'
         write(*,'(a,a)') '  Mesh file    : ', trim(ucas)
         write(*,'(a,a)') '  Control file : ', trim(uctl)
         write(*,'(a)') '=========================================================='
      end if
      call run_uns( MPI_COMM_WORLD, ucas, uctl )

   case default
      ! detection already reported the precise failure on rank 0
      call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
   end select

   call MPI_Finalize( ierr )

contains

   !---------------------------------------------------------------------------
   ! detect_in_cwd -- decide the mode from cwd files.  Deterministic on every
   ! rank (no MPI comm).  On failure prints the candidates and leaves mode
   ! unset so the caller aborts.
   !---------------------------------------------------------------------------
   subroutine detect_in_cwd( dmode, dscas, dsctl, ducas, ductl, dn_pair, dstems, drank )
      implicit none
      integer, intent(out) :: dmode, dn_pair
      character(len=512), intent(out) :: dscas, dsctl, ducas, ductl
      character(len=512), intent(out) :: dstems(:)
      integer, intent(in) :: drank

      logical :: has_mix, has_cec, has_mesh
      logical :: struct_present

      dn_pair = 0
      dstems  = ''

      inquire( file='mix.control', exist=has_mix )
      if ( has_mix ) then
         dmode = MODE_COUPLED
         dscas = 'Mesh3d.x'          ! struct reads Mesh3d.x internally; informational
         dsctl = 'control.ec'
         ducas = 'unMesh.cas'
         ductl = 'unMesh.control'
         return
      end if

      ! no mix.control: try the unstructured pair, then struct
      call find_uns_pair( ducas, ductl, dn_pair, dstems )

      inquire( file='control.ec', exist=has_cec )
      inquire( file='Mesh3d.x', exist=has_mesh )
      struct_present = ( has_cec .and. has_mesh )

      if ( dn_pair > 0 .and. struct_present ) then
         dmode = MODE_UNSET
         if ( drank == 0 ) then
            write(*,'(a)') '[mixnsolver] ERROR: ambiguous cwd -- both an uns .cas/.control'
            write(*,'(a)') '    pair and control.ec+Mesh3d.x are present.'
            write(*,'(a)') '    Add a mix.control (=> COUPLED) or move files to disambiguate.'
            call list_uns_pairs( dn_pair, dstems )
         end if
         return
      end if

      if ( dn_pair > 0 ) then
         if ( dn_pair > 1 ) then
            dmode = MODE_UNSET
            if ( drank == 0 ) then
               write(*,'(a,i0,a)') '[mixnsolver] ERROR: ', dn_pair, &
                     ' unstructured .cas/.control pairs found -- ambiguous.'
               call list_uns_pairs( dn_pair, dstems )
               write(*,'(a)') '    Use the canonical names unMesh.cas + unMesh.control,'
               write(*,'(a)') '    or reduce the directory to a single pair.'
            end if
            return
         end if
         dmode = MODE_UNS
         ! ducas/ductl already set by find_uns_pair for dn_pair==1
         return
      end if

      if ( struct_present ) then
         dmode = MODE_STRUCT
         dscas = 'Mesh3d.x'
         dsctl = 'control.ec'
         ducas = ''
         ductl = ''
         return
      end if

      dmode = MODE_UNSET
      if ( drank == 0 ) then
         write(*,'(a)') '[mixnsolver] ERROR: cannot auto-detect a case in this directory.'
         write(*,'(a,l1)') '    mix.control present       : ', has_mix
         write(*,'(a,l1)') '    control.ec present        : ', has_cec
         write(*,'(a,l1)') '    Mesh3d.x present          : ', has_mesh
         write(*,'(a,l1)') '    uns .cas/.control pair    : ', .false.
         write(*,'(a)') '    Need one of:'
         write(*,'(a)') '      COUPLED: mix.control [+ Mesh3d.x control.ec unMesh.cas unMesh.control]'
         write(*,'(a)') '      UNS    : a unique <stem>.cas + <stem>.control'
         write(*,'(a)') '      STRUCT : control.ec + Mesh3d.x'
      end if
   end subroutine detect_in_cwd

   !---------------------------------------------------------------------------
   ! find_uns_pair -- look for a unique <stem>.cas + <stem>.control pair in cwd.
   ! Prefers the canonical names unMesh.cas / unMesh.control.  Falls back to a
   ! directory scan ('ls *.cas') to collect every candidate; an ambiguous
   ! directory yields n_pair > 1 and the caller errors out.
   ! Sets ucas/uctl to the first pair when n_pair == 1.
   !---------------------------------------------------------------------------
   subroutine find_uns_pair( oucas, ouctl, on_pair, ostems )
      implicit none
      character(len=512), intent(out) :: oucas, ouctl
      integer, intent(out) :: on_pair
      character(len=512), intent(out) :: ostems(:)

      logical :: c1, c2
      integer :: u, ios, nl, j, dotpos
      character(len=512) :: line, stem

      on_pair = 0
      ostems  = ''

      ! canonical pair first
      inquire( file='unMesh.cas', exist=c1 )
      inquire( file='unMesh.control', exist=c2 )
      if ( c1 .and. c2 ) then
         oucas = 'unMesh.cas'
         ouctl = 'unMesh.control'
         on_pair = 1
         ostems(1) = 'unMesh'
         return
      end if

      ! generic scan; a ".cas" whose same-stem ".control" also exists is a pair
      call execute_command_line( 'ls *.cas > .mixnsolver_caslist 2>/dev/null', wait=.true. )
      open( newunit=u, file='.mixnsolver_caslist', status='old', action='read', iostat=ios )
      if ( ios == 0 ) then
         nl = 0
         do
            read( u, '(a)', iostat=ios ) line
            if ( ios /= 0 ) exit
            nl = nl + 1
            if ( nl > MAXCAND ) exit
            ! strip trailing carriage-return/blank, drop the ".cas" extension
            line = trim(line)
            dotpos = index( line, '.cas' )
            if ( dotpos > 1 ) then
               stem = line(1:dotpos-1)
               inquire( file=trim(stem)//'.control', exist=c2 )
               if ( c2 ) then
                  ! register if not already listed
                  do j = 1, on_pair
                     if ( trim(ostems(j)) == trim(stem) ) exit
                     if ( j == on_pair ) then
                        on_pair = on_pair + 1
                        ostems(on_pair) = stem
                     end if
                  end do
                  if ( on_pair == 0 ) then
                     on_pair = 1
                     ostems(1) = stem
                  end if
               end if
            end if
         end do
         close( u )
      end if
      call execute_command_line( 'rm -f .mixnsolver_caslist', wait=.true. )

      ! de-duplicate is inherent above (loop guard); enforce the 1:1 result
      if ( on_pair == 1 ) then
         oucas = trim(ostems(1))//'.cas'
         ouctl = trim(ostems(1))//'.control'
      end if
   end subroutine find_uns_pair

   subroutine list_uns_pairs( on_pair, ostems )
      implicit none
      integer, intent(in) :: on_pair
      character(len=512), intent(in) :: ostems(:)
      integer :: j
      do j = 1, on_pair
         write(*,'(a,a)') '      - ', trim(ostems(j))//'.cas + '//trim(ostems(j))//'.control'
      end do
   end subroutine list_uns_pairs

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