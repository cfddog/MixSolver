!===============================================================================
! test_nan_guard_mpi.f90 -- unit test for the partitioned fail-fast non-finite
! guard (mod_uns_simple_mpi::nan_guard_mpi) of the unstructured solver.
!
! Run:  make uns_nan_guard_mpi_test && mpirun -np 2 bin/uns_nan_guard_mpi_test
!                                      [mesh.cas]
!       (default mesh: cases/porous_plug/plug.cas; needs >= 2 ranks)
!
! EXPECTED exit code: 60.  The last check injects a NaN on the LAST rank only
! and then aborts with that rank's ier (NAN_IER = 60) -- i.e. the run stops here
! exactly like a diverged partitioned solver does.  So:
!     exit status 60 + "ALL MPI NAN-GUARD CHECKS PASSED" in the log  = pass
!     exit status  1                                                 = FAILED
! A NaN injected on the last rank is the interesting case: rank 0 is finite, so
! only the MPI_MIN reduction can select the reporter, and the printed global
! cell id can only come from the rank-local -> global map (lm%cell_l2g).
!
! Groups:
!   (1) finite field set            -> ier = 0 on EVERY rank (no false abort)
!   (2) NaN confined to a halo cell -> ier = 0 (owned cells only: halos are
!                                      refreshed by the exchange)
!   (3) NaN on the last rank        -> ier = NAN_IER on EVERY rank, reported by
!                                      the last rank, naming that rank's global
!                                      cell, then abort(60)
!===============================================================================
program test_nan_guard_mpi
   use mpi
   use mod_precision, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mod_uns_mpi_core, only: mpi_bootstrap, mpi_check, myrank, nprocs, mpi_comm
   use mod_uns_mesh, only: mesh_t
   use mod_uns_cas_reader, only: read_cas
   use mod_uns_connectivity, only: conn_t, build_connectivity
   use mod_uns_geometry, only: geom_t, compute_geometry
   use mod_uns_local_mesh, only: local_mesh_t, extract_local_mesh
   use mod_uns_fields, only: fields_t, init_fields
   use mod_uns_simple_mpi, only: nan_guard_mpi
   use mod_uns_simple, only: NAN_IER
   implicit none

   type(mesh_t)       :: m_g
   type(conn_t)       :: c_g
   type(geom_t)       :: g_g
   type(local_mesh_t) :: lm
   type(fields_t)     :: fld
   character(len=256) :: casfile
   integer, allocatable :: part(:)
   real(dp) :: nanv
   integer  :: i, ierr, ier, nfail, nfail_all, narg, last_rank

   ! mpi_bootstrap() calls MPI_Init itself (only one MPI_Init is legal)
   call mpi_bootstrap()

   if ( nprocs < 2 ) then
      if ( myrank == 0 ) write(*,'(a)') &
         'FATAL: this test needs >= 2 ranks (mpirun -np 2 ...)'
      call MPI_Abort( MPI_COMM_WORLD, 2, ierr )
   end if

   nfail     = 0
   last_rank = nprocs - 1
   casfile   = 'cases/porous_plug/plug.cas'
   narg      = command_argument_count()
   if ( narg >= 1 ) call get_command_argument( 1, casfile )

   if ( myrank == 0 ) write(*,'(a,a,a,i0)') 'mesh: ', trim(casfile), &
                                            '   ranks: ', nprocs

   ! ---- global mesh, exactly as the partitioned driver builds it -------------
   ! (resolve_cell_zones is skipped on purpose: the guard must survive a mesh
   !  whose m%cztype is not allocated yet -- is_solid_cell is guarded for that)
   call read_cas( trim(casfile), m_g, ierr )
   if ( ierr /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: read_cas failed, ier=', ierr
      call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
   end if
   call build_connectivity( m_g, c_g, ierr )
   call compute_geometry( m_g, g_g, ierr )
   if ( myrank == 0 ) write(*,'(a,i0,a,i0)') 'global mesh: ncells=', &
      m_g%ncells, '  nfaces=', m_g%nfaces
   if ( m_g%ncells < nprocs ) then
      if ( myrank == 0 ) write(*,'(a)') 'FATAL: need at least one cell per rank'
      call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
   end if

   ! deterministic round-robin partition (rank ids in part are 0-based, as
   ! extract_local_mesh expects): cell i -> rank mod(i-1,nprocs)
   allocate( part(m_g%ncells) )
   do i = 1, m_g%ncells
      part(i) = mod( i - 1, nprocs )
   end do
   call extract_local_mesh( m_g, c_g, g_g, part, lm, ierr )
   if ( ierr /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: extract_local_mesh ier=', ierr
      call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
   end if
   call init_fields( lm%m, lm%g, fld )

   ! round-robin layout => the first owned cell of rank r is global cell r+1
   call check( lm%nowned >= 1, 1, 'rank owns at least one cell', nfail )
   call check( lm%cell_l2g(1) == myrank + 1, 1, &
               'rank-local cell 1 maps to global cell myrank+1', nfail )

   nanv = ieee_value( 1.0_dp, ieee_quiet_nan )

   ! ---- (1) clean field set: every rank must agree there is nothing wrong ----
   call nan_guard_mpi( lm, fld, 'mpi-clean', 0, ier )
   call check( ier == 0, 1, 'clean field set: ier=0 on this rank', nfail )

   ! ---- (2) NaN confined to a HALO cell: not owned, so not an abort ----------
   if ( lm%m%ncells > lm%nowned ) then
      fld%p(lm%nowned + 1) = nanv
      call nan_guard_mpi( lm, fld, 'mpi-halo', 0, ier )
      call check( ier == 0, 2, &
                  'NaN in a halo cell does not abort (owned cells only)', nfail )
      fld%p(lm%nowned + 1) = 0.0_dp
   else
      fld%p(1) = nanv                        ! no halo: only rank 0 has one
      call nan_guard_mpi( lm, fld, 'mpi-halo', 0, ier )
      if ( myrank == 0 ) then
         call check( ier == NAN_IER, 2, 'no halo here: owned NaN still caught', &
                     nfail )
      end if
      fld%p(1) = 0.0_dp
   end if
   call nan_guard_mpi( lm, fld, 'mpi-clean-2', 0, ier )
   call check( ier == 0, 2, 'guard is stateless: ier=0 after cleanup', nfail )

   ! ---- (3) NaN on the LAST rank only: the reduction must find it ------------
   if ( myrank == last_rank ) fld%p(1) = nanv      ! global cell nprocs
   call nan_guard_mpi( lm, fld, 'mpi-nan-last-rank', 42, ier )
   call check( ier == NAN_IER, 3, &
               'NaN on the last rank: ier=NAN_IER on every rank', nfail )
   call check( NAN_IER == 60, 3, 'NAN_IER is 60 (documented abort code)', nfail )

   ! the diagnosis for group (3) is printed by rank nprocs-1 itself (nan_report)
   flush(6)

   call MPI_Allreduce( nfail, nfail_all, 1, MPI_INTEGER, MPI_SUM, mpi_comm, ierr )
   call mpi_check( ierr, 'test summary' )
   if ( myrank == 0 .and. nfail_all == 0 ) then
      write(*,'(a,i0,a)') '=== 3 groups on ', nprocs, ' ranks, 0 failures ==='
      write(*,'(a)') 'ALL MPI NAN-GUARD CHECKS PASSED'
   else if ( myrank == 0 ) then
      write(*,'(a,i0,a)') '=== FAILED: ', nfail_all, ' assertion(s) ==='
   end if
   flush(6)

   ! always abort: 60 when the guard behaved (the runner expects 60), else 1
   if ( nfail_all == 0 ) then
      call MPI_Abort( mpi_comm, NAN_IER, ierr )
   else
      call MPI_Abort( mpi_comm, 1, ierr )
   end if

contains

   !----------------------------------------------------------------------------
   ! check -- count a failed assertion (the group number only labels the text)
   !----------------------------------------------------------------------------
   subroutine check( ok, group, what, nfail )
      logical,          intent(in)    :: ok
      integer,          intent(in)    :: group
      character(len=*), intent(in)    :: what
      integer,          intent(inout) :: nfail

      if ( ok ) then
         write(*,'(a,i0,a,i0,a,a)') 'PASS (', group, ', rank ', myrank, '): ', &
                                    trim(what)
      else
         write(*,'(a,i0,a,i0,a,a)') 'FAIL (', group, ', rank ', myrank, '): ', &
                                    trim(what)
         nfail = nfail + 1
      end if
   end subroutine check

end program test_nan_guard_mpi
