!===============================================================================
! test_nan_guard.f90 -- unit test for the fail-fast non-finite guard of the
! unstructured solver (mod_uns_simple: nan_scan / nan_report / nan_guard,
! NAN_IER = 60).
!
! Run:  make uns_nan_guard_test && bin/uns_nan_guard_test [mesh.cas]
!       (default mesh: cases/porous_plug/plug.cas -- any .cas works, the test
!        only needs the mesh, its geometry and an initialised field set)
!
! Why this test exists (2026-10-13): divergence used to be SILENT -- the steady
! cavity probe printed it=500 du_max=2.79e71 and it=550 all-NaN, and the solver
! kept iterating, hiding both the first bad iteration and the cell that caused
! it.  The guard now stops the run at the first NaN/Inf (ier = 60 -> abort code
! 60).  The guard is pure bookkeeping over the field arrays (no physics), so it
! is unit-tested here: read a real mesh, init the fields, inject NaN / Inf and
! check the scan outcome (field, owner cell, kind), the printed report and the
! returned ier.  The scan is also checked to be *deterministic* (first offender
! in cell-index order, then in field order) because the MPI guard reports the
! lowest offending rank and relies on that ordering.
! The MPI counterpart (mod_uns_simple_mpi::nan_guard_mpi) adds a collective
! reduction on top of this core and is unit-tested separately by
! test_nan_guard_mpi.f90 (mpirun -np 2, expected exit code 60).
!===============================================================================
program test_nan_guard
   use mod_precision, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                           ieee_positive_inf
   use mod_uns_mesh, only: mesh_t
   use mod_uns_cas_reader, only: read_cas
   use mod_uns_connectivity, only: conn_t, build_connectivity
   use mod_uns_geometry, only: geom_t, compute_geometry
   use mod_uns_fields, only: fields_t, init_fields
   use mod_uns_simple, only: nan_scan, nan_guard, NAN_IER
   implicit none

   type(mesh_t)      :: m
   type(conn_t)      :: c
   type(geom_t)      :: g
   type(fields_t)    :: fld
   character(len=256) :: casfile
   character(len=32) :: bfld
   character(len=8)  :: bkind
   real(dp) :: bval, nanv, infv
   integer  :: bcell, ier, ifin, kc1, kc2, nfail, narg

   nfail = 0
   casfile = 'cases/porous_plug/plug.cas'
   narg = command_argument_count()
   if ( narg >= 1 ) call get_command_argument( 1, casfile )

   write(*,'(a,a)') 'mesh: ', trim(casfile)
   call read_cas( trim(casfile), m, ier )
   if ( ier /= 0 ) then
      write(*,'(a,i0)') 'FATAL: read_cas failed, ier=', ier
      stop 1
   end if
   call build_connectivity( m, c, ier )
   call compute_geometry( m, g, ier )
   call init_fields( m, g, fld )
   write(*,'(a,i0,a,i0)') 'mesh loaded: ncells=', m%ncells, '  nfaces=', m%nfaces
   if ( m%ncells < 4 ) then
      write(*,'(a,i0)') 'FATAL: this test needs at least 4 cells, got ', m%ncells
      stop 1
   end if
   kc1 = m%ncells / 2 + 1
   kc2 = max( 1, m%ncells / 4 )
   if ( kc1 == kc2 ) kc1 = kc1 + 1

   nanv = ieee_value( 1.0_dp, ieee_quiet_nan )
   infv = ieee_value( 1.0_dp, ieee_positive_inf )   ! +Inf, portably

   ! ---- (1) finite field set: no false positive ------------------------------
   call nan_scan( m, fld, m%ncells, bfld, bcell, bval, bkind, ifin )
   call check( ifin == 0, 1, 'clean field set scans finite (ifin=0)', nfail )
   call nan_guard( m, g, fld, m%ncells, 'unit-clean', 0, ier )
   call check( ier == 0, 1, 'clean nan_guard returns ier=0', nfail )

   ! ---- (2) NaN in u_y at kc1 ------------------------------------------------
   fld%u(2,kc1) = nanv
   call nan_scan( m, fld, m%ncells, bfld, bcell, bval, bkind, ifin )
   call check( ifin == 1, 2, 'NaN is detected (ifin=1)', nfail )
   call check( trim(bfld) == 'u_y', 2, 'NaN scan names the field u_y', nfail )
   call check( bcell == kc1, 2, 'NaN scan reports the owner cell', nfail )
   call check( trim(bkind) == 'NaN', 2, 'NaN is labelled NaN (not Inf)', nfail )
   call nan_guard( m, g, fld, m%ncells, 'unit-nan', 77, ier )
   call check( ier == NAN_IER, 2, 'nan_guard returns NAN_IER', nfail )
   call check( NAN_IER == 60, 2, 'NAN_IER is 60 (documented abort code)', nfail )
   fld%u(2,kc1) = 0.0_dp

   ! ---- (3) Inf in p at kc2 --------------------------------------------------
   fld%p(kc2) = infv
   call nan_scan( m, fld, m%ncells, bfld, bcell, bval, bkind, ifin )
   call check( ifin == 1 .and. trim(bfld) == 'p', 3, 'Inf in p is detected', nfail )
   call check( bcell == kc2, 3, 'Inf scan reports the owner cell', nfail )
   call check( trim(bkind) == 'Inf', 3, 'Inf is labelled Inf (not NaN)', nfail )
   call nan_guard( m, g, fld, m%ncells, 'unit-inf', 5, ier )
   call check( ier == NAN_IER, 3, 'nan_guard on Inf returns NAN_IER', nfail )
   fld%p(kc2) = 0.0_dp

   ! ---- (4) deterministic choice of the first offender -----------------------
   ! NaN in T at the SMALLER cell index must win over Inf in u_x at the larger
   ! one, even though u_x comes first in the field order (cell index dominates)
   fld%T(kc2)   = nanv
   fld%u(1,kc1) = infv
   call nan_scan( m, fld, m%ncells, bfld, bcell, bval, bkind, ifin )
   call check( ifin == 1 .and. bcell == kc2 .and. trim(bfld) == 'T', 4, &
               'first offender = lowest cell index (T at kc2, not u_x at kc1)', &
               nfail )
   fld%T(kc2)   = 0.0_dp
   fld%u(1,kc1) = 0.0_dp

   ! ---- (5) LTNE solid temperature is scanned too ----------------------------
   fld%T_s(kc1) = nanv
   call nan_scan( m, fld, m%ncells, bfld, bcell, bval, bkind, ifin )
   call check( ifin == 1 .and. trim(bfld) == 'T_s', 5, &
               'T_s (LTNE solid temperature) is scanned', nfail )
   fld%T_s(kc1) = 0.0_dp

   ! ---- (6) no state left behind: the guard recovers -------------------------
   call nan_guard( m, g, fld, m%ncells, 'unit-clean-2', 999, ier )
   call check( ier == 0, 6, 'guard is stateless: ier=0 after cleanup', nfail )

   write(*,'(a,i0,a,i0,a)') '=== 6 test groups, ', nfail, ' failures ==='
   if ( nfail > 0 ) then
      write(*,'(a)') 'FAILED'
      stop 1
   end if
   write(*,'(a)') 'ALL NaN-GUARD CHECKS PASSED'

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
         write(*,'(a,i0,a,a)') 'PASS (', group, '): ', trim(what)
      else
         write(*,'(a,i0,a,a)') 'FAIL (', group, '): ', trim(what)
         nfail = nfail + 1
      end if
   end subroutine check

end program test_nan_guard
