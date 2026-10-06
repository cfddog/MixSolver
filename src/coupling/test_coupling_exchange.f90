!===============================================================================
! test_coupling_exchange.f90 -- phase-6 MPI end-to-end tests for
! mod_coupling_exchange (the cross-group count+payload protocol).
!
! Run with: mpirun -np 2 bin/coupling_test
!   rank 0 plays the struct coupling root, rank 1 the uns coupling root.
!
! Tests:
!   1. constant field struct->uns: every received value equals the expected
!      SI conversion of the non-dim constant (protocol + units in the chain)
!   2. ramp field uns->struct: per-face distinct values arrive in the same
!      order (index preservation through the payload buffer)
!   3. zero faces on the struct side (count=0 round, no deadlock, peer count ok)
!   4. zero faces on the uns side   (mirror of 3)
!   5. ramp field struct->uns: per-face SI values checked against hand-computed
!      expectations (conversion + ordering), stored for the echo round
!   6. unit round-trip over the wire: the SI state received in round 5 is
!      echoed back; the struct side converts SI->non-dim and compares with
!      the ramp it originally sent (nd -> SI -> MPI -> MPI -> SI -> nd)
!===============================================================================
program test_coupling_exchange
   use mpi
   use mod_precision, only: dp
   use mod_reference_state, only: init_reference_state, get_ref_state, &
                                  reference_state_t
   use mod_interface_units, only: SI_to_struct
   use mod_coupling_exchange, only: exchange_struct_to_uns, &
                                    exchange_uns_to_struct, recv_uns_iface_state
   implicit none

   integer, parameter :: NTEST = 6
   integer, parameter :: STRUCT_ROOT = 0, UNS_ROOT = 1
   integer :: ierr, rank, nproc, nfail, nfail_all

   ! echo storage: uns keeps the SI state received in round 4 and sends it
   ! back verbatim in round 5 (mimics the real coupling iteration offset)
   real(dp), allocatable :: echo_rho(:), echo_u(:,:), echo_T(:), echo_p(:)

   call MPI_Init( ierr )
   call MPI_Comm_rank( MPI_COMM_WORLD, rank, ierr )
   call MPI_Comm_size( MPI_COMM_WORLD, nproc, ierr )

   if ( nproc /= 2 ) then
      if ( rank == 0 ) write(*,'(a)') 'ERROR: this test requires mpirun -np 2'
      call MPI_Finalize( ierr )
      stop 1
   end if

   call init_reference_state()
   nfail = 0

   ! ---- round 1: tests 1 (const s->u) and 2 (ramp u->s) ----------------------
   if ( rank == STRUCT_ROOT ) then
      call round1_struct( nfail )
   else
      call round1_uns( nfail )
   end if

   ! ---- round 2: test 3 (zero struct faces) ----------------------------------
   if ( rank == STRUCT_ROOT ) then
      call round2_struct( nfail )
   else
      call round2_uns( nfail )
   end if

   ! ---- round 3: test 4 (zero uns faces) --------------------------------------
   if ( rank == STRUCT_ROOT ) then
      call round3_struct( nfail )
   else
      call round3_uns( nfail )
   end if

   ! ---- rounds 4-5: tests 5 (ramp s->u) and 6 (SI echo round-trip) ------------
   if ( rank == STRUCT_ROOT ) then
      call round4_struct( nfail )
      call round5_struct( nfail )
   else
      call round4_uns( nfail )
      call round5_uns( nfail )
   end if

   ! ---- verdict ----------------------------------------------------------------
   call MPI_Allreduce( nfail, nfail_all, 1, MPI_INTEGER, MPI_SUM, &
                       MPI_COMM_WORLD, ierr )
   if ( rank == 0 ) then
      write(*,'(a,i0,a,i0,a)') '=== ', NTEST, ' tests, ', nfail_all, ' failures ==='
      if ( nfail_all > 0 ) then
         write(*,'(a)') 'FAIL'
      else
         write(*,'(a)') 'PASS'
      end if
   end if
   call MPI_Finalize( ierr )
   if ( nfail_all > 0 ) stop 1

contains

   !---------------------------------------------------------------------------
   ! expect -- scalar comparison with combined abs/rel tolerance
   !---------------------------------------------------------------------------
   subroutine expect( got, want, name, nfail )
      real(dp),         intent(in)    :: got, want
      character(len=*), intent(in)    :: name
      integer,          intent(inout) :: nfail
      real(dp), parameter :: tol = 1.0e-12_dp
      if ( abs(got - want) > tol * (1.0_dp + abs(want)) ) then
         write(*,'(a,i0,a,a,a,es16.8,a,es16.8)') '[rank ', rank, '] FAIL ', &
               name, ': got ', got, ' want ', want
         nfail = nfail + 1
      end if
   end subroutine expect

   !---------------------------------------------------------------------------
   ! Round 1, struct side: send 3-face constant non-dim field; receive and
   ! check the 5-face uns ramp.
   !---------------------------------------------------------------------------
   subroutine round1_struct( nfail )
      integer, intent(inout) :: nfail
      integer, parameter :: NS = 3
      real(dp) :: rho_nd(NS), u_nd(NS), v_nd(NS), w_nd(NS), T_nd(NS), p_nd(NS)
      real(dp), allocatable :: urho(:), uu(:,:), uT(:), up(:)
      integer :: nu, i, ie

      rho_nd = 1.234_dp; u_nd =  0.111_dp; v_nd = -0.222_dp
      w_nd   = 0.333_dp; T_nd =  0.999_dp; p_nd =  0.00123_dp

      call exchange_struct_to_uns( rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                   NS, UNS_ROOT, ie )
      call expect( real(ie,dp), 0.0_dp, 'r1 struct send ierr', nfail )

      call recv_uns_iface_state( nu, urho, uu, uT, up, UNS_ROOT, ie )
      call expect( real(nu,dp), 5.0_dp, 'r1 struct recv nfaces', nfail )
      if ( nu == 5 ) then
         do i = 1, nu
            call expect( urho(i),  1.1_dp + 0.01_dp*i,            'r1 urho', nfail )
            call expect( uu(1,i),  10.0_dp + i,                   'r1 uu_x', nfail )
            call expect( uu(2,i), -20.0_dp - i,                   'r1 uu_y', nfail )
            call expect( uu(3,i),  30.0_dp + 2.0_dp*i,            'r1 uu_z', nfail )
            call expect( uT(i),    300.0_dp + i,                  'r1 uT',   nfail )
            call expect( up(i),    101325.0_dp + 100.0_dp*i,      'r1 up',   nfail )
         end do
      end if
   end subroutine round1_struct

   !---------------------------------------------------------------------------
   ! Round 1, uns side: receive the 3-face constant field and check against
   ! the expected SI conversion; send back a 5-face ramp.
   !---------------------------------------------------------------------------
   subroutine round1_uns( nfail )
      integer, intent(inout) :: nfail
      integer, parameter :: NU = 5
      real(dp) :: irho(NU), iu(3,NU), iT(NU), ip(NU)
      real(dp), allocatable :: srho(:), su(:,:), sT(:), sp(:)
      type(reference_state_t) :: r
      integer :: i, ie

      do i = 1, NU
         irho(i)   =  1.1_dp + 0.01_dp*i
         iu(:,i)   = [ 10.0_dp + i, -20.0_dp - i, 30.0_dp + 2.0_dp*i ]
         iT(i)     =  300.0_dp + i
         ip(i)     =  101325.0_dp + 100.0_dp*i
      end do

      call exchange_uns_to_struct( irho, iu, iT, ip, NU, STRUCT_ROOT, &
                                   srho, su, sT, sp, ie )
      call expect( real(ie,dp), 0.0_dp, 'r1 uns exch ierr', nfail )

      r = get_ref_state()
      call expect( real(size(srho),dp), 3.0_dp, 'r1 uns recv nfaces', nfail )
      if ( size(srho) == 3 ) then
         do i = 1, 3
            call expect( srho(i), 1.234_dp  * r%rho_ref, 'r1 srho', nfail )
            call expect( su(1,i), 0.111_dp  * r%a_ref,   'r1 su_x', nfail )
            call expect( su(2,i), -0.222_dp * r%a_ref,   'r1 su_y', nfail )
            call expect( su(3,i), 0.333_dp  * r%a_ref,   'r1 su_z', nfail )
            call expect( sT(i),   0.999_dp  * r%T_ref,   'r1 sT',   nfail )
            call expect( sp(i),   0.00123_dp * r%p_scale, 'r1 sp',  nfail )
         end do
      end if
   end subroutine round1_uns

   !---------------------------------------------------------------------------
   ! Round 2, struct side: send zero faces; receive and check a 4-face ramp.
   !---------------------------------------------------------------------------
   subroutine round2_struct( nfail )
      integer, intent(inout) :: nfail
      real(dp), allocatable :: empty(:)
      real(dp), allocatable :: urho(:), uu(:,:), uT(:), up(:)
      integer :: nu, i, ie

      allocate( empty(0) )
      call exchange_struct_to_uns( empty, empty, empty, empty, empty, empty, &
                                   0, UNS_ROOT, ie )
      call expect( real(ie,dp), 0.0_dp, 'r2 struct send0 ierr', nfail )

      call recv_uns_iface_state( nu, urho, uu, uT, up, UNS_ROOT, ie )
      call expect( real(nu,dp), 4.0_dp, 'r2 struct recv nfaces', nfail )
      if ( nu == 4 ) then
         do i = 1, nu
            call expect( urho(i), 2.0_dp + i,    'r2 urho', nfail )
            call expect( uT(i),   310.0_dp + i,  'r2 uT',   nfail )
         end do
      end if
   end subroutine round2_struct

   !---------------------------------------------------------------------------
   ! Round 2, uns side: receive zero struct faces (count must be 0); send a
   ! 4-face ramp.
   !---------------------------------------------------------------------------
   subroutine round2_uns( nfail )
      integer, intent(inout) :: nfail
      integer, parameter :: NU = 4
      real(dp) :: irho(NU), iu(3,NU), iT(NU), ip(NU)
      real(dp), allocatable :: srho(:), su(:,:), sT(:), sp(:)
      integer :: i, ie

      do i = 1, NU
         irho(i)   =  2.0_dp + i
         iu(:,i)   = [ 1.0_dp, 2.0_dp, 3.0_dp ] * i
         iT(i)     =  310.0_dp + i
         ip(i)     =  100000.0_dp + i
      end do

      call exchange_uns_to_struct( irho, iu, iT, ip, NU, STRUCT_ROOT, &
                                   srho, su, sT, sp, ie )
      call expect( real(ie,dp), 0.0_dp, 'r2 uns exch ierr', nfail )
      call expect( real(size(srho),dp), 0.0_dp, 'r2 uns recv nfaces=0', nfail )
   end subroutine round2_uns

   !---------------------------------------------------------------------------
   ! Round 3, struct side: send a 2-face constant variant; receive zero uns
   ! faces (count must be 0).
   !---------------------------------------------------------------------------
   subroutine round3_struct( nfail )
      integer, intent(inout) :: nfail
      integer, parameter :: NS = 2
      real(dp) :: rho_nd(NS), u_nd(NS), v_nd(NS), w_nd(NS), T_nd(NS), p_nd(NS)
      real(dp), allocatable :: urho(:), uu(:,:), uT(:), up(:)
      integer :: nu, ie

      rho_nd = 2.0_dp; u_nd = 0.5_dp; v_nd = 0.0_dp
      w_nd   = 0.0_dp;  T_nd = 1.1_dp; p_nd = 0.9_dp

      call exchange_struct_to_uns( rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                   NS, UNS_ROOT, ie )
      call expect( real(ie,dp), 0.0_dp, 'r3 struct send ierr', nfail )

      call recv_uns_iface_state( nu, urho, uu, uT, up, UNS_ROOT, ie )
      call expect( real(nu,dp), 0.0_dp, 'r3 struct recv nfaces=0', nfail )
   end subroutine round3_struct

   !---------------------------------------------------------------------------
   ! Round 3, uns side: receive the 2-face constant variant; send zero faces.
   !---------------------------------------------------------------------------
   subroutine round3_uns( nfail )
      integer, intent(inout) :: nfail
      real(dp), allocatable :: empty(:), empty3(:,:)
      real(dp), allocatable :: srho(:), su(:,:), sT(:), sp(:)
      type(reference_state_t) :: r
      integer :: i, ie

      allocate( empty(0), empty3(3,0) )
      call exchange_uns_to_struct( empty, empty3, empty, empty, 0, &
                                   STRUCT_ROOT, srho, su, sT, sp, ie )
      call expect( real(ie,dp), 0.0_dp, 'r3 uns exch ierr', nfail )

      r = get_ref_state()
      call expect( real(size(srho),dp), 2.0_dp, 'r3 uns recv nfaces', nfail )
      if ( size(srho) == 2 ) then
         do i = 1, 2
            call expect( srho(i), 2.0_dp * r%rho_ref,  'r3 srho', nfail )
            call expect( su(1,i), 0.5_dp * r%a_ref,    'r3 su_x', nfail )
            call expect( sT(i),   1.1_dp * r%T_ref,    'r3 sT',   nfail )
            call expect( sp(i),   0.9_dp * r%p_scale,  'r3 sp',   nfail )
         end do
      end if
   end subroutine round3_uns

   !---------------------------------------------------------------------------
   ! ramp_nd -- the non-dim ramp sent by struct in round 4 (shared definition)
   !---------------------------------------------------------------------------
   pure subroutine ramp_nd( i, rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd )
      integer,  intent(in)  :: i
      real(dp), intent(out) :: rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd
      rho_nd = 1.0_dp + 0.1_dp  * i
      u_nd   = 0.2_dp + 0.01_dp * i
      v_nd   = -0.1_dp - 0.02_dp * i
      w_nd   = 0.05_dp * i
      T_nd   = 0.9_dp + 0.01_dp * i
      p_nd   = 0.8_dp + 0.05_dp * i
   end subroutine ramp_nd

   !---------------------------------------------------------------------------
   ! Round 4, struct side: send a 4-face non-dim ramp; receive the dummy
   ! 1-face uns payload (protocol symmetry only, values unchecked).
   !---------------------------------------------------------------------------
   subroutine round4_struct( nfail )
      integer, intent(inout) :: nfail
      integer, parameter :: NS = 4
      real(dp) :: rho_nd(NS), u_nd(NS), v_nd(NS), w_nd(NS), T_nd(NS), p_nd(NS)
      real(dp), allocatable :: urho(:), uu(:,:), uT(:), up(:)
      integer :: nu, i, ie

      do i = 1, NS
         call ramp_nd( i, rho_nd(i), u_nd(i), v_nd(i), w_nd(i), T_nd(i), p_nd(i) )
      end do

      call exchange_struct_to_uns( rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                   NS, UNS_ROOT, ie )
      call expect( real(ie,dp), 0.0_dp, 'r4 struct send ierr', nfail )

      call recv_uns_iface_state( nu, urho, uu, uT, up, UNS_ROOT, ie )
      call expect( real(nu,dp), 1.0_dp, 'r4 struct recv nfaces', nfail )
   end subroutine round4_struct

   !---------------------------------------------------------------------------
   ! Round 4, uns side: receive the 4-face ramp, check each SI value against
   ! hand-computed expectations, and store it for the round-5 echo.  Send a
   ! dummy 1-face payload back.
   !---------------------------------------------------------------------------
   subroutine round4_uns( nfail )
      integer, intent(inout) :: nfail
      real(dp) :: irho(1), iu(3,1), iT(1), ip(1)
      type(reference_state_t) :: r
      real(dp) :: r_nd, u_nd, v_nd, w_nd, T_nd, p_nd
      integer :: i, ie

      irho = 1.0_dp; iu = 0.0_dp; iT = 300.0_dp; ip = 101325.0_dp

      call exchange_uns_to_struct( irho, iu, iT, ip, 1, STRUCT_ROOT, &
                                   echo_rho, echo_u, echo_T, echo_p, ie )
      call expect( real(ie,dp), 0.0_dp, 'r4 uns exch ierr', nfail )

      r = get_ref_state()
      call expect( real(size(echo_rho),dp), 4.0_dp, 'r4 uns recv nfaces', nfail )
      if ( size(echo_rho) == 4 ) then
         do i = 1, 4
            call ramp_nd( i, r_nd, u_nd, v_nd, w_nd, T_nd, p_nd )
            call expect( echo_rho(i),   r_nd * r%rho_ref, 'r4 srho', nfail )
            call expect( echo_u(1,i),   u_nd * r%a_ref,   'r4 su_x', nfail )
            call expect( echo_u(2,i),   v_nd * r%a_ref,   'r4 su_y', nfail )
            call expect( echo_u(3,i),   w_nd * r%a_ref,   'r4 su_z', nfail )
            call expect( echo_T(i),     T_nd * r%T_ref,   'r4 sT',   nfail )
            call expect( echo_p(i),     p_nd * r%p_scale, 'r4 sp',   nfail )
         end do
      end if
   end subroutine round4_uns

   !---------------------------------------------------------------------------
   ! Round 5, struct side: send a dummy 1-face payload; receive the echoed SI
   ! ramp, convert back to non-dim and compare with the original round-4 ramp.
   !---------------------------------------------------------------------------
   subroutine round5_struct( nfail )
      integer, intent(inout) :: nfail
      real(dp) :: rho_nd(1), u_nd(1), v_nd(1), w_nd(1), T_nd(1), p_nd(1)
      real(dp), allocatable :: urho(:), uu(:,:), uT(:), up(:)
      real(dp) :: r_nd, u_nd0, v_nd0, w_nd0, T_nd0, p_nd0
      real(dp) :: r2, u2, v2, w2, T2, p2
      integer :: nu, i, ie

      rho_nd = 1.0_dp; u_nd = 0.0_dp; v_nd = 0.0_dp
      w_nd   = 0.0_dp;  T_nd = 1.0_dp; p_nd = 1.0_dp

      call exchange_struct_to_uns( rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                   1, UNS_ROOT, ie )
      call expect( real(ie,dp), 0.0_dp, 'r5 struct send ierr', nfail )

      call recv_uns_iface_state( nu, urho, uu, uT, up, UNS_ROOT, ie )
      call expect( real(nu,dp), 4.0_dp, 'r5 struct recv nfaces', nfail )
      if ( nu == 4 ) then
         do i = 1, nu
            call SI_to_struct( urho(i), uu(1,i), uu(2,i), uu(3,i), uT(i), up(i), &
                               r2, u2, v2, w2, T2, p2 )
            call ramp_nd( i, r_nd, u_nd0, v_nd0, w_nd0, T_nd0, p_nd0 )
            call expect( r2, r_nd,  'r5 roundtrip rho', nfail )
            call expect( u2, u_nd0, 'r5 roundtrip u',   nfail )
            call expect( v2, v_nd0, 'r5 roundtrip v',   nfail )
            call expect( w2, w_nd0, 'r5 roundtrip w',   nfail )
            call expect( T2, T_nd0, 'r5 roundtrip T',   nfail )
            call expect( p2, p_nd0, 'r5 roundtrip p',   nfail )
         end do
      end if
   end subroutine round5_struct

   !---------------------------------------------------------------------------
   ! Round 5, uns side: echo the round-4 SI state back verbatim; receive the
   ! dummy 1-face struct payload (unchecked).
   !---------------------------------------------------------------------------
   subroutine round5_uns( nfail )
      integer, intent(inout) :: nfail
      real(dp), allocatable :: srho(:), su(:,:), sT(:), sp(:)
      integer :: ie

      call exchange_uns_to_struct( echo_rho, echo_u, echo_T, echo_p, &
                                   size(echo_rho), STRUCT_ROOT, &
                                   srho, su, sT, sp, ie )
      call expect( real(ie,dp), 0.0_dp, 'r5 uns exch ierr', nfail )
      call expect( real(size(srho),dp), 1.0_dp, 'r5 uns recv nfaces', nfail )
   end subroutine round5_uns

end program test_coupling_exchange
