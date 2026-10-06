!===============================================================================
! test_units_exchange.f90 -- phase-5 unit tests:
!   (1) struct <-> SI unit conversion is its own inverse (round-trip)
!   (2) struct->uns bilinear interpolation is exact for constant fields
!   (3) uns->struct area-weighted average preserves total flux
!===============================================================================
program test_units_exchange
   use mod_precision, only: dp
   use mod_reference_state, only: init_reference_state, get_ref_state
   use mod_interface_units, only: struct_to_SI, SI_to_struct
   use mod_interface, only: Interface_FACE_TYPE, Interface_List, Num_Interface, &
                            PEER_STRUCT, PEER_UNS, MATCH_MATCHED
   use mod_interface_exchange, only: iface_state_t, alloc_iface_state, &
                                     struct_to_uns_interpolate, &
                                     uns_to_struct_interpolate
   implicit none

   integer, parameter :: NTEST = 3
   integer :: nfail
   nfail = 0

   ! ---- test 1: unit conversion round-trip ----
   call test_conversion_roundtrip( nfail )

   ! ---- tests 2 & 3: interpolation on a synthetic 1x2 interface ----
   call test_interpolation( nfail )

   write(*,'(a,i0,a,i0,a)') '=== ', NTEST, ' tests, ', nfail, ' failures ==='
   if ( nfail > 0 ) then
      write(*,'(a)') 'FAIL'
      stop 1
   end if
   write(*,'(a)') 'PASS'

contains

   subroutine test_conversion_roundtrip( nfail )
      integer, intent(inout) :: nfail
      real(dp) :: rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd
      real(dp) :: rho, u, v, w, T, p
      real(dp) :: rho2, u2, v2, w2, T2, p2
      real(dp) :: tol = 1.0e-14_dp

      call init_reference_state()

      rho_nd = 1.234_dp; u_nd = 0.111_dp; v_nd = -0.222_dp
      w_nd   = 0.333_dp; T_nd = 0.999_dp; p_nd = 0.00123_dp

      call struct_to_SI( rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                         rho, u, v, w, T, p )
      call SI_to_struct( rho, u, v, w, T, p, &
                         rho2, u2, v2, w2, T2, p2 )

      if ( abs(rho2-rho_nd) > tol ) nfail = nfail + 1
      if ( abs(u2-u_nd)     > tol ) nfail = nfail + 1
      if ( abs(v2-v_nd)     > tol ) nfail = nfail + 1
      if ( abs(w2-w_nd)     > tol ) nfail = nfail + 1
      if ( abs(T2-T_nd)     > tol ) nfail = nfail + 1
      if ( abs(p2-p_nd)     > tol ) nfail = nfail + 1
      if ( nfail > 0 ) then
         write(*,'(a)') 'FAIL: conversion round-trip'
      else
         write(*,'(a)') 'PASS: conversion round-trip (tol 1e-14)'
      end if
   end subroutine test_conversion_roundtrip

   !---------------------------------------------------------------------------
   ! Synthetic interface: 1 structured quad (global idx 1) and 2 unstructured
   ! faces (global idx 2,3) both pointing at the struct face.  The struct quad
   ! has 4 vertices; uns face 2 has bilinear weights giving vertex A, uns face 3
   ! giving vertex C.  Constant field => both interpolations reproduce it.
   !---------------------------------------------------------------------------
   subroutine test_interpolation( nfail )
      integer, intent(inout) :: nfail
      type(iface_state_t) :: st
      real(dp) :: tol = 1.0e-14_dp
      real(dp) :: const_rho = 1.2_dp, const_u = 10.0_dp, const_T = 300.0_dp, &
                  const_p = 1.0e5_dp

      ! build Interface_List manually
      Num_Interface = 3
      allocate( Interface_List(3) )

      ! struct face (idx 1)
      Interface_List(1)%solver = PEER_STRUCT
      Interface_List(1)%nv = 4
      Interface_List(1)%match_state = MATCH_MATCHED
      Interface_List(1)%peer_id = 2    ! points to first uns face (any is fine)
      Interface_List(1)%area = 1.0_dp

      ! uns face 2 (idx 2): weights = vertex A of struct quad
      Interface_List(2)%solver = PEER_UNS
      Interface_List(2)%match_state = MATCH_MATCHED
      Interface_List(2)%peer_id = 1
      allocate( Interface_List(2)%peer_w(4) )
      Interface_List(2)%peer_w = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      Interface_List(2)%area = 0.5_dp

      ! uns face 3 (idx 3): weights = vertex C of struct quad
      Interface_List(3)%solver = PEER_UNS
      Interface_List(3)%match_state = MATCH_MATCHED
      Interface_List(3)%peer_id = 1
      allocate( Interface_List(3)%peer_w(4) )
      Interface_List(3)%peer_w = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
      Interface_List(3)%area = 0.5_dp

      call alloc_iface_state( st )

      ! set struct vertex state to a constant field
      st%v_rho(1,:) = const_rho
      st%v_T(1,:)   = const_T
      st%v_p(1,:)   = const_p
      st%v_u(1,1,:) = const_u
      st%v_u(1,2,:) = 0.0_dp
      st%v_u(1,3,:) = 0.0_dp

      ! struct -> uns
      call struct_to_uns_interpolate( st )

      if ( abs(st%rho(2)-const_rho) > tol .or. &
           abs(st%rho(3)-const_rho) > tol .or. &
           abs(st%u(1,2)-const_u)  > tol .or. &
           abs(st%T(2)-const_T)    > tol .or. &
           abs(st%p(2)-const_p)    > tol ) then
         nfail = nfail + 1
         write(*,'(a)') 'FAIL: struct->uns constant interpolation'
      else
         write(*,'(a)') 'PASS: struct->uns constant interpolation'
      end if

      ! now perturb uns faces and check uns->struct area-weighted average
      st%rho(2) = 2.0_dp; st%rho(3) = 4.0_dp
      st%u(1,2) = 0.0_dp; st%u(1,3) = 2.0_dp
      call uns_to_struct_interpolate( st )
      ! areas are 0.5 and 0.5 -> average = (0.5*2 + 0.5*4)/1.0 = 3.0
      if ( abs(st%rho(1)-3.0_dp) > tol .or. &
           abs(st%u(1,1)-1.0_dp) > tol ) then
         nfail = nfail + 1
         write(*,'(a,es12.4,a,es12.4)') 'FAIL: uns->struct area-weighted avg; got rho=', &
               st%rho(1), ' u=', st%u(1,1)
      else
         write(*,'(a)') 'PASS: uns->struct area-weighted average'
      end if
   end subroutine test_interpolation

end program test_units_exchange
