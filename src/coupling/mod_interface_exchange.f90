!===============================================================================
! mod_interface_exchange.f90 -- conservative interface data exchange (phase 5).
!
! Builds on phase-4 geometric matching (mod_interface_match), which stored
! peer_id / peer_w on each matched interface face in mod_interface's
! Interface_List.  This module owns the *flow state* on those faces (SI units)
! and performs the two directional interpolations:
!
!   STRUCT -> UNS : each unstructured face centroid state is the bilinear
!                   interpolation of its peer structured quad's 4 vertex
!                   states, using the weights peer_w (sum = 1 by construction).
!
!   UNS -> STRUCT : each structured face state is the area-weighted average
!                   of all unstructured faces whose peer_id points here.
!                   (Conservative: sum(area_uns * q_uns) / sum(area_uns).)
!
! Both directions operate on SI quantities (mod_interface_units handles the
! non-dimensional <-> SI conversion on each solver side; this module stays
! unit-agnostic).
!
! MPI architecture (phase 5 skeleton, exercised in phase 6):
!   - struct and uns solvers share MPI_COMM_WORLD; split into struct/uns
!     sub-communicators.
!   - The interface state is assembled on a "coupling root" rank (rank 0 of
!     each sub-comm), interpolated, then scattered.  For the serial build the
!     exchange routines operate directly on the in-memory state.
!
! The interpolations are pure array routines so they are unit-testable without
! MPI; the MPI wrappers are thin.
!===============================================================================
module mod_interface_exchange
   use mod_precision, only: dp
   use mod_interface, only: Interface_FACE_TYPE, Interface_List, Num_Interface, &
                            PEER_STRUCT, PEER_UNS, MATCH_MATCHED
   implicit none
   private

   public :: iface_state_t
   public :: alloc_iface_state, free_iface_state
   public :: struct_to_uns_interpolate
   public :: uns_to_struct_interpolate

   ! ---------------------------------------------------------------------------
   ! iface_state_t -- SI flow state on interface faces, parallel to
   ! Interface_List(1:Num_Interface).
   !
   !   rho(:), u(3,:), T(:), p(:)   -- per-face state (both solvers)
   !   v_rho/ v_u / v_T / v_p       -- per-struct-face vertex state, indexed
   !                                    by the struct-face local index in
   !                                    sidx (:) (see alloc_iface_state).
   ! ---------------------------------------------------------------------------
   type :: iface_state_t
      real(dp), allocatable :: rho(:), T(:), p(:)
      real(dp), allocatable :: u(:,:)          ! (3, Num_Interface)
      ! struct-vertex state (only faces with solver==PEER_STRUCT)
      integer,  allocatable :: sidx(:)         ! local->global struct face map
      real(dp), allocatable :: v_rho(:,:)      ! (ns, 4)
      real(dp), allocatable :: v_u(:,:,:)      ! (ns, 3, 4)
      real(dp), allocatable :: v_T(:,:)        ! (ns, 4)
      real(dp), allocatable :: v_p(:,:)        ! (ns, 4)
   end type iface_state_t

contains

   !---------------------------------------------------------------------------
   ! Allocate state arrays sized to the current Interface_List.  Struct-vertex
   ! arrays are sized to the number of PEER_STRUCT faces.
   !---------------------------------------------------------------------------
   subroutine alloc_iface_state( st )
      type(iface_state_t), intent(out) :: st
      integer :: i, ns
      if ( allocated(st%rho) ) call free_iface_state( st )
      allocate( st%rho(Num_Interface), st%T(Num_Interface), &
                st%p(Num_Interface), st%u(3, Num_Interface) )
      st%rho = 0.0_dp; st%T = 0.0_dp; st%p = 0.0_dp; st%u = 0.0_dp
      ns = count( Interface_List(1:Num_Interface)%solver == PEER_STRUCT )
      allocate( st%sidx(ns) )
      ns = 0
      do i = 1, Num_Interface
         if ( Interface_List(i)%solver == PEER_STRUCT ) then
            ns = ns + 1
            st%sidx(ns) = i
         end if
      end do
      allocate( st%v_rho(ns,4), st%v_u(ns,3,4), st%v_T(ns,4), st%v_p(ns,4) )
      st%v_rho = 0.0_dp; st%v_u = 0.0_dp
      st%v_T = 0.0_dp;   st%v_p = 0.0_dp
   end subroutine alloc_iface_state

   subroutine free_iface_state( st )
      type(iface_state_t), intent(inout) :: st
      if ( allocated(st%rho) )  deallocate( st%rho )
      if ( allocated(st%T) )    deallocate( st%T )
      if ( allocated(st%p) )    deallocate( st%p )
      if ( allocated(st%u) )    deallocate( st%u )
      if ( allocated(st%sidx) ) deallocate( st%sidx )
      if ( allocated(st%v_rho) ) deallocate( st%v_rho )
      if ( allocated(st%v_u) )   deallocate( st%v_u )
      if ( allocated(st%v_T) )   deallocate( st%v_T )
      if ( allocated(st%v_p) )   deallocate( st%v_p )
   end subroutine free_iface_state

   !---------------------------------------------------------------------------
   ! STRUCT -> UNS interpolation.
   ! For every matched unstructured face, compute its SI state as the bilinear
   ! weighted sum of its peer structured quad's 4 vertex states:
   !     q_uns = sum_k  peer_w(k) * q_struct_vertex(k)
   ! Struct vertex states are read from st%v_* ; results written to st%rho/u/T/p
   ! at the unstructured face's global index.
   !---------------------------------------------------------------------------
   subroutine struct_to_uns_interpolate( st )
      type(iface_state_t), intent(inout) :: st
      integer :: j, k, istruct, peer
      real(dp) :: w(4)
      do j = 1, Num_Interface
         if ( Interface_List(j)%solver /= PEER_UNS ) cycle
         if ( Interface_List(j)%match_state /= MATCH_MATCHED ) cycle
         peer = Interface_List(j)%peer_id          ! global struct face index
         ! find local struct index in sidx
         istruct = 0
         do k = 1, size(st%sidx)
            if ( st%sidx(k) == peer ) then
               istruct = k
               exit
            end if
         end do
         if ( istruct == 0 ) cycle
         w = Interface_List(j)%peer_w(1:4)
         st%rho(j) = dot_product( w, st%v_rho(istruct,:) )
         st%T(j)   = dot_product( w, st%v_T(istruct,:) )
         st%p(j)   = dot_product( w, st%v_p(istruct,:) )
         do k = 1, 3
            st%u(k,j) = dot_product( w, st%v_u(istruct,k,:) )
         end do
      end do
   end subroutine struct_to_uns_interpolate

   !---------------------------------------------------------------------------
   ! UNS -> STRUCT interpolation (area-weighted, conservative).
   ! For every matched structured face, average the states of all unstructured
   ! faces whose peer_id points here, weighted by face area:
   !     q_struct = sum( area_uns * q_uns ) / sum( area_uns )
   ! Result written to st%rho/u/T/p at the structured face's global index.
   !---------------------------------------------------------------------------
   subroutine uns_to_struct_interpolate( st )
      type(iface_state_t), intent(inout) :: st
      integer :: i, j
      real(dp) :: wsum, qsum_rho, qsum_T, qsum_p, qsum_u(3), ai
      do i = 1, Num_Interface
         if ( Interface_List(i)%solver /= PEER_STRUCT ) cycle
         if ( Interface_List(i)%match_state /= MATCH_MATCHED ) cycle
         wsum = 0.0_dp
         qsum_rho = 0.0_dp; qsum_T = 0.0_dp; qsum_p = 0.0_dp
         qsum_u = 0.0_dp
         do j = 1, Num_Interface
            if ( Interface_List(j)%solver /= PEER_UNS ) cycle
            if ( Interface_List(j)%match_state /= MATCH_MATCHED ) cycle
            if ( Interface_List(j)%peer_id /= i ) cycle
            ai = Interface_List(j)%area
            wsum     = wsum + ai
            qsum_rho = qsum_rho + ai * st%rho(j)
            qsum_T   = qsum_T   + ai * st%T(j)
            qsum_p   = qsum_p   + ai * st%p(j)
            qsum_u   = qsum_u   + ai * st%u(:,j)
         end do
         if ( wsum > 0.0_dp ) then
            st%rho(i) = qsum_rho / wsum
            st%T(i)   = qsum_T   / wsum
            st%p(i)   = qsum_p   / wsum
            st%u(:,i) = qsum_u   / wsum
         end if
      end do
   end subroutine uns_to_struct_interpolate

end module mod_interface_exchange
