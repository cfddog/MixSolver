!===============================================================================
! mod_coupling_exchange.f90 -- MPI cross-group data exchange (phase 6)
!
! Orchestrates the weak-coupling data exchange between the structured and
! unstructured solver groups.  Both groups share MPI_COMM_WORLD; rank 0 of
! each group acts as the "coupling root" for its side.
!
! Exchange protocol (per coupling iteration):
!   Phase 1 (struct -> uns):
!     struct root sends: (a) nfaces_struct (int)
!                        (b) buf(7, nfaces_struct) SI state at struct cells +
!                            struct-side wall heat flux q_f [W/m^2] in row 7
!                            (phase 17: positive = heat entering the peer; the
!                            solid side applies it as a Neumann surface flux)
!     uns root receives and stores as "struct SI state at struct faces"
!     (the interpolation struct->uns happens in the coupling driver after
!     receive, using the phase-5 mod_interface_exchange routines).
!
!   Phase 2 (uns -> struct):
!     uns root sends:    (a) nfaces_uns (int)
!                        (b) buf(7, nfaces_uns) SI state at uns faces +
!                            per-face cell-zone class in row 7
!                            (IFACE_CZ_FLUID / IFACE_CZ_POROUS / IFACE_CZ_SOLID,
!                            phase 17: lets the struct root dispatch the peer
!                            type without re-deriving it from geometry)
!     struct root receives and stores as "uns SI state at uns faces"
!     (the interpolation uns->struct happens in the coupling driver).
!
! The two-phase send-then-send pattern with matched recv avoids deadlock
! because the two groups are disjoint: each side sends first, then receives.
! Blocking Send is used; MPI buffer (attached at MPI_Init) absorbs the
! message so both sides can proceed to the Recv without waiting.
!
! NOTE on count=0: when one side has no interface faces (nfaces=0) it still
! sends the count and a zero-length payload so the peer's Recv matches.
!===============================================================================
module mod_coupling_exchange
   use mpi
   use mod_precision, only: dp
   use mod_interface_units, only: struct_to_SI, SI_to_struct
   use mod_reference_state, only: reference_state_t, get_ref_state
   implicit none
   private

   public :: exchange_struct_to_uns
   public :: exchange_uns_to_struct
   public :: recv_uns_iface_state

   ! MPI tags
   integer, parameter :: TAG_CNT_S2U = 110, TAG_DAT_S2U = 111
   integer, parameter :: TAG_CNT_U2S = 210, TAG_DAT_U2S = 211

contains

   !---------------------------------------------------------------------------
   ! exchange_struct_to_uns -- struct root: convert non-dim interface state to
   ! SI and send (count + payload) to the uns root.
   !
   ! nfaces: number of struct interface faces (may be 0).
   ! rho_nd/u_nd/v_nd/w_nd/T_nd/p_nd: non-dim state at struct interface cells
   !   (dummy args may be unallocated when nfaces==0; only size 0 is sent).
   ! q_nd: struct-side wall heat flux at the interface faces, in SI [W/m^2]
   !   (phase 17: already SI, no conversion; positive = heat entering the peer
   !   solid).  May be unallocated / size 0 when nfaces==0.
   ! uns_root_rank: global rank of the uns coupling root.
   !---------------------------------------------------------------------------
   subroutine exchange_struct_to_uns(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                      q_nd, nfaces, uns_root_rank, ierr)
      integer,  intent(in)  :: nfaces, uns_root_rank
      real(dp), intent(in)  :: rho_nd(*), u_nd(*), v_nd(*), w_nd(*), &
                                T_nd(*), p_nd(*), q_nd(*)
      integer,  intent(out) :: ierr

      real(dp), allocatable :: buf(:,:)
      integer :: i
      type(reference_state_t) :: r

      ! ---- send the face count first ----
      call MPI_Send(nfaces, 1, MPI_INTEGER, uns_root_rank, TAG_CNT_S2U, &
                    MPI_COMM_WORLD, ierr)

      if (nfaces <= 0) return

      ! ---- convert to SI and send payload ----
      r = get_ref_state()
      allocate(buf(7, nfaces))
      do i = 1, nfaces
         call struct_to_SI(rho_nd(i), u_nd(i), v_nd(i), w_nd(i), &
                           T_nd(i), p_nd(i), &
                           buf(1,i), buf(2,i), buf(3,i), buf(4,i), &
                           buf(5,i), buf(6,i), r)
         buf(7,i) = q_nd(i)          ! already SI [W/m^2]
      end do
      call MPI_Send(buf, 7*nfaces, MPI_DOUBLE_PRECISION, uns_root_rank, &
                    TAG_DAT_S2U, MPI_COMM_WORLD, ierr)
      deallocate(buf)
   end subroutine exchange_struct_to_uns

   !---------------------------------------------------------------------------
   ! exchange_uns_to_struct -- uns root:
   !   1. receive struct SI state (count + payload)
   !   2. send uns interface SI state + per-face cell-zone class
   !      (count + 7-row payload) back to struct root
   !
   ! iface_rho/iface_u/iface_T/iface_p: uns SI state at uns interface faces.
   ! iface_cz: per-face cell-zone class (IFACE_CZ_*).
   ! n_uns_faces: number of uns interface faces (may be 0).
   ! struct_root_rank: global rank of the struct coupling root.
   !
   ! On return, srho/su/sT/sp hold the struct SI state at struct faces
   ! (size = nfaces_struct received from the struct side; may be 0), and sq
   ! holds the struct-side wall heat flux q_f [W/m^2] (payload row 7).
   !---------------------------------------------------------------------------
   subroutine exchange_uns_to_struct(iface_rho, iface_u, iface_T, iface_p, &
                                      iface_cz, n_uns_faces, struct_root_rank, &
                                      srho, su, sT, sp, sq, ierr)
      integer,  intent(in)  :: n_uns_faces, struct_root_rank
      real(dp), intent(in)  :: iface_rho(*), iface_u(3,*), iface_T(*), iface_p(*)
      integer,  intent(in)  :: iface_cz(*)
      real(dp), allocatable, intent(out) :: srho(:), su(:,:), sT(:), sp(:)
      real(dp), allocatable, intent(out) :: sq(:)
      integer,  intent(out) :: ierr

      real(dp), allocatable :: buf(:,:)
      integer :: i, n_struct_faces

      ! ---- Phase 1: receive struct face count ----
      call MPI_Recv(n_struct_faces, 1, MPI_INTEGER, struct_root_rank, &
                    TAG_CNT_S2U, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      allocate(srho(n_struct_faces), su(3,n_struct_faces), &
               sT(n_struct_faces), sp(n_struct_faces), sq(n_struct_faces))

      if (n_struct_faces > 0) then
         allocate(buf(7, n_struct_faces))
         call MPI_Recv(buf, 7*n_struct_faces, MPI_DOUBLE_PRECISION, &
                       struct_root_rank, TAG_DAT_S2U, MPI_COMM_WORLD, &
                       MPI_STATUS_IGNORE, ierr)
         do i = 1, n_struct_faces
            srho(i)   = buf(1,i)
            su(:,i)   = buf(2:4,i)
            sT(i)     = buf(5,i)
            sp(i)     = buf(6,i)
            sq(i)     = buf(7,i)     ! struct-side wall heat flux [W/m^2]
         end do
         deallocate(buf)
      end if

      ! ---- Phase 2: send uns face count + payload (7 rows) to struct root ----
      call MPI_Send(n_uns_faces, 1, MPI_INTEGER, struct_root_rank, &
                    TAG_CNT_U2S, MPI_COMM_WORLD, ierr)

      if (n_uns_faces > 0) then
         allocate(buf(7, n_uns_faces))
         do i = 1, n_uns_faces
            buf(1,i)   = iface_rho(i)
            buf(2:4,i) = iface_u(:,i)
            buf(5,i)   = iface_T(i)
            buf(6,i)   = iface_p(i)
            buf(7,i)   = real(iface_cz(i), dp)
         end do
         call MPI_Send(buf, 7*n_uns_faces, MPI_DOUBLE_PRECISION, &
                       struct_root_rank, TAG_DAT_U2S, MPI_COMM_WORLD, ierr)
         deallocate(buf)
      end if
   end subroutine exchange_uns_to_struct

   !---------------------------------------------------------------------------
   ! recv_uns_iface_state -- struct root: receive the uns interface SI state.
   !
   ! n_uns_faces: set to the received uns face count.
   ! urho/uu/uT/up: uns SI state at uns faces (allocated on return).
   ! ucz: per-face cell-zone class received in payload row 7 (phase 17).
   !---------------------------------------------------------------------------
   subroutine recv_uns_iface_state(n_uns_faces, urho, uu, uT, up, ucz, &
                                    uns_root_rank, ierr)
      integer,  intent(out) :: n_uns_faces, ierr
      integer,  intent(in)  :: uns_root_rank
      real(dp), allocatable, intent(out) :: urho(:), uu(:,:), uT(:), up(:)
      integer,  allocatable, intent(out) :: ucz(:)

      real(dp), allocatable :: buf(:,:)
      integer :: i

      call MPI_Recv(n_uns_faces, 1, MPI_INTEGER, uns_root_rank, TAG_CNT_U2S, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      allocate(urho(n_uns_faces), uu(3,n_uns_faces), &
               uT(n_uns_faces), up(n_uns_faces), ucz(n_uns_faces))

      if (n_uns_faces > 0) then
         allocate(buf(7, n_uns_faces))
         call MPI_Recv(buf, 7*n_uns_faces, MPI_DOUBLE_PRECISION, &
                       uns_root_rank, TAG_DAT_U2S, MPI_COMM_WORLD, &
                       MPI_STATUS_IGNORE, ierr)
         do i = 1, n_uns_faces
            urho(i)   = buf(1,i)
            uu(:,i)   = buf(2:4,i)
            uT(i)     = buf(5,i)
            up(i)     = buf(6,i)
            ucz(i)    = nint(buf(7,i))
         end do
         deallocate(buf)
      end if
   end subroutine recv_uns_iface_state

end module mod_coupling_exchange
