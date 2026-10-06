!===============================================================================
! mod_halo.f90 -- Halo (ghost-cell) communication for parallel UNSSolver
!                  (Step 4)
!
! halo_setup builds, for the local domain described by local_mesh_t, the
! symmetric neighbour schedule:
!   - recv side: which of MY halo cells I need, and from which owner rank
!     (discovered locally from the global partition vector)
!   - send side: which of MY owned cells other ranks need (discovered by a
!     two-phase negotiation: exchange request counts with MPI_Alltoall, then
!     exchange the requested global cell IDs with MPI_Alltoallv)
!
! halo_exchange_scalar / halo_exchange_vector then update the halo slots of a
! field in one non-blocking exchange per call (Irecv + Isend + Waitall).
!
! Storage inside halo_info_t is flat: per-peer ids are concatenated in peers()
! order with per-peer counts and offsets, so message buffers are contiguous.
!===============================================================================
module mod_uns_halo
   use mod_precision, only: dp, ip, pi
   use mod_uns_mpi_core, only: myrank, nprocs, mpi_check
   use mpi
   use mod_uns_local_mesh, only: local_mesh_t
   implicit none
   private
   public :: halo_info_t, halo_setup, halo_exchange_scalar, halo_exchange_vector

   type :: halo_info_t
      integer :: npeers = 0
      integer, allocatable :: peers(:)    ! neighbour ranks (size npeers)

      ! send schedule (my owned cells requested by each peer)
      integer, allocatable :: scnt(:)    ! cells sent per peer (size npeers)
      integer, allocatable :: soff(:)    ! flat offsets, 0-based (npeers+1)
      integer, allocatable :: sids(:)    ! local owned cell IDs (total send)

      ! recv schedule (my halo cells owned by each peer)
      integer, allocatable :: rcnt(:)
      integer, allocatable :: roff(:)
      integer, allocatable :: rids(:)    ! local halo cell IDs (total recv)
   end type halo_info_t

   integer, parameter :: HALO_TAG = 200

contains

   !----------------------------------------------------------------------------
   ! Build the halo schedule. part(:) is the global partition vector.
   !----------------------------------------------------------------------------
   subroutine halo_setup( lm, part, hi, ierr )
      type(local_mesh_t), intent(in)  :: lm
      integer,            intent(in)  :: part(:)
      type(halo_info_t),  intent(out) :: hi
      integer,            intent(out) :: ierr

      integer  :: npr, r, j, h, g, total_in, p
      integer, allocatable :: want_cnt(:)     ! cells I want from each rank
      integer, allocatable :: send_cnt(:)     ! cells each rank wants from me
      integer, allocatable :: want_global(:)  ! global IDs I want (grouped)
      integer, allocatable :: want_cursor(:)
      integer, allocatable :: in_global(:)    ! global IDs requested of me
      integer, allocatable :: sdispl(:), rdispl(:)

      ierr = 0
      npr = nprocs

      ! ---- 1. recv side: group halo cells by owner rank ---------------------
      allocate( want_cnt(0:npr-1), want_cursor(0:npr-1) )
      want_cnt = 0
      do h = 1, lm%nghost
         r = part(lm%halo_global(h))
         want_cnt(r) = want_cnt(r) + 1
      end do
      allocate( want_global(sum(want_cnt)) )
      ! offsets per rank
      want_cursor = 0
      do r = 1, npr-1
         want_cursor(r) = want_cursor(r-1) + want_cnt(r-1)
      end do
      ! fill grouped requested global IDs in halo order
      block
         integer, allocatable :: cur(:)
         allocate( cur(0:npr-1) ); cur = want_cursor
         do h = 1, lm%nghost
            r = part(lm%halo_global(h))
            cur(r) = cur(r) + 1
            want_global(cur(r)) = lm%halo_global(h)
         end do
         deallocate( cur )
      end block

      ! ---- 2. exchange request counts --------------------------------------
      allocate( send_cnt(0:npr-1) )
      call MPI_Alltoall( want_cnt, 1, MPI_INTEGER, &
                         send_cnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
      call mpi_check( ierr, 'halo_setup alltoall counts' )

      ! ---- 3. exchange requested global IDs --------------------------------
      allocate( sdispl(0:npr), rdispl(0:npr) )
      sdispl(0) = 0; rdispl(0) = 0
      do r = 0, npr-1
         sdispl(r+1) = sdispl(r) + want_cnt(r)
         rdispl(r+1) = rdispl(r) + send_cnt(r)
      end do
      total_in = rdispl(npr)
      allocate( in_global(total_in) )
      call MPI_Alltoallv( want_global, want_cnt, sdispl, MPI_INTEGER, &
                          in_global,   send_cnt, rdispl, MPI_INTEGER, &
                          MPI_COMM_WORLD, ierr )
      call mpi_check( ierr, 'halo_setup alltoallv ids' )

      ! ---- 4. compact schedule, dropping zero pairs ------------------------
      hi%npeers = count( want_cnt(0:npr-1) + send_cnt(0:npr-1) > 0 )
      allocate( hi%peers(hi%npeers), hi%scnt(hi%npeers), hi%rcnt(hi%npeers), &
                hi%soff(hi%npeers+1), hi%roff(hi%npeers+1) )
      j = 0
      do r = 0, npr-1
         if ( want_cnt(r) + send_cnt(r) == 0 ) cycle
         j = j + 1
         hi%peers(j) = r
         hi%scnt(j)  = send_cnt(r)
         hi%rcnt(j)  = want_cnt(r)
      end do

      hi%soff(1) = 0; hi%roff(1) = 0
      do j = 1, hi%npeers
         hi%soff(j+1) = hi%soff(j) + hi%scnt(j)
         hi%roff(j+1) = hi%roff(j) + hi%rcnt(j)
      end do
      allocate( hi%sids(hi%soff(hi%npeers+1)), &
                hi%rids(hi%roff(hi%npeers+1)) )

      ! send ids: map incoming requested global IDs to my owned local IDs
      do j = 1, hi%npeers
         r = hi%peers(j)
         do p = 1, hi%scnt(j)
            g = in_global( rdispl(r) + p )
            hi%sids(hi%soff(j)+p) = lm%cell_g2l(g)
         end do
      end do
      ! recv ids: local halo IDs in want order
      do j = 1, hi%npeers
         r = hi%peers(j)
         do p = 1, hi%rcnt(j)
            g = want_global( sdispl(r) + p )
            hi%rids(hi%roff(j)+p) = lm%cell_g2l(g)
         end do
      end do

      deallocate( want_cnt, send_cnt, want_global, want_cursor, in_global, &
                  sdispl, rdispl )

   end subroutine halo_setup

   !----------------------------------------------------------------------------
   ! Exchange one scalar field: update field() at halo ids.
   !----------------------------------------------------------------------------
   subroutine halo_exchange_scalar( hi, field, ierr )
      type(halo_info_t), intent(in)    :: hi
      real(dp),          intent(inout) :: field(:)
      integer,           intent(out)   :: ierr

      integer :: ntot_s, ntot_r
      real(dp), allocatable :: sbuf(:), rbuf(:)
      integer :: k

      ntot_s = hi%soff(hi%npeers+1)
      ntot_r = hi%roff(hi%npeers+1)
      allocate( sbuf(ntot_s), rbuf(ntot_r) )

      do k = 1, ntot_s
         sbuf(k) = field(hi%sids(k))
      end do

      call halo_comm( hi, 1, sbuf, rbuf, ierr )
      if ( ierr /= 0 ) return

      do k = 1, ntot_r
         field(hi%rids(k)) = rbuf(k)
      end do
   end subroutine halo_exchange_scalar

   !----------------------------------------------------------------------------
   ! Exchange a vector field with shape (ncomp, nlocal): update halo columns.
   !----------------------------------------------------------------------------
   subroutine halo_exchange_vector( hi, field, ierr )
      type(halo_info_t), intent(in)    :: hi
      real(dp),          intent(inout) :: field(:,:)
      integer,           intent(out)   :: ierr

      integer :: ncomp, ntot_s, ntot_r
      real(dp), allocatable :: sbuf(:), rbuf(:)
      integer :: k, c, m

      ncomp  = size(field,1)
      ntot_s = hi%soff(hi%npeers+1)
      ntot_r = hi%roff(hi%npeers+1)
      allocate( sbuf(ncomp*ntot_s), rbuf(ncomp*ntot_r) )

      ! pack: components of each cell contiguous
      do k = 1, ntot_s
         c = hi%sids(k)
         do m = 1, ncomp
            sbuf(ncomp*(k-1)+m) = field(m,c)
         end do
      end do

      call halo_comm( hi, ncomp, sbuf, rbuf, ierr )
      if ( ierr /= 0 ) return

      do k = 1, ntot_r
         c = hi%rids(k)
         do m = 1, ncomp
            field(m,c) = rbuf(ncomp*(k-1)+m)
         end do
      end do
   end subroutine halo_exchange_vector

   !----------------------------------------------------------------------------
   ! Core message exchange over already-packed contiguous buffers.
   ! sbuf / rbuf layout: per-peer blocks in peers() order, each block holds
   ! ncomp*ncnt values.
   !----------------------------------------------------------------------------
   subroutine halo_comm( hi, ncomp, sbuf, rbuf, ierr )
      type(halo_info_t), intent(in)  :: hi
      integer,           intent(in)  :: ncomp
      real(dp),          intent(in)  :: sbuf(:)
      real(dp),          intent(out) :: rbuf(:)
      integer,           intent(out) :: ierr

      integer :: j, nmsg, peer, cnt
      integer, allocatable :: req(:), stat(:,:)

      nmsg = 2 * hi%npeers
      if ( nmsg == 0 ) then
         ierr = 0
         return
      end if
      allocate( req(nmsg), stat(MPI_STATUS_SIZE,nmsg) )

      ! post all receives first
      do j = 1, hi%npeers
         peer = hi%peers(j)
         cnt  = ncomp * hi%rcnt(j)
         if ( cnt > 0 ) then
            call MPI_Irecv( rbuf(hi%roff(j)*ncomp+1), cnt, &
                            MPI_DOUBLE_PRECISION, peer, HALO_TAG, &
                            MPI_COMM_WORLD, req(j), ierr )
            call mpi_check( ierr, 'halo_comm irecv' )
         else
            req(j) = MPI_REQUEST_NULL
         end if
      end do
      ! post all sends
      do j = 1, hi%npeers
         peer = hi%peers(j)
         cnt  = ncomp * hi%scnt(j)
         if ( cnt > 0 ) then
            call MPI_Isend( sbuf(hi%soff(j)*ncomp+1), cnt, &
                            MPI_DOUBLE_PRECISION, peer, HALO_TAG, &
                            MPI_COMM_WORLD, req(hi%npeers+j), ierr )
            call mpi_check( ierr, 'halo_comm isend' )
         else
            req(hi%npeers+j) = MPI_REQUEST_NULL
         end if
      end do

      call MPI_Waitall( nmsg, req, stat, ierr )
      call mpi_check( ierr, 'halo_comm waitall' )
      deallocate( req, stat )
   end subroutine halo_comm

end module mod_uns_halo
