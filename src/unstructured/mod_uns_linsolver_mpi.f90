!===============================================================================
! mod_linsolver_mpi.f90 -- Distributed sparse linear solvers (Step 4)
!
! Local matrix convention:
!   A is a csr_t with nrows >= nowned. Rows 1..nowned are the OWNED rows and
!   are fully populated; their columns may reference
!     - owned cells : 1..nowned
!     - halo  cells : nowned+1..nowned+nghost
!   Halo rows (nowned+1..) are never read and may be empty.
!
! Strategy:
!   - Krylov vectors are OWNED-sized (exactly the rows this rank solves). For a
!     matrix product, the owned input is copied into a full-length work array,
!     halo-exchanged, and the owned output rows are evaluated with the full
!     stencil. Thus no owned-sized array is ever indexed with a halo id.
!   - preconditioning: extract the owned-owned block B (halo columns dropped)
!     and reuse the serial ilu0/ic0 factor + apply on B. The residual already
!     contains the halo contributions, so dropping them from the
!     preconditioner is the standard distributed block treatment.
!   - all inner products / norms are global (MPI_Allreduce over owned entries)
!     so every rank follows the identical iteration path.
!===============================================================================
module mod_uns_linsolver_mpi
   use mpi
   use mod_precision, only: dp, ip, pi
   use mod_uns_mpi_core, only: mpi_check
   use mod_uns_halo
   use mod_uns_linsolver, only: csr_t, ilu0_factor, ilu0_apply, ic0_factor, &
                            ic0_apply
   implicit none
   private
   public :: matvec_owned, owned_block, &
             bicgstab_ilu0_mpi, cg_jacobi_mpi, cg_ic0_mpi

contains

   !----------------------------------------------------------------------------
   ! Global dot product over owned entries.
   !----------------------------------------------------------------------------
   real(dp) function gdot( a, b )
      real(dp), intent(in) :: a(:), b(:)
      integer :: ierr
      real(dp) :: s
      s = dot_product( a, b )
      call MPI_Allreduce( MPI_IN_PLACE, s, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                          MPI_COMM_WORLD, ierr )
      call mpi_check( ierr, 'gdot' )
      gdot = s
   end function gdot

   !----------------------------------------------------------------------------
   ! Distributed matrix-vector product on FULL vectors (used by callers that
   ! already hold owned + halo storage): x is halo-exchanged, then the owned
   ! output rows are evaluated.
   !----------------------------------------------------------------------------
   subroutine matvec_owned( A, x, y, hi, nowned )
      type(csr_t),       intent(in)    :: A
      real(dp),          intent(inout) :: x(:)      ! halo entries refreshed
      real(dp),          intent(out)   :: y(:)
      type(halo_info_t), intent(in)    :: hi
      integer,           intent(in)    :: nowned
      integer :: i, k, ierr

      call halo_exchange_scalar( hi, x, ierr )
      call mpi_check( ierr, 'matvec_owned exchange' )

      do i = 1, nowned
         y(i) = 0.0_dp
         do k = A%row_ptr(i), A%row_ptr(i+1)-1
            y(i) = y(i) + A%val(k) * x(A%col_idx(k))
         end do
      end do
   end subroutine matvec_owned

   !----------------------------------------------------------------------------
   ! Matrix-vector product for OWNED-sized Krylov vectors. The owned input xo
   ! is copied into the full-length work array xfull and halo-exchanged before
   ! the owned output rows yo are evaluated.
   !----------------------------------------------------------------------------
   subroutine matvec_krylov( A, xo, yo, hi, nowned, xfull )
      type(csr_t),       intent(in)  :: A
      real(dp),          intent(in)  :: xo(:)    ! size nowned
      real(dp),          intent(out) :: yo(:)    ! size nowned
      type(halo_info_t), intent(in)  :: hi
      integer,           intent(in)  :: nowned
      real(dp),          intent(inout) :: xfull(:)   ! full-length work
      integer :: i, k, ierr

      xfull(1:nowned) = xo(1:nowned)
      call halo_exchange_scalar( hi, xfull, ierr )
      call mpi_check( ierr, 'matvec_krylov exchange' )

      do i = 1, nowned
         yo(i) = 0.0_dp
         do k = A%row_ptr(i), A%row_ptr(i+1)-1
            yo(i) = yo(i) + A%val(k) * xfull(A%col_idx(k))
         end do
      end do
   end subroutine matvec_krylov

   !----------------------------------------------------------------------------
   ! Extract the owned-owned block: nowned rows, only columns <= nowned.
   ! Column numbering is unchanged (stays 1..nowned); columns stay sorted.
   !----------------------------------------------------------------------------
   subroutine owned_block( A, nowned, B )
      type(csr_t), intent(in)  :: A
      integer,     intent(in)  :: nowned
      type(csr_t), intent(out) :: B

      integer :: i, k, nkeep
      integer, allocatable :: keep(:)
      integer :: maxdeg

      B%nrows = nowned
      allocate( B%row_ptr(nowned+1) )

      ! Pass 1: per-row kept counts (builds row_ptr) and max kept degree.
      maxdeg = 0
      B%row_ptr(1) = 1
      do i = 1, nowned
         nkeep = 0
         do k = A%row_ptr(i), A%row_ptr(i+1)-1
            if ( A%col_idx(k) <= nowned ) nkeep = nkeep + 1
         end do
         maxdeg = max( maxdeg, nkeep )
         B%row_ptr(i+1) = B%row_ptr(i) + nkeep
      end do
      allocate( keep(maxdeg) )

      allocate( B%col_idx(B%row_ptr(nowned+1)-1), &
                B%val(B%row_ptr(nowned+1)-1) )

      ! Pass 2: per row, record the kept source positions and copy immediately
      ! (keep must be refilled for every row -- it is not valid across rows).
      do i = 1, nowned
         nkeep = 0
         do k = A%row_ptr(i), A%row_ptr(i+1)-1
            if ( A%col_idx(k) <= nowned ) then
               nkeep = nkeep + 1
               keep(nkeep) = k
            end if
         end do
         do k = 1, nkeep
            B%col_idx(B%row_ptr(i)+k-1) = A%col_idx(keep(k))
            B%val(B%row_ptr(i)+k-1)     = A%val(keep(k))
         end do
      end do
      deallocate( keep )
   end subroutine owned_block

   !----------------------------------------------------------------------------
   ! Distributed BiCGSTAB with ILU(0) on the owned block.
   ! ierr: 0 converged, 1 max iterations reached, 2 breakdown
   !----------------------------------------------------------------------------
   subroutine bicgstab_ilu0_mpi( A, b, x, nowned, hi, tol, maxit, it, res, ierr )
      type(csr_t),       intent(in)    :: A
      real(dp),          intent(in)    :: b(:)    ! size nowned
      real(dp),          intent(inout) :: x(:)    ! size nowned
      integer,           intent(in)    :: nowned
      type(halo_info_t), intent(in)    :: hi
      real(dp),          intent(in)    :: tol
      integer,           intent(in)    :: maxit
      integer,           intent(out)   :: it
      real(dp),          intent(out)   :: res
      integer,           intent(out)   :: ierr

      type(csr_t) :: Bblk, lu
      real(dp), allocatable :: r0(:), r(:), p(:), v(:), s(:), t(:), &
                               ph(:), sh(:), xfull(:)
      real(dp) :: bnrm, rho, rho_old, alpha, omega, beta, tt, snrm

      it = 0; ierr = 0; res = 0.0_dp
      bnrm = sqrt( gdot(b, b) )
      if ( bnrm == 0.0_dp ) then
         x = 0.0_dp
         return
      end if

      ! ILU(0) factor the owned-owned block
      call owned_block( A, nowned, Bblk )
      lu = Bblk
      call ilu0_factor( lu, ierr )
      if ( ierr /= 0 ) then; ierr = ierr + 10; return; end if

      allocate( r0(nowned), r(nowned), p(nowned), v(nowned), s(nowned), &
                t(nowned), ph(nowned), sh(nowned), xfull(A%nrows) )

      call matvec_krylov( A, x, v, hi, nowned, xfull )
      r0 = b - v
      r  = r0
      p = 0.0_dp; v = 0.0_dp
      rho_old = 1.0_dp; alpha = 1.0_dp; omega = 1.0_dp
      res = sqrt( gdot(r0,r0) ) / bnrm
      if ( res <= tol ) return

      do it = 1, maxit
         rho = gdot( r, r0 )
         if ( abs(rho) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         beta = ( rho / rho_old ) * ( alpha / omega )
         p = r + beta * ( p - omega * v )
         call ilu0_apply( lu, p, ph )
         call matvec_krylov( A, ph, v, hi, nowned, xfull )
         alpha = rho / gdot( r0, v )
         s = r - alpha * v
         snrm = sqrt( gdot(s,s) )
         if ( snrm / bnrm <= tol ) then
            x = x + alpha * ph
            res = snrm / bnrm
            return
         end if
         call ilu0_apply( lu, s, sh )
         call matvec_krylov( A, sh, t, hi, nowned, xfull )
         tt = gdot( t, t )
         if ( tt < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         omega = gdot( t, s ) / tt
         if ( abs(omega) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         x = x + alpha*ph + omega*sh
         r = s - omega * t
         res = sqrt( gdot(r,r) ) / bnrm
         if ( res <= tol ) return
         rho_old = rho
      end do

      if ( ierr == 0 ) then
         ierr = 1
         it = maxit
      end if
      deallocate( r0, r, p, v, s, t, ph, sh, xfull )
   end subroutine bicgstab_ilu0_mpi

   !----------------------------------------------------------------------------
   ! Distributed CG with Jacobi (diagonal) preconditioning.
   ! ierr: 0 converged, 1 max iterations reached, 2 breakdown
   !----------------------------------------------------------------------------
   subroutine cg_jacobi_mpi( A, b, x, nowned, hi, tol, maxit, it, res, ierr )
      type(csr_t),       intent(in)    :: A
      real(dp),          intent(in)    :: b(:)
      real(dp),          intent(inout) :: x(:)
      integer,           intent(in)    :: nowned
      type(halo_info_t), intent(in)    :: hi
      real(dp),          intent(in)    :: tol
      integer,           intent(in)    :: maxit
      integer,           intent(out)   :: it
      real(dp),          intent(out)   :: res
      integer,           intent(out)   :: ierr

      real(dp), allocatable :: r(:), z(:), p(:), apv(:), diag(:), xfull(:)
      real(dp) :: bnrm, rz, rz_new, pap, alpha, beta
      integer  :: i, k

      it = 0; ierr = 0; res = 0.0_dp
      bnrm = sqrt( gdot(b,b) )
      if ( bnrm == 0.0_dp ) then
         x = 0.0_dp
         return
      end if

      allocate( r(nowned), z(nowned), p(nowned), apv(nowned), diag(nowned), &
                xfull(A%nrows) )
      do i = 1, nowned
         diag(i) = 0.0_dp
         do k = A%row_ptr(i), A%row_ptr(i+1)-1
            if ( A%col_idx(k) == i ) then
               diag(i) = A%val(k); exit
            end if
         end do
         if ( abs(diag(i)) < tiny(1.0_dp) ) then
            write(*,'(a,i0)') 'CG_MPI ERROR: missing/bad diagonal row ', i
            ierr = 2; return
         end if
      end do

      call matvec_krylov( A, x, apv, hi, nowned, xfull )
      r = b - apv
      res = sqrt( gdot(r,r) ) / bnrm
      if ( res <= tol ) return
      z = r / diag
      p = z
      rz = gdot( r, z )

      do it = 1, maxit
         call matvec_krylov( A, p, apv, hi, nowned, xfull )
         pap = gdot( p, apv )
         if ( abs(pap) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         alpha = rz / pap
         x = x + alpha*p
         r = r - alpha*apv
         res = sqrt( gdot(r,r) ) / bnrm
         if ( res <= tol ) return
         z = r / diag
         rz_new = gdot( r, z )
         if ( abs(rz) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         beta = rz_new / rz
         p = z + beta*p
         rz = rz_new
      end do

      if ( ierr == 0 ) ierr = 1
      deallocate( r, z, p, apv, diag, xfull )
   end subroutine cg_jacobi_mpi

   !----------------------------------------------------------------------------
   ! Distributed CG with ICC(0) on the owned block.
   ! ierr: 0 converged, 1 max iterations reached, 2 breakdown
   !----------------------------------------------------------------------------
   subroutine cg_ic0_mpi( A, b, x, nowned, hi, tol, maxit, it, res, ierr )
      type(csr_t),       intent(in)    :: A
      real(dp),          intent(in)    :: b(:)
      real(dp),          intent(inout) :: x(:)
      integer,           intent(in)    :: nowned
      type(halo_info_t), intent(in)    :: hi
      real(dp),          intent(in)    :: tol
      integer,           intent(in)    :: maxit
      integer,           intent(out)   :: it
      real(dp),          intent(out)   :: res
      integer,           intent(out)   :: ierr

      type(csr_t) :: Bblk, L
      real(dp), allocatable :: r(:), z(:), p(:), apv(:), xfull(:)
      real(dp) :: bnrm, rz, rz_new, pap, alpha, beta

      it = 0; ierr = 0; res = 0.0_dp
      bnrm = sqrt( gdot(b,b) )
      if ( bnrm == 0.0_dp ) then
         x = 0.0_dp
         return
      end if

      call owned_block( A, nowned, Bblk )
      L = Bblk
      call ic0_factor( L, ierr )
      if ( ierr /= 0 ) then; ierr = ierr + 10; return; end if

      allocate( r(nowned), z(nowned), p(nowned), apv(nowned), xfull(A%nrows) )

      call matvec_krylov( A, x, apv, hi, nowned, xfull )
      r = b - apv
      res = sqrt( gdot(r,r) ) / bnrm
      if ( res <= tol ) return
      call ic0_apply( L, r, z )
      p = z
      rz = gdot( r, z )

      do it = 1, maxit
         call matvec_krylov( A, p, apv, hi, nowned, xfull )
         pap = gdot( p, apv )
         if ( abs(pap) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         alpha = rz / pap
         x = x + alpha*p
         r = r - alpha*apv
         res = sqrt( gdot(r,r) ) / bnrm
         if ( res <= tol ) return
         call ic0_apply( L, r, z )
         rz_new = gdot( r, z )
         if ( abs(rz) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         beta = rz_new / rz
         p = z + beta*p
         rz = rz_new
      end do

      if ( ierr == 0 ) ierr = 1
      deallocate( r, z, p, apv, xfull )
   end subroutine cg_ic0_mpi

end module mod_uns_linsolver_mpi
