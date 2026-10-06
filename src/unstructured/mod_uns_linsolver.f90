!===============================================================================
! mod_linsolver.f90 -- Built-in sparse linear solvers (phase 2)
!
! Storage : CSR (compressed sparse row). Column indices within a row must
!           be sorted in increasing order (required by ILU(0)).
! Solvers :
!   bicgstab_ilu0 -- right-preconditioned BiCGSTAB with ILU(0); general
!                    nonsymmetric matrices (workhorse for SIMPLE/PISO).
!   cg_jacobi     -- preconditioned Conjugate Gradient with Jacobi (diagonal)
!                    preconditioning; for symmetric positive definite systems.
!   cg_ic0        -- preconditioned CG with incomplete Cholesky (zero fill,
!                    ICC(0)); much stronger than Jacobi on Laplacian-type
!                    pressure-correction matrices.
! Notes  :
!   - bicgstab_ilu0 / cg_ic0 factor a private copy of A on each call. The
!     pressure-correction matrix changes every outer iteration, so caching
!     the factorization is not beneficial yet.
!   - Convergence measure: ||r||_2 / max(||b||_2, tiny) <= tol.
!===============================================================================
module mod_uns_linsolver
   use mod_precision, only: dp, ip, pi
   implicit none
   private
   public :: csr_t, csr_matvec, ilu0_factor, ilu0_apply, bicgstab_ilu0, &
             cg_jacobi, ic0_factor, ic0_apply, cg_ic0

   type :: csr_t
      integer              :: nrows = 0
      integer, allocatable :: row_ptr(:)    ! size nrows+1
      integer, allocatable :: col_idx(:)    ! size nnz, sorted per row
      real(dp), allocatable :: val(:)       ! size nnz
      integer, allocatable :: diag_pos(:)   ! position of A(i,i) in val
   end type csr_t

contains

   !----------------------------------------------------------------------------
   ! y = A * x
   !----------------------------------------------------------------------------
   subroutine csr_matvec( A, x, y )
      type(csr_t), intent(in)  :: A
      real(dp),    intent(in)  :: x(:)
      real(dp),    intent(out) :: y(:)
      integer :: i, k

      do i = 1, A%nrows
         y(i) = 0.0_dp
         do k = A%row_ptr(i), A%row_ptr(i+1)-1
            y(i) = y(i) + A%val(k) * x(A%col_idx(k))
         end do
      end do
   end subroutine csr_matvec

   !----------------------------------------------------------------------------
   ! Find the position of column j in row i; 0 if absent.
   !----------------------------------------------------------------------------
   integer function csr_find( A, i, j )
      type(csr_t), intent(in) :: A
      integer,     intent(in) :: i, j
      integer :: k

      csr_find = 0
      do k = A%row_ptr(i), A%row_ptr(i+1)-1
         if ( A%col_idx(k) == j ) then
            csr_find = k
            return
         else if ( A%col_idx(k) > j ) then
            return                            ! sorted: cannot appear later
         end if
      end do
   end function csr_find

   !----------------------------------------------------------------------------
   ! ILU(0) factorization in the IKJ variant; done in place.
   ! L is unit lower triangular, U upper; both stored in the original
   ! sparsity pattern (no fill-in). Requires a full diagonal and sorted
   ! column indices.
   !----------------------------------------------------------------------------
   subroutine ilu0_factor( A, ierr )
      type(csr_t), intent(inout) :: A
      integer,     intent(out)   :: ierr

      integer  :: i, k, kk, jj, pos, j
      real(dp) :: aik

      ierr = 0
      if ( .not. allocated(A%diag_pos) ) allocate( A%diag_pos(A%nrows) )
      do i = 1, A%nrows
         A%diag_pos(i) = csr_find( A, i, i )
         if ( A%diag_pos(i) == 0 ) then
            write(*,'(a,i0)') 'ILU0 ERROR: missing diagonal in row ', i
            ierr = 1
            return
         end if
      end do

      do i = 2, A%nrows
         do kk = A%row_ptr(i), A%diag_pos(i)-1      ! strictly lower part
            k = A%col_idx(kk)
            if ( abs( A%val(A%diag_pos(k)) ) < tiny(1.0_dp) ) then
               write(*,'(a,i0)') 'ILU0 ERROR: zero pivot in row ', k
               ierr = 2
               return
            end if
            aik = A%val(kk) / A%val(A%diag_pos(k))
            A%val(kk) = aik
            do jj = A%diag_pos(k), A%row_ptr(k+1)-1  ! row k, cols >= k
               j = A%col_idx(jj)
               if ( j <= k ) cycle
               pos = csr_find( A, i, j )
               if ( pos > 0 ) A%val(pos) = A%val(pos) - aik * A%val(jj)
            end do
         end do
      end do
   end subroutine ilu0_factor

   !----------------------------------------------------------------------------
   ! Apply z = (LU)^{-1} r with unit-diagonal L (forward) then U (backward)
   !----------------------------------------------------------------------------
   subroutine ilu0_apply( A, r, z )
      type(csr_t), intent(in)  :: A
      real(dp),    intent(in)  :: r(:)
      real(dp),    intent(out) :: z(:)
      integer :: i, k
      real(dp) :: s

      z = r
      do i = 1, A%nrows                             ! forward: L z = r
         s = z(i)
         do k = A%row_ptr(i), A%diag_pos(i)-1
            s = s - A%val(k) * z(A%col_idx(k))
         end do
         z(i) = s
      end do
      do i = A%nrows, 1, -1                         ! backward: U z = z
         s = z(i)
         do k = A%diag_pos(i)+1, A%row_ptr(i+1)-1
            s = s - A%val(k) * z(A%col_idx(k))
         end do
         z(i) = s / A%val(A%diag_pos(i))
      end do
   end subroutine ilu0_apply

   !----------------------------------------------------------------------------
   ! Preconditioned BiCGSTAB with ILU(0)
   ! ierr: 0 converged, 1 max iterations reached, 2 breakdown
   !----------------------------------------------------------------------------
   subroutine bicgstab_ilu0( A, b, x, tol, maxit, it, res, ierr, verb )
      type(csr_t), intent(in)    :: A
      real(dp),    intent(in)    :: b(:)
      real(dp),    intent(inout) :: x(:)          ! initial guess on entry
      real(dp),    intent(in)    :: tol
      integer,     intent(in)    :: maxit
      integer,     intent(out)   :: it            ! iterations used
      real(dp),    intent(out)   :: res           ! final relative residual
      integer,     intent(out)   :: ierr
      integer,     intent(in), optional :: verb   ! print residual every verb

      type(csr_t) :: lu
      real(dp), allocatable :: r0(:), r(:), p(:), v(:), s(:), t(:), &
                               ph(:), sh(:), tmp(:)
      real(dp) :: bnrm, rho, rho_old, alpha, omega, beta, tt, snrm

      it = 0; ierr = 0; res = 0.0_dp
      bnrm = norm2( b )
      if ( bnrm == 0.0_dp ) then                   ! homogeneous system
         x = 0.0_dp
         return
      end if

      ! private ILU(0) factorization of A
      lu = A
      call ilu0_factor( lu, ierr )
      if ( ierr /= 0 ) then; ierr = ierr + 10; return; end if

      allocate( r0(A%nrows), r(A%nrows), p(A%nrows), v(A%nrows), &
                s(A%nrows), t(A%nrows), ph(A%nrows), sh(A%nrows), &
                tmp(A%nrows) )

      call csr_matvec( A, x, tmp )
      r0 = b - tmp
      r = r0
      p = 0.0_dp; v = 0.0_dp
      rho_old = 1.0_dp; alpha = 1.0_dp; omega = 1.0_dp
      res = norm2( r0 ) / bnrm
      if ( res <= tol ) return

      do it = 1, maxit
         rho = dot_product( r, r0 )
         if ( abs(rho) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         beta = ( rho / rho_old ) * ( alpha / omega )
         p = r + beta * ( p - omega * v )
         call ilu0_apply( lu, p, ph )
         call csr_matvec( A, ph, v )
         alpha = rho / dot_product( r0, v )   ! NB: inner product with r0!
         s = r - alpha * v
         snrm = norm2( s )
         if ( snrm / bnrm <= tol ) then
            x = x + alpha * ph
            res = snrm / bnrm
            return
         end if
         call ilu0_apply( lu, s, sh )
         call csr_matvec( A, sh, t )
         tt = dot_product( t, t )
         if ( tt < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         omega = dot_product( t, s ) / tt
         if ( abs(omega) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         x = x + alpha * ph + omega * sh
         r = s - omega * t
         res = norm2( r ) / bnrm
         if ( present(verb) ) then
            if ( verb > 0 .and. mod(it,verb) == 0 ) &
               write(*,'(a,i0,a,es10.3)') '  bicg it=', it, '  res=', res
         end if
         if ( res <= tol ) return
         rho_old = rho
      end do

      if ( ierr == 0 ) then
         ierr = 1                                ! no convergence
         it = maxit
      end if
   end subroutine bicgstab_ilu0

   !----------------------------------------------------------------------------
   ! Conjugate Gradient with Jacobi (diagonal) preconditioning, for SPD A
   ! ierr: 0 converged, 1 max iterations reached, 2 breakdown
   !----------------------------------------------------------------------------
   subroutine cg_jacobi( A, b, x, tol, maxit, it, res, ierr )
      type(csr_t), intent(in)    :: A
      real(dp),    intent(in)    :: b(:)
      real(dp),    intent(inout) :: x(:)
      real(dp),    intent(in)    :: tol
      integer,     intent(in)    :: maxit
      integer,     intent(out)   :: it
      real(dp),    intent(out)   :: res
      integer,     intent(out)   :: ierr

      real(dp), allocatable :: r(:), z(:), p(:), ap(:), diag(:)
      real(dp) :: bnrm, rz, rz_new, pap, alpha, beta
      integer  :: i, k

      it = 0; ierr = 0; res = 0.0_dp
      bnrm = norm2( b )
      if ( bnrm == 0.0_dp ) then
         x = 0.0_dp
         return
      end if

      allocate( r(A%nrows), z(A%nrows), p(A%nrows), ap(A%nrows), &
                diag(A%nrows) )
      do i = 1, A%nrows                          ! extract the diagonal
         k = csr_find( A, i, i )
         if ( k == 0 ) then
            write(*,'(a,i0)') 'CG ERROR: missing diagonal in row ', i
            ierr = 2
            return
         end if
         diag(i) = A%val(k)
      end do

      call csr_matvec( A, x, ap )
      r = b - ap
      res = norm2( r ) / bnrm
      if ( res <= tol ) return
      z = r / diag
      p = z
      rz = dot_product( r, z )

      do it = 1, maxit
         call csr_matvec( A, p, ap )
         pap = dot_product( p, ap )
         if ( abs(pap) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         alpha = rz / pap
         x = x + alpha * p
         r = r - alpha * ap
         res = norm2( r ) / bnrm
         if ( res <= tol ) return
         z = r / diag
         rz_new = dot_product( r, z )
         if ( abs(rz) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         beta = rz_new / rz
         p = z + beta * p
         rz = rz_new
      end do

      if ( ierr == 0 ) ierr = 1

      deallocate( r, z, p, ap, diag )
   end subroutine cg_jacobi

   !----------------------------------------------------------------------------
   ! Dot product of values at matching column indices in two CSR slices
   ! (intersection by two-pointer merge; both slices must be column-sorted)
   !----------------------------------------------------------------------------
   real(dp) function ic0_dot( A, p1, p2, q1, q2 )
      type(csr_t), intent(in) :: A
      integer,     intent(in) :: p1, p2, q1, q2
      integer :: p, q

      p = p1; q = q1
      ic0_dot = 0.0_dp
      do while ( p <= p2 .and. q <= q2 )
         if ( A%col_idx(p) < A%col_idx(q) ) then
            p = p + 1
         else if ( A%col_idx(p) > A%col_idx(q) ) then
            q = q + 1
         else
            ic0_dot = ic0_dot + A%val(p) * A%val(q)
            p = p + 1; q = q + 1
         end if
      end do
   end function ic0_dot

   !----------------------------------------------------------------------------
   ! Incomplete Cholesky factorization without fill-in (ICC(0)), done in
   ! place: overwrites the lower triangle and diagonal with L such that the
   ! stored A ~ L L^T; the upper triangle is left untouched and unused.
   !
   ! For each row i:
   !   L(i,k) = ( A(i,k) - sum_{j in pat(i) intersect pat(k), j<k}
   !                       L(i,j) L(k,j) ) / L(k,k)        for k < i
   !   L(i,i) = sqrt( A(i,i) - sum_{k in pat(i), k<i} L(i,k)^2 )
   ! Exists without breakdown for symmetric M-matrices (e.g. the pressure-
   ! correction matrix). Requires a full diagonal and sorted columns.
   !----------------------------------------------------------------------------
   subroutine ic0_factor( A, ierr )
      type(csr_t), intent(inout) :: A
      integer,     intent(out)   :: ierr

      integer  :: i, kp, k, dp_i
      real(dp) :: s

      ierr = 0
      if ( .not. allocated(A%diag_pos) ) allocate( A%diag_pos(A%nrows) )
      do i = 1, A%nrows
         A%diag_pos(i) = csr_find( A, i, i )
         if ( A%diag_pos(i) == 0 ) then
            write(*,'(a,i0)') 'ICC0 ERROR: missing diagonal in row ', i
            ierr = 1
            return
         end if
      end do

      do i = 1, A%nrows
         dp_i = A%diag_pos(i)
         do kp = A%row_ptr(i), dp_i - 1
            k = A%col_idx(kp)
            s = A%val(kp) - ic0_dot( A, A%row_ptr(i), kp-1, &
                                    A%row_ptr(k), A%diag_pos(k)-1 )
            A%val(kp) = s / A%val(A%diag_pos(k))
         end do
         s = A%val(dp_i)
         do kp = A%row_ptr(i), dp_i - 1
            s = s - A%val(kp)**2
         end do
         if ( s <= 0.0_dp ) then
            write(*,'(a,i0)') 'ICC0 ERROR: non-positive pivot in row ', i
            ierr = 2
            return
         end if
         A%val(dp_i) = sqrt( s )
      end do
   end subroutine ic0_factor

   !----------------------------------------------------------------------------
   ! Apply z = (L L^T)^{-1} r: forward solve L y = r, then backward solve
   ! L^T z = y. Only the lower triangle and diagonal of L are used.
   !----------------------------------------------------------------------------
   subroutine ic0_apply( L, r, z )
      type(csr_t), intent(in)  :: L
      real(dp),    intent(in)  :: r(:)
      real(dp),    intent(out) :: z(:)

      integer  :: i, kp, k
      real(dp) :: s
      real(dp), allocatable :: w(:)

      allocate( w(L%nrows) )
      ! forward: L y = r (y stored in z)
      do i = 1, L%nrows
         s = r(i)
         do kp = L%row_ptr(i), L%diag_pos(i)-1
            s = s - L%val(kp) * z(L%col_idx(kp))
         end do
         z(i) = s / L%val(L%diag_pos(i))
      end do
      ! backward: L^T z = y (w evolves as the right-hand side)
      w = z
      do i = L%nrows, 1, -1
         z(i) = w(i) / L%val(L%diag_pos(i))
         do kp = L%row_ptr(i), L%diag_pos(i)-1
            k = L%col_idx(kp)
            w(k) = w(k) - L%val(kp) * z(i)
         end do
      end do
      deallocate( w )
   end subroutine ic0_apply

   !----------------------------------------------------------------------------
   ! Preconditioned Conjugate Gradient with ICC(0), for SPD A
   ! ierr: 0 converged, 1 max iterations reached, 2 breakdown
   !----------------------------------------------------------------------------
   subroutine cg_ic0( A, b, x, tol, maxit, it, res, ierr )
      type(csr_t), intent(in)    :: A
      real(dp),    intent(in)    :: b(:)
      real(dp),    intent(inout) :: x(:)
      real(dp),    intent(in)    :: tol
      integer,     intent(in)    :: maxit
      integer,     intent(out)   :: it
      real(dp),    intent(out)   :: res
      integer,     intent(out)   :: ierr

      type(csr_t) :: L
      real(dp), allocatable :: r(:), z(:), p(:), ap(:)
      real(dp) :: bnrm, rz, rz_new, pap, alpha, beta

      it = 0; ierr = 0; res = 0.0_dp
      bnrm = norm2( b )
      if ( bnrm == 0.0_dp ) then
         x = 0.0_dp
         return
      end if

      ! private ICC(0) factorization of A
      L = A
      call ic0_factor( L, ierr )
      if ( ierr /= 0 ) then; ierr = ierr + 10; return; end if

      allocate( r(A%nrows), z(A%nrows), p(A%nrows), ap(A%nrows) )

      call csr_matvec( A, x, ap )
      r = b - ap
      res = norm2( r ) / bnrm
      if ( res <= tol ) return
      call ic0_apply( L, r, z )
      p = z
      rz = dot_product( r, z )

      do it = 1, maxit
         call csr_matvec( A, p, ap )
         pap = dot_product( p, ap )
         if ( abs(pap) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         alpha = rz / pap
         x = x + alpha * p
         r = r - alpha * ap
         res = norm2( r ) / bnrm
         if ( res <= tol ) return
         call ic0_apply( L, r, z )
         rz_new = dot_product( r, z )
         if ( abs(rz) < tiny(1.0_dp) ) then; ierr = 2; exit; end if
         beta = rz_new / rz
         p = z + beta * p
         rz = rz_new
      end do

      if ( ierr == 0 ) ierr = 1

      deallocate( r, z, p, ap )
   end subroutine cg_ic0

end module mod_uns_linsolver
