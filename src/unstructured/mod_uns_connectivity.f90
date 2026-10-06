!===============================================================================
! mod_connectivity.f90 -- Build cell-based connectivity from face lists
!
! Produces:
!   cf_ptr / cf     : CSR-style cell -> face lists
!   c2c_ptr / c2c   : CSR-style cell -> neighbor cell lists (boundary faces
!                     are skipped)
!   nbf / bfaces    : number of boundary faces and their index list
!===============================================================================
module mod_uns_connectivity
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   implicit none
   private
   public :: conn_t, build_connectivity

   type :: conn_t
      integer, allocatable :: cf_ptr(:)     ! size ncells+1
      integer, allocatable :: cf(:)         ! face indices per cell
      integer, allocatable :: c2c_ptr(:)    ! size ncells+1
      integer, allocatable :: c2c(:)        ! neighbor cells per cell
      integer              :: nbf = 0
      integer, allocatable :: bfaces(:)     ! boundary face indices
   end type conn_t

contains

   subroutine build_connectivity( m, c, ier )
      type(mesh_t),  intent(in)  :: m
      type(conn_t),  intent(out) :: c
      integer,       intent(out) :: ier

      integer :: i, k, n, cell

      ier = 0

      ! ---- pass 1: count faces per cell ------------------------------------
      allocate( c%cf_ptr(m%ncells+1) )
      c%cf_ptr = 0
      do i = 1, m%nfaces
         c%cf_ptr(m%f(i)%c0+1) = c%cf_ptr(m%f(i)%c0+1) + 1
         if ( m%f(i)%c1 > 0 ) &
            c%cf_ptr(m%f(i)%c1+1) = c%cf_ptr(m%f(i)%c1+1) + 1
      end do
      do k = 1, m%ncells
         c%cf_ptr(k+1) = c%cf_ptr(k+1) + c%cf_ptr(k)
      end do

      allocate( c%cf(c%cf_ptr(m%ncells+1)) )

      ! ---- pass 2: fill ------------------------------------------------------
      faces_fill: block
         integer, allocatable :: fill(:)
         allocate( fill(m%ncells) )
         fill = c%cf_ptr(1:m%ncells) + 1
         do i = 1, m%nfaces
            cell = m%f(i)%c0
            c%cf(fill(cell)) = i; fill(cell) = fill(cell) + 1
            if ( m%f(i)%c1 > 0 ) then
               cell = m%f(i)%c1
               c%cf(fill(cell)) = i; fill(cell) = fill(cell) + 1
            end if
         end do
      end block faces_fill

      ! ---- boundary face list ------------------------------------------------
      c%nbf = 0
      do i = 1, m%nfaces
         if ( m%f(i)%c1 == 0 ) c%nbf = c%nbf + 1
      end do
      allocate( c%bfaces(c%nbf) )
      n = 0
      do i = 1, m%nfaces
         if ( m%f(i)%c1 == 0 ) then
            n = n + 1
            c%bfaces(n) = i
         end if
      end do

      ! ---- cell -> neighbor cells (interior faces only) ----------------------
      allocate( c%c2c_ptr(m%ncells+1) )
      c%c2c_ptr = 0
      do i = 1, m%nfaces
         if ( m%f(i)%c1 > 0 ) then
            c%c2c_ptr(m%f(i)%c0+1) = c%c2c_ptr(m%f(i)%c0+1) + 1
            c%c2c_ptr(m%f(i)%c1+1) = c%c2c_ptr(m%f(i)%c1+1) + 1
         end if
      end do
      do k = 1, m%ncells
         c%c2c_ptr(k+1) = c%c2c_ptr(k+1) + c%c2c_ptr(k)
      end do
      allocate( c%c2c(c%c2c_ptr(m%ncells+1)) )

      neigh_fill: block
         integer, allocatable :: fill(:)
         allocate( fill(m%ncells) )
         fill = c%c2c_ptr(1:m%ncells) + 1
         do i = 1, m%nfaces
            if ( m%f(i)%c1 > 0 ) then
               cell = m%f(i)%c0
               c%c2c(fill(cell)) = m%f(i)%c1; fill(cell) = fill(cell) + 1
               cell = m%f(i)%c1
               c%c2c(fill(cell)) = m%f(i)%c0; fill(cell) = fill(cell) + 1
            end if
         end do
      end block neigh_fill

   end subroutine build_connectivity

end module mod_uns_connectivity
