!===============================================================================
! mod_gather.f90 -- Field gathering and parallel result output (Step 6)
!
! The local solver operates on owned + halo cells; output routines (vtk_write,
! tecplot_write, ghia_compare) need the global field on rank 0.  This module
! provides two helpers:
!
!   gather_fields_to_root:
!       Collects each rank's owned slice of u/p/T into rank 0's global
!       fields_t.  Each rank contributes exactly nowned contiguous entries
!       (its owned cells in local 1..nowned order).  Rank 0 places them
!       into the global field at the cell IDs indicated by the partition
!       vector (part(c) = owner rank), which reproduces the original
!       extract_local_mesh owned-global ordering deterministically.
!
!   write_snapshot_mpi:
!       Convenience wrapper: gathers fields, then rank 0 dispatches the
!       existing serial vtk_write or tecplot_write on the global mesh.
!       Non-root ranks participate in the gather but do no file I/O.
!
! NB: rank 0 must have already built m_g / c_g / g_g (the global mesh,
! connectivity and geometry) and have fld_g allocated to global size before
! the first call to either routine (see main_mpi).
!===============================================================================
module mod_uns_gather
   use mod_precision, only: dp, ip, pi
   use mod_uns_mpi_core, only: myrank, nprocs, mpi_check, mpi_comm
   use mpi
   use mod_uns_mesh
   use mod_uns_connectivity
   use mod_uns_geometry
   use mod_uns_control
   use mod_uns_fields
   use mod_uns_output
   use mod_uns_local_mesh, only: local_mesh_t
   implicit none
   private
   public :: gather_fields_to_root, write_snapshot_mpi

contains

   !----------------------------------------------------------------------------
   ! Gather u/p/T from each rank's owned cells to rank 0's global fld_g.
   ! All ranks call this collective routine.  fld_g is only meaningful on
   ! rank 0 after the call; non-root ranks may pass an unallocated fld_g.
   !----------------------------------------------------------------------------
   subroutine gather_fields_to_root( lm, fld, part, fld_g, ier )
      type(local_mesh_t), intent(in)    :: lm
      type(fields_t),     intent(in)    :: fld
      integer,            intent(in)    :: part(:)
      type(fields_t),     intent(inout) :: fld_g
      integer,            intent(out)   :: ier

      integer :: r, k, c, ierr
      integer :: nown, ntot
      integer, allocatable :: recvcounts(:), displs(:)
      real(dp), allocatable :: sbuf_u(:), sbuf_p(:), sbuf_T(:), sbuf_Ts(:)
      real(dp), allocatable :: rbuf_u(:), rbuf_p(:), rbuf_T(:), rbuf_Ts(:)

      ier = 0

      ! ---- per-rank owned counts and displacement table ----------------------
      ! (rank 0 owns the global partition vector; the k-th owned local cell
      ! of rank r is the k-th global cell c with part(c) == r, scanned in
      ! ascending c order -- this matches extract_local_mesh exactly)
      allocate( recvcounts(0:nprocs-1), displs(0:nprocs) )
      do r = 0, nprocs-1
         recvcounts(r) = count( part(1:size(part)) == r )
      end do
      displs(0) = 0
      do r = 1, nprocs
         displs(r) = displs(r-1) + recvcounts(r-1)
      end do
      nown = lm%nowned
      ntot = displs(nprocs)

      ! ---- pack owned slices (1..nowned) into contiguous send buffers -------
      allocate( sbuf_u(3*nown), sbuf_p(nown), sbuf_T(nown), sbuf_Ts(nown) )
      do k = 1, nown
         sbuf_u(3*k-2:3*k) = fld%u(:, k)
         sbuf_p(k)          = fld%p(k)
         sbuf_T(k)          = fld%T(k)
         sbuf_Ts(k)         = fld%T_s(k)
      end do

      ! ---- rank 0 allocates receive buffers; others pass null-size ---------
      if ( myrank == 0 ) then
         allocate( rbuf_u(3*ntot), rbuf_p(ntot), rbuf_T(ntot), rbuf_Ts(ntot) )
      else
         allocate( rbuf_u(0), rbuf_p(0), rbuf_T(0), rbuf_Ts(0) )
      end if

      call MPI_Gatherv( sbuf_u, 3*nown, MPI_DOUBLE_PRECISION, &
                        rbuf_u, 3*recvcounts, 3*displs, MPI_DOUBLE_PRECISION, &
                        0, mpi_comm, ierr )
      call mpi_check( ierr, 'gather u' )
      call MPI_Gatherv( sbuf_p, nown, MPI_DOUBLE_PRECISION, &
                        rbuf_p, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                        0, mpi_comm, ierr )
      call mpi_check( ierr, 'gather p' )
      call MPI_Gatherv( sbuf_T, nown, MPI_DOUBLE_PRECISION, &
                        rbuf_T, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                        0, mpi_comm, ierr )
      call mpi_check( ierr, 'gather T' )
      call MPI_Gatherv( sbuf_Ts, nown, MPI_DOUBLE_PRECISION, &
                        rbuf_Ts, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                        0, mpi_comm, ierr )
      call mpi_check( ierr, 'gather T_s' )

      ! ---- rank 0: scatter receive buffer into global fields by partition --
      if ( myrank == 0 ) then
         do r = 0, nprocs-1
            k = displs(r)            ! 0-based offset into rbuf
            do c = 1, size(part)     ! part is global, sized ncells_g
               if ( part(c) == r ) then
                  k = k + 1
                  fld_g%u(:, c) = rbuf_u(3*k-2:3*k)
                  fld_g%p(c)    = rbuf_p(k)
                  fld_g%T(c)    = rbuf_T(k)
                  fld_g%T_s(c)  = rbuf_Ts(k)
               end if
            end do
         end do
      end if

      deallocate( sbuf_u, sbuf_p, sbuf_T, sbuf_Ts, rbuf_u, rbuf_p, rbuf_T, &
                  rbuf_Ts, recvcounts, displs )
   end subroutine gather_fields_to_root

   !----------------------------------------------------------------------------
   ! Collective: gather fields and rank 0 writes VTU or PLT on the global mesh.
   ! fld_g is a persistent work array pre-allocated on rank 0 (unallocated on
   ! other ranks is fine -- they only participate in the gather).
   !----------------------------------------------------------------------------
   subroutine write_snapshot_mpi( lm, fld, part, m_g, c_g, g_g, ctrl, &
                                   filename, fld_g, ier )
      type(local_mesh_t), intent(in)    :: lm
      type(fields_t),     intent(in)    :: fld
      integer,            intent(in)    :: part(:)
      type(mesh_t),       intent(in)    :: m_g
      type(conn_t),       intent(in)    :: c_g
      type(geom_t),       intent(in)    :: g_g
      type(ctrl_t),       intent(in)    :: ctrl
      character(len=*),   intent(in)    :: filename
      type(fields_t),     intent(inout) :: fld_g
      integer,            intent(out)   :: ier

      integer :: ios

      ier = 0
      call gather_fields_to_root( lm, fld, part, fld_g, ier )
      if ( ier /= 0 ) return

      if ( myrank == 0 ) then
         if ( ctrl%out_format == 'tecplot' ) then
            call tecplot_write( m_g, c_g, g_g, ctrl, fld_g, trim(filename), ios )
         else
            call vtk_write( m_g, c_g, g_g, fld_g, trim(filename), ios )
         end if
         if ( ios /= 0 ) then
            write(*,'(a,a,a,i0)') 'WARNING: result write failed for ', &
               trim(filename), ' (continuing), ier=', ios
         else
            write(*,'(a,a)') '  wrote ', trim(filename)
         end if
      end if
   end subroutine write_snapshot_mpi

end module mod_uns_gather
