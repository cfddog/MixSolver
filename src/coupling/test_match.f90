!===============================================================================
! test_match.f90 -- phase 4 verification driver.
!
! Loads the structured mesh (control.ec + Mesh3d.x + bc3d.inp) and the
! unstructured mesh (unMesh.cas), registers both sides' coupling interface
! faces into the shared Interface_List, then runs match_interfaces and
! prints the matching report.  Does NOT advance either solver.
!
! Build: see Makefile target 'match_test'.
! Run:   cd grid_BC && mpirun -np 1 ../bin/match_test
!===============================================================================
program test_match
   use mpi
   use mod_struct_init, only: read_parameter, init
   use mod_uns_cas_reader, only: read_cas
   use mod_uns_connectivity, only: build_connectivity
   use mod_uns_geometry, only: compute_geometry, register_interface_zones
   use mod_interface_match, only: match_interfaces
   use mod_interface, only: Interface_List, Num_Interface
   use mod_uns_mesh, only: mesh_t
   use mod_uns_connectivity, only: conn_t
   use mod_uns_geometry, only: geom_t
   implicit none

   integer :: ierr, myrank
   type(mesh_t) :: m
   type(conn_t) :: conn
   type(geom_t) :: g

   call mpi_init(ierr)
   call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)

   ! ---- structured side: read control.ec + Mesh3d.x + bc3d.inp, register ----
   if (myrank == 0) write(*,'(a)') '=== phase-4 match test ==='
   call read_parameter
   call init
   if (myrank == 0) &
      write(*,'(a,i0)') '  structured interface entries registered: ', &
                         count(Interface_List(1:Num_Interface)%solver == 1)

   ! ---- unstructured side: read unMesh.cas, register ----
   call read_cas('unMesh.cas', m, ierr)
   if (ierr /= 0) then
      write(*,'(a,i0)') 'FATAL: read_cas failed, ier = ', ierr
      call mpi_abort(MPI_COMM_WORLD, 1, ierr)
   end if
   call build_connectivity(m, conn, ierr)
   call compute_geometry(m, g, ierr)
   call register_interface_zones(m, g)
   if (myrank == 0) &
      write(*,'(a,i0)') '  unstructured interface entries registered: ', &
                         count(Interface_List(1:Num_Interface)%solver == 2)

   ! ---- geometric matching ----
   if (myrank == 0) call match_interfaces

   call mpi_finalize(ierr)
end program test_match
