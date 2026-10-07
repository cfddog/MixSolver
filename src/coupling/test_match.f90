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
   use mod_uns_geometry, only: compute_geometry, register_interface_zones, &
                               tag_interface_cell_zones
   use mod_uns_control, only: read_control, resolve_cell_zones, ctrl_t
   use mod_interface_match, only: match_interfaces, dispatch_interfaces
   use mod_interface, only: Interface_List, Num_Interface, IFACE_COMP_FLUID_FLUID, &
                            IFACE_COMP_FLUID_POROUS
   use mod_uns_mesh, only: mesh_t
   use mod_uns_connectivity, only: conn_t
   use mod_uns_geometry, only: geom_t
   implicit none

   integer :: ierr, myrank, n_ff, n_fp
   type(mesh_t) :: m
   type(conn_t) :: conn
   type(geom_t) :: g
   type(ctrl_t) :: ctrl

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

   ! ---- cell-zone tagging (phase 11 dispatch needs m%cztype) ----
   call read_control('unMesh.control', ctrl, ierr)
   if (ierr /= 0) then
      write(*,'(a,i0)') 'FATAL: read_control failed, ier = ', ierr
      call mpi_abort(MPI_COMM_WORLD, 1, ierr)
   end if
   call resolve_cell_zones(m, ctrl, ierr)
   if (ierr /= 0) then
      write(*,'(a,i0)') 'FATAL: resolve_cell_zones failed, ier = ', ierr
      call mpi_abort(MPI_COMM_WORLD, 1, ierr)
   end if
   call tag_interface_cell_zones(m)

   ! ---- geometric matching + phase-11 dispatch ----
   if (myrank == 0) call match_interfaces
   if (myrank == 0) call dispatch_interfaces

   ! ---- assertion: grid_BC is an all-fluid case -> never classified porous ----
   n_ff = count(Interface_List(1:Num_Interface)%iface_type == IFACE_COMP_FLUID_FLUID)
   n_fp = count(Interface_List(1:Num_Interface)%iface_type == IFACE_COMP_FLUID_POROUS)
   if (myrank == 0) then
      write(*,'(a,i0,a,i0)') '  classified comp-fluid<->lowspeed-fluid: ', &
                             n_ff, ' / ', Num_Interface
      write(*,'(a,i0)')      '  classified comp-fluid<->lowspeed-porous: ', n_fp
   end if
   if (n_ff == 0 .or. n_fp /= 0) then
      if (myrank == 0) &
         write(*,'(a)') 'ASSERT FAILED: grid_BC must be all-fluid flow B'
      call mpi_abort(MPI_COMM_WORLD, 2, ierr)
   end if
   if (myrank == 0) write(*,'(a)') '  dispatch assertion PASSED'

   call mpi_finalize(ierr)
end program test_match
