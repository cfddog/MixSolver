!===============================================================================
! main_mpi.f90 -- Parallel UNSSolver driver (Steps 5 + 6)
!
! All ranks read the global mesh, control, and build connectivity/geometry.
! The mesh is partitioned via METIS (rank 0) and the partition is broadcast.
! Each rank extracts its local mesh, sets up halo exchange, builds local BCs
! and fields, then calls the parallel SIMPLE or PISO driver.
!
! After the solve, owned cells of each rank are gathered to rank 0, which
! holds the global mesh and writes the result file (VTU / Tecplot PLT) plus
! the Ghia et al. (1982) benchmark comparison.  PISO snapshots are gathered
! and written inside piso_run_mpi at every n_out_every steps.
!
! Usage: mpirun -np N ./bin/unsolver_mpi [mesh.cas] [control]
!===============================================================================
program unsolver_mpi
   use mod_precision, only: dp, ip, pi
   use mod_uns_mpi_core
   use mpi
   use mod_uns_mesh
   use mod_uns_cas_reader
   use mod_uns_connectivity
   use mod_uns_geometry
   use mod_uns_control
   use mod_uns_bc
   use mod_uns_fields
   use mod_uns_output
   use mod_uns_partition
   use mod_uns_local_mesh
   use mod_uns_halo
   use mod_uns_linsolver_mpi
   use mod_uns_simple_mpi
   use mod_uns_gather
   use mod_uns_restart
   use mod_uns_forces
   implicit none

   ! Global mesh data (read by all ranks)
   type(mesh_t)   :: m_g
   type(conn_t)   :: c_g
   type(geom_t)   :: g_g
   type(ctrl_t)   :: ctrl
   ! Local mesh + halo
   type(local_mesh_t) :: lm
   type(halo_info_t)  :: hi
   type(bc_t)         :: bcs_loc
   type(fields_t)      :: fld
   ! Global BC (rank 0 only, for force post-processing)
   type(bc_t)         :: bcs_glb
   ! Global fields on rank 0 (for output gathering)
   type(fields_t)      :: fld_g
   ! Partition
   integer, allocatable :: part(:)
   integer              :: edgecut
   ! Command line / output file
   character(len=512) :: casfile, ctlfile, arg, mapfile, stem, vtufile
   character(len=8)   :: out_ext
   character(len=256) :: src_dump
   integer :: ier, j, k, ios
   real(dp) :: umax_local, umax, pmax_local, pmax, pmin_local, pmin, mscale
   ! Force post-processing (rank 0 only)
   real(dp) :: Fx, Fy, Fz, CD, CL, u_inf, A_ref, p_inf
   character(len=512) :: cpfile

   ! ---- MPI init ----------------------------------------------------------
   call mpi_bootstrap()

   ! ---- command line -----------------------------------------------------
   casfile = 'Grid_BC/cavity.cas'
   ctlfile = 'cases/cavity.control'
   if ( command_argument_count() >= 1 ) then
      call get_command_argument( 1, arg ); casfile = trim(arg)
   end if
   if ( command_argument_count() >= 2 ) then
      call get_command_argument( 2, arg ); ctlfile = trim(arg)
   end if

   ! derive map file name from casfile basename
   j = scan( casfile, '/', back=.true. )
   stem = casfile(j+1:)
   k = scan( stem, '.', back=.true. )
   if ( k > 0 ) stem = stem(:k-1)
   mapfile = trim(stem)//'.part.map'

   if ( myrank == 0 ) then
      write(*,'(a)') '=============================================================='
      write(*,'(a)') '  UNSSolver-MPI -- parallel CFD solver'
      write(*,'(a)') '=============================================================='
      write(*,'(a,a)') '  Mesh file    : ', trim(casfile)
      write(*,'(a,a)') '  Control file : ', trim(ctlfile)
      write(*,'(a,i0,a,i0)') '  MPI ranks    : ', nprocs, '  (rank 0 reporting)'
      write(*,'(a,a)') '  Part map     : ', trim(mapfile)
   end if

   ! ---- all ranks: read mesh + connectivity + geometry --------------------
   call read_cas( trim(casfile), m_g, ier )
   if ( ier /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: mesh reading failed, ier = ', ier
      call mpi_shutdown(); stop 1
   end if

   ! optional unit scaling (e.g. mm -> m) before connectivity/geometry
   call read_mesh_scale( trim(ctlfile), mscale )
   if ( mscale /= 1.0_dp ) then
      m_g%x = m_g%x * mscale
      if ( myrank == 0 ) write(*,'(a,es10.3)') '  mesh coordinates scaled by ', mscale
   end if

   call build_connectivity( m_g, c_g, ier )
   if ( ier /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: connectivity failed, ier = ', ier
      call mpi_shutdown(); stop 1
   end if

   call compute_geometry( m_g, g_g, ier )
   if ( ier /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: geometry failed, ier = ', ier
      call mpi_shutdown(); stop 1
   end if

   ! ---- register CAS interface zones for coupling (phase 3 step B) -----------
   ! Only rank 0 has the global mesh and geometry, so registration happens
   ! on rank 0 only.  The Interface_List is populated on rank 0 for later
   ! use by the coupling layer.
   if ( myrank == 0 ) then
      write(*,'(a)') ''
      write(*,'(a)') '--- Coupling interface registration ---'
      call register_interface_zones( m_g, g_g )
   end if

   if ( myrank == 0 ) then
      write(*,'(a)') ''
      write(*,'(a)') '--- Global mesh summary ---'
      write(*,'(a,i10)') '  nodes : ', m_g%nnodes
      write(*,'(a,i10)') '  faces : ', m_g%nfaces
      write(*,'(a,i10)') '  cells : ', m_g%ncells
   end if

   ! ---- all ranks: read control -------------------------------------------
   call read_control( trim(ctlfile), ctrl, ier )
   if ( ier /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: control reading failed, ier = ', ier
      call mpi_shutdown(); stop 1
   end if

   ! ---- derive final output filename from stem + out_format ----------------
   if ( ctrl%out_format == 'tecplot' ) then
      out_ext = '.plt'
   else
      out_ext = '.vtu'
   end if
   vtufile = trim(stem)//trim(out_ext)
   if ( myrank == 0 ) then
      write(*,'(a,a)') '  Output file  : ', trim(vtufile)
   end if

   ! ---- phase 3 step C: resolve cell (volume) zones on every rank ---------
   ! (each rank holds the global mesh; only rank 0 prints the table)
   call resolve_cell_zones( m_g, ctrl, ier, verbose = (myrank == 0) )
   if ( ier /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: cell-zone resolution failed, ier = ', ier
      call mpi_shutdown(); stop 1
   end if

   ! ---- partition --------------------------------------------------------
   ! Restart with user-supplied partition_file: try to reuse if nprocs/ncells
   ! match, else re-partition from scratch.
   if ( ctrl%restart .and. len_trim(ctrl%partition_file) > 0 ) then
      call load_or_partition( c_g, nprocs, part, edgecut, &
                              trim(ctrl%partition_file), trim(mapfile), &
                              trim(stem), ier )
   else
      call partition_mesh_distributed( c_g, nprocs, part, edgecut, &
                                       trim(mapfile), trim(stem), ier )
   end if
   if ( ier /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: partition failed, ier = ', ier
      call mpi_shutdown(); stop 1
   end if

   if ( myrank == 0 ) then
      write(*,'(a,i0,a,i0,a)') '  Partitioned into ', nprocs, &
         ' zones, edgecut = ', edgecut, ' (METIS)'
   end if

   ! ---- extract local mesh + setup halo ----------------------------------
   call extract_local_mesh( m_g, c_g, g_g, part, lm, ier )
   if ( ier /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: local mesh extraction failed, ier = ', ier
      call mpi_shutdown(); stop 1
   end if

   call halo_setup( lm, part, hi, ier )
   if ( ier /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: halo setup failed, ier = ', ier
      call mpi_shutdown(); stop 1
   end if

   if ( myrank == 0 ) then
      write(*,'(a)') ''
      write(*,'(a)') '--- Local mesh info (rank 0) ---'
      write(*,'(a,i10)') '  local cells  : ', lm%m%ncells
      write(*,'(a,i10)') '    owned      : ', lm%nowned
      write(*,'(a,i10)') '    halo       ', lm%m%ncells - lm%nowned
      write(*,'(a,i10)') '  local faces  : ', lm%m%nfaces
   end if

   ! ---- build local BCs + fields -----------------------------------------
   ! allow_empty: a bc zone may have no face on this rank (owned by peers)
   call build_bc( lm%m, ctrl, bcs_loc, ier, allow_empty = .true. )
   if ( ier /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: local BC setup failed, ier = ', ier
      call mpi_shutdown(); stop 1
   end if

   call init_fields( lm%m, lm%g, fld )

   ! populate per-cell porous-medium coefficients from ctrl%cz
   call setup_porous_fields( lm%m, ctrl, fld )

   ! ---- initial temperature field (uniform + x-asymmetric perturbation) ----
   ! Global node extents (m_g%x) give the correct global mode phase on every
   ! rank; local cell centroids (lm%g%xc) position each local cell.
   if ( ctrl%init_T /= 0.0_dp .or. ctrl%t_pert /= 0.0_dp ) then
      block
         real(dp) :: xmin, xmax, lx, pi_v
         integer  :: ic
         pi_v = acos( -1.0_dp )
         xmin = minval( m_g%x(1,:) )
         xmax = maxval( m_g%x(1,:) )
         lx   = max( xmax - xmin, epsilon(1.0_dp) )
         do ic = 1, lm%m%ncells
            fld%T(ic) = ctrl%init_T + ctrl%t_pert &
                        * sin( pi_v * ( lm%g%xc(1,ic) - xmin ) / lx )
            fld%T_s(ic) = fld%T(ic)
         end do
         fld%T_old      = fld%T
         fld%T_old_old  = fld%T
         fld%T_s_old    = fld%T_s
         fld%T_s_old_old = fld%T_s
      end block
   end if

   ! ---- initialize velocity to free-stream (from far-field BC, if any) ----
   ! Otherwise SIMPLE starts from u=0 and the far-field inflow/outflow test
   ! (which keys on uP) sees zero normal velocity and returns uf=0, so the
   ! zero state is a trivial fixed point and the solver "converges" at it=1.
   block
      integer :: igrp, kk
      real(dp) :: u_inf_vec(3)
      u_inf_vec = 0.0_dp
      do igrp = 1, ctrl%nbc
         if ( ctrl%bc(igrp)%btype == BC_FARFIELD ) then
            u_inf_vec = ctrl%bc(igrp)%uvel
            exit
         end if
      end do
      if ( norm2(u_inf_vec) > 0.0_dp ) then
         do kk = 1, lm%nowned
            fld%u(:,kk) = u_inf_vec
         end do
         if ( myrank == 0 ) then
            write(*,'(a,3(es10.3,1x))') '  Initialized u to free-stream: ', u_inf_vec
         end if
      end if
   end block

   ! ---- allocate global fields container fld_g ------------------------------
   ! In restart mode, fld_g is allocated on ALL ranks by read_field_dump
   ! (so each rank can scatter_global_to_local into its local fld).
   ! Otherwise only rank 0 needs fld_g for output gathering.
   if ( ctrl%restart ) then
      if ( myrank == 0 ) then
         write(*,'(a)') ''
         write(*,'(a)') '--- Restart mode ---'
         write(*,'(a,a)') '  Dump file    : ', trim(ctrl%restart_file)
      end if
      call read_field_dump( trim(ctrl%restart_file), m_g%ncells, fld_g, &
                            src_dump, ier )
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: dump read failed, ier = ', ier
         call mpi_shutdown(); stop 1
      end if
      if ( myrank == 0 .and. len_trim(src_dump) > 0 ) then
         write(*,'(a,a)') '  Dump src mesh: ', trim(src_dump)
      end if
      call scatter_global_to_local( lm, fld_g, fld, ier )
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: field scatter failed, ier = ', ier
         call mpi_shutdown(); stop 1
      end if
      if ( myrank == 0 ) then
         write(*,'(a,i0)') '  Restart: scattered global fld -> local fld (nowned=', &
                           lm%nowned, ')'
      end if
   else
      if ( myrank == 0 ) call init_fields( m_g, g_g, fld_g )
   end if

   ! Populate porous coefficients of the rank-0 global work field so the
   ! VTU writer's LTNE gate (maxval(h_sf*a_sf) > 0) emits the solid-
   ! temperature block for the gathered field.
   if ( myrank == 0 ) call setup_porous_fields( m_g, ctrl, fld_g )

   ! ---- dispatch: steady SIMPLE, transient PISO or PIMPLE -----------------
   if ( ctrl%transient ) then
      if ( ctrl%pimple ) then
         call pimple_run_mpi( lm, hi, ctrl, bcs_loc, fld, trim(vtufile), &
                               m_g, c_g, g_g, part, fld_g, ier )
      else
         call piso_run_mpi( lm, hi, ctrl, bcs_loc, fld, trim(vtufile), &
                             m_g, c_g, g_g, part, fld_g, ier )
      end if
   else
      call simple_run_mpi( lm, hi, ctrl, bcs_loc, fld, ier )
   end if

   if ( ier /= 0 ) then
      if ( myrank == 0 ) write(*,'(a,i0)') 'FATAL: solver failed, ier = ', ier
      call mpi_shutdown(); stop 1
   end if

   ! ---- solution summary (global) -----------------------------------------
   umax_local = maxval( norm2( fld%u(:,1:lm%nowned), dim=1 ) )
   pmax_local = maxval( fld%p(1:lm%nowned) )
   pmin_local = minval( fld%p(1:lm%nowned) )
   call MPI_Allreduce( MPI_IN_PLACE, umax_local, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, mpi_comm, ier )
   call MPI_Allreduce( MPI_IN_PLACE, pmax_local, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, mpi_comm, ier )
   call MPI_Allreduce( MPI_IN_PLACE, pmin_local, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MIN, mpi_comm, ier )
   umax = umax_local; pmax = pmax_local; pmin = pmin_local

   if ( myrank == 0 ) then
      write(*,'(a)') ''
      write(*,'(a)') '--- Solution summary (MPI) ---'
      write(*,'(a,es12.5)') '  max |u|  : ', umax
      write(*,'(a,es14.6)') '  max p    : ', pmax
      write(*,'(a,es14.6)') '  min p    : ', pmin
      write(*,'(a,es12.5)') '  Tmin     : ', minval( fld%T(1:lm%nowned) )
      write(*,'(a,es12.5)') '  Tmax     : ', maxval( fld%T(1:lm%nowned) )
   end if

   ! ---- Step 6: gather to rank 0 + write final VTU/PLT + Ghia compare ----
   ! NB: PISO with n_out_every > 0 has already written per-step snapshots
   ! inside piso_run_mpi.  The final-state file <stem>.<ext> is always
   ! written here (mirroring the serial main driver).
   call write_snapshot_mpi( lm, fld, part, m_g, c_g, g_g, ctrl, &
                             trim(vtufile), fld_g, ios )
   if ( ios /= 0 .and. myrank == 0 ) then
      write(*,'(a,i0)') 'WARNING: final output write failed, ier=', ios
   end if

   if ( myrank == 0 ) then
      call ghia_compare( m_g, g_g, fld_g, &
                         ctrl%rho / max( ctrl%mu, tiny(1.0_dp) ), ier )

      ! ---- Force post-processing: build global BC, integrate pressure on
      ! wall zones, report CD/CL and write Cp distribution. Skips silently
      ! if no wall zone exists.
      call build_bc( m_g, ctrl, bcs_glb, ier )
      if ( ier == 0 ) then
         cpfile = trim(stem)//'_cp.dat'
         call compute_forces( m_g, g_g, ctrl, bcs_glb, fld_g, &
                              Fx, Fy, Fz, CD, CL, u_inf, A_ref, p_inf )
         call write_cp_distribution( m_g, g_g, ctrl, bcs_glb, fld_g, &
                                     trim(cpfile) )
         call compute_separation( m_g, g_g, ctrl, bcs_glb, fld_g )
      else
         write(*,'(a,i0)') 'WARNING: global BC build failed (force skip), ier=', ier
      end if
   end if

   ! ---- Step 7: write field dump for future restart ----------------------
   ! fld_g has been filled by write_snapshot_mpi's gather on rank 0, so the
   ! global state is already in place.  write_field_dump is a rank-0 no-op
   ! on other ranks.
   if ( len_trim(ctrl%dump_file) > 0 ) then
      if ( myrank == 0 ) then
         call write_field_dump( trim(ctrl%dump_file), m_g, fld_g, &
                                trim(casfile), ios )
         if ( ios /= 0 ) write(*,'(a,i0)') 'WARNING: dump write failed, ier=', ios
      end if
   end if

   call mpi_shutdown()

end program unsolver_mpi
