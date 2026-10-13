!===============================================================================
! mod_uns_driver_mpi.f90 -- METIS-partitioned PARALLEL unstructured driver.
!
! Entry point
! -----------
!   uns_run_mpi( comm, casfile, ctlfile )
!
! This restores the parallel unstructured path that used to live in the
! (phase-15 removed) src/unstructured/main_uns_mpi.f90, but callable from the
! single self-dispatching entry point:
!
!   mpirun -np N bin/mixnsolver_mpi        (mix.control:  mode = uns)
!
! Wiring (all ranks execute this routine; one partition per rank)
! ---------------------------------------------------------------
!   1. every rank reads the GLOBAL mesh (unMesh.cas) + control (unMesh.control),
!      builds connectivity/geometry and resolves the cell zones;
!   2. rank 0 partitions the CELLS with METIS (nparts = nprocs) and the
!      partition is broadcast; the map is cached in <stem>.part.map and (on a
!      restart, ctrl%partition_file) reloaded when it matches;
!   3. every rank extracts its LOCAL mesh (owned cells + one halo layer), sets
!      up halo exchange and builds LOCAL boundary conditions;
!   4. the parallel SIMPLE / PISO / PIMPLE drivers advance the local domain;
!   5. rank 0 gathers the owned cells, writes <stem>.vtu and runs the
!      post-processing (Ghia comparison, force / Cp integration, field dump).
!
! The file name ends in `_mpi.f90`, so the Makefile keeps it out of the serial
! (non -DHAVE_MPI) tree and links it into the MPI tree only.
!
! NOTE: mod_uns_partition / mod_uns_halo / mod_uns_linsolver_mpi /
! mod_uns_gather issue their collectives on MPI_COMM_WORLD, so this driver is
! restricted to the world communicator -- it is wired for the standalone
! `mode = uns` run only (there the uns group owns every rank).  The coupled
! mode keeps the full-mesh replication path (mod_mix_driver::uns_group_driver).
!===============================================================================
module mod_uns_driver_mpi
   use mpi
   use mod_precision, only: dp
   use mod_uns_mpi_core, only: myrank, nprocs, mpi_set_comm
   use mod_uns_mesh, only: mesh_t
   use mod_uns_cas_reader, only: read_cas
   use mod_uns_connectivity, only: conn_t, build_connectivity
   use mod_uns_geometry, only: geom_t, compute_geometry, register_interface_zones, &
                               tag_interface_cell_zones
   use mod_uns_control, only: ctrl_t, read_control, read_mesh_scale, &
                              resolve_cell_zones, BC_INTERFACE
   use mod_uns_bc, only: bc_t, build_bc
   use mod_uns_fields, only: fields_t, init_fields, setup_porous_fields, is_solid_cell
   use mod_uns_partition, only: partition_mesh_distributed
   use mod_uns_local_mesh, only: local_mesh_t, extract_local_mesh
   use mod_uns_halo, only: halo_info_t, halo_setup
   use mod_uns_simple_mpi, only: simple_run_mpi, piso_run_mpi, pimple_run_mpi
   use mod_uns_simple, only: NAN_IER
   use mod_uns_gather, only: write_snapshot_mpi
   use mod_uns_restart, only: read_field_dump, scatter_global_to_local, &
                              load_or_partition, write_field_dump
   use mod_uns_output, only: ghia_compare
   use mod_uns_forces, only: compute_forces, write_cp_distribution, compute_separation
   use mod_interface, only: report_interface_dispatch, PEER_STRUCT, IFACE_CZ_FLUID
   implicit none
   private
   public :: uns_run_mpi

contains

   !---------------------------------------------------------------------------
   ! uns_run_mpi -- parallel unstructured solve (see module header).
   !
   ! comm    : communicator holding ALL unstructured ranks (MPI_COMM_WORLD in
   !           the standalone `mode = uns` run)
   ! casfile : mesh file, e.g. 'unMesh.cas'
   ! ctlfile : control file, e.g. 'unMesh.control'
   !---------------------------------------------------------------------------
   subroutine uns_run_mpi( comm, casfile, ctlfile )
      integer,          intent(in) :: comm
      character(len=*), intent(in) :: casfile, ctlfile

      ! --- global mesh/geometry (every rank holds the full mesh) -------------
      type(mesh_t)   :: m_g
      type(conn_t)   :: c_g
      type(geom_t)   :: g_g
      type(ctrl_t)   :: ctrl
      ! --- local domain ------------------------------------------------------
      type(local_mesh_t) :: lm
      type(halo_info_t)  :: hi
      type(bc_t)         :: bcs_loc      ! local boundary conditions
      type(bc_t)         :: bcs_glb      ! rank-0 global BC (force post-proc)
      type(fields_t)     :: fld          ! local fields
      type(fields_t)     :: fld_g        ! global work fields (rank 0)
      ! --- partition ---------------------------------------------------------
      integer, allocatable :: part(:)
      integer :: edgecut
      ! --- misc --------------------------------------------------------------
      integer :: ier, ios, j, k
      real(dp) :: mscale
      real(dp) :: umax_l, pmax_l, pmin_l, umax, pmax, pmin
      real(dp) :: Fx, Fy, Fz, CD, CL, u_inf, A_ref, p_inf
      character(len=512) :: stem, mapfile, vtufile, src_dump, cpfile

      ier = 0

      ! ---- adopt the caller's already-initialised communicator --------------
      ! (main_dispatch called MPI_Init; a second MPI_Init would be fatal)
      if ( comm /= MPI_COMM_WORLD ) then
         write(*,'(a)') '[uns-mpi] FATAL: the partitioned uns driver only ' // &
                        'supports MPI_COMM_WORLD and needs >= 2 ranks'
         call MPI_Abort( comm, 6, ios )
         return
      end if
      call mpi_set_comm( comm )

      ! ---- 1. global mesh: every rank reads the same files ------------------
      call read_cas( trim(casfile), m_g, ier )
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: mesh read failed, ier=', ier
         call MPI_Abort( comm, 10, ios ); return
      end if

      ! optional unit scaling (e.g. mm -> m) before connectivity/geometry
      call read_mesh_scale( trim(ctlfile), mscale )
      if ( mscale /= 1.0_dp ) then
         m_g%x = m_g%x * mscale
         if ( myrank == 0 ) write(*,'(a,es10.3)') '[uns-mpi] mesh coordinates scaled by ', mscale
      end if

      call build_connectivity( m_g, c_g, ier )
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: connectivity failed, ier=', ier
         call MPI_Abort( comm, 11, ios ); return
      end if

      call compute_geometry( m_g, g_g, ier )
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: geometry failed, ier=', ier
         call MPI_Abort( comm, 12, ios ); return
      end if

      ! coupling-interface registration needs the global mesh + geometry, so it
      ! is done on rank 0 only (as in the removed main_uns_mpi)
      if ( myrank == 0 ) call register_interface_zones( m_g, g_g )

      call read_control( trim(ctlfile), ctrl, ier )
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: control read failed, ier=', ier
         call MPI_Abort( comm, 13, ios ); return
      end if

      ! phase 3 step C: resolve cell (volume) zones; only rank 0 prints the table
      call resolve_cell_zones( m_g, ctrl, ier, verbose = (myrank == 0) )
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: cell-zone resolution failed, ier=', ier
         call MPI_Abort( comm, 14, ios ); return
      end if

      if ( myrank == 0 ) then
         ! phase-11 dispatch: tag interface faces with the owning cell class and
         ! print the interface dispatch table (peer = structured solver)
         call tag_interface_cell_zones( m_g )
         call report_interface_dispatch( PEER_STRUCT, IFACE_CZ_FLUID )
      end if


      ! ---- output / partition file names from the mesh basename --------------
      stem = 'result'
      j = scan( casfile, '/', back=.true. )
      k = scan( casfile(j+1:), '.', back=.true. )
      if ( k > 1 ) stem = casfile(j+1:j+k-1)
      mapfile = trim(stem)//'.part.map'
      vtufile = trim(stem)//'.vtu'
      if ( myrank == 0 ) then
         write(*,'(a)')    '----------------------------------------------------------'
         write(*,'(a,i0,a)') '[uns-mpi] partitioned unstructured run on ', nprocs, ' ranks'
         write(*,'(a,a)')    '[uns-mpi] mesh          : ', trim(casfile)
         write(*,'(a,a)')    '[uns-mpi] control       : ', trim(ctlfile)
         write(*,'(a,a)')    '[uns-mpi] partition map : ', trim(mapfile)
         write(*,'(a,a)')    '[uns-mpi] output        : ', trim(vtufile)
         write(*,'(a,i0,a,i0)') '[uns-mpi] global mesh   : ', m_g%ncells, &
                              ' cells / ', m_g%nfaces
         write(*,'(a)')    '----------------------------------------------------------'
      end if

      ! ---- 2. partition the cells (rank 0 METIS, then broadcast) ------------
      ! Restart with a user-supplied partition_file: reuse the map if its rank
      ! count / mesh size still match, else re-partition from scratch.
      if ( ctrl%restart .and. len_trim(ctrl%partition_file) > 0 ) then
         call load_or_partition( c_g, nprocs, part, edgecut, &
                                 trim(ctrl%partition_file), trim(mapfile), &
                                 trim(stem), ier )
      else
         call partition_mesh_distributed( c_g, nprocs, part, edgecut, &
                                          trim(mapfile), trim(stem), ier )
      end if
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: partition failed, ier=', ier
         call MPI_Abort( comm, 15, ios ); return
      end if
      if ( myrank == 0 ) write(*,'(a,i0,a,i0)') '[uns-mpi] METIS parts = ', nprocs, &
                        ', edgecut = ', edgecut

      ! ---- 3. local mesh (owned + halo) + halo exchange ---------------------
      call extract_local_mesh( m_g, c_g, g_g, part, lm, ier )
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: local mesh extraction failed, ier=', ier
         call MPI_Abort( comm, 16, ios ); return
      end if

      call halo_setup( lm, part, hi, ier )
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: halo setup failed, ier=', ier
         call MPI_Abort( comm, 17, ios ); return
      end if

      if ( myrank == 0 ) write(*,'(a,i0,a,i0,a,i0,a)') '[uns-mpi] rank 0 local: ', &
                        lm%nowned, ' owned + ', lm%m%ncells - lm%nowned, &
                        ' halo cells / ', lm%m%nfaces, ' faces'

      ! ---- local BCs + local fields (mirrors uns_solver_init) ---------------
      ! allow_empty: a bc zone may own no face on this rank (its faces live on
      ! a peer rank), which is normal in the partitioned run.
      call build_bc( lm%m, ctrl, bcs_loc, ier, allow_empty = .true. )
      if ( ier /= 0 ) then
         if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: local BC setup failed, ier=', ier
         call MPI_Abort( comm, 18, ios ); return
      end if

      call init_fields( lm%m, lm%g, fld )
      call setup_porous_fields( lm%m, ctrl, fld )

      ! phase 17: a coupling-interface face whose owner is a PURE SOLID cell has
      ! no peer temperature to advect -> default it to an ADIABATIC (mode=1,
      ! q=0) patch, exactly as uns_solver_init does (a coupled run overwrites
      ! iface_q every coupling iteration via set_interface_q).
      block
         integer :: gi2, i2, c02
         logical :: any_solid
         any_solid = .false.
         do gi2 = 1, bcs_loc%nb
            if ( bcs_loc%gb(gi2)%btype /= BC_INTERFACE ) cycle
            do i2 = 1, bcs_loc%gb(gi2)%nf
               c02 = lm%m%f( bcs_loc%gb(gi2)%faces(i2) )%c0
               if ( is_solid_cell( lm%m, fld, c02 ) ) then
                  bcs_loc%iface_T_mode( bcs_loc%gb(gi2)%faces(i2) ) = 1
                  any_solid = .true.
               end if
            end do
         end do
         if ( any_solid .and. myrank == 0 ) write(*,'(a)') &
            '[uns-mpi] CHT interface: pure-solid faces default to adiabatic (Neumann)'
      end block

      ! optional uniform initial temperature (+ x-asymmetric perturbation).
      ! The global node extent gives the correct global mode phase on every
      ! rank; the local cell centroids position each local cell.
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
            fld%T_old       = fld%T
            fld%T_old_old   = fld%T
            fld%T_s_old     = fld%T_s
            fld%T_s_old_old = fld%T_s
         end block
      end if

      ! ---- global work field fld_g: rank-0 gather target (or restart state) --
      if ( ctrl%restart .and. len_trim(ctrl%restart_file) > 0 ) then
         ! read_field_dump is a collective (rank 0 reads + Bcast), so every rank
         ! takes part; fld_g is then allocated on ALL ranks and each rank seeds
         ! its local fld from it.
         if ( myrank == 0 ) write(*,'(a,a)') '[uns-mpi] restart from ', trim(ctrl%restart_file)
         call read_field_dump( trim(ctrl%restart_file), m_g%ncells, fld_g, &
                               src_dump, ier )
         if ( ier /= 0 ) then
            if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: dump read failed, ier=', ier
            call MPI_Abort( comm, 19, ios ); return
         end if
         call scatter_global_to_local( lm, fld_g, fld, ier )
         if ( ier /= 0 ) then
            if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: field scatter failed, ier=', ier
            call MPI_Abort( comm, 20, ios ); return
         end if
         if ( myrank == 0 ) write(*,'(a,i0,a)') '[uns-mpi] scattered global -> ', &
                           lm%nowned, ' owned local cells per rank'
      else
         if ( ctrl%restart .and. myrank == 0 ) write(*,'(a)') &
            '[uns-mpi] WARNING: restart=.true. but restart_file is empty -> cold start'
         ! rank 0 only needs fld_g as the gather target; the porous coefficients
         ! must be populated so the VTU writer's LTNE gate emits the solid block.
         if ( myrank == 0 ) then
            call init_fields( m_g, g_g, fld_g )
            call setup_porous_fields( m_g, ctrl, fld_g )
         end if
      end if

      ! ---- 4. solve: steady SIMPLE, transient PISO or PIMPLE ----------------
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
         if ( myrank == 0 ) write(*,'(a,i0)') '[uns-mpi] FATAL: solver failed, ier=', ier
         if ( ier == NAN_IER ) then
            if ( myrank == 0 ) write(*,'(a)') &
               '[uns-mpi] non-finite solution value -- see the nan-guard ' // &
               'report above (field / cell / iteration); the run is stopped here.'
            flush(6)
            call MPI_Abort( comm, NAN_IER, ios ); return
         end if
         flush(6)
         call MPI_Abort( comm, 21, ios ); return
      end if

      ! ---- global solution summary (owned cells only) -----------------------
      umax_l = maxval( norm2( fld%u(:,1:lm%nowned), dim = 1 ) )
      pmax_l = maxval( fld%p(1:lm%nowned) )
      pmin_l = minval( fld%p(1:lm%nowned) )
      call MPI_Allreduce( MPI_IN_PLACE, umax_l, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, comm, ier )
      call MPI_Allreduce( MPI_IN_PLACE, pmax_l, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, comm, ier )
      call MPI_Allreduce( MPI_IN_PLACE, pmin_l, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MIN, comm, ier )
      umax = umax_l; pmax = pmax_l; pmin = pmin_l

      if ( myrank == 0 ) then
         write(*,'(a)')        '[uns-mpi] --- solution summary (global) ---'
         write(*,'(a,es12.5)') '[uns-mpi]   max |u| : ', umax
         write(*,'(a,es14.6)') '[uns-mpi]   max p   : ', pmax
         write(*,'(a,es14.6)') '[uns-mpi]   min p   : ', pmin
         write(*,'(a,es12.5)') '[uns-mpi]   Tmin    : ', minval( fld%T(1:lm%nowned) )
         write(*,'(a,es12.5)') '[uns-mpi]   Tmax    : ', maxval( fld%T(1:lm%nowned) )
      end if

      ! ---- 5. gather to rank 0, write the result file -----------------------
      ! (transient PISO/PIMPLE with n_out_every > 0 has already written its
      ! per-step snapshots inside the solver; the final state is always written
      ! here, mirroring the serial driver.)
      call write_snapshot_mpi( lm, fld, part, m_g, c_g, g_g, trim(vtufile), &
                               fld_g, ios )
      if ( ios /= 0 .and. myrank == 0 ) write(*,'(a,i0)') &
         '[uns-mpi] WARNING: final output write failed, ier=', ios

      ! ---- rank-0 post-processing on the gathered global state --------------
      if ( myrank == 0 ) then
         call ghia_compare( m_g, g_g, fld_g, &
                            ctrl%rho / max( ctrl%mu, tiny(1.0_dp) ), ier )

         ! force post-processing: global BC, pressure integration on wall zones,
         ! CD/CL report + Cp distribution + separation detection
         call build_bc( m_g, ctrl, bcs_glb, ier )
         if ( ier == 0 ) then
            cpfile = trim(stem)//'_cp.dat'
            call compute_forces( m_g, g_g, ctrl, bcs_glb, fld_g, &
                                 Fx, Fy, Fz, CD, CL, u_inf, A_ref, p_inf )
            call write_cp_distribution( m_g, g_g, ctrl, bcs_glb, fld_g, &
                                        trim(cpfile) )
            call compute_separation( m_g, g_g, ctrl, bcs_glb, fld_g )
            if ( u_inf /= 0.0_dp ) write(*,'(a,es12.5,a,es12.5)') &
               '[uns-mpi]   Cd    : ', CD, '   Cl : ', CL
         else
            write(*,'(a,i0)') '[uns-mpi] WARNING: global BC build failed (force skip), ier=', ier
         end if
      end if

      ! ---- field dump for future restarts (rank-0 no-op on the peers) -------
      if ( len_trim(ctrl%dump_file) > 0 ) then
         if ( myrank == 0 ) then
            call write_field_dump( trim(ctrl%dump_file), m_g, fld_g, &
                                   trim(casfile), ios )
            if ( ios /= 0 ) write(*,'(a,i0)') &
               '[uns-mpi] WARNING: dump write failed, ier=', ios
         end if
      end if

      if ( myrank == 0 ) write(*,'(a)') '[uns-mpi] partitioned run finished'
   end subroutine uns_run_mpi

end module mod_uns_driver_mpi

