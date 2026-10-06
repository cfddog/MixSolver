!===============================================================================
! mod_simple_mpi.f90 -- Parallel SIMPLE/PISO drivers (Step 5)
!
! Reuses the assembly routines from mod_simple (momentum_assembly,
! correct_fields, temperature_assembly, flux_rhiechow, build_cell_csr,
! make_step_filename, pos_of) on the LOCAL mesh (owned + halo cells).
!
! Key differences from the serial drivers:
!   - Linear solves use mod_linsolver_mpi (bicgstab_ilu0_mpi / cg_jacobi_mpi /
!     cg_ic0_mpi) which operate on owned-sized vectors with halo exchange
!     inside the matvec.
!   - Cell-based fields (u, p, T, gu, gp, gt, apc, pp) are halo-exchanged at
!     six points each outer iteration so that the assembly routines see
!     up-to-date ghost-cell values.
!   - Convergence criteria (mass imbalance, du_max) use MPI_Allreduce to get
!     the global maximum.
!   - PPE pin cell is applied only on rank 0 (one global pin suffices to make
!     the pure-Neumann problem non-singular).
!   - Only rank 0 prints iteration history.
!===============================================================================
module mod_uns_simple_mpi
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   use mod_uns_connectivity
   use mod_uns_geometry
   use mod_uns_control
   use mod_uns_bc
   use mod_uns_fields
   use mod_uns_linsolver
   use mod_uns_output
   use mod_uns_mpi_core, only: myrank, nprocs, mpi_check, mpi_comm
   use mpi
   use mod_uns_local_mesh, only: local_mesh_t
   use mod_uns_halo, only: halo_info_t, halo_setup, halo_exchange_scalar, &
                       halo_exchange_vector
   use mod_uns_linsolver_mpi, only: bicgstab_ilu0_mpi, cg_jacobi_mpi, cg_ic0_mpi
   use mod_uns_gather, only: write_snapshot_mpi
   use mod_uns_simple, only: build_cell_csr, momentum_assembly, correct_fields, &
                            temperature_assembly, solid_temperature_assembly, &
                            flux_rhiechow, make_step_filename, &
                            outflow_mass_sums, outflow_mass_scale, &
                         pos_of
   implicit none
   private
   public :: simple_run_mpi, piso_run_mpi, pimple_run_mpi

contains

   !----------------------------------------------------------------------------
   ! Parallel SIMPLE driver
   !----------------------------------------------------------------------------
   subroutine simple_run_mpi( lm, hi, ctrl, bcs, fld, ier )
      type(local_mesh_t), intent(in)    :: lm
      type(halo_info_t),  intent(in)    :: hi
      type(ctrl_t),       intent(in)    :: ctrl
      type(bc_t),         intent(in)    :: bcs
      type(fields_t),     intent(inout) :: fld
      integer,            intent(out)   :: ier

      type(csr_t) :: amat, pmat
      real(dp), allocatable :: rhs(:), uo(:,:), pp(:), dpp(:,:), phif(:)
      real(dp), allocatable :: ap(:), imbalc(:), Tprev(:)
      real(dp) :: uscale, area_btot, flux_ref
      real(dp) :: imbal = 0.0_dp, du_max = 0.0_dp, dT_max = 0.0_dp
      real(dp) :: resl, ressum, gtmp
      integer  :: it, comp, itl, ierr, itl_max, i, k, c0
      integer, parameter :: PRINT_EVERY = 50

      ier = 0

      ! ---- matrix pattern (shared by momentum and pressure equations) --------
      call build_cell_csr( lm%m, lm%c, amat )
      pmat = amat
      allocate( rhs(lm%m%ncells), ap(lm%m%ncells), uo(3,lm%m%ncells) )
      allocate( pp(lm%m%ncells), dpp(3,lm%m%ncells), phif(lm%m%nfaces) )
      allocate( imbalc(lm%m%ncells) )
      allocate( Tprev(lm%m%ncells) )

      ! ---- reference scales (global) ------------------------------------------
      uscale = 0.0_dp
      area_btot = 0.0_dp
      do i = 1, bcs%nb
         uscale = max( uscale, norm2( bcs%gb(i)%lid_vel ), &
                               norm2( bcs%gb(i)%uvel ), &
                               bcs%gb(i)%uspeed )
         do k = 1, bcs%gb(i)%nf
            area_btot = area_btot + lm%g%area( bcs%gb(i)%faces(k) )
         end do
      end do
      if ( uscale <= 0.0_dp ) uscale = 1.0_dp
      if ( area_btot <= 0.0_dp ) area_btot = 1.0_dp
      call MPI_Allreduce( MPI_IN_PLACE, uscale, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, mpi_comm, ierr )
      call MPI_Allreduce( MPI_IN_PLACE, area_btot, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, mpi_comm, ierr )
      flux_ref = ctrl%rho * uscale * area_btot

      if ( myrank == 0 ) then
         write(*,'(a)') ''
         write(*,'(a)') '--- SIMPLE (MPI) iteration ---'
         write(*,'(a,es10.3,a,es10.3)') '  reference: uscale =', uscale, &
                                        '   flux_ref =', flux_ref
      end if

      fld%apc = 1.0_dp
      fld%ts_order = 0
      call compute_gradients( lm%m, lm%g, bcs, fld )

      ! ---- outer iterations ---------------------------------------------------
      it = 0
      do it = 1, ctrl%outer_max
         uo = fld%u
         itl_max = 0
         ressum = 0.0_dp

         ! mass-flow-inlet cold-start ramp (inlet_ramp<=1 disables it)
         if ( ctrl%inlet_ramp > 1 ) &
            call set_inlet_ramp_factor( real(it,dp) / real(ctrl%inlet_ramp,dp) )

         ! 0. halo exchange: u, p, T, gu, gp, gt before assembly ---------------
         call exchange_all_fields( hi, fld, ierr )

         ! 1. momentum equations -------------------------------------------------
         do comp = 1, 3
            call momentum_assembly( lm%m, lm%g, ctrl, bcs, fld, comp, amat, rhs, ap )
            call bicgstab_ilu0_mpi( amat, rhs(1:lm%nowned), fld%u(comp,:), &
                                    lm%nowned, hi, ctrl%lin_tol, ctrl%lin_max, &
                                    itl, resl, ierr )
            if ( ierr == 2 ) then
               if ( myrank == 0 ) &
                  write(*,'(a,i0,a,i0)') 'FATAL: BiCGSTAB breakdown, comp=', &
                                         comp, ', outer it=', it
               ier = 20; return
            end if
            itl_max = max( itl_max, itl )
            ressum = ressum + resl
         end do
         fld%apc = ap / ctrl%alpha_u

         ! 1b. exchange u + apc --------------------------------------------------
         call halo_exchange_vector( hi, fld%u, ierr )
         call mpi_check( ierr, 'simple u exchange' )
         call halo_exchange_scalar( hi, fld%apc, ierr )
         call mpi_check( ierr, 'simple apc exchange' )

         ! 2. gradients ----------------------------------------------------------
         call compute_gradients( lm%m, lm%g, bcs, fld )

         ! 2b. exchange gradients ------------------------------------------------
         call exchange_gradients( hi, fld, ierr )

         ! 3. Rhie-Chow face fluxes ---------------------------------------------
         call flux_rhiechow( lm%m, lm%g, ctrl, bcs, fld )
         call outflow_mass_rescale_mpi( lm, bcs, fld )

         ! 4. pressure-correction equation --------------------------------------
         call ppe_assembly_mpi( lm, ctrl, bcs, fld, pmat, rhs )
         pp = 0.0_dp
         select case ( trim(ctrl%ppe_precond) )
         case ( 'ic0' )
            call cg_ic0_mpi( pmat, rhs(1:lm%nowned), pp, lm%nowned, hi, &
                             ctrl%lin_tol, ctrl%lin_max, itl, resl, ierr )
         case default
            call cg_jacobi_mpi( pmat, rhs(1:lm%nowned), pp, lm%nowned, hi, &
                                ctrl%lin_tol, ctrl%lin_max, itl, resl, ierr )
         end select
         if ( ierr == 2 ) then
            if ( myrank == 0 ) &
               write(*,'(a,i0)') 'FATAL: CG breakdown in pressure correction, it=', it
            ier = 21; return
         end if
         itl_max = max( itl_max, itl )

         ! 4b. exchange pp -------------------------------------------------------
         call halo_exchange_scalar( hi, pp, ierr )
         call mpi_check( ierr, 'simple pp exchange' )

         ! 5. corrections --------------------------------------------------------
         call correct_fields( lm%m, lm%g, ctrl, bcs, fld, pp, phif, dpp )

         ! 5b. exchange u, p -----------------------------------------------------
         call halo_exchange_vector( hi, fld%u, ierr )
         call mpi_check( ierr, 'simple u2 exchange' )
         call halo_exchange_scalar( hi, fld%p, ierr )
         call mpi_check( ierr, 'simple p exchange' )

         ! 5c. temperature equation ---------------------------------------------
         Tprev = fld%T
         call temperature_assembly( lm%m, lm%g, ctrl, bcs, fld, amat, rhs, ap )
         call bicgstab_ilu0_mpi( amat, rhs(1:lm%nowned), fld%T, &
                                 lm%nowned, hi, ctrl%lin_tol, ctrl%lin_max, &
                                 itl, resl, ierr )
         if ( ierr == 2 ) then
            if ( myrank == 0 ) &
               write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in temperature, it=', it
            ier = 22; return
         end if
         itl_max = max( itl_max, itl )

         ! 5d. exchange T --------------------------------------------------------
         call halo_exchange_scalar( hi, fld%T, ierr )
         call mpi_check( ierr, 'simple T exchange' )

         ! 5e. solid temperature equation (LTNE only) ----------------------------
         if ( ctrl%thermal_model == 'ltne' ) then
            call solid_temperature_assembly( lm%m, lm%g, ctrl, bcs, fld, amat, rhs, ap )
            call bicgstab_ilu0_mpi( amat, rhs(1:lm%nowned), fld%T_s, &
                                    lm%nowned, hi, ctrl%lin_tol, ctrl%lin_max, &
                                    itl, resl, ierr )
            if ( ierr == 2 ) then
               if ( myrank == 0 ) &
                  write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in solid T, it=', it
               ier = 23; return
            end if
            itl_max = max( itl_max, itl )
            call halo_exchange_scalar( hi, fld%T_s, ierr )
            call mpi_check( ierr, 'simple Ts exchange' )
         end if


         ! 6. convergence measures (global) --------------------------------------
         imbalc = 0.0_dp
         do i = 1, lm%m%nfaces
            c0 = lm%m%f(i)%c0
            imbalc(c0) = imbalc(c0) + fld%flux(i)
            if ( lm%m%f(i)%c1 > 0 ) &
               imbalc(lm%m%f(i)%c1) = imbalc(lm%m%f(i)%c1) - fld%flux(i)
         end do
         imbal = maxval( abs( imbalc(1:lm%nowned) ) ) / flux_ref
         du_max = maxval( abs( fld%u(:,1:lm%nowned) - uo(:,1:lm%nowned) ) )
         ! temperature change per sweep: required for buoyancy-coupled runs
         dT_max = maxval( abs( fld%T(1:lm%nowned) - Tprev(1:lm%nowned) ) )
         call MPI_Allreduce( MPI_IN_PLACE, imbal, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )
         call MPI_Allreduce( MPI_IN_PLACE, du_max, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )
         call MPI_Allreduce( MPI_IN_PLACE, dT_max, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )

         if ( myrank == 0 ) then
            if ( mod(it,PRINT_EVERY) == 0 .or. it == 1 .or. &
                 ( imbal < ctrl%outer_tol .and. du_max < ctrl%outer_tol .and. &
                   dT_max < ctrl%outer_tol ) ) then
               gtmp = ressum / 3.0_dp
               write(*,'(a,i6,a,es10.3,a,es10.3,a,es10.3,a,i5,a,es9.2,a,es9.2,a,es9.2)') &
                  '  it=', it, '  mass-imbal=', imbal, '  du_max=', du_max, &
                  '  dT_max=', dT_max, &
                  '  lin-it=', itl_max, '  lin-res=', gtmp, &
                  '  Tmin=', minval(fld%T(1:lm%nowned)), &
                  '  Tmax=', maxval(fld%T(1:lm%nowned))
            end if
         end if

         if ( imbal < ctrl%outer_tol .and. du_max < ctrl%outer_tol .and. &
              dT_max < ctrl%outer_tol ) then
            if ( myrank == 0 ) &
               write(*,'(a,i0,a)') 'CONVERGED after ', it, ' outer iterations'
            exit
         end if
      end do

      if ( it >= ctrl%outer_max .and. &
           .not. ( imbal < ctrl%outer_tol .and. du_max < ctrl%outer_tol .and. &
                   dT_max < ctrl%outer_tol ) ) then
         if ( myrank == 0 ) &
            write(*,'(a)') 'WARNING: SIMPLE reached outer_max without converging'
      end if

      deallocate( rhs, ap, uo, pp, dpp, phif, imbalc, Tprev )

   end subroutine simple_run_mpi

   !----------------------------------------------------------------------------
   ! Parallel PISO driver (transient)
   !
   ! Snapshots are written every n_out_every steps and at the final step:
   ! fields are gathered from all ranks to rank 0, which then calls the
   ! existing serial vtk_write / tecplot_write on the global mesh.
   !----------------------------------------------------------------------------
   subroutine piso_run_mpi( lm, hi, ctrl, bcs, fld, vtu_prefix, &
                             m_g, c_g, g_g, part, fld_g, ier )
      type(local_mesh_t), intent(in)    :: lm
      type(halo_info_t),  intent(in)    :: hi
      type(ctrl_t),       intent(in)    :: ctrl
      type(bc_t),         intent(in)    :: bcs
      type(fields_t),     intent(inout) :: fld
      character(len=*),   intent(in)    :: vtu_prefix
      type(mesh_t),       intent(in)    :: m_g
      type(conn_t),       intent(in)    :: c_g
      type(geom_t),       intent(in)    :: g_g
      integer,            intent(in)    :: part(:)
      type(fields_t),     intent(inout) :: fld_g  ! pre-allocated on rank 0
      integer,            intent(out)   :: ier

      type(csr_t) :: amat, pmat
      real(dp), allocatable :: rhs(:), pp(:), dpp(:,:), phif(:)
      real(dp), allocatable :: ap(:), imbalc(:)
      real(dp) :: uscale, area_btot, flux_ref
      real(dp) :: imbal = 0.0_dp, du_max = 0.0_dp
      real(dp) :: resl, ressum, t_now, umax, cfl_max, dt_dx, gtmp
      integer  :: it, comp, itl, ierr, itl_max, icorr, i, k, c0
      integer  :: n_steps, n_out, ios
      real(dp) :: umax_local
      character(len=8)   :: step_ext
      character(len=512) :: vtufile
      integer, parameter :: PRINT_EVERY = 1

      ier = 0

      if ( .not. ctrl%transient ) then
         if ( myrank == 0 ) &
            write(*,'(a)') 'FATAL: piso_run_mpi invoked with transient=false'
         ier = 30; return
      end if
      if ( ctrl%dt <= 0.0_dp ) then
         if ( myrank == 0 ) &
            write(*,'(a,es12.4)') 'FATAL: dt must be > 0, got ', ctrl%dt
         ier = 31; return
      end if
      n_steps = max( 1, ctrl%n_time_max )
      n_out   = max( 0, ctrl%n_out_every )
      if ( ctrl%out_format == 'tecplot' ) then
         step_ext = '.plt'
      else
         step_ext = '.vtu'
      end if

      ! ---- matrix pattern ----------------------------------------------------
      call build_cell_csr( lm%m, lm%c, amat )
      pmat = amat
      allocate( rhs(lm%m%ncells), ap(lm%m%ncells), pp(lm%m%ncells) )
      allocate( dpp(3,lm%m%ncells), phif(lm%m%nfaces), imbalc(lm%m%ncells) )

      ! ---- reference scales (global) ------------------------------------------
      uscale = 0.0_dp
      area_btot = 0.0_dp
      do i = 1, bcs%nb
         uscale = max( uscale, norm2( bcs%gb(i)%lid_vel ), &
                               norm2( bcs%gb(i)%uvel ), &
                               bcs%gb(i)%uspeed )
         do k = 1, bcs%gb(i)%nf
            area_btot = area_btot + lm%g%area( bcs%gb(i)%faces(k) )
         end do
      end do
      if ( uscale <= 0.0_dp ) uscale = 1.0_dp
      if ( area_btot <= 0.0_dp ) area_btot = 1.0_dp
      call MPI_Allreduce( MPI_IN_PLACE, uscale, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, mpi_comm, ierr )
      call MPI_Allreduce( MPI_IN_PLACE, area_btot, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, mpi_comm, ierr )
      flux_ref = ctrl%rho * uscale * area_btot

      if ( myrank == 0 ) then
         write(*,'(a)') ''
         write(*,'(a)') '--- PISO (MPI) time stepping ---'
         write(*,'(a,es10.3,a,i0,a,i0)') '  dt=', ctrl%dt, '  n_steps=', &
                                          n_steps, '  n_correct=', ctrl%n_correct
         write(*,'(a,es10.3,a,es10.3)') '  reference: uscale =', uscale, &
                                        '   flux_ref =', flux_ref
      end if

      fld%apc = 1.0_dp
      call compute_gradients( lm%m, lm%g, bcs, fld )
      call flux_rhiechow( lm%m, lm%g, ctrl, bcs, fld )

      t_now = 0.0_dp

      ! ---- time stepping ------------------------------------------------------
      do it = 1, n_steps

         ! 1. advance history + select time scheme ------------------------------
         fld%u_old_old = fld%u_old
         fld%u_old     = fld%u
         fld%T_old_old = fld%T_old
         fld%T_old     = fld%T
         fld%T_s_old_old = fld%T_s_old
         fld%T_s_old     = fld%T_s
         if ( it == 1 .or. ctrl%time_scheme == 1 ) then
            fld%ts_order = 1
         else
            fld%ts_order = 2
         end if

         ! 0. halo exchange before momentum assembly ---------------------------
         call exchange_all_fields( hi, fld, ierr )

         ! 2. predict momentum --------------------------------------------------
         itl_max = 0
         ressum = 0.0_dp
         do comp = 1, 3
            call momentum_assembly( lm%m, lm%g, ctrl, bcs, fld, comp, amat, rhs, ap )
            call bicgstab_ilu0_mpi( amat, rhs(1:lm%nowned), fld%u(comp,:), &
                                    lm%nowned, hi, ctrl%lin_tol, ctrl%lin_max, &
                                    itl, resl, ierr )
            if ( ierr == 2 ) then
               if ( myrank == 0 ) &
                  write(*,'(a,i0,a,i0)') 'FATAL: BiCGSTAB breakdown, comp=', &
                                         comp, ', step=', it
               ier = 33; return
            end if
            itl_max = max( itl_max, itl )
            ressum = ressum + resl
         end do
         fld%apc = ap

         ! 2b. exchange u + apc -------------------------------------------------
         call halo_exchange_vector( hi, fld%u, ierr )
         call mpi_check( ierr, 'piso u exchange' )
         call halo_exchange_scalar( hi, fld%apc, ierr )
         call mpi_check( ierr, 'piso apc exchange' )

         ! 3. gradients + Rhie-Chow fluxes --------------------------------------
         call compute_gradients( lm%m, lm%g, bcs, fld )
         call exchange_gradients( hi, fld, ierr )
         call flux_rhiechow( lm%m, lm%g, ctrl, bcs, fld )
         call outflow_mass_rescale_mpi( lm, bcs, fld )

         ! 4. n_correct pressure-correction sweeps ------------------------------
         do icorr = 1, ctrl%n_correct
            call ppe_assembly_mpi( lm, ctrl, bcs, fld, pmat, rhs )
            pp = 0.0_dp
            select case ( trim(ctrl%ppe_precond) )
            case ( 'ic0' )
               call cg_ic0_mpi( pmat, rhs(1:lm%nowned), pp, lm%nowned, hi, &
                                ctrl%lin_tol, ctrl%lin_max, itl, resl, ierr )
            case default
               call cg_jacobi_mpi( pmat, rhs(1:lm%nowned), pp, lm%nowned, hi, &
                                   ctrl%lin_tol, ctrl%lin_max, itl, resl, ierr )
            end select
            if ( ierr == 2 ) then
               if ( myrank == 0 ) &
                  write(*,'(a,i0,a,i0)') 'FATAL: CG breakdown in PPE, sweep=', &
                                         icorr, ', step=', it
               ier = 34; return
            end if
            itl_max = max( itl_max, itl )

            ! exchange pp before correct_fields
            call halo_exchange_scalar( hi, pp, ierr )
            call mpi_check( ierr, 'piso pp exchange' )

            call correct_fields( lm%m, lm%g, ctrl, bcs, fld, pp, phif, dpp )

            ! exchange u, p after corrections
            call halo_exchange_vector( hi, fld%u, ierr )
            call mpi_check( ierr, 'piso u2 exchange' )
            call halo_exchange_scalar( hi, fld%p, ierr )
            call mpi_check( ierr, 'piso p2 exchange' )
         end do

         ! 4b. temperature equation ---------------------------------------------
         call temperature_assembly( lm%m, lm%g, ctrl, bcs, fld, amat, rhs, ap )
         call bicgstab_ilu0_mpi( amat, rhs(1:lm%nowned), fld%T, &
                                 lm%nowned, hi, ctrl%lin_tol, ctrl%lin_max, &
                                 itl, resl, ierr )
         if ( ierr == 2 ) then
            if ( myrank == 0 ) &
               write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in temperature, step=', it
            ier = 35; return
         end if
         itl_max = max( itl_max, itl )

         ! exchange T
         call halo_exchange_scalar( hi, fld%T, ierr )
         call mpi_check( ierr, 'piso T exchange' )

         ! solid temperature (LTNE only)
         if ( ctrl%thermal_model == 'ltne' ) then
            call solid_temperature_assembly( lm%m, lm%g, ctrl, bcs, fld, amat, rhs, ap )
            call bicgstab_ilu0_mpi( amat, rhs(1:lm%nowned), fld%T_s, &
                                    lm%nowned, hi, ctrl%lin_tol, ctrl%lin_max, &
                                    itl, resl, ierr )
            if ( ierr == 2 ) then
               if ( myrank == 0 ) &
                  write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in solid T, step=', it
               ier = 36; return
            end if
            itl_max = max( itl_max, itl )
            call halo_exchange_scalar( hi, fld%T_s, ierr )
            call mpi_check( ierr, 'piso Ts exchange' )
         end if

         ! 5. advance time -------------------------------------------------------
         t_now = t_now + ctrl%dt

         ! ---- diagnostics (global) --------------------------------------------
         imbalc = 0.0_dp
         do i = 1, lm%m%nfaces
            c0 = lm%m%f(i)%c0
            imbalc(c0) = imbalc(c0) + fld%flux(i)
            if ( lm%m%f(i)%c1 > 0 ) &
               imbalc(lm%m%f(i)%c1) = imbalc(lm%m%f(i)%c1) - fld%flux(i)
         end do
         imbal = maxval( abs( imbalc(1:lm%nowned) ) ) / flux_ref
         du_max = maxval( abs( fld%u(:,1:lm%nowned) - fld%u_old(:,1:lm%nowned) ) ) / ctrl%dt
         umax_local = maxval( norm2( fld%u(:,1:lm%nowned), dim=1 ) )
         call MPI_Allreduce( MPI_IN_PLACE, imbal, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )
         call MPI_Allreduce( MPI_IN_PLACE, du_max, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )
         call MPI_Allreduce( MPI_IN_PLACE, umax_local, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )
         umax = umax_local

         cfl_max = 0.0_dp
         do i = 1, lm%nowned
            dt_dx = umax * ctrl%dt / ( lm%g%vol(i) ** (1.0_dp/3.0_dp) )
            if ( dt_dx > cfl_max ) cfl_max = dt_dx
         end do
         call MPI_Allreduce( MPI_IN_PLACE, cfl_max, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )

         if ( myrank == 0 ) then
            if ( mod(it-1,PRINT_EVERY) == 0 .or. it == 1 .or. it == n_steps ) then
               gtmp = ressum / 3.0_dp
               write(*,'(a,i6,a,es10.3,a,es9.2,a,es9.2,a,es9.2,a,es9.2,a,i5,a,es9.2)') &
                  '  step=', it, '  t=', t_now, '  mass-imbal=', imbal, &
                  '  du/dt=', du_max, '  max|u|=', umax, '  CFL~=', cfl_max, &
                  '  lin-it=', itl_max, '  lin-res=', gtmp
               write(*,'(a,es9.2,a,es9.2)') '    Tmin=', &
                  minval(fld%T(1:lm%nowned)), '  Tmax=', &
                  maxval(fld%T(1:lm%nowned))
            end if
         end if

         ! ---- optional snapshot (Step 6): gather to rank 0 + global write ----
         if ( n_out > 0 .and. ( mod(it,n_out) == 0 .or. it == n_steps ) ) then
            call make_step_filename( vtu_prefix, it, step_ext, vtufile )
            call write_snapshot_mpi( lm, fld, part, m_g, c_g, g_g, ctrl, &
                                     trim(vtufile), fld_g, ios )
            if ( ios /= 0 .and. myrank == 0 ) then
               write(*,'(a,i0,a)') 'WARNING: snapshot write failed at step ', &
                                  it, ' (continuing)'
            end if
         end if

      end do

      if ( myrank == 0 ) &
         write(*,'(a,i0,a,es12.4)') 'PISO (MPI) finished: ', n_steps, &
                                     ' time steps, t_end =', t_now

      deallocate( rhs, ap, pp, dpp, phif, imbalc )

   end subroutine piso_run_mpi

   !----------------------------------------------------------------------------
   ! PIMPLE (MPI) driver -- phase 11.
   ! Same algorithm as the serial pimple_run: outer iterations per time step,
   ! momentum under-relaxation (alpha_u) on top of the implicit time term,
   ! full (alpha_p=1) pressure corrections.  Halo exchanges follow the same
   ! pattern as piso_run_mpi, repeated inside the outer loop.
   !----------------------------------------------------------------------------
   subroutine pimple_run_mpi( lm, hi, ctrl, bcs, fld, vtu_prefix, &
                              m_g, c_g, g_g, part, fld_g, ier )
      type(local_mesh_t), intent(in)    :: lm
      type(halo_info_t),  intent(in)    :: hi
      type(ctrl_t),       intent(in)    :: ctrl
      type(bc_t),         intent(in)    :: bcs
      type(fields_t),     intent(inout) :: fld
      character(len=*),   intent(in)    :: vtu_prefix
      type(mesh_t),       intent(in)    :: m_g
      type(conn_t),       intent(in)    :: c_g
      type(geom_t),       intent(in)    :: g_g
      integer,            intent(in)    :: part(:)
      type(fields_t),     intent(inout) :: fld_g
      integer,            intent(out)   :: ier

      type(csr_t) :: amat, pmat
      real(dp), allocatable :: rhs(:), pp(:), dpp(:,:), phif(:)
      real(dp), allocatable :: ap(:), imbalc(:)
      real(dp) :: uscale, area_btot, flux_ref
      real(dp) :: imbal = 0.0_dp, du_max = 0.0_dp
      real(dp) :: resl, ressum, t_now, umax, cfl_max, dt_dx, gtmp
      integer  :: it, comp, itl, ierr, itl_max, icorr, iouter, i, k, c0
      integer  :: n_steps, n_out, n_outer, ios
      real(dp) :: umax_local
      character(len=8)   :: step_ext
      character(len=512) :: vtufile
      integer, parameter :: PRINT_EVERY = 1

      ier = 0

      if ( .not. ctrl%transient ) then
         if ( myrank == 0 ) &
            write(*,'(a)') 'FATAL: pimple_run_mpi invoked with transient=false'
         ier = 50; return
      end if
      if ( ctrl%dt <= 0.0_dp ) then
         if ( myrank == 0 ) &
            write(*,'(a,es12.4)') 'FATAL: dt must be > 0, got ', ctrl%dt
         ier = 51; return
      end if
      n_steps = max( 1, ctrl%n_time_max )
      n_outer = max( 1, ctrl%n_outer_iter )
      n_out   = max( 0, ctrl%n_out_every )
      if ( ctrl%out_format == 'tecplot' ) then
         step_ext = '.plt'
      else
         step_ext = '.vtu'
      end if

      call build_cell_csr( lm%m, lm%c, amat )
      pmat = amat
      allocate( rhs(lm%m%ncells), ap(lm%m%ncells), pp(lm%m%ncells) )
      allocate( dpp(3,lm%m%ncells), phif(lm%m%nfaces), imbalc(lm%m%ncells) )

      uscale = 0.0_dp
      area_btot = 0.0_dp
      do i = 1, bcs%nb
         uscale = max( uscale, norm2( bcs%gb(i)%lid_vel ), &
                               norm2( bcs%gb(i)%uvel ), &
                               bcs%gb(i)%uspeed )
         do k = 1, bcs%gb(i)%nf
            area_btot = area_btot + lm%g%area( bcs%gb(i)%faces(k) )
         end do
      end do
      if ( uscale <= 0.0_dp ) uscale = 1.0_dp
      if ( area_btot <= 0.0_dp ) area_btot = 1.0_dp
      call MPI_Allreduce( MPI_IN_PLACE, uscale, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, mpi_comm, ierr )
      call MPI_Allreduce( MPI_IN_PLACE, area_btot, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, mpi_comm, ierr )
      flux_ref = ctrl%rho * uscale * area_btot

      if ( myrank == 0 ) then
         write(*,'(a)') ''
         write(*,'(a)') '--- PIMPLE (MPI) time stepping ---'
         write(*,'(a,es10.3,a,i0,a,i0,a,i0)') '  dt=', ctrl%dt, '  n_steps=', &
               n_steps, '  n_outer_iter=', n_outer, '  n_correct=', ctrl%n_correct
         write(*,'(a,es10.3,a,es10.3)') '  alpha_u=', ctrl%alpha_u, &
               '   flux_ref =', flux_ref
      end if

      fld%apc = 1.0_dp
      call compute_gradients( lm%m, lm%g, bcs, fld )
      call flux_rhiechow( lm%m, lm%g, ctrl, bcs, fld )

      t_now = 0.0_dp

      do it = 1, n_steps

         ! history advance once per step (outside outer loop)
         fld%u_old_old = fld%u_old
         fld%u_old     = fld%u
         fld%T_old_old = fld%T_old
         fld%T_old     = fld%T
         fld%T_s_old_old = fld%T_s_old
         fld%T_s_old     = fld%T_s
         if ( it == 1 .or. ctrl%time_scheme == 1 ) then
            fld%ts_order = 1
         else
            fld%ts_order = 2
         end if

         call exchange_all_fields( hi, fld, ierr )

         itl_max = 0
         ressum  = 0.0_dp

         ! ---- outer iterations ------------------------------------------------
         do iouter = 1, n_outer

            do comp = 1, 3
               call momentum_assembly( lm%m, lm%g, ctrl, bcs, fld, comp, amat, rhs, ap )
               call bicgstab_ilu0_mpi( amat, rhs(1:lm%nowned), fld%u(comp,:), &
                                       lm%nowned, hi, ctrl%lin_tol, ctrl%lin_max, &
                                       itl, resl, ierr )
               if ( ierr == 2 ) then
                  if ( myrank == 0 ) &
                     write(*,'(a,i0,a,i0)') 'FATAL: BiCGSTAB breakdown, comp=', &
                                            comp, ', step=', it
                  ier = 53; return
               end if
               itl_max = max( itl_max, itl )
               ressum = ressum + resl
            end do
            fld%apc = ap / ctrl%alpha_u   ! relaxed diag for Rhie-Chow

            call halo_exchange_vector( hi, fld%u, ierr )
            call mpi_check( ierr, 'pimple u exchange' )
            call halo_exchange_scalar( hi, fld%apc, ierr )
            call mpi_check( ierr, 'pimple apc exchange' )

            call compute_gradients( lm%m, lm%g, bcs, fld )
            call exchange_gradients( hi, fld, ierr )
            call flux_rhiechow( lm%m, lm%g, ctrl, bcs, fld )
            call outflow_mass_rescale_mpi( lm, bcs, fld )

            do icorr = 1, ctrl%n_correct
               call ppe_assembly_mpi( lm, ctrl, bcs, fld, pmat, rhs )
               pp = 0.0_dp
               select case ( trim(ctrl%ppe_precond) )
               case ( 'ic0' )
                  call cg_ic0_mpi( pmat, rhs(1:lm%nowned), pp, lm%nowned, hi, &
                                   ctrl%lin_tol, ctrl%lin_max, itl, resl, ierr )
               case default
                  call cg_jacobi_mpi( pmat, rhs(1:lm%nowned), pp, lm%nowned, hi, &
                                      ctrl%lin_tol, ctrl%lin_max, itl, resl, ierr )
               end select
               if ( ierr == 2 ) then
                  if ( myrank == 0 ) &
                     write(*,'(a,i0,a,i0)') 'FATAL: CG breakdown in PPE, sweep=', &
                                            icorr, ', step=', it
                  ier = 54; return
               end if
               itl_max = max( itl_max, itl )

               call halo_exchange_scalar( hi, pp, ierr )
               call mpi_check( ierr, 'pimple pp exchange' )
               call correct_fields( lm%m, lm%g, ctrl, bcs, fld, pp, phif, dpp )

               call halo_exchange_vector( hi, fld%u, ierr )
               call mpi_check( ierr, 'pimple u2 exchange' )
               call halo_exchange_scalar( hi, fld%p, ierr )
               call mpi_check( ierr, 'pimple p2 exchange' )
            end do

         end do   ! outer

         ! temperature (once per step, with final corrected velocity)
         call temperature_assembly( lm%m, lm%g, ctrl, bcs, fld, amat, rhs, ap )
         call bicgstab_ilu0_mpi( amat, rhs(1:lm%nowned), fld%T, &
                                 lm%nowned, hi, ctrl%lin_tol, ctrl%lin_max, &
                                 itl, resl, ierr )
         if ( ierr == 2 ) then
            if ( myrank == 0 ) &
               write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in T, step=', it
            ier = 55; return
         end if
         itl_max = max( itl_max, itl )
         call halo_exchange_scalar( hi, fld%T, ierr )
         call mpi_check( ierr, 'pimple T exchange' )

         ! solid temperature (LTNE only)
         if ( ctrl%thermal_model == 'ltne' ) then
            call solid_temperature_assembly( lm%m, lm%g, ctrl, bcs, fld, amat, rhs, ap )
            call bicgstab_ilu0_mpi( amat, rhs(1:lm%nowned), fld%T_s, &
                                    lm%nowned, hi, ctrl%lin_tol, ctrl%lin_max, &
                                    itl, resl, ierr )
            if ( ierr == 2 ) then
               if ( myrank == 0 ) &
                  write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in solid T, step=', it
               ier = 56; return
            end if
            itl_max = max( itl_max, itl )
            call halo_exchange_scalar( hi, fld%T_s, ierr )
            call mpi_check( ierr, 'pimple Ts exchange' )
         end if

         t_now = t_now + ctrl%dt

         ! diagnostics (global)
         imbalc = 0.0_dp
         do i = 1, lm%m%nfaces
            c0 = lm%m%f(i)%c0
            imbalc(c0) = imbalc(c0) + fld%flux(i)
            if ( lm%m%f(i)%c1 > 0 ) &
               imbalc(lm%m%f(i)%c1) = imbalc(lm%m%f(i)%c1) - fld%flux(i)
         end do
         imbal = maxval( abs( imbalc(1:lm%nowned) ) ) / flux_ref
         du_max = maxval( abs( fld%u(:,1:lm%nowned) - fld%u_old(:,1:lm%nowned) ) ) / ctrl%dt
         umax_local = maxval( norm2( fld%u(:,1:lm%nowned), dim=1 ) )
         call MPI_Allreduce( MPI_IN_PLACE, imbal, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )
         call MPI_Allreduce( MPI_IN_PLACE, du_max, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )
         call MPI_Allreduce( MPI_IN_PLACE, umax_local, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )
         umax = umax_local

         cfl_max = 0.0_dp
         do i = 1, lm%nowned
            dt_dx = umax * ctrl%dt / ( lm%g%vol(i) ** (1.0_dp/3.0_dp) )
            if ( dt_dx > cfl_max ) cfl_max = dt_dx
         end do
         call MPI_Allreduce( MPI_IN_PLACE, cfl_max, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MAX, mpi_comm, ierr )

         if ( myrank == 0 ) then
            if ( mod(it-1,PRINT_EVERY) == 0 .or. it == 1 .or. it == n_steps ) then
               gtmp = ressum / ( 3.0_dp * n_outer )
               write(*,'(a,i6,a,es10.3,a,es9.2,a,es9.2,a,es9.2,a,es9.2,a,i5,a,es9.2)') &
                  '  step=', it, '  t=', t_now, '  mass-imbal=', imbal, &
                  '  du/dt=', du_max, '  max|u|=', umax, '  CFL~=', cfl_max, &
                  '  lin-it=', itl_max, '  lin-res=', gtmp
               write(*,'(a,es9.2,a,es9.2)') '    Tmin=', &
                  minval(fld%T(1:lm%nowned)), '  Tmax=', &
                  maxval(fld%T(1:lm%nowned))
            end if
         end if

         if ( n_out > 0 .and. ( mod(it,n_out) == 0 .or. it == n_steps ) ) then
            call make_step_filename( vtu_prefix, it, step_ext, vtufile )
            call write_snapshot_mpi( lm, fld, part, m_g, c_g, g_g, ctrl, &
                                     trim(vtufile), fld_g, ios )
            if ( ios /= 0 .and. myrank == 0 ) then
               write(*,'(a,i0,a)') 'WARNING: snapshot write failed at step ', &
                                  it, ' (continuing)'
            end if
         end if
      end do

      if ( myrank == 0 ) &
         write(*,'(a,i0,a,es12.4)') 'PIMPLE (MPI) finished: ', n_steps, &
                                     ' time steps, t_end =', t_now

      deallocate( rhs, ap, pp, dpp, phif, imbalc )

   end subroutine pimple_run_mpi

   !----------------------------------------------------------------------------
   ! MPI wrapper for the outflow (fully-developed) global mass scaling: local
   ! partial sums on each rank, allreduce to global totals, then the purely
   ! local scaling application.  See outflow_mass_sums / outflow_mass_scale.
   !----------------------------------------------------------------------------
   subroutine outflow_mass_rescale_mpi( lm, bcs, fld )
      type(local_mesh_t), intent(in)    :: lm
      type(bc_t),         intent(in)    :: bcs
      type(fields_t),     intent(inout) :: fld

      integer  :: ierr
      real(dp) :: m_req, m_out, area_out, beta
      real(dp) :: buf(3)

      call outflow_mass_sums( lm%m, lm%g, bcs, fld, m_req, m_out, area_out )
      buf = [ m_req, m_out, area_out ]
      call mpi_allreduce( MPI_IN_PLACE, buf, 3, MPI_DOUBLE_PRECISION, MPI_SUM, &
                          mpi_comm, ierr )
      call mpi_check( ierr, 'outflow_mass_rescale_mpi: allreduce' )
      m_req = buf(1); m_out = buf(2); area_out = buf(3)
      call outflow_mass_scale( lm%m, lm%g, bcs, fld, m_req, m_out, area_out, &
                               beta )

   end subroutine outflow_mass_rescale_mpi

   !----------------------------------------------------------------------------
   ! PPE assembly with rank-0-only pin. Same as the serial ppe_assembly but
   ! the pin (p'=0 at cell 1) is applied only on rank 0 to avoid
   ! over-constraining the distributed system.
   !----------------------------------------------------------------------------
   subroutine ppe_assembly_mpi( lm, ctrl, bcs, fld, A, rhs )
      type(local_mesh_t), intent(in)    :: lm
      type(ctrl_t),       intent(in)    :: ctrl
      type(bc_t),          intent(in)    :: bcs
      type(fields_t),     intent(in)    :: fld
      type(csr_t),         intent(inout) :: A
      real(dp),            intent(out)   :: rhs(lm%m%ncells)

      integer  :: i, k, kk, c0, c1, kd, gi, ierr
      real(dp) :: dbf, dn, af
      logical  :: has_pdir

      A%val = 0.0_dp
      rhs   = 0.0_dp

      do i = 1, lm%m%nfaces
         c1 = lm%m%f(i)%c1
         c0 = lm%m%f(i)%c0

         if ( c1 == 0 ) then
            ! boundary face: add flux to RHS; Dirichlet-p faces add af to diagonal
            gi = bcs%fgrp(i)
            if ( gi > 0 ) then
               if ( bcs%gb(gi)%btype == BC_POUTLET .or. &
                    bcs%gb(gi)%btype == BC_FARFIELD ) then
                  dbf = lm%g%vol(c0) / fld%apc(c0)
                  dn  = norm2( lm%g%xf(:,i) - lm%g%xc(:,c0) )
                  af  = ctrl%rho * dbf * lm%g%area(i) / dn
                  kd  = pos_of( A, c0, c0 )
                  A%val(kd) = A%val(kd) + af
               end if
            end if
            rhs(c0) = rhs(c0) - fld%flux(i)
            cycle
         end if

         dbf = fld%lf(i)           * lm%g%vol(c0) / fld%apc(c0) &
             + (1.0_dp-fld%lf(i)) * lm%g%vol(c1) / fld%apc(c1)
         dn = norm2( lm%g%xc(:,c1) - lm%g%xc(:,c0) )
         af = ctrl%rho * dbf * lm%g%area(i) / dn

         kd = pos_of( A, c0, c0 ); A%val(kd) = A%val(kd) + af
         kd = pos_of( A, c0, c1 ); A%val(kd) = A%val(kd) - af
         kd = pos_of( A, c1, c1 ); A%val(kd) = A%val(kd) + af
         kd = pos_of( A, c1, c0 ); A%val(kd) = A%val(kd) - af

         rhs(c0) = rhs(c0) - fld%flux(i)
         rhs(c1) = rhs(c1) + fld%flux(i)
      end do

      ! ---- pin cell 1 only on rank 0, all-Neumann domains only ----------------
      ! Same rationale as the serial ppe_assembly: with any pressure-Dirichlet
      ! boundary face (POUTLET / FARFIELD, on ANY rank) the system is already
      ! nonsingular and the pin would discard the pinned cell's continuity row,
      ! turning it into a permanent mass source/sink.  The flag is reduced over
      ! all ranks; the pin (needed only for the closed-domain case) is applied
      ! on rank 0 alone.
      has_pdir = .false.
      do i = 1, lm%m%nfaces
         if ( lm%m%f(i)%c1 == 0 ) then
            gi = bcs%fgrp(i)
            if ( gi > 0 ) then
               if ( bcs%gb(gi)%btype == BC_POUTLET .or. &
                    bcs%gb(gi)%btype == BC_FARFIELD ) has_pdir = .true.
            end if
         end if
      end do
      call mpi_allreduce( MPI_IN_PLACE, has_pdir, 1, MPI_LOGICAL, MPI_LOR, &
                          mpi_comm, ierr )
      call mpi_check( ierr, 'ppe_assembly_mpi: allreduce has_pdir' )

      if ( myrank == 0 .and. .not. has_pdir ) then
         do k = A%row_ptr(1), A%row_ptr(2) - 1
            if ( A%col_idx(k) == 1 ) then
               A%val(k) = 1.0_dp
            else
               A%val(k) = 0.0_dp
            end if
         end do
         do kk = 2, A%nrows
            kd = pos_of( A, kk, 1 )
            if ( kd > 0 ) A%val(kd) = 0.0_dp
         end do
         rhs(1) = 0.0_dp
      end if

   end subroutine ppe_assembly_mpi

   !----------------------------------------------------------------------------
   ! Exchange all cell-based fields: u, p, T, gu, gp, gt
   !----------------------------------------------------------------------------
   subroutine exchange_all_fields( hi, fld, ierr )
      type(halo_info_t), intent(in)    :: hi
      type(fields_t),   intent(inout) :: fld
      integer,          intent(out)   :: ierr
      integer :: d

      call halo_exchange_vector( hi, fld%u, ierr )
      call mpi_check( ierr, 'ex u' )
      call halo_exchange_scalar( hi, fld%p, ierr )
      call mpi_check( ierr, 'ex p' )
      call halo_exchange_scalar( hi, fld%T, ierr )
      call mpi_check( ierr, 'ex T' )
      ! gu is (3,3,ncells): exchange each direction
      do d = 1, 3
         call halo_exchange_vector( hi, fld%gu(d,:,:), ierr )
         call mpi_check( ierr, 'ex gu' )
      end do
      call halo_exchange_vector( hi, fld%gp, ierr )
      call mpi_check( ierr, 'ex gp' )
      call halo_exchange_vector( hi, fld%gt, ierr )
      call mpi_check( ierr, 'ex gt' )
   end subroutine exchange_all_fields

   !----------------------------------------------------------------------------
   ! Exchange only gradient fields: gu, gp, gt
   !----------------------------------------------------------------------------
   subroutine exchange_gradients( hi, fld, ierr )
      type(halo_info_t), intent(in)    :: hi
      type(fields_t),   intent(inout) :: fld
      integer,          intent(out)   :: ierr
      integer :: d

      do d = 1, 3
         call halo_exchange_vector( hi, fld%gu(d,:,:), ierr )
         call mpi_check( ierr, 'ex gu post' )
      end do
      call halo_exchange_vector( hi, fld%gp, ierr )
      call mpi_check( ierr, 'ex gp post' )
      call halo_exchange_vector( hi, fld%gt, ierr )
      call mpi_check( ierr, 'ex gt post' )
   end subroutine exchange_gradients

end module mod_uns_simple_mpi
