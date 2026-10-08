!===============================================================================
! mod_simple.f90 -- Incompressible SIMPLE (steady) and PISO (transient) solvers
!
! Collocated arrangement with Rhie-Chow momentum interpolation.
!
! Steady SIMPLE (simple_run, phase 3), per outer iteration:
!   1. assemble + solve the three momentum components (upwind convection,
!      optionally blended with a limited second-order upwind deferred
!      correction, diffusion with lagged non-orthogonal correction,
!      BiCGSTAB/ILU(0)), with under-relaxation
!   2. Green-Gauss gradients (needed by Rhie-Chow)
!   3. face mass fluxes via Rhie-Chow interpolation
!   4. assemble + solve the pressure-correction (Poisson) equation (CG)
!   5. correct velocity, pressure and fluxes
! Convergence: maximum cell mass imbalance and velocity change.
!
! Transient PISO (piso_run, phase 8), per time step (Issa 1986, collocated):
!   1. predict: solve momentum with an implicit time term added to ap/rhs;
!      no under-relaxation. The time scheme is set by ctrl%time_scheme and the
!      per-step flag fld%ts_order (phase 9):
!        ts_order=1  implicit Euler  rho*V/dt*(u^{n+1} - u^n)
!        ts_order=2  BDF2           (3u^{n+1}-4u^n+u^{n-1})/(2dt)
!      Step 1 bootstraps with Euler (no u^{n-1}); step 2+ uses BDF2 when
!      time_scheme=2. History shift u_old_old<-u_old<-u happens each step.
!   2. gradients + Rhie-Chow fluxes from predicted velocity
!   3. n_correct non-iterative pressure-correction sweeps (default 2):
!      a. assemble + solve PPE for p' from current flux imbalance
!      b. correct u, p, flux (alpha_p = 1; full correction)
!      c. between sweeps re-evaluate gradients + Rhie-Chow fluxes
!   4. advance: t <- t + dt, write VTU per n_out_every steps
!
! Boundary treatment:
!   wall            : no-slip (optionally moving lid), zero-gradient p
!   symmetry        : specular reflection of velocity, zero-gradient p
!   velocity-inlet  : fixed velocity, zero-gradient p
!   pressure-outlet : fixed pressure, extrapolated velocity
!   outflow         : fully-developed -- zero normal gradient for u, p, T;
!                     outflow face fluxes rescaled globally every iteration
!                     (outflow_mass_rescale) so total outflow = total inflow,
!                     keeping the pure-Neumann PPE consistent (pinned cell 1)
!===============================================================================
module mod_uns_simple
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   use mod_uns_connectivity
   use mod_uns_geometry
   use mod_uns_control
   use mod_uns_bc
   use mod_uns_fields
   use mod_uns_linsolver
   use mod_uns_output
   implicit none
   private
   public :: simple_run, piso_run, pimple_run
   public :: build_cell_csr, momentum_assembly, correct_fields, &
             temperature_assembly, solid_temperature_assembly, &
             flux_rhiechow, make_step_filename, pos_of, &
             outflow_mass_sums, outflow_mass_scale, outflow_mass_rescale

contains

   !----------------------------------------------------------------------------
   ! Main SIMPLE driver
   !----------------------------------------------------------------------------
   subroutine simple_run( m, conn, g, ctrl, bcs, fld, ier, nsteps )
      type(mesh_t),    intent(in)    :: m
      type(conn_t),    intent(in)    :: conn
      type(geom_t),    intent(in)    :: g
      type(ctrl_t),    intent(in)    :: ctrl
      type(bc_t),      intent(in)    :: bcs
      type(fields_t),  intent(inout) :: fld
      integer,         intent(out)   :: ier
      integer,         intent(in), optional :: nsteps  ! cap outer iterations
                                                     ! (coupling driver: run
                                                     ! N steps then exchange)

      type(csr_t) :: amat, pmat
      real(dp), allocatable :: rhs(:), uo(:,:), pp(:), dpp(:,:), phif(:)
      real(dp), allocatable :: ap(:), imbalc(:), Tprev(:)
      real(dp) :: uscale, area_btot, flux_ref, imbal, du_max, dT_max
      real(dp) :: resl, ressum
      integer  :: it, comp, itl, ierr, itl_max, i, k, c0, it_max
      integer, parameter :: PRINT_EVERY = 50

      ier = 0

      ! ---- matrix pattern (shared by momentum and pressure equations) --------
      call build_cell_csr( m, conn, amat )
      pmat = amat                        ! same pattern, values set per solve
      allocate( rhs(m%ncells), ap(m%ncells), uo(3,m%ncells), pp(m%ncells) )
      allocate( dpp(3,m%ncells), phif(m%nfaces), imbalc(m%ncells), Tprev(m%ncells) )
      dT_max = huge(1.0_dp)

      ! ---- reference scales for normalization --------------------------------
      uscale = 0.0_dp
      area_btot = 0.0_dp
      do i = 1, bcs%nb
         uscale = max( uscale, norm2( bcs%gb(i)%lid_vel ), &
                               norm2( bcs%gb(i)%uvel ), &
                               bcs%gb(i)%uspeed )
         do k = 1, bcs%gb(i)%nf
            area_btot = area_btot + g%area( bcs%gb(i)%faces(k) )
         end do
      end do
      if ( uscale <= 0.0_dp ) uscale = 1.0_dp
      if ( area_btot <= 0.0_dp ) area_btot = 1.0_dp
      flux_ref = ctrl%rho * uscale * area_btot

      write(*,'(a)') ''
      write(*,'(a)') '--- SIMPLE iteration ---'
      write(*,'(a,es10.3,a,es10.3)') '  reference: uscale =', uscale, &
                                     '   flux_ref =', flux_ref

      fld%apc = 1.0_dp
      fld%ts_order = 0                  ! steady: no time term, use under-relaxation
      call compute_gradients( m, g, bcs, fld )

      ! ---- outer iterations ---------------------------------------------------
      it = 0
      it_max = ctrl%outer_max
      if ( present(nsteps) ) it_max = min( it_max, nsteps )
      do it = 1, it_max
         uo = fld%u
         itl_max = 0
         ressum = 0.0_dp

         ! cold-start ramp: scale the mass-flow-inlet speed AND the
         ! pressure-Dirichlet (outlet) value together (inlet_ramp<=1 disables)
         if ( ctrl%inlet_ramp > 1 ) then
            call set_inlet_ramp_factor( real(it,dp) / real(ctrl%inlet_ramp,dp) )
            call set_pval_ramp_factor( real(it,dp) / real(ctrl%inlet_ramp,dp) )
         end if

         ! 1. momentum equations -------------------------------------------------
         do comp = 1, 3
            call momentum_assembly( m, g, ctrl, bcs, fld, comp, amat, rhs, ap )
            call bicgstab_ilu0( amat, rhs, fld%u(comp,:), ctrl%lin_tol, &
                                ctrl%lin_max, itl, resl, ierr )
            if ( ierr == 2 ) then
               write(*,'(a,i0,a,i0)') 'FATAL: BiCGSTAB breakdown, comp=', &
                                      comp, ', outer it=', it
               ier = 20
               return
            end if
            itl_max = max( itl_max, itl )
            ressum = ressum + resl
         end do
         fld%apc = ap / ctrl%alpha_u    ! relaxed diagonal for Rhie-Chow

         ! 2. gradients (current state) -------------------------------------------
         call compute_gradients( m, g, bcs, fld )

         ! 3. Rhie-Chow face fluxes -------------------------------------------------
         call flux_rhiechow( m, g, ctrl, bcs, fld )
         call outflow_mass_rescale( m, g, bcs, fld )

         ! 4. pressure-correction equation -------------------------------------------
         call ppe_assembly( m, g, ctrl, bcs, fld, pmat, rhs )
         pp = 0.0_dp
         select case ( trim(ctrl%ppe_precond) )
         case ( 'ic0' )
            call cg_ic0( pmat, rhs, pp, ctrl%lin_tol, ctrl%lin_max, &
                         itl, resl, ierr )
         case default
            call cg_jacobi( pmat, rhs, pp, ctrl%lin_tol, ctrl%lin_max, &
                            itl, resl, ierr )
         end select
         if ( ierr == 2 ) then
            write(*,'(a,i0)') 'FATAL: CG breakdown in pressure correction, it=', it
            ier = 21
            return
         end if
         itl_max = max( itl_max, itl )

         ! 5. corrections --------------------------------------------------------------
         call correct_fields( m, g, ctrl, bcs, fld, pp, phif, dpp )

         ! 5b. temperature equation (passive scalar; gradients lagged) -------------
         Tprev = fld%T
         call temperature_assembly( m, g, ctrl, bcs, fld, amat, rhs, ap )
         call bicgstab_ilu0( amat, rhs, fld%T, ctrl%lin_tol, &
                             ctrl%lin_max, itl, resl, ierr )
         if ( ierr == 2 ) then
            write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in temperature, it=', it
            ier = 22; return
         end if
         itl_max = max( itl_max, itl )

         ! 5c. solid temperature equation (LTNE only; no convection) -------------
         if ( ctrl%thermal_model == 'ltne' ) then
            call solid_temperature_assembly( m, g, ctrl, bcs, fld, amat, rhs, ap )
            call bicgstab_ilu0( amat, rhs, fld%T_s, ctrl%lin_tol, &
                                ctrl%lin_max, itl, resl, ierr )
            if ( ierr == 2 ) then
               write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in solid T, it=', it
               ier = 23; return
            end if
            itl_max = max( itl_max, itl )
         end if

         ! 6. convergence measures -------------------------------------------------------
         imbalc = 0.0_dp
         do i = 1, m%nfaces
            c0 = m%f(i)%c0
            imbalc(c0) = imbalc(c0) + fld%flux(i)
            if ( m%f(i)%c1 > 0 ) &
               imbalc(m%f(i)%c1) = imbalc(m%f(i)%c1) - fld%flux(i)
         end do
         imbal = maxval( abs( imbalc ) ) / flux_ref
         du_max = maxval( abs( fld%u - uo ) )
         ! temperature change per sweep: required for buoyancy-coupled runs
         ! where the velocity can stall while T is still evolving
         dT_max = maxval( abs( fld%T - Tprev ) )

         if ( mod(it,PRINT_EVERY) == 0 .or. it == 1 .or. &
              ( imbal < ctrl%outer_tol .and. du_max < ctrl%outer_tol .and. &
                dT_max < ctrl%outer_tol ) ) then
            if ( ctrl%thermal_model == 'ltne' ) then
               write(*,'(a,i6,a,es10.3,a,es10.3,a,es10.3,a,i5,a,es9.2,a,es9.2,a,es9.2,a,es9.2,a,es9.2,a)') &
                  '  it=', it, '  mass-imbal=', imbal, '  du_max=', du_max, &
                  '  dT_max=', dT_max, &
                  '  lin-it=', itl_max, '  lin-res=', ressum/3.0_dp, &
                  '  Tf[', minval(fld%T), ',', maxval(fld%T), &
                  '] Ts[', minval(fld%T_s), ',', maxval(fld%T_s), ']'
            else
               write(*,'(a,i6,a,es10.3,a,es10.3,a,es10.3,a,i5,a,es9.2,a,es9.2,a,es9.2)') &
                  '  it=', it, '  mass-imbal=', imbal, '  du_max=', du_max, &
                  '  dT_max=', dT_max, &
                  '  lin-it=', itl_max, '  lin-res=', ressum/3.0_dp, &
                  '  Tmin=', minval(fld%T), '  Tmax=', maxval(fld%T)
            end if
         end if

         if ( imbal < ctrl%outer_tol .and. du_max < ctrl%outer_tol .and. &
              dT_max < ctrl%outer_tol ) then
            write(*,'(a,i0,a)') 'CONVERGED after ', it, ' outer iterations'
            exit
         end if
      end do

      if ( it >= it_max .and. &
           .not. ( imbal < ctrl%outer_tol .and. du_max < ctrl%outer_tol .and. &
                   dT_max < ctrl%outer_tol ) ) then
         if ( present(nsteps) ) then
            write(*,'(a,i0,a)') 'SIMPLE: stopped after ', it, &
                                 ' outer iterations (nsteps cap)'
         else
            write(*,'(a)') 'WARNING: SIMPLE reached outer_max without converging'
         end if
      end if

      deallocate( rhs, ap, uo, pp, dpp, phif, imbalc, Tprev )

   end subroutine simple_run

   !----------------------------------------------------------------------------
   ! Main PISO driver (transient, phase 8 + BDF2 phase 9)
   !
   ! Per time step:
   !   1. advance history (u_old_old<-u_old<-u) and set ts_order
   !      (step 1 / time_scheme=1 => implicit Euler; step 2+ / time_scheme=2
   !       => BDF2)
   !   2. predict momentum (selected time term, no under-relaxation)
   !   3. gradients + Rhie-Chow fluxes
   !   4. n_correct pressure-correction sweeps (default 2); between sweeps
   !      re-evaluate gradients + Rhie-Chow to keep F consistent with u
   !   5. report mass imbalance, advance time, optional VTU output
   !----------------------------------------------------------------------------
   subroutine piso_run( m, conn, g, ctrl, bcs, fld, vtu_prefix, ier )
      type(mesh_t),    intent(in)    :: m
      type(conn_t),    intent(in)    :: conn
      type(geom_t),    intent(in)    :: g
      type(ctrl_t),    intent(in)    :: ctrl
      type(bc_t),      intent(in)    :: bcs
      type(fields_t),  intent(inout) :: fld
      character(len=*),intent(in)    :: vtu_prefix   ! VTU filename prefix
      integer,         intent(out)   :: ier

      type(csr_t) :: amat, pmat
      real(dp), allocatable :: rhs(:), pp(:), dpp(:,:), phif(:)
      real(dp), allocatable :: ap(:), imbalc(:)
      real(dp) :: uscale, area_btot, flux_ref, imbal, du_max
      real(dp) :: resl, ressum, t_now, umax, cfl_max, dt_dx
      integer  :: it, comp, itl, ierr, itl_max, icorr, i, k, c0
      integer  :: n_steps, n_out, ios
      character(len=512) :: vtufile
      character(len=8)   :: step_ext
      integer, parameter :: PRINT_EVERY = 1

      ier = 0

      ! ---- sanity checks ------------------------------------------------------
      if ( .not. ctrl%transient ) then
         write(*,'(a)') 'FATAL: piso_run invoked with transient=false'
         ier = 30; return
      end if
      if ( ctrl%dt <= 0.0_dp ) then
         write(*,'(a,es12.4)') 'FATAL: dt must be > 0 for transient, got ', &
                               ctrl%dt
         ier = 31; return
      end if
      if ( ctrl%n_time_max <= 0 ) then
         write(*,'(a,i0)') 'FATAL: n_time_max must be > 0, got ', ctrl%n_time_max
         ier = 32; return
      end if
      if ( ctrl%n_correct < 1 ) then
         write(*,'(a,i0)') 'WARNING: n_correct < 1, forcing to 1'
      end if
      n_steps  = max( 1, ctrl%n_time_max )
      n_out    = max( 0, ctrl%n_out_every )
      if ( ctrl%out_format == 'tecplot' ) then
         step_ext = '.plt'
      else
         step_ext = '.vtu'
      end if

      ! ---- matrix pattern (shared by momentum and pressure equations) --------
      call build_cell_csr( m, conn, amat )
      pmat = amat
      allocate( rhs(m%ncells), ap(m%ncells), pp(m%ncells) )
      allocate( dpp(3,m%ncells), phif(m%nfaces), imbalc(m%ncells) )

      ! ---- reference scales for normalization --------------------------------
      uscale = 0.0_dp
      area_btot = 0.0_dp
      do i = 1, bcs%nb
         uscale = max( uscale, norm2( bcs%gb(i)%lid_vel ), &
                               norm2( bcs%gb(i)%uvel ), &
                               bcs%gb(i)%uspeed )
         do k = 1, bcs%gb(i)%nf
            area_btot = area_btot + g%area( bcs%gb(i)%faces(k) )
         end do
      end do
      if ( uscale <= 0.0_dp ) uscale = 1.0_dp
      if ( area_btot <= 0.0_dp ) area_btot = 1.0_dp
      flux_ref = ctrl%rho * uscale * area_btot

      write(*,'(a)') ''
      write(*,'(a)') '--- PISO time stepping ---'
      write(*,'(a,es10.3,a,i0,a,i0)') '  dt=', ctrl%dt, '  n_steps=', &
                                       n_steps, '  n_correct=', ctrl%n_correct
      write(*,'(a,es10.3,a,es10.3)') '  reference: uscale =', uscale, &
                                     '   flux_ref =', flux_ref

      ! ---- initial state: gradients + flux at t=0 ----------------------------
      fld%apc = 1.0_dp
      call compute_gradients( m, g, bcs, fld )
      call flux_rhiechow( m, g, ctrl, bcs, fld )

      t_now = 0.0_dp

      ! ---- time stepping ------------------------------------------------------
      do it = 1, n_steps

         ! 1. advance history + select time scheme ------------------------------
         ! BDF2 needs u^{n} (u_old) and u^{n-1} (u_old_old): shift the history
         ! back BEFORE overwriting u_old with the current solution. Step 1 has
         ! no u^{n-1}, so it bootstraps with implicit Euler (ts_order=1); from
         ! step 2 on (when time_scheme=2) the full BDF2 formula is used.
         fld%u_old_old = fld%u_old
         fld%u_old     = fld%u
         fld%T_old_old = fld%T_old
         fld%T_old     = fld%T
         fld%T_s_old_old = fld%T_s_old
         fld%T_s_old     = fld%T_s
         if ( it == 1 .or. ctrl%time_scheme == 1 ) then
            fld%ts_order = 1   ! implicit Euler (bootstrap or user choice)
         else
            fld%ts_order = 2   ! BDF2 (2nd order)
         end if

         ! 2. predict momentum with the selected time term ----------------------
         itl_max = 0
         ressum = 0.0_dp
         do comp = 1, 3
            call momentum_assembly( m, g, ctrl, bcs, fld, comp, amat, rhs, ap )
            call bicgstab_ilu0( amat, rhs, fld%u(comp,:), ctrl%lin_tol, &
                                ctrl%lin_max, itl, resl, ierr )
            if ( ierr == 2 ) then
               write(*,'(a,i0,a,i0)') 'FATAL: BiCGSTAB breakdown, comp=', &
                                      comp, ', time step=', it
               ier = 33; return
            end if
            itl_max = max( itl_max, itl )
            ressum = ressum + resl
         end do
         fld%apc = ap                      ! un-relaxed diag for Rhie-Chow

         ! 3. gradients + Rhie-Chow fluxes from predicted velocity --------------
         call compute_gradients( m, g, bcs, fld )
         call flux_rhiechow( m, g, ctrl, bcs, fld )
         call outflow_mass_rescale( m, g, bcs, fld )

         ! 4. n_correct pressure-correction sweeps ------------------------------
         ! Standard collocated PISO: each sweep assembles the PPE from the
         ! current fld%flux (which correct_fields updates in-place via the
         ! PISO flux correction F -= A_f*(p'_c1 - p'_c0)). The Rhie-Chow is
         ! NOT re-evaluated between sweeps -- the flux correction is the
         ! consistent PISO update. Gradients of p' used for the velocity
         ! correction are computed inside correct_fields.
         do icorr = 1, ctrl%n_correct
            call ppe_assembly( m, g, ctrl, bcs, fld, pmat, rhs )
            pp = 0.0_dp
            select case ( trim(ctrl%ppe_precond) )
            case ( 'ic0' )
               call cg_ic0( pmat, rhs, pp, ctrl%lin_tol, ctrl%lin_max, &
                            itl, resl, ierr )
            case default
               call cg_jacobi( pmat, rhs, pp, ctrl%lin_tol, ctrl%lin_max, &
                               itl, resl, ierr )
            end select
            if ( ierr == 2 ) then
               write(*,'(a,i0,a,i0)') 'FATAL: CG breakdown in PPE, sweep=', &
                                      icorr, ', time step=', it
               ier = 34; return
            end if
            itl_max = max( itl_max, itl )

            call correct_fields( m, g, ctrl, bcs, fld, pp, phif, dpp )
         end do

         ! 4b. temperature equation (passive scalar; uses corrected flux + lagged gt)
         call temperature_assembly( m, g, ctrl, bcs, fld, amat, rhs, ap )
         call bicgstab_ilu0( amat, rhs, fld%T, ctrl%lin_tol, &
                             ctrl%lin_max, itl, resl, ierr )
         if ( ierr == 2 ) then
            write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in temperature, step=', it
            ier = 35; return
         end if
         itl_max = max( itl_max, itl )

         ! 4c. solid temperature equation (LTNE only)
         if ( ctrl%thermal_model == 'ltne' ) then
            call solid_temperature_assembly( m, g, ctrl, bcs, fld, amat, rhs, ap )
            call bicgstab_ilu0( amat, rhs, fld%T_s, ctrl%lin_tol, &
                                ctrl%lin_max, itl, resl, ierr )
            if ( ierr == 2 ) then
               write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in solid T, step=', it
               ier = 36; return
            end if
            itl_max = max( itl_max, itl )
         end if

         ! 5. advance time -------------------------------------------------------
         t_now = t_now + ctrl%dt

         ! ---- diagnostics: max mass imbalance, |du|/dt, CFL, max|u| ----------
         imbalc = 0.0_dp
         do i = 1, m%nfaces
            c0 = m%f(i)%c0
            imbalc(c0) = imbalc(c0) + fld%flux(i)
            if ( m%f(i)%c1 > 0 ) &
               imbalc(m%f(i)%c1) = imbalc(m%f(i)%c1) - fld%flux(i)
         end do
         imbal = maxval( abs( imbalc ) ) / flux_ref
         du_max = maxval( abs( fld%u - fld%u_old ) ) / ctrl%dt
         umax   = maxval( norm2( fld%u, dim=1 ) )

         ! crude CFL estimate from umax * dt / min cell size (bbox-based)
         cfl_max = 0.0_dp
         do i = 1, m%ncells
            dt_dx = umax * ctrl%dt / ( g%vol(i) ** (1.0_dp/3.0_dp) )
            if ( dt_dx > cfl_max ) cfl_max = dt_dx
         end do

         if ( mod(it-1,PRINT_EVERY) == 0 .or. it == 1 .or. it == n_steps ) then
            write(*,'(a,i6,a,es10.3,a,es9.2,a,es9.2,a,es9.2,a,es9.2,a,i5,a,es9.2)') &
               '  step=', it, '  t=', t_now, '  mass-imbal=', imbal, &
               '  du/dt=', du_max, '  max|u|=', umax, '  CFL~=', cfl_max, &
               '  lin-it=', itl_max, '  lin-res=', ressum/3.0_dp
            write(*,'(a,es9.2,a,es9.2)') '    Tmin=', minval(fld%T), '  Tmax=', &
                                         maxval(fld%T)
         end if

         ! ---- optional result snapshot (VTU or Tecplot PLT) ------------------
         if ( n_out > 0 .and. ( mod(it,n_out) == 0 .or. it == n_steps ) ) then
            call make_step_filename( vtu_prefix, it, step_ext, vtufile )
            if ( ctrl%out_format == 'tecplot' ) then
               call tecplot_write( m, conn, g, ctrl, fld, trim(vtufile), ios )
            else
               call vtk_write( m, conn, g, fld, trim(vtufile), ios )
            end if
            if ( ios /= 0 ) then
               write(*,'(a,i0,a)') 'WARNING: result write failed at step ', &
                                  it, ' (continuing)'
            else
               write(*,'(a,a)') '    wrote ', trim(vtufile)
            end if
         end if
      end do

      write(*,'(a,i0,a,es12.4)') 'PISO finished: ', n_steps, &
                                 ' time steps, t_end =', t_now

      deallocate( rhs, ap, pp, dpp, phif, imbalc )

   end subroutine piso_run

   !----------------------------------------------------------------------------
   ! Main PIMPLE driver (transient, phase 11)
   !
   ! PIMPLE = PISO + SIMPLE-style momentum under-relaxation + optional outer
   ! iterations per time step.  It targets the stability gap of plain PISO:
   ! for large dt (CFL > 1) the lagged convective flux can make the non-
   ! iterative PISO correctors diverge.  PIMPLE recovers stability by:
   !   1. solving the momentum equation with under-relaxation alpha_u
   !      (the time term stays fully implicit, so time accuracy is retained);
   !   2. running n_outer_iter outer sweeps per time step, each sweeping
   !      momentum -> gradients -> Rhie-Chow -> n_correct pressure
   !      corrections (full alpha_p=1, like PISO).
   ! With n_outer_iter=1 and alpha_u<1 this is "PISO with momentum
   ! relaxation", the most common PIMPLE form in production codes.
   !
   ! History (u_old, u_old_old) is advanced once per time step, OUTSIDE the
   ! outer loop, so the time discretisation always references the start-of-
   ! step solution; the outer iterations only converge the nonlinearities
   ! within the step.
   !----------------------------------------------------------------------------
   subroutine pimple_run( m, conn, g, ctrl, bcs, fld, vtu_prefix, ier )
      type(mesh_t),    intent(in)    :: m
      type(conn_t),    intent(in)    :: conn
      type(geom_t),    intent(in)    :: g
      type(ctrl_t),    intent(in)    :: ctrl
      type(bc_t),      intent(in)    :: bcs
      type(fields_t),  intent(inout) :: fld
      character(len=*),intent(in)    :: vtu_prefix
      integer,         intent(out)   :: ier

      type(csr_t) :: amat, pmat
      real(dp), allocatable :: rhs(:), pp(:), dpp(:,:), phif(:)
      real(dp), allocatable :: ap(:), imbalc(:)
      real(dp) :: uscale, area_btot, flux_ref, imbal, du_max
      real(dp) :: resl, ressum, t_now, umax, cfl_max, dt_dx
      integer  :: it, comp, itl, ierr, itl_max, icorr, iouter, i, k, c0
      integer  :: n_steps, n_out, n_outer, ios
      character(len=512) :: vtufile
      character(len=8)   :: step_ext
      integer, parameter :: PRINT_EVERY = 1

      ier = 0

      if ( .not. ctrl%transient ) then
         write(*,'(a)') 'FATAL: pimple_run invoked with transient=false'
         ier = 40; return
      end if
      if ( ctrl%dt <= 0.0_dp ) then
         write(*,'(a,es12.4)') 'FATAL: dt must be > 0 for PIMPLE, got ', ctrl%dt
         ier = 41; return
      end if
      if ( ctrl%n_time_max <= 0 ) then
         write(*,'(a,i0)') 'FATAL: n_time_max must be > 0, got ', ctrl%n_time_max
         ier = 42; return
      end if
      if ( ctrl%n_correct < 1 ) then
         write(*,'(a,i0)') 'WARNING: n_correct < 1, forcing to 1'
      end if
      n_steps = max( 1, ctrl%n_time_max )
      n_outer = max( 1, ctrl%n_outer_iter )
      n_out   = max( 0, ctrl%n_out_every )
      if ( ctrl%out_format == 'tecplot' ) then
         step_ext = '.plt'
      else
         step_ext = '.vtu'
      end if

      call build_cell_csr( m, conn, amat )
      pmat = amat
      allocate( rhs(m%ncells), ap(m%ncells), pp(m%ncells) )
      allocate( dpp(3,m%ncells), phif(m%nfaces), imbalc(m%ncells) )

      uscale = 0.0_dp
      area_btot = 0.0_dp
      do i = 1, bcs%nb
         uscale = max( uscale, norm2( bcs%gb(i)%lid_vel ), &
                               norm2( bcs%gb(i)%uvel ), &
                               bcs%gb(i)%uspeed )
         do k = 1, bcs%gb(i)%nf
            area_btot = area_btot + g%area( bcs%gb(i)%faces(k) )
         end do
      end do
      if ( uscale <= 0.0_dp ) uscale = 1.0_dp
      if ( area_btot <= 0.0_dp ) area_btot = 1.0_dp
      flux_ref = ctrl%rho * uscale * area_btot

      write(*,'(a)') ''
      write(*,'(a)') '--- PIMPLE time stepping ---'
      write(*,'(a,es10.3,a,i0,a,i0,a,i0)') '  dt=', ctrl%dt, '  n_steps=', &
            n_steps, '  n_outer_iter=', n_outer, '  n_correct=', ctrl%n_correct
      write(*,'(a,es10.3,a,es10.3)') '  alpha_u=', ctrl%alpha_u, &
            '  reference flux_ref =', flux_ref

      fld%apc = 1.0_dp
      call compute_gradients( m, g, bcs, fld )
      call flux_rhiechow( m, g, ctrl, bcs, fld )

      t_now = 0.0_dp

      do it = 1, n_steps

         ! advance history once per time step (start-of-step values are fixed
         ! for the time-term; outer iterations only converge the nonlinear
         ! part within the step)
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

         itl_max = 0
         ressum  = 0.0_dp

         ! ---- outer iterations (PIMPLE stabilisation loop) -----------------
         do iouter = 1, n_outer

            ! 1. momentum with under-relaxation (ctrl%pimple=.true. triggers
            !    alpha_u relaxation inside momentum_assembly even though
            !    ts_order > 0)
            do comp = 1, 3
               call momentum_assembly( m, g, ctrl, bcs, fld, comp, amat, rhs, ap )
               call bicgstab_ilu0( amat, rhs, fld%u(comp,:), ctrl%lin_tol, &
                                   ctrl%lin_max, itl, resl, ierr )
               if ( ierr == 2 ) then
                  write(*,'(a,i0,a,i0)') 'FATAL: BiCGSTAB breakdown, comp=', &
                                         comp, ', step=', it
                  ier = 43; return
               end if
               itl_max = max( itl_max, itl )
               ressum = ressum + resl
            end do
            fld%apc = ap / ctrl%alpha_u   ! relaxed diag for Rhie-Chow

            ! 2. gradients + Rhie-Chow fluxes
            call compute_gradients( m, g, bcs, fld )
            call flux_rhiechow( m, g, ctrl, bcs, fld )
            call outflow_mass_rescale( m, g, bcs, fld )

            ! 3. n_correct pressure-correction sweeps (full alpha_p=1)
            do icorr = 1, ctrl%n_correct
               call ppe_assembly( m, g, ctrl, bcs, fld, pmat, rhs )
               pp = 0.0_dp
               select case ( trim(ctrl%ppe_precond) )
               case ( 'ic0' )
                  call cg_ic0( pmat, rhs, pp, ctrl%lin_tol, ctrl%lin_max, &
                               itl, resl, ierr )
               case default
                  call cg_jacobi( pmat, rhs, pp, ctrl%lin_tol, ctrl%lin_max, &
                                  itl, resl, ierr )
               end select
               if ( ierr == 2 ) then
                  write(*,'(a,i0,a,i0)') 'FATAL: CG breakdown in PPE, sweep=', &
                                         icorr, ', step=', it
                  ier = 44; return
               end if
               itl_max = max( itl_max, itl )
               call correct_fields( m, g, ctrl, bcs, fld, pp, phif, dpp )
            end do

         end do   ! outer iterations

         ! 4. temperature equation (passive scalar; solved once per step with
         !    the final corrected velocity field)
         call temperature_assembly( m, g, ctrl, bcs, fld, amat, rhs, ap )
         call bicgstab_ilu0( amat, rhs, fld%T, ctrl%lin_tol, &
                             ctrl%lin_max, itl, resl, ierr )
         if ( ierr == 2 ) then
            write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in T, step=', it
            ier = 45; return
         end if
         itl_max = max( itl_max, itl )

         ! 4b. solid temperature equation (LTNE only)
         if ( ctrl%thermal_model == 'ltne' ) then
            call solid_temperature_assembly( m, g, ctrl, bcs, fld, amat, rhs, ap )
            call bicgstab_ilu0( amat, rhs, fld%T_s, ctrl%lin_tol, &
                                ctrl%lin_max, itl, resl, ierr )
            if ( ierr == 2 ) then
               write(*,'(a,i0)') 'FATAL: BiCGSTAB breakdown in solid T, step=', it
               ier = 46; return
            end if
            itl_max = max( itl_max, itl )
         end if

         ! 5. advance time
         t_now = t_now + ctrl%dt

         ! diagnostics
         imbalc = 0.0_dp
         do i = 1, m%nfaces
            c0 = m%f(i)%c0
            imbalc(c0) = imbalc(c0) + fld%flux(i)
            if ( m%f(i)%c1 > 0 ) &
               imbalc(m%f(i)%c1) = imbalc(m%f(i)%c1) - fld%flux(i)
         end do
         imbal = maxval( abs( imbalc ) ) / flux_ref
         du_max = maxval( abs( fld%u - fld%u_old ) ) / ctrl%dt
         umax   = maxval( norm2( fld%u, dim=1 ) )

         cfl_max = 0.0_dp
         do i = 1, m%ncells
            dt_dx = umax * ctrl%dt / ( g%vol(i) ** (1.0_dp/3.0_dp) )
            if ( dt_dx > cfl_max ) cfl_max = dt_dx
         end do

         if ( mod(it-1,PRINT_EVERY) == 0 .or. it == 1 .or. it == n_steps ) then
            write(*,'(a,i6,a,es10.3,a,es9.2,a,es9.2,a,es9.2,a,es9.2,a,i5,a,es9.2)') &
               '  step=', it, '  t=', t_now, '  mass-imbal=', imbal, &
               '  du/dt=', du_max, '  max|u|=', umax, '  CFL~=', cfl_max, &
               '  lin-it=', itl_max, '  lin-res=', ressum/(3.0_dp*n_outer)
            write(*,'(a,es9.2,a,es9.2)') '    Tmin=', minval(fld%T), '  Tmax=', &
                                         maxval(fld%T)
         end if

         if ( n_out > 0 .and. ( mod(it,n_out) == 0 .or. it == n_steps ) ) then
            call make_step_filename( vtu_prefix, it, step_ext, vtufile )
            if ( ctrl%out_format == 'tecplot' ) then
               call tecplot_write( m, conn, g, ctrl, fld, trim(vtufile), ios )
            else
               call vtk_write( m, conn, g, fld, trim(vtufile), ios )
            end if
            if ( ios /= 0 ) then
               write(*,'(a,i0,a)') 'WARNING: result write failed at step ', &
                                  it, ' (continuing)'
            else
               write(*,'(a,a)') '    wrote ', trim(vtufile)
            end if
         end if
      end do

      write(*,'(a,i0,a,es12.4)') 'PIMPLE finished: ', n_steps, &
                                 ' time steps, t_end =', t_now

      deallocate( rhs, ap, pp, dpp, phif, imbalc )

   end subroutine pimple_run

   !----------------------------------------------------------------------------
   ! Build "<prefix>_stepNNNN<ext>" filename (zero-padded to 4 digits).
   ! ext carries the leading dot, e.g. '.vtu' or '.plt'.
   !----------------------------------------------------------------------------
   subroutine make_step_filename( prefix, step, ext, fname )
      character(len=*), intent(in)  :: prefix
      integer,          intent(in)  :: step
      character(len=*), intent(in)  :: ext
      character(len=*), intent(out) :: fname

      character(len=8) :: num
      integer :: ns, npref, slen, next, i

      write( num, '(i0)' ) step
      ns   = len_trim(num)
      npref = len_trim(prefix)
      next = len_trim(ext)

      ! zero-pad to 4 digits
      slen = max( 4 - ns, 0 )
      fname = ''
      fname(1:npref) = prefix(1:npref)
      fname(npref+1:npref+5) = '_step'
      do i = 1, slen
         fname(npref+5+i:npref+5+i) = '0'
      end do
      fname(npref+5+slen+1:npref+5+slen+ns) = num(1:ns)
      fname(npref+5+slen+ns+1:npref+5+slen+ns+next) = ext(1:next)

   end subroutine make_step_filename

   !----------------------------------------------------------------------------
   ! Cell-based CSR pattern: row = cell, columns = self + unique c2c
   ! neighbours, sorted ascending (mod_linsolver contract)
   !----------------------------------------------------------------------------
   subroutine build_cell_csr( m, conn, A )
      type(mesh_t),   intent(in)  :: m
      type(conn_t),   intent(in)  :: conn
      type(csr_t),    intent(out) :: A

      integer, allocatable :: tmp(:), row(:)
      integer :: k, j, n, cnt, nuniq, pass

      A%nrows = m%ncells
      allocate( A%row_ptr(m%ncells+1) )
      A%row_ptr(1) = 1

      n = 0
      do k = 1, m%ncells
         n = max( n, conn%c2c_ptr(k+1) - conn%c2c_ptr(k) )
      end do
      allocate( tmp(n+1), row(n+1) )

      cnt = 0
      do pass = 1, 2
         cnt = 0
         do k = 1, m%ncells
            ! candidates: neighbours + self
            n = conn%c2c_ptr(k+1) - conn%c2c_ptr(k)
            tmp(1:n+1) = (/ conn%c2c(conn%c2c_ptr(k)+1:conn%c2c_ptr(k+1)), k /)
            call sort_int( n+1, tmp )
            ! unique (already sorted)
            nuniq = 0
            do j = 1, n+1
               if ( nuniq == 0 ) then
                  nuniq = nuniq + 1
                  row(nuniq) = tmp(j)
               else if ( tmp(j) /= row(nuniq) ) then
                  nuniq = nuniq + 1
                  row(nuniq) = tmp(j)
               end if
            end do
            if ( pass == 1 ) then
               A%row_ptr(k+1) = A%row_ptr(k) + nuniq
            else
               do j = 1, nuniq
                  cnt = cnt + 1
                  A%col_idx(cnt) = row(j)
               end do
            end if
         end do
         if ( pass == 1 ) allocate( A%col_idx(A%row_ptr(m%ncells+1)-1), &
                                    A%val(A%row_ptr(m%ncells+1)-1) )
      end do

      deallocate( tmp, row )

   end subroutine build_cell_csr

   ! insertion sort, ascending, in place
   pure subroutine sort_int( n, a )
      integer, intent(in)    :: n
      integer, intent(inout) :: a(n)
      integer :: i, j, key
      do i = 2, n
         key = a(i)
         j = i - 1
         do while ( j >= 1 )
            if ( a(j) <= key ) exit
            a(j+1) = a(j)
            j = j - 1
         end do
         a(j+1) = key
      end do
   end subroutine sort_int

   !----------------------------------------------------------------------------
   ! Position of column j in row i of the CSR matrix (0 if absent)
   !----------------------------------------------------------------------------
   integer function pos_of( A, i, j )
      type(csr_t), intent(in) :: A
      integer,     intent(in) :: i, j
      integer :: k

      pos_of = 0
      do k = A%row_ptr(i), A%row_ptr(i+1) - 1
         if ( A%col_idx(k) == j ) then
            pos_of = k
            return
         else if ( A%col_idx(k) > j ) then
            return
         end if
      end do
   end function pos_of

   !----------------------------------------------------------------------------
   ! Assemble one momentum component into A%val / rhs; ap = unrelaxed diagonal
   !----------------------------------------------------------------------------
   subroutine momentum_assembly( m, g, ctrl, bcs, fld, comp, A, rhs, ap )
      type(mesh_t),   intent(in)    :: m
      type(geom_t),   intent(in)    :: g
      type(ctrl_t),   intent(in)    :: ctrl
      type(bc_t),     intent(in)    :: bcs
      type(fields_t), intent(in)    :: fld
      integer,        intent(in)    :: comp
      type(csr_t),    intent(inout) :: A
      real(dp),       intent(out)   :: rhs(m%ncells)
      real(dp),       intent(out)   :: ap(m%ncells)

      integer  :: i, kk, c0, c1, kd, gi
      real(dp) :: F, D, dn, pf, dw, Scomp, phi_ud, dc, noc, sd, at
      real(dp) :: uf(3), gf(3), dvec(3)
      real(dp) :: mu_eff0, mu_eff1, mu_f, umag
      ! Beavers-Joseph interface-face work arrays
      integer  :: cp, cf
      real(dp) :: D_face, bja, Knn, dPf, Cbj, nh(3), rhs_x
      logical  :: is_bj
      real(dp), allocatable :: psi(:)

      A%val = 0.0_dp
      rhs   = 0.0_dp
      ap    = 0.0_dp
      sd    = 0.0_dp

      ! Barth-Jespersen limiter for the high-order face value (lagged)
      if ( ctrl%conv_blend > 0.0_dp ) then
         allocate( psi(m%ncells) )
         call velocity_limiter( m, g, bcs, fld, comp, psi )
      end if

      ! ---- interior faces ------------------------------------------------------
      do i = 1, m%nfaces
         c1 = m%f(i)%c1
         if ( c1 == 0 ) cycle
         c0 = m%f(i)%c0

         ! Face mass flux (kg/s).  fld%flux already carries rho (see
         ! flux_rhiechow and the temperature assembly comment); multiplying
         ! by rho again here over-weights interior convection by a factor rho
         ! and is inconsistent with the boundary convective coefficient below
         ! (F = rho*(u_f.S_f)), which would leave a spurious O((rho-1)*rho*u^2)
         ! pressure offset stamped by the open boundaries.
         F    = fld%flux(i)
         dvec = g%xc(:,c1) - g%xc(:,c0)
         dn   = norm2( dvec )
         ! effective (Brinkman) viscosity: mu_eff = mu / porosity
         mu_eff0 = ctrl%mu / fld%porosity(c0)
         mu_eff1 = ctrl%mu / fld%porosity(c1)
         mu_f    = fld%lf(i) * mu_eff0 + (1.0_dp - fld%lf(i)) * mu_eff1
         if ( ctrl%nonorth_corr > 0.0_dp ) then
            ! over-relaxed decomposition: implicit normal part
            ! D = mu_eff*|Sf|^2/(Sf.d) >= mu_eff*|Sf|/|d| (Sf.d > 0 for convex cells)
            sd = dot_product( g%sf(:,i), dvec )
            D  = mu_f * g%area(i)**2 / sd
         else
            D  = mu_f * g%area(i) / dn
         end if

         ! ---- Beavers-Joseph slip at fluid/porous interface faces ------------
         ! Face with exactly one porous neighbour and bj_alpha > 0: the BJ
         ! stress jump  du_t/dn|fluid = (alpha/sqrt(K_nn)) (u_face - u_D)
         ! is imposed as a face Robin condition.  Eliminating u_face between
         ! the one-sided fluid gradient (u_face-u_Pf)/d_Pf and the jump
         ! condition gives the effective tangential diffusion coefficient
         !   C = mu*|Sf| * (alpha/sqrt(K_nn)) / (1 + alpha*d_Pf/sqrt(K_nn))
         ! with the wall-at-face limit C -> mu*|Sf|/d_Pf as K -> 0 and free
         ! slip as K -> inf.  The face-normal component keeps the standard
         ! two-point D:  D_face = C + (D-C)*n_comp^2 (exact for axis-aligned
         ! faces); the (D-C)*n_comp*n_k cross part is deferred to the rhs
         ! (lagged, equal and opposite on both cells -> momentum conserving).
         ! Default bj_alpha = 0 keeps this branch bit-identically inactive.
         is_bj  = .false.
         D_face = D
         if ( m%cztype(c0) /= m%cztype(c1) ) then
            if ( m%cztype(c0) == CZ_POROUS ) then
               cp = c0 ; cf = c1
            else
               cp = c1 ; cf = c0
            end if
            bja = fld%bj_alpha(cp)
            if ( bja > 0.0_dp ) then
               nh  = g%sf(:,i) / g%area(i)
               Knn = dot_product( fld%perm_dir(:,cp), nh**2 )
               if ( Knn > 0.0_dp ) then
                  dPf = norm2( g%xf(:,i) - g%xc(:,cf) )
                  Cbj = ctrl%mu * g%area(i) * ( bja / sqrt(Knn) ) &
                        / ( 1.0_dp + bja * dPf / sqrt(Knn) )
                  D_face = Cbj + ( D - Cbj ) * nh(comp)**2
                  rhs_x  = ( D - Cbj ) * nh(comp) * ( &
                           ( fld%u(1,c0) - fld%u(1,c1) ) * nh(1) &
                         + ( fld%u(2,c0) - fld%u(2,c1) ) * nh(2) &
                         + ( fld%u(3,c0) - fld%u(3,c1) ) * nh(3) &
                         - ( fld%u(comp,c0) - fld%u(comp,c1) ) * nh(comp) )
                  rhs(c0) = rhs(c0) + rhs_x
                  rhs(c1) = rhs(c1) - rhs_x
                  is_bj   = .true.
               end if
            end if
         end if

         ! upwind convection + implicit diffusion
         ap(c0) = ap(c0) + max( F, 0.0_dp ) + D_face
         ap(c1) = ap(c1) + max( -F, 0.0_dp ) + D_face
         kd = pos_of( A, c0, c1 )
         A%val(kd) = A%val(kd) - ( D_face + max( -F, 0.0_dp ) )
         kd = pos_of( A, c1, c0 )
         A%val(kd) = A%val(kd) - ( D_face + max( F, 0.0_dp ) )

         ! deferred correction towards a limited second-order upwind face
         ! value phi_f = phi_U + psi_U * grad(phi)_U . (x_f - x_U): the
         ! matrix keeps the (diagonally dominant) first-order upwind part,
         ! while blend * F * (phi_ho - phi_ud) goes explicitly into the rhs,
         ! lagged with the previous outer iterate held in fld%u / fld%gu
         if ( ctrl%conv_blend > 0.0_dp ) then
            if ( F > 0.0_dp ) then
               phi_ud = fld%u(comp,c0)
               dc = psi(c0) * dot_product( fld%gu(comp,:,c0), &
                                           g%xf(:,i) - g%xc(:,c0) )
            else
               phi_ud = fld%u(comp,c1)
               dc = psi(c1) * dot_product( fld%gu(comp,:,c1), &
                                           g%xf(:,i) - g%xc(:,c1) )
            end if
            dc = ctrl%conv_blend * F * dc
            rhs(c0) = rhs(c0) - dc
            rhs(c1) = rhs(c1) + dc
         end if

         ! non-orthogonal diffusion correction (deferred, lagged): with the
         ! over-relaxed split the matrix keeps the two-point normal gradient
         ! mu*|Sf|^2/(Sf.d)*(phi_c1-phi_c0); only the tangential cross term
         ! mu*grad(phi)_f.(Sf - |Sf|^2 d/(Sf.d)) goes explicitly into the
         ! rhs, with the face gradient interpolated from the cell gradients
         if ( ctrl%nonorth_corr > 0.0_dp .and. .not. is_bj ) then
            gf  = fld%lf(i) * fld%gu(comp,:,c0) &
                + ( 1.0_dp - fld%lf(i) ) * fld%gu(comp,:,c1)
            noc = ctrl%nonorth_corr * mu_f * dot_product( gf, &
                     g%sf(:,i) - dvec * ( g%area(i)**2 / sd ) )
            rhs(c0) = rhs(c0) + noc
            rhs(c1) = rhs(c1) - noc
         end if

         ! pressure force -p_f S_f.
         ! At fluid/porous interface faces the plain distance-weighted
         ! interpolation misses the (kinked) interface pressure by
         ! O((s_por-s_flu)*d0*d1/(d0+d1)), which acts on the two adjacent cells
         ! as an equal-and-opposite force dipole and drives the odd-even
         ! velocity jitter seen next to the interface.  Use the kink-consistent
         ! face pressure there (see kink_face_pressure); it collapses to the
         ! plain interpolation for a locally linear field, so single-phase
         ! cases are untouched.
         if ( is_porous_cell(fld,c0) .neqv. is_porous_cell(fld,c1) ) then
            call kink_face_pressure( m, g, fld, i, pf )
         else
            pf = fld%lf(i) * fld%p(c0) + (1.0_dp - fld%lf(i)) * fld%p(c1)
         end if
         Scomp = g%sf(comp,i)
         rhs(c0) = rhs(c0) - pf * Scomp
         rhs(c1) = rhs(c1) + pf * Scomp
      end do

      ! ---- boundary faces --------------------------------------------------------
      do i = 1, m%nfaces
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         c0 = m%f(i)%c0
         call bc_face_vel( bcs, i, g%xf(:,i), g%sf(:,i), fld%u(:,c0), uf )
         call bc_face_p( bcs, i, fld%p(c0), pf )
         dw = norm2( g%xf(:,i) - g%xc(:,c0) )
         D  = ( ctrl%mu / fld%porosity(c0) ) * g%area(i) / dw
         F  = ctrl%rho * dot_product( uf, g%sf(:,i) )

         select case ( bcs%gb(gi)%btype )
         case ( BC_WALL, BC_SYMMETRY, BC_VINLET, BC_VINLET_PARAB, BC_SLIPWALL, &
                BC_MASSINLET, BC_INTERFACE )
            ! fixed face velocity: implicit diffusion (+outflow safety);
            ! (VINLET_PARAB: same Dirichlet treatment, face value = parabolic
            !  profile along the inward normal)
            ! (MASSINLET: uf = -(mdot/rho) n, prescribed normal inflow)
            ! (INTERFACE: velocity is Dirichlet-set by the coupling peer;
            !  without this branch the prescribed velocity never enters the
            !  momentum matrix -- the first cell layer sees only the PPE
            !  flux and settles to ~half the face velocity, producing a
            !  spurious factor-of-two across the coupling interface)
            ! NB: the diagonal is written from ap() in the relaxation loop
            ! below, so boundary contributions must go into ap(c0)
            ap(c0) = ap(c0) + D + max( F, 0.0_dp )
            rhs(c0) = rhs(c0) + D*uf(comp) - F*uf(comp)
         case ( BC_POUTLET, BC_OUTFLOW )
            ! extrapolated face velocity u_f = u_P (first order), upwind split:
            ! outflow (F>0) is implicit; reversed inflow (F<0) must stay
            ! explicit -- adding F to the diagonal then would WEAKEN it and
            ! can drive the momentum (and hence PPE) matrix indefinite
            ap(c0)  = ap(c0) + max( F, 0.0_dp )
            rhs(c0) = rhs(c0) - min( F, 0.0_dp ) * fld%u(comp,c0)
         case ( BC_FARFIELD )
            ! Inflow (F<0, flow into domain): Dirichlet velocity = u_far.
            ! Outflow (F>0): zeroth-order extrapolation (like POUTLET).
            if ( F < 0.0_dp ) then
               ap(c0)  = ap(c0)  + D + max( F, 0.0_dp )
               rhs(c0) = rhs(c0) + D*uf(comp) - F*uf(comp)
            else
               ap(c0)  = ap(c0)  + F
            end if
         end select
         rhs(c0) = rhs(c0) - pf * g%sf(comp,i)
      end do

      ! ---- Boussinesq buoyancy source (phase 10) ------------------------------
      ! Buoyancy force per volume = -rho0*beta*(T-Tref)*g_vec; added to the
      ! momentum rhs. gravity points toward earth; T>Tref gives an upward
      ! (anti-gravity) force. Only active when ctrl%boussinesq = .true.
      if ( ctrl%boussinesq ) then
         do kk = 1, m%ncells
            rhs(kk) = rhs(kk) - ctrl%rho * ctrl%beta &
                      * ( fld%T(kk) - ctrl%tref ) * ctrl%gravity(comp) &
                      * g%vol(kk)
         end do
      end if

      ! ---- uniform body force ------------------------------------------------
      ! Constant force per unit volume f (N/m^3), e.g. a mean pressure-gradient
      ! surrogate in open channels (in a closed box it would be cancelled
      ! exactly by the pressure gradient).  rhs += f_comp * V.
      if ( any( ctrl%body_force /= 0.0_dp ) ) then
         do kk = 1, m%ncells
            rhs(kk) = rhs(kk) + ctrl%body_force(comp) * g%vol(kk)
         end do
      end if

      ! ---- Darcy-Forchheimer porous source (phase 12) -------------------------
      ! Momentum sink per volume:  S = -(mu/K) u - rho*inertial*|u| u
      !   - linear Darcy term -> implicit (added to diagonal ap)
      !   - nonlinear Forchheimer term -> explicit (lagged velocity, into rhs)
      ! K is the diagonal permeability tensor component along the current
      ! momentum direction (axis-aligned anisotropy; off-diagonal terms would
      ! need a block-coupled momentum assembly and are not supported).
      ! Fluid cells have perm_dir=0 and inertial=0, so they are skipped.
      ! NB the Darcy coefficient carries NO porosity: this is the
      ! "divided-by-porosity" form of the volume-averaged momentum equation
      !   rho/eps d(u)/dt + rho/eps^2 div(u u) = -grad p + mu/eps lap(u)
      !                                           - mu/K u - rho cE/sqrt(K)|u|u
      ! i.e. the sink acts on the SUPERFICIAL (seepage) velocity with mu/K
      ! (Brinkman/Vafai-Kim convention, mu_eff = mu/eps for the viscous term
      ! above and in the faces).  Scaling it by 1/eps as well multiplies the
      ! plug pressure drop by 1/eps (1.43 here) and is inconsistent with the
      ! Betchen 2006 reference solution (dp/dx = (mu/K) U/Da-form).
      do kk = 1, m%ncells
         if ( fld%perm_dir(comp,kk) > 0.0_dp ) then
            ap(kk) = ap(kk) + ctrl%mu / fld%perm_dir(comp,kk) * g%vol(kk)
         end if
         if ( fld%inertial(kk) > 0.0_dp ) then
            umag = sqrt( fld%u(1,kk)**2 + fld%u(2,kk)**2 + fld%u(3,kk)**2 )
            rhs(kk) = rhs(kk) - ctrl%rho * fld%inertial(kk) * umag &
                      * fld%u(comp,kk) * g%vol(kk)
         end if
      end do

      ! ---- transient time term ------------------------------------------------
      ! ts_order = 1: implicit Euler    rho*V/dt*(u^{n+1} - u^n)
      !            => at = rho*V/dt;   ap += at;            rhs += at*u_old
      ! ts_order = 2: BDF2  (3u^{n+1} - 4u^n + u^{n-1})/(2dt) = RHS(u^{n+1})
      !            => at = rho*V/(2dt); ap += 3*at;
      !               rhs += at*(4*u_old - u_old_old)
      ! ts_order = 0: steady (no time term; steady under-relaxation below)
      ! Added BEFORE the diagonal under-relaxation loop so that, in transient
      ! mode (alpha_u forced to 1 internally), the time term stays fully
      ! implicit and provides the diagonal dominance that the under-relaxation
      ! would otherwise supply.
      select case ( fld%ts_order )
      case ( 1 )   ! implicit Euler
         do kk = 1, m%ncells
            at = ctrl%rho * g%vol(kk) / ctrl%dt
            ap(kk)  = ap(kk)  + at
            rhs(kk) = rhs(kk) + at * fld%u_old(comp,kk)
         end do
      case ( 2 )   ! BDF2 (2nd order)
         do kk = 1, m%ncells
            at = ctrl%rho * g%vol(kk) / ( 2.0_dp * ctrl%dt )
            ap(kk)  = ap(kk)  + 3.0_dp * at
            rhs(kk) = rhs(kk) + at * ( 4.0_dp * fld%u_old(comp,kk) &
                                       - fld%u_old_old(comp,kk) )
         end do
      end select

      ! ---- diagonal + under-relaxation --------------------------------------------
      ! Three modes:
      !   * steady SIMPLE (ts_order=0):  alpha_u under-relaxation, no time term.
      !   * transient PISO (ts_order>0, pimple=.false.): no relaxation; the
      !     time term supplies diagonal dominance.
      !   * transient PIMPLE (ts_order>0, pimple=.true.): alpha_u under-
      !     relaxation ON TOP of the time term -- the time term keeps the
      !     scheme time-accurate while alpha_u damps the explicit/lagged
      !     convective instability for large dt.
      if ( fld%ts_order > 0 .and. .not. ctrl%pimple ) then
         do kk = 1, m%ncells
            kd = pos_of( A, kk, kk )
            A%val(kd) = ap(kk)
         end do
      else
         do kk = 1, m%ncells
            kd = pos_of( A, kk, kk )
            A%val(kd) = ap(kk) / ctrl%alpha_u
            rhs(kk) = rhs(kk) &
                    + ( 1.0_dp - ctrl%alpha_u ) / ctrl%alpha_u * ap(kk) * fld%u(comp,kk)
         end do
      end if

      if ( allocated(psi) ) deallocate( psi )

   end subroutine momentum_assembly

   !----------------------------------------------------------------------------
   ! Barth-Jespersen slope limiter for one velocity component:
   !   psi(c) = min over faces of the factor keeping the reconstructed face
   !   value phi_c + grad(c).(x_f - x_c) within [phi_min, phi_max] of the
   !   cell neighbourhood (face-neighbours + boundary face values)
   !----------------------------------------------------------------------------
   subroutine velocity_limiter( m, g, bcs, fld, comp, psi )
      type(mesh_t),   intent(in)  :: m
      type(geom_t),   intent(in)  :: g
      type(bc_t),     intent(in)  :: bcs
      type(fields_t), intent(in)  :: fld
      integer,        intent(in)  :: comp
      real(dp),       intent(out) :: psi(m%ncells)

      integer  :: i, c0, c1, gi
      real(dp) :: pmin(m%ncells), pmax(m%ncells)
      real(dp) :: uc, fr, df, uf(3)

      ! neighbourhood extrema (interior neighbours + boundary face values)
      do i = 1, m%ncells
         pmin(i) = fld%u(comp,i)
         pmax(i) = fld%u(comp,i)
      end do
      do i = 1, m%nfaces
         c0 = m%f(i)%c0
         c1 = m%f(i)%c1
         if ( c1 > 0 ) then
            pmin(c0) = min( pmin(c0), fld%u(comp,c1) )
            pmax(c0) = max( pmax(c0), fld%u(comp,c1) )
            pmin(c1) = min( pmin(c1), fld%u(comp,c0) )
            pmax(c1) = max( pmax(c1), fld%u(comp,c0) )
         else
            gi = bcs%fgrp(i)
            if ( gi > 0 ) then
               call bc_face_vel( bcs, i, g%xf(:,i), g%sf(:,i), &
                                 fld%u(:,c0), uf )
               pmin(c0) = min( pmin(c0), uf(comp) )
               pmax(c0) = max( pmax(c0), uf(comp) )
            end if
         end if
      end do

      ! limiter factor from the reconstructed value at every face of the cell
      psi = 1.0_dp
      do i = 1, m%nfaces
         c0 = m%f(i)%c0
         uc = fld%u(comp,c0)
         fr = uc + dot_product( fld%gu(comp,:,c0), g%xf(:,i) - g%xc(:,c0) )
         df = fr - uc
         if ( df > 0.0_dp ) then
            psi(c0) = min( psi(c0), ( pmax(c0) - uc ) / df )
         else if ( df < 0.0_dp ) then
            psi(c0) = min( psi(c0), ( pmin(c0) - uc ) / df )
         end if
         c1 = m%f(i)%c1
         if ( c1 > 0 ) then
            uc = fld%u(comp,c1)
            fr = uc + dot_product( fld%gu(comp,:,c1), g%xf(:,i) - g%xc(:,c1) )
            df = fr - uc
            if ( df > 0.0_dp ) then
               psi(c1) = min( psi(c1), ( pmax(c1) - uc ) / df )
            else if ( df < 0.0_dp ) then
               psi(c1) = min( psi(c1), ( pmin(c1) - uc ) / df )
            end if
         end if
      end do
      psi = max( 0.0_dp, min( 1.0_dp, psi ) )

   end subroutine velocity_limiter

   !----------------------------------------------------------------------------
   ! Face mass fluxes via Rhie-Chow momentum interpolation
   !----------------------------------------------------------------------------
   subroutine flux_rhiechow( m, g, ctrl, bcs, fld )
      type(mesh_t),   intent(in)    :: m
      type(geom_t),   intent(in)    :: g
      type(ctrl_t),   intent(in)    :: ctrl
      type(bc_t),     intent(in)    :: bcs
      type(fields_t), intent(inout) :: fld

      integer  :: i, c0, c1, gi
      real(dp) :: dbf, dn, pfgrad_n
      real(dp) :: ufl(3), gpl(3), dvec(3), uf(3)

      do i = 1, m%nfaces
         c1 = m%f(i)%c1
         if ( c1 == 0 ) then
            gi = bcs%fgrp(i)
            c0 = m%f(i)%c0
            call bc_face_vel( bcs, i, g%xf(:,i), g%sf(:,i), fld%u(:,c0), uf )
            fld%flux(i) = ctrl%rho * dot_product( uf, g%sf(:,i) )
            cycle
         end if
         c0 = m%f(i)%c0

         ufl = fld%lf(i) * fld%u(:,c0) + (1.0_dp - fld%lf(i)) * fld%u(:,c1)
         gpl = fld%lf(i) * fld%gp(:,c0) + (1.0_dp - fld%lf(i)) * fld%gp(:,c1)
         dvec = g%xc(:,c1) - g%xc(:,c0)
         dn = norm2( dvec )
         ! compact (orthogonal) face pressure gradient along d
         pfgrad_n = ( fld%p(c1) - fld%p(c0) ) / dn

         dbf = fld%lf(i)         * g%vol(c0) / fld%apc(c0) &
             + (1.0_dp-fld%lf(i)) * g%vol(c1) / fld%apc(c1)

         fld%flux(i) = ctrl%rho * ( dot_product( ufl, g%sf(:,i) ) &
            + dbf * ( dot_product( gpl, g%sf(:,i) ) &
                      - pfgrad_n * dot_product( dvec, g%sf(:,i) ) / dn ) )
      end do

   end subroutine flux_rhiechow

   !----------------------------------------------------------------------------
   ! Partial sums for the outflow (fully-developed) global mass scaling.
   !   m_req    = total mass flux that MUST leave through the outflow faces,
   !              i.e. -(sum of flux over fixed-flux boundaries: inlets, walls).
   !              Pressure-Dirichlet faces (POUTLET / FARFIELD) are excluded:
   !              they self-adjust through their p' flux correction.
   !   m_out    = current total flux through the outflow faces
   !   area_out = total outflow face area
   ! Purely local: on MPI ranks call this on the local mesh and allreduce the
   ! three sums before outflow_mass_scale (see outflow_mass_rescale_mpi).
   !----------------------------------------------------------------------------
   subroutine outflow_mass_sums( m, g, bcs, fld, m_req, m_out, area_out )
      type(mesh_t),   intent(in)  :: m
      type(geom_t),   intent(in)  :: g
      type(bc_t),     intent(in)  :: bcs
      type(fields_t), intent(in)  :: fld
      real(dp),       intent(out) :: m_req, m_out, area_out

      integer :: i, gi, bt

      m_req = 0.0_dp; m_out = 0.0_dp; area_out = 0.0_dp
      do i = 1, m%nfaces
         if ( m%f(i)%c1 > 0 ) cycle          ! interior face
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         bt = bcs%gb(gi)%btype
         if ( bt == BC_OUTFLOW ) then
            m_out    = m_out + fld%flux(i)
            area_out = area_out + g%area(i)
         else if ( bt /= BC_POUTLET .and. bt /= BC_FARFIELD ) then
            m_req = m_req - fld%flux(i)      ! inflow flux is negative (Sf out)
         end if
      end do

   end subroutine outflow_mass_sums

   !----------------------------------------------------------------------------
   ! Apply the outflow global mass scaling (purely local, given GLOBAL sums).
   ! Normal case: multiply every outflow face flux by beta = m_req/m_out so
   ! the total outflow exactly matches the required inflow -- this makes the
   ! pure-Neumann PPE consistent (sum rhs = 0) and lets the pinned cell 1 be
   ! harmless instead of acting as a fake mass sink.
   ! Startup fallback: with no/negligible outflow yet (m_out ~ 0 after a
   ! from-rest momentum solve) the scaled shape would be amplified noise;
   ! distribute the required flux uniformly by face area instead.
   !----------------------------------------------------------------------------
   subroutine outflow_mass_scale( m, g, bcs, fld, m_req, m_out, area_out, beta )
      type(mesh_t),   intent(in)    :: m
      type(geom_t),   intent(in)    :: g
      type(bc_t),     intent(in)    :: bcs
      type(fields_t), intent(inout) :: fld
      real(dp),       intent(in)    :: m_req, m_out, area_out
      real(dp),       intent(out)   :: beta

      integer :: i, gi

      beta = 1.0_dp
      if ( area_out <= 0.0_dp .or. m_req <= 0.0_dp ) return   ! no outflow work

      if ( m_out > 1.0e-2_dp * m_req ) then
         beta = m_req / m_out
         do i = 1, m%nfaces
            if ( m%f(i)%c1 > 0 ) cycle
            gi = bcs%fgrp(i)
            if ( gi == 0 ) cycle
            if ( bcs%gb(gi)%btype == BC_OUTFLOW ) &
               fld%flux(i) = beta * fld%flux(i)
         end do
      else
         beta = 0.0_dp   ! signals the uniform fallback (diagnostics)
         do i = 1, m%nfaces
            if ( m%f(i)%c1 > 0 ) cycle
            gi = bcs%fgrp(i)
            if ( gi == 0 ) cycle
            if ( bcs%gb(gi)%btype == BC_OUTFLOW ) &
               fld%flux(i) = m_req * g%area(i) / area_out
         end do
      end if

   end subroutine outflow_mass_scale

   !----------------------------------------------------------------------------
   ! Serial convenience wrapper: local sums ARE the global sums.
   !----------------------------------------------------------------------------
   subroutine outflow_mass_rescale( m, g, bcs, fld )
      type(mesh_t),   intent(in)    :: m
      type(geom_t),   intent(in)    :: g
      type(bc_t),     intent(in)    :: bcs
      type(fields_t), intent(inout) :: fld

      real(dp) :: m_req, m_out, area_out, beta

      call outflow_mass_sums( m, g, bcs, fld, m_req, m_out, area_out )
      call outflow_mass_scale( m, g, bcs, fld, m_req, m_out, area_out, beta )

   end subroutine outflow_mass_rescale

   !----------------------------------------------------------------------------
   ! Assemble the pressure-correction equation (pure Neumann + pinned cell 1)
   !----------------------------------------------------------------------------
   subroutine ppe_assembly( m, g, ctrl, bcs, fld, A, rhs )
      type(mesh_t),   intent(in)    :: m
      type(geom_t),   intent(in)    :: g
      type(ctrl_t),   intent(in)    :: ctrl
      type(bc_t),     intent(in)    :: bcs
      type(fields_t), intent(in)    :: fld
      type(csr_t),    intent(inout) :: A
      real(dp),       intent(out)   :: rhs(m%ncells)

      integer  :: i, k, kk, c0, c1, kd, gi
      real(dp) :: dbf, dn, af
      logical  :: has_pdir

      A%val = 0.0_dp
      rhs   = 0.0_dp

      do i = 1, m%nfaces
         c1 = m%f(i)%c1
         c0 = m%f(i)%c0

         if ( c1 == 0 ) then
            ! ---- boundary face -------------------------------------------------
            ! The boundary mass flux must enter the PPE RHS, otherwise the
            ! pressure correction cannot balance inflow/outflow through open
            ! boundaries (far-field, outlet). For Dirichlet-pressure BCs
            ! (POUTLET, FARFIELD) p'=0 on the boundary, so the face acts like
            ! a neighbour with p'=0: add af to the owner diagonal and the
            ! boundary flux to the RHS. For zero-gradient pressure BCs
            ! (WALL, SYMMETRY, VINLET) dp'/dn=0, so no matrix contribution;
            ! the boundary flux still goes to the RHS so the PPE sees the
            ! full cell mass imbalance.
            gi = bcs%fgrp(i)
            if ( gi > 0 ) then
               if ( bcs%gb(gi)%btype == BC_POUTLET .or. &
                    bcs%gb(gi)%btype == BC_FARFIELD ) then
                  dbf = g%vol(c0) / fld%apc(c0)
                  dn  = norm2( g%xf(:,i) - g%xc(:,c0) )
                  af  = ctrl%rho * dbf * g%area(i) / dn
                  kd  = pos_of( A, c0, c0 )
                  A%val(kd) = A%val(kd) + af
               end if
            end if
            rhs(c0) = rhs(c0) - fld%flux(i)
            cycle
         end if

         ! ---- interior face ---------------------------------------------------
         dbf = fld%lf(i)         * g%vol(c0) / fld%apc(c0) &
             + (1.0_dp-fld%lf(i)) * g%vol(c1) / fld%apc(c1)
         dn = norm2( g%xc(:,c1) - g%xc(:,c0) )
         af = ctrl%rho * dbf * g%area(i) / dn

         kd = pos_of( A, c0, c0 ); A%val(kd) = A%val(kd) + af
         kd = pos_of( A, c0, c1 ); A%val(kd) = A%val(kd) - af
         kd = pos_of( A, c1, c1 ); A%val(kd) = A%val(kd) + af
         kd = pos_of( A, c1, c0 ); A%val(kd) = A%val(kd) - af

         rhs(c0) = rhs(c0) - fld%flux(i)
         rhs(c1) = rhs(c1) + fld%flux(i)
      end do

      ! ---- pin cell 1 (all-Neumann / closed domains only) ----------------------
      ! With at least one pressure-Dirichlet boundary face (POUTLET / FARFIELD)
      ! the matrix is already nonsingular and pinning cell 1 OVER-CONSTRAINS the
      ! system: the pinned cell's continuity row is discarded, so it acts as a
      ! permanent mass source/sink, freezing a spurious mass imbalance and
      ! corrupting the flow (and the energy field) around it.  Pin only when no
      ! pressure-Dirichlet face exists anywhere (pure Neumann problem).
      has_pdir = .false.
      do i = 1, m%nfaces
         if ( m%f(i)%c1 == 0 ) then
            gi = bcs%fgrp(i)
            if ( gi > 0 ) then
               if ( bcs%gb(gi)%btype == BC_POUTLET .or. &
                    bcs%gb(gi)%btype == BC_FARFIELD ) has_pdir = .true.
            end if
         end if
      end do

      if ( .not. has_pdir ) then
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

   end subroutine ppe_assembly

   !----------------------------------------------------------------------------
   ! Apply velocity / pressure / flux corrections from p'
   !----------------------------------------------------------------------------
   subroutine correct_fields( m, g, ctrl, bcs, fld, pp, phif, dpp )
      type(mesh_t),   intent(in)    :: m
      type(geom_t),   intent(in)    :: g
      type(ctrl_t),   intent(in)    :: ctrl
      type(bc_t),     intent(in)    :: bcs
      type(fields_t), intent(inout) :: fld
      real(dp),       intent(in)    :: pp(m%ncells)
      real(dp),       intent(inout) :: phif(m%nfaces)   ! work array
      real(dp),       intent(inout) :: dpp(3,m%ncells)   ! work array

      integer  :: i, c0, c1, gi
      real(dp) :: dbf, dn, af

      ! gradient of p'. Interior faces: linear interpolation. Boundary
      ! faces: zero-gradient extrapolation (owner value) EXCEPT at fixed-
      ! pressure boundaries (POUTLET, FARFIELD), where the correction is a
      ! homogeneous Dirichlet condition p'_f = 0 -- treating those faces as
      ! zero-gradient underestimates dpp in the owner cell and leaves the
      ! boundary column under-corrected (observed as a spurious multi-cell
      ! exit adjustment in developed channel flow).
      do i = 1, m%nfaces
         c1 = m%f(i)%c1
         if ( c1 > 0 ) then
            phif(i) = fld%lf(i) * pp(m%f(i)%c0) &
                    + (1.0_dp - fld%lf(i)) * pp(c1)
         else
            gi = bcs%fgrp(i)
            if ( gi > 0 ) then
               if ( bcs%gb(gi)%btype == BC_POUTLET .or. &
                    bcs%gb(gi)%btype == BC_FARFIELD ) then
                  phif(i) = 0.0_dp
                  cycle
               end if
            end if
            phif(i) = pp(m%f(i)%c0)
         end if
      end do
      call grad_scalar( m, g, phif, dpp )

      ! velocity: u += -(V/a_P) grad(p')
      do i = 1, m%ncells
         fld%u(:,i) = fld%u(:,i) - ( g%vol(i) / fld%apc(i) ) * dpp(:,i)
      end do

      ! pressure: in transient PISO the correction is applied in full
      ! (alpha_p = 1) since each corrector step is non-iterative; the time
      ! term supplies the stability that steady SIMPLE draws from under-
      ! relaxation. In steady SIMPLE the standard alpha_p is applied.
      if ( fld%ts_order > 0 ) then
         fld%p = fld%p + pp
      else
         fld%p = fld%p + ctrl%alpha_p * pp
      end if

      ! flux
      do i = 1, m%nfaces
         c1 = m%f(i)%c1
         if ( c1 == 0 ) then
            ! Boundary face. At a fixed-pressure boundary (POUTLET,
            ! FARFIELD) p'_f = 0, so -- exactly as on the PPE row -- the
            ! face mass flux receives the correction
            !   F'_f = -rho*d_b*S/dn*(p'_f - p'_P) = +af*p'_P.
            ! Without it the outlet column could only rebalance through the
            ! (also mis-set) cell-velocity correction, stalling a few cells
            ! upstream of the exit on a spurious non-parallel state.
            ! Zero-gradient-pressure boundaries (walls, symmetry, fixed-
            ! velocity inlets) keep F'_f = 0.
            gi = bcs%fgrp(i)
            if ( gi > 0 ) then
               if ( bcs%gb(gi)%btype == BC_POUTLET .or. &
                    bcs%gb(gi)%btype == BC_FARFIELD ) then
                  c0  = m%f(i)%c0
                  dbf = g%vol(c0) / fld%apc(c0)
                  dn  = norm2( g%xf(:,i) - g%xc(:,c0) )
                  af  = ctrl%rho * dbf * g%area(i) / dn
                  fld%flux(i) = fld%flux(i) + af * pp(c0)
               end if
            end if
            cycle
         end if

         c0 = m%f(i)%c0
         dbf = fld%lf(i)         * g%vol(c0) / fld%apc(c0) &
             + (1.0_dp-fld%lf(i)) * g%vol(c1) / fld%apc(c1)
         dn = norm2( g%xc(:,c1) - g%xc(:,c0) )
         af = ctrl%rho * dbf * g%area(i) / dn
         fld%flux(i) = fld%flux(i) - af * ( pp(c1) - pp(c0) )
      end do

   end subroutine correct_fields

   !----------------------------------------------------------------------------
   ! Assemble the energy (temperature) equation into A%val / rhs.
   ! Mirrors momentum_assembly for a single scalar but: (a) has no pressure
   ! force, (b) uses cp / k_cond, (c) the convective mass flux is cp*fld%flux
   ! (fld%flux already carries rho, so no extra rho factor).
   ! ap_T = unrelaxed diagonal (handed back for reuse / diagnostics).
   !----------------------------------------------------------------------------
   subroutine temperature_assembly( m, g, ctrl, bcs, fld, A, rhs, ap_T )
      type(mesh_t),   intent(in)    :: m
      type(geom_t),   intent(in)    :: g
      type(ctrl_t),   intent(in)    :: ctrl
      type(bc_t),     intent(in)    :: bcs
      type(fields_t), intent(in)    :: fld
      type(csr_t),    intent(inout) :: A
      real(dp),       intent(out)   :: rhs(m%ncells)
      real(dp),       intent(out)   :: ap_T(m%ncells)

      integer  :: i, kk, c0, c1, kd, gi
      real(dp) :: FT, DT, dn, dw, sd, at, T_face, q_face, phi_ud, dc, noc
      real(dp) :: gf(3), dvec(3), ndsf(3)
      real(dp) :: k_eff0, k_eff1, kf, hasf, kft, kst, kd_f, umag, umag2
      logical  :: is_neumann, is_ltne
      real(dp), allocatable :: psi_T(:)
      real(dp), allocatable :: kdisp(:,:)   ! (3,ncells) dispersion conductivity

      is_ltne = ( ctrl%thermal_model == 'ltne' )

      A%val = 0.0_dp
      rhs   = 0.0_dp
      ap_T  = 0.0_dp
      sd    = 0.0_dp

      ! Barth-Jespersen limiter for the high-order T face value (lagged)
      if ( ctrl%conv_blend > 0.0_dp ) then
         allocate( psi_T(m%ncells) )
         call temperature_limiter( m, g, bcs, fld, psi_T )
      end if

      ! ---- thermal dispersion (Bear component-wise model) ----------------------
      ! Extra fluid-phase conductivity from velocity-driven pore-scale mixing:
      !   D_dd = rho*cp/|u| * ( disp_l*u_d^2 + disp_t*(|u|^2-u_d^2) )  [W/m/K]
      ! i.e. rho*cp*( disp_l|u| along the flow, disp_t|u| transverse ); the
      ! 1/|u| factor is shared by both terms.  Zero in stagnant cells and in
      ! cells without dispersivities (pure fluid).
      allocate( kdisp(3,m%ncells) )
      kdisp = 0.0_dp
      do kk = 1, m%ncells
         if ( fld%disp_l(kk) <= 0.0_dp .and. fld%disp_t(kk) <= 0.0_dp ) cycle
         umag = sqrt( dot_product( fld%u(:,kk), fld%u(:,kk) ) )
         if ( umag < 1.0e-14_dp ) cycle
         umag2 = umag * umag
         do kd = 1, 3
            kdisp(kd,kk) = ctrl%rho * ctrl%cp / umag * ( &
               fld%disp_l(kk) * fld%u(kd,kk)**2 &
               + fld%disp_t(kk) * ( umag2 - fld%u(kd,kk)**2 ) )
         end do
      end do

      ! ---- interior faces ------------------------------------------------------
      do i = 1, m%nfaces
         c1 = m%f(i)%c1
         if ( c1 == 0 ) cycle
         c0 = m%f(i)%c0

         ! convective mass flux for T: cp * (rho*u.S) -- flux already has rho.
         ! Convection is carried by the fluid phase, so it always uses cp.
         FT   = ctrl%cp * fld%flux(i)
         dvec = g%xc(:,c1) - g%xc(:,c0)
         dn   = norm2( dvec )
         ! effective diffusion coefficient:
         !   LTE  : k_eff = eps*k_f + (1-eps)*k_s  (volume-weighted mixture)
         !   LTNE : k_f_diff = eps*k_f             (fluid-phase conduction only)
         if ( is_ltne ) then
            k_eff0 = fld%porosity(c0) * ctrl%k_cond
            k_eff1 = fld%porosity(c1) * ctrl%k_cond
         else
            k_eff0 = fld%porosity(c0)*ctrl%k_cond &
                   + (1.0_dp - fld%porosity(c0)) * fld%k_s(c0)
            k_eff1 = fld%porosity(c1)*ctrl%k_cond &
                   + (1.0_dp - fld%porosity(c1)) * fld%k_s(c1)
         end if
         kf = fld%lf(i) * k_eff0 + (1.0_dp - fld%lf(i)) * k_eff1
         ! dispersion: face-projected diagonal tensor n.D.n (face-averaged);
         ! the non-orthogonal correction below reuses kf, which is exact for
         ! orthogonal meshes and consistent otherwise
         ndsf = ( g%sf(:,i) / g%area(i) )**2
         kd_f = fld%lf(i) * dot_product( ndsf, kdisp(:,c0) ) &
              + ( 1.0_dp - fld%lf(i) ) * dot_product( ndsf, kdisp(:,c1) )
         kf = kf + kd_f
         if ( ctrl%nonorth_corr > 0.0_dp ) then
            sd = dot_product( g%sf(:,i), dvec )
            DT = kf * g%area(i)**2 / sd
         else
            DT = kf * g%area(i) / dn
         end if

         ! upwind convection + implicit diffusion
         ap_T(c0) = ap_T(c0) + max( FT, 0.0_dp ) + DT
         ap_T(c1) = ap_T(c1) + max( -FT, 0.0_dp ) + DT
         kd = pos_of( A, c0, c1 )
         A%val(kd) = A%val(kd) - ( DT + max( -FT, 0.0_dp ) )
         kd = pos_of( A, c1, c0 )
         A%val(kd) = A%val(kd) - ( DT + max( FT, 0.0_dp ) )

         ! deferred correction towards a limited second-order upwind T face
         if ( ctrl%conv_blend > 0.0_dp ) then
            if ( FT > 0.0_dp ) then
               phi_ud = fld%T(c0)
               dc = psi_T(c0) * dot_product( fld%gt(:,c0), &
                                             g%xf(:,i) - g%xc(:,c0) )
            else
               phi_ud = fld%T(c1)
               dc = psi_T(c1) * dot_product( fld%gt(:,c1), &
                                             g%xf(:,i) - g%xc(:,c1) )
            end if
            dc = ctrl%conv_blend * FT * dc
            rhs(c0) = rhs(c0) - dc
            rhs(c1) = rhs(c1) + dc
         end if

         ! non-orthogonal diffusion correction (deferred, lagged)
         if ( ctrl%nonorth_corr > 0.0_dp ) then
            gf  = fld%lf(i) * fld%gt(:,c0) &
                + ( 1.0_dp - fld%lf(i) ) * fld%gt(:,c1)
            noc = ctrl%nonorth_corr * kf * dot_product( gf, &
                     g%sf(:,i) - dvec * ( g%area(i)**2 / sd ) )
            rhs(c0) = rhs(c0) + noc
            rhs(c1) = rhs(c1) - noc
         end if
      end do

      ! ---- boundary faces ------------------------------------------------------
      do i = 1, m%nfaces
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         c0 = m%f(i)%c0
         call bc_face_T( bcs, i, g%xf(:,i), fld%T(c0), T_face, q_face, is_neumann )
         dw = norm2( g%xf(:,i) - g%xc(:,c0) )
         ! owner-cell effective diffusion coefficient (same model as interior)
         if ( is_ltne ) then
            kf = fld%porosity(c0) * ctrl%k_cond
         else
            kf = fld%porosity(c0)*ctrl%k_cond &
               + (1.0_dp - fld%porosity(c0)) * fld%k_s(c0)
         end if
         ! dispersion contribution of the owner cell along the face normal
         ndsf = ( g%sf(:,i) / g%area(i) )**2
         kf = kf + dot_product( ndsf, kdisp(:,c0) )
         DT = kf * g%area(i) / dw
         FT = ctrl%cp * fld%flux(i)              ! rho*u.S already in flux

         if ( .not. is_neumann ) then
            ! Dirichlet (fixed-T wall / inlet)
            ap_T(c0) = ap_T(c0) + DT + max( FT, 0.0_dp )
            rhs(c0)  = rhs(c0)  + DT*T_face - FT*T_face
         else
            ! Neumann (fixed-flux / adiabatic / symmetry / outlet)
            if ( abs(q_face) > 0.0_dp ) then
               if ( is_ltne ) then
                  ! LTNE: split the boundary heat flux between the phases by
                  ! effective conductivity fractions (Nield-style wall split):
                  ! the fluid receives q*eps*k_f/k_tot here and the solid
                  ! q*(1-eps)*k_s/k_tot in solid_temperature_assembly, so the
                  ! sum stays exactly q and no phase is counted twice.
                  kft = fld%porosity(c0) * ctrl%k_cond
                  kst = (1.0_dp - fld%porosity(c0)) * fld%k_s(c0)
                  rhs(c0) = rhs(c0) + q_face * g%area(i) * kft / ( kft + kst )
               else
                  rhs(c0) = rhs(c0) + q_face * g%area(i)  ! q>0 heats domain
               end if
            end if
            ! outflow convection (implicit; zero for walls/symmetry where F=0)
            if ( bcs%gb(gi)%btype == BC_POUTLET .or. &
                 bcs%gb(gi)%btype == BC_OUTFLOW ) &
               ap_T(c0) = ap_T(c0) + max( FT, 0.0_dp )
         end if
      end do

      ! ---- LTNE fluid-solid interfacial heat source ----------------------------
      ! Source term +h_sf*a_sf*(T_s - T_f): the T_f part is treated implicitly
      ! (ap_T += h_sf*a_sf*V) and the T_s part explicitly (rhs += h_sf*a_sf*T_s*V).
      if ( is_ltne ) then
         do kk = 1, m%ncells
            hasf = fld%h_sf(kk) * fld%a_sf(kk)
            if ( hasf > 0.0_dp ) then
               ap_T(kk) = ap_T(kk) + hasf * g%vol(kk)
               rhs(kk)  = rhs(kk)  + hasf * fld%T_s(kk) * g%vol(kk)
            end if
         end do
      end if

      ! ---- transient time term ------------------------------------------------
      ! LTE : at = (rho*cp)_eff * V / dt, with (rho*cp)_eff = eps*(rho*cp)_f
      !                                              + (1-eps)*(rho*cp)_s
      ! LTNE: at = eps*(rho*cp)_f * V / dt   (fluid-phase heat capacity only)
      select case ( fld%ts_order )
      case ( 1 )
         do kk = 1, m%ncells
            if ( is_ltne ) then
               at = fld%porosity(kk) * ctrl%rho * ctrl%cp * g%vol(kk) / ctrl%dt
            else
               at = ( fld%porosity(kk) * ctrl%rho * ctrl%cp &
                    + (1.0_dp - fld%porosity(kk)) * fld%rho_s(kk) * fld%cp_s(kk) ) &
                    * g%vol(kk) / ctrl%dt
            end if
            ap_T(kk) = ap_T(kk) + at
            rhs(kk)  = rhs(kk)  + at * fld%T_old(kk)
         end do
      case ( 2 )
         do kk = 1, m%ncells
            if ( is_ltne ) then
               at = fld%porosity(kk) * ctrl%rho * ctrl%cp * g%vol(kk) &
                    / ( 2.0_dp * ctrl%dt )
            else
               at = ( fld%porosity(kk) * ctrl%rho * ctrl%cp &
                    + (1.0_dp - fld%porosity(kk)) * fld%rho_s(kk) * fld%cp_s(kk) ) &
                    * g%vol(kk) / ( 2.0_dp * ctrl%dt )
            end if
            ap_T(kk) = ap_T(kk) + 3.0_dp * at
            rhs(kk)  = rhs(kk)  + at * ( 4.0_dp * fld%T_old(kk) &
                                         - fld%T_old_old(kk) )
         end do
      end select

      ! ---- diagonal + under-relaxation --------------------------------------------
      ! Same three-mode logic as momentum_assembly: PISO skips relaxation,
      ! SIMPLE and PIMPLE apply alpha_u (reused for the scalar equation to
      ! keep the parameter set lean).
      if ( fld%ts_order > 0 .and. .not. ctrl%pimple ) then
         do kk = 1, m%ncells
            kd = pos_of( A, kk, kk )
            A%val(kd) = ap_T(kk)
         end do
      else
         do kk = 1, m%ncells
            kd = pos_of( A, kk, kk )
            A%val(kd) = ap_T(kk) / ctrl%alpha_u
            rhs(kk) = rhs(kk) &
                    + ( 1.0_dp - ctrl%alpha_u ) / ctrl%alpha_u * ap_T(kk) * fld%T(kk)
         end do
      end if

      if ( allocated(psi_T) ) deallocate( psi_T )

   end subroutine temperature_assembly

   !----------------------------------------------------------------------------
   ! Solid-phase energy equation (LTNE, phase 12).
   !
   !   (1-eps)*rho_s*cp_s * dT_s/dt = div( (1-eps)*k_s grad T_s )
   !                                 + h_sf*a_sf*(T_f - T_s)
   !
   ! No convection in the solid.  The wall thermal BC (fixed-T or fixed-flux)
   ! is shared with the fluid phase via bc_face_T.  The interfacial heat
   ! exchange is split implicitly (T_s -> diagonal) and explicitly (T_f -> rhs).
   !----------------------------------------------------------------------------
   subroutine solid_temperature_assembly( m, g, ctrl, bcs, fld, A, rhs, ap_Ts )
      type(mesh_t),   intent(in)    :: m
      type(geom_t),   intent(in)    :: g
      type(ctrl_t),   intent(in)    :: ctrl
      type(bc_t),     intent(in)    :: bcs
      type(fields_t), intent(in)    :: fld
      type(csr_t),    intent(inout) :: A
      real(dp),       intent(out)   :: rhs(m%ncells)
      real(dp),       intent(out)   :: ap_Ts(m%ncells)

      integer  :: i, kk, c0, c1, kd, gi
      real(dp) :: Ds, dn, dw, sd, at, T_face, q_face, noc, hasf
      real(dp) :: gf(3), dvec(3), k_s0, k_s1, ksf, kft
      logical  :: is_neumann

      A%val = 0.0_dp
      rhs   = 0.0_dp
      ap_Ts = 0.0_dp
      sd    = 0.0_dp

      ! ---- interior faces ------------------------------------------------------
      do i = 1, m%nfaces
         c1 = m%f(i)%c1
         if ( c1 == 0 ) cycle
         c0 = m%f(i)%c0

         dvec = g%xc(:,c1) - g%xc(:,c0)
         dn   = norm2( dvec )
         ! solid-phase diffusion coefficient: (1-eps)*k_s
         k_s0 = (1.0_dp - fld%porosity(c0)) * fld%k_s(c0)
         k_s1 = (1.0_dp - fld%porosity(c1)) * fld%k_s(c1)
         ksf  = fld%lf(i) * k_s0 + (1.0_dp - fld%lf(i)) * k_s1
         if ( ctrl%nonorth_corr > 0.0_dp ) then
            sd = dot_product( g%sf(:,i), dvec )
            Ds = ksf * g%area(i)**2 / sd
         else
            Ds = ksf * g%area(i) / dn
         end if

         ap_Ts(c0) = ap_Ts(c0) + Ds
         ap_Ts(c1) = ap_Ts(c1) + Ds
         kd = pos_of( A, c0, c1 )
         A%val(kd) = A%val(kd) - Ds
         kd = pos_of( A, c1, c0 )
         A%val(kd) = A%val(kd) - Ds

         ! non-orthogonal diffusion correction (deferred, lagged)
         if ( ctrl%nonorth_corr > 0.0_dp ) then
            gf  = fld%lf(i) * fld%gts(:,c0) &
                + ( 1.0_dp - fld%lf(i) ) * fld%gts(:,c1)
            noc = ctrl%nonorth_corr * ksf * dot_product( gf, &
                     g%sf(:,i) - dvec * ( g%area(i)**2 / sd ) )
            rhs(c0) = rhs(c0) + noc
            rhs(c1) = rhs(c1) - noc
         end if
      end do

      ! ---- boundary faces ------------------------------------------------------
      do i = 1, m%nfaces
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         c0 = m%f(i)%c0
         call bc_face_T( bcs, i, g%xf(:,i), fld%T_s(c0), T_face, q_face, is_neumann )
         ! The solid phase has no through-flow: at flow-inlet faces the shared
         ! thermal BC does not apply -- use an adiabatic (zero-flux) condition.
         if ( bcs%gb(gi)%btype == BC_VINLET .or. &
              bcs%gb(gi)%btype == BC_VINLET_PARAB .or. &
              bcs%gb(gi)%btype == BC_MASSINLET ) then
            is_neumann = .true.
            q_face     = 0.0_dp
         end if
         dw = norm2( g%xf(:,i) - g%xc(:,c0) )
         ksf = (1.0_dp - fld%porosity(c0)) * fld%k_s(c0)
         kft = fld%porosity(c0) * ctrl%k_cond
         Ds  = ksf * g%area(i) / dw

         if ( .not. is_neumann ) then
            ! Dirichlet (fixed-T wall)
            ap_Ts(c0) = ap_Ts(c0) + Ds
            rhs(c0)   = rhs(c0)  + Ds * T_face
         else
            ! Neumann (fixed-flux / adiabatic / symmetry / outlet): the solid
            ! receives its conductivity fraction of the boundary heat flux
            ! (complement of the fluid share in temperature_assembly)
            if ( abs(q_face) > 0.0_dp ) &
               rhs(c0) = rhs(c0) + q_face * g%area(i) * ksf / ( kft + ksf )
         end if
      end do

      ! ---- fluid-solid interfacial heat source --------------------------------
      ! +h_sf*a_sf*(T_f - T_s): T_s implicit (diagonal), T_f explicit (rhs).
      do kk = 1, m%ncells
         hasf = fld%h_sf(kk) * fld%a_sf(kk)
         if ( hasf > 0.0_dp ) then
            ap_Ts(kk) = ap_Ts(kk) + hasf * g%vol(kk)
            rhs(kk)   = rhs(kk)  + hasf * fld%T(kk) * g%vol(kk)
         end if
      end do

      ! ---- transient time term ------------------------------------------------
      ! at = (1-eps)*rho_s*cp_s * V / dt
      select case ( fld%ts_order )
      case ( 1 )
         do kk = 1, m%ncells
            at = (1.0_dp - fld%porosity(kk)) * fld%rho_s(kk) * fld%cp_s(kk) &
                 * g%vol(kk) / ctrl%dt
            ap_Ts(kk) = ap_Ts(kk) + at
            rhs(kk)   = rhs(kk)  + at * fld%T_s_old(kk)
         end do
      case ( 2 )
         do kk = 1, m%ncells
            at = (1.0_dp - fld%porosity(kk)) * fld%rho_s(kk) * fld%cp_s(kk) &
                 * g%vol(kk) / ( 2.0_dp * ctrl%dt )
            ap_Ts(kk) = ap_Ts(kk) + 3.0_dp * at
            rhs(kk)   = rhs(kk)  + at * ( 4.0_dp * fld%T_s_old(kk) &
                                          - fld%T_s_old_old(kk) )
         end do
      end select

      ! ---- diagonal + under-relaxation -----------------------------------------
      if ( fld%ts_order > 0 .and. .not. ctrl%pimple ) then
         do kk = 1, m%ncells
            kd = pos_of( A, kk, kk )
            A%val(kd) = ap_Ts(kk)
         end do
      else
         do kk = 1, m%ncells
            kd = pos_of( A, kk, kk )
            A%val(kd) = ap_Ts(kk) / ctrl%alpha_u
            rhs(kk) = rhs(kk) &
                    + ( 1.0_dp - ctrl%alpha_u ) / ctrl%alpha_u * ap_Ts(kk) * fld%T_s(kk)
         end do
      end if

   end subroutine solid_temperature_assembly

   !----------------------------------------------------------------------------
   ! Barth-Jespersen slope limiter for the temperature field (mirrors
   ! velocity_limiter with fld%T / fld%gt and bc_face_T for boundary values).
   !----------------------------------------------------------------------------
   subroutine temperature_limiter( m, g, bcs, fld, psi_T )
      type(mesh_t),   intent(in)  :: m
      type(geom_t),   intent(in)  :: g
      type(bc_t),     intent(in)  :: bcs
      type(fields_t), intent(in)  :: fld
      real(dp),       intent(out) :: psi_T(m%ncells)

      integer  :: i, c0, c1, gi
      real(dp) :: pmin(m%ncells), pmax(m%ncells)
      real(dp) :: tc, fr, df, tf, qf
      logical  :: neum

      do i = 1, m%ncells
         pmin(i) = fld%T(i)
         pmax(i) = fld%T(i)
      end do
      do i = 1, m%nfaces
         c0 = m%f(i)%c0
         c1 = m%f(i)%c1
         if ( c1 > 0 ) then
            pmin(c0) = min( pmin(c0), fld%T(c1) )
            pmax(c0) = max( pmax(c0), fld%T(c1) )
            pmin(c1) = min( pmin(c1), fld%T(c0) )
            pmax(c1) = max( pmax(c1), fld%T(c0) )
         else
            gi = bcs%fgrp(i)
            if ( gi > 0 ) then
               call bc_face_T( bcs, i, g%xf(:,i), fld%T(c0), tf, qf, neum )
               pmin(c0) = min( pmin(c0), tf )
               pmax(c0) = max( pmax(c0), tf )
            end if
         end if
      end do

      psi_T = 1.0_dp
      do i = 1, m%nfaces
         c0 = m%f(i)%c0
         tc = fld%T(c0)
         fr = tc + dot_product( fld%gt(:,c0), g%xf(:,i) - g%xc(:,c0) )
         df = fr - tc
         if ( df > 0.0_dp ) then
            psi_T(c0) = min( psi_T(c0), ( pmax(c0) - tc ) / df )
         else if ( df < 0.0_dp ) then
            psi_T(c0) = min( psi_T(c0), ( pmin(c0) - tc ) / df )
         end if
         c1 = m%f(i)%c1
         if ( c1 > 0 ) then
            tc = fld%T(c1)
            fr = tc + dot_product( fld%gt(:,c1), g%xf(:,i) - g%xc(:,c1) )
            df = fr - tc
            if ( df > 0.0_dp ) then
               psi_T(c1) = min( psi_T(c1), ( pmax(c1) - tc ) / df )
            else if ( df < 0.0_dp ) then
               psi_T(c1) = min( psi_T(c1), ( pmin(c1) - tc ) / df )
            end if
         end if
      end do
      psi_T = max( 0.0_dp, min( 1.0_dp, psi_T ) )

   end subroutine temperature_limiter

end module mod_uns_simple
