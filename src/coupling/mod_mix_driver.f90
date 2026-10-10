!===============================================================================
! mod_mix_driver.f90 -- reusable weak-coupling group drivers (phase "mixnsolver").
!
! Extracts the two group drivers (and their helpers) that were internal
! `contains` subprograms of `program mixsolver` (src/main.f90) so that both the
! legacy mixed driver (bin/mixsolver) and the new single self-dispatching
! executable (bin/mixnsolver) can reuse the exact same coupling loop without
! duplicating it.
!
! All 6 procedures are moved VERBATIM from src/main.f90 -- the legacy
! `bin/mixsolver` behaviour is bit-identical (verified by the couple_porous
! regression).  Only two are public: struct_group_driver / uns_group_driver.
!===============================================================================
module mod_mix_driver
   use mpi
   use mod_precision, only: dp
   implicit none
   private

   public :: struct_group_driver
   public :: uns_group_driver

contains

   !===========================================================================
   ! Structured-group driver.
   ! Initialises the structured solver, then runs the coupling loop:
   !   1. Extract interface state (non-dim), convert to SI, send to uns root
   !   2. Receive uns SI state from uns root, convert to non-dim, set ghost cells
   !   3. Advance structured solver by one time step
   !===========================================================================
   subroutine struct_group_driver( comm, casfile, ctlfile, rank, nproc )
      use mod_struct_driver, only: struct_solver_init, struct_solver_step, &
                                    struct_extract_iface, struct_set_iface_bc, &
                                    struct_solver_save
      use mod_coupling_exchange, only: exchange_struct_to_uns, recv_uns_iface_state
      use mod_interface_units, only: struct_to_SI, SI_to_struct
      use mod_reference_state, only: get_ref_state, reference_state_t, &
                                     get_coupling_params, scheduled_steps
      implicit none

      integer,          intent(in) :: comm, rank, nproc
      character(len=*), intent(in) :: casfile, ctlfile

      integer :: iter, iter0, n_couple, n_uns_root, ierr2
      integer :: nfaces, n_uns_recv, n_uns_steps_d, iface_ramp_d
      integer :: n_struct_steps_d, isub, save_interval, couple_restart
      integer :: n_struct_steps_start_d, step_decay_d, n_sub
      integer :: iface_p_model_d
      real(dp) :: iface_relax_d
      real(dp), allocatable :: rho_nd(:), u_nd(:), v_nd(:), w_nd(:), T_nd(:), p_nd(:)
      real(dp), allocatable :: urho(:), uu(:,:), uT(:), up(:)
      real(dp), allocatable :: rho_nd2(:), u_nd2(:), v_nd2(:), w_nd2(:), T_nd2(:), p_nd2(:)
      real(dp), allocatable :: rho_si(:), u_si(:), v_si(:), w_si(:), T_si(:), p_si(:)
      type(reference_state_t) :: r

      call get_coupling_params( n_couple, n_uns_steps_d, iface_relax_d, iface_ramp_d, &
                                n_struct_steps = n_struct_steps_d, &
                                save_interval = save_interval, &
                                couple_restart = couple_restart, &
                                n_struct_steps_start = n_struct_steps_start_d, &
                                step_decay_every = step_decay_d, &
                                iface_p_model = iface_p_model_d )
      n_uns_root = 1   ! uns root is global rank 1 (struct is rank 0)
      if ( n_struct_steps_start_d > n_struct_steps_d ) then
         write(*,'(a,i0,a,i0,a)') '[struct rank 0] struct substeps/coupling iter: ', &
               n_struct_steps_start_d, ' (halving every ', step_decay_d, &
               ' coupling iter) -> steady '
         write(*,'(a,i0)') '[struct rank 0]   steady struct substeps: ', n_struct_steps_d
      else
         write(*,'(a,i0)') '[struct rank 0] struct pseudo-time substeps/coupling iter: ', &
                           n_struct_steps_d
      end if

      ! joint restart: recover the last saved coupling iteration count
      iter0 = 0
      if ( couple_restart == 1 ) call read_couple_state( iter0 )

      write(*,'(a,i0,a)') '[struct rank ', rank, '] initialising structured solver'
      call struct_solver_init(comm, ctlfile, force_restart = (iter0 > 0) )
      write(*,'(a,i0,a)') '[struct rank ', rank, '] init OK'
      if ( iter0 > 0 ) &
         write(*,'(a,i0,a,i0)') '[struct rank ', rank, '] restart from coupling iter ', &
                                iter0

      ! ---- coupling loop ------------------------------------------------------
      do iter = iter0+1, n_couple
         ! extract interface state (non-dim)
         call struct_extract_iface(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, nfaces)
         if (iter == 1) &
            write(*,'(a,i0,a,i0,a)') '[struct rank ', rank, '] interface faces: ', nfaces

         ! send struct SI state to uns root (count + payload, works for nfaces=0)
         if (nfaces > 0) then
            call exchange_struct_to_uns(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                        nfaces, n_uns_root, ierr2)
         else
            call exchange_struct_to_uns(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                        0, n_uns_root, ierr2)
         end if

         ! receive uns SI state from uns root
         call recv_uns_iface_state(n_uns_recv, urho, uu, uT, up, n_uns_root, ierr2)
         if (iter == 1) &
            write(*,'(a,i0,a,i0,a)') '[struct rank ', rank, '] received uns faces: ', n_uns_recv

         ! periodic SI diagnostics (interface means, both peers)
         if ( mod(iter,25)==0 .or. iter==1 .or. iter==n_couple ) then
            if (nfaces > 0) then
               allocate(rho_si(nfaces), u_si(nfaces), v_si(nfaces), &
                        w_si(nfaces), T_si(nfaces), p_si(nfaces))
               call convert_to_si(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                  rho_si, u_si, v_si, w_si, T_si, p_si, nfaces)
               write(*,'(a,i0,a,i0,a,3es12.4)') &
                  '[struct rank ', rank, '] diag iter ', iter, &
                  ' struct SI: mean u_x, mean p(abs), mean T =', &
                  sum(u_si)/nfaces, sum(p_si)/nfaces, sum(T_si)/nfaces
               deallocate(rho_si, u_si, v_si, w_si, T_si, p_si)
            end if
            if (n_uns_recv > 0) then
               write(*,'(a,i0,a,i0,a,3es12.4)') &
                  '[struct rank ', rank, '] diag iter ', iter, &
                  ' uns    SI: mean u_x, mean p(gauge), mean T =', &
                  sum(uu(1,:))/n_uns_recv, sum(up)/n_uns_recv, &
                  sum(uT)/n_uns_recv
            end if
         end if

         ! advance structured solver (subcycled pseudo-time steps to keep
         ! pace with the uns side, which takes several SIMPLE outer steps).
         ! Phase-13 start-up schedule: the sub-step count decays geometrically
         ! (e.g. 1000, 500, 250, ...) until it reaches the steady cadence, so
         ! the two sides exchange more and more frequently and finally at the
         ! same frequency.
         n_sub = scheduled_steps( iter, n_struct_steps_start_d, n_struct_steps_d, &
                                  step_decay_d )
         if ( n_struct_steps_start_d > n_struct_steps_d ) &
            write(*,'(a,i0,a,i0)') '[struct rank ', rank, &
                  '] coupling iter ', iter, ' -> struct substeps: ', n_sub
         do isub = 1, n_sub
            call struct_solver_step()
         end do

         ! set ghost cells from the received uns first-cell state.  Done
         ! after the substeps.  The characteristic back-pressure strength is
         ! alpha = ALPHA_MAX * min(1, iter/iface_ramp): the ramp cushions the
         ! cold start, and ALPHA_MAX<1 permanently damps the partitioned
         ! coupling loop (full alpha=1 was empirically unstable at Ma=0.1:
         ! the acoustic gain 1/(rho*c) makes the exchange loop gain > 1).
         ! At convergence p1 -> p_peer for any alpha>0 because the struct
         ! inner pressure relaxes toward pb every coupling iteration.
         if (nfaces > 0 .and. n_uns_recv > 0) then
            block
               real(dp), parameter :: ALPHA_MAX = 0.3_dp
               real(dp) :: alpha_p
            allocate(rho_nd2(nfaces), u_nd2(nfaces), v_nd2(nfaces), &
                     w_nd2(nfaces), T_nd2(nfaces), p_nd2(nfaces))
            r = get_ref_state()
            call convert_uns_to_struct_nd(urho, uu, uT, up, n_uns_recv, &
                                          rho_nd2, u_nd2, v_nd2, w_nd2, T_nd2, p_nd2, &
                                          nfaces, r)
            alpha_p = ALPHA_MAX * min(1.0_dp, &
                        real(iter,dp)/real(max(iface_ramp_d,1),dp))
            ! Phase 14 (iface_p_model = dirichlet|momentum): the interface
            ! pressure is now a DIRICHLET condition on the uns side, taken from
            ! THIS side -- the struct is the pressure master.  Relaxing the
            ! struct interface pressure towards the slab pressure as well would
            ! close a two-way pressure loop (struct -> slab -> struct) whose
            ! measured gain in C_P_test is ~5-10x per coupling iteration, i.e.
            ! divergent.  With alpha = 0 the struct keeps extrapolating its own
            ! interface pressure (consistent with the mainstream) while still
            ! taking the slab mass flux through the ghost VELOCITY, which is
            ! the Dirichlet-Neumann partition the literature prescribes.
            if ( iface_p_model_d /= 0 ) alpha_p = 0.0_dp
            call struct_set_iface_bc(rho_nd2, u_nd2, v_nd2, w_nd2, T_nd2, p_nd2, nfaces, &
                                     alpha = alpha_p)
            deallocate(rho_nd2, u_nd2, v_nd2, w_nd2, T_nd2, p_nd2)
            end block
         else if (nfaces > 0) then
            ! no uns data received; keep current state as BC (extrapolation)
            call struct_set_iface_bc(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, nfaces, &
                                     alpha = 0.0_dp)
         end if
         if ( mod(iter,100)==0 ) &
            write(*,'(a,i0,a,i0,a)') '[struct rank ', rank, '] coupling iter ', iter, ' done'

         ! phase-10 auto-save: aligned with the uns side at the same coupling
         ! iter (after both sides have exchanged and stepped).
         if ( save_interval > 0 ) then
            if ( mod(iter, save_interval) == 0 ) then
               call struct_solver_save
               if ( rank == 0 ) call write_couple_state( iter )
               write(*,'(a,i0,a,i0,a)') '[struct rank ', rank, &
                     '] saved flow3d.dat at coupling iter ', iter
            end if
         end if

         if (allocated(rho_nd)) deallocate(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd)
         if (allocated(urho)) deallocate(urho, uu, uT, up)
      end do

      write(*,'(a,i0,a)') '[struct rank ', rank, '] coupling loop done'
   end subroutine struct_group_driver

   !===========================================================================
   ! Unstructured-group driver.
   ! Initialises the unstructured solver, then runs the coupling loop:
   !   1. Receive struct SI state from struct root, set interface BC
   !   2. Advance SIMPLE by n_uns_steps
   !   3. Extract interface SI state, send to struct root
   !===========================================================================
   subroutine uns_group_driver( comm, casfile, ctlfile, rank, nproc )
      use mod_uns_driver, only: uns_solver_init, uns_solver_step, &
                                uns_solver_extract_iface
      use mod_uns_mesh, only: mesh_t
      use mod_uns_connectivity, only: conn_t
      use mod_uns_geometry, only: geom_t
      use mod_uns_control, only: ctrl_t, BC_INTERFACE, BC_POUTLET, BC_FARFIELD
      use mod_uns_bc, only: bc_t, set_interface_vel, &
                            set_interface_T, bc_face_vel, set_interface_p, &
                            set_interface_slip, set_iface_p_dirichlet, &
                            set_iface_slip_mode, set_interface_T_ltne, &
                            set_iface_zhang_vel, &
                            set_iface_t_ltne_mode, iface_t_ltne_is_on
      use mod_iface_law, only: iface_slip_conductance, iface_p_momentum, &
                               iface_p_blend, iface_p_flux_response, &
                               iface_p_gain, iface_zhang_velocity, &
                               iface_zhang_T
      use mod_uns_fields, only: fields_t
#ifdef HAVE_MPI
      ! The uns field dump/restart lives in the MPI-only mod_uns_restart stack
      ! (-> mod_uns_mpi_core / mod_uns_partition / mod_uns_local_mesh).  The
      ! serial 'single-process MPI' build (make all) omits that stack, so the
      ! save path below is compiled out there; the MPI build enables it.
      use mod_uns_restart, only: write_field_dump
#endif
      use mod_coupling_exchange, only: exchange_uns_to_struct
      use mod_reference_state, only: get_coupling_params, scheduled_steps, &
                                     get_ref_state, reference_state_t
      implicit none

      integer,          intent(in) :: comm, rank, nproc
      character(len=*), intent(in) :: casfile, ctlfile

      type(mesh_t)   :: m
      type(conn_t)   :: c
      type(geom_t)   :: g
      type(ctrl_t)   :: ctrl
      type(bc_t)     :: bcs
      type(fields_t) :: fld
      integer :: ier, ierr, iter, iter0, n_couple, n_uns_steps, n_struct_root
      integer, allocatable :: ifaces(:)
      real(dp), allocatable :: irho(:), iT(:), ip(:), iu(:,:)
      real(dp), allocatable :: srho(:), su(:,:), sT(:), sp(:)
      real(dp), allocatable :: su_bc(:,:)
      integer :: n_uns_faces, iface_ramp, save_interval, couple_restart
      integer :: n_uns_steps_start, step_decay, n_uns_cur, iface_vel_mode
      integer :: iface_p_anchor
      integer :: iface_p_model, iface_slip, iface_pm_sub
      real(dp) :: iface_p_w, iface_slip_alpha
      integer :: iface_t_model
      real(dp) :: iface_beta, iface_beta1, iface_df_ratio
      logical :: iface_p_dir
      real(dp) :: iface_relax, omega
      logical :: has_struct
      character(len=*), parameter :: uns_dump = 'unMesh_restart.dat'

      call get_coupling_params( n_couple, n_uns_steps, iface_relax, iface_ramp, &
                                save_interval = save_interval, &
                                couple_restart = couple_restart, &
                                n_uns_steps_start = n_uns_steps_start, &
                                step_decay_every = step_decay, &
                                iface_vel_mode = iface_vel_mode, &
                                iface_p_anchor = iface_p_anchor, &
                                iface_p_model = iface_p_model, &
                                iface_slip = iface_slip, &
                                iface_pm_subiter = iface_pm_sub, &
                                iface_p_blend = iface_p_w, &
                                iface_slip_alpha = iface_slip_alpha, &
                                 iface_t_model = iface_t_model, &
                                 iface_beta = iface_beta, &
                                 iface_beta1 = iface_beta1, &
                                 iface_df_ratio = iface_df_ratio )
      n_struct_root = 0   ! struct root is global rank 0
      ! Phase-14 interface closure switches (module-level flags in mod_uns_bc;
      ! the standalone uns solver leaves them at .false., so its behaviour is
      ! unchanged).  iface_p_model /= 0 makes the coupling interface a
      ! pressure-Dirichlet patch, which also switches the interface velocity
      ! treatment inside the momentum assembly.
      iface_p_dir = ( iface_p_model /= 0 )
      call set_iface_p_dirichlet( iface_p_dir )
      call set_iface_slip_mode( iface_slip == 1 )
      ! iface_velocity = zhang: the interface velocity is a viscous closure
      ! value (Zhang Eqs.22/25/26), so the interface face must not convect the
      ! prescribed tangential slip into the first slab cell (see mod_uns_bc).
      call set_iface_zhang_vel( iface_vel_mode == 3 )
      if ( rank == 0 ) then
         write(*,'(a,i0,a,i0,a,i0)') '[uns] iface_p_model = ', iface_p_model, &
            '  iface_slip = ', iface_slip, '  pm_subiter = ', iface_pm_sub
      end if
      ! When there is no struct group (nproc==1, n_struct_ranks==0) the
      ! exchange is skipped entirely.
      has_struct = ( nproc > 1 )

      ! joint restart: recover the last saved coupling iteration count
      iter0 = 0
      if ( couple_restart == 1 ) call read_couple_state( iter0 )

      write(*,'(a,i0,a)') '[uns rank ', rank, '] initialising unstructured solver'
      if ( iter0 > 0 ) then
         call uns_solver_init( casfile, ctlfile, m, c, g, ctrl, bcs, fld, ier, &
                               restart_file = uns_dump )
      else
         call uns_solver_init( casfile, ctlfile, m, c, g, ctrl, bcs, fld, ier )
      end if
      if ( ier /= 0 ) then
         write(*,'(a,i0,a,i0)') '[uns rank ', rank, '] init FAILED ier=', ier
         call MPI_Abort( MPI_COMM_WORLD, ier, ierr )
         return
      end if
      ! ---- phase-14 Zhang 2011 interface closures ------------------------------
      ! Velocity (Eqs.22/25/26) is handled per face in the coupling loop; the
      ! temperature side (Eqs.27-29) needs the thermal model to decide whether
      ! there is a flux to split (LTNE) or the volume average degenerates to
      ! the exchanged value (LTE, T_f = T_s).
      call set_iface_t_ltne_mode( iface_t_model == 1 .and. &
                                  trim(ctrl%thermal_model) == 'ltne' )
      if ( rank == 0 ) then
         if ( iface_vel_mode == 3 .and. iface_p_dir ) &
            write(*,'(a)') '[uns] WARNING: iface_velocity=zhang is overridden ' // &
               'by iface_p_model=dirichlet|momentum (the pressure-Dirichlet ' // &
               'interface extrapolates the velocity; use iface_p_model=grad0)'
         if ( iface_t_model == 1 ) &
            write(*,'(a,l1)') '[uns] iface_t_model=zhang, LTNE split on = ', &
                              iface_t_ltne_is_on()
      end if
      write(*,'(a,i0,a,i0)') '[uns rank ', rank, '] init OK, ncells=', m%ncells
      if ( iter0 > 0 ) &
         write(*,'(a,i0,a,i0)') '[uns rank ', rank, '] restart from coupling iter ', iter0

      ! ---- extract initial interface state (for exchange sizing) ---------------
      call uns_solver_extract_iface( m, g, bcs, fld, ifaces, irho, iu, iT, ip, ier, ctrl%rho )
      n_uns_faces = size(ifaces)
      write(*,'(a,i0,a,i0)') '[uns rank ', rank, '] interface faces: ', n_uns_faces

      ! ---- coupling loop ------------------------------------------------------
      do iter = iter0+1, n_couple

         ! receive struct SI state, send uns SI state back
         ! (skipped when there is no struct group, e.g. nproc==1)
         if ( has_struct ) then
            call exchange_uns_to_struct(irho, iu, iT, ip, n_uns_faces, &
                                        n_struct_root, srho, su, sT, sp, ierr)

            ! Optional interface pressure datum anchoring (iface_p_anchor=1).
            ! The uns PPE is pure-Neumann (no pressure BC: the slab is closed
            ! apart from the coolant inlet and the interface), so its gauge
            ! level is a free parameter that drifts.  The struct consumes
            ! (gauge + p_ref) as an ABSOLUTE back pressure, so an un-anchored
            ! datum lets the two sides' interface pressures drift apart until
            ! the loop goes unstable.  Shift the whole uns gauge field so the
            ! mean interface absolute pressure equals the struct's; a constant
            ! shift leaves every gradient (i.e. the uns solution) untouched.
            !
            ! Phase 14 note (iface_p_model = dirichlet|momentum): the anchor is
            ! then NOT redundant -- and is in fact essential.  The interface
            ! pressure becomes a Dirichlet value taken from the peer, so a
            ! cold-start peer interface pressure that is O(10 kPa) off would be
            ! imposed on the slab and the loop diverges.  The anchor first
            ! re-datums the uns field onto the peer's level (a pure shift, no
            ! solution change), after which the Dirichlet update is only the
            ! small physical difference (the Eq.44 momentum term plus the
            ! slab/face extrapolation difference) instead of an absolute value.
            if ( iface_p_anchor == 1 .and. n_uns_faces > 0 ) then
               block
                  type(reference_state_t) :: rr
                  real(dp) :: p_uns_abs, p_str_abs, dp_anchor
                  rr        = get_ref_state()
                  p_uns_abs = sum(ip) / real(n_uns_faces, dp) + rr%p_ref
                  p_str_abs = sum(sp) / real(n_uns_faces, dp)
                  dp_anchor = p_str_abs - p_uns_abs
                  fld%p     = fld%p + dp_anchor
                  if ( iter <= 8 .or. mod(iter,25) == 0 ) &
                     write(*,'(a,i0,a,3es12.4)') '[uns rank ', rank, &
                        '] p-anchor: p_uns_abs, p_str_abs, shift =', &
                        p_uns_abs, p_str_abs, dp_anchor
               end block
            end if

            ! set interface BC from struct state with under-relaxation and
            ! linear ramp over the first iface_ramp coupling iterations:
            !   omega  = iface_relax * min(1, iter/iface_ramp)
            !   BC     = omega*struct + (1-omega)*uns_current
            ! (iu/ip hold the uns interface state from the previous extract)
            if (n_uns_faces > 0) then
               omega = iface_relax * min( 1.0_dp, &
                            real(iter,dp) / real(max(iface_ramp,1),dp) )
               allocate( su_bc(3,n_uns_faces) )
               if ( iface_vel_mode == 2 .or. iface_p_dir ) then
                  ! 'balance' (iface_velocity = balance): the peer velocity is
                  ! ignored; the interface is the slab's own fully-developed
                  ! outflow (zero gradient, i.e. the previous uns interface
                  ! state).  Its total flux is then rescaled by the mass-balance
                  ! correction below so that exactly the injected coolant leaves
                  ! through it.
                  ! Phase 14 (iface_p_model = dirichlet|momentum): the interface
                  ! velocity must NOT be Dirichlet-ised either -- the pressure is
                  ! the interface constraint -- so the same own extrapolation is
                  ! used for the normal part (it also keeps the energy equation's
                  ! interface enthalpy flux, which is built from iface_vel,
                  ! equal to the actual transpiration flux instead of the peer's
                  ! ~100 m/s tangential slip).
                  su_bc = iu
               else
                  su_bc = omega * su + (1.0_dp - omega) * iu
               end if
               ! Optional normal-only interface velocity (iface_velocity =
               ! normal in mix.control).  Only the interface-normal component
               ! is taken from the peer (that is the exchanged mass flux, the
               ! Dirichlet half of the Dirichlet-Neumann partition); the
               ! tangential components keep the local uns interface state
               ! (zero-gradient).  Forcing the peer free-stream tangential
               ! slip (~100 m/s) into the porous slab would require
               ! dp/dx = mu*u/K ~ 5e5 Pa/m and drives the coupled loop to NaN.
               if ( iface_vel_mode == 1 .and. .not. iface_p_dir ) then
                  block
                     integer  :: fi2, kf2
                     real(dp) :: nvec2(3), un2, nrm2
                     do fi2 = 1, n_uns_faces
                        kf2   = ifaces(fi2)
                        nvec2 = g%sf(:,kf2)
                        nrm2  = norm2( nvec2 )
                        if ( nrm2 > 0.0_dp ) then
                           nvec2 = nvec2 / nrm2
                           un2   = dot_product( su_bc(:,fi2), nvec2 )
                           su_bc(:,fi2) = ( iu(:,fi2) &
                                          - dot_product( iu(:,fi2), nvec2 )*nvec2 ) &
                                          + un2 * nvec2
                        end if
                     end do
                  end block
               end if
               ! ---- Zhang 2011 Eqs.22/25/26: literature interface velocity ----
               ! ONE interface velocity is shared by the two domains (Eq.22);
               ! its normal part solves the two-sided normal stress balance
               ! (Eq.25, a conductance/reciprocal-distance blend of the two
               ! near-interface cell velocities) and its tangential part the
               ! stress-jump balance (Eq.26, Ochoa-Tapia & Whitaker with the
               ! excess viscous coefficient beta and the inertial beta1).  The
               ! result is imposed as a full-vector Dirichlet on the slab, so
               ! the interface has NO unconstrained tangential mode (the
               ! divergence seen with 'balance'/'full', where the peer-free
               ! tangential slip was either ignored or copied as ~100 m/s);
               ! the peer velocity enters the blend instead of replacing the
               ! local value.  beta defaults to eps*bj_alpha, the
               ! correspondence that makes this law identical to the validated
               ! internal fluid/porous face treatment of cases/beavers_joseph.
               if ( iface_vel_mode == 3 .and. .not. iface_p_dir ) then
                  block
                     integer  :: fi3, kf3, c03
                     real(dp) :: nv3(3), ar3, dp3, eps3, lam3, beta3, Vz(3)
                     do fi3 = 1, n_uns_faces
                        kf3  = ifaces(fi3)
                        c03  = m%f(kf3)%c0
                        ar3  = g%area(kf3)
                        nv3  = g%sf(:,kf3) / ar3
                        dp3  = norm2( g%xf(:,kf3) - g%xc(:,c03) )
                        eps3 = max( fld%porosity(c03), 1.0e-6_dp )
                        ! sqrt(K) along the interface normal (same reduction as
                        ! the internal bj_alpha branch)
                        lam3 = sqrt( max( dot_product( fld%perm_dir(:,c03), &
                                                       nv3**2 ), 0.0_dp ) )
                        beta3 = iface_beta
                        if ( beta3 <= 0.0_dp ) then
                           if ( iface_slip_alpha > 0.0_dp ) then
                              beta3 = eps3 * iface_slip_alpha
                           else if ( fld%bj_alpha(c03) > 0.0_dp ) then
                              beta3 = eps3 * fld%bj_alpha(c03)
                           else
                              beta3 = 1.0_dp        ! Ochoa-Tapia: O(1)
                           end if
                        end if
                        Vz = iface_zhang_velocity( &
                           iu(:,fi3), su(:,fi3), nv3, dp3, &
                           iface_df_ratio * dp3, ctrl%mu, ctrl%mu / eps3, &
                           eps3, lam3, beta3, iface_beta1, ctrl%rho )
                        ! Same cold-start ramp / under-relaxation as the legacy
                        ! exchange: the law output is a *fixed-point* update of
                        ! the interface state, so on a cold field it must be
                        ! phased in against the slab's own value -- applying the
                        ! full peer-driven value at iter 1 slams the interface
                        ! (measured 21 m/s on C_P_test) and the slab's thin
                        ! first cell, being strongly Dirichlet-pinned, then
                        ! amplifies it toward the free stream.
                        su_bc(:,fi3) = omega * Vz + (1.0_dp - omega) * iu(:,fi3)
                     end do
                  end block
               end if
               ! the uns side imposes only the peer face VELOCITY here;
               ! interface pressure is zero-gradient (bc_face_p for
               ! BC_INTERFACE returns pP), while the struct side imposes the
               ! uns back pressure via its subsonic-outflow characteristic
               ! ghost state.  Setting pressure as well over-constrains both
               ! sides and freezes a kPa-scale pressure jump.  The received
               ! sp (absolute Pa) is retained for diagnostics only.
               ! mass-conservation fix: with all-Dirichlet velocity
               ! boundaries (mass-flow-inlet + interface) the net boundary
               ! volume flux must vanish or the pure-Neumann PPE has no
               ! solution.  Sum the flux over all non-interface boundaries
               ! with the CURRENT bc state (includes the inlet ramp factor),
               ! then add a uniform normal correction to the interface face
               ! velocities so the total boundary flux is zero.
               block
                  integer  :: kf, gi, fi
                  logical  :: has_p_bc
                  real(dp) :: uf(3), F_other, F_iface, A_iface, du_n, nvec(3)
                  ! Uniform-normal flux correction is required only when every
                  ! non-interface boundary is velocity-Dirichlet (pure-Neumann
                  ! PPE).  With a pressure-specified boundary (pressure-outlet
                  ! / far-field) the PPE is well posed and the interface inflow
                  ! must NOT be cancelled to balance the not-yet-developed
                  ! outlet flux: that blocks the inlet at startup, the struct
                  ! side compresses against a closed wall, and the interface
                  ! pressure gradient explodes.
                  has_p_bc = any( bcs%gb(:)%btype == BC_POUTLET ) .or. &
                             any( bcs%gb(:)%btype == BC_FARFIELD )
                  ! phase 14: with iface_p_model /= grad0 the interface is a
                  ! PRESSURE-Dirichlet patch and its velocity is extrapolated
                  ! (bc_face_vel -> uP), exactly like pressure-outlet, so the
                  ! PPE closes continuity through the af*p' term and no uniform
                  ! flux correction is needed or wanted here: correcting an
                  ! extrapolated flux just fights the PPE (measured divergent
                  ! interface flux in C_P_test).
                  if ( iface_p_dir .and. n_uns_faces > 0 ) has_p_bc = .true.
                  F_other = 0.0_dp
                  do kf = 1, m%nfaces
                     if ( m%f(kf)%c1 /= 0 ) cycle
                     gi = bcs%fgrp(kf)
                     if ( gi == 0 ) cycle
                     if ( bcs%gb(gi)%btype == BC_INTERFACE ) cycle
                     call bc_face_vel( bcs, kf, g%xf(:,kf), g%sf(:,kf), &
                                       fld%u(:,m%f(kf)%c0), uf )
                     F_other = F_other + dot_product( uf, g%sf(:,kf) )
                  end do
                  F_iface = 0.0_dp
                  A_iface = 0.0_dp
                  do fi = 1, n_uns_faces
                     kf = ifaces(fi)
                     F_iface = F_iface + dot_product( su_bc(:,fi), g%sf(:,kf) )
                     A_iface = A_iface + g%area(kf)
                  end do
                  du_n = 0.0_dp
                  if ( .not. has_p_bc .and. A_iface > 0.0_dp ) then
                     du_n = -( F_other + F_iface ) / A_iface
                     do fi = 1, n_uns_faces
                        kf = ifaces(fi)
                        nvec = g%sf(:,kf) / g%area(kf)
                        su_bc(:,fi) = su_bc(:,fi) + du_n * nvec
                     end do
                  end if
                  if ( iter <= 8 .or. mod(iter,25)==0 ) then
                     write(*,'(a,i0,a,l1,a,3es11.3)') '[uns rank ', rank, &
                           '] iface has_p_bc=', has_p_bc, &
                           ' F_other/F_iface/du_n =', F_other, F_iface, du_n
                     write(*,'(a,i0,a,6es11.3)') '[uns rank ', rank, &
                           '] recv su x min/max, sT min/max, sp min/max =', &
                           minval(su(1,:)), maxval(su(1,:)), &
                           minval(sT), maxval(sT), minval(sp), maxval(sp)
                     write(*,'(a,i0,a,4es11.3)') '[uns rank ', rank, &
                           '] BC: su_bc x min/max =', &
                           minval(su_bc(1,:)), maxval(su_bc(1,:))
                  end if
               end block
               if ( mod(iter,100)==0 .or. iter==1 ) &
                  write(*,'(a,i0,a,i0,a,f6.3)') '[uns rank ', rank, &
                        '] iter ', iter, ' iface omega = ', omega
               ! ---- phase-14 interface pressure model --------------------------
               ! iface_p_model = dirichlet|momentum: the peer (struct) interface
               ! pressure is imposed per face as a DIRICHLET pressure.  The
               ! slab's PPE then stops being pure-Neumann: its gauge datum is
               ! fixed by the mainstream, the interface mass flux follows from
               ! continuity, and iface_p_anchor / the du_n correction are no
               ! longer needed (both bypassed above).  'momentum' adds the
               ! Betchen 2006 Eq.43/44 normal-momentum correction across the
               ! flow-area change, the two-sided (porous-side extrapolation)
               ! blend of his Sec.4.2 and the deferred p-mdot sub-iteration.
               if ( iface_p_dir ) then
                  block
                     type(reference_state_t) :: rr
                     integer  :: fi, kf, c0, ksub
                     real(dp) :: nhat(3), area_f, dp_g, eps_c, lam_c, bja_c
                     real(dp) :: mdot_f, mdot_it, un_f, af_f, p_fl, p_f, p_p
                     real(dp) :: p_cur, p_new, p_applied, p_old_face
                     real(dp) :: dp_corr
                     real(dp) :: p_sum, fl_sum, c_mean, c_max, gain_f, g_max
                     real(dp) :: om_p
                     real(dp), allocatable :: p_dir_bc(:), c_slip(:), u_peer(:,:)

                     allocate( p_dir_bc(n_uns_faces), c_slip(n_uns_faces), &
                               u_peer(3,n_uns_faces) )
                     rr     = get_ref_state()
                     p_sum  = 0.0_dp
                     dp_corr = 0.0_dp
                     fl_sum = 0.0_dp
                     c_mean = 0.0_dp
                     c_max  = 0.0_dp
                     g_max  = 0.0_dp

                     do fi = 1, n_uns_faces
                        kf     = ifaces(fi)
                        c0     = m%f(kf)%c0
                        area_f = g%area(kf)
                        nhat   = g%sf(:,kf) / area_f   ! outward (out of slab)
                        dp_g   = norm2( g%xf(:,kf) - g%xc(:,c0) )
                        eps_c  = fld%porosity(c0)
                        mdot_f = fld%flux(kf)          ! kg/s, + = out of slab
                        un_f   = mdot_f / ( ctrl%rho * area_f )

                        ! fluid-side estimate: peer interface pressure (gauge)
                        p_fl = sp(fi) - rr%p_ref
                        p_f  = p_fl
                        if ( iface_p_model == 2 ) &          ! Betchen Eq.43/44
                           p_f = iface_p_momentum( p_fl, mdot_f, un_f, &
                                                   area_f, eps_c )
                        ! porous-side deferred estimate: one-sided linear
                        ! extrapolation of the slab pressure to the face
                        ! (cf. mod_uns_fields:kink_face_pressure, which does the
                        ! same two-sidedly on internal fluid/porous faces)
                        p_p = fld%p(c0) + dp_g * &
                              dot_product( fld%gp(:,c0), nhat )
                        if ( iface_p_model == 2 ) then
                           p_cur = iface_p_blend( p_f, p_p, iface_p_w )
                        else
                           p_cur = p_f
                        end if

                        ! ---- p-mdot coupling (Betchen Sec.4.2) -----------------
                        ! Eq.43/44 makes the interface pressure depend on the
                        ! interface mass flow rate, which in turn depends on the
                        ! pressure: Betchen notes that "a small number of
                        ! iterations is required" and that the refined estimate
                        ! is "best implemented in a deferred fashion".
                        ! The local loop gain of the pair is
                        !   G = 2*(1-eps)*af*|mdot| / (eps*rho*A^2)
                        ! (mod_iface_law:iface_p_gain).  G is reported so the
                        ! stiffness of the pair is visible; when
                        ! iface_pm_subiter >= 1 the pair is then refined with an
                        ! unconditionally stable UNDER-RELAXED iteration
                        ! (omega = 1/(1+G) <= 1, the classic damped fixed point:
                        ! the undamped form oscillates with growth G per step and
                        ! is unusable when G >> 1).  Default 0 = pure deferred:
                        ! p_cur is Eq.44 evaluated on the flux the solver
                        ! actually produced, and the outer coupling loop closes
                        ! the p-mdot loop.
                        p_old_face = bcs%iface_p(kf)
                        if ( iter == iter0 + 1 ) &
                           p_old_face = fld%p(c0)   ! 1st iter: the interface was
                                                    ! zero-gradient, so the flux
                                                    ! in hand was produced under
                                                    ! the cell pressure
                        af_f = ctrl%rho * ( g%vol(c0) / fld%apc(c0) ) &
                               * area_f / dp_g
                        if ( iface_p_model == 2 ) then
                           gain_f = iface_p_gain( af_f, mdot_f, eps_c, &
                                                  ctrl%rho, area_f )
                           g_max  = max( g_max, gain_f )
                           ! The refined sub-iteration is admissible only in the
                           ! contracting regime G < 1.  For G > 1 the composed
                           ! map is monotone increasing with slope G, so its
                           ! fixed point REPELS: neither the raw iteration nor
                           ! any under-relaxation with omega in (0,1) can
                           ! converge (verified in src/coupling/test_iface_law
                           ! test 7).  Keep the deferred estimate there and let
                           ! the outer coupling iteration close the p-mdot loop.
                           ! Admissibility test for the local linearisation.
                           ! af is a SINGLE-CELL PPE response, not a steady
                           ! sensitivity of the converged interface flux (with a
                           ! mass-flow-inlet BC that flux is pinned by the inlet
                           ! and a uniform interface-pressure shift changes it by
                           ! ~zero).  The sub-iteration is therefore only run
                           ! when the implied flux change is small compared
                           ! with the physical flux; on C_P_test
                           ! af*(p_cur-p_prev) ~ 0.2 kg/s against mdot ~ 3.5e-3
                           ! kg/s, i.e. 55x too large -> skipped, and the
                           ! DEFERRED estimate (one Eq.44 evaluation on the
                           ! flux the solver produced) is used instead --
                           ! exactly the "implemented in a deferred fashion"
                           ! route of Betchen Sec.4.2.
                           if ( iface_pm_sub > 0 .and. gain_f < 1.0_dp .and. &
                                abs( af_f * ( p_cur - p_old_face ) ) &
                                < 0.1_dp * abs( mdot_f ) ) then
                              mdot_it = mdot_f
                              do ksub = 1, iface_pm_sub
                                 p_applied = p_cur
                                 mdot_it = iface_p_flux_response( mdot_it, af_f, &
                                               p_cur, p_old_face )
                                 un_f  = mdot_it / ( ctrl%rho * area_f )
                                 p_f   = iface_p_momentum( p_fl, mdot_it, un_f, &
                                                           area_f, eps_c )
                                 p_new = iface_p_blend( p_f, p_p, iface_p_w )
                                 om_p  = 1.0_dp / ( 1.0_dp + gain_f )
                                 p_cur = p_cur + om_p * ( p_new - p_cur )
                                    p_old_face = p_applied
                              end do
                           end if
                        end if

                        ! cold-start ramp / under-relaxation: reuse the interface
                        ! velocity omega (iface_relax, iface_ramp).  The
                        ! fallback (omega -> 0) is the CURRENT cell pressure,
                        ! i.e. the legacy zero-gradient treatment, so the
                        ! Dirichlet constraint is phased in gradually instead of
                        ! being applied at full strength on a cold field.
                        p_dir_bc(fi) = omega * p_cur &
                                     + (1.0_dp - omega) * fld%p(c0)

                        ! ---- iface_slip = bj: tangential stress-jump friction
                        ! (Beavers-Joseph / Ochoa-Tapia-Whitaker; the same
                        ! series-resistance law the solver applies to internal
                        ! fluid/porous faces via bj_alpha, see mod_iface_law and
                        ! cases/beavers_joseph).  The momentum assembly turns
                        ! this conductance into a tangential Robin term; only
                        ! the tangential part of the peer velocity is used.
                        if ( iface_slip == 1 ) then
                           lam_c = sqrt( max( dot_product( fld%perm_dir(:,c0), &
                                                           nhat**2 ), 0.0_dp ) )
                           if ( iface_slip_alpha > 0.0_dp ) then
                              bja_c = iface_slip_alpha
                           else
                              bja_c = fld%bj_alpha(c0)
                              if ( bja_c <= 0.0_dp ) bja_c = 1.0_dp
                           end if
                           ! The peer-side gap d_f is not part of the exchange
                           ! protocol; assume the peer resolves the interface
                           ! with a gap comparable to the slab's (d_f = d_p).
                           if ( iface_slip == 1 ) then
                              ! stress jump (Beavers-Joseph / Ochoa-Tapia)
                              c_slip(fi) = iface_slip_conductance( area_f, dp_g, &
                                               dp_g, lam_c, ctrl%mu, bja_c, eps_c )
                           else
                              ! tangential velocity continuity: flush (no-slip)
                              ! porous-side conductance C = mu_e*A/d_p, which is
                              ! the alpha -> infinity limit of the same law
                              c_slip(fi) = ( ctrl%mu / eps_c ) * area_f / dp_g
                           end if
                           u_peer(:,fi) = su(:,fi)
                           c_mean = c_mean + c_slip(fi)
                           c_max  = max( c_max, c_slip(fi) )
                        end if

                        dp_corr = max( dp_corr, abs( p_cur - fld%p(c0) ) )
                        p_sum  = p_sum + p_dir_bc(fi)
                        fl_sum = fl_sum + mdot_f
                     end do

                     c_mean = c_mean / real( max(n_uns_faces,1), dp )
                     call set_interface_p( bcs, ifaces, p_dir_bc )
                     if ( iface_slip == 1 ) &
                        call set_interface_slip( bcs, ifaces, c_slip, u_peer )

                     if ( iter <= 8 .or. mod(iter,25) == 0 ) then
                        write(*,'(a,i0,a,5es12.4)') '[uns rank ', rank, &
                           '] iface p-dir: p_bc(gauge) p_fl dp_corr mdot_sum area =', &
                           p_sum / real(max(n_uns_faces,1),dp), &
                           sum(sp)/real(max(n_uns_faces,1),dp) - rr%p_ref, &
                           dp_corr, fl_sum, sum( g%area(ifaces) )
                        if ( iface_p_model == 2 ) &
                           write(*,'(a,i0,a,2es12.4)') '[uns rank ', rank, &
                              '] iface p-mdot: G_max, 1/(1+G) =', g_max, &
                              1.0_dp/(1.0_dp+g_max)
                        if ( iface_slip == 1 ) &
                           write(*,'(a,i0,a,2es12.4)') '[uns rank ', rank, &
                              '] iface slip: mean/max C [kg/s] =', c_mean, c_max
                     end if
                     deallocate( p_dir_bc, c_slip, u_peer )
                  end block
               end if
               call set_interface_vel( bcs, ifaces, su_bc )
               ! ---- Zhang 2011 Eqs.27-29: per-phase interface temperature ----
               ! The peer's near-interface temperature T_flP (sT) plus the slab
               ! cell values (T_f, T_s) are solved simultaneously for the three
               ! interface temperatures: T_fl (the clear-fluid value = the
               ! volume average <T>^p), T_fi and T_si (the porous phases).
               ! Imposing T_fi on the fluid-phase equation and T_si on the
               ! solid-phase one makes the phase fluxes eps*F and (1-eps)*F,
               ! i.e. the Zhang Eq.28/29 porosity split, with no extra source
               ! term.  Under LTE the closure degenerates to the legacy single
               ! Dirichlet value, so nothing changes there.
               if ( iface_t_model == 1 .and. iface_t_ltne_is_on() .and. &
                    n_uns_faces > 0 ) then
                  block
                     integer  :: fi4, kf4, c04
                     real(dp) :: dp4, df4, e4, kf4v, kfe4, kse4
                     real(dp) :: Tfl4, Tfi4, Tsi4, F4
                     real(dp), allocatable :: Tf_bc(:), Ts_bc(:)
                     allocate( Tf_bc(n_uns_faces), Ts_bc(n_uns_faces) )
                     do fi4 = 1, n_uns_faces
                        kf4  = ifaces(fi4)
                        c04  = m%f(kf4)%c0
                        dp4  = norm2( g%xf(:,kf4) - g%xc(:,c04) )
                        df4  = iface_df_ratio * dp4
                        e4   = max( fld%porosity(c04), 1.0e-6_dp )
                        kf4v = max( ctrl%k_cond, 0.0_dp )
                        kfe4 = e4 * kf4v
                        kse4 = ( 1.0_dp - e4 ) * fld%k_s(c04)
                        call iface_zhang_T( sT(fi4), fld%T(c04), fld%T_s(c04), &
                                            e4, kf4v, kfe4, kse4, df4, dp4, &
                                            Tfl4, Tfi4, Tsi4, F4 )
                        Tf_bc(fi4) = Tfi4
                        Ts_bc(fi4) = Tsi4
                     end do
                     call set_interface_T_ltne( bcs, ifaces, Tf_bc, Ts_bc )
                     deallocate( Tf_bc, Ts_bc )
                  end block
               else
                  call set_interface_T( bcs, ifaces, sT )
               end if
               deallocate( su_bc )
            end if
         end if

         ! advance unstructured solver (phase-13 start-up schedule may decay
         ! n_uns_cur from n_uns_steps_start down to the steady n_uns_steps)
         n_uns_cur = scheduled_steps( iter, n_uns_steps_start, n_uns_steps, &
                                      step_decay )
         if ( n_uns_steps_start > n_uns_steps ) &
            write(*,'(a,i0,a,i0)') '[uns rank ', rank, &
                  '] coupling iter ', iter, ' -> uns steps: ', n_uns_cur
         call uns_solver_step( m, c, g, ctrl, bcs, fld, n_uns_cur, ier )
         if ( ier /= 0 ) then
            write(*,'(a,i0,a,i0)') '[uns rank ', rank, '] step FAILED ier=', ier
            call MPI_Abort( MPI_COMM_WORLD, ier, ierr )
            return
         end if

         ! extract interface state for next exchange
         call uns_solver_extract_iface( m, g, bcs, fld, ifaces, irho, iu, iT, ip, ier, ctrl%rho )
         ! ---- Zhang 2011 Eq.27: what the clear-fluid side must see ----------
         ! The fluid side may NOT be handed the porous FLUID-phase temperature:
         ! Eq.27 makes the interface temperature it sees the porosity-weighted
         ! VOLUME AVERAGE of the two porous phases,
         !    T_fl = <T>^p = eps <T_f>^f + (1-eps) <T_s>^s,
         ! with the flux split eps*F / (1-eps)*F solved from Eqs.28-29.  The
         ! same three-equation solve as on the BC side is therefore repeated
         ! here (same inputs: the peer value sT, the slab cell state, the local
         ! conductivities) and its T_fl replaces the extracted value.  Under
         ! LTE this is the identity T_f = T_s = T, i.e. bit-identical legacy.
         if ( iface_t_model == 1 .and. iface_t_ltne_is_on() .and. &
              size(ifaces) > 0 ) then
            block
               integer  :: fi5, kf5, c05
               real(dp) :: dp5, df5, e5, kf5v
               real(dp) :: Tfl5, Tfi5, Tsi5, F5
               do fi5 = 1, size(ifaces)
                  kf5  = ifaces(fi5)
                  c05  = m%f(kf5)%c0
                  dp5  = norm2( g%xf(:,kf5) - g%xc(:,c05) )
                  df5  = iface_df_ratio * dp5
                  e5   = max( fld%porosity(c05), 1.0e-6_dp )
                  kf5v = max( ctrl%k_cond, 0.0_dp )
                  call iface_zhang_T( sT(fi5), fld%T(c05), fld%T_s(c05), &
                                      e5, kf5v, e5*kf5v, &
                                      (1.0_dp-e5)*fld%k_s(c05), df5, dp5, &
                                      Tfl5, Tfi5, Tsi5, F5 )
                  iT(fi5) = Tfl5          ! Eq.27 value for the fluid side
               end do
            end block
         end if
         if ( size(ifaces) > 0 .and. (mod(iter,25)==0 .or. iter==1 .or. iter==n_couple) ) then
            write(*,'(a,i0,a,i0,a,3es12.4)') '[uns rank ', rank, '] diag iter ', iter, &
                  ' uns SI: mean u_x, mean p(gauge), |u|max =', &
                  sum(iu(1,:))/size(ifaces), sum(ip)/size(ifaces), &
                  maxval( sqrt( sum(iu*iu, dim=1) ) )
         end if

         ! No reflection: pass the raw uns cell-centre state back to the
         ! struct solver.  The struct side will compute a self-consistent
         ! ghost cell (ghost = 2*desired_face - inner) so that its own
         ! extract_iface arithmetic-mean returns exactly the received value.

         ! phase-10 auto-save: write the unstructured dump aligned with the
         ! struct side.  Only the uns root writes (global rank n_struct_ranks:
         ! 1 for a coupled run, 0 for a standalone nproc==1 run); the dump is
         ! a global field, serial-read on restart.
         if ( save_interval > 0 ) then
            if ( mod(iter, save_interval) == 0 ) then
               if ( rank == merge(0, 1, nproc == 1) ) then
#ifdef HAVE_MPI
                  call write_field_dump( uns_dump, m, fld, trim(casfile), ier )
                  if ( ier == 0 ) write(*,'(a,a,a,i0)') '[uns] saved dump: ', &
                        trim(uns_dump), ' at coupling iter ', iter
#endif
                  ! standalone run (no struct group): the uns root owns
                  ! couple_state.dat so a restart can recover the iter count.
                  if ( .not. has_struct ) call write_couple_state( iter )
               end if
            end if
         end if

         if (allocated(srho)) deallocate(srho, su, sT, sp)
      end do

      write(*,'(a,i0,a)') '[uns rank ', rank, '] coupling loop done'

      ! final field snapshot for archiving (ASCII VTU) + dump
      if ( .not. has_struct ) then
         write(*,'(a)') '[uns] standalone run: skipping coupled VTU snapshot'
      else
         block
            use mod_uns_output, only: vtk_write
            integer :: iov
            call vtk_write( m, c, g, fld, 'unMesh_coupled.vtu', iov )
            if ( iov == 0 ) write(*,'(a)') '[uns] wrote unMesh_coupled.vtu'
         end block
      end if
   end subroutine uns_group_driver

   !---------------------------------------------------------------------------
   ! Joint restart state helpers (phase 10).  couple_state.dat stores the last
   ! saved coupling iteration so that both sides restart from the same iter.
   ! Only rank 0 (struct root) reads/writes; the other side reads the file
   ! independently after a barrier-free filesystem sync.
   !---------------------------------------------------------------------------
   subroutine write_couple_state( iter )
      integer, intent(in) :: iter
      integer :: u, ios
      open( newunit=u, file='couple_state.dat', status='replace', action='write', &
            iostat=ios )
      if ( ios == 0 ) then
         write(u,*) iter
         close(u)
      end if
   end subroutine write_couple_state

   subroutine read_couple_state( iter )
      integer, intent(out) :: iter
      logical :: ex
      integer :: u, ios
      iter = 0
      inquire( file='couple_state.dat', exist=ex )
      if ( .not. ex ) return
      open( newunit=u, file='couple_state.dat', status='old', action='read', &
            iostat=ios )
      if ( ios == 0 ) then
         read(u,*,iostat=ios) iter
         close(u)
      end if
   end subroutine read_couple_state

   !===========================================================================
   ! convert_to_si -- array wrapper for struct_to_SI (non-dim -> SI).
   !===========================================================================
   subroutine convert_to_si(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                            rho, u, v, w, T, p, n)
      use mod_interface_units, only: struct_to_SI
      use mod_reference_state, only: get_ref_state, reference_state_t
      implicit none
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: rho_nd(n), u_nd(n), v_nd(n), w_nd(n), T_nd(n), p_nd(n)
      real(dp), intent(out) :: rho(n), u(n), v(n), w(n), T(n), p(n)
      integer :: i
      type(reference_state_t) :: r
      r = get_ref_state()
      do i = 1, n
         call struct_to_SI(rho_nd(i), u_nd(i), v_nd(i), w_nd(i), T_nd(i), p_nd(i), &
                           rho(i), u(i), v(i), w(i), T(i), p(i), r)
      end do
   end subroutine convert_to_si

   !===========================================================================
   ! convert_uns_to_struct_nd -- map the uns SI interface state onto the
   ! struct interface faces, then SI -> struct non-dim.
   !
   ! Two mappings are supported:
   !   * n_struct == n_uns: direct 1:1 mapping in face-list order.  Both
   !     patches are generated as matching (j outer, k inner) quad arrays,
   !     so no geometric interpolation is needed in this case.
   !   * otherwise: uniform area-weighted average as a placeholder until
   !     the full uns->struct area-weighted interpolation
   !     (mod_interface_exchange) is wired into the driver.
   ! NB: the uns pressure is GAUGE; the struct non-dim pressure is absolute
   ! (p_inf* = 1/(gamma*Ma^2)), so p_ref is added before conversion.
   !===========================================================================
   subroutine convert_uns_to_struct_nd(urho, uu, uT, up, n_uns, &
                                       rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                       n_struct, r)
      use mod_interface_units, only: SI_to_struct
      use mod_reference_state, only: reference_state_t
      implicit none
      integer,  intent(in)  :: n_uns, n_struct
      real(dp), intent(in)  :: urho(n_uns), uu(3,n_uns), uT(n_uns), up(n_uns)
      real(dp), intent(out) :: rho_nd(n_struct), u_nd(n_struct), v_nd(n_struct), &
                               w_nd(n_struct), T_nd(n_struct), p_nd(n_struct)
      type(reference_state_t), intent(in) :: r

      integer :: i, j
      real(dp) :: avg_rho, avg_u(3), avg_T, avg_p
      real(dp) :: rho_s, u_s, v_s, w_s, T_s, p_s, p_abs

      if ( n_struct == n_uns ) then
         ! direct 1:1 mapping (matched patches)
         do i = 1, n_struct
            p_abs = up(i) + r%p_ref
            call SI_to_struct(urho(i), uu(1,i), uu(2,i), uu(3,i), uT(i), p_abs, &
                              rho_s, u_s, v_s, w_s, T_s, p_s, r)
            rho_nd(i) = rho_s
            u_nd(i)   = u_s
            v_nd(i)   = v_s
            w_nd(i)   = w_s
            T_nd(i)   = T_s
            p_nd(i)   = p_s
         end do
      else
         ! uniform area-weighted average (equal weights for now)
         avg_rho = sum(urho) / real(n_uns, dp)
         avg_u   = sum(uu, dim=2) / real(n_uns, dp)
         avg_T   = sum(uT) / real(n_uns, dp)
         avg_p   = sum(up) / real(n_uns, dp) + r%p_ref

         do i = 1, n_struct
            call SI_to_struct(avg_rho, avg_u(1), avg_u(2), avg_u(3), avg_T, avg_p, &
                              rho_s, u_s, v_s, w_s, T_s, p_s, r)
            rho_nd(i) = rho_s
            u_nd(i)   = u_s
            v_nd(i)   = v_s
            w_nd(i)   = w_s
            T_nd(i)   = T_s
            p_nd(i)   = p_s
         end do
      end if
   end subroutine convert_uns_to_struct_nd

end module mod_mix_driver