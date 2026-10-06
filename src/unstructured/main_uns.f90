!===============================================================================
! main.f90 -- UNSSolver driver
!
! Phase 1: read a Fluent ASCII .cas file, build connectivity and geometric
!          data, print a verification report (zone table, volumes, GCL).
! Phase 3: read the control file, apply boundary conditions and run the
!          steady incompressible SIMPLE solver.
!
! Usage: unsolver [mesh.cas] [control]     (defaults: Grid_BC/cavity.cas,
!                                           cases/cavity.control)
!===============================================================================
program unsolver
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   use mod_uns_cas_reader
   use mod_uns_connectivity
   use mod_uns_geometry
   use mod_uns_control
   use mod_uns_bc
   use mod_uns_fields
   use mod_uns_simple
   use mod_uns_output
   implicit none

   type(mesh_t)   :: m
   type(conn_t)   :: c
   type(geom_t)   :: g
   type(ctrl_t)   :: ctrl
   type(bc_t)     :: bcs
   type(fields_t) :: fld
   character(len=512) :: casfile, ctlfile, arg, vtufile, outstem
   character(len=8)   :: outext
   logical            :: out_user
   integer        :: ier, i, j, k, zid, nint, nquad, ntri, npoly, ct
   integer        :: nb_int, nb_bnd, bad_self
   real(dp)       :: vol_sum, vol_min, vol_max, gcl_max
   real(dp)       :: bbox_lo(3), bbox_hi(3), dV, umax, mscale
   integer        :: type_cnt(0:7)
   integer        :: idbg, idf   ! TEMP DEBUG face flux dump

   ! ---- command line -------------------------------------------------------
   casfile = 'Grid_BC/cavity.cas'
   ctlfile = 'cases/cavity.control'
   if ( command_argument_count() >= 1 ) then
      call get_command_argument( 1, arg )
      casfile = trim(arg)
   end if
   if ( command_argument_count() >= 2 ) then
      call get_command_argument( 2, arg )
      ctlfile = trim(arg)
   end if

   ! output file stem: derived from the mesh basename unless given as arg 3.
   ! The extension (.vtu / .plt) is reconciled with ctrl%out_format after the
   ! control file is read, since the format choice is not known until then.
   outstem = 'result'
   j = scan( casfile, '/', back=.true. )
   k = scan( casfile(j+1:), '.', back=.true. )
   if ( k > 1 ) outstem = casfile(j+1:j+k-1)
   out_user = .false.
   vtufile = trim(outstem)//'.vtu'
   if ( command_argument_count() >= 3 ) then
      call get_command_argument( 3, arg )
      vtufile = trim(arg)
      out_user = .true.
   end if

   write(*,'(a)') '=============================================================='
   write(*,'(a)') '  UNSSolver -- unstructured multiphysics CFD (pre + SIMPLE)'
   write(*,'(a)') '=============================================================='
   write(*,'(a,a)') '  Mesh file    : ', trim(casfile)
   write(*,'(a,a)') '  Control file : ', trim(ctlfile)
   ! NB: the output-file line is printed after read_control reconciles the
   ! extension with ctrl%out_format (see below).

   ! ---- read mesh ------------------------------------------------------------
   call read_cas( trim(casfile), m, ier )
   if ( ier /= 0 ) then
      write(*,'(a,i0)') 'FATAL: mesh reading failed, ier = ', ier
      stop 1
   end if

   ! optional unit scaling (e.g. mm -> m) before connectivity/geometry
   call read_mesh_scale( trim(ctlfile), mscale )
   if ( mscale /= 1.0_dp ) then
      m%x = m%x * mscale
      write(*,'(a,es10.3)') '  mesh coordinates scaled by ', mscale
   end if

   write(*,'(a)') ''
   write(*,'(a)') '--- Mesh summary ---'
   write(*,'(a,i10)') '  dimension       : ', m%ndim
   write(*,'(a,i10)') '  nodes           : ', m%nnodes
   write(*,'(a,i10)') '  faces           : ', m%nfaces
   write(*,'(a,i10)') '  cells           : ', m%ncells

   ! cell type histogram
   type_cnt = 0
   do k = 1, m%ncells
      ct = 0
      if ( m%ctype(k) >= 0 .and. m%ctype(k) <= 7 ) ct = m%ctype(k)
      type_cnt(ct) = type_cnt(ct) + 1
   end do
   write(*,'(a)') '  cell types      : ' // &
      ' mixed=  tri=  tet=  quad=  hex=  pyr=  wedge=  poly='
   write(*,'(a,8(i8,2x))') '                    ', type_cnt(0:7)

   ! ---- connectivity ----------------------------------------------------------
   call build_connectivity( m, c, ier )
   if ( ier /= 0 ) then
      write(*,'(a,i0)') 'FATAL: connectivity failed, ier = ', ier
      stop 1
   end if

   ! topology sanity checks
   nint = m%nfaces - c%nbf
   nb_int = 0; nb_bnd = 0; bad_self = 0
   do i = 1, m%nfaces
      if ( m%f(i)%c1 == 0 ) then
         nb_bnd = nb_bnd + 1
      else
         nb_int = nb_int + 1
         if ( m%f(i)%c0 == m%f(i)%c1 ) bad_self = bad_self + 1
      end if
   end do
   write(*,'(a)') ''
   write(*,'(a)') '--- Topology checks ---'
   write(*,'(a,i10)') '  interior faces  : ', nb_int
   write(*,'(a,i10)') '  boundary faces  : ', nb_bnd
   write(*,'(a,i10)') '  self-connected  : ', bad_self
   write(*,'(a,i10)') '  max faces/cell  : ', maxval( c%cf_ptr(2:m%ncells+1) - c%cf_ptr(1:m%ncells) )

   ! ---- geometry ---------------------------------------------------------------
   call compute_geometry( m, g, ier )
   if ( ier /= 0 ) then
      write(*,'(a,i0)') 'FATAL: geometry failed, ier = ', ier
      stop 1
   end if
   call geom_stats( m, g, vol_sum, vol_min, vol_max, gcl_max, ier )

   ! ---- register CAS interface zones for coupling (phase 3 step B) -----------
   write(*,'(a)') ''
   write(*,'(a)') '--- Coupling interface registration ---'
   call register_interface_zones( m, g )

   ! bounding box from nodes
   bbox_lo = minval( m%x, dim=2 )
   bbox_hi = maxval( m%x, dim=2 )
   dV = (bbox_hi(1)-bbox_lo(1)) * (bbox_hi(2)-bbox_lo(2)) * &
        (bbox_hi(3)-bbox_lo(3))

   write(*,'(a)') ''
   write(*,'(a)') '--- Geometry checks ---'
   write(*,'(a,3(es12.5,1x))') '  bbox min        : ', bbox_lo
   write(*,'(a,3(es12.5,1x))') '  bbox max        : ', bbox_hi
   write(*,'(a,es14.6)') '  sum(vol)        : ', vol_sum
   write(*,'(a,es14.6)') '  bbox volume     : ', dV
   write(*,'(a,es14.6)') '  vol min         : ', vol_min
   write(*,'(a,es14.6)') '  vol max         : ', vol_max
   write(*,'(a,es14.6)') '  GCL max residual: ', gcl_max

   ! ---- zone table ----------------------------------------------------------------
   write(*,'(a)') ''
   write(*,'(a)') '--- Face zone table ---'
   write(*,'(a)') '   id  condition          user name           nf      tri    quad   polygon'
   do i = 1, m%nzone
      ntri = 0; nquad = 0; npoly = 0
      do k = 1, m%nfaces
         if ( m%f(k)%zone == m%zone(i)%id ) then
            select case ( m%f(k)%nn )
            case ( 2 ); ntri = ntri          ! line (should not occur in 3D)
            case ( 3 ); ntri = ntri + 1
            case ( 4 ); nquad = nquad + 1
            case default; npoly = npoly + 1
            end select
         end if
      end do
      zid = m%zone(i)%id
      write(*,'(1x,i5,2x,a18,2x,a18,2x,i8,2x,i6,2x,i6,2x,i6)') &
         zid, m%zone(i)%cond_name, m%zone(i)%user_name, m%zone(i)%nf, &
         ntri, nquad, npoly
   end do

   ! ---- overall verdict ---------------------------------------------------------
   write(*,'(a)') ''
   if ( ier /= 0 .or. bad_self > 0 .or. gcl_max >= 1.0e-10_dp ) then
      write(*,'(a)') 'RESULT: FAIL -- see messages above.'
      stop 1
   end if
   write(*,'(a)') 'RESULT: PASS -- mesh, connectivity and geometry verified.'

   ! ---- phase 3: boundary conditions + SIMPLE solver -----------------------------
   call read_control( trim(ctlfile), ctrl, ier )
   if ( ier /= 0 ) then
      write(*,'(a,i0)') 'FATAL: control file reading failed, ier = ', ier
      stop 1
   end if

   ! reconcile the output-file extension with the requested format unless the
   ! user explicitly supplied the output path as command-line argument 3
   if ( .not. out_user ) then
      if ( ctrl%out_format == 'tecplot' ) then
         outext = '.plt'
      else
         outext = '.vtu'
      end if
      vtufile = trim(outstem)//trim(outext)
   end if
   write(*,'(a,a)') '  Output file  : ', trim(vtufile)

   ! ---- phase 3 step C: resolve cell (volume) zones from the control file ----
   call resolve_cell_zones( m, ctrl, ier )
   if ( ier /= 0 ) then
      write(*,'(a,i0)') 'FATAL: cell-zone resolution failed, ier = ', ier
      stop 1
   end if

   call build_bc( m, ctrl, bcs, ier )
   if ( ier /= 0 ) then
      write(*,'(a,i0)') 'FATAL: boundary condition setup failed, ier = ', ier
      stop 1
   end if

   call init_fields( m, g, fld )

   ! populate per-cell porous-medium coefficients from ctrl%cz
   call setup_porous_fields( m, ctrl, fld )

   ! ---- initial temperature field (uniform + optional x-asymmetric mode) ----
   ! Used for natural-convection (Rayleigh-Benard) cases where the field must
   ! start near the wall mean and a small perturbation breaks the symmetry of
   ! the pure-conduction state.
   if ( ctrl%init_T /= 0.0_dp .or. ctrl%t_pert /= 0.0_dp ) then
      block
         real(dp) :: xmin, xmax, lx, pi_v
         integer  :: ic
         pi_v = acos( -1.0_dp )
         xmin = minval( m%x(1,:) )
         xmax = maxval( m%x(1,:) )
         lx   = max( xmax - xmin, epsilon(1.0_dp) )
         do ic = 1, m%ncells
            fld%T(ic) = ctrl%init_T + ctrl%t_pert &
                        * sin( pi_v * ( g%xc(1,ic) - xmin ) / lx )
            fld%T_s(ic) = fld%T(ic)
         end do
         fld%T_old      = fld%T
         fld%T_old_old  = fld%T
         fld%T_s_old    = fld%T_s
         fld%T_s_old_old = fld%T_s
      end block
   end if

   ! ---- dispatch: steady SIMPLE, transient PISO or PIMPLE -------------------
   if ( ctrl%transient ) then
      if ( ctrl%pimple ) then
         call pimple_run( m, c, g, ctrl, bcs, fld, trim(vtufile), ier )
         if ( ier /= 0 ) then
            write(*,'(a,i0)') 'FATAL: PIMPLE solver failed, ier = ', ier
            stop 1
         end if
         write(*,'(a)') 'RESULT: PIMPLE run finished.'
      else
         ! PISO writes per-step VTU files as <vtufile-stem>_stepNNNN.vtu
         ! when n_out_every > 0; the final state is also written to <vtufile>.
         call piso_run( m, c, g, ctrl, bcs, fld, trim(vtufile), ier )
         if ( ier /= 0 ) then
            write(*,'(a,i0)') 'FATAL: PISO solver failed, ier = ', ier
            stop 1
         end if
         write(*,'(a)') 'RESULT: PISO run finished.'
      end if
   else
      call simple_run( m, c, g, ctrl, bcs, fld, ier )
      if ( ier /= 0 ) then
         write(*,'(a,i0)') 'FATAL: SIMPLE solver failed, ier = ', ier
         stop 1
      end if
      write(*,'(a)') 'RESULT: SIMPLE run finished.'
   end if

   ! ---- solution summary ----------------------------------------------------------
   umax = maxval( norm2( fld%u, dim=1 ) )
   k = minloc( norm2( g%xc - 0.5_dp, dim=1 ), dim=1 )
   write(*,'(a)') ''
   write(*,'(a)') '--- Solution summary ---'
   write(*,'(a,es12.5)') '  max |u|         : ', umax
   write(*,'(a,3(es11.4,1x),a,es12.5)') '  cell nearest (0.5,0.5,0.5): ', &
      g%xc(:,k), ' u = ', norm2( fld%u(:,k) )
   write(*,'(a,3(es11.4,1x))') '                        u_vec = ', fld%u(:,k)
   write(*,'(a,es14.6)') '  max p           : ', maxval( fld%p )
   write(*,'(a,es14.6)') '  min p           : ', minval( fld%p )

   ! ---- wall heat balance / Nusselt number (natural convection) -------------
   if ( ctrl%boussinesq ) call thermal_wall_report( m, g, ctrl, bcs, fld )

   ! TEMP DEBUG: dump face fluxes for transverse-mode diagnosis (remove later)
   open( newunit = idbg, file = 'flux_debug.txt', status = 'replace' )
   do idf = 1, m%nfaces
      write(idbg,'(i8,1x,i8,1x,i8,1x,4(es14.6,1x))') &
         idf, m%f(idf)%c0, m%f(idf)%c1, g%sf(:,idf), fld%flux(idf)
   end do
   close(idbg)

   ! ---- output: final result file + benchmark comparison -----------------------
   ! Format selected by ctrl%out_format ('vtu' default, or 'tecplot' PLT).
   ! PISO has already written intermediate snapshots as <stem>_stepNNNN.<ext>.
   if ( ctrl%out_format == 'tecplot' ) then
      call tecplot_write( m, c, g, ctrl, fld, trim(vtufile), ier )
   else
      call vtk_write( m, c, g, fld, trim(vtufile), ier )
   end if
   if ( ier /= 0 ) then
      write(*,'(a,i0)') 'FATAL: result output failed, ier = ', ier
      stop 1
   end if
   call ghia_compare( m, g, fld, ctrl%rho / max( ctrl%mu, tiny(1.0_dp) ), ier )

end program unsolver
