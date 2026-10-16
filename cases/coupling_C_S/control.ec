$control_ec
! coupling_C_S -- structured (compressible) side: Ma=3 laminar flat plate
! convecting over a 100 mm steel slab (conjugate heat transfer, phase 17).
! Mesh coupling_C_S/Mesh3d.x (mm units, 81 x 101 x 2 nodes, Plot3D
! unformatted): x in [0,1000], y in [0,500], z in [0,100]; the first wall
! layer is 1 um.  The unstructured peer (unMesh.cas, same mm units) is the
! SOLID slab x in [0,1000], y in [-100,0], z in [0,100], so the two domains
! butt at y=0 and the shared face (80 quads, 12.5 mm x 100 mm each) is
! registered as gridgen generic:8 (bc3d.inp line 5, bc3d.fvbnd code 8).
!
! Boundary conditions (bc3d.inp / bc3d.fvbnd, Gridgen codes):
!   i=1   (x=0)     code 5 = Inflow  -> boundary_Farfield(FLAG_INLET) forces
!                    the ghost state to the free stream (all four
!                    characteristics enter): a SUPERSONIC inlet, as the case
!                    readme asks for ("来流采用超声速入口").
!   i=81  (x=1000)  code 6 = Outflow -> boundary_Farfield(FLAG_OUTLET) is
!                    pure interior extrapolation: a SUPERSONIC outlet.
!   k=1,2 (z faces) code 3 = Symmetry
!   j=1   (y=0)     code 8 = coupling interface -> isothermal NO-SLIP wall
!                    whose temperature is the solid surface temperature and
!                    which returns its wall heat flux to the solid side.
!   j=101 (y=500)   code 6 = Outflow (supersonic; free stream at the top)
! NOTE there is no ordinary wall in this mesh: the ONLY wall is the coupling
! interface, which is exactly what a CHT plate needs.
!
! Reference state (mix.control, SI): rho_ref = 0.44124 kg/m^3,
! T_ref = 800 K, L_ref = 1e-3 m (one grid unit = 1 mm),
! u_ref = Ma*a_ref = 3*sqrt(1.4*287.058*800) = 1700.95 m/s.
! OpenCFD-EC nondimensionalisation: u* = u/U_inf (= 1 in the free stream),
! p* = p/(rho_ref u_ref^2) with p_inf* = 1/(gamma Ma^2) so that
! p_inf = rho_ref R T_ref = 101331 Pa (~1 atm), T* = T/T_ref, mu* = 1/Re.
! The coupling layer converts back with mu_SI = rho_ref*u_ref*L_ref/Re,
! i.e. Re = 1000 per grid unit (1 mm) <-> mu = 7.506e-4 Pa.s (air at 800 K
! would be ~3.7e-5, so this is a viscosity-enlarged model fluid: the case
! readme prescribes Ma = 3, Re = 1000 and we keep its Re per grid unit).
  Ma=3.0d0
  Re=1000.d0              ! rho_ref*u_ref*L_ref/mu (L_ref = 1 mm); mu=7.506e-4 Pa.s
  gamma=1.4d0
  PrL=0.7d0
  PrT=0.9d0
  AoA=0.d0
  AoS=0.d0
  P_outlet=-1.d0          ! <0: extrapolation; unused (in/outflow are forced)
  t_end=1.0d12            ! steady local-dt march; do not stop by physical time
  Iflag_local_dt=1
  dt_global=0.01d0
  Time_Method=0           ! LU-SGS
  CFL=0.5d0               ! case readme: "cfd=0.5"
  dtmax=10.d0
  dtmin=1.d-9
  Iflag_init=0            ! uniform free-stream initialisation
  If_viscous=1
  Iflag_turbulence_model=0
  Iflag_Scheme=5          ! MUSCL3
  Iflag_Flux=5            ! Van_Leer
  IFlag_Reconstruction=0
  Bound_Scheme=4
  Kstep_save=5000
  Kstep_show=100
  Kstep_average=0
  Kstep_smooth=-1
  Kstep_init_smooth=0
  If_Residual_smoothing=0
  If_dtime_mesh=1
  w_LU=1.d0
  Mesh_File_Format=0
  Num_Mesh=1
  T_inf=800.d0            ! free-stream static temperature [K] (Ma=3 -> T_aw~2079 K)
  Twall=-1.d0             ! no ordinary wall in this mesh; interface is a wall
  Kt_inf=1.d-5
  Wt_inf=1.d-2
  Step_Inner_Limit=20
  MUT_MAX=-1.d0
  IF_Debug=0
  Ref_S=1.d0
  Ref_L=1.d0              ! grid unit 1 mm (physical length comes from mix.control)
  Centroid=0.d0, 0.d0, 0.d0
  Cood_Y_UP=1
  Periodic_dX=0.d0
  Periodic_dY=0.d0
  Periodic_dZ=0.d0
  IFLAG_LIMIT_FLOW=1
  Ldmin=1.d-6
  Ldmax=1000.d0
  Lpmin=1.d-6
  Lpmax=1000.d0
  Lumax=1000.d0
  CP1_NSA=0.2d0
  CP2_NSA=100.d0
  IF_Scheme_Positivity=1
  Iflag_savefile=0
  NUM_THREADS=1
/
