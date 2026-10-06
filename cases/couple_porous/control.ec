$control_ec
! couple_porous -- structured side (OpenCFD-EC), low-Ma nearly-incompressible
! strictly 1D through-flow (all four lateral faces symmetry), i- forced
! uniform inlet, i+ coupling interface (generic:8).  Matched to the
! unstructured SI side:
!   T=300 K, rho=1.177 kg/m^3, Ma=0.1 -> u=34.72 m/s
!   physical air mu=1.846e-5 Pa.s -> Re=2214 per mm (U_inf based)
  Ma=0.1d0
  Re=2214.d0              ! U_inf based (rho*U_inf*1mm/mu)
  gamma=1.4d0
  PrL=0.7d0
  PrT=0.9d0
  AoA=0.d0
  AoS=0.d0
  P_outlet=-1.d0
  t_end=1.0d12            ! local-dt steady march; do not stop by physical time
  Iflag_local_dt=1
  dt_global=0.01d0
  Time_Method=0           ! LU-SGS
  CFL=5.d0
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
  T_inf=300.d0
  Twall=-1.d0             ! adiabatic (uniform-T field)
  Kt_inf=1.d-5
  Wt_inf=1.d-2
  Step_Inner_Limit=20
  MUT_MAX=-1.d0
  IF_Debug=0
  Ref_S=1.d0
  Ref_L=1.d0              ! grid unit 1 mm
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
