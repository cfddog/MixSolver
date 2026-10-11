$control_ec
! C_P_test -- structured (compressible) side, transpiration-cooling plate.
! Domain x in [0,1], y in [0,1] (blocks 2-1-3 along x); interface at the
! block-1 j=1 face (y=0, x in [0.2,0.8]) via gridgen generic:8.
  Ma=0.3d0
  Re=5.0d4                ! based on Ref_L (grid unit)
  gamma=1.4d0
  PrL=0.7d0
  PrT=0.9d0
  AoA=0.d0
  AoS=0.d0
  P_outlet=-1.d0          ! <0: extrapolation outlet
  t_end=1.0d12            ! steady local-dt march
  Iflag_local_dt=1
  dt_global=0.01d0
  Time_Method=0           ! LU-SGS
  CFL=5.d0
  dtmax=10.d0
  dtmin=1.d-9
  Iflag_init=0
  If_viscous=1
  Iflag_turbulence_model=0
  Iflag_Scheme=5          ! MUSCL3
  Iflag_Flux=5            ! Van_Leer
  IFlag_Reconstruction=0
  Bound_Scheme=4
  Kstep_save=5000
  Kstep_show=1
  Kstep_average=0
  Kstep_smooth=-1
  Kstep_init_smooth=0
  If_Residual_smoothing=0
  If_dtime_mesh=1
  w_LU=1.d0
  Mesh_File_Format=0
  Num_Mesh=1
  T_inf=309.d0            ! mainstream (inflow/farfield) temperature [K]
  Twall=-1.d0             ! adiabatic walls
  Kt_inf=1.d-5
  Wt_inf=1.d-2
  Step_Inner_Limit=20
  MUT_MAX=-1.d0
  IF_Debug=0
  Ref_S=1.d0
  Ref_L=1.d0              ! grid unit = 1 (mesh and uns mesh share units)
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
