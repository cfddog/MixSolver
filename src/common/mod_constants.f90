!===============================================================================
! mod_constants.f90 -- Physical and numerical constants for MixNSSolver
!
! This module defines:
!   (1) Mathematical constants shared by both solvers.
!   (2) Thermodynamic / gas-dynamic constants (ideal-gas assumptions).
!   (3) Reference-state parameters used for non-dimensionalisation and
!       for converting between the structured solver's internal non-dimensional
!       variables and SI units at the coupling interface.
!   (4) Small numbers to guard against division by zero.
!
! Solver-specific boundary-condition type codes and scheme identifiers are
! NOT placed here; they remain in each solver's own module (mod_struct_bc,
! mod_uns_bc, etc.) to avoid polluting the common layer with solver details.
!===============================================================================
module mod_constants
   use mod_precision, only: dp
   implicit none
   private

   ! ---------------------------------------------------------------------------
   ! Mathematical constants
   ! ---------------------------------------------------------------------------
   !> Threshold below which a quantity is treated as zero
   real(dp), parameter, public :: LIM_ZERO = 1.0e-20_dp

   ! ---------------------------------------------------------------------------
   ! Thermodynamic constants (ideal gas, calorically perfect)
   ! ---------------------------------------------------------------------------
   !> Ratio of specific heats  cp/cv  for air
   real(dp), parameter, public :: GAMMA = 1.4_dp

   !> Specific gas constant for air  R = R_universal / M_air  [J/(kg K)]
   real(dp), parameter, public :: R_GAS = 287.058_dp

   !> Specific heat at constant pressure  cp = gamma*R/(gamma-1)  [J/(kg K)]
   real(dp), parameter, public :: CP_AIR = GAMMA * R_GAS / (GAMMA - 1.0_dp)

   !> Specific heat at constant volume  cv = R/(gamma-1)  [J/(kg K)]
   real(dp), parameter, public :: CV_AIR = R_GAS / (GAMMA - 1.0_dp)

   !> Prandtl number for air (molecular)
   real(dp), parameter, public :: PR_AIR = 0.72_dp

   ! ---------------------------------------------------------------------------
   ! Reference-state parameters for non-dimensionalisation
   ! ---------------------------------------------------------------------------
   ! These reference values define the mapping between the structured solver's
   ! non-dimensional variables and SI units.  The unstructured solver works
   ! entirely in SI, so these are only needed at the coupling interface.
   !
   ! Non-dimensionalisation convention (OpenCFD-EC style):
   !   rho* = rho / rho_ref
   !   u*   = u   / a_ref        (velocity scaled by reference speed of sound)
   !   T*   = T   / T_ref
   !   p*   = p   / (rho_ref * a_ref^2)
   !   e*   = e   / (rho_ref * a_ref^2)    (total energy per unit volume)
   !
   ! The reference speed of sound is  a_ref = sqrt(gamma * R_GAS * T_ref).
   !
   ! Default values correspond to standard sea-level conditions so that the
   ! solver can run without a configuration file; the user should override
   ! them via the coupling control file for each case.
   ! ---------------------------------------------------------------------------

   !> Reference temperature [K]
   real(dp), parameter, public :: T_REF  = 288.15_dp

   !> Reference density [kg/m^3]
   real(dp), parameter, public :: RHO_REF = 1.225_dp

   !> Reference speed of sound [m/s]  = sqrt(gamma * R_GAS * T_REF)
   real(dp), parameter, public :: A_REF  = 340.29702908_dp   ! sqrt(1.4*287.058*288.15)

   !> Reference velocity [m/s]  (free-stream, for convective scaling)
   real(dp), parameter, public :: U_REF  = A_REF        ! default: sonic

   !> Reference pressure [Pa]  = rho_ref * R_GAS * T_REF
   real(dp), parameter, public :: P_REF  = RHO_REF * R_GAS * T_REF

   !> Reference length [m]  (user-defined per case; default 1 m)
   real(dp), parameter, public :: L_REF  = 1.0_dp

   ! ---------------------------------------------------------------------------
   ! Derived reference quantities
   ! ---------------------------------------------------------------------------
   !> Reference dynamic pressure  q_ref = 0.5 * rho_ref * U_REF^2  [Pa]
   real(dp), parameter, public :: Q_REF = 0.5_dp * RHO_REF * U_REF * U_REF

   ! ---------------------------------------------------------------------------
   ! Unit conversion helpers (compile-time constants)
   ! ---------------------------------------------------------------------------
   ! Conversion factors for the structured solver <-> SI interface.
   ! These are the denominators used when converting FROM SI TO non-dimensional:
   !   rho_nd = rho_SI / RHO_REF
   !   u_nd   = u_SI   / A_REF
   !   T_nd   = T_SI   / T_REF
   !   p_nd   = p_SI   / (RHO_REF * A_REF**2)
   !
   ! And the reciprocals for converting FROM non-dimensional TO SI:
   !   rho_SI = rho_nd * RHO_REF
   !   u_SI   = u_nd   * A_REF
   !   T_SI   = T_nd   * T_REF
   !   p_SI   = p_nd   * (RHO_REF * A_REF**2)

   !> p_ref_dynamic = rho_ref * a_ref^2  (pressure scaling for struct solver)
   real(dp), parameter, public :: P_SCALE = RHO_REF * A_REF * A_REF

   ! ---------------------------------------------------------------------------
   ! Standard atmosphere (for convenience / initialisation)
   ! ---------------------------------------------------------------------------
   !> Standard sea-level temperature [K]
   real(dp), parameter, public :: T_SL  = 288.15_dp

   !> Standard sea-level pressure [Pa]
   real(dp), parameter, public :: P_SL  = 101325.0_dp

   !> Standard sea-level density [kg/m^3]
   real(dp), parameter, public :: RHO_SL = P_SL / (R_GAS * T_SL)

end module mod_constants
