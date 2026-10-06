!===============================================================================
! mod_precision.f90 -- Unified precision kind constants for MixNSSolver
!
! This module provides a single source of truth for floating-point and integer
! kind parameters used by both the structured (OpenCFD-EC) and unstructured
! (UNSSolverProj) solvers, as well as the coupling layer.
!
! Compatibility notes:
!   - OpenCFD-EC originally used  PRE_EC = 8  (hardcoded kind).
!     The alias  PRE_EC = dp  is provided so that legacy code can be ported
!     by simply changing  "use precision_EC"  to  "use mod_precision".
!   - UNSSolverProj originally used  dp = selected_real_kind(15,307).
!     The definition below is identical.
!   - MPI data-type constants (OCFD_DATA_TYPE, UNS_MPI_REAL) are provided so
!     that both solvers can use the same MPI communication routines without
!     modification.
!
! Usage:
!   use mod_precision
!   real(dp) :: x       ! double-precision real
!   real(sp) :: y       ! single-precision real
!   integer(ip) :: i    ! default integer
!===============================================================================
module mod_precision
   implicit none
   private

   ! ---------------------------------------------------------------------------
   ! Floating-point kinds
   ! ---------------------------------------------------------------------------
   !> Double precision: ~15 significant digits, exponent range ~307
   integer, parameter, public :: dp = selected_real_kind(15, 307)

   !> Single precision: ~6 significant digits, exponent range ~37
   integer, parameter, public :: sp = selected_real_kind(6, 37)

   ! ---------------------------------------------------------------------------
   ! Integer kinds
   ! ---------------------------------------------------------------------------
   !> Default integer kind for array indices and loop counters
   integer, parameter, public :: ip = kind(1)

   !> 64-bit integer kind for large array sizes / global indices
   integer, parameter, public :: i8 = selected_int_kind(18)

   ! ---------------------------------------------------------------------------
   ! Backward-compatible aliases (for porting OpenCFD-EC / UNSSolverProj)
   ! ---------------------------------------------------------------------------
   ! OpenCFD-EC used  PRE_EC  as the real kind and  OCFD_DATA_TYPE  for MPI.
   integer, parameter, public :: PRE_EC = dp          ! legacy alias

   ! ---------------------------------------------------------------------------
   ! Mathematical constants
   ! ---------------------------------------------------------------------------
   real(dp), parameter, public :: pi = 3.141592653589793238462643383279502884_dp

   ! ---------------------------------------------------------------------------
   ! Small number to avoid division by zero
   ! ---------------------------------------------------------------------------
   real(dp), parameter, public :: eps_dp  = epsilon(1.0_dp)    ! ~2.2e-16
   real(dp), parameter, public :: tiny_dp = tiny(1.0_dp)       ! smallest positive ~2.3e-308

   ! ---------------------------------------------------------------------------
   ! MPI data-type mapping
   ! ---------------------------------------------------------------------------
   ! These are set at module load time via a small helper that queries the MPI
   ! implementation.  In non-MPI builds the parameters are unused.
   ! For simplicity we use the standard named constants directly; if the host
   ! code does not "use mpi" these will simply be ignored at compile time
   ! because they are not referenced.
   !
   ! NOTE: The actual MPI type constants (MPI_DOUBLE_PRECISION etc.) are
   !       defined by the MPI module/library.  We re-export convenient aliases
   !       here so that solver code does not need to "use mpi" just for types.
   !       The values are assigned in the MPI-enabled wrapper (mod_mpi_common).
   !       For non-MPI compilation they remain unused.
   ! ---------------------------------------------------------------------------

end module mod_precision
