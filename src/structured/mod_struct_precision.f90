!===============================================================================
! mod_struct_precision.f90 -- precision/MPI bridge for ported OpenCFD-EC code
!
! The pristine OpenCFD-EC declares its own module 'precision_EC' holding PRE_EC
! and (via the in-module 'include "mpif.h"') ALL MPI named constants; those are
! re-exported through the use-chain const_var -> Global_Var to every legacy
! file, which is how routines without their own 'use mpi' see MPI_COMM_WORLD
! etc.  During the phase-2a lift-and-shift port this behaviour is preserved:
!   * 'use mpi' (Fortran 2008) replaces 'include "mpif.h"'; MPI entities stay
!     public, so no legacy file needs an added 'use mpi';
!   * PRE_EC is defined with the same selected_real_kind() expression as the
!     project-wide 'dp' in common/mod_precision.f90 (identical kind value,
!     verified at link time), but 'dp' itself is NOT imported: several legacy
!     routines declare local variables named 'dp' (density increment) which
!     would clash with a use-associated 'dp'.
!===============================================================================
module precision_EC
   use mpi
   implicit none

   ! Double precision: identical kind parameter to common/mod_precision 'dp'
   integer, parameter :: PRE_EC = selected_real_kind(15, 307)

   ! MPI element data type matching PRE_EC
   integer, parameter :: OCFD_DATA_TYPE = MPI_DOUBLE_PRECISION
end module precision_EC
