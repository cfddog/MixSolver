!===============================================================================
! mod_interface_units.f90 -- unit conversion at the coupling interface (phase 5).
!
! IRON RULE (constraints.md section 2): at the interface each solver converts
! its own variables TO SI, the SI data are exchanged, and each solver converts
! the received SI data BACK to its own internal form.  All conversion factors
! come from mod_reference_state (runtime configurable via mix.control).
!
! Structured solver (OpenCFD-EC) internal variables are non-dimensional by
! free-stream dynamic scales (verified against mod_struct_init/mod_struct_bc:
! u_inf* = 1, p_inf* = 1/(gamma*Ma^2)):
!   rho* = rho / rho_ref
!   u*   = u   / U_inf        (U_inf = Ma*a_ref, key u_ref in mix.control)
!   T*   = T   / T_ref
!   p*   = p   / (rho_ref * U_inf^2)
! With u_ref unset (0) U_inf = a_ref, recovering the Ma=1 special case.
!
! Unstructured solver works entirely in SI, so no conversion is needed on its
! side -- the SI values are used directly.
!
! This module provides elementwise and array conversion routines.  Velocity is
! handled as a 3-vector (u,v,w) scaled uniformly by u_ref.
!===============================================================================
module mod_interface_units
   use mod_precision, only: dp
   use mod_reference_state, only: reference_state_t, get_ref_state
   implicit none
   private

   public :: struct_to_SI, SI_to_struct
   public :: struct_vel_to_SI, SI_vel_to_struct

contains

   !---------------------------------------------------------------------------
   ! scalar conversions (elementwise)
   !---------------------------------------------------------------------------
   subroutine struct_to_SI( rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                            rho, u, v, w, T, p, st )
      real(dp),           intent(in)            :: rho_nd, u_nd, v_nd, w_nd
      real(dp),           intent(in)            :: T_nd, p_nd
      real(dp),           intent(out)           :: rho, u, v, w, T, p
      type(reference_state_t), intent(in), optional :: st
      type(reference_state_t) :: r
      r = get_ref_state(); if ( present(st) ) r = st
      rho = rho_nd * r%rho_ref
      u   = u_nd   * r%u_ref
      v   = v_nd   * r%u_ref
      w   = w_nd   * r%u_ref
      T   = T_nd   * r%T_ref
      p   = p_nd   * r%p_scale
   end subroutine struct_to_SI

   subroutine SI_to_struct( rho, u, v, w, T, p, &
                            rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, st )
      real(dp),           intent(in)            :: rho, u, v, w, T, p
      real(dp),           intent(out)           :: rho_nd, u_nd, v_nd, w_nd
      real(dp),           intent(out)           :: T_nd, p_nd
      type(reference_state_t), intent(in), optional :: st
      type(reference_state_t) :: r
      r = get_ref_state(); if ( present(st) ) r = st
      rho_nd = rho / r%rho_ref
      u_nd   = u   / r%u_ref
      v_nd   = v   / r%u_ref
      w_nd   = w   / r%u_ref
      T_nd   = T   / r%T_ref
      p_nd   = p   / r%p_scale
   end subroutine SI_to_struct

   !---------------------------------------------------------------------------
   ! velocity-vector only (for convenience when assembling flux BCs)
   !---------------------------------------------------------------------------
   subroutine struct_vel_to_SI( u_nd, v_nd, w_nd, u, v, w, st )
      real(dp),           intent(in)            :: u_nd, v_nd, w_nd
      real(dp),           intent(out)           :: u, v, w
      type(reference_state_t), intent(in), optional :: st
      type(reference_state_t) :: r
      r = get_ref_state(); if ( present(st) ) r = st
      u = u_nd * r%u_ref
      v = v_nd * r%u_ref
      w = w_nd * r%u_ref
   end subroutine struct_vel_to_SI

   subroutine SI_vel_to_struct( u, v, w, u_nd, v_nd, w_nd, st )
      real(dp),           intent(in)            :: u, v, w
      real(dp),           intent(out)           :: u_nd, v_nd, w_nd
      type(reference_state_t), intent(in), optional :: st
      type(reference_state_t) :: r
      r = get_ref_state(); if ( present(st) ) r = st
      u_nd = u / r%u_ref
      v_nd = v / r%u_ref
      w_nd = w / r%u_ref
   end subroutine SI_vel_to_struct

end module mod_interface_units
