!===============================================================================
! mod_struct_flowvar.f90 -- per-block scratch flow-field arrays
! Extracted (phase 2a) from pristine OpenCFD-EC opencfd_ec3d_v1.16a.f90 so it
! can be compiled independently before program main (auto dependencies).
! Memory is allocated for the block currently being computed and freed after.
!===============================================================================
  module Flow_Var
   use precision_EC
   real(PRE_EC), save,pointer,dimension(:,:,:)::  d,uu,v,w,T,p,cc ! 密度、x-速度、y-速度、z-速度、压力、声速
   real(PRE_EC), save,pointer,dimension(:,:,:,:):: Flux                 ! i- ,j-及k-方向的通量
   real(PRE_EC), save,pointer,dimension(:,:,:):: Lvi,Lvj,Lvk,Lci,Lcj,Lck   ! 无粘项及粘性项Jocabian的谱半径 

  end module Flow_Var
