!===============================================================================
! mod_interface.f90 -- mixed structured/unstructured coupling interface type
!
! Defines BC_Interface, the registration type for gridgen "generic: 8"
! coupling faces.  This is a pure marker/registration type: it records which
! faces are coupling interfaces plus enough geometry (centroid, unit normal,
! area, bounding box) for the phase-3 geometric auto-matching against the
! unstructured side (Fluent CAS interface zones).
!
! Scope decisions (phase 2b s4b):
!   * pure marker / registration -- no coupling algorithm, no SI conversion
!     (SI conversion happens in the separate coupling layer at exchange time)
!   * coexists with the legacy BC_MSG_TYPE -- the bc=8 face is additionally
!     registered here, the original bc_msg array is untouched
!   * geometry is computed at registration (init time, coordinates ready)
!
! Cross-solver: the unstructured side (phase 3) will register its CAS
! interface zones into the same Interface_List so a geometric match can pair
! each structured face with its unstructured peer.
!===============================================================================
module mod_interface
   use mod_precision, only: dp
   implicit none

   ! gridgen "generic" coupling-face code (undefined in pristine OpenCFD-EC,
   ! introduced by MixNS for the structured/unstructured coupling surface).
   integer, parameter :: BC_INTERFACE = 8

   ! peer-solver identification
   integer, parameter :: PEER_NONE = 0, PEER_STRUCT = 1, PEER_UNS = 2

   ! matching state
   integer, parameter :: MATCH_UNMATCHED = 0, MATCH_MATCHED = 1

   type :: Interface_FACE_TYPE
      ! ---- registration identity (written by the owning side) ----
      integer :: solver                 ! PEER_STRUCT / PEER_UNS
      integer :: block_no               ! struct: global block no; uns: zone id
      integer :: face                   ! struct: face dir 1-6; uns: 0 (n/a)
      integer :: f_no                   ! struct: sub-face no; uns: 0
      integer :: ib, ie, jb, je, kb, ke ! struct: sub-face ijk range; uns: 0
      ! ---- geometry (computed at registration; basis for auto-matching) ----
      real(dp) :: centroid(3)           ! area-weighted face centre
      real(dp) :: normal(3)             ! area-weighted mean unit normal
      real(dp) :: area                  ! total face area
      real(dp) :: bbox_min(3)           ! bounding box lower corner
      real(dp) :: bbox_max(3)           ! bounding box upper corner
      ! ---- face vertices (phase 4: needed for projection + bilinear weights) ----
      ! For structured quads nv=4 ordered A,B,C,D CCW around the face normal.
      ! For unstructured polygonal faces nv = m%f(fi)%nn.
      integer  :: nv = 0
      real(dp), allocatable :: verts(:,:)  ! (3, nv) vertex coordinates
      ! ---- pairing state (filled by the phase-4 geometric matcher) ----
      integer :: match_state = MATCH_UNMATCHED
      integer :: peer_id = 0            ! index of the peer face in Interface_List
      ! ---- interpolation weights (phase 4: peer vertices contributing to this
      !      face centroid; size = peer%nv, sums to ~1 for a bilinear fit) ----
      real(dp), allocatable :: peer_w(:)
      ! ---- structured-side ghost-cell indices (phase 6) ----
      ! (ic,jc,kc) = inner cell next to the face, (ig,jg,kg) = ghost cell to write.
      integer :: ic, jc, kc
      integer :: ig, jg, kg
   end type Interface_FACE_TYPE

   ! registered coupling faces (structured side now; unstructured side in
   ! phase 3).  Geometry is stored in each side's local units; the coupling
   ! layer converts to SI at exchange time.
   type(Interface_FACE_TYPE), allocatable :: Interface_List(:)
   integer :: Num_Interface = 0

end module mod_interface
