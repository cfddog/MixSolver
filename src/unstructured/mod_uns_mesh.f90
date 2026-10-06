!===============================================================================
! mod_mesh.f90 -- Mesh data types for unstructured Fluent-format meshes
!
! Conventions:
!   - Faces store owner cell (c0) and neighbor cell (c1); c1 = 0 means the
!     face is a boundary face.
!   - Face zone condition types follow the Fluent zone type codes, but zone
!     condition names read from the (45 ... ) records take precedence.
!===============================================================================
module mod_uns_mesh
   use mod_precision, only: dp, ip, pi
   implicit none
   private
   public :: face_t, zone_t, mesh_t

   ! Face connectivity record
   type :: face_t
      integer              :: nn = 0        ! number of nodes of the face
      integer, allocatable :: nodes(:)      ! node indices (1-based)
      integer              :: c0 = 0        ! owner cell index
      integer              :: c1 = 0        ! neighbor cell index (0 = boundary)
      integer              :: zone = 0      ! face zone id
   end type face_t

   ! Zone (condition) information
   type :: zone_t
      integer              :: id = 0            ! zone id
      integer              :: cond_code = 0     ! numeric condition type (13 rec)
      integer              :: nf = 0            ! number of faces in zone
      ! Names are kept long: VC tags (e.g. "Zone 2 ..., VC: porous Fluid = 1")
      ! may sit well beyond the first 32 characters of a zone name.
      character(len=128)   :: cond_name = ''    ! condition name from (45 rec)
      character(len=128)   :: user_name = ''    ! user zone name from (45 rec)
   end type zone_t

   ! Mesh container
   type :: mesh_t
      integer                 :: nnodes = 0
      integer                 :: nfaces = 0
      integer                 :: ncells = 0
      real(dp), allocatable   :: x(:,:)        ! node coordinates, x(1:3,1:nnodes)
      type(face_t), allocatable :: f(:)        ! faces, f(1:nfaces)
      integer, allocatable    :: ctype(:)      ! per-cell Fluent element type
      integer                 :: nzone = 0
      type(zone_t), allocatable :: zone(:)     ! face zone table
      ! ---- cell (volume) zones, phase 3 step C ----
      integer, allocatable    :: czone(:)      ! per-cell cell-zone id (1:ncells)
      integer                 :: nczone = 0    ! number of cell zones
      type(zone_t), allocatable :: czt(:)      ! cell zone table (nf = cell count)
      ! per-cell resolved block type from .control "cell_zone" lines:
      ! 0 = unspecified, 1 = fluid, 2 = porous (filled by resolve_cell_zones)
      integer, allocatable    :: cztype(:)
      integer                 :: ndim = 3      ! problem dimension
   end type mesh_t

end module mod_uns_mesh
