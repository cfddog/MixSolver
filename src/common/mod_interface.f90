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

   !---------------------------------------------------------------------------
   ! interface-type classification (phase 11 dispatch table)
   !
   ! The type is a pure function of (own solver, own cell-zone class, peer
   ! solver, peer cell-zone class).  The structured side is always
   ! compressible + fluid; only the unstructured side distinguishes fluid vs
   ! porous cells, so a cross-solver interface type is decided by the
   ! unstructured cell-zone class.
   !---------------------------------------------------------------------------
   integer, parameter :: IFACE_UNKNOWN           = 0
   integer, parameter :: IFACE_COMP_FLUID_FLUID  = 1  ! struct(comp|fluid) <-> uns(low-speed|fluid)  : flow B
   integer, parameter :: IFACE_COMP_FLUID_POROUS = 2  ! struct(comp|fluid) <-> uns(low-speed|porous) : C2
   integer, parameter :: IFACE_UNS_FLUID_POROUS  = 3  ! uns(fluid)        <-> uns(porous)            : C1 (internal)

   ! own-side cell-zone class stored on each registered face.  Mirrors the
   ! unstructured CZ_FLUID / CZ_POROUS constants without creating a
   ! common -> unstructured module dependency.
   integer, parameter :: IFACE_CZ_FLUID  = 1
   integer, parameter :: IFACE_CZ_POROUS = 2

   !---- exchange-quantity list vocabulary (phase 11) -----------------------
   ! exchanged quantity id
   integer, parameter :: Q_RHO = 1, Q_U = 2, Q_T = 3, Q_P = 4
   ! role of an exchanged quantity at the interface (bit flags)
   integer, parameter :: XQ_NONE      = 0
   integer, parameter :: XQ_DIRICHLET = 1   ! state imposed as Dirichlet on the receiver
   integer, parameter :: XQ_FLUX      = 2   ! flux-type (Neumann / zero-gradient)
   integer, parameter :: XQ_JUMP      = 4   ! jump condition (slip / stress)
   integer, parameter :: XQ_CHARACTER = 8   ! characteristic / linearised-Riemann state

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
      ! ---- interface dispatch (phase 11) ----
      ! own-side cell-zone class of the cell touching this interface face:
      !   IFACE_CZ_FLUID / IFACE_CZ_POROUS  (structured side is always fluid)
      integer :: cz_type    = IFACE_CZ_FLUID
      ! own-side local face index (unstructured: m%f() index; structured: 0).
      ! Lets the dispatch recover the owning cell's cell-zone type.
      integer :: loc_face   = 0
      ! classified interface type; set by dispatch_interfaces /
      ! report_interface_dispatch after registration / pairing.
      integer :: iface_type = IFACE_UNKNOWN
   end type Interface_FACE_TYPE

   ! registered coupling faces (structured side now; unstructured side in
   ! phase 3).  Geometry is stored in each side's local units; the coupling
   ! layer converts to SI at exchange time.
   type(Interface_FACE_TYPE), allocatable :: Interface_List(:)
   integer :: Num_Interface = 0

contains

   !---------------------------------------------------------------------------
   ! Classify one coupling interface from the solver + cell-zone class of both
   ! sides.  Cross-solver (struct <-> uns): the structured side is always
   ! compressible + fluid, so the type is decided by the unstructured peer's
   ! cell-zone class.  uns <-> uns is a purely internal fluid/porous interface.
   !---------------------------------------------------------------------------
   pure integer function classify_interface( solver_own, cz_own, &
                                             solver_peer, cz_peer ) result( it )
      integer, intent(in) :: solver_own, cz_own, solver_peer, cz_peer
      integer :: uns_porous

      if ( (solver_own == PEER_STRUCT .and. solver_peer == PEER_UNS) .or. &
           (solver_own == PEER_UNS .and. solver_peer == PEER_STRUCT) ) then
         ! cross-solver: the unstructured cell-zone class picks fluid vs porous
         uns_porous = 0
         if ( solver_own  == PEER_UNS .and. cz_own  == IFACE_CZ_POROUS ) uns_porous = 1
         if ( solver_peer == PEER_UNS .and. cz_peer == IFACE_CZ_POROUS ) uns_porous = 1
         if ( uns_porous == 1 ) then
            it = IFACE_COMP_FLUID_POROUS
         else
            it = IFACE_COMP_FLUID_FLUID
         end if
      else if ( solver_own == PEER_UNS .and. solver_peer == PEER_UNS ) then
         ! internal uns fluid/porous interface (C1) -- not exchanged cross-solver
         it = IFACE_UNS_FLUID_POROUS
      else
         it = IFACE_UNKNOWN
      end if
   end function classify_interface

   pure function iface_type_name( it ) result( s )
      integer, intent(in) :: it
      character(len=32) :: s
      select case ( it )
      case ( IFACE_COMP_FLUID_FLUID );  s = 'comp-fluid<->lowspeed-fluid'
      case ( IFACE_COMP_FLUID_POROUS ); s = 'comp-fluid<->lowspeed-porous'
      case ( IFACE_UNS_FLUID_POROUS );  s = 'fluid<->porous (internal)'
      case default;                     s = 'unknown'
      end select
   end function iface_type_name

   pure function quantity_name( q ) result( s )
      integer, intent(in) :: q
      character(len=8) :: s
      select case ( q )
      case ( Q_RHO ); s = 'rho'
      case ( Q_U );   s = 'u'
      case ( Q_T );   s = 'T'
      case ( Q_P );   s = 'p'
      case default;   s = '?'
      end select
   end function quantity_name

   pure function role_name( r ) result( s )
      integer, intent(in) :: r
      character(len=16) :: s
      select case ( r )
      case ( XQ_DIRICHLET ); s = 'dirichlet'
      case ( XQ_FLUX );      s = 'flux'
      case ( XQ_JUMP );      s = 'jump'
      case ( XQ_CHARACTER ); s = 'characteristic'
      case default;          s = '-'
      end select
   end function role_name

   !---------------------------------------------------------------------------
   ! Exchange-quantity list for one interface type.  Fills qid(:) / role(:) with
   ! the exchanged quantities and their role, viewed from this side's
   ! perspective (a quantity flagged XQ_DIRICHLET is imposed by the peer on this
   ! side; XQ_CHARACTER is a characteristic / back-pressure state).  Returns the
   ! number of entries; qid/role must hold at least 4 elements.
   !---------------------------------------------------------------------------
   pure subroutine iface_exchange_recipe( it, qid, role, nq )
      integer, intent(in)  :: it
      integer, intent(out) :: qid(4), role(4)
      integer, intent(out) :: nq
      nq = 0
      select case ( it )
      case ( IFACE_COMP_FLUID_FLUID, IFACE_COMP_FLUID_POROUS )
         ! flow B / C2: velocity + T Dirichlet on the low-speed side, pressure and
         ! density passed by characteristic to the compressible ghost cell.
         nq = 4
         qid(1) = Q_U;   role(1) = XQ_DIRICHLET
         qid(2) = Q_T;   role(2) = XQ_DIRICHLET
         qid(3) = Q_P;   role(3) = XQ_CHARACTER
         qid(4) = Q_RHO; role(4) = XQ_CHARACTER
      case ( IFACE_UNS_FLUID_POROUS )
         ! internal C1: Beavers-Joseph slip + stress jump, no cross-solver exchange
         nq = 1
         qid(1) = Q_U;   role(1) = XQ_JUMP
      case default
         nq = 0
      end select
   end subroutine iface_exchange_recipe

   ! human-readable recipe string, e.g. "u:dirichlet T:dirichlet p:characteristic"
   pure function iface_recipe_string( it ) result( s )
      integer, intent(in) :: it
      character(len=96) :: s
      integer :: n, nq, qid(4), role(4)
      s = ''
      call iface_exchange_recipe( it, qid, role, nq )
      do n = 1, nq
         s = trim(s) // ' ' // trim(quantity_name(qid(n))) &
                    // ':' // trim(role_name(role(n)))
      end do
      s = adjustl( s )
   end function iface_recipe_string

   !---------------------------------------------------------------------------
   ! Per-side dispatch: classify every registered face of THIS side against an
   ! assumed peer solver / cell-zone class and print the dispatch table.  Used
   ! by the coupled driver where each rank group only holds its own entries.
   !---------------------------------------------------------------------------
   subroutine report_interface_dispatch( peer_solver, peer_cz )
      integer, intent(in) :: peer_solver   ! PEER_STRUCT / PEER_UNS
      integer, intent(in) :: peer_cz       ! assumed peer cell-zone class
      integer :: i, it, cnt(3)
      character(len=40) :: peer_lbl

      if ( Num_Interface <= 0 ) return
      cnt = 0
      do i = 1, Num_Interface
         it = classify_interface( Interface_List(i)%solver, &
                                  Interface_List(i)%cz_type, &
                                  peer_solver, peer_cz )
         Interface_List(i)%iface_type = it
         select case ( it )
         case ( IFACE_COMP_FLUID_FLUID );  cnt(1) = cnt(1) + 1
         case ( IFACE_COMP_FLUID_POROUS ); cnt(2) = cnt(2) + 1
         case ( IFACE_UNS_FLUID_POROUS );  cnt(3) = cnt(3) + 1
         end select
      end do

      if ( peer_solver == PEER_STRUCT ) then
         peer_lbl = 'struct (compressible|fluid)'
      else
         peer_lbl = 'uns (low-speed)'
      end if

      write(*,'(a)') '--- Interface dispatch table (phase 11) ---'
      write(*,'(a,a)')  '  assumed peer : ', trim(peer_lbl)
      write(*,'(a,i0)') '  faces        : ', Num_Interface
      if ( cnt(1) > 0 ) then
         write(*,'(a,i0,a,a)') '  faces ', cnt(1), '  ', &
            trim(iface_type_name(IFACE_COMP_FLUID_FLUID))
         write(*,'(a,a)')      '        recipe:', &
            trim(iface_recipe_string(IFACE_COMP_FLUID_FLUID))
      end if
      if ( cnt(2) > 0 ) then
         write(*,'(a,i0,a,a)') '  faces ', cnt(2), '  ', &
            trim(iface_type_name(IFACE_COMP_FLUID_POROUS))
         write(*,'(a,a)')      '        recipe:', &
            trim(iface_recipe_string(IFACE_COMP_FLUID_POROUS))
      end if
      if ( cnt(3) > 0 ) then
         write(*,'(a,i0,a,a)') '  faces ', cnt(3), '  ', &
            trim(iface_type_name(IFACE_UNS_FLUID_POROUS))
         write(*,'(a,a)')      '        recipe:', &
            trim(iface_recipe_string(IFACE_UNS_FLUID_POROUS))
      end if
      write(*,'(a)') '--- end dispatch ---'
   end subroutine report_interface_dispatch

end module mod_interface
