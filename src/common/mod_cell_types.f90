!===============================================================================
! mod_cell_types.f90 -- Fluent cell element code -> VTK cell type (stage 17)
!
! The .cas reader stores the per-cell Fluent element code in mesh_t%ctype
! (section `(12 (zid first last act ctype)(...))`: the 5th header field of a
! uniform cell zone, or the per-cell token of a "mixed" zone, whose header
! carries ctype = 0).  This module is the single place that translates those
! codes into VTK cell types for the native VTU output path
! (`out_cells = native`, see docs/plan.md stage 17 and mod_uns_output).
!
! Scope decision (user, 2026-10-15; see docs/plan.md stage 17 and
! memory-bank/worklog.md):
!   * 4 = hexahedron and 6 = wedge are the only cell shapes present in this
!     repository's .cas files and the only ones mapped here:
!         4 -> VTK_HEXAHEDRON (12),  6 -> VTK_WEDGE (13)
!   * 1 (point), 2 (unspecified), 3 (line) and 5 (triangle) are 0-D/1-D/2-D
!     entities.  This program is three-dimensional and does not deal with
!     them; degenerate cells (zero volume) are not tolerated either.  Such
!     codes are therefore reported as NOT writable and the native writer stops
!     with an error (with the cell id and code in the message) instead of
!     silently emitting a broken or zero-volume cell.
!   * 7 = polyhedron has no usable VTK encoding (the legacy
!     VTK_POLYHEDRON face-stream is dropped by VTK 9 / ParaView) -> the caller
!     falls back to the tetrahedral fan for the whole file, with one warning.
!   * 0 is NOT a cell type: it only marks a "mixed" cell zone; the real
!     per-cell codes live in the record payload (the reader expands them).
!
! Evidence for the code values (see docs/plan.md stage 17):
!   * 4 = hex  : 55 zone records `(12 (zid 1 <n> 1 4))` in the tracked .cas
!                files, plus the in-repo generators cases/porous_plug/
!                gen_plug.py:74 and cases/heat_conduction/gen_annulus.py:99,
!                which write exactly that record with the comment "Hex cells";
!   * 6 = wedge: 33164 per-cell tokens in the payload of the mixed cell zones
!                (cases/cylinder_cas/cylinder.cas, cases/ltne/LNTE1D.cas,
!                grid_BC/unMesh.cas, cluster_cases/cylinder_Re{40,200});
!   * codes 1/2/3/5/7 never occur in this repository's inputs.
!===============================================================================
module mod_cell_types
   use mod_vtu_writer, only: VTK_HEXAHEDRON, VTK_WEDGE
   implicit none
   private
   public :: fluent_to_vtk, fluent_code_name, cell_type_supported, &
             cell_type_is_polyhedron, FLU_MIXED_ZONE, FLU_POLYHEDRON

   ! Fluent codes that the native writer treats specially
   integer, parameter :: FLU_MIXED_ZONE = 0   ! zone-level "mixed" marker
   integer, parameter :: FLU_POLYHEDRON = 7   ! no usable VTK encoding

contains

   !----------------------------------------------------------------------------
   ! VTK cell type for a Fluent element code; 0 = cannot be written natively.
   !----------------------------------------------------------------------------
   pure integer function fluent_to_vtk( code )

      integer, intent(in) :: code

      select case ( code )
      case ( 4 )
         fluent_to_vtk = VTK_HEXAHEDRON
      case ( 6 )
         fluent_to_vtk = VTK_WEDGE
      case default
         fluent_to_vtk = 0
      end select

   end function fluent_to_vtk

   !----------------------------------------------------------------------------
   ! True when the native writer can emit this Fluent code as one VTK cell.
   !----------------------------------------------------------------------------
   pure logical function cell_type_supported( code )

      integer, intent(in) :: code

      cell_type_supported = ( fluent_to_vtk(code) > 0 )

   end function cell_type_supported

   !----------------------------------------------------------------------------
   ! True for the polyhedron code (fan fallback, not an error).
   !----------------------------------------------------------------------------
   pure logical function cell_type_is_polyhedron( code )

      integer, intent(in) :: code

      cell_type_is_polyhedron = ( code == FLU_POLYHEDRON )

   end function cell_type_is_polyhedron

   !----------------------------------------------------------------------------
   ! Human-readable name of a Fluent element code, for error messages.
   !----------------------------------------------------------------------------
   pure function fluent_code_name( code ) result( name )

      integer, intent(in) :: code
      character(len=56)   :: name

      select case ( code )
      case ( 0 )
         name = 'mixed cell zone (no per-cell type)'
      case ( 1 )
         name = 'point/vertex (0-D: not used by this 3-D solver)'
      case ( 2 )
         name = 'unspecified 0-D/1-D/2-D entity (not supported)'
      case ( 3 )
         name = 'line/edge (1-D: not used by this 3-D solver)'
      case ( 4 )
         name = 'hexahedron'
      case ( 5 )
         name = 'triangle (2-D: not used by this 3-D solver)'
      case ( 6 )
         name = 'wedge (triangular prism)'
      case ( 7 )
         name = 'polyhedron (no usable VTK encoding)'
      case default
         name = 'unknown Fluent element code'
      end select

   end function fluent_code_name

end module mod_cell_types
