!===============================================================================
! mod_local_mesh.f90 -- Local mesh extraction for parallel UNSSolver (Step 3)
!
! Given the global mesh/conn/geom and a partition vector part(1:ncells_global)
! (part(c) = owner rank), extract_local_mesh builds a self-contained local
! domain consisting of:
!   - owned cells : global cells with part(c) == myrank
!   - halo cells  : distinct 1-layer neighbours of owned cells (via c2c) owned
!                   by other ranks
! and the local faces/nodes touching owned cells, all renumbered locally.
!
! Layout convention of the local domain:
!   local cell IDs 1..nowned           -> owned cells
!                    nowned+1..nowned+nghost -> halo cells
!
! Geometry strategy: remap from the GLOBAL geom_t, do not recompute locally.
! Rationale: halo cells are missing some of their faces in the extracted mesh
! (halo-halo and halo-foreign faces are dropped), so a local compute_geometry
! would give wrong halo volumes/centroids and, critically, wrong rc1 vectors
! on owned-halo faces.  Face node ordering and the c0/c1 roles are preserved
! during extraction, hence the already-oriented global Sf/Xf and global cell
! V/Xc can be copied verbatim and remain mutually consistent.
!
! Local connectivity (cell->face, cell->cell CSR) is rebuilt from the local
! mesh with build_connectivity, so it follows the local numbering.
!===============================================================================
module mod_uns_local_mesh
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   use mod_uns_connectivity
   use mod_uns_geometry
   use mod_uns_mpi_core, only: myrank
   use mpi
   implicit none
   private
   public :: local_mesh_t, extract_local_mesh

   type :: local_mesh_t
      type(mesh_t) :: m          ! local mesh in local numbering
      type(conn_t) :: c          ! local connectivity (local numbering)
      type(geom_t) :: g          ! local geometry, remapped from global geom

      integer :: nowned = 0      ! number of owned cells
      integer :: nghost = 0      ! number of halo (ghost) cells

      ! Global <-> local maps. g2l arrays are sized by GLOBAL entity counts;
      ! entries are 0 for entities absent from the local domain.
      integer, allocatable :: owned_global(:)  ! global IDs of owned cells (nowned)
      integer, allocatable :: halo_global(:)   ! global IDs of halo cells  (nghost)
      integer, allocatable :: cell_l2g(:)      ! local -> global cells   (nowned+nghost)
      integer, allocatable :: cell_g2l(:)      ! global -> local cells  (ncells_global)
      integer, allocatable :: face_l2g(:)      ! local -> global faces  (nfaces_local)
      integer, allocatable :: face_g2l(:)      ! global -> local faces  (nfaces_global)
      integer, allocatable :: node_l2g(:)      ! local -> global nodes  (nnodes_local)
      integer, allocatable :: node_g2l(:)      ! global -> local nodes  (nnodes_global)
   end type local_mesh_t

contains

   !----------------------------------------------------------------------------
   ! Extract the local domain for rank myrank.
   !----------------------------------------------------------------------------
   subroutine extract_local_mesh( m_g, c_g, g_g, part, lm, ierr )
      type(mesh_t),        intent(in)  :: m_g
      type(conn_t),        intent(in)  :: c_g
      type(geom_t),        intent(in)  :: g_g
      integer,             intent(in)  :: part(:)             ! size ncells_global
      type(local_mesh_t),  intent(out) :: lm
      integer,             intent(out) :: ierr

      integer :: nc_g, nf_g, nn_g, i, k, lid, nf_local, nn_local
      integer, allocatable :: halo_mark(:)   ! marker: global cell -> halo flag
      integer, allocatable :: face_keep(:)   ! global face -> keep flag
      integer, allocatable :: node_keep(:)   ! global node -> keep flag
      ierr = 0
      nc_g = m_g%ncells
      nf_g = m_g%nfaces
      nn_g = m_g%nnodes

      ! ---- 1. owned cells ---------------------------------------------------
      lm%nowned = count( part(1:nc_g) == myrank )
      allocate( lm%owned_global(lm%nowned), lm%cell_g2l(nc_g) )
      lm%cell_g2l = 0
      k = 0
      do i = 1, nc_g
         if ( part(i) == myrank ) then
            k = k + 1
            lm%owned_global(k) = i
            lm%cell_g2l(i) = k               ! owned local IDs 1..nowned
         end if
      end do

      ! ---- 2. halo cells: distinct foreign c2c neighbours of owned cells ---
      allocate( halo_mark(nc_g) ); halo_mark = 0
      do k = 1, lm%nowned
         i = lm%owned_global(k)
         call add_halo_neighbors( c_g, i, part, halo_mark )
      end do
      ! count distinct halo candidates, then allocate to the exact size
      lm%nghost = count( halo_mark(1:nc_g) /= 0 )
      allocate( lm%halo_global(lm%nghost), lm%cell_l2g(lm%nowned+lm%nghost) )
      ! assign halo local IDs deterministically (scan halo_mark in global ID
      ! order so results are reproducible independent of c2c traversal order)
      lm%nghost = 0
      do i = 1, nc_g
         if ( halo_mark(i) /= 0 .and. part(i) /= myrank ) then
            lm%nghost = lm%nghost + 1
            lm%halo_global(lm%nghost) = i
            lm%cell_g2l(i) = lm%nowned + lm%nghost
         end if
      end do
      lm%cell_l2g(1:lm%nowned) = lm%owned_global
      lm%cell_l2g(lm%nowned+1:lm%nowned+lm%nghost) = lm%halo_global
      deallocate( halo_mark )

      ! ---- 3. local faces: any face touching an OWNED cell -----------------
      allocate( face_keep(nf_g), lm%face_g2l(nf_g) )
      face_keep = 0; lm%face_g2l = 0
      nf_local = 0
      do i = 1, nf_g
         if ( face_touches_owned( m_g, i, lm%cell_g2l, lm%nowned ) ) then
            nf_local = nf_local + 1
            face_keep(i) = 1
            lm%face_g2l(i) = nf_local
         end if
      end do
      allocate( lm%face_l2g(nf_local) )
      k = 0
      do i = 1, nf_g
         if ( face_keep(i) == 1 ) then
            k = k + 1
            lm%face_l2g(k) = i
         end if
      end do

      ! ---- 4. local nodes: union of nodes of kept faces --------------------
      allocate( node_keep(nn_g), lm%node_g2l(nn_g) )
      node_keep = 0; lm%node_g2l = 0
      nn_local = 0
      do k = 1, nf_local
         i = lm%face_l2g(k)
         do lid = 1, m_g%f(i)%nn
            if ( node_keep(m_g%f(i)%nodes(lid)) == 0 ) then
               node_keep(m_g%f(i)%nodes(lid)) = 1
               nn_local = nn_local + 1
               lm%node_g2l(m_g%f(i)%nodes(lid)) = nn_local
            end if
         end do
      end do
      allocate( lm%node_l2g(nn_local) )
      do i = 1, nn_g
         if ( lm%node_g2l(i) > 0 ) lm%node_l2g(lm%node_g2l(i)) = i
      end do

      ! ---- 5. build local mesh_t -------------------------------------------
      call build_local_mesh( m_g, lm )

      deallocate( face_keep, node_keep )

      ! ---- 6. local connectivity (CSR in local numbering) ------------------
      call build_connectivity( lm%m, lm%c, ierr )
      if ( ierr /= 0 ) return

      ! ---- 7. remap geometry from global geom_t ----------------------------
      call remap_geometry( g_g, lm )

   end subroutine extract_local_mesh

   !----------------------------------------------------------------------------
   ! Mark foreign c2c neighbours of owned cell 'cell' as halo candidates.
   !----------------------------------------------------------------------------
   subroutine add_halo_neighbors( c_g, cell, part, halo_mark )
      type(conn_t), intent(in)    :: c_g
      integer,      intent(in)    :: cell
      integer,      intent(in)    :: part(:)
      integer,      intent(inout) :: halo_mark(:)
      integer :: p, neigh
      do p = c_g%c2c_ptr(cell)+1, c_g%c2c_ptr(cell+1)
         neigh = c_g%c2c(p)
         if ( part(neigh) /= myrank .and. halo_mark(neigh) == 0 ) then
            halo_mark(neigh) = 1
         end if
      end do
   end subroutine add_halo_neighbors

   !----------------------------------------------------------------------------
   ! True if global face 'i' has an OWNED cell as c0 or c1.
   ! Boundary faces (c1 = 0) are kept iff c0 is owned, avoiding duplication
   ! across ranks.
   !----------------------------------------------------------------------------
   logical function face_touches_owned( m_g, i, cell_g2l, nowned )
      type(mesh_t), intent(in) :: m_g
      integer,      intent(in) :: i, cell_g2l(:), nowned
      integer :: l0, l1
      face_touches_owned = .false.
      l0 = cell_g2l(m_g%f(i)%c0)
      if ( l0 >= 1 .and. l0 <= nowned ) then
         face_touches_owned = .true.; return
      end if
      if ( m_g%f(i)%c1 > 0 ) then
         l1 = cell_g2l(m_g%f(i)%c1)
         if ( l1 >= 1 .and. l1 <= nowned ) face_touches_owned = .true.
      end if
   end function face_touches_owned

   !----------------------------------------------------------------------------
   ! Populate lm%m (mesh_t) with local numbering. c0/c1 roles and face node
   ! order are preserved; only indices are remapped. Zone table copied as
   ! global metadata (faces reference zone ids).
   !----------------------------------------------------------------------------
   subroutine build_local_mesh( m_g, lm )
      type(mesh_t),       intent(in)    :: m_g
      type(local_mesh_t), intent(inout) :: lm
      integer :: nloc, nf_local, nn_local, k, i, ln, z
      nf_local = size(lm%face_l2g)
      nn_local = size(lm%node_l2g)
      nloc     = lm%nowned + lm%nghost

      lm%m%nnodes = nn_local
      lm%m%nfaces = nf_local
      lm%m%ncells = nloc
      lm%m%nzone  = m_g%nzone
      lm%m%ndim   = m_g%ndim

      allocate( lm%m%x(3,nn_local), lm%m%f(nf_local), lm%m%ctype(nloc), &
                lm%m%zone(m_g%nzone) )

      ! node coordinates
      do i = 1, nn_local
         lm%m%x(:,i) = m_g%x(:,lm%node_l2g(i))
      end do

      ! faces: remap node list + cell ids
      do k = 1, nf_local
         i = lm%face_l2g(k)
         ln = m_g%f(i)%nn
         allocate( lm%m%f(k)%nodes(ln) )
         lm%m%f(k)%nn   = ln
         lm%m%f(k)%zone = m_g%f(i)%zone
         do ln = 1, size(m_g%f(i)%nodes)
            lm%m%f(k)%nodes(ln) = lm%node_g2l(m_g%f(i)%nodes(ln))
         end do
         lm%m%f(k)%c0 = lm%cell_g2l(m_g%f(i)%c0)
         if ( m_g%f(i)%c1 > 0 ) then
            lm%m%f(k)%c1 = lm%cell_g2l(m_g%f(i)%c1)
         else
            lm%m%f(k)%c1 = 0
         end if
      end do

      ! cell types
      do i = 1, nloc
         lm%m%ctype(i) = m_g%ctype(lm%cell_l2g(i))
      end do

      ! cell (volume) zones: per-cell ids must follow the cells into the
      ! local mesh, otherwise setup_porous_fields finds no porous zone on
      ! any rank and MPI runs silently lose the Darcy/LTNE physics
      if ( allocated(m_g%czone) ) then
         allocate( lm%m%czone(nloc) )
         do i = 1, nloc
            lm%m%czone(i) = m_g%czone(lm%cell_l2g(i))
         end do
      end if
      if ( allocated(m_g%cztype) ) then
         allocate( lm%m%cztype(nloc) )
         do i = 1, nloc
            lm%m%cztype(i) = m_g%cztype(lm%cell_l2g(i))
         end do
      end if
      lm%m%nczone = m_g%nczone
      if ( allocated(m_g%czt) ) then
         allocate( lm%m%czt(m_g%nczone) )
         do z = 1, m_g%nczone
            lm%m%czt(z) = m_g%czt(z)
         end do
      end if

      ! zone table copied verbatim (global nf counts kept as metadata)
      do z = 1, m_g%nzone
         lm%m%zone(z) = m_g%zone(z)
      end do
   end subroutine build_local_mesh

   !----------------------------------------------------------------------------
   ! Remap all geometric quantities from the global geom_t, preserving the
   ! face orientations (c0 -> c1 roles preserved in build_local_mesh).
   !----------------------------------------------------------------------------
   subroutine remap_geometry( g_g, lm )
      type(geom_t),       intent(in)    :: g_g
      type(local_mesh_t), intent(inout) :: lm
      integer :: nf_local, nloc, k, i
      nf_local = size(lm%face_l2g)
      nloc     = lm%nowned + lm%nghost

      allocate( lm%g%xf(3,nf_local), lm%g%sf(3,nf_local), lm%g%area(nf_local) )
      allocate( lm%g%vol(nloc), lm%g%xc(3,nloc) )
      allocate( lm%g%rc0(3,nf_local), lm%g%rc1(3,nf_local) )

      do k = 1, nf_local
         i = lm%face_l2g(k)
         lm%g%xf(:,k)  = g_g%xf(:,i)
         lm%g%sf(:,k)  = g_g%sf(:,i)
         lm%g%area(k)  = g_g%area(i)
         lm%g%rc0(:,k) = g_g%rc0(:,i)
         lm%g%rc1(:,k) = g_g%rc1(:,i)
      end do
      do i = 1, nloc
         k = lm%cell_l2g(i)
         lm%g%vol(i) = g_g%vol(k)
         lm%g%xc(:,i) = g_g%xc(:,k)
      end do
   end subroutine remap_geometry

end module mod_uns_local_mesh
