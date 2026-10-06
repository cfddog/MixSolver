!===============================================================================
! mod_partition.f90 -- Mesh partitioning for parallel UNSSolver (Step 2)
!
! Workflow:
!   1. build_dual_graph(conn) -- extract CSR dual graph from conn_t%c2c.
!      conn_t%c2c_ptr is already 0-based CSR row pointers; conn_t%c2c holds
!      1-based cell IDs. METIS v5 expects 0-based, so adjncy = c2c - 1.
!   2. partition_mesh_serial(conn, nparts) -- rank 0 only: build graph,
!      call METIS_PartGraphKway(nparts=nprocs). Returns part(1..ncells) with
!      0-based owner rank for each cell, plus edgecut.
!   3. write_partition_map -- ASCII map file (cell_id -> owner_rank) for
!      restart with different nproc. Format documented in file header.
!   4. read_partition_map -- read existing map for restart.
!   5. partition_mesh_distributed -- orchestrator: rank 0 does serial
!      partition + writes map, then MPI_Bcast distributes part to all ranks.
!
! Map file format (ASCII):
!   # UNSSolver partition map v1
!   # source_mesh: <casfile>
!   # ncells: <N>
!   # nprocs: <P>
!   # edgecut: <E>
!   # generated: <YYYY-MM-DD HH:MM:SS>
!   # columns: cell_id(1-based) owner_rank(0-based)
!   1 0
!   2 0
!   ...
!===============================================================================
module mod_uns_partition
   use mod_precision, only: dp, ip, pi
   use mod_uns_connectivity, only: conn_t
   use iso_c_binding,    only: c_int, c_int32_t, c_null_ptr, c_loc
   use mod_uns_metis_iface,  only: metis_partkway, metis_setdefaultoptions, metis_ok
   implicit none
   private
   public :: build_dual_graph, partition_mesh_serial, &
             write_partition_map, read_partition_map, &
             partition_mesh_distributed, partition_balance

   integer, parameter :: MAP_VERSION = 1

contains

   !----------------------------------------------------------------------------
   ! Build CSR dual graph from conn_t%c2c, METIS 0-based ready.
   !   xadj(1:ncells+1) : row pointers (copy of c2c_ptr, already 0-based)
   !   adjncy(1:nnz)    : column indices (c2c with 1-based -> 0-based conversion)
   !----------------------------------------------------------------------------
   subroutine build_dual_graph( conn, xadj, adjncy, ierr )
      type(conn_t),        intent(in)  :: conn
      integer(c_int32_t), allocatable, intent(out) :: xadj(:), adjncy(:)
      integer,             intent(out) :: ierr
      integer :: n, nnz, i
      ierr = 0
      if ( .not. allocated(conn%c2c_ptr) .or. .not. allocated(conn%c2c) ) then
         ierr = 1; return
      end if
      n    = size(conn%c2c_ptr) - 1
      nnz  = conn%c2c_ptr(n+1)
      allocate( xadj(n+1), adjncy(nnz) )
      do i = 1, n+1
         xadj(i) = int( conn%c2c_ptr(i), c_int32_t )
      end do
      do i = 1, nnz
         adjncy(i) = int( conn%c2c(i) - 1, c_int32_t )   ! 1-based -> 0-based
      end do
   end subroutine build_dual_graph

   !----------------------------------------------------------------------------
   ! Rank 0 only: build dual graph and call METIS_PartGraphKway(nparts).
   ! Allocates part(1:ncells) with 0-based owner rank per cell.
   !----------------------------------------------------------------------------
   subroutine partition_mesh_serial( conn, nparts, part, edgecut, ierr )
      type(conn_t),        intent(in)  :: conn
      integer,             intent(in)  :: nparts
      integer, allocatable, intent(out) :: part(:)
      integer,             intent(out) :: edgecut, ierr
      integer(c_int32_t), allocatable :: xadj(:), adjncy(:), part_c(:)
      integer(c_int32_t), target :: options(40)
      integer(c_int32_t) :: nvtxs, ncon, nparts_c, objval
      integer(c_int) :: ret
      integer :: n, i
      ierr = 0; edgecut = -1
      n = size(conn%c2c_ptr) - 1

      ! Degenerate single-partition case: METIS k-way divides by zero when
      ! nparts == 1, so assign all cells to partition 0 directly.
      if ( nparts <= 1 ) then
         allocate( part(n) )
         part = 0
         edgecut = 0
         return
      end if

      call build_dual_graph( conn, xadj, adjncy, ierr )
      if ( ierr /= 0 ) return

      allocate( part(n), part_c(n) )
      part   = -1
      part_c = -1
      nvtxs   = int( n,       c_int32_t )
      ncon    = 1_c_int32_t
      nparts_c = int( nparts, c_int32_t )
      objval  = -1_c_int32_t

      ret = metis_setdefaultoptions( options )
      if ( ret /= metis_ok ) then; ierr = 2; return; end if

      ret = metis_partkway( nvtxs, ncon, xadj, adjncy, &
           c_null_ptr, c_null_ptr, c_null_ptr, &   ! vwgt, vsize, adjwgt
           nparts_c, c_null_ptr, c_null_ptr, &     ! tpwgts, ubvec (NULL = equal/1.05)
           c_loc(options), objval, part_c )
      if ( ret /= metis_ok ) then; ierr = 3; return; end if

      do i = 1, n
         part(i) = int( part_c(i) )
      end do
      edgecut = int( objval )
      deallocate( xadj, adjncy, part_c )
   end subroutine partition_mesh_serial

   !----------------------------------------------------------------------------
   ! Write ASCII partition map file (for restart with different nproc).
   !----------------------------------------------------------------------------
   subroutine write_partition_map( filename, part, ncells, nparts, edgecut, &
                                    src_file, ierr )
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: part(:), ncells, nparts, edgecut
      character(len=*), intent(in)  :: src_file
      integer,          intent(out) :: ierr
      integer :: unit, i
      character(len=8)  :: date_str
      character(len=10) :: time_str
      ierr = 0
      open( newunit=unit, file=trim(filename), status='replace', &
            action='write', form='formatted', iostat=ierr )
      if ( ierr /= 0 ) return
      call date_and_time( date_str, time_str )
      write(unit,'(a,i0)')     '# UNSSolver partition map v', MAP_VERSION
      write(unit,'(a,a)')     '# source_mesh: ', trim(src_file)
      write(unit,'(a,i0)')    '# ncells: ', ncells
      write(unit,'(a,i0)')    '# nprocs: ', nparts
      write(unit,'(a,i0)')    '# edgecut: ', edgecut
      write(unit,'(a,a,a,a)') '# generated: ', &
           date_str(1:4)//'-'//date_str(5:6)//'-'//date_str(7:8), ' ', &
           time_str(1:2)//':'//time_str(3:4)//':'//time_str(5:6)
      write(unit,'(a)')       '# columns: cell_id(1-based) owner_rank(0-based)'
      do i = 1, ncells
         write(unit,'(i0,1x,i0)') i, part(i)
      end do
      close(unit)
   end subroutine write_partition_map

   !----------------------------------------------------------------------------
   ! Read ASCII partition map file (for restart).
   ! Returns: part(1:ncells), ncells, nparts, edgecut.
   !----------------------------------------------------------------------------
   subroutine read_partition_map( filename, part, ncells, nparts, edgecut, ierr )
      character(len=*),   intent(in)  :: filename
      integer, allocatable, intent(out) :: part(:)
      integer,            intent(out) :: ncells, nparts, edgecut, ierr
      integer :: unit, cell_id, owner, n_parsed, ios
      character(len=512) :: line
      ierr = 0; ncells = 0; nparts = 0; edgecut = -1
      open( newunit=unit, file=trim(filename), status='old', action='read', &
            form='formatted', iostat=ierr )
      if ( ierr /= 0 ) return
      n_parsed = 0
      do
         read( unit, '(a)', iostat=ios ) line
         if ( ios /= 0 ) exit
         line = adjustl( line )
         if ( len_trim(line) == 0 ) cycle
         if ( line(1:1) == '#' ) then
            if ( index(line, 'ncells:') > 0 ) then
               read( line(index(line,':')+1:), * ) ncells
            else if ( index(line, 'nprocs:') > 0 ) then
               read( line(index(line,':')+1:), * ) nparts
            else if ( index(line, 'edgecut:') > 0 ) then
               read( line(index(line,':')+1:), * ) edgecut
            end if
         else
            if ( ncells > 0 .and. .not. allocated(part) ) then
               allocate( part(ncells) ); part = -1
            end if
            if ( allocated(part) ) then
               read( line, * ) cell_id, owner
               if ( cell_id >= 1 .and. cell_id <= ncells ) then
                  part(cell_id) = owner
                  n_parsed = n_parsed + 1
               end if
            end if
         end if
      end do
      close( unit )
      if ( .not. allocated(part) ) then; ierr = 1; return; end if
      if ( n_parsed /= ncells ) ierr = 2
      if ( any(part < 0) ) ierr = 3
   end subroutine read_partition_map

   !----------------------------------------------------------------------------
   ! Compute balance statistics: count cells per partition.
   ! Returns counts(0:nparts-1) and imbalance = max(counts)/avg - 1.
   !----------------------------------------------------------------------------
   subroutine partition_balance( part, nparts, counts, imbalance )
      integer, intent(in)  :: part(:), nparts
      integer, allocatable, intent(out) :: counts(:)
      real(dp), intent(out) :: imbalance
      integer :: i, p
      real(dp) :: avg
      allocate( counts(0:nparts-1) ); counts = 0
      do i = 1, size(part)
         p = part(i)
         if ( p >= 0 .and. p < nparts ) counts(p) = counts(p) + 1
      end do
      avg = real( size(part), dp ) / real( nparts, dp )
      if ( avg > 0.0_dp ) then
         imbalance = real( maxval(counts), dp ) / avg - 1.0_dp
      else
         imbalance = 0.0_dp
      end if
   end subroutine partition_balance

   !----------------------------------------------------------------------------
   ! Orchestrator: rank 0 does serial partition + writes map, then broadcasts
   ! part[] to all ranks.  conn is only accessed on rank 0; non-root ranks may
   ! pass an unallocated conn_t.  All ranks receive part(1:ncells).
   ! Caller must have called mpi_bootstrap() first.
   !----------------------------------------------------------------------------
   subroutine partition_mesh_distributed( conn, nparts, part, edgecut, mapfile, &
                                          src_file, ierr )
#ifdef HAVE_MPI
      use mod_uns_mpi_core, only: myrank
   use mpi
      use mpi
#endif
      type(conn_t),     intent(in)  :: conn
      integer,           intent(in)  :: nparts
      integer, allocatable, intent(out) :: part(:)
      integer,           intent(out) :: edgecut, ierr
      character(len=*),  intent(in)  :: mapfile, src_file
      integer :: n, bcast_ierr
      ierr = 0; edgecut = -1
      ! Rank 0 extracts ncells from conn; broadcast to all ranks
      if ( myrank == 0 ) then
         if ( .not. allocated(conn%c2c_ptr) ) then
            write(*,'(a)') 'FATAL: rank 0 conn%c2c_ptr not allocated'
            ierr = 1
         else
            n = size(conn%c2c_ptr) - 1
         end if
      end if
      call MPI_Bcast( ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, bcast_ierr )
      if ( ierr /= 0 ) then
         call MPI_Abort( MPI_COMM_WORLD, ierr, bcast_ierr )
      end if
      call MPI_Bcast( n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, bcast_ierr )
      ! Allocate part on all ranks
      allocate( part(n) ); part = -1
      ! Rank 0: build graph, call METIS, write map
      if ( myrank == 0 ) then
         call partition_mesh_serial( conn, nparts, part, edgecut, ierr )
         if ( ierr /= 0 ) then
            write(*,'(a,i0)') 'FATAL: partition_mesh_serial failed, ierr=', ierr
         else if ( mapfile /= '' ) then
            call write_partition_map( mapfile, part, n, nparts, edgecut, &
                                      src_file, bcast_ierr )
            if ( bcast_ierr /= 0 ) then
               write(*,'(a,i0,a)') 'WARN: write map failed, ierr=', bcast_ierr, &
                    ' (continuing)'
               bcast_ierr = 0
            else
               write(*,'(a,a)') '  Partition map written to: ', trim(mapfile)
            end if
         end if
      end if
      ! Broadcast partition + edgecut + ierr to all ranks
      call MPI_Bcast( part,    n, MPI_INTEGER, 0, MPI_COMM_WORLD, bcast_ierr )
      call MPI_Bcast( edgecut, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, bcast_ierr )
      call MPI_Bcast( ierr,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, bcast_ierr )
      if ( ierr /= 0 ) then
         if ( myrank == 0 ) write(*,'(a)') 'FATAL: partitioning failed, aborting'
         call MPI_Abort( MPI_COMM_WORLD, ierr, bcast_ierr )
      end if
   end subroutine partition_mesh_distributed

end module mod_uns_partition
