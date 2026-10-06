!===============================================================================
! mod_restart.f90 -- Field dump / restart for parallel UNSSolver (Step 7)
!
! Provides:
!   write_field_dump:
!       Rank 0 writes the GLOBAL field state (u / p / T + time-history
!       u_old / u_old_old / T_old / T_old_old + ts_order) to an ASCII dump
!       file.  Header is human-readable; field data uses es24.16 to
!       preserve all significant digits.  Non-root ranks are no-ops (the
!       global field only exists on rank 0).
!
!   read_field_dump:
!       Rank 0 reads the dump and MPI_Bcast-s the global field arrays to
!       all ranks.  On return fld_g is allocated on every rank with the
!       global mesh size and contains the dumped state.
!
!   scatter_global_to_local:
!       Copies fld_g into a local fld using lm%owned_global (maps local
!       owned 1..nowned -> global cell ID).  Halo cells are left at their
!       init zero -- the first halo exchange in the solver fills them.
!
!   load_or_partition:
!       Restart-aware partition orchestrator.  Tries to read the user-
!       supplied partition map (ctrl%partition_file); if it matches the
!       current nprocs and mesh size, reuses it (and broadcasts).  Other-
!       wise, falls back to partition_mesh_distributed (which re-runs
!       METIS on rank 0 and overwrites the default map file).
!
! File format (ASCII, line-oriented):
!   # UNSolver field dump v1
!   # source_mesh: <casfile>
!   # ncells: <N>
!   # ts_order: <0|1|2>
!   # generated: <YYYY-MM-DD HH:MM:SS>
!   # fields: u p T u_old u_old_old T_old T_old_old
!   <u(1,1)> <u(2,1)> <u(3,1)>          (one line per cell, 3 doubles)
!   ...
!   # end u
!   <p(1)>                                (one line per cell, 1 double)
!   ...
!   # end p
!   <T(1)>
!   ...
!   # end T
!   <u_old(1,1)> <u_old(2,1)> <u_old(3,1)>
!   ...
!   # end u_old
!   <u_old_old(1,1)> ...
!   ...
!   # end u_old_old
!   <T_old(1)>
!   ...
!   # end T_old
!   <T_old_old(1)>
!   ...
!   # end T_old_old
!===============================================================================
module mod_uns_restart
   use mod_precision, only: dp, ip, pi
   use mod_uns_mpi_core, only: myrank, nprocs, mpi_check, mpi_comm
   use mpi
   use mod_uns_mesh
   use mod_uns_connectivity, only: conn_t
   use mod_uns_fields
   use mod_uns_local_mesh, only: local_mesh_t
   use mod_uns_partition,   only: read_partition_map, partition_mesh_distributed
   implicit none
   private
   public :: write_field_dump, read_field_dump, scatter_global_to_local, &
             load_or_partition

   integer, parameter :: DUMP_VERSION = 1

contains

   !----------------------------------------------------------------------------
   ! Rank 0: write global fields to ASCII dump file.
   ! Non-root ranks are no-ops (caller may invoke collectively).
   !----------------------------------------------------------------------------
   subroutine write_field_dump( fname, m_g, fld_g, src_file, ier )
      character(len=*), intent(in)  :: fname
      type(mesh_t),     intent(in)  :: m_g
      type(fields_t),   intent(in)  :: fld_g
      character(len=*), intent(in)  :: src_file
      integer,          intent(out) :: ier
      integer :: u, ios, c
      character(len=8)  :: dstr
      character(len=10) :: tstr

      ier = 0
      if ( myrank /= 0 ) return   ! only rank 0 holds the global field

      open( newunit=u, file=trim(fname), status='replace', action='write', &
            form='formatted', iostat=ios )
      if ( ios /= 0 ) then
         write(*,'(a,a,a,i0)') 'ERROR: cannot open dump file ', trim(fname), &
                               ' for writing, iostat=', ios
         ier = 1; return
      end if

      call date_and_time( dstr, tstr )
      write(u,'(a,i0)')           '# UNSolver field dump v', DUMP_VERSION
      write(u,'(a,a)')           '# source_mesh: ', trim(src_file)
      write(u,'(a,i0)')          '# ncells: ', m_g%ncells
      write(u,'(a,i0)')          '# ts_order: ', fld_g%ts_order
      write(u,'(a,a,a,a,a,a,a,a)') '# generated: ', &
           dstr(1:4)//'-'//dstr(5:6)//'-'//dstr(7:8), ' ', &
           tstr(1:2)//':'//tstr(3:4)//':'//tstr(5:6)
      write(u,'(a)')             '# fields: u p T u_old u_old_old T_old T_old_old'

      ! ---- u (3 doubles per line per cell) -----------------------------------
      do c = 1, m_g%ncells
         write(u,'(3(es24.16,1x))') fld_g%u(:,c)
      end do
      write(u,'(a)') '# end u'

      ! ---- p -----------------------------------------------------------------
      do c = 1, m_g%ncells
         write(u,'(es24.16)') fld_g%p(c)
      end do
      write(u,'(a)') '# end p'

      ! ---- T -----------------------------------------------------------------
      do c = 1, m_g%ncells
         write(u,'(es24.16)') fld_g%T(c)
      end do
      write(u,'(a)') '# end T'

      ! ---- u_old -------------------------------------------------------------
      do c = 1, m_g%ncells
         write(u,'(3(es24.16,1x))') fld_g%u_old(:,c)
      end do
      write(u,'(a)') '# end u_old'

      ! ---- u_old_old --------------------------------------------------------
      do c = 1, m_g%ncells
         write(u,'(3(es24.16,1x))') fld_g%u_old_old(:,c)
      end do
      write(u,'(a)') '# end u_old_old'

      ! ---- T_old -------------------------------------------------------------
      do c = 1, m_g%ncells
         write(u,'(es24.16)') fld_g%T_old(c)
      end do
      write(u,'(a)') '# end T_old'

      ! ---- T_old_old --------------------------------------------------------
      do c = 1, m_g%ncells
         write(u,'(es24.16)') fld_g%T_old_old(c)
      end do
      write(u,'(a)') '# end T_old_old'

      close(u)
      write(*,'(a,a)') '  Field dump written to: ', trim(fname)
   end subroutine write_field_dump

   !----------------------------------------------------------------------------
   ! Rank 0 reads dump, then MPI_Bcast distributes the global field arrays to
   ! all ranks.  fld_g is allocated on every rank with global size ncells_g.
   ! src_file_out returns the source_mesh header field (for caller sanity
   ! check); it is broadcast to all ranks.
   !----------------------------------------------------------------------------
   subroutine read_field_dump( fname, ncells_g, fld_g, src_file_out, ier, serial )
      character(len=*),   intent(in)    :: fname
      integer,            intent(in)    :: ncells_g
      type(fields_t),     intent(inout) :: fld_g
      character(len=*),   intent(out)   :: src_file_out
      integer,            intent(out)   :: ier
      logical, optional,  intent(in)    :: serial
      integer :: u, ios, ver, n_in, ts_in, bcast_ierr
      character(len=512) :: line
      logical :: ser

      ser = .false.
      if ( present(serial) ) ser = serial

      ier = 0
      src_file_out = ''

      ! ---- allocate fld_g on all ranks (global size) -------------------------
      if ( .not. allocated(fld_g%u) ) then
         allocate( fld_g%u(3,ncells_g), fld_g%p(ncells_g) )
         allocate( fld_g%gp(3,ncells_g), fld_g%gu(3,3,ncells_g) )
         allocate( fld_g%flux(1), fld_g%lf(1), fld_g%apc(ncells_g) )
         allocate( fld_g%u_old(3,ncells_g), fld_g%u_old_old(3,ncells_g) )
         allocate( fld_g%T(ncells_g), fld_g%T_old(ncells_g), fld_g%T_old_old(ncells_g) )
         allocate( fld_g%gt(3,ncells_g) )
         fld_g%u = 0.0_dp;        fld_g%p = 0.0_dp
         fld_g%gp = 0.0_dp;       fld_g%gu = 0.0_dp
         fld_g%flux = 0.0_dp;     fld_g%apc = 1.0_dp
         fld_g%u_old = 0.0_dp;    fld_g%u_old_old = 0.0_dp
         fld_g%ts_order = 0
         fld_g%T = 0.0_dp;        fld_g%T_old = 0.0_dp
         fld_g%T_old_old = 0.0_dp
         fld_g%gt = 0.0_dp
      end if

      ! ---- rank 0 reads file (or every rank in serial mode) -------------------
      if ( myrank == 0 .or. ser ) then
         open( newunit=u, file=trim(fname), status='old', action='read', &
               form='formatted', iostat=ios )
         if ( ios /= 0 ) then
            write(*,'(a,a,a,i0)') 'ERROR: cannot open dump file ', trim(fname), &
                                  ' for reading, iostat=', ios
            ier = 1
         else
            ! ---- parse header lines starting with '#'; stop at '# fields:' -
            ver = 0; n_in = 0; ts_in = 0
            do
               read(u,'(a)',iostat=ios) line
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: dump file ended in header'
                  ier = 2; exit
               end if
               line = adjustl(line)
               if ( len_trim(line) == 0 ) cycle
               if ( line(1:1) /= '#' ) then
                  ! Not a header line -- file is malformed (we expect '# fields:'
                  ! as the last header line, then data follows).
                  write(*,'(a,a)') 'ERROR: dump file malformed header at: ', &
                                  trim(line)
                  ier = 2; exit
               end if
               if ( index(line,'# UNSolver field dump v') > 0 ) then
                  ! marker '# UNSolver field dump v' followed by version int;
                  ! do NOT use index(line,'v') -- 'v' also appears in 'UNSolver'.
                  read( line(len('# UNSolver field dump v')+1:), *, iostat=ios ) ver
               else if ( index(line,'# source_mesh:') > 0 ) then
                  read( line(index(line,':')+1:), '(a)') src_file_out
                  src_file_out = trim(adjustl(src_file_out))
               else if ( index(line,'# ncells:') > 0 ) then
                  read( line(index(line,':')+1:), *, iostat=ios ) n_in
               else if ( index(line,'# ts_order:') > 0 ) then
                  read( line(index(line,':')+1:), *, iostat=ios ) ts_in
               else if ( index(line,'# fields:') > 0 ) then
                  exit   ! header done, data follows
               end if
            end do

            if ( ier == 0 ) then
               if ( ver /= DUMP_VERSION ) then
                  write(*,'(a,i0,a,i0)') 'ERROR: dump version mismatch, file=', &
                        ver, ' expected=', DUMP_VERSION
                  ier = 3
               else if ( n_in /= ncells_g ) then
                  write(*,'(a,i0,a,i0)') 'ERROR: dump ncells=', n_in, &
                        ' but global mesh ncells=', ncells_g
                  ier = 4
               end if
            end if

            ! ---- read field data list-directed (crosses line boundaries) ----
            if ( ier == 0 ) then
               read(u,*,iostat=ios) fld_g%u
               if ( ios /= 0 ) then
                  write(*,'(a,i0)') 'ERROR: reading u failed, ios=', ios
                  ier = 5
               end if
            end if
            if ( ier == 0 ) then
               read(u,'(a)',iostat=ios) line   ! skip '# end u'
               read(u,*,iostat=ios) fld_g%p
               if ( ios /= 0 ) then; ier = 5; end if
            end if
            if ( ier == 0 ) then
               read(u,'(a)',iostat=ios) line   ! skip '# end p'
               read(u,*,iostat=ios) fld_g%T
               if ( ios /= 0 ) then; ier = 5; end if
            end if
            if ( ier == 0 ) then
               read(u,'(a)',iostat=ios) line   ! skip '# end T'
               read(u,*,iostat=ios) fld_g%u_old
               if ( ios /= 0 ) then; ier = 5; end if
            end if
            if ( ier == 0 ) then
               read(u,'(a)',iostat=ios) line   ! skip '# end u_old'
               read(u,*,iostat=ios) fld_g%u_old_old
               if ( ios /= 0 ) then; ier = 5; end if
            end if
            if ( ier == 0 ) then
               read(u,'(a)',iostat=ios) line   ! skip '# end u_old_old'
               read(u,*,iostat=ios) fld_g%T_old
               if ( ios /= 0 ) then; ier = 5; end if
            end if
            if ( ier == 0 ) then
               read(u,'(a)',iostat=ios) line   ! skip '# end T_old'
               read(u,*,iostat=ios) fld_g%T_old_old
               if ( ios /= 0 ) then; ier = 5; end if
            end if
            close(u)
            if ( ier == 0 ) then
               fld_g%ts_order = ts_in
               write(*,'(a,a)') '  Field dump read from: ', trim(fname)
               write(*,'(a,i0,a,i0,a,i0)') '    ncells=', n_in, &
                  '  ts_order=', ts_in, '  ver=', ver
            end if
         end if
      end if

      ! ---- broadcast error status to all ranks -------------------------------
      ! Skipped in serial mode: mpi_comm is MPI_COMM_WORLD in the coupled
      ! driver and the struct ranks never call this routine.
      if ( .not. ser ) then
         call MPI_Bcast( ier, 1, MPI_INTEGER, 0, mpi_comm, bcast_ierr )
      end if
      if ( ier /= 0 ) return

      ! ---- broadcast ts_order and field arrays to all ranks ------------------
      if ( .not. ser ) then
         call MPI_Bcast( fld_g%ts_order, 1, MPI_INTEGER, 0, mpi_comm, bcast_ierr )
         call MPI_Bcast( fld_g%u,         size(fld_g%u),         MPI_DOUBLE_PRECISION, 0, mpi_comm, bcast_ierr )
         call MPI_Bcast( fld_g%p,          size(fld_g%p),         MPI_DOUBLE_PRECISION, 0, mpi_comm, bcast_ierr )
         call MPI_Bcast( fld_g%T,          size(fld_g%T),         MPI_DOUBLE_PRECISION, 0, mpi_comm, bcast_ierr )
         call MPI_Bcast( fld_g%u_old,      size(fld_g%u_old),    MPI_DOUBLE_PRECISION, 0, mpi_comm, bcast_ierr )
         call MPI_Bcast( fld_g%u_old_old,  size(fld_g%u_old_old), MPI_DOUBLE_PRECISION, 0, mpi_comm, bcast_ierr )
         call MPI_Bcast( fld_g%T_old,      size(fld_g%T_old),    MPI_DOUBLE_PRECISION, 0, mpi_comm, bcast_ierr )
         call MPI_Bcast( fld_g%T_old_old,  size(fld_g%T_old_old), MPI_DOUBLE_PRECISION, 0, mpi_comm, bcast_ierr )

         ! ---- broadcast src_file_out (so caller can compare with mesh file) ----
         call MPI_Bcast( src_file_out, len(src_file_out), MPI_CHARACTER, 0, &
                         mpi_comm, bcast_ierr )
      end if

   end subroutine read_field_dump

   !----------------------------------------------------------------------------
   ! Scatter global fld_g into local fld using lm%owned_global mapping.
   ! Local owned cell k (1..nowned) corresponds to global cell owned_global(k).
   ! Halo cells (k > nowned) are left at their init zero -- the first halo
   ! exchange in the solver will fill them from neighbouring ranks.
   !----------------------------------------------------------------------------
   subroutine scatter_global_to_local( lm, fld_g, fld, ier )
      type(local_mesh_t), intent(in)    :: lm
      type(fields_t),     intent(in)    :: fld_g
      type(fields_t),     intent(inout) :: fld
      integer,            intent(out)   :: ier
      integer :: k, g
      ier = 0
      do k = 1, lm%nowned
         g = lm%owned_global(k)
         if ( g < 1 .or. g > size(fld_g%p) ) then
            write(*,'(a,i0,a,i0)') 'ERROR: scatter out-of-range, local=', k, &
                                  ' global=', g
            ier = 1; return
         end if
         fld%u(:,k)         = fld_g%u(:,g)
         fld%p(k)           = fld_g%p(g)
         fld%T(k)           = fld_g%T(g)
         fld%u_old(:,k)     = fld_g%u_old(:,g)
         fld%u_old_old(:,k) = fld_g%u_old_old(:,g)
         fld%T_old(k)       = fld_g%T_old(g)
         fld%T_old_old(k)   = fld_g%T_old_old(g)
      end do
      fld%ts_order = fld_g%ts_order
   end subroutine scatter_global_to_local

   !----------------------------------------------------------------------------
   ! Restart-aware partition orchestrator.
   !   1. Rank 0 inquires user_mapfile; if it exists, read it.
   !   2. If nprocs and ncells match, reuse the loaded part (broadcast).
   !   3. Otherwise, fall back to partition_mesh_distributed (which re-runs
   !      METIS on rank 0 and overwrites default_mapfile).
   !----------------------------------------------------------------------------
   subroutine load_or_partition( c_g, nprocs_in, part, edgecut, &
                                  user_mapfile, default_mapfile, &
                                  src_file, ier )
      type(conn_t), intent(in) :: c_g
      integer,      intent(in) :: nprocs_in
      integer, allocatable, intent(out) :: part(:)
      integer,      intent(out) :: edgecut, ier
      character(len=*), intent(in) :: user_mapfile, default_mapfile, src_file
      integer :: n_map, p_map, e_map, bcast_ierr, my_ierr
      logical :: exists
      integer, allocatable :: part_tmp(:)

      ier = 0; my_ierr = 0; edgecut = -1

      if ( myrank == 0 ) then
         inquire( file=trim(user_mapfile), exist=exists )
         if ( exists ) then
            call read_partition_map( trim(user_mapfile), part_tmp, &
                                      n_map, p_map, e_map, my_ierr )
            if ( my_ierr == 0 .and. p_map == nprocs_in .and. &
                 n_map == size(c_g%c2c_ptr)-1 ) then
               call move_alloc( part_tmp, part )
               edgecut = e_map
               write(*,'(a,a)') '  Reusing partition map: ', trim(user_mapfile)
               my_ierr = 0    ! success
            else
               if ( my_ierr /= 0 ) then
                  write(*,'(a)') '  Partition map unreadable, will re-partition'
               else
                  write(*,'(a,i0,a,i0,a,i0)') '  Map nprocs=', p_map, &
                       ' (current=', nprocs_in, ', ncells_map=', n_map, &
                       ') mismatch, re-partitioning'
               end if
               if ( allocated(part_tmp) ) deallocate(part_tmp)
               my_ierr = -1   ! signal re-partition
            end if
         else
            write(*,'(a,a)') '  Partition map not found, will partition: ', &
                             trim(user_mapfile)
            my_ierr = -1
         end if
      end if

      ! broadcast the decision so all ranks agree on the path
      call MPI_Bcast( my_ierr, 1, MPI_INTEGER, 0, mpi_comm, bcast_ierr )

      if ( my_ierr == 0 ) then
         ! ---- reuse path: broadcast loaded part + edgecut --------------------
         if ( .not. allocated(part) ) then
             ! non-root ranks need to allocate part to the global size
             allocate( part(size(c_g%c2c_ptr)-1) )
         end if
         call MPI_Bcast( part, size(part), MPI_INTEGER, 0, mpi_comm, bcast_ierr )
         call MPI_Bcast( edgecut, 1, MPI_INTEGER, 0, mpi_comm, bcast_ierr )
         ier = 0
      else
         ! ---- re-partition from scratch (overwrites default_mapfile) --------
         call partition_mesh_distributed( c_g, nprocs_in, part, edgecut, &
                                          trim(default_mapfile), &
                                          trim(src_file), ier )
      end if
   end subroutine load_or_partition

end module mod_uns_restart
