!===============================================================================
! main_dispatch.f90 -- single self-dispatching executable (bin/mixnsolver).
!
! Takes NO command-line arguments.  The coupling / solver mode is declared in a
! REQUIRED co-located control file `mix.control` via the key `mode`:
!
!     mode = coupled | struct | uns
!
!   * COUPLED : weak-coupling; reads Mesh3d.x + control.ec (struct) and
!              unMesh.cas + unMesh.control (uns), all relative to cwd, and runs
!              the same split + group drivers as bin/mixsolver.
!   * STRUCT  : structured (compressible) solver; reads control.ec + Mesh3d.x.
!   * UNS     : unstructured (incompressible) solver; reads unMesh.cas +
!              unMesh.control.
!
! The mode is authoritative -- the dispatcher does NOT guess from which files
! happen to be present.  `mix.control` must exist and contain a valid `mode`
! key, otherwise the run aborts with a clear error.
!
! Every rank shares the cwd and parses the same mix.control, so the decision is
! deterministic on all ranks with NO MPI communication (and it completes before
! any MPI_Comm_split, because struct vs uns membership depends on the mode).
!
! Run modes:
!   - COUPLED : same split + group drivers as bin/mixsolver (mod_mix_driver).
!   - STRUCT  : mod_struct_driver::struct_solver_run on the full MPI_COMM_WORLD.
!   - UNS     : unstructured (incompressible) solver on the full communicator.
!               With '-np N > 1' in the MPI build the cells are METIS-
!               partitioned and each rank solves its local subdomain
!               (mod_uns_driver_mpi::uns_run_mpi, partition map <casstem>.part.map);
!               with '-np 1' (or in the serial build) the legacy full-mesh path
!               (mod_uns_driver::uns_solver_init/step) is kept, which holds the
!               FULL mesh on every rank.
!
! Optional key in mix.control selecting the unstructured path:
!
!     uns_parallel = auto | always | never        (default: auto)
!
!   auto   : partitioned when nproc > 1, full mesh otherwise (default)
!   always : always partitioned (MPI build only; also valid with -np 1, which
!             is useful to compare the partitioned solution against the legacy
!             one at the same rank count)
!   never  : never partitioned (legacy full-mesh replication on every rank)
!
! Optional rank-count declaration in mix.control (CHECKED against the launch):
!
!     n_struct_ranks = <int>      ranks for the structured solver
!     n_uns_ranks    = <int>      ranks for the unstructured solver
!     (sum = total ranks, which MUST equal MPI_Comm_size(MPI_COMM_WORLD))
!
!   Declaring one key without the other, or declaring a sum that differs from
!   the launched rank count, prints a diagnosis on rank 0 and ABORTS the job
!   (MPI_Abort) -- the guard against submitting a job whose -np / job-script
!   task count does not match the split the case was prepared for.  When
!   neither key is present the legacy implicit split is kept (COUPLED: 1
!   struct rank + the remaining ranks uns; STRUCT/UNS: every rank that side)
!   and only an informational line is printed, so old case directories keep
!   working unchanged.
!
!   Mode-specific rules:
!     COUPLED : n_struct_ranks = 1, n_uns_ranks = 1 (sum 2).  The cross-group
!               exchange connects the two group roots directly
!               (mod_coupling_exchange) and the uns group has no intra-group
!               collectives, so extra ranks block in the exchange (verified:
!               -np 3 hangs); struct-side interface extraction is single-rank
!               too (global Interface_List vs per-rank local blocks).
!     STRUCT  : n_struct_ranks = total ranks (n_uns_ranks must be 0).
!     UNS     : n_uns_ranks    = total ranks (n_struct_ranks must be 0).
!===============================================================================
program mixnsolver
   use mpi
   use mod_precision, only: dp
   use mod_mix_driver, only: struct_group_driver, uns_group_driver
   use mod_struct_driver, only: struct_solver_run
   use mod_uns_driver, only: uns_solver_init, uns_solver_step
#ifdef HAVE_MPI
   ! Partitioned parallel unstructured driver.  The module lives in
   ! mod_uns_driver_mpi.f90, which the serial build excludes by file name, so the
   ! USE is compiled for the MPI tree only (see the MODE_UNS branch below).
   use mod_uns_driver_mpi, only: uns_run_mpi
#endif
   use mod_reference_state, only: read_mix_control
   implicit none

   integer, parameter :: MODE_COUPLED = 1, MODE_UNS = 2, MODE_STRUCT = 3
   ! unstructured path selection (optional `uns_parallel` key in mix.control)
   integer, parameter :: UNS_AUTO = 0, UNS_ALWAYS = 1, UNS_NEVER = 2
   ! sentinel of the optional rank-count keys: value = key absent
   integer, parameter :: NS_UNSET = -1

   integer :: ierr, rank, nproc, colour, key, mode
   integer :: comm_sub, n_struct_ranks, uns_parallel
   integer :: n_struct_decl, n_uns_decl, n_struct_left

   ! ---- MPI init -------------------------------------------------------------
   call MPI_Init( ierr )
   call MPI_Comm_rank( MPI_COMM_WORLD, rank, ierr )
   call MPI_Comm_size( MPI_COMM_WORLD, nproc, ierr )

   ! ---- mode: REQUIRED from mix.control (no file-presence guessing) ----------
   call parse_mix_mode( mode, rank )

   ! ---- rank split declared in mix.control (n_struct_ranks / n_uns_ranks) ----
   ! Every rank shares the cwd and parses the same file, so the decision needs
   ! no MPI communication and completes before any MPI_Comm_split.  A
   ! declaration that cannot be honoured aborts the whole job (rank 0 prints
   ! the diagnosis).  `n_struct_left` = struct ranks resolved for the split.
   call parse_mix_int_key( 'n_struct_ranks', n_struct_decl, rank )
   call parse_mix_int_key( 'n_uns_ranks',    n_uns_decl,    rank )
   call resolve_rank_split( mode, nproc, rank, n_struct_decl, n_uns_decl, &
                            n_struct_left )

   ! ---- dispatch -------------------------------------------------------------
   select case ( mode )
   case ( MODE_COUPLED )
      ! reference state + coupling parameters from the same mix.control
      call read_mix_control( 'mix.control', ierr )

      ! struct ranks come from the mix.control split (resolve_rank_split):
      ! legacy 1 struct rank (or none at -np 1), or the declared value.
      n_struct_ranks = n_struct_left
      if ( rank < n_struct_ranks ) then
         colour = 0
      else
         colour = 1
      end if
      key = rank
      call MPI_Comm_split( MPI_COMM_WORLD, colour, key, comm_sub, ierr )

      if ( rank == 0 ) then
         write(*,'(a)') '=========================================================='
         write(*,'(a)') '  MixNSSolver [mixnsolver] -- mode: COUPLED'
         write(*,'(a,i0,a,i0,a,i0)') '  Ranks        : ', nproc, &
            ' total = struct ', n_struct_ranks, ' + uns ', nproc - n_struct_ranks
         write(*,'(a)') '=========================================================='
      end if
      if ( rank < n_struct_ranks ) then
         call struct_group_driver( comm_sub, 'Mesh3d.x', 'control.ec', rank, nproc )
      else
         call uns_group_driver( comm_sub, 'unMesh.cas', 'unMesh.control', rank, nproc )
      end if

   case ( MODE_STRUCT )
      if ( rank == 0 ) then
         write(*,'(a)') '=========================================================='
         write(*,'(a)') '  MixNSSolver [mixnsolver] -- mode: STRUCT (compressible)'
         write(*,'(a)') '  (reads control.ec + Mesh3d.x)'
         write(*,'(a,i0)') '  MPI ranks    : ', nproc
         write(*,'(a)') '=========================================================='
      end if
      call struct_solver_run( MPI_COMM_WORLD )

   case ( MODE_UNS )
      ! optional mix.control key selecting the unstructured run path
      call parse_mix_uns_parallel( uns_parallel, rank )
      if ( rank == 0 ) then
         write(*,'(a)') '=========================================================='
         write(*,'(a)') '  MixNSSolver [mixnsolver] -- mode: UNS (incompressible)'
         write(*,'(a)') '  Mesh file    : unMesh.cas'
         write(*,'(a)') '  Control file : unMesh.control'
         write(*,'(a,i0)') '  MPI ranks    : ', nproc
      end if
#ifdef HAVE_MPI
      if ( uns_parallel == UNS_ALWAYS .or. &
           ( uns_parallel == UNS_AUTO .and. nproc > 1 ) ) then
         if ( rank == 0 ) write(*,'(a)') &
            '  Path         : METIS-partitioned parallel unstructured solve'
         call uns_run_mpi( MPI_COMM_WORLD, 'unMesh.cas', 'unMesh.control' )
      else
         if ( rank == 0 ) write(*,'(a)') &
            '  Path         : full-mesh unstructured solve (no partition)'
         call run_uns( MPI_COMM_WORLD, 'unMesh.cas', 'unMesh.control' )
      end if
#else
      ! Serial build: mod_uns_driver_mpi / METIS are not linked here, so the
      ! full-mesh path is the only option (uns_parallel=always is ignored).
      if ( uns_parallel == UNS_ALWAYS .and. rank == 0 ) write(*,'(a)') &
         '[mixnsolver] WARNING: uns_parallel=always needs the MPI build ' // &
         '(bin/mixnsolver_mpi); using the serial full-mesh path.'
      if ( rank == 0 ) write(*,'(a)') &
         '  Path         : full-mesh unstructured solve (serial build)'
      call run_uns( MPI_COMM_WORLD, 'unMesh.cas', 'unMesh.control' )
#endif
      if ( rank == 0 ) write(*,'(a)') &
         '=========================================================='

   case default
      ! parse_mix_mode already printed the precise failure on rank 0
      call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
   end select

   call MPI_Finalize( ierr )

contains

   !---------------------------------------------------------------------------
   ! parse_mix_mode -- read `mode = <val>` from the REQUIRED co-located
   ! mix.control and set the dispatch mode.  Deterministic on every rank (no MPI
   ! comm).  Aborts with a clear message if mix.control is missing or the mode
   ! key is absent / invalid.
   !---------------------------------------------------------------------------
   subroutine parse_mix_mode( o_mode, drank )
      implicit none
      integer,          intent(out) :: o_mode
      integer,          intent(in)  :: drank

      logical :: ex
      integer :: u, ios, ie
      character(len=512) :: line, ckey, val
      character(len=16)  :: lv

      o_mode = 0
      inquire( file='mix.control', exist=ex )
      if ( .not. ex ) then
         if ( drank == 0 ) then
            write(*,'(a)') '[mixnsolver] ERROR: mix.control is REQUIRED in the'
            write(*,'(a)') '    current directory to declare the solve mode.'
            write(*,'(a)') '    Add:  mode = coupled | struct | uns'
         end if
         call MPI_Abort( MPI_COMM_WORLD, 2, ios )
         return
      end if

      open( newunit=u, file='mix.control', status='old', action='read', iostat=ios )
      if ( ios /= 0 ) then
         if ( drank == 0 ) write(*,'(a)') '[mixnsolver] ERROR: cannot open mix.control'
         call MPI_Abort( MPI_COMM_WORLD, 3, ios )
         return
      end if

      do
         read( u, '(a)', iostat=ios ) line
         if ( ios /= 0 ) exit
         ie = scan( line, '#' )
         if ( ie > 0 ) line = line(:ie-1)
         line = adjustl( line )
         if ( len_trim(line) == 0 ) cycle
         ie = scan( line, '=' )
         if ( ie == 0 ) cycle
         ckey = trim( lowercase_str( adjustl( line(:ie-1) ) ) )
         if ( trim(ckey) /= 'mode' ) cycle
         val = adjustl( line(ie+1:) )
         ! lowercase value
         lv = lowercase_str( val )
         select case ( trim(lv) )
         case ( 'coupled', 'coupling', 'mix', '1' )
            o_mode = MODE_COUPLED
         case ( 'struct', 'structured', 'compressible', '2' )
            o_mode = MODE_STRUCT
         case ( 'uns', 'unstruct', 'unstructured', 'incompressible', '3' )
            o_mode = MODE_UNS
         case default
            if ( drank == 0 ) write(*,'(a,a,a)') '[mixnsolver] ERROR: invalid mode "', &
                  trim(val), '" in mix.control (use coupled | struct | uns)'
            call MPI_Abort( MPI_COMM_WORLD, 4, ios )
            return
         end select
      end do
      close( u )

      if ( o_mode == 0 ) then
         if ( drank == 0 ) then
            write(*,'(a)') '[mixnsolver] ERROR: mix.control has no "mode" key.'
            write(*,'(a)') '    Add:  mode = coupled | struct | uns'
         end if
         call MPI_Abort( MPI_COMM_WORLD, 5, ios )
      end if
   end subroutine parse_mix_mode

   !---------------------------------------------------------------------------
   ! parse_mix_uns_parallel -- optional `uns_parallel = auto | always | never`
   ! key in mix.control, selecting the unstructured run path in MODE_UNS:
   !
   !   auto   (default) : METIS-partitioned when nproc > 1, legacy full mesh at
   !                      -np 1 (so a single-rank run reproduces the serial
   !                      solution bit for bit)
   !   always           : always partitioned (MPI build only; also under -np 1,
   !                      which lets one compare partitioned vs legacy results
   !                      at the same rank count)
   !   never            : never partitioned (legacy full-mesh replication)
   !
   ! Deterministic on every rank (no MPI communication); a missing key or an
   ! unknown value keeps `auto` (rank 0 warns about the latter).
   !---------------------------------------------------------------------------
   subroutine parse_mix_uns_parallel( o_par, drank )
      implicit none
      integer, intent(out) :: o_par
      integer, intent(in)  :: drank

      logical :: ex
      integer :: u, ios, ie
      character(len=512) :: line, ckey, val
      character(len=16)  :: lv

      o_par = UNS_AUTO
      inquire( file='mix.control', exist=ex )
      if ( .not. ex ) return

      open( newunit=u, file='mix.control', status='old', action='read', &
            iostat=ios )
      if ( ios /= 0 ) return

      do
         read( u, '(a)', iostat=ios ) line
         if ( ios /= 0 ) exit
         ie = scan( line, '#' )
         if ( ie > 0 ) line = line(:ie-1)
         line = adjustl( line )
         if ( len_trim(line) == 0 ) cycle
         ie = scan( line, '=' )
         if ( ie == 0 ) cycle
         ckey = trim( lowercase_str( adjustl( line(:ie-1) ) ) )
         if ( trim(ckey) /= 'uns_parallel' ) cycle
         val = adjustl( line(ie+1:) )
         ! lowercase value
         lv = lowercase_str( val )
         select case ( trim(lv) )
         case ( 'auto', '0' )
            o_par = UNS_AUTO
         case ( 'always', 'true', 'yes', 'on', '1' )
            o_par = UNS_ALWAYS
         case ( 'never', 'false', 'no', 'off', '2' )
            o_par = UNS_NEVER
         case default
            if ( drank == 0 ) write(*,'(a,a,a)') &
               '[mixnsolver] WARNING: unknown uns_parallel value "', trim(val), &
               '" (use auto | always | never); using auto'
            o_par = UNS_AUTO
         end select
         exit
      end do
      close( u )
   end subroutine parse_mix_uns_parallel

   !---------------------------------------------------------------------------
   ! itoa / keyval -- tiny message-formatting helpers.  keyval prints the
   ! NS_UNSET sentinel as "(absent)" so a missing key is visible in the error.
   ! Deferred-length results so the concatenated messages carry no padding.
   !---------------------------------------------------------------------------
   function itoa( n ) result( s )
      implicit none
      integer, intent(in) :: n
      character(len=:), allocatable :: s
      character(len=24) :: tmp
      write( tmp, '(i0)' ) n
      s = trim( tmp )
   end function itoa

   function keyval( n ) result( s )
      implicit none
      integer, intent(in) :: n
      character(len=:), allocatable :: s
      character(len=24) :: tmp
      if ( n == NS_UNSET ) then
         s = '(absent)'
      else
         write( tmp, '(i0)' ) n
         s = trim( tmp )
      end if
   end function keyval

   !---------------------------------------------------------------------------
   ! lowercase_str -- ASCII lowercase copy of a string (case-insensitive keys).
   !---------------------------------------------------------------------------
   function lowercase_str( input ) result( output )
      implicit none
      character(len=*), intent(in) :: input
      character(len=len(input))    :: output
      integer :: i, ic
      output = input
      do i = 1, len(input)
         ic = iachar( input(i:i) )
         if ( ic >= iachar('A') .and. ic <= iachar('Z') ) &
            output(i:i) = achar( ic + 32 )
      end do
   end function lowercase_str

   !---------------------------------------------------------------------------
   ! parse_mix_int_key -- read an optional `<keyname> = <int>` key from
   ! mix.control (case-insensitive key, `#` comments, the value may carry a
   ! trailing comment).  Returns NS_UNSET when the file or the key is missing;
   ! an unparsable value is reported once by rank 0 and treated as absent.
   !---------------------------------------------------------------------------
   subroutine parse_mix_int_key( keyname, o_val, drank )
      implicit none
      character(len=*), intent(in)  :: keyname
      integer,          intent(out) :: o_val
      integer,          intent(in)  :: drank

      logical :: ex
      integer :: u, ios, ie, iv
      character(len=512) :: line, ckey, val
      character(len=64)  :: lkey

      o_val = NS_UNSET
      inquire( file='mix.control', exist=ex )
      if ( .not. ex ) return
      open( newunit=u, file='mix.control', status='old', action='read', &
            iostat=ios )
      if ( ios /= 0 ) return

      lkey = trim( lowercase_str( keyname ) )
      do
         read( u, '(a)', iostat=ios ) line
         if ( ios /= 0 ) exit
         ie = scan( line, '#' )
         if ( ie > 0 ) line = line(:ie-1)
         line = adjustl( line )
         if ( len_trim(line) == 0 ) cycle
         ie = scan( line, '=' )
         if ( ie == 0 ) cycle
         ckey = trim( lowercase_str( adjustl( line(:ie-1) ) ) )
         if ( trim(ckey) /= trim(lkey) ) cycle
         val = adjustl( line(ie+1:) )
         read( val, *, iostat=ios ) iv
         if ( ios /= 0 ) then
            if ( drank == 0 ) write(*,'(a,a,a)') '[mixnsolver] WARNING: ', &
               trim(keyname), ' has a non-integer value; ignoring the key'
            o_val = NS_UNSET
         else
            o_val = iv
         end if
         exit
      end do
      close( u )
   end subroutine parse_mix_int_key

   !---------------------------------------------------------------------------
   ! rank_split_error -- print the rank-count diagnosis (rank 0) and abort.
   !---------------------------------------------------------------------------
   subroutine rank_split_error( drank, code, why, ns_decl, nu_decl, nproc, total )
      implicit none
      integer,          intent(in) :: drank, code, ns_decl, nu_decl, nproc, total
      character(len=*), intent(in) :: why
      integer :: ios

      if ( drank == 0 ) then
         write(*,'(a)') '=========================================================='
         write(*,'(a)') '[mixnsolver] ERROR: ' // trim(why)
         write(*,'(a)') '----------------------------------------------------------'
         write(*,'(a)') '  mix.control   : n_struct_ranks = ' // keyval(ns_decl) // &
                        ' , n_uns_ranks = ' // keyval(nu_decl)
         write(*,'(a)') '  declared sum  : ' // keyval(total) // ' rank(s)'
         write(*,'(a)') '  launched with : ' // itoa(nproc) // &
                        ' rank(s) (MPI_COMM_WORLD)'
         write(*,'(a)') '  Fix mix.control (n_struct_ranks / n_uns_ranks) or the'
         write(*,'(a)') '  -np / job-script task count so that they agree, or run'
         write(*,'(a)') '  without these keys (legacy implicit split).'
         write(*,'(a)') '=========================================================='
      end if
      call MPI_Abort( MPI_COMM_WORLD, code, ios )
   end subroutine rank_split_error

   !---------------------------------------------------------------------------
   ! resolve_rank_split -- validate the optional mix.control rank-count keys
   ! against the launched rank count, and return the number of ranks owned by
   ! the STRUCTURED group (the MPI_Comm_split split used in COUPLED mode).
   !
   !   n_struct_ranks = <int>   ranks for the structured solver
   !   n_uns_ranks    = <int>   ranks for the unstructured solver
   !   declared sum = n_struct_ranks + n_uns_ranks, MUST equal nproc
   !
   ! Aborts (rank 0 prints the diagnosis) when only one of the two keys is
   ! present, when a count is negative, when the declared sum differs from the
   ! launched rank count, or when the declaration conflicts with the mode:
   !
   !   COUPLED : only n_struct_ranks = 1 + n_uns_ranks = 1 is supported (the
   !             group exchange connects the two group roots directly and the
   !             uns group has no intra-group collectives, so extra ranks
   !             block in the exchange); for the same reason the legacy
   !             undeclared default is only allowed up to 2 ranks.
   !   STRUCT  : n_uns_ranks must be 0 and n_struct_ranks = all ranks;
   !   UNS     : n_struct_ranks must be 0 and n_uns_ranks = all ranks.
   !
   ! With neither key present the legacy implicit split is kept (COUPLED: 1
   ! struct rank + the rest uns; STRUCT/UNS: all ranks) and an informational
   ! line is printed instead of a check.
   !---------------------------------------------------------------------------
   subroutine resolve_rank_split( mode, nproc, drank, ns_decl, nu_decl, o_ns )
      implicit none
      integer, intent(in)  :: mode, nproc, drank, ns_decl, nu_decl
      integer, intent(out) :: o_ns

      logical :: any_decl, both_decl
      integer :: total

      any_decl  = ( ns_decl /= NS_UNSET ) .or. ( nu_decl /= NS_UNSET )
      both_decl = ( ns_decl /= NS_UNSET ) .and. ( nu_decl /= NS_UNSET )

      ! ---- undeclared: legacy implicit split -------------------------------
      if ( .not. any_decl ) then
         o_ns = 1
         if ( nproc == 1 ) o_ns = 0   ! all ranks -> uns (no struct side)
         if ( mode == MODE_COUPLED .and. nproc > 2 ) call rank_split_error( &
            drank, 8, 'mode = coupled is limited to 1 struct rank + 1 uns ' // &
            'rank (total 2), which the launched rank count exceeds and no ' // &
            'n_struct_ranks / n_uns_ranks declaration was found in mix.control.', &
            ns_decl, nu_decl, nproc, NS_UNSET )
         if ( drank == 0 ) then
            select case ( mode )
            case ( MODE_COUPLED )
               write(*,'(a)') '[mixnsolver] rank split : undeclared -> legacy ' // &
                  'struct ' // itoa(o_ns) // ' + uns ' // itoa(nproc - o_ns) // &
                  '   (add n_struct_ranks + n_uns_ranks to have it checked)'
            case ( MODE_STRUCT )
               write(*,'(a)') '[mixnsolver] rank split : undeclared -> ' // &
                  'structured solver on all ' // itoa(nproc) // ' rank(s)'
            case default
               write(*,'(a)') '[mixnsolver] rank split : undeclared -> ' // &
                  'unstructured solver on all ' // itoa(nproc) // ' rank(s)'
            end select
         end if
         return
      end if

      ! ---- declared: both keys must be present together ---------------------
      if ( .not. both_decl ) call rank_split_error( drank, 6, &
         'incomplete rank-count declaration: declare BOTH n_struct_ranks and ' // &
         'n_uns_ranks (or neither).', ns_decl, nu_decl, nproc, NS_UNSET )
      if ( ns_decl < 0 .or. nu_decl < 0 ) call rank_split_error( drank, 6, &
         'negative rank count: both counts must be >= 0.', ns_decl, nu_decl, &
         nproc, NS_UNSET )

      total = ns_decl + nu_decl

      select case ( mode )
      case ( MODE_COUPLED )
         ! implementation limit first: a matching total alone cannot run
         if ( ns_decl /= 1 .or. nu_decl /= 1 ) call rank_split_error( drank, 8, &
            'unsupported split for mode = coupled: exactly 1 struct rank + 1 ' // &
            'uns rank (total 2) is supported today; for a single-solver run ' // &
            'use mode = struct / uns.', ns_decl, nu_decl, nproc, total )
         if ( total /= nproc ) call rank_split_error( drank, 7, &
            'rank-count mismatch: the declared total is not the launched ' // &
            'rank count.', ns_decl, nu_decl, nproc, total )
         o_ns = ns_decl
      case ( MODE_STRUCT )
         if ( nu_decl /= 0 ) call rank_split_error( drank, 8, &
            'mode = struct runs every rank as the structured solver, so ' // &
            'n_uns_ranks must be 0 (use mode = coupled for a struct + uns ' // &
            'split).', ns_decl, nu_decl, nproc, total )
         if ( ns_decl == 0 ) call rank_split_error( drank, 6, &
            'mode = struct needs at least 1 struct rank.', ns_decl, nu_decl, &
            nproc, total )
         if ( ns_decl /= nproc ) call rank_split_error( drank, 7, &
            'rank-count mismatch: the declared total is not the launched ' // &
            'rank count.', ns_decl, nu_decl, nproc, total )
         o_ns = ns_decl
      case ( MODE_UNS )
         if ( ns_decl /= 0 ) call rank_split_error( drank, 8, &
            'mode = uns runs every rank as the unstructured solver, so ' // &
            'n_struct_ranks must be 0 (use mode = coupled for a struct + uns ' // &
            'split).', ns_decl, nu_decl, nproc, total )
         if ( nu_decl == 0 ) call rank_split_error( drank, 6, &
            'mode = uns needs at least 1 uns rank.', ns_decl, nu_decl, nproc, &
            total )
         if ( nu_decl /= nproc ) call rank_split_error( drank, 7, &
            'rank-count mismatch: the declared total is not the launched ' // &
            'rank count.', ns_decl, nu_decl, nproc, total )
         o_ns = 0
      case default
         o_ns = 1   ! unreachable: parse_mix_mode aborts on an unknown mode
         if ( nproc == 1 ) o_ns = 0
      end select

      if ( drank == 0 ) write(*,'(a)') '[mixnsolver] rank split : ' // &
         'n_struct_ranks = ' // keyval(ns_decl) // ' + n_uns_ranks = ' // &
         keyval(nu_decl) // ' = ' // itoa(total) // ' rank(s) == launched ' // &
         itoa(nproc) // ' rank(s)   [checked]'
   end subroutine resolve_rank_split

   !---------------------------------------------------------------------------
   ! run_uns -- unstructured-only solve on the full communicator.
   ! uns_solver_init runs the FULL mesh on every rank (same as the coupling
   ! uns side, no METIS partition), so all ranks step in lockstep and the field
   ! equals the serial solve.  Only the global rank 0 writes the result file.
   !---------------------------------------------------------------------------
   subroutine run_uns( comm, casfile, ctlfile )
      use mod_uns_driver, only: uns_solver_init, uns_solver_step
      use mod_uns_mesh, only: mesh_t
      use mod_uns_connectivity, only: conn_t
      use mod_uns_geometry, only: geom_t
      use mod_uns_control, only: ctrl_t, out_is_binary, out_native_cells
      use mod_uns_bc, only: bc_t
      use mod_uns_fields, only: fields_t
      use mod_uns_output, only: vtk_write, ghia_compare, thermal_wall_report
      implicit none
      integer, intent(in) :: comm
      character(len=*), intent(in) :: casfile, ctlfile

      type(mesh_t)   :: m
      type(conn_t)   :: c
      type(geom_t)   :: g
      type(ctrl_t)   :: ctrl
      type(bc_t)     :: bcs
      type(fields_t) :: fld
      integer :: ier, ier2, j, k, kk
      real(dp) :: umax
      character(len=512) :: outstem, vtufile
      character(len=8)   :: outext

      call uns_solver_init( casfile, ctlfile, m, c, g, ctrl, bcs, fld, ier )
      if ( ier /= 0 ) then
         write(*,'(a,i0)') '[mixnsolver/uns] FATAL: init failed ier=', ier
         if ( ier == 40 ) write(*,'(a)') '  (transient controls on the ' // &
            'steady-only full-mesh path -- see the report above: use ' // &
            'bin/mixnsolver_mpi with uns_parallel = always)'
         flush(6)
         call MPI_Abort( comm, ier, ier2 )
         return
      end if

      call uns_solver_step( m, c, g, ctrl, bcs, fld, ctrl%outer_max, ier )
      if ( ier /= 0 ) then
         write(*,'(a,i0)') '[mixnsolver/uns] FATAL: solver failed ier=', ier
         if ( ier == 60 ) write(*,'(a)') '  (non-finite solution value -- ' // &
            'see the nan-guard report above; the fields printed above are the ' // &
            'last finite state)'
         flush(6)
         call MPI_Abort( comm, ier, ier2 )
         return
      end if

      call MPI_Barrier( comm, ier2 )

      if ( rank == 0 ) then
         ! ---- solution summary and thermal diagnostics ------------------------
         ! The pre-phase-15 entry point (src/unstructured/main_uns.f90) printed
         ! this block; it moved here when that file was deleted (commit
         ! da255a0), and the natural-convection heat-balance / Nusselt report
         ! (thermal_wall_report, otherwise unreferenced dead code) was lost with
         ! it.  Restored so the full-mesh path again documents the quantities
         ! the case READMEs quote.
         umax = maxval( norm2( fld%u, dim=1 ) )
         kk   = minloc( norm2( g%xc - 0.5_dp, dim=1 ), dim=1 )
         write(*,'(a)') ''
         write(*,'(a)') '--- Solution summary ---'
         write(*,'(a,es12.5)') '  max |u|         : ', umax
         write(*,'(a,3(es11.4,1x),a,es12.5)') '  cell nearest (0.5,0.5,0.5): ', &
            g%xc(:,kk), ' u = ', norm2( fld%u(:,kk) )
         write(*,'(a,3(es11.4,1x))') '                        u_vec = ', fld%u(:,kk)
         write(*,'(a,es14.6)') '  max p           : ', maxval( fld%p )
         write(*,'(a,es14.6)') '  min p           : ', minval( fld%p )
         if ( ctrl%boussinesq ) call thermal_wall_report( m, g, ctrl, bcs, fld )

         ! result file name: <cas-basename>.vtu
         outstem = 'result'
         j = scan( casfile, '/', back=.true. )
         k = scan( casfile(j+1:), '.', back=.true. )
         if ( k > 1 ) outstem = casfile(j+1:j+k-1)
         outext  = '.vtu'
         vtufile = trim(outstem)//trim(outext)
         call vtk_write( m, c, g, fld, trim(vtufile), ier, &
                         out_is_binary(ctrl), out_native_cells(ctrl) )
         if ( ier /= 0 ) then
            write(*,'(a,i0)') '[mixnsolver/uns] FATAL: result output failed ier=', ier
            call MPI_Abort( comm, ier, ier2 )
            return
         end if
         write(*,'(a,a)') '[mixnsolver/uns] wrote ', trim(vtufile)

         ! ---- benchmark comparison on the printed centreline ------------------
         ! Re = rho*U*L/mu with L = 1 (same convention as the MPI path); prints a
         ! table only for the tabulated Reynolds numbers, otherwise a skip line.
         call ghia_compare( m, g, fld, ctrl%rho / max( ctrl%mu, tiny(1.0_dp) ), ier )
      end if
   end subroutine run_uns

end program mixnsolver