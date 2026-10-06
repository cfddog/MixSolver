!===============================================================================
! mod_cas_reader.f90 -- ASCII Fluent .cas/.msh mesh reader
!
! Handles the sections required for face-based unstructured FVM:
!   (2  ndim)                                    problem dimension
!   (10 (zid first last act ndim)(coords))       node coordinates
!   (12 (zid first last act ctype)(payload))     cells
!   (13 (zid first last cond ftype)(payload))    faces with owner/neighbor
!   (45 (zid condname username)())               zone condition names
!
! IMPORTANT: Fluent writes all integer indices in these sections in
! HEXADECIMAL notation. All index parsing therefore uses the Z format.
!===============================================================================
module mod_uns_cas_reader
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   implicit none
   private
   public :: read_cas

   ! --- tokenizer state ---
   character(len=1), allocatable :: buf(:)   ! whole file as characters
   integer :: ipos = 0                       ! current read position
   integer :: ilen = 0                       ! buffer length
   integer, parameter :: TOKLEN = 96         ! max token length

   ! --- fluent zone condition code -> name lookup (fallback only) ---
   integer, parameter :: NCOD = 16
   integer :: cod( NCOD ) = (/ 2, 3, 5, 7, 9, 12, 14, 15, 17, 18, 20, &
                               24, 31, 36, 37, 38 /)
   character(len=18) :: codname( NCOD ) = (/ &
      'interior          ', 'wall              ', 'pressure-outlet   ', &
      'symmetry          ', 'velocity-inlet    ', 'pressure-inlet    ', &
      'inlet-vent        ', 'intake-fan        ', 'outlet-vent       ', &
      'exhaust-fan       ', 'mass-flow-inlet   ', 'outflow           ', &
      'axis              ', 'pressure-far-field', 'fluid             ', &
      'solid             ' /)

contains

   !----------------------------------------------------------------------------
   ! Main entry: read a Fluent ASCII .cas/.msh file into mesh_t
   !----------------------------------------------------------------------------
   subroutine read_cas( filename, m, ier )
      character(len=*), intent(in)  :: filename
      type(mesh_t),     intent(out) :: m
      integer,          intent(out) :: ier

      character(len=TOKLEN) :: t
      integer               :: ios

      ier = 0
      call slurp( filename, ier )
      if ( ier /= 0 ) return

      ! main token loop: sections are identified by '(' followed by the
      ! section number
      do
         call next_token( t, ios )
         if ( ios /= 0 ) exit              ! EOF -> done
         if ( trim(t) /= '(' ) cycle       ! ignore stray tokens
         call next_token( t, ios )
         if ( ios /= 0 ) exit
         select case ( trim(t) )
         case ( '2' )
            call sec_dim( m )
         case ( '10' )
            call sec_nodes( m, ier )
         case ( '12' )
            call sec_cells( m, ier )
         case ( '13' )
            call sec_faces( m, ier )
         case ( '45' )
            call sec_zonedef( m )
         case default
            call skip_record()             ! unknown/comment sections
         end select
         if ( ier /= 0 ) return
      end do

      call finalize_zones( m )
      call validate_mesh( m, ier )

   end subroutine read_cas

   !----------------------------------------------------------------------------
   ! Slurp the whole file into the module buffer
   !----------------------------------------------------------------------------
   subroutine slurp( filename, ier )
      character(len=*), intent(in)  :: filename
      integer,          intent(out) :: ier
      integer :: u, sz, ios

      ier = 0
      inquire( file = trim(filename), size = sz )
      if ( sz <= 0 ) then
         write(*,'(a)') 'ERROR: cannot read file (or empty): ' // trim(filename)
         ier = 1
         return
      end if
      allocate( buf(sz) )
      open( newunit = u, file = trim(filename), access = 'stream', &
            status = 'old', action = 'read', iostat = ios )
      if ( ios /= 0 ) then
         write(*,'(a)') 'ERROR: cannot open file: ' // trim(filename)
         ier = 1
         return
      end if
      read( u, iostat = ios ) buf
      close( u )
      ilen = sz
      ipos = 1
   end subroutine slurp

   !----------------------------------------------------------------------------
   ! Get the next whitespace-delimited token; '(' and ')' are always
   ! returned as separate single-character tokens.
   ! ier: 0 = ok, -1 = EOF
   !----------------------------------------------------------------------------
   subroutine next_token( t, ier )
      character(len=TOKLEN), intent(out) :: t
      integer,               intent(out) :: ier
      character(len=1) :: c
      integer :: n

      t = repeat( ' ', TOKLEN ); ier = 0

      do                                        ! skip whitespace
         if ( ipos > ilen ) then; ier = -1; return; end if
         c = buf(ipos)
         if ( c == ' ' .or. c == achar(9) .or. c == achar(10) .or. &
              c == achar(13) ) then
            ipos = ipos + 1
         else
            exit
         end if
      end do

      n = 0
      do while ( ipos <= ilen )
         c = buf(ipos)
         if ( c == ' ' .or. c == achar(9) .or. c == achar(10) .or. &
              c == achar(13) ) exit
         if ( c == '(' .or. c == ')' ) then
            if ( n == 0 ) then                 ! paren starts its own token
               t(1:1) = c
               n = 1
               ipos = ipos + 1
            end if
            exit                               ! paren also ends a token
         end if
         n = n + 1
         if ( n > TOKLEN ) then; ier = 2; return; end if
         t(n:n) = c
         ipos = ipos + 1
      end do

   end subroutine next_token

   !----------------------------------------------------------------------------
   ! Read a balanced "(...)" group; the current token must be '('.
   ! Returns the tokens strictly inside the group.
   ! ier: 0 = ok, 3 = current token is not '(' (token consumed!), -1 = EOF
   !----------------------------------------------------------------------------
   subroutine read_group( g, ng, ier )
      character(len=TOKLEN), allocatable, intent(out) :: g(:)
      integer,               intent(out) :: ng
      integer,               intent(out) :: ier
      character(len=TOKLEN) :: t
      integer :: ios, depth, cap

      ng = 0; ier = 0
      allocate( g(64) ); cap = 64

      call next_token( t, ios )                  ! must be '('
      if ( ios /= 0 ) then; ier = -1; return; end if
      if ( trim(t) /= '(' ) then; ier = 3; return; end if

      depth = 1
      do
         call next_token( t, ios )
         if ( ios /= 0 ) then; ier = -1; return; end if
         if ( trim(t) == '(' ) then
            depth = depth + 1
         else if ( trim(t) == ')' ) then
            depth = depth - 1
            if ( depth == 0 ) exit
         end if
         ng = ng + 1
         if ( ng > cap ) call grow( g, cap )
         g(ng) = t
      end do

   end subroutine read_group

   subroutine grow( g, cap )
      character(len=TOKLEN), allocatable, intent(inout) :: g(:)
      integer, intent(inout) :: cap
      character(len=TOKLEN), allocatable :: tmp(:)
      allocate( tmp(cap*2) )
      tmp(1:cap) = g(1:cap)
      call move_alloc( tmp, g )
      cap = cap * 2
   end subroutine grow

   ! Skip tokens until the currently open record is balanced again.
   ! Called right after a group has been consumed: depth starts at 1.
   subroutine skip_record()
      character(len=TOKLEN) :: t
      integer :: ios, depth
      depth = 1
      do
         call next_token( t, ios )
         if ( ios /= 0 ) return
         if ( trim(t) == '(' ) depth = depth + 1
         if ( trim(t) == ')' ) then
            depth = depth - 1
            if ( depth == 0 ) return
         end if
      end do
   end subroutine skip_record

   !----------------------------------------------------------------------------
   ! Integer parse: Fluent indices are hexadecimal (Z edit descriptor),
   ! with a plain-decimal fallback.
   !----------------------------------------------------------------------------
   integer function tok_int( t, ier )
      character(len=*), intent(in)  :: t
      integer,          intent(out) :: ier
      integer :: lt, ios
      character(len=16) :: fmt

      lt = len_trim(t)
      write( fmt, '(a,i0,a)' ) '(Z', lt, ')'
      read( t(1:lt), fmt, iostat = ios ) tok_int
      ier = 0
      if ( ios /= 0 ) then
         read( t(1:lt), *, iostat = ios ) tok_int
         if ( ios /= 0 ) ier = 4
      end if
   end function tok_int

   real(dp) function tok_real( t, ier )
      character(len=*), intent(in)  :: t
      integer,          intent(out) :: ier
      integer :: ios
      read( t(1:len_trim(t)), *, iostat = ios ) tok_real
      ier = 0
      if ( ios /= 0 ) ier = 5
   end function tok_real

   !----------------------------------------------------------------------------
   ! Section (2 ndim)
   !----------------------------------------------------------------------------
   subroutine sec_dim( m )
      type(mesh_t), intent(inout) :: m
      character(len=TOKLEN) :: t
      integer :: ios
      call next_token( t, ios )
      if ( ios == 0 ) read( t(1:len_trim(t)), *, iostat = ios ) m%ndim
   end subroutine sec_dim

   !----------------------------------------------------------------------------
   ! Section (10 (zid first last act [ndim])(coords))
   !----------------------------------------------------------------------------
   subroutine sec_nodes( m, ier )
      type(mesh_t), intent(inout) :: m
      integer,      intent(out)   :: ier

      character(len=TOKLEN), allocatable :: g(:), c(:)
      integer :: ng, nh(8), k, n, i, ier2
      integer :: zid, first, last, act, nd

      ier = 0
      call read_group( g, ng, ier )
      if ( ier /= 0 ) return
      if ( ng < 4 ) then; ier = 6; return; end if

      do k = 1, min(ng, 8)
         nh(k) = tok_int( g(k), ier2 )
         if ( ier2 /= 0 ) then; ier = 6; return; end if
      end do
      zid = nh(1); first = nh(2); last = nh(3); act = nh(4)
      nd = 2
      if ( ng >= 5 ) nd = nh(5)
      m%nnodes = max( m%nnodes, last )

      if ( zid == 0 .or. act == 0 ) then         ! header record
         call skip_record()
         return
      end if

      if ( .not. allocated(m%x) ) allocate( m%x(3, m%nnodes) )

      call read_group( c, ng, ier )              ! coordinate payload
      if ( ier /= 0 ) return
      n = last - first + 1
      if ( ng < n*nd ) then; ier = 7; return; end if

      do i = 1, n
         do k = 1, nd
            m%x(k, first+i-1) = tok_real( c((i-1)*nd + k), ier2 )
            if ( ier2 /= 0 ) then; ier = 7; return; end if
         end do
      end do
      call skip_record()                         ! consume closing ')'
   end subroutine sec_nodes

   !----------------------------------------------------------------------------
   ! Section (12 (zid first last act ctype)(payload))
   ! Data records: act /= 0. For mixed cells (ctype=0) the payload holds one
   ! element type per cell; polyhedra (7) list face counts and face indices.
   ! Definition records: (zid first last cond name) -- no payload.
   !----------------------------------------------------------------------------
   subroutine sec_cells( m, ier )
      type(mesh_t), intent(inout) :: m
      integer,      intent(out)   :: ier

      character(len=TOKLEN), allocatable :: g(:)
      integer :: ng, nh(8), k, ier2, p, ct, nfa, ic
      integer :: zid, first, last, act

      ier = 0
      call read_group( g, ng, ier )
      if ( ier /= 0 ) then
         if ( ier == 3 ) ier = 0                 ! record ended: nothing more
         return
      end if
      if ( ng < 4 ) then; ier = 8; return; end if

      do k = 1, min(ng, 8)
         nh(k) = tok_int( g(k), ier2 )
         if ( ier2 /= 0 ) then
            ! non-numeric header field -> definition record with name
            call skip_record()
            return
         end if
      end do
      zid = nh(1); first = nh(2); last = nh(3); act = nh(4)

      if ( zid == 0 ) then                       ! header record
         call skip_record()
         return
      end if

      if ( ng == 4 .or. ng == 5 ) then
         ! could be a definition record (cond + name): distinguish by
         ! checking whether the 5th field parses as an integer (element type)
         if ( ng == 4 ) then                     ! definition without elem type
            call skip_record()
            return
         end if
         ct = tok_int( g(5), ier2 )
         if ( ier2 /= 0 ) then                   ! name -> definition record
            call skip_record()
            return
         end if
      else
         ct = -1
      end if

      m%ncells = max( m%ncells, last )
      call grow_iarray( m%ctype, m%ncells )
      if ( ct > 0 ) m%ctype(first:last) = ct

      ! phase 3 step C: record the owning cell (volume) zone of every cell.
      ! Cell-zone names arrive later via the (45 records; finalize_zones
      ! splits cell zones out of the common zone table.
      call grow_iarray( m%czone, m%ncells )
      m%czone(first:last) = zid

      ! payload group (present for mixed cells); if next token is ')' the
      ! record is already complete
      call read_group( g, ng, ier )
      if ( ier == 3 ) then                       ! no payload
         ier = 0
         return
      else if ( ier /= 0 ) then
         return
      end if

      if ( ct > 0 ) then
         ! fixed element type with (unexpected) payload: ignore payload
         call skip_record()
         return
      end if

      ! mixed cells: parse per-cell records
      ic = first
      p = 1
      do while ( p <= ng )
         ct = tok_int( g(p), ier2 ); p = p + 1
         if ( ier2 /= 0 ) then; ier = 8; return; end if
         if ( ct == 7 ) then                     ! polyhedron: [7 nf f1..fnf]
            if ( p > ng ) then; ier = 8; return; end if
            nfa = tok_int( g(p), ier2 ); p = p + 1
            if ( ier2 /= 0 ) then; ier = 8; return; end if
            p = p + nfa
            if ( p-1 > ng ) then; ier = 8; return; end if
         end if
         if ( ic > m%ncells ) then; ier = 8; return; end if
         m%ctype(ic) = ct
         ic = ic + 1
      end do
      call skip_record()                         ! consume closing ')'
   end subroutine sec_cells

   !----------------------------------------------------------------------------
   ! Section (13 (zid first last cond ftype)(payload))
   ! Data records carry ftype as 5th field; definition records carry a name.
   !----------------------------------------------------------------------------
   subroutine sec_faces( m, ier )
      type(mesh_t), intent(inout) :: m
      integer,      intent(out)   :: ier

      character(len=TOKLEN), allocatable :: g(:)
      type(face_t), allocatable :: fl(:)
      type(zone_t) :: z
      integer :: ng, nh(8), k, ier2, p, ft, nn, j, zi
      integer :: zid, first, last, cond, ftype, nfl

      ier = 0
      call read_group( g, ng, ier )
      if ( ier /= 0 ) return
      if ( ng < 4 ) then; ier = 9; return; end if

      do k = 1, min(ng, 8)
         nh(k) = tok_int( g(k), ier2 )
         if ( ier2 /= 0 ) then                   ! name field -> definition
            call skip_record()
            return
         end if
      end do
      zid = nh(1); first = nh(2); last = nh(3); cond = nh(4)

      if ( zid == 0 .or. cond == 0 ) then        ! header record
         call skip_record()
         return
      end if

      if ( ng == 4 ) then                        ! definition (no ftype here
         z%id = zid; z%cond_code = cond          ! in this file layout)
         z%nf = last - first + 1
         call zone_put( m, z )
         call skip_record()
         return
      end if

      ftype = tok_int( g(5), ier2 )
      if ( ier2 /= 0 ) then                      ! 5th field is a name
         z%id = zid; z%cond_code = cond
         z%nf = last - first + 1
         call zone_put( m, z )
         call skip_record()
         return
      end if

      nfl = last - first + 1
      call read_group( g, ng, ier )              ! face payload
      if ( ier /= 0 ) return

      allocate( fl(nfl) )
      p = 0
      do k = 1, nfl
         if ( ftype == 0 ) then                  ! mixed: leading face type
            p = p + 1
            if ( p > ng ) then; ier = 10; return; end if
            ft = tok_int( g(p), ier2 )
            if ( ier2 /= 0 ) then; ier = 10; return; end if
         else
            ft = ftype
         end if

         select case ( ft )
         case ( 2, 3, 4 )
            nn = ft                              ! code == node count:
                                                 ! line(2)/tri(3)/quad(4)
         case ( 5 )                              ! polygon: count follows
            p = p + 1
            if ( p > ng ) then; ier = 10; return; end if
            nn = tok_int( g(p), ier2 )
            if ( ier2 /= 0 ) then; ier = 10; return; end if
         case default
            ier = 11
            write(*,'(a,i0,a,i0,a,i0,a)') 'ERROR: unknown face type ', ft, &
               ' at face #', k, ' (token pos ', p, '): ' // trim(g(p))
            write(*,'(a,20(1x,a))') 'DEBUG first tokens:', &
               ('"'//trim(g(j))//'"', j = 1, min(ng,20))
            return
         end select

         allocate( fl(k)%nodes(nn) )
         fl(k)%nn = nn
         do j = 1, nn
            p = p + 1
            if ( p > ng ) then; ier = 10; return; end if
            fl(k)%nodes(j) = tok_int( g(p), ier2 )
            if ( ier2 /= 0 ) then; ier = 10; return; end if
         end do
         p = p + 1
         if ( p+1 > ng ) then; ier = 10; return; end if
         fl(k)%c0 = tok_int( g(p), ier2 )
         fl(k)%c1 = tok_int( g(p+1), ier2 )
         if ( ier2 /= 0 ) then; ier = 10; return; end if
         p = p + 1
         fl(k)%zone = zid
      end do

      call mesh_append_faces( m, fl, first, nfl )

      z%id = zid; z%cond_code = cond; z%nf = nfl
      zi = zone_find( m, zid )
      if ( zi > 0 ) then
         m%zone(zi)%cond_code = cond
         m%zone(zi)%nf = nfl
      else
         call zone_put( m, z )
      end if

      deallocate( fl )
      call skip_record()                         ! consume closing ')'
   end subroutine sec_faces

   !----------------------------------------------------------------------------
   ! Section (45 (zid condname username)())   -- zone condition names
   ! (also used for polyhedral cell-face tables; numeric records are skipped)
   !----------------------------------------------------------------------------
   subroutine sec_zonedef( m )
      type(mesh_t), intent(inout) :: m

      character(len=TOKLEN), allocatable :: g(:)
      type(zone_t) :: z
      integer :: ng, ier2, zi, zid

      call read_group( g, ng, ier2 )
      if ( ier2 /= 0 ) then
         if ( ier2 == 3 ) return                 ! record already closed
         call skip_record()
         return
      end if

      if ( ng >= 3 ) then
         zid = tok_int( g(1), ier2 )
         if ( ier2 == 0 ) then
            block
               integer :: probe
               ! 2nd/3rd fields must be non-numeric for a zone definition
               probe = tok_int( g(2), ier2 )
               if ( ier2 /= 0 ) then             ! g(2) is a name -> definition
                  z%id = zid
                  z%cond_name = adjustl( g(2) )
                  z%user_name = adjustl( g(3) )
                  zi = zone_find( m, zid )
                  if ( zi > 0 ) then
                     m%zone(zi)%cond_name = z%cond_name
                     m%zone(zi)%user_name = z%user_name
                  else
                     call zone_put( m, z )
                  end if
               end if
            end block
         end if
      end if
      call skip_record()                         ! consume remaining groups
   end subroutine sec_zonedef

   ! --- zone table helpers ------------------------------------------------------
   integer function zone_find( m, zid )
      type(mesh_t), intent(in) :: m
      integer, intent(in) :: zid
      integer :: i
      zone_find = 0
      if ( .not. allocated(m%zone) ) return
      do i = 1, m%nzone
         if ( m%zone(i)%id == zid ) then
            zone_find = i
            return
         end if
      end do
   end function zone_find

   subroutine zone_put( m, z )
      type(mesh_t), intent(inout) :: m
      type(zone_t), intent(in)    :: z
      type(zone_t), allocatable :: tmp(:)
      if ( .not. allocated(m%zone) ) allocate( m%zone(8) )
      if ( m%nzone == size(m%zone) ) then
         allocate( tmp(m%nzone*2) )
         tmp(1:m%nzone) = m%zone(1:m%nzone)
         call move_alloc( tmp, m%zone )
      end if
      m%nzone = m%nzone + 1
      m%zone(m%nzone) = z
   end subroutine zone_put

   ! fill missing condition names from the numeric code table, then split
   ! the common zone table into face zones (m%zone) and cell zones (m%czt).
   ! A zone id referenced by faces is a face zone; one referenced by cells
   ! is a cell zone.  If an id appears on both (unusual), it is kept in
   ! both tables.
   subroutine finalize_zones( m )
      type(mesh_t), intent(inout) :: m
      integer :: i, j, id, nc
      type(zone_t), allocatable :: fz(:)
      integer :: nfz

      ! ---- numeric condition-name fallback (existing behaviour) --------------
      do i = 1, m%nzone
         if ( len_trim(m%zone(i)%cond_name) == 0 ) then
            do j = 1, NCOD
               if ( m%zone(i)%cond_code == cod(j) ) then
                  m%zone(i)%cond_name = trim( codname(j) )
                  exit
               end if
            end do
         end if
      end do

      ! ---- make sure the per-cell zone map exists and covers every cell ------
      if ( .not. allocated(m%czone) ) then
         allocate( m%czone( max(m%ncells,1) ), source = 0 )
      else if ( size(m%czone) < m%ncells ) then
         call grow_iarray( m%czone, m%ncells )
      end if

      ! ---- build the cell-zone table from distinct czone ids -----------------
      do id = minval( m%czone(1:m%ncells) ), maxval( m%czone(1:m%ncells) )
         if ( id < 1 ) cycle                  ! 0 = uncovered cells
         nc = count( m%czone(1:m%ncells) == id )
         if ( nc == 0 ) cycle
         call czone_put( m, id, nc )
      end do

      ! ---- keep only face-referenced zones in m%zone -------------------------
      allocate( fz( max(m%nzone,1) ) )
      nfz = 0
      do i = 1, m%nzone
         if ( any( m%f(1:m%nfaces)%zone == m%zone(i)%id ) ) then
            nfz = nfz + 1
            fz(nfz) = m%zone(i)
         end if
      end do
      if ( allocated(m%zone) ) deallocate( m%zone )
      m%nzone = nfz
      if ( nfz > 0 ) then
         allocate( m%zone(nfz) )
         m%zone(1:nfz) = fz(1:nfz)
      end if
   end subroutine finalize_zones

   ! append/complete one cell-zone entry; names are taken from the (45 record
   ! if present in the temporary common zone table, otherwise left blank
   subroutine czone_put( m, id, nc )
      type(mesh_t), intent(inout) :: m
      integer,      intent(in)    :: id, nc
      type(zone_t), allocatable :: tmp(:)
      integer :: k

      do k = 1, m%nczone
         if ( m%czt(k)%id == id ) then
            m%czt(k)%nf = nc
            return
         end if
      end do

      if ( .not. allocated(m%czt) ) allocate( m%czt(4) )
      if ( m%nczone == size(m%czt) ) then
         allocate( tmp(m%nczone*2) )
         tmp(1:m%nczone) = m%czt(1:m%nczone)
         call move_alloc( tmp, m%czt )
      end if
      m%nczone = m%nczone + 1
      m%czt(m%nczone)%id   = id
      m%czt(m%nczone)%nf   = nc
      ! pick up condition/user names supplied by a (45 record
      do k = 1, m%nzone
         if ( m%zone(k)%id == id ) then
            m%czt(m%nczone)%cond_code = m%zone(k)%cond_code
            m%czt(m%nczone)%cond_name = m%zone(k)%cond_name
            m%czt(m%nczone)%user_name = m%zone(k)%user_name
            exit
         end if
      end do
      if ( len_trim(m%czt(m%nczone)%cond_name) == 0 ) then
         do k = 1, NCOD
            if ( m%czt(m%nczone)%cond_code == cod(k) ) then
               m%czt(m%nczone)%cond_name = trim( codname(k) )
               exit
            end if
         end do
      end if
   end subroutine czone_put

   ! grow an allocatable 1-D integer array to at least n elements; the newly
   ! added tail is zero-filled so index ranges written by the caller start
   ! from a defined state
   subroutine grow_iarray( a, n )
      integer, allocatable, intent(inout) :: a(:)
      integer,              intent(in)    :: n
      integer, allocatable :: t(:)
      if ( allocated(a) ) then
         if ( size(a) >= n ) return
         allocate( t(n), source = 0 )
         t(1:size(a)) = a
         call move_alloc( t, a )
      else
         allocate( a(n), source = 0 )
      end if
   end subroutine grow_iarray

   ! --- mesh face storage helper ------------------------------------------------
   subroutine mesh_append_faces( m, fl, first, nfl )
      type(mesh_t),  intent(inout) :: m
      type(face_t),  intent(in)    :: fl(:)
      integer,       intent(in)    :: first, nfl
      type(face_t), allocatable :: tmp(:)

      if ( .not. allocated(m%f) ) then
         allocate( m%f( max(nfl,1024) ) )
      else if ( size(m%f) < first + nfl - 1 ) then
         allocate( tmp( max( 2*size(m%f), first+nfl-1 ) ) )
         tmp(1:size(m%f)) = m%f
         call move_alloc( tmp, m%f )
      end if
      m%f(first:first+nfl-1) = fl(1:nfl)
      m%nfaces = max( m%nfaces, first + nfl - 1 )
   end subroutine mesh_append_faces

   !----------------------------------------------------------------------------
   ! Final consistency validation of the parsed mesh
   !----------------------------------------------------------------------------
   subroutine validate_mesh( m, ier )
      type(mesh_t), intent(in)  :: m
      integer,      intent(out) :: ier
      integer :: i, j, badn, badc, badz

      ier = 0
      if ( m%nnodes == 0 .or. m%ncells == 0 .or. m%nfaces == 0 ) then
         write(*,'(a)') 'ERROR: incomplete mesh sections in .cas file'
         ier = 20
         return
      end if

      badn = 0; badc = 0; badz = 0
      do i = 1, m%nfaces
         do j = 1, m%f(i)%nn
            if ( m%f(i)%nodes(j) < 1 .or. m%f(i)%nodes(j) > m%nnodes ) &
               badn = badn + 1
         end do
         if ( m%f(i)%c0 < 1 .or. m%f(i)%c0 > m%ncells ) badc = badc + 1
         if ( m%f(i)%c1 < 0 .or. m%f(i)%c1 > m%ncells ) badc = badc + 1
         if ( zone_find( m, m%f(i)%zone ) == 0 ) badz = badz + 1
      end do

      if ( badn > 0 .or. badc > 0 ) then
         write(*,'(a,i0,a,i0)') 'ERROR: invalid node refs: ', badn, &
                                ', invalid cell refs: ', badc
         write(*,'(a)') '       -> index base mismatch (hex/decimal)?'
         ier = 21
      end if
      if ( badz > 0 ) then
         write(*,'(a,i0)') 'WARNING: faces referencing unknown zones: ', badz
      end if
   end subroutine validate_mesh

end module mod_uns_cas_reader
