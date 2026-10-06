!===============================================================================
! mod_metis_iface.f90 -- Fortran 2008 iso_c_binding interface for METIS v5
!
! Serial METIS graph/mesh partitioning.  Rank 0 builds the dual graph (from
! conn_t%c2c CSR), calls METIS_PartGraphKway with nparts = nprocs, then
! broadcasts the partition map.  This is the domain-decomposition entry point
! for the parallel solver (Step 2 of the MPI plan).
!
! Type assumptions (validated by tests/test_metis.f90):
!   idx_t  = c_int32_t  (32-bit signed int)
!   real_t = c_float     (32-bit IEEE float)
!   If METIS is ever rebuilt with IDXTYPEWIDTH=64 or REALTYPEWIDTH=64,
!   change the idx_t / real_t parameter below and rerun the smoke test.
!
! METIS v5 returns METIS_OK = 1 on success (not 0), METIS_ERROR = -2 on
! failure.
!===============================================================================
module mod_uns_metis_iface
   use iso_c_binding, only: c_int, c_float, c_int32_t, c_ptr, c_null_ptr
   implicit none
   private
   public :: metis_partkway, metis_setdefaultoptions, metis_ok, metis_error

   integer, parameter :: idx_t  = c_int32_t    ! METIS index type
   integer, parameter :: real_t = c_float      ! METIS real type

   ! METIS return codes (from metis.h)
   integer(c_int), parameter :: metis_ok    =  1_c_int
   integer(c_int), parameter :: metis_error = -2_c_int

   ! METIS option indices (v5 metis.h enum) -- documented for reference.
   ! PTYPE=0, OBJTYPE=1, CTYPE=2, IATYPE=3, RTYPE=4, DBGLEVEL=5,
   ! DBGFILE=6, NIPARTS=7, NITER=8, NCUTS=9, SEED=10, ONDISK=11,
   ! COMPRESS=12, CCORDER=13, CITYPE=14, CONTIG=15, MINCONN=16,
   ! NUMBERING=17, DROPEDGES=18, NOOUTPUT=19, NOPTIONS=40.

   interface
      !-----------------------------------------------------------------------
      ! int METIS_PartGraphKway(
      !   idx_t *nvtxs,  idx_t *ncon,  idx_t *xadj,  idx_t *adjncy,
      !   idx_t *vwgt,   idx_t *vsize, idx_t *adjwgt,
      !   idx_t *nparts, real_t *tpwgts, real_t *ubvec,
      !   idx_t *options, idx_t *objval, idx_t *part);
      !
      ! vwgt, vsize, adjwgt, tpwgts, ubvec, options may be NULL in C.
      ! In Fortran we pass c_null_ptr for the unused ones.
      !-----------------------------------------------------------------------
      integer(c_int) function metis_partkway( &
           nvtxs, ncon, xadj, adjncy, &
           vwgt, vsize, adjwgt, &
           nparts, tpwgts, ubvec, &
           options, objval, part ) &
           bind(c, name="METIS_PartGraphKway")
         import :: c_int, c_int32_t, c_float, c_ptr
         integer(c_int32_t), intent(in)    :: nvtxs
         integer(c_int32_t), intent(in)    :: ncon
         integer(c_int32_t), intent(in)    :: xadj(*)
         integer(c_int32_t), intent(in)    :: adjncy(*)
         type(c_ptr),        value         :: vwgt
         type(c_ptr),        value         :: vsize
         type(c_ptr),        value         :: adjwgt
         integer(c_int32_t), intent(in)    :: nparts
         type(c_ptr),        value         :: tpwgts
         type(c_ptr),        value         :: ubvec
         type(c_ptr),        value         :: options
         integer(c_int32_t), intent(out)   :: objval
         integer(c_int32_t), intent(out)   :: part(*)
      end function metis_partkway

      !-----------------------------------------------------------------------
      ! int METIS_SetDefaultOptions(idx_t *options);
      ! Fills options[0..METIS_NOPTIONS-1] with default values.  options must
      ! be at least 40 elements.  Returns METIS_OK on success.
      !-----------------------------------------------------------------------
      integer(c_int) function metis_setdefaultoptions( options ) &
           bind(c, name="METIS_SetDefaultOptions")
         import :: c_int, c_int32_t
         integer(c_int32_t), intent(out) :: options(*)
      end function metis_setdefaultoptions
   end interface

end module mod_uns_metis_iface
