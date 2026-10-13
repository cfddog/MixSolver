#!/usr/bin/env bash
# =========================================================================
# Cylinder scene runner (2D laminar cylinder, Re = 200, UNSTEADY / PIMPLE)
#   ./run.sh [np] [repo-root]
#     np        number of MPI ranks.  Default: the rank split declared in
#               mix.control (n_struct_ranks + n_uns_ranks), or 1 when the
#               scene declares neither.
#     repo-root path holding bin/ (default $MIXROOT or ../..)
#
# Unlike the steady scenes this one ALWAYS uses bin/mixnsolver_mpi: PIMPLE is
# implemented only in the MPI driver (mod_uns_driver_mpi -> pimple_run_mpi),
# and `uns_parallel = always` in mix.control makes the METIS-partitioned path
# active even at np = 1.  Build it with:  make mpi -j1
#
# Artifacts land in ./logs and ./out, so the scene stays pristine:
#   out/<tag>.vtu               final field (== last snapshot)
#   out/<tag>.step<NNNN>.vtu    one snapshot per n_out_every time steps
#   out/<tag>.cp.dat            wall Cp distribution at the final state
#   logs/<tag>.log              full solver log
# =========================================================================
set -euo pipefail
cd "$(dirname "$0")"
# np default: the rank split declared in mix.control (n_struct_ranks + n_uns_ranks),
# 1 when the scene declares neither.  A DECLARED split is authoritative: the
# solver itself re-checks it against the launched rank count and aborts (exit
# code 6/7/8) on any mismatch, so a wrong np fails loudly instead of silently.
MIX_NP=$(awk -F= '/^[[:space:]]*n_(struct|uns)_ranks[[:space:]]*=/ { v=$2; gsub(/[^0-9]/,"",v); s+=v } END { print s+0 }' mix.control 2>/dev/null || true)
MIX_NP=${MIX_NP:-0}
NP=${1:-${MIX_NP:-1}}
if [ "$NP" -lt 1 ]; then NP=1; fi
if [ "$MIX_NP" -gt 0 ] && [ "$MIX_NP" -ne "$NP" ]; then
   echo "!! mix.control declares $MIX_NP rank(s) (n_struct_ranks+n_uns_ranks) but np=$NP was requested"
   echo "!! the solver will abort with a rank-count mismatch -- use ./run.sh $MIX_NP"
fi
ROOT=${2:-${MIXROOT:-../..}}
BIN=$ROOT/bin/mixnsolver_mpi
if [ ! -x "$BIN" ]; then
   echo "!! $BIN not found or not executable -- build it with: make mpi -j1" >&2
   exit 1
fi
mkdir -p logs out
# Unbuffered solver output: libgfortran does its own buffering, so a redirected
# solver log stays EMPTY for minutes (piping through stdbuf does NOT help).
# This keeps logs/*.log usable with `tail -f` while the run proceeds.
export GFORTRAN_UNBUFFERED_ALL=y
TAG=cylinder_Re200_np$NP

echo "== cylinder_Re200: U_inf=1, D=1, Re=200 (UNSTEADY PIMPLE), np=$NP =="
echo "   bin      : $BIN"
mpirun ${MPIRUN_EXTRA:-} -np "$NP" "$BIN" 2>&1 | tee logs/$TAG.log

# --- archive the transient output ------------------------------------------
# PIMPLE writes the final field to unMesh.vtu and one snapshot every
# n_out_every steps to unMesh.vtu_step<7-digit>.vtu; move both into out/.
mv -f unMesh.vtu "out/$TAG.vtu"
n_snap=0
for f in unMesh.vtu_step*.vtu; do
   [ -e "$f" ] || continue
   # strip the leading "unMesh.vtu_" -> "<TAG>.step<NNNNNNN>.vtu"
   mv -f "$f" "out/${TAG}.${f#unMesh.vtu_}"
   n_snap=$((n_snap + 1))
done
if [ -f unMesh_cp.dat ]; then mv -f unMesh_cp.dat "out/$TAG.cp.dat"; fi

echo
echo "wrote  out/$TAG.vtu              (final field)"
echo "       out/$TAG.step*.vtu        ($n_snap snapshots)"
echo "wrote  logs/$TAG.log"
echo "verify/plot with:"
echo "  python3 check_cylinder_unsteady.py 0.005 out/$TAG.step*.vtu"
