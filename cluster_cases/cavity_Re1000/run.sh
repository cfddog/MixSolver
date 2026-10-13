#!/usr/bin/env bash
# =========================================================================
# Lid-driven cavity scene runner (2D cavity, Re = 1000, steady SIMPLE)
#   ./run.sh [np] [repo-root]
#     np        number of MPI ranks.  Default: the rank split declared in
#               mix.control (n_struct_ranks + n_uns_ranks), or 1 when the
#               scene declares neither.
#     repo-root path holding bin/ (default $MIXROOT or ../..)
#
# This scene uses bin/mixnsolver_mpi; `uns_parallel = always` in mix.control
# makes the METIS-partitioned SIMPLE path active even at np = 1.
#
# The MPI path also runs the rank-0 post-processing on the gathered field, so
# the log ends with the Ghia et al. (1982) Re = 1000 comparison table; this
# scene is judged by that table (expect relative L2 ~ 3e-2 on this mesh).
#
# Artifacts land in ./logs and ./out, so the scene stays pristine:
#   out/<tag>.vtu               final field (= the steady state)
#   out/<tag>.step<NNNN>.vtu    one snapshot per n_out_every time steps
#   logs/<tag>.log              full solver log incl. the Ghia table
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
TAG=cavity_Re1000_SIMPLE10_np$NP

echo "== cavity_Re1000: lid = 1, L = 1, Re = 1000 (SIMPLE, max 10 iterations), np=$NP =="
echo "   bin      : $BIN"
mpirun ${MPIRUN_EXTRA:-} -np "$NP" "$BIN" 2>&1 | tee logs/$TAG.log

# --- archive the solver output ---------------------------------------------
# SIMPLE writes the final field to unMesh.vtu; retain any snapshots if enabled.
n_snap=0
if [ -f unMesh.vtu ]; then mv -f unMesh.vtu "out/$TAG.vtu"; fi
for f in unMesh.vtu_step*.vtu; do
   [ -e "$f" ] || continue
   # strip the leading "unMesh.vtu_" -> "<TAG>.step<NNNNNNN>.vtu"
   mv -f "$f" "out/${TAG}.${f#unMesh.vtu_}"
   n_snap=$((n_snap + 1))
done
if [ -f unMesh_cp.dat ]; then mv -f unMesh_cp.dat "out/$TAG.cp.dat"; fi

echo
echo "wrote  out/$TAG.vtu              (field after at most 10 SIMPLE iterations)"
echo "       out/$TAG.step*.vtu        ($n_snap snapshots)"
echo "wrote  logs/$TAG.log             (SIMPLE iteration log and Ghia comparison)"
echo "verify/plot with:"
echo "  python3 check_cavity.py out/$TAG.vtu"
