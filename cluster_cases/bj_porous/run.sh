#!/usr/bin/env bash
# =========================================================================
# Beavers-Joseph (fluid / porous interface) scene runner
#   ./run.sh [np] [repo-root]
#     np        number of MPI ranks.  Default: the rank split declared in
#               mix.control (n_struct_ranks + n_uns_ranks), or 1 when the scene
#               declares neither.  np>1 uses bin/mixnsolver_mpi; in the
#               unstructured mode np>1 runs the METIS-partitioned parallel path
#               (mix.control: uns_parallel = auto, the default), i.e. the cells
#               really are split over the ranks.  Add "uns_parallel = never" to
#               mix.control for the old full-mesh-per-rank replication
#               (determinism check only, no speed-up, bit-identical to serial).
#               NOTE: when mix.control declares n_struct_ranks / n_uns_ranks,
#               np MUST equal their sum -- the solver checks this before the
#               MPI rank split and aborts with a printed diagnosis (exit code
#               6/7/8) instead of running on the wrong number of ranks.
#     repo-root path holding bin/ (default $MIXROOT or ../..)
#
# Runs the three slip coefficients alpha = 0 / 1 / 2 by swapping the active
# control out of ./variants (the working copy is always unMesh.control).
# Artifacts land in ./logs and ./out, so the scene stays pristine.
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
BIN=$ROOT/bin/mixnsolver
[ "$NP" -gt 1 ] && BIN=$ROOT/bin/mixnsolver_mpi
mkdir -p logs out
# Unbuffered solver output: libgfortran does its own buffering, so a redirected
# solver log stays EMPTY for minutes (piping through stdbuf does NOT help).
# This keeps logs/*.log usable with `tail -f` while the run proceeds.
export GFORTRAN_UNBUFFERED_ALL=y

echo "== bj_porous: open channel + Darcy bed, alpha = 0 / 1 / 2, np=$NP =="
echo "   bin      : $BIN"
for A in a0 a1 a2; do
   cp -f "variants/bj_${A}.control" unMesh.control
   echo "-- alpha case bj_${A} (control = variants/bj_${A}.control)"
   if [ "$NP" -le 1 ]; then
      "$BIN" 2>&1 | tee "logs/bj_${A}_np$NP.log"
   else
      mpirun ${MPIRUN_EXTRA:-} -np "$NP" "$BIN" 2>&1 | tee "logs/bj_${A}_np$NP.log"
   fi
   mv -f unMesh.vtu "out/bj_${A}_np$NP.vtu"
done
# leave the alpha = 1 control active (matches the archived bj_a1 case)
cp -f variants/bj_a1.control unMesh.control

echo
echo "wrote  out/bj_a{0,1,2}_np$NP.vtu   (logs: logs/bj_a*_np$NP.log)"
echo "verify with:"
echo "  BJ_PLOT_OUT=out/bj_check_np$NP.png python3 plot_bj.py \\"
echo "      0.0=out/bj_a0_np$NP.vtu 1.0=out/bj_a1_np$NP.vtu 2.0=out/bj_a2_np$NP.vtu"
