#!/bin/bash
# A/B test of the fluid/porous interface-pressure fix (kink_face_pressure) and the
# Darcy-sink mu/K convention: build a pre-fix binary (source git-stashed) and a
# post-fix one, run the SAME cases with both, then compare column/field data.
#
#   hir   Betchen PLUG Re_H=1000 (the only case whose pre-fix output was not kept)
#   dae2  PLUG Da=1e-2          bj2/bj3  Beavers-Joseph Da=1e-2 / 1e-3
#   dae3  PLUG Da=1e-3          fluid    same mesh with NO porous zone (no-op test)
#   ppd   1-D porous plug, Darcy (isolates mu/K vs (mu/eps)/K)
#
# usage: bash abtest_interface_pressure.sh [scratch-dir]     (default /tmp/plugtest)
#
# Requirements: the fix must still be an *uncommitted* worktree change of
#   src/unstructured/mod_uns_fields.f90 + src/unstructured/mod_uns_simple.f90
# (that is what gets stashed to produce the pre-fix binary), and a clean build
# tree (make unstructured -j4 must work).  Aborts (and pops the stash) on failure.
set -u
REPO=$(cd "$(dirname "$0")/../.." && pwd)
PT=${1:-/tmp/plugtest}
FILES="src/unstructured/mod_uns_fields.f90 src/unstructured/mod_uns_simple.f90"
mkdir -p "$PT"
exec > "$PT/abtest.log" 2>&1
set -x
cd "$REPO" || exit 1

git diff --quiet -- $FILES && {
  echo "ABTEST_ABORT: $FILES have no uncommitted changes -> nothing to stash."; exit 1; }
git stash push -- $FILES || { echo PRE_STASH_FAIL; exit 1; }
make unstructured -j4 || { echo PRE_BUILD_FAIL; git stash pop; exit 1; }
cp bin/uns_solver "$PT/uns_solver_PREFIX"; md5sum "$PT/uns_solver_PREFIX"

# ---------------- case scratch dirs ----------------
for k in hir dae2 dae3 bj2 bj3 fluid ppd; do
  for TAG in PRE POST; do mkdir -p "$PT/${k}${TAG}"; done
done
cp cases/betchen/plug_hir.cas  cases/betchen/plug_hir.control  "$PT/hirPRE/"
cp cases/betchen/plug_dae2.cas cases/betchen/plug_dae2.control "$PT/dae2PRE/"
cp cases/betchen/plug_dae3.cas cases/betchen/plug_dae3.control "$PT/dae3PRE/"
cp cases/betchen/bj_dae2.cas   cases/betchen/bj_dae2.control   "$PT/bj2PRE/"
cp cases/betchen/bj_dae3.cas   cases/betchen/bj_dae3.control   "$PT/bj3PRE/"
cp cases/porous_plug/plug.cas  cases/porous_plug/plug_darcy.control "$PT/ppdPRE/"
# fluid = the PLUG Da=1e-3 mesh with the cell_zone line removed (zone then
# defaults to fluid): identical geometry, no porous cells -> the fix must be a
# bit-level no-op here (VTU and solver log md5 must match).
cp cases/betchen/plug_dae3.cas "$PT/fluidPRE/"
grep -v '^[[:space:]]*cell_zone' cases/betchen/plug_dae3.control \
  > "$PT/fluidPRE/plug_dae3.control"
for k in hir dae2 dae3 bj2 bj3 fluid ppd; do
  cp "$PT/${k}PRE"/* "$PT/${k}POST/"
done

declare -A CAS CTL
CAS[hir]=plug_hir.cas;    CTL[hir]=plug_hir.control
CAS[dae2]=plug_dae2.cas;  CTL[dae2]=plug_dae2.control
CAS[dae3]=plug_dae3.cas;  CTL[dae3]=plug_dae3.control
CAS[bj2]=bj_dae2.cas;     CTL[bj2]=bj_dae2.control
CAS[bj3]=bj_dae3.cas;     CTL[bj3]=bj_dae3.control
CAS[fluid]=plug_dae3.cas; CTL[fluid]=plug_dae3.control
CAS[ppd]=plug.cas;        CTL[ppd]=plug_darcy.control

runblock() {   # runblock TAG EXE
  local TAG=$1 EXE=$2
  for k in hir dae2 dae3 bj2 bj3 fluid ppd; do
    echo "### RUN $TAG $k"
    ( cd "$PT/${k}${TAG}" && date && $EXE "${CAS[$k]}" "${CTL[$k]}" > run.log 2>&1 \
        && date ) || echo "RUN_FAIL $TAG $k"
  done
}
runblock PRE "$PT/uns_solver_PREFIX"

# ---------------- restore + rebuild ----------------
cd "$REPO" || exit 1
git stash pop || { echo POST_POP_FAIL; exit 1; }
make unstructured -j4 || { echo POST_BUILD_FAIL; exit 1; }
cp bin/uns_solver "$PT/uns_solver_FIXED"
md5sum "$PT/uns_solver_FIXED"

runblock POST "$PT/uns_solver_FIXED"

echo ABTEST_DONE
git stash list
git status --porcelain

# ---------------- quick report ----------------
echo '--- fluid no-op check (PRE vs POST: all four md5 must match) ---'
md5sum "$PT"/fluid{PRE,POST}/plug_dae3.vtu "$PT"/fluid{PRE,POST}/run.log
echo '--- finer analysis: ---'
echo "  python3 cmp_centerline.py $PT/dae3PRE/plug_dae3.vtu=PRE  ..."
echo "  python3 vs_ref.py <vtu>=<ucsv>=<pcsv>=<label> ..."
echo "  python3 fdiff_vtu.py $PT/dae2PRE/plug_dae2.vtu $PT/dae2POST/plug_dae2.vtu label"
echo "  ../porous_plug/fit_pp.py $PT/ppdPOST/plug.vtu"
