#!/usr/bin/env bash
#===============================================================================
# M6-wing bit-level regression for the structured solver.
#
# Compares the in-tree structured path (bin/mixnsolver, ser tree, mode=struct)
# against the pristine OpenCFD-EC 1.16a reference (external/OpenCFD-EC-1.16a),
# run np1 on the 4-block M6-wing case (t_end=0.501, Kstep_save=50 -> 51 steps,
# save at 50).
#
# The pristine reference reproduces the historical /tmp/m6reg baseline
# bit-for-bit: flow3d.dat md5 = dc134a2d196422043ecad7c86ac8f898.  Prebuilt
# baseline outputs live in baseline/, so by default the new solver is checked
# against those files without rebuilding the reference.  Set M6_REBUILD_REF=1 to
# instead rebuild + rerun the pristine reference from source (iconv recipe).
#
# Env knobs:
#   M6_WORK        scratch dir           (default /tmp/m6wing_regress)
#   M6_NP          ranks                 (default 1)
#   M6_OPT         reference build flags (default "-O2 -std=legacy -ffree-line-length-none")
#   M6_REBUILD_REF 1 = rebuild pristine reference instead of using baseline/
#===============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
EXT_SRC="$ROOT/external/OpenCFD-EC-1.16a"
CASE_SRC="$EXT_SRC/cases/M6-wing"
CTL="${M6_CTL:-$HERE/control.ec}"
WORK="${M6_WORK:-/tmp/m6wing_regress}"
NP="${M6_NP:-1}"
OPT="${M6_OPT:--O2 -std=legacy -ffree-line-length-none}"

BINS=(flow3d.dat SA3d.dat wall_dist.dat partation-auto.dat part_grid.dat)
TEXTS=(Step_mess.dat bc3d.inc mesh-quality.dat)

echo "== [1/4] build in-tree structured solver =="
make -C "$ROOT" all -j"${M6_J:-4}"

echo "== [2/4] prepare case dir (mixnsolver reads Mesh3d.x, mode=struct) =="
rm -rf "$WORK"; mkdir -p "$WORK/new"
cp "$CASE_SRC/bc3d.inp" "$CASE_SRC/Mesh3d.dat" "$CTL" "$WORK/new/"
mv "$WORK/new/Mesh3d.dat" "$WORK/new/Mesh3d.x"
# control.ec shipped here is already the t_end=0.501 / Kstep_save=50 regression
# control.  mixnsolver requires a mix.control declaring the zero-arg solve mode.
printf 'mode = struct\n' > "$WORK/new/mix.control"

if [[ "${M6_REBUILD_REF:-0}" == "1" ]]; then
   echo "== [2b] build pristine reference (GBK->UTF-8, fullwidth ! -> !) =="
   REF="$WORK/ref_src"
   cp -r "$EXT_SRC" "$REF"
   ( cd "$REF"
     for f in *.f90; do
        if iconv -f GBK -t UTF-8 "$f" > "$f.u8" 2>/dev/null; then mv "$f.u8" "$f"; else rm -f "$f.u8"; fi
        sed -i 's/\xEF\xBC\x81/!/g' "$f"
     done
     make f77=mpif90 opt="$OPT" > build.log 2>&1 )
   REF_BIN="$REF/opencfd-ec1.16a.out"
   mkdir -p "$WORK/ref"
   cp "$CASE_SRC/bc3d.inp" "$CASE_SRC/Mesh3d.dat" "$CTL" "$WORK/ref/"
   echo "== [3/4] run pristine reference (np$NP) =="
   ( cd "$WORK/ref" && mpirun --oversubscribe -np "$NP" "$REF_BIN" > run.log 2>&1 || true )
   BASE="$WORK/ref"
else
   echo "== [3/4] compare against prebuilt baseline/ (set M6_REBUILD_REF=1 to rebuild) =="
   BASE="$HERE/baseline"
fi

echo "== [3b] run in-tree solver (np$NP) =="
( cd "$WORK/new" && mpirun --oversubscribe -np "$NP" "$ROOT/bin/mixnsolver" > run.log 2>&1 )

echo "== [4/4] compare =="
fail=0
for f in "${BINS[@]}"; do
   if cmp -s "$BASE/$f" "$WORK/new/$f"; then echo "  IDENTICAL  $f"; else echo "  DIFFER     $f"; fail=1; fi
done
for f in "${TEXTS[@]}"; do
   if diff -q "$BASE/$f" "$WORK/new/$f" >/dev/null; then echo "  IDENTICAL  $f"; else echo "  DIFFER     $f"; fail=1; fi
done

echo "== flow3d.dat md5 =="
( cd "$WORK/new" && md5sum flow3d.dat )
echo "   expected: dc134a2d196422043ecad7c86ac8f898"

if [[ "$fail" == 0 ]]; then echo "M6-WING REGRESSION: PASS"; else echo "M6-WING REGRESSION: FAIL"; fi
exit "$fail"
