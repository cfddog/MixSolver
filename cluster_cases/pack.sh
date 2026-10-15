#!/usr/bin/env bash
# =========================================================================
# Pack the cluster scenes for transfer to a cluster login node.
#
#   ./pack.sh [out.tar.gz] [--with-repo]
#     out.tar.gz   default: cluster_cases_<YYYYmmdd>.tar.gz
#     --with-repo  also write <out>_repo.tar.gz: a WHITELIST of the source
#                  tree (Makefile, src/, lib/, docs/ + the scenes) so the
#                  cluster needs only a Fortran toolchain -- no git clone and
#                  no 1.2 GB cases/ checkout. Measured ~4.0 MB / 126 entries.
#
# The build reads nothing outside src/ and lib/metis/libmetis.a (see Makefile
# rules 'src/%.f90' and METIS_LIB), so the whitelist is sufficient to produce
# bin/mixnsolver and bin/mixnsolver_mpi.
#
# Run artifacts (logs/, out/, *.vtu) and build debris are excluded either way,
# so the archives always contain pristine scenes plus a sha256 for transfer.
# =========================================================================
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)

WITH_REPO=0
ARGS=()
for a in "$@"; do
   case "$a" in
      --with-repo) WITH_REPO=1 ;;
      *)           ARGS+=("$a") ;;
   esac
done
OUT=${ARGS[0]:-$HERE/cluster_cases_$(date +%Y%m%d).tar.gz}

EXCL=('*/logs' '*/out' '*.vtu' '*/__pycache__' '*.pyc')

echo "== scenes   : $HERE"
tar "${EXCL[@]/#/--exclude=}" -czf "$OUT" -C "$ROOT" cluster_cases
printf '   %-8s %s  (%s entries)\n' "$(du -h "$OUT" | cut -f1)" "$OUT" \
       "$(tar -tzf "$OUT" | wc -l)"

if [ "$WITH_REPO" -eq 1 ]; then
   REPO_OUT="${OUT%.tar.gz}_repo.tar.gz"
   echo "== repo tree: $ROOT (whitelist: Makefile src lib docs cluster_cases)"
   WHITE=(Makefile src lib docs cluster_cases)
   KEEP=()
   for w in "${WHITE[@]}"; do if [ -e "$ROOT/$w" ]; then KEEP+=("$w"); fi; done
   if [ ${#KEEP[@]} -eq 0 ]; then
      echo "!! no source-tree entries found under $ROOT -- skipping" >&2
   else
      tar "${EXCL[@]/#/--exclude=}" --transform "s,^,$(basename "$ROOT")/," \
          -czf "$REPO_OUT" -C "$ROOT" "${KEEP[@]}"
      printf '   %-8s %s\n' "$(du -h "$REPO_OUT" | cut -f1)" "$REPO_OUT"
   fi
fi

echo "== sha256"
sha256sum "$OUT"; [ "$WITH_REPO" -eq 1 ] && sha256sum "${OUT%.tar.gz}_repo.tar.gz"

cat <<'EOT'

on the cluster:
  scenes only :  tar -xzf <out>.tar.gz       -> cluster_cases/<scene>/
                 (drop it into an existing repo, or use your own bin/)
  with repo   :  tar -xzf <out>_repo.tar.gz  -> MixNSSolver/
                 cd MixNSSolver && make all -j1                 # -> bin/mixnsolver
                 make mpi -j1                                   # -> bin/mixnsolver_mpi
                 cd cluster_cases/<scene> && MIXROOT=$PWD/../.. ./run.sh

verify before/after transfer:  sha256sum -c <(sha256sum <out>.tar.gz)
EOT
