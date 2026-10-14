#!/usr/bin/env bash
# ============================================================================
# run_all.sh -- 串行运行 cases/heat_conduction 的 5 个纯导热验证算例
#
# 遵守「本机资源保护」硬规则（见 .trae/rules/project_rules.md）：
#   * 同一时刻只跑 1 个算例（本脚本严格串行，不并行提交）；
#   * 每个算例 wall-clock 上限 600 s（timeout 600，可用 TIMEOUT= 覆盖）；
#   * 运行现场（unMesh.cas / unMesh.control / *.vtu / run.log）全部落在
#     .dbg/heat_conduction/<case>/ 下 —— 均为可重算的运行产物，不入库。
#
# 用法：
#   bash cases/heat_conduction/run_all.sh            # 跑全部 5 个
#   bash cases/heat_conduction/run_all.sh c1 c3      # 只跑指定算例
#   bash cases/heat_conduction/run_all.sh c3f        # 附加：c3 的径向加密网格
#                                                    #（Nr=40，网格现场生成、不入库）
#   TIMEOUT=90 bash cases/heat_conduction/run_all.sh c1   # 缩短单例上限
#
# 结果：每个算例打印 run.log 尾部，并提示 .vtu 路径；对拍用
#       python3 check_heat_conduction.py [--refine]
#
# 网格（*.cas）已入库，如需重新生成（在 cases/heat_conduction/ 下）：
#   python3 gen_1d.py      c1_1d.cas
#   python3 gen_1d.py      c5_1d_src.cas 40 --src --ny 1   # --src 必须 NY=1
#   python3 gen_annulus.py c3_annulus.cas                  # 默认 Nr=20 Nt=72
# ============================================================================
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
BIN=$ROOT/bin/mixnsolver
RUN=$ROOT/.dbg/heat_conduction
TIMEOUT=${TIMEOUT:-600}

CASES=${*:-"c1 c2 c3 c4 c5"}

declare -A MESH=(
  [c1]=c1_1d.cas       [c2]=c1_1d.cas
  [c3]=c3_annulus.cas  [c4]=c3_annulus.cas
  [c5]=c5_1d_src.cas )
declare -A CTL=(
  [c1]=c1_1d_fixedT.control
  [c2]=c2_1d_adiabatic.control
  [c3]=c3_annulus_fixedT.control
  [c4]=c4_annulus_adiabatic.control
  [c5]=c5_1d_source.control
  [c3f]=c3_annulus_fixedT.control )


if [ ! -x "$BIN" ]; then
  echo "ERROR: $BIN not found -- build first (make -j)" >&2
  exit 1
fi

# ---- 网格：缺失则按 README 记录的命令重新生成（可复现） --------------------
if [ ! -f "$HERE/c1_1d.cas" ] || [ ! -f "$HERE/c5_1d_src.cas" ]; then
  ( cd "$HERE" && python3 gen_1d.py c1_1d.cas \
             && python3 gen_1d.py c5_1d_src.cas 40 --src --ny 1 ) || exit 1
fi
if [ ! -f "$HERE/c3_annulus.cas" ]; then
  ( cd "$HERE" && python3 gen_annulus.py c3_annulus.cas ) || exit 1
fi

FAIL=0
for c in $CASES; do
  d=$RUN/$c
  mkdir -p "$d"
  case $c in
    # c3f = c3 的径向加密版（Nr=40, Nt=72）。mesh 在 .dbg 里现场生成（不入库），
    # 控制文件复用 c3 的（两者 zone id 完全一致：2 固体 / 4 symmetry z /
    # 5 内壁 / 6 外壁）。用于网格收敛性检查：python3 check_heat_conduction.py --refine
    c3f) meshshow="annulus Nr=${C3F_NR:-40} Nt=72 (现场生成, 不入库)"
         python3 "$HERE/gen_annulus.py" "$d/unMesh.cas" "${C3F_NR:-40}" 72 ;;
    *)   meshshow="${MESH[$c]}"
         cp "$HERE/${MESH[$c]}" "$d/unMesh.cas" ;;
  esac
  cp "$HERE/${CTL[$c]}"  "$d/unMesh.control"
  cp "$HERE/mix.control" "$d/mix.control"
  echo "=================================================================="
  echo "[$c] $(date '+%F %T')  mesh=${meshshow}  ctl=${CTL[$c]}  timeout=${TIMEOUT}s"
  ( cd "$d" && timeout "$TIMEOUT" "$BIN" > run.log 2>&1 )
  rc=$?
  if [ $rc -eq 124 ]; then
    echo "[$c] *** TIMEOUT after ${TIMEOUT}s (本机资源保护上限) ***"
    FAIL=1
  elif [ $rc -ne 0 ]; then
    echo "[$c] *** exit code $rc ***"
    FAIL=1
  fi
  echo "----- run.log tail -----"
  tail -n 25 "$d/run.log"
  if [ -f "$d/unMesh.vtu" ]; then
    echo "[$c] OK  -> $d/unMesh.vtu"
  else
    echo "[$c] *** no unMesh.vtu produced ***"
    FAIL=1
  fi
done

echo "=================================================================="
if [ $FAIL -eq 0 ]; then
  echo "all cases finished; now: python3 $HERE/check_heat_conduction.py"
else
  echo "some cases failed or timed out (see above)"
fi
exit $FAIL
