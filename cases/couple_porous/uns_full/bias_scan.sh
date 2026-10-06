#!/bin/bash
# bias_scan.sh -- C2 uns 绝对压力偏置（≈ -0.176943*rho*u_in^2）的差分实验矩阵。
#
# 用法:
#   bash bias_scan.sh <workdir>                     # 收敛子集（~15 s）
#   SKIP_PROBES=0 bash bias_scan.sh <workdir>       # 加稳定性探针（预期多数发散）
#   PROBE_TIMEOUT=30 SKIP_PROBES=0 bash ...         # 探针超时（默认 60 s）
#
# 背景/结论见 ../README.md 第 5 节（§5.1 机理定位、§5.2 本脚本）。
# 所有算例都是同一"全流体"网格（gen_full.py 生成，zone 8 改 fluid），
# 只改 出口/入口 BC、inlet_ramp、以及 x 分辨率；逐平面 p/u 由
# plane_profile.py 的 load() 读取（uns 写出器把每个单元质心点放在单元块的第 5
# 个节点，故不需要任何平均，也不会把流体/多孔界面抹平）。
#
# 产物（*.vtu/*.log）都落在 <workdir>；仓库里不跟踪运行产物。
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
BIN=$REPO/bin/uns_solver
WORK=${1:-.}
SKIP_PROBES=${SKIP_PROBES:-1}
PROBE_TIMEOUT=${PROBE_TIMEOUT:-60}
BASE=$HERE/unMesh.control

[ -x "$BIN" ] || { echo "ERROR: missing $BIN -- run 'make unstructured' first"; exit 1; }
[ -f "$BASE" ] || { echo "ERROR: missing $BASE"; exit 1; }

mkdir -p "$WORK" && cd "$WORK"
cp "$HERE/gen_full.py" "$HERE/plane_profile.py" .
if [ ! -f unMesh.cas ]; then
   python3 gen_full.py > gen.log 2>&1 || { echo "ERROR: gen_full.py failed"; exit 1; }
fi

# 基准 = 全流体 + pressure-outlet 0.0 + mass-flow-inlet（其余取 unMesh.control）
sed 's/^cell_zone = 8 .*/cell_zone = 8 fluid/' "$BASE" > base.control

run() {   # name controlfile [timeout]
   local n=$1 c=$2 t=${3:-90}
   printf '  %-10s ' "$n"
   timeout "$t" "$BIN" unMesh.cas "$c" "$n.vtu" > "$n.log" 2>&1
   local rc=$?
   if grep -q CONVERGED "$n.log"; then
      printf 'converged   %s\n' "$(grep -E 'max p|min p' "$n.log" | tr '\n' ' ' | tr -s ' ')"
   else
      printf 'NOT converged (rc=%d)\n' "$rc"
   fi
}

echo "--- 收敛子集（Δx = 2 mm, 全流体）---"
cp base.control c1.control                                    # baseline
sed 's/^inlet_ramp = .*/inlet_ramp = 10/' base.control > c2.control
sed 's/^bc = 7 .*/bc = 7 mass-flow-inlet 20.434 300.0/' base.control > c3.control
run c1 c1.control
run c2 c2.control
run c3 c3.control

echo "--- x 分辨率（Δx = 1 / 4 mm）---"
for tag in 200 50; do
   d="nx$tag"
   mkdir -p "$d"
   if [ "$tag" = "50" ]; then
      sed -e 's/^NX, NY, NZ = 100, 20, 2/NX, NY, NZ = 50, 20, 2/' \
          -e 's/^NF_HALF = 50/NF_HALF = 25/' gen_full.py > "$d/gen_full.py"
   else
      sed 's/^NX, NY, NZ = 100, 20, 2/NX, NY, NZ = 200, 20, 2/' gen_full.py > "$d/gen_full.py"
   fi
   cp plane_profile.py base.control "$d/"
   ( cd "$d" && { [ -f unMesh.cas ] || python3 gen_full.py > gen.log 2>&1; } )
   printf '  %-10s ' "$d"
   ( cd "$d" && timeout 150 "$BIN" unMesh.cas base.control out.vtu > out.log 2>&1 )
   if grep -q CONVERGED "$d/out.log"; then
      printf 'converged   %s\n' "$(grep -E 'max p|min p' "$d/out.log" | tr '\n' ' ' | tr -s ' ')"
   else
      printf 'NOT converged\n'
   fi
done

if [ "$SKIP_PROBES" = "0" ]; then
   echo "--- 稳定性探针（见 README §5.1；预期多数不收敛，每条上限 ${PROBE_TIMEOUT}s）---"
   sed 's/^bc = 6 .*/bc = 6 pressure-outlet 10.0/'  base.control > p1.control
   sed 's/^bc = 6 .*/bc = 6 pressure-outlet -10.0/' base.control > p2.control
   sed 's/^inlet_ramp = .*/inlet_ramp = 1000/'      base.control > p3.control
   sed 's/^bc = 7 .*/bc = 7 velocity-inlet 34.7224 0.0 0.0 300.0/' base.control > p4.control
   run p1 p1.control "$PROBE_TIMEOUT"
   run p2 p2.control "$PROBE_TIMEOUT"
   run p3 p3.control "$PROBE_TIMEOUT"
   run p4 p4.control "$PROBE_TIMEOUT"
fi

# ---------------------------------------------------------------------------
# 汇总工具（生成到工作目录，随用随看）
# ---------------------------------------------------------------------------
cat > sum.py <<'PY'
#!/usr/bin/env python3
"""逐平面平台值 / 入口邻格 / 凹陷汇总：sum.py <vtu> ..."""
import sys, numpy as np, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plane_profile import load
hdr = ('%-16s %5s %8s %15s %11s %17s %11s %17s %16s'
       % ('case', 'npl', 'dx[mm]', 'plateau p[Pa]', 'std', 'plateau u[m/s]',
          'p_1[Pa]', 'u_1/u_in', 'dip p[Pa]@mm'))
print(hdr)
for f in sys.argv[1:]:
    if not os.path.exists(f):
        continue
    x, p, u = load(f)
    xr = np.round(x, 9)
    pl = np.array(sorted(set(xr)))
    mp = np.array([p[xr == v].mean() for v in pl])
    mu = np.array([u[xr == v].mean() for v in pl])
    mid = (pl > 0.30 * pl[-1]) & (pl < 0.85 * pl[-1])
    print('%-16s %5d %8.2f %15.5f %11.2e %17.5f %11.5f %17.6f %12.3f @%.1f'
          % (f, pl.size, (pl[1] - pl[0]) * 1e3, mp[mid].mean(), mp[mid].std(),
             mu[mid].mean(), mp[0], mu[0] / mu[mid].mean(),
             mp.min(), pl[mp.argmin()] * 1e3))
PY

cat > shift.py <<'PY'
#!/usr/bin/env python3
"""床组（组 1）整场平移：measured p - 物理解(流体半段 302.13 平、床半段 3021.7 Pa/m)。
用法: python3 shift.py g1.vtu   （g1.vtu = 用 unMesh.control 原样跑出的组 1）"""
import sys, numpy as np, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plane_profile import load
x, p, u = load(sys.argv[1])
xr = np.round(x, 9)
pl = np.array(sorted(set(xr)))
mp = np.array([p[xr == v].mean() for v in pl])
mu = np.array([u[xr == v].mean() for v in pl])
Pb, G = 302.13, 3021.7                       # 解析：界面 302.13 Pa，床梯度 3021.7 Pa/m
phys = np.where(pl < 0.1, Pb, Pb - G * (pl - 0.1))
sh = mp - phys
fl = (pl > 0.02) & (pl < 0.098)
bd = (pl > 0.102) & (pl < 0.198)
print('%s: 流体半段 mean p = %.4f (物理解 %.2f)  床半段 mean p = %.4f'
      % (sys.argv[1], mp[fl].mean(), Pb, mp[bd].mean()))
print('  整场平移（剔除两端单元）: mean %.4f  std %.2e Pa' % (sh[fl | bd].mean(), sh[fl | bd].std()))
for k in list(range(0, 4)) + [24, 49, 50, 74] + list(range(pl.size - 3, pl.size)):
    print('    x=%8.3f mm  measured=%11.4f  physical=%10.4f  shift=%10.4f  u=%9.5f'
          % (pl[k] * 1e3, mp[k], phys[k], sh[k], mu[k]))
PY

echo "--- 平台/入口邻格汇总（c1 baseline / c2 ramp=10 / c3 u-halved / nx200 / nx50）---"
python3 sum.py c1.vtu c2.vtu c3.vtu nx200/out.vtu nx50/out.vtu
if [ -f g1.vtu ]; then
   echo "--- 床组整场平移（g1.vtu）---"
   python3 shift.py g1.vtu
else
   echo "（床组平移：cp $BASE g1.control && $BIN unMesh.cas g1.control g1.vtu && python3 shift.py g1.vtu）"
fi
