# C2 可压缩（struct）– 多孔（uns）跨组界面验证

> 日期：2026-10-06　求解器：`bin/uns_solver`、`bin/mixsolver(_mpi)`
> 主题：结构侧流体段（x=0–100 mm）+ 非结构侧多孔床（x=100–200 mm）弱耦合，
> 界面在 x=100 mm；用"单求解器全流域算例 `uns_full`"作参考对拍。
> 本轮修了两个缺陷：① `uns_full` 网格生成器的 **cell-id 排序**（把串联写成了
> 并联，梯度只有解析值 ~60%）；② `compare_iface.py` 的**取窗**（跨两个平面，
> 幻影 112 Pa 跳变）。

## 1. 布局

| 角色 | 文件 | 说明 |
|---|---|---|
| 结构侧（耦合） | `Mesh3d.x`、`bc3d.inp`、`control.ec`、`mix.control` | x=0–100 mm 纯流体；四侧 symmetry（严格 1D）；i− 强制入口、i+ 界面 |
| 非结构侧（耦合） | `unMesh.cas`、`unMesh.control` | x=100–200 mm，**全域 = 床**（单 cell zone `2`）；x=200 pressure-outlet、x=100 界面 |
| 参考（单求解器） | `uns_full/{gen_full.py, unMesh.cas, unMesh.control, unMesh.vtu, full.log}` | x=0–200 mm 用**一个** uns 求解器：zone 2 fluid（0–100 mm）＋ zone 8 `VC:porous`（100–200 mm） |
| 界面比较 | `compare_iface.py` | 界面处 struct 面／uns 首排／参考三方对比 |
| 床梯度校验 | `uns_full/check_bed_gradient.py` | 由 VTU 逐格心还原 p(x)，校验床梯度＝解析值（**同时是回归守卫**） |
| 逐平面剖析 | `uns_full/plane_profile.py` | 打印 pmax/pmin 格位置、逐 x 平面 p/u 均值与拟合 dP/dx（12 子四面体 / 5 节点两种 VTU 布局通用）——§5 四组对照的取证工具 |
| 失效网格留档 | `scratch/unMesh.cas.zsplit_bug` | 修复前的 z-分片网格（复现旧结果用） |

物理：300 K 空气 ρ=1.177 kg/m³、μ=1.846e-5 Pa·s；入口 u=34.7224 m/s
（`mass flux m'` = ρu = 40.868 kg/m²/s）。

## 2. 解析参考（床 x=100–200 mm）

`unMesh.control`：`cell_zone = 8 porous perm=1e-6 inertial=1.0 porosity=0.4`，故

    dP/dx = μu/(εK) + ρβu² = 1602.2 + 1419.5 = 3021.7 Pa/m = 3.0217 Pa/mm
    Δp(床) = 302.1 Pa  (L = 0.1 m)

流体半段无阻力、四侧 symmetry ⇒ 严格 1D、无粘无损失 ⇒ 应有 dP/dx = 0。

## 3. 缺陷 ①（已修）：`gen_full.py` 的 cell-id 排序把"串联"写成"并联"

**症状**：`uns_full` 复现出 dP/dx ≈ **1.78–1.83 Pa/mm**——流体半段 −1829.2 Pa/m、
床半段 −1775.5 Pa/m，只有解析值的 ~60%；`max|u|` = 42.87 m/s（＞入口 34.72）、
`max p` = 398.06 Pa（＞床压降 302 Pa）。把 zone 2／zone 8 互换后结果**逐位相同**。

**根因**：Fluent 的 cell-zone 记录是**一段连续 id 区间**，而 `gen_full.py` 原式

    cid = 1 + i + j*NX + k*NX*NY        # 旧：k(z) 最慢

让 k 成为最慢索引，于是 `zone 2 = cells 1..2000` 落在 **z=2.5 mm 的整长条**、
`zone 8 = 2001..4000` 落在 **z=7.5 mm 的整长条**——两个 zone 都是 x=0–200 mm
**全长**的 z 薄片，在 z 方向**并联**，而不是 x 方向**串联**的
"流体 100 mm ＋ 床 100 mm"。

压力场是单值 p(x)，两个并联片只能共享同一个 dP/dx：流体片"被拖下水"、
床片"被抬起来"，量出 ~1.8 Pa/mm 的折中梯度。zone 号互换只是 z 镜像
（几何对称）⇒ 逐位相同，这曾是"换 zone 无影响"这一假象的来源。

**修法**（`gen_full.py`，`cid()` 上方已留说明注释）：

    cid = 1 + k + j*NZ + i*NY*NZ        # 新：x 最慢，再 z，再 y

于是 zone 2 → 1..2000（i=0..49，x=0–100 mm）、zone 8 → 2001..4000
（i=50..99，x=100–200 mm）各为一段连续区间。

**验证**（`check_bed_gradient.py`，逐格心 p(x)）：

| 量 | 修复前（z 并联） | 修复后（x 串联） | 解析 |
|---|---|---|---|
| 床 dP/dx | −1775.5 Pa/m | **−3021.3 Pa/m**（≈−3.0213 Pa/mm） | −3021.7 Pa/m |
| 流体半段 dP/dx | −1829.2 Pa/m | 14.5 Pa/m | 0 |
| 床压降（拟合×0.1 m） | — | **302.13 Pa（偏差 0.01%）** | 302.1 Pa |
| `max p` | 398.06 Pa | 301.88 Pa | ≈302 |
| `max\|u\|` | 42.87 m/s | 35.23 m/s | 34.72（入口） |

`check_bed_gradient.py` 现 **PASS**；在 z-分片旧网格上它会**主动报错**
（`AssertionError: zone 2 is not the x<100 mm half!`），即该脚本兼作回归守卫。

## 4. 缺陷 ②（已修）：`compare_iface.py` 的取窗跨越了两个平面

**症状**：界面处 struct 与 uns 的**压力差 112 Pa**（struct 295.0 / uns 182.6），
而两侧物理上应当连续。

**根因**：取窗原用"中心 ±1.25 mm"（宽 2.5 mm），而两侧 dx 均为 **2 mm**，
窗口同时套住了 x=101 mm（p=+295 Pa，界面首排）与 x=103 mm（p=+70 Pa，
已进入床内 2 mm），把两者**平均**成 182.6 Pa ⇒ 凭空出现"112 Pa 跳变"。

**修法**：改用**绝对**窗 `[X_IFACE, X_IFACE + DX) = [100, 102) mm`，恰只含一个
平面（格心 x=101 mm）。两张网格的首个多孔侧格心**都**在 x=0.101 m，可直接对比。

**验证**（`compare_iface.py`，收敛态）：

| 量 | struct 面 | uns 首排 | 参考 @x=100 |
|---|---|---|---|
| u_x (m/s) | 34.233 | 31.351 | 34.708 |
| p (Pa gauge) | 295.0 | 295.0 | 48.5 |

**界面压力跳变 = −0.08 Pa**（判据 (a) 通过）。速度两条结论：

- **uns 首排 −9.67% 不是耦合缺陷，是求解器入口邻格的固有特征**：同一相对亏损
  在三个独立算例里重复——参考 x=1 mm 的 31.832 / 入口 34.7224 = **0.91676**、
  耦合 31.351 / 界面 34.233 = **0.91581**、参考改 `outflow` 后仍为 31.8319；
  且 2–4 层内即恢复（参考 x=3/5/9 mm：34.80/35.23/34.67 m/s）⇒ 只影响首排取值。
- **耦合侧床内梯度自洽（由界面速度直接预测，无需拟合）**：
  `dP/dx = μu/(εK) + ρβu²` 代入 struct 界面面速度 34.233 m/s 得 **2959.2 Pa/m**，
  与耦合 uns 床内实测（拟合 130–190 mm）**−2959.2 Pa/m** 逐位吻合
  （110–190 mm 窗为 −2964.5）。相对解析值的 −2.1% 完全来自界面速度本身的
  −1.37%（0.9863 ⇒ 线性项 ×0.9863、平方项 ×0.9728），即 struct 侧近壁/交换层的
  既有残余（plan.md 流 B 节：struct 近壁 u=13.1 vs uns 9.05 @1.25 mm），
  **不是**界面通量缺陷。

## 5. 已知限制：单求解器参考的**绝对压力水平**含 ≈−250 Pa 内部偏置

`uns_full` 的**梯度**精确（床 −3021.5 Pa/m，相对解析 0.01%），但**绝对水平**整体
偏低 ≈250 Pa。为定位它做了四组对照（同一 4000 格网格，只改 `cell_zone` 与出口 BC）：

| 算例 | zone 8 | 出口 BC | 入口邻格 (x=1 mm) | 内部平台 | 内部最低点 |
|---|---|---|---|---|---|
| 参考（标准对拍） | porous | pressure-outlet | +301.88 Pa（物理解 302 ✓） | 物理解 −251 | −248.93 @x=197（物理解 2.9 −251 ✓） |
| 全流体对照 | fluid | pressure-outlet | −0.05 Pa（物理解 0 ✓） | **−251.17**（整段严格平） | −285.76 @x=5 mm |
| 全流体 ＋ outflow | fluid | outflow | +3.02 Pa | **−248.1** | −282.69 @x=5 mm |
| 床 ＋ outflow | porous | outflow | （被床压降掩盖） | ≈−244 | −546.71 ≈ (−244) − 302 ✓ |

物理解基线：流体半段 dP/dx=0 ⇒ p=0；床半段 dP/dx=−3021.7 Pa/m ⇒ 界面 302 Pa、
出口 0。四组只差两行配置（`cell_zone = 8 fluid` ↔ `perm=1e-6 inertial=1.0
porosity=0.4`；`bc = 6 pressure-outlet` ↔ `outflow`），其余完全相同；取证工具
`python3 plane_profile.py <vtu> <x0> <x1>`（逐平面 p/u 均值 ＋ 拟合 dP/dx）。

四条要点：

1. **与多孔/床无关**：全流体对照内部本应恒为 0，实测**整段严格平**在
   −251.17 Pa（同一 x 平面 40 格的散布 ~1e−3 Pa）⇒ 偏置纯属求解器侧。
2. **与出口 BC 类型无关**：换成 `outflow`（纯 Neumann PPE ＋ 全局质量缩放）后
   内部平台 −248.1，与 pressure-outlet 的 −251.17 只差 3 Pa（1.2%）⇒ 既不是
   `cases/channel/README.md` §3.3 记录的 pressure-outlet 出口段伪调整（那条记的是
   **出口末列**的局部抬高/横向汇聚，与这里的**内部整段水平下移**是两个现象），
   也不能靠换出口 BC 消掉。附带：`outflow` 会把出口末列速度从 34.72 抬到
   36.79 m/s（+5.9%，与 §3.3 的"末列 +3.5% 体速度外观偏差（1/β 不动点）"同族），
   而内部场仍精确 34.72 ⇒ 质量守恒无损，但**不要**拿末列数据、也不要把 `outflow`
   当作该偏置的解药。
3. **严格 ∝ u²**：入口后凹陷深度 −285.76 Pa，在 u/2 → −71.44（×0.250）、
   2u → −1143.06（×4.000）都精确按 u² 缩放 ⇒ 与 ρu²=1419.5 Pa 同阶的**离散
   动量／压力修正相容性项**，不是舍入或迭代历史残留（也解释了逐位可复现）。
4. **空间结构**：内部整体常数下移 ＋ **紧邻入口边界的那一层保持物理解水平**
   （两组独立验证：+301.88≈302、−0.05≈0）＋ 出口末列向上回收。故**梯度不受
   影响**，但凡**绝对压力**对比都会被污染。

**对 C2 判据的直接影响**：参考在 x=101 mm（床内第一层）读到 48.5 Pa
＝ 物理解 302 − 253；而**耦合侧**在 x=101 mm 恰是它自己的**入口邻格**（界面即
入口），读到 295 Pa ≈ 物理解。故 `compare_iface.py` 里"参考 48.5 vs 耦合 295"
那 ~246 Pa 差**不是耦合误差**，而是两侧"域边界位置不同"（一个在域内部、一个是
入口邻格）叠加该偏置所致。

**结论（C2 判据）**：(a) 界面**连续性**——压力跳变 −0.08 Pa，struct 面速度
−1.37%（见 §4）；(b) 床**梯度**——−3021.5 vs 解析 −3021.7（0.01%）。
**不**与参考做绝对压力对齐。耦合侧自洽：界面 295 − 床压降 302 ≈ 出口 0 ✓。

（同理，`check_bed_gradient.py` 的拟合窗取 `0.005–0.095` 与 `0.105–0.195` m，
即剔除首/末几层——那里是 mass-flow-inlet 渐启残留、入口邻格异常与出口回收区；
中间区段梯度精确：`check_bed_gradient.py` 自带窗 `0.105–0.195 m` 给 −3021.3 Pa/m，
另取 `0.110/0.130/0.150/0.170–0.190 m` 四组窗都给 −3021.5 Pa/m——差异只来自窗
端点，两者都 ≈0.01% 于解析 −3021.7 Pa/m。`plane_profile.py` 用同一取法，故其
默认窗的输出与 `check_bed_gradient.py` 逐位一致。）

**开放项**（建议单独立项，本任务不修）：该偏置的机理。候选：(i) 给定通量入口
（`mass-flow-inlet`）在 collocated ＋ Rhie-Chow 下与投影步的相容性源项；
(ii) 压力水平锚点——`cases/channel` §3.3 记录 `outflow` 的 PPE 为纯 Neumann 并
pin cell 1，而本网格 cell 1 正是入口邻格（与"入口邻格保持物理解水平"吻合），
但 pressure-outlet 下（PPE 有 Dirichlet）偏置仍在 ⇒ 锚点说不足以解释全部。
影响面：**所有** uns 单求解器算例的**绝对压力**（既有验证都只看梯度/型线，故一
直未暴露）；凡要与其它求解器或边界值对绝对压力，需先在 `src/unstructured/` 的
投影步与 BC 层定位。

## 6. 复现

```bash
# 参考（单求解器，全流域）
cd uns_full
python3 gen_full.py                                    # 重新生成 unMesh.cas
../../bin/uns_solver unMesh.cas unMesh.control unMesh.vtu > full.log 2>&1
python3 check_bed_gradient.py unMesh.vtu               # 期望 PASS（3.0213 vs 3.0217）

# 耦合算例（在算例目录里跑；产出界面比较所需的两份文件）
#   unMesh_coupled.vtu ← uns 侧终场快照（src/main.f90:475）
#   flow3d.dat         ← struct 侧周期存档；**须 save_interval ≤ n_couple**
#                        （仓库里的 mix.control 是 n_couple=1000、
#                        save_interval=1000，即跑满 1000 迭代才落盘一次）
#   本轮取证：n_couple=400、save_interval=100，mpirun -np 2 约 4 min
cd ..
mpirun -np 2 ../../bin/mixsolver_mpi mix.control Mesh3d.x control.ec \
        unMesh.cas unMesh.control > mix.log 2>&1

# 界面比较（读 unMesh_coupled.vtu + flow3d.dat + uns_full/unMesh.vtu）
python3 compare_iface.py
```

产物落在**运行目录**；本轮取证是在 `/tmp` 的拷贝里跑的（保持算例目录干净），
故仓库里没有这两个文件。400 个耦合迭代已足够：界面量从 ~iter 100 起就不再变
（iter 300/325/350/375/400 的平均 u_x、p 逐位相同），界面压力跳变 −0.08 Pa、
uns 首排 −9.67% 等都已在 §4 定值。

预期关键行：

```
  bed   half gradient :   -3021.3 Pa/m  (-3.0213 Pa/mm, expect 3.0217)
  bed drop (analytic gradient over L=0.1 m) :   302.13 Pa
  RESULT: PASS (within 1% of analytic)
```

`compare_iface.py` 预期关键行（400 迭代、绝对窗 `[100,102) mm`）：

```
  struct face: max +1.38%  RMS 1.37%  mean-ux -1.37%
  uns  cell0 : max +9.68%  RMS 9.67%  mean-ux -9.67%
  pressure: struct 294.9 Pa, uns 295.0 Pa, ref 48.5 Pa (gauge)
  interface pressure jump (struct face - uns cell0): -0.08 Pa
  analytic bed drop: 302.1 Pa
```

### 6.1 §5 绝对压力偏置的四组对照（复现）

```bash
cd uns_full
BIN=../../bin/uns_solver
# 组 1（参考：床 + pressure-outlet）：unMesh.control 原样
$BIN unMesh.cas unMesh.control   ref.vtu  > ref.log  2>&1
# 组 2（全流体 + pressure-outlet）：只把 zone 8 改 fluid
sed 's/^cell_zone = 8 .*/cell_zone = 8 fluid/' unMesh.control > af.control
$BIN unMesh.cas af.control       af.vtu   > af.log   2>&1
# 组 3（全流体 + outflow）
sed -e 's/^cell_zone = 8 .*/cell_zone = 8 fluid/' \
    -e 's/^bc = 6 .*/bc = 6 outflow/' unMesh.control > afo.control
$BIN unMesh.cas afo.control      afo.vtu  > afo.log  2>&1
# 组 4（床 + outflow）：只改出口
sed 's/^bc = 6 .*/bc = 6 outflow/' unMesh.control > bo.control
$BIN unMesh.cas bo.control       bo.vtu   > bo.log   2>&1

# 取证（逐 x 平面 p/u 均值 + 拟合 dP/dx）
python3 plane_profile.py af.vtu  0.110 0.190   # 组 2 → 整段严格平 -251.17 Pa
python3 plane_profile.py afo.vtu 0.110 0.190   # 组 3 → -248.1 Pa（仅差 +3 Pa）
python3 plane_profile.py bo.vtu  0.110 0.190   # 组 4 → 最低点 ≈ -546.7 Pa
python3 plane_profile.py ref.vtu 0.105 0.195   # 组 1 → 床梯度 -3021.3 Pa/m
```

对照逻辑：组 1/2 差"床 vs 全流体" ⇒ 偏置与多孔无关；组 2/3 差"出口 BC 类型"
⇒ 与出口 BC 无关；组 2 把 `bc = 7 mass-flow-inlet 40.868` 改成 `20.434`（u/2）
或 `81.736`（2u）⇒ 凹陷深度 −71.44 / −1143.06 Pa，严格 ∝ u²。

已知无害副作用：uns 单求解器每次运行都会在当前目录落一个 `flux_debug.txt`
（`src/unstructured/main_uns.f90:297` 无条件写出），可删。
