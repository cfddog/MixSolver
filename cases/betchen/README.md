# Betchen 2006 多孔-流体界面算例（BJ / PLUG / HT）

复现 Betchen, Straatman & Khandelwal (2006), *Int. J. Heat Mass Transfer* 49: 4596-4610
的多孔-流体界面算例，用论文给出的解析/参考解验证非结构求解器的:

- **PLUG（算例 2，垂直界面，Fig.7/8）**：管道中部一段多孔塞，界面处压力梯度不连续
  （kink），是"界面条件"最敏感的算例；
- **BJ（算例 2，平行界面，Fig.5）**：下多孔 / 上流体，验证 Beavers-Joseph 界面剪切；
- **HT（算例 3，Fig.9）**：局部热非平衡换热（`ht_a` / `ht_b`，另见 `plot_ht.py`）。

参考数据（用户提供的数字化结果）读自 `/mnt/c/temp/validate_case/`：
`BJ_1.csv`/`BJ_2.csv`（Da=1e-2/1e-3 的 u/U-y 剖面），
`plug_1_u.csv`/`plug_1_p.csv`（Da=1e-2）、`plug_2_u.csv`/`plug_2_p.csv`（Da=1e-3）
（中心线 u/U 与 p/ρU² 沿 x），另见 `validation_targets_summary.csv`、`数据说明.md`。

## 1. 算例与文件

| 算例 | 网格/脚本 | 控制文件 | 关键设置 | 状态 |
|---|---|---|---|---|
| PLUG Da=1e-2 (Fig.7a) | `gen_plug.py` → `plug_dae2.cas` | `plug_dae2.control` | Re_H=1, K=1e-6, ε=0.7, 3H\|2H\|3H | ✅ 已对拍 |
| PLUG Da=1e-3 (Fig.7b) | `plug_dae3.cas` | `plug_dae3.control` | Re_H=1, K=1e-7 | ✅ 已对拍 |
| PLUG 高 Re (Fig.8) | `plug_hir.cas` | `plug_hir.control` | Re_H=1000, K=1e-6, inertial=244, 5H\|5H\|50H | 定性（无数字化参考） |
| BJ Da=1e-2 (Fig.5) | `gen_bj.py` → `bj_dae2.cas` | `bj_dae2.control` | Re_H=1, K=1e-6, 8H×2H, 100×40 | ✅ 形状对拍 |
| BJ Da=1e-3 (Fig.5) | `bj_dae3.cas` | `bj_dae3.control` | K=1e-7 | ✅ 形状对拍 |
| HT (Fig.9) | `gen_ht.py` → `ht_a/b.cas` | `ht_a/b.control` | ε=0.9118, q=315 W/m | 另见 `plot_ht.py` |

`gen_plug.py` 的第 3 个参数支持 `n1,n2,n3`（每段 x 方向网格数，缺省 20/20/20），
用于界面附近加密网格的网格收敛性研究。

## 2. 数值设置与变量约定

**体积平均（Brinkman + Darcy-Forchheimer）动量方程**，采用 Vafai-Kim 的"按孔隙率
除过"形式（速度取表观/seepage 速度 u）：

```
ρ/ε·∂u/∂t + ρ/ε²·∇·(uu) = -∇p + μ/ε·∇²u - μ/K·u - ρ·cE/√K·|u|u
```

映射到 `.control`：

| 方程项 | 代码 | 说明 |
|---|---|---|
| 粘性项 μ/ε ∇²u | 面扩散系数 D = (μ/ε)·A/dw | 两格取线性平均（`compute_gradients` 路径） |
| Darcy 汇 μ/K·u | `ap += μ/K·V`（**不含 1/ε**） | `cell_zone ... perm=K` |
| Forchheimer 汇 ρ·cE/√K·\|u\|u | `rhs -= ρ·inertial·\|u\|·u·V`，`inertial = cE/√K` | `cell_zone ... inertial=...` |
| 孔隙率 | `porosity=ε` | 只进粘性项与对流项（`ρ/ε`, `ρ/ε²`） |

无量纲化（与参考数据一致）：`U0 = μ/(ρH)`（Re_H=1）、`u/U0`、`p/(ρU0²)`；
高 Re 算例 `U0 = 1000μ/(ρH) = 1.5684 m/s`。多孔区若为 **bevel-free 均匀床**，则
`cE = 1.75ε/(150ε³)^0.5`（ε=0.7 → 1.4166；代码里用 `inertial = cE/√K`）。

入口：论文用充分发展的抛物剖面 `u(y)=6Uy/H(1-y/H)`（峰值 1.5U0、均值 U0）。
求解器自 2026-10-08 起支持

```
bc = <zone> velocity-inlet-parabolic <Umean> <span> <axis> <origin> [<T>]
```

即在面法向施加 `|u| = 6·Umean·ξ(1−ξ)`，`ξ = (x_axis − origin)/span`（详见 §4.4）。
PLUG 三个算例已改用它，例如 `bc = 9 velocity-inlet-parabolic 1.5684e-3 0.01 2 0.0 300.0`
（H=0.01 m、y 向、下壁 y=0）。3H（高 Re 5H）上游段仍保留：Re_H=1 时剖面在
1–2 个单元内就已与解析抛物一致，上游段的作用是把界面推离出口边界的影响。

**界面处理**：CAS 中流体/多孔分块为内部面，默认走"速度 + 应力连续"；压力面值
在流体/多孔界面处使用 §4.1 的一致化重构。

## 3. 验证结果（2026-10-08 界面压力一致化修复后）

A/B 对拍：`abtest_interface_pressure.sh` 把修复 stash 掉编出 pre-fix 二进制，两版跑
同一批算例（每格的 `.cas/.control` 复制到临时目录，串行跑），再逐列/全场对比。
抖动度取界面邻域（x/H∈[2.3,5.8]）相邻列速度的二阶差分幅值
`alt_amp = max|u_i-(u_{i-1}+u_{i+1})/2|`，速度以 U0=μ/(ρH) 无量纲化
（高 Re 算例 U0=1000μ/(ρH)）；参考解为 `/mnt/c/temp/validate_case/plug_{1,2}_{u,p}.csv`
的数字化曲线（`vs_ref.py` 取最近截面插值比较）。

### 3.1 PLUG（垂直界面，Fig.7a/7b）

| 工况 | alt_amp 前 → 后 | 界面邻域 u/U0 前 → 后 | u L2(x≥2H) 前 → 后 |
|---|---|---|---|
| Da=1e-2 | 0.147 → **0.018** | 1.225–1.497 → 1.266–1.494 | 2.94% → **0.58%** |
| Da=1e-3 | 0.781 → **0.074** | 0.853–1.735（±40%）→ 1.069–1.491 | 13.72% → **2.34%** |
| Re_H=1000 | 0.142 → **0.022** | 1.024–1.238 → 1.032–1.195 | —（图 8，定性） |

修复前是"锯齿"解：Da=1e-3 中心线在 x/H=2.9/3.0/4.95/5.08 处为
`… 1.437 1.518 0.858 1.724 …`（±35% 逐列振荡），修复后同一段单调
`… 1.272 1.126 1.080 1.081 …`。Re_H=1000 有符号二阶差分从
`+0.0142 −0.0472 +0.1417 −0.1082`（正负交替）变为 `+0.0058 +0.0071 +0.0215 −0.0196`
（单调过渡层）。

压力：plug 段（x/H∈[3.25,4.75]）线性拟合 dp/dx 与入口压力（p/ρU0²）：

| 工况 | 量 | 前 → 后 | 数字化参考 | 相对误差 |
|---|---|---|---|---|
| Da=1e-2 | dp/dx (/H) | 177.5 → **130.7** | 140.3 | −6.8% |
| | p_inlet/(ρU0²) | 423.6 → **330.9** | 345.8 | −4.3% |
| | p L2 | 15.81% → **2.57%** | — | — |
| Da=1e-3 | dp/dx (/H) | 1506.1 → **1069.3** | 1091.6 | −2.0% |
| | p_inlet/(ρU0²) | 3035.4 → **2216.9** | 2264.4 | −2.1% |
| | p L2 | 21.81% → **1.83%** | — | — |

界面压力 jump（p_inlet 与 plug 段斜率）在修复后与数字化参考一致到 2–7%，剩余的
偏差主要是入口用均匀速度而非抛物线（见 §3.4 与 §5）。

### 3.2 BJ（平行界面，Fig.9a/9b）

| 工况 | 与参考解 L2(u) 前 → 后 | 收敛 |
|---|---|---|
| Da=1e-2 | 4.09% → **3.54%** | mass-imbal 1.2e-19, du_max 8.7e-19（机器零） |
| Da=1e-3 | 2.88% → **2.87%** | 同上 |

BJ 的界面与流向平行，Darcy 汇造成的斜率折转只在横向起作用，故收益远小于 PLUG。
参考曲线用约 2×U0 的归一化速度，脚本按常数比例对齐峰后只比形状（对齐后
my-fluid-mean-peak 1.394/1.455 vs ref-peak 1.387/1.457，量级 <1%）。

### 3.3 位级无副作用（no-op）与 Darcy 汇约定

| 算例 | 检查 | 前 → 后 |
|---|---|---|
| fluid（PLUG 网格，无多孔块） | VTU + 求解日志 md5；全字段差分 | **完全相同**（md5 一致，max\|du\|=max\|dp\|=0） |
| ppd（1-D Darcy 塞） | plug 段 dp/dx 拟合 | −461.5 Pa/m（+150.0%）→ **−184.6 Pa/m（−0.0%）** |

fluid 算例的两版二进制给出逐位相同的速度/压力场与求解日志 ⇒ §4.1/§4.2 的改动
只作用于多孔单元或流体/多孔界面，单相路径完全不受影响。

### 3.4 抛物入口剖面（第二轮：UNI → PAR，同一二进制）

只改 `.control` 里入口那一行（`velocity-inlet` → §4.4 的
`velocity-inlet-parabolic`，均值仍取 U0，故流量不变）：

| 工况 | u L2（全中心线） | max\|du\| | u L2（x≥2H） | p L2 | p_inlet/(ρU0²) |
|---|---|---|---|---|---|
| Da=1e-2 | 8.06% → **0.63%** | 0.401 → **0.021** | 0.579% → 0.494% | 2.57% → **2.37%** | 330.9 → **333.8**（参考 345.8） |
| Da=1e-3 | 7.02% → **2.23%** | 0.402 → **0.097** | 2.34% → 2.35% | 1.83% → **1.76%** | 2216.9 → **2221.9**（参考 2264.4） |

（`vs_ref.py /tmp/plugparab/{dae2,dae3}UNI|PAR/plug_*.vtu=<u.csv>=<p.csv>=<标签>`；
max\|du\| 修复前都落在入口采样点 x/H≈0.03–0.06，修复后降到 1 个单元以内的
有限体积滞后。）

入流形状（Da=1e-3，z 中心面）：修复前第一列（x/H=0.075）平顶 `u/U0≈1.10`
（壁面格 0.33，列内 τ 积分 0.944·U0H），修复后第一列即抛物（峰值 1.483）、
第 4 列（x/H=0.525）与解析 `6ξ(1−ξ)` 逐点差 ≤0.004（≈0.3% 峰值 1.5U0），
`∫u dy/(U0H)` = 1.001；越往下游越准（第 5 列 0.675 处 max 差 0.0027）：

```
ξ=y/H : 0.02 0.07 0.12 0.17 0.21 0.26 0.31 0.36 0.40 0.45 0.50
解析  : .139 .398 .629 .833 1.010 1.160 1.282 1.378 1.446 1.486 1.500
我们  : .143 .401 .632 .835 1.011 1.161 1.283 1.378 1.446 1.486 1.500
```

- **流量不变**：剖面均值恰为 Umean，plug 段 dp/dx 与 p_inlet 只变 0.2%
  ⇒ 收益完全来自入流**形状**，与「压力整体偏差 ~1.2×」的旧猜测无关。
- **残余**：Da=1e-3 的 max\|du\|=0.097 已移到多孔段界面过渡层
  （x/H≈2.96：参考 1.140 vs 我们 1.237），与入口无关（见 §5 第 2 条）。
- **无回归**：均匀入口路径逐位不变 —— 用新二进制 + 旧入口行重跑 `plug_dae3`，
  `plug_dae3.vtu` 与仓库内（旧二进制产出的）文件 **md5 完全一致**。

## 4. 本次修复（2026-10-08）

### 4.1 界面面压力的一致化重构（kink-consistent face pressure）

**症状**：PLUG（压力梯度在界面处折转）算例在界面两侧出现逐列速度"锯齿"，
Da=1e-3 达 ±40%、Re_H=1000 达 ±10%，且跨界面压力/速度误差随网格细化不收敛
（odd-even 模态），与参考解对比 u L2 13.7%。

**机理**：跨流体/多孔界面压力只有 C0 连续、斜率跳变（多孔侧多出 Darcy 汇
μ/K·u，PLUG 算例局部约 30–60 倍）。距离加权插值

```
pf_lin = lf*p(c0) + (1-lf)*p(c1)
```

在界面面上的误差为 `du = d0*d1*(s_por-s_flu)/(d0+d1)`（d = 单元中心到面的距离），
对 PLUG 是 O(30 ρU0²) 量级，远大于局部粘性压力变化。这个 pf 以等值反号进入两侧
单元的动量方程（压力力 −p_f·S_f），构成**力偶极子**：不改变总动量，却驱动一个
odd-even 速度模态，而该模态唯一的阻尼是界面附近很弱的流向粘性/对流耦合 —— 这正是
观察到的锯齿。它同时污染 `compute_gradients` 的 Green-Gauss 压力梯度与 Rhie-Chow
通量。

**修复**：新增 `kink_face_pressure`（`mod_uns_fields.f90`），对界面面改用两侧
单侧二次重构的平均

```
pf_quad = 0.5*( p(c0) + d0*g(c0).n + p(c1) + d1*g(c1).n )
```

- 对分段线性（kink）场，两侧重构在界面上精确相等 ⇒ pf_quad 精确等界面压力；
- 对局部线性场退化为线性插值 ⇒ 光滑压力无副作用（§4.3 位级验证）；
- 该面值同时被 `compute_gradients`（压力面值循环）与 `momentum_assembly`
  （−p_f·S_f）使用，两者只由单元类型判定 `is_porous_cell(...) .neqv. ...` 触发，
  因此（动量压力力, 单元梯度）的自洽不动点恰为两侧各自的真实斜率。

### 4.2 Darcy 汇的 1/ε 因子（一维多孔塞算例暴露）

`momentum_assembly` 的 Darcy 汇原为 `ap += (μ/ε)/K·V`，多乘了 1/ε。体积平均动量
方程（Vafai-Kim "除孔隙率"形式）

```
ρ/ε·∂u/∂t + ρ/ε²·∇·(uu) = -∇p + μ/ε·∇²u - μ/K·u - ρ cE/√K·|u|u
```

中 1/ε 只出现在时间/对流/粘性项，**Darcy 汇本身是 μ/K·u（表观速度）**。多乘 1/ε
使一维塞压降放大 1/ε：`cases/porous_plug`（K=1e-8, ε=0.4）解析 dp/dx=−184.6 Pa/m，
修复前算得 −461.5 Pa/m（+150%）；Betchen 算例 ε=0.7（1.43×），入口压力相应偏高
30–40%。修复后 porous_plug 四个工况误差
0.000% / 0.001% / 0.317% / 0.000%，详见 `cases/porous_plug/README.md` 第 7 节。

### 4.3 单相算例的 no-op 验证

`fluid` 工况（同一 PLUG 网格与边界，但无多孔块 ⇒ `is_porous_cell` 恒假）在两版
二进制下的 VTU md5 完全相同、求解日志逐字节一致 ⇒ 修复对单相路径无影响。

### 4.4 抛物（充分发展）入口边界 `velocity-inlet-parabolic`

新增 BC 类型 `BC_VINLET_PARAB = 10`，语法

```
bc = <zone> velocity-inlet-parabolic <Umean> <span> <axis> <origin> [<T>]
```

- `bc_face_vel`（`mod_uns_bc.f90`）把面速度取为沿**面内法向**的一维抛物剖面
  `|u_f| = 6·Umean·ξ(1−ξ)`，`ξ = clamp((x_axis − origin)/span, 0, 1)`，
  `u_f = −|u_f|·n_outward`（与 `mass-flow-inlet` 一样指向域内，故 `span/origin`
  可直接写成"壁面到壁面"的几何量）；`clamp` 保证面心略微落在名义区间外时仍
  `|u_f| ≥ 0`。
- `axis`（1=x / 2=y / 3=z）选定剖面方向；2D 直接用，3D 时剖面沿 `axis`、
  其余方向均匀（面法向仅由 `axis` 分量构成 ⇒ 无横向入流）。
- **流量等价**：`6∫₀¹ξ(1−ξ)dξ = 1`，剖面连续平均恰为 `Umean` ⇒ 与
  `velocity-inlet Umean 0 0 0` 的总流量相同（离散面心采样误差 O(dx²)，PLUG 网格
  实测相对 ~1e-3 量级）。
- 接入面：动量装配（定速度 Dirichlet，与 `velocity-inlet` 同一个 case 列表）、
  温度（`tval` 作入口 Dirichlet，与 `velocity-inlet` 同）、LTNE 固相入口识别、
  BC 报告（打印 `Umean`/`span`/`axis`/`origin`）。`parse_bc` 校验
  `span > 0`、`axis ∈ {1,2,3}`，缺 token 或非法参数即报错退出。

**验证**：见 §3.4（PLUG Da=1e-2/1e-3 的 u L2 8.06%/7.02% → 0.63%/2.23%，
第 4 列入流剖面与解析 `6ξ(1−ξ)` 逐点差 ≤0.004 = 0.3% 峰值、`∫u dy/(U0H)`=1.001）；
均匀入口路径逐位不变（§3.4 末条）；MPI np2 与串行该网格 max|du|=2.8e-9。

## 5. 已知局限与后续

1. ~~**入口剖面**~~（**已解决**，2026-10-08 第二轮）：`velocity-inlet-parabolic`
   （§4.4）直接施加论文的抛物剖面，Da=1e-2 的 u L2 由 8.06% 降到 0.63%（§3.4）。
   ⇐ **残余 2–7% 压力偏差不再由入口解释**：抛物剖面均值同为 U0，压降/入口压力只变
   0.2%，故偏差更可能来自 (2) 的界面过渡层分辨率或多孔区有效黏度/惯性项约定。
   Re_H=1000 仍只能定性：发展长度 ~60H，5H 上游段远不足以在界面处给出充分发展剖面
   （`plug_hir` 也已改用抛物入口，但界面处的剖面形状与论文的充分发展直方图不同）。
2. **界面加密**：Da=1e-3 的界面过渡层只有 2–3 层网格宽。`gen_plug.py` 现支持
   第 3 个参数 `n1,n2,n3`（分段网格数），可做界面加密与网格收敛性研究（待办）。
3. **Beavers-Joseph 滑移**：当前界面条件只有"速度连续 + 应力连续"，未实现 BJ 滑移
   系数；BJ 算例（平行界面）修复收益相对小（L2 4.09%→3.54%），残余差异可能与其
   有关。
4. **BJ 参考解约定**：参考曲线用约 2×U0 的归一化尺度，现按常数比例对齐后只比较
   形状；近界面剖面差异尚未逐项分解。

## 6. 复现命令

```bash
# 网格 + 求解
python3 gen_plug.py plug_dae2.cas          && ../../bin/uns_solver plug_dae2.cas plug_dae2.control > plug_dae2.log 2>&1
python3 gen_plug.py plug_dae3.cas          && ../../bin/uns_solver plug_dae3.cas plug_dae3.control > plug_dae3.log 2>&1
python3 gen_plug.py plug_hir.cas hir       && ../../bin/uns_solver plug_hir.cas plug_hir.control > plug_hir.log 2>&1
python3 gen_plug.py plug_ref.cas 40,40,40  # 分段网格数（界面加密用）
python3 compare_plug.py                    # -> images/plug_compare_ref.png

# 分析脚本（本次新增）
python3 cmp_centerline.py plug_dae3.vtu=post            # 抖动度 alt_amp / 压降 / 入口压力
python3 vs_ref.py plug_dae3.vtu=plug_2_u.csv=plug_2_p.csv=post   # 与数字化参考解 L2
python3 fdiff_vtu.py a.vtu b.vtu label                  # 全字段 VTU 差分（任意网格/宽高比）
../porous_plug/fit_pp.py ../porous_plug/plug_darcy.vtu  # 1-D 压降拟合 vs 解析

# 入口 BC 的 UNI→PAR 对拍（同一二进制，只改 .control 的一行）
python3 inlet_profile.py plug_dae3.vtu=PAR   # 前几列 u(y)/U0 profile vs 解析抛物

# 完整 A/B（pre-fix vs post-fix 二进制；会临时 git stash 源码，仅验证用）
#   要求修复尚未提交（脚本靠 stash 出 pre-fix 版本）
bash abtest_interface_pressure.sh
```

A/B 脚本做的事：把 `mod_uns_fields.f90`/`mod_uns_simple.f90` 的改动 `git stash`
⇒ 编译出 pre-fix 二进制 ⇒ 用同一批 `.cas/.control` 在 `/tmp` 里串行跑 PRE ⇒
`git stash pop` + 重编译（md5 与 A/B 前完全一致，构建可复现）⇒ 跑 POST ⇒
打印 `git stash list`/`git status`。两版结果都用上面的分析脚本对比。

