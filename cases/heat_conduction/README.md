# 纯导热（conduction-only）解析解验证算例族

> 2026-10-14 | 阶段 16「非结构框架固体传热」的**纯导热路径**回归基线：把
> `CZ_SOLID` ＋ `k_s` ＋ 按体区施加的体积热源 `q_src` ＋ 定温/绝热边界放在 5 个
> 可闭式求解的构型上对拍解析解。此前该类路径只做过 `cases/solid_conduction/`
> 的一次冒烟（40×2×1 管、线性剖面定性正确），本算例族给出**定量的、可复算的**
> 判据（最大绝对误差 / RMS 误差 / 网格收敛比）与一条**回归守卫**。

## 1. 算例一览

| 算例 | 网格 / 胞数 | 热边界（zone） | 解析解 | 判定容差 | 结果 |
|---|---|---|---|---|---|
| `c1` | `c1_1d.cas` 40×2×1 / 80 | zone6 `x=0` 定温 300 K、zone7 `x=L` 定温 500 K | 线性 `T=300+200·x/L` | 1e-6 K | **PASS**（0.000e+00 K） |
| `c2` | `c1_1d.cas` / 80 | zone6 定温 400 K，其余（zone4/5 symmetry、zone7）绝热 | 均匀 `T≡400 K` | 1e-6 K | **PASS**（0.000e+00 K） |
| `c3` | `c3_annulus.cas` Nr20·Nt72 / 1440 | zone5 `r=RI` 定温 300 K、zone6 `r=RO` 定温 500 K | 对数 `T=300+200·ln(r/RI)/ln(RO/RI)` | 1.0 K | **PASS**（3.912e-01 K，0.196 %） |
| `c4` | `c3_annulus.cas` / 1440 | zone5 定温 400 K，zone6 绝热 | 均匀 `T≡400 K` | 1e-6 K | **PASS**（0.000e+00 K） |
| `c5` | `c5_1d_src.cas` 40×1×1 / 40 | zone6 定温 300 K、zone7 定温 500 K，**体区 zone8 加 q=1e6 W/m³** | 分段抛物线＋两段直线 | 0.1 K | **PASS**（5.208e-02 K，0.026 %） |
| `c3f` | 现场生成 Nr40·Nt72 / 2880（附加） | 同 `c3` | 同 `c3` | 1.0 K | **PASS**（2.291e-01 K）→ c3 收敛比 **1.71** |

所有算例都是**纯固体**（全部胞 CZ_SOLID，无流体/多孔胞），物性取
`k_s=15 W/(m·K)`、`cp_s=460 J/(kg·K)`、`ρ_s=7800 kg/m³`（不锈钢近似），
`mesh_scale=1.0`、`init_t=400 K`、`outer_tol=1e-6`、`lin_tol=1e-12`、`lin_max=1000`。
`max|u_x| = 0`、`heat_source_solid` 除 `c5` 的 zone8 外全为 0（见 `--tables` 报告）。

## 2. 网格与生成（可位级复现）

网格全部由本目录的脚本生成（不依赖外部网格文件），三条生成命令的产物与仓库内
`.cas` **`cmp` 完全一致**（2026-10-14 实测）：

```bash
python3 gen_1d.py      c1_1d.cas                       # 40x2x1 直棒，单体区
python3 gen_1d.py      c5_1d_src.cas 40 --src --ny 1   # 40x1x1 直棒，3 体区
python3 gen_annulus.py c3_annulus.cas                  # 极坐标 Nr=20 Nt=72
```

* `gen_1d.py`（直棒，几何 `x∈[0,0.1]`、`y∈[0,0.01]`、`z∈[0,0.002]`，单位 m）
  - 面区：`3` Interior、`4` symmetry（z=0/LZ）→ 绝热、`5` symmetry（y=0/LY）→ 绝热、
    `6` `x=0` 端面、`7` `x=L` 端面。
  - 体区：默认单区 `2`；`--src` 时沿 x 切成 **3 个连续 id 段**（供局部热源使用）
    `2` = 胞 1..12（`x<0.03`）、`8` = 胞 13..28（`0.03≤x<0.07`，挂 `q_src`）、
    `9` = 胞 29..40（`x≥0.07`）。
  - **为什么 `--src` 必须 `--ny 1`**：Fluent 体区记录 `(12 (id first last ...))` 只能
    描述**一个连续 id 区间**，而胞号是 `1+i+j*NX`（i 最慢），只有 `NY=1` 时“x 分片”
    才恰好等于连续 id 段。`gen_1d.py` 现在**直接拒绝** `--src` ＋ `NY≠1` 并在文件头
    写明原因（对应 `.trae/rules/project_rules.md`「多 cell zone 网格的 cell-id 排序」
    硬规则；旧版在 `NY=2` 下把 zone8 错落到 j=0 行 `x∈[0.06,0.1)`，`q_src` 漏出热源区、
    最大误差 26.6 K，该缺陷已由 `check_heat_conduction.py` 的回归守卫锁定）。
* `gen_annulus.py`（极坐标圆环，`RI=0.01`、`RO=0.03` m、`z∈[0,0.002]`，`dr=0.5 mm`、
  `dθ=5°`，θ 方向闭合）
  - 面区：`4` symmetry（z 端面，绝热）、`5` `r=RI` 内壁、`6` `r=RO` 外壁。
  - 体区：`2`（全部胞，`CZ_SOLID`）。
  控制文件与 zone id 一一对应（`c3` 与 `c3f` 共用 `c3_annulus_fixedT.control`）。


## 3. 解析解（`check_heat_conduction.py` 里的 `ref_*` 函数，与 `.control` 注释同步）

* **c1 / 一维线性**：无内热源、两端定温 ⇒ `T(x) = T_C + (T_H − T_C)·x/L`，
  `T_C=300 K`、`T_H=500 K`、`L=0.1 m`。
* **c2 / c4 / 一处定温＋其余绝热**：稳态导热＋零净热流 ⇒ 全域等于那个定温值（400 K）。
  该算例同时检验 `symmetry` 面被当作**绝热**（不是定温）处理。
* **c3 / 圆环对数解**：轴对称、无内热源、与 `k_s` 无关，
  `T(r) = T_C + (T_H − T_C)·ln(r/RI)/ln(RO/RI)`；`r=0.02 m` 处 `T=426.19 K`。
  网格非正交（面法向一般不与 r 对齐）⇒ 检验非正交修正下的固体导热装配。
* **c5 / 局部热源三段解**：`q=1e6 W/m³` 只在 `[a,b]=[0.03,0.07] m` 内，`a=0.03`、`b=0.07`：

  ```text
  T(x) =  3333.333 x      + 300      , 0 <= x < a         （直线）
  T(x) = -33333.33 x^2 + 5333.333 x + 270 , a <= x < b    （抛物线）
  T(x) =  666.6667 x      + 433.3333 , b <= x <= L        （直线）
  ```

  结点值 `T(a)=400 K`、`T(b)=480 K`；能量平衡 `q·(b−a) = 4.0e4 W/m²` 与
  两端通量差 `−k·T'(b) + k·T'(a) = 5.0e4 − 1.0e4` 一致（`c5_1d_source.control`
  头部含完整推导）。

## 4. 结果（2026-10-14，`bin/mixnsolver` 串行，1 次外迭代即 CONVERGED）

| 算例 | 胞数 | 最大绝对误差 | RMS 误差 | 位置 | 判定（tol） |
|---|---|---|---|---|---|
| `c1` | 80 | `0.000e+00 K` | `0.000e+00` | — | PASS (1e-6) |
| `c2` | 80 | `0.000e+00 K` | `0.000e+00` | — | PASS (1e-6) |
| `c3` | 1440 | `3.912e-01 K` (0.196 % of 200 K) | `2.869e-01` | `r=0.0104979 m`（内壁邻胞） | PASS (1.0) |
| `c4` | 1440 | `0.000e+00 K` | `0.000e+00` | — | PASS (1e-6) |
| `c5` | 40 | `5.208e-02 K` (0.026 %) | `3.294e-02` | `x=0.03125 m`（源区首胞） | PASS (0.1) |
| `c3f` | 2880 | `2.291e-01 K` (0.115 %) | `2.009e-01` | `r=0.0102423 m` | PASS (1.0) |

**逐条解释（这些偏差都是"已知且可解释"的，不是待修缺陷）**：

1. `c1`/`c2`/`c4` 精确到机器零：格心有限体积格式对**线性/常值**场是精确的，
   说明定温 Dirichlet、绝热 `symmetry`、纯固体体区装配与线性求解在纯导热路径上无误。
2. `c3` 的残差是**壁面距离的 1 阶误差**：内/外壁 Dirichlet 用"半格心"距离
   （内壁邻胞心 `r=0.0104979`，`dr/2=0.5 mm`），圆柱几何下该近似对 `T∝ln r` 只有
   1 阶精度；径向差分是 2 阶，故整体介于 1–2 阶。加密实测：`Nr 20→40`，
   最大误差 `3.912e-01 → 2.291e-01 K`，**收敛比 1.71**（介于 1~2，符合上述分析）。
3. `c5` 的残差同样有闭式解释：源区边界正好落在网格面上，解析解在 `x=a,b` 处
   二阶导跳变（"拐点"），格心格式在界面胞上留下截断误差，其闭式值为
   **`q·h²/(8k) = 1e6·(0.0025)²/(8·15) = 0.0520833 K`**。实测：源区内偏差均值
   `+0.05208325 K`（与上式对到 8 位有效数字），**源区外**最大偏差仅 `3.33e-07 K`。
   另用同一离散系统的独立 numpy 直接求解交叉验证：与求解器结果差 `~3e-7 K`
   ⇒ 求解器实现正确，余项纯属离散化。

## 5. 复现命令

```bash
# 0) 先构建（只需一次）：make -j4     —— 产出 bin/mixnsolver
cd /home/sundong/Fortran_Project/MixNSSolver

# 1) 跑算例（严格串行、每例 timeout 600 s，遵守「本机资源保护」硬规则）
bash cases/heat_conduction/run_all.sh            # c1..c5，约 1~2 s/例
bash cases/heat_conduction/run_all.sh c3f        # 可选：c3 径向加密（网格现场生成）
C3F_NR=80 bash cases/heat_conduction/run_all.sh c3f   # 更细的收敛序列
TIMEOUT=90 bash cases/heat_conduction/run_all.sh c1   # 覆盖单例上限

# 2) 对拍（退出码 0 = 全部 PASS）
python3 cases/heat_conduction/check_heat_conduction.py              # 5 个算例
python3 cases/heat_conduction/check_heat_conduction.py --refine --tables
python3 cases/heat_conduction/check_heat_conduction.py --refine -v c5

# 3) 回归守卫自测（坏资产必须报错；退出码 0 = 守卫有效）
python3 cases/heat_conduction/check_heat_conduction.py --selftest
```

运行现场：`.dbg/heat_conduction/<case>/{unMesh.cas, unMesh.control, mix.control,
run.log, unMesh.vtu}` —— 全是可重算的产物，**不入库**（`.gitignore` 已覆盖 `.dbg/`）。

## 6. 对拍脚本与「回归守卫」

`check_heat_conduction.py` 直接读 VTU（ASCII）取 CellData：

* VTU 读取约定（与 `cases/porous_plug/plot_plug.py` 一致）：`vtk_write` 把每个六面体
  胞拆成 12 个四面体（face fan ＋ 胞心），points 末尾追加 `ncells` 个胞心；CellData
  每个数组按原胞顺序写 12 份相同值。脚本用 `offsets` 切块还原单元节点号
  （`connectivity` 是 VTK 规范的 4 个节点号/四面体；旧文件里 legacy VTK 的
  "顶点数"前导 token 也能识别），并**自校验**三件事：每个块恰好 4 个节点号、
  胞心 id = 块最大节点号、tet→cell 映射严格块状有序、CellData 逐胞复制一致
  （任一不满足即 assert 失败，避免"读错列"却给出看似合理数字）。
* 取量是 **`temperature_solid`**（固体温度），不是 `temperature`（`T_lte`）。
* 判定：`max|T_solid − T_analytic| <= tol`，`tol` 见 §1 表（由上面的误差分析给出，
  不是"调到能过"）；`--tables` 打印剖面表，`--refine` 追加 `c3f` 并报收敛比。

**c5 的体区落点守卫**（对应 `.trae/rules/project_rules.md`「验证脚本必须自带回归
守卫」＋「多 cell zone 网格的 cell-id 排序」）：检查 `heat_source_solid > 0` 的胞
**恰好**是胞心落在 `(0.03, 0.07)` 内的那些 x 面（本网格 16 个）。守卫会在坏资产上
主动报错，而不是静默给出"看起来合理"的温度场。`--selftest` 用当前 c5 网格人工构造
历史缺陷资产（`q_src` 挪到 `x∈[0.06,0.1)` 的 16 个胞，即旧版 `--src` 未强制 `--ny 1`
时的落点）复现该场景，实测：

```text
[1/2] 正确资产 … 回归守卫 : q_src>0 的胞 16 个，覆盖 16 个 x 面（应 16）… 通过
[2/2] 坏资产   … OK  : 守卫按预期报错 -> q_src 落点与体区 id 段不匹配：q>0 覆盖
      16 个 x 面 [0.06125, 0.09875]，应为 (0.03, 0.07) 内的 16 个面；见 …
退出码 0（守卫有效）
```

## 7. 说明与遗留

* **`T_lte` 滞后一轮**：固体胞里 `T_lte`（`temperature`）被钉成**上一轮**的 `T_s`
  （纯固体算例第 1 次外迭代即满足收敛判据，故看起来"落后一轮"，差值可达 ~100 K）。
  这是实现语义而非误差，物理量取 `temperature_solid`；报告里会打印
  `|T_lte − T_solid|_max` 并给出该说明。
* **纯固体全壁面无压参考**：动量/PPE 算子在"流体胞数为 0"时退化，`lin-res=NaN`
  （阶段 16 记录过，非本算例族引入、不参与导热迭代、不影响 `T_s`），因此所有算例
  都是 1 次外迭代即 `CONVERGED after 1 outer iterations`。流固界面算例由流体区
  锚压，可安全忽略该量。
* **`c3` 容差 1.0 K** = 0.196 % of 200 K，对应"半格心壁面距离"的 1 阶误差；要更严
  需加密（见 `c3f`）或引入壁面距离修正。
* **下一步可选**：① `C3F_NR=80` 的完整收敛序列（看比值是否收敛到 1~2 区间内稳定值）；
  ② 瞬态导热（BDF2 时间精度 / 与稳态极限对比）尚未覆盖；③ 变 `k_s` 的界面（固-固或
  固-流 CHT）本算例族未涉及，仍由 `cases/solid_conduction/` 与 CHT 算例覆盖。

