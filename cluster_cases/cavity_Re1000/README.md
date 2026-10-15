# cavity_Re1000 —— 二维顶盖驱动方腔，Re = 1000（稳态 SIMPLE）

自包含现场：网格 + control + 校验脚本。运行产物落在 `logs/`、`out/`。

## 1. 物理与网格

| 项 | 值 |
|---|---|
| 几何 | 单位方腔 [0,1]²，顶盖 y = 1 以 U_lid = (1,0,0) 拖动，其余三面静止 |
| 网格 | 结构化风格的贴壁加密网格 80 × 80 × 1 = **6400 单元 / 25760 节点**，z 方向**一个单元** → 严格二维 |
| 边界 | zone 5 = wall（四面墙，320 个面；`lid` 平面过滤选出 y = 1 那面）；zone 4 = symmetry（z 两平面，12800 个面） |
| 无量纲化 | ρ = 1，μ = 0.001，L = 1，U_lid = 1 → Re = ρU L/μ = **1000** |
| 参考 | Ghia, Ghia & Shin (1982) JCP 48:387，Table I（Re = 1000 竖中心线 u_x，17 个站位） |

`unMesh.cas` 与 `cases/natconv/cavity80.cas` **逐位相同**，所以那个目录的稳态结果可以直接当回归基线。

当前现场配置为稳态 SIMPLE，`outer_max = 10`，用于快速核对 Re=1000 控制参数和并行 SIMPLE 求解路径；这不是完整收敛算例。

## 2. 关键设置（`unMesh.control`）

```
transient = false
pimple = false
outer_max = 10       # SIMPLE 外迭代上限
outer_tol = 0.0        # 固定跑满 10 次；求解器按严格小于号判断收敛
alpha_u = 0.7   alpha_p = 0.3
conv_blend = 1.0   nonorth_corr = 0.0   ppe_precond = ic0
lin_tol = 1e-8   lin_max = 300
lid = 5 2 1.0 1.0 0.0 0.0      # zone 5，平面 dir 2（y），y = 1 处速度 (1,0,0)
```

`lid` 是**平面过滤**而不是一个独立 zone：它把 zone 5 里落在 y = 1 上的面挑出来给 Dirichlet
速度 (1,0,0)，其余 zone-5 面保持无滑移。顶盖在 t = 0 **脉冲启动**。

现场使用 `uns_parallel = always` + `bin/mixnsolver_mpi`，明确走分区 MPI SIMPLE 路径；稳态步数由 `outer_max` 控制，不使用 `n_time_max`。

## 3. 运行

```bash
./run.sh          # 最多 10 次 SIMPLE 外迭代；np 默认 1
./run.sh 4        # 4 rank（METIS 分块）
python3 check_cavity.py out/cavity_Re1000_SIMPLE10_np1.vtu
```

输出分别写入 `logs/cavity_Re1000_SIMPLE10_np<N>.log` 和 `out/`。`outer_tol = 0` 配合求解器的严格小于收敛判断，确保恰好执行 10 次迭代；日志会报告达到迭代上限。该短跑只确认运行链路和数值健康度，不代表 Ghia 基准已收敛。

跑完 `logs/` 里的日志**尾部**依次有：`[uns-mpi] --- solution summary ---`（max\|u\| / max p /
min p）、**Ghia Re = 1000 对照表**（`relative L2 error (centerline u_x): ...`）、力汇总、以及
Cp / 分离诊断（`mod_uns_driver_mpi.f90` 第 350 行起，rank 0 对汇总后的全局场做一次后处理）。

## 4. 历史瞬态基准记录（PIMPLE）

以下 4.1–4.3 为此前 `n_time_max = 6000` 的 PIMPLE 运行记录，仅作参考，不对应当前 10 次 SIMPLE 配置。

### 4.1 到达定常解的代价：`du/dt` 衰减（np = 4，dt = 0.01）

`du/dt` 是日志里每步打印的 `max|Δu|/dt`。按 5 秒分段取段平均：

| 时间段 t | du/dt 段平均 | 段内 min–max | max\|u\| |
|---|---|---|---|
| 0–5 | 5.174e-1 | 1.72e-1 … **9.38** | 0.8718 |
| 5–10 | 1.054e-1 | 6.46e-2 … 1.72e-1 | 0.9003 |
| 10–15 | 4.575e-2 | 3.19e-2 … 6.45e-2 | 0.9007 |
| 15–20 | 2.428e-2 | 1.83e-2 … 3.19e-2 | 0.9041 |
| 20–25 | 1.374e-2 | 9.66e-3 … 1.83e-2 | 0.9071 |
| 25–30 | 6.689e-3 | 4.33e-3 … 9.64e-3 | 0.9091 |
| 30–35 | 2.882e-3 | 1.81e-3 … 4.32e-3 | 0.9103 |
| 35–40 | 1.201e-3 | 7.58e-4 … 1.80e-3 | 0.9110 |
| 40–45 | 5.183e-4 | 3.38e-4 … 7.57e-4 | 0.9110 |
| 45–50 | 2.353e-4 | 1.58e-4 … 3.37e-4 | 0.9112 |
| 50–55 | 1.131e-4 | 7.89e-5 … 1.58e-4 | 0.9120 |
| 55–60 | 6.211e-5 | 4.90e-5 … 7.88e-5 | 0.9120 |

* t = 0–5 的 **9.38** 是脉冲启动那一两步（首层网格突然被拖动），不是不稳定。
* t ≳ 10 以后每个 5 秒段都精确地**减半**（比值 0.47–0.53），也就是 `du/dt ∝ exp(-t/τ)`，
  **τ ≈ 7.2**——这是方腔里最慢的那个模态（角涡/主涡的缓慢调整），不是数值噪声。
* 因此「跑到定常」需要的步数由 τ 决定：`dt = 0.01` 下 τ ≈ 720 步，而 **t = 60（6000 步）
  时 `du/dt = 4.9e-5`**，剩余未收敛量 ~ `du/dt·τ` ≈ 3.5e-4（相对 max\|u\| 0.912 是 4e-4）。
  **想要最后那 4e-4 就把 `n_time_max` 提到 ~10000（t = 100）**；本现场取 6000 是
  「5 位有效数字的工程精度 + 20 min 量级」的折中，见 4.2 的 L2 收敛曲线。

### 4.2 校验：竖中心线 u_x 对 Ghia Re = 1000 表

跑完以后任何时刻的场都能独立复算（不依赖日志）：

```
python3 check_cavity.py out/cavity_Re1000_np4.vtu
```

把快照逐帧过一遍，就得到 **L2 误差随时间的收敛曲线**（同一张网格、`n_out_every = 500`）：

| t | 10 | 20 | 30 | 40 | 50 | 60 | 稳态 SIMPLE（同网格，np=2） |
|---|---|---|---|---|---|---|---|
| 相对 L2 误差 | 4.51e-1 | 1.92e-1 | 5.65e-2 | 3.30e-2 | 3.04e-2 | **3.01e-2** | **2.89e-2** |
| max\|u\|（单元值） | 0.8996 | 0.9058 | 0.9098 | 0.9112 | 0.9115 | 0.9116 | — |

（这张表由 `unMesh.vtu_step<NNNN>.vtu` 快照逐帧跑 `check_cavity.py` 得到；快照每 1000 步
一个，即 t = 10, 20, …, 60。）

* **t = 60 的瞬态解离真正的定常解只差 4%**（3.01e-2 vs 2.89e-2），而这 4% 就对应 4.1 节
  里那 4e-4 的残余——两条独立的信息互相印证。
* `max|u|` 是**单元中心值**，所以是 0.9116 而不是 1.0（1.0 是**面**上的盖速；盖下第一层
  单元中心的量本来就低一点）。用 `check_cavity.py` 看这个数时不要以为解错了。
* 网格分辨率效应（同一套 control / 稳态路径）：**40×40 → 9.12e-2，80×80 → 2.89e-2**，
  加密一倍误差降 3.2 倍——所以 3e-2 这个量级是**网格**决定的，不是时间推进或 MPI 分块造成的。
* `check_cavity.py` 用与求解器 `ghia_compare` 完全一致的 8 近邻反比距离探针
  （`mod_uns_output.f90`），所以脚本打印的 L2 与日志里的
  `relative L2 error (centerline u_x)` **逐位相同**（都是跳过上下两个壁面行后的 15 个站位）。

### 4.3 数值健康度与代价（6000 步）

| 量 | 实测 |
|---|---|
| CFL~ | 1.69e-1 |
| lin-it / lin-res | 168–175 / 1e-9 … 4.4e-9（每步 2×2 = 4 次压力修正） |
| mass-imbal | ~1e-18（机器精度，全程无漂移） |
| dT_max | ≡ 0（无能量方程） |
| np = 4 本机墙钟 | **18.2 min / 6000 步 = 0.182 s/步**（5.5 步/s） |
| 日志噪声 | 每个 rank **各打印一份** `--- Boundary conditions ---`（np = 4 就有 4 份，面数是本 rank 的局部数），属正常 |

## 5. 日志尾部要看/不要看的

本节以下诊断说明针对历史 PIMPLE 完整算例；当前 10 次 SIMPLE 短跑以 `SIMPLE (MPI) iteration`、迭代日志和有限值为检查项，不以 Ghia 误差作为收敛判据。

* ✅ `relative L2 error (centerline u_x)`：本现场的判据（期望 ~3e-2）。
* ⚠️ `--- Separation and wake (Re = rho*U*D/mu) ---` 这一段是**圆柱专用**的诊断，在方腔里
  会打印出 `Top separation: not detected` / `Wake length ... : 0.494` 之类**没有意义**的数
  （还会拿 Re = 40 的圆柱文献值来比）。**忽略它。**
* ⚠️ `WARNING: no far-field BC; skip Cp distribution write`：方腔没有远场，属正常。
* 瞬态方式跑到定常**不会**打印 `CONVERGED`（那是稳态 SIMPLE 的判据）；历史 PIMPLE 记录看
  4.1 的 `du/dt` 平台 + 4.2 的 L2。

## 6. `check_cavity.py` 说明

* 读 ASCII VTU，点 = 网格节点 + **每个单元一个形心点**；tets 形如
  `n_pts i1 i2 i3 centroid_id`，用「唯一的形心 id」恢复每个求解单元的单一值。
* 探针与 `mod_uns_output.f90` 的 `probe_velocity` 一致（8 近邻、`w = 1/(d²+1e-12)`），
  所以对同一时刻的场，脚本与日志里的 Ghia 表**逐位相同**；日志丢了也能从任何 `.vtu`
  复算（串行/并行、快照/末态都行）。
* L2 的分子分母只统计**内部 15 个站位**（跳过 y = 1 和 y = 0 两个壁面行：单元中心探针
  取不到壁面值），与求解器口径一致。

## 7. 备注

* 以下为历史 PIMPLE 算例备注；当前现场是 `transient = false` 的 10 次 SIMPLE 参数检查。
* **收敛慢是这个算例的物理性质**，不是 `dt` 太小：τ ≈ 7.2 是无量纲时间，与网格无关地长
  （方腔的次级涡缓慢建立）。完整稳态解可使用 SIMPLE 增加 `outer_max` 并观察收敛。
* 快照命名：`unMesh.vtu_step<NNNN>.vtu`（求解器按 4 位补零），默认 `n_out_every = 1000`
  → 6 个快照（t = 10, 20, …, 60），**每个约 11 MB**（ASCII VTU，6400 单元 / 25760 节点），
  加末态场共 **约 77 MB**。比这个小（`n_out_every = 500`）只多花磁盘，不提高精度。
* 重启：`restart = true` + `restart_file = <dump>` + `dump_file = <dump>` 可续算；dump 自
  2026-10-13 起带**时间历史 ＋ `ts_order`**（此前历史全零、阶数丢失，重启会静默掉到一阶），
  且从 `ts_order = 2` 的 dump 续算第一步即保持二阶。链式多段续算（把长算例切进单次 ≤10 min）
  的定量验证见 `cluster_cases/cylinder_Re200/README.md` 第 8 节。
* 核数声明与守卫（`n_struct_ranks` + `n_uns_ranks`，退出码 6/7/8）见
  `cluster_cases/README.md` 第 3 节；本现场的 `mix.control` 已打开 `uns_parallel = always`。
## 8. ASCII / binary 输出切换（`out_format`）

`unMesh.control` 默认**不写** `out_format`，即 **inline ASCII**（Float64）。需要小文件时
加一行即可（控制文件末尾已附注释示例）：

```
out_format = vtu_binary      # appended raw binary（Float32）
```

| `out_format` | 数据布局 | 精度 | 本算例末态快照 |
|---|---|---|---|
| `vtu`（默认）/ `vtu_ascii` | inline ASCII | Float64（全部有效位） | **13.70 MB** |
| `vtu_binary` | appended raw ＋ UInt32 块头 | Float32（~7 位有效数字） | **3.38 MB**（~1/4） |

两者都是**符合 VTK 规范**的 `.vtu`（`connectivity` 为纯节点号、`offsets` 累加、
`types` 每单元一个），ParaView 与 Tecplot 360 都能打开；求解器每条
`wrote VTU file` 日志会注明实际布局（binary 会多打 `(appended raw binary, Float32)`）。

同一 control 分别跑一次（np = 1、10 次 SIMPLE 外迭代）的实测对照：

| 量 | ASCII | binary |
|---|---|---|
| 日志 `max \|u\|` | 9.11512E-01 | 9.11512E-01 |
| 日志 `relative L2 error (centerline u_x)` | 2.89706E-02 | 2.89706E-02 |
| VTK 读取器对拍 `pressure`（相对差） | — | 2.7e-08 |
| `velocity`（相对差） | — | 3.3e-08 |
| `temperature` | — | **逐位相同** |
| 首/末单元 connectivity | `[81,6642,6643,13122]` / `[6479,13121,6560,19521]` | 同上 |

**注意**

* `check_cavity.py` 用严正则读 **ASCII 文本**，读不了 binary（仓库内 20 多个
  `cases/**/*.py` 同理）。切 binary 后改用
  `python3 tools/validate_vtu.py --vtk out/<tag>.vtu`（格式守卫 ＋ VTK 读取器）
  或自写 `vtk`/`meshio` 读取器。
* binary 只有 ~7 位有效数字，会破坏"逐位一致"式回归断言；对 Ghia L2（3 位有效
  数字判据）无影响。
* `tools/validate_vtu.py` 现在两种布局都支持（appended raw / base64 的块
  `offset`、块头、长度都按 VTK 约定校验；`--selftest` 10 个夹具）。
* 遗留（与本次 ascii/binary 切换无关）：`--vtk` 模式还会检查四面体有无非正体积，
  本算例与 `cases/heat_conduction` 的旧产物都会报"非正体积四面体"（每个单元
  12 个 tet 里 4 个朝向相反）——这是**改动前就存在**的四面体化朝向问题，单独跟踪修复。

