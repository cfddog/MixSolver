# cluster_cases —— 五个可直接上集群的「干净现场」

本目录把五类物理整理成**自包含、零参数、可整体搬走**的算例现场：

| 现场 | 物理 | 求解路径 | 规模 |
|---|---|---|---|
| `channel_lowspeed/` | 平面 Poiseuille 槽道，Re_H = 100（低速） | 稳态 SIMPLE | 30000 单元 |
| `cylinder_Re40/` | 二维圆柱绕流稳态尾迹，Re = 40 | 稳态 SIMPLE | 3928 单元 |
| `bj_porous/` | Beavers-Joseph 流体/多孔界面，α = 0 / 1 / 2 | 稳态 SIMPLE | 4000 单元 |
| `cylinder_Re200/` | 二维圆柱绕流**非稳态**卡门涡街，Re = 200 | **瞬态 PIMPLE** | 3928 单元 |
| `cavity_Re1000/` | 二维顶盖驱动方腔，Re = 1000 | **瞬态 PIMPLE** | 6400 单元 |

后两个现场走的是**瞬态 PIMPLE** 路径，而该路径**只存在于 MPI 驱动里**
（`mod_uns_driver_mpi -> pimple_run_mpi`）：串行 `bin/mixnsolver` 只会跑稳态 SIMPLE 并
**静默忽略** `transient = true`。因此它们在自己的 `mix.control` 里写了
`uns_parallel = always`，`run.sh` 也**总是**调 `bin/mixnsolver_mpi`（即使 `np = 1`）。
踩过的坑见第 6 节第 8–9 条，数字见第 5 节，细节（边界条件、判据、代价）在各现场自己的
`README.md`。

每个现场只放**输入件**（网格 + control + 运行/校验脚本），不放运行产物：产物一律落到
各自的 `out/` 与 `logs/`（`*.vtu`、`*.log` 已在 `.gitignore` 里屏蔽）。

```
cluster_cases/
├── README.md          本文件
├── pack.sh            打包脚本（tar.gz；--with-repo 可连带源码树）
├── channel_lowspeed/  低速槽道     → README.md, run.sh, plot_channel.py, variants/
├── cylinder_Re40/     圆柱 Re=40   → README.md, run.sh, check_cylinder.py, reference_cylinder_cp.dat
├── bj_porous/         BJ 界面      → README.md, run.sh, plot_bj.py, gen_bj.py, variants/
├── cylinder_Re200/    圆柱 Re=200  → README.md, run.sh, check_cylinder_unsteady.py（瞬态）
└── cavity_Re1000/     方腔 Re=1000 → README.md, run.sh, check_cavity.py（瞬态）
```

---

## 1. 二进制怎么跑（零参数自调度）

`src/main_dispatch.f90` 让可执行文件**不吃任何命令行参数**，只读当前工作目录下的三个文件：

| 文件 | 作用 | 缺失后果 |
|---|---|---|
| `mix.control` | 自调度模式：`mode = uns` / `struct` / `coupled` | 报错退出 |
| `unMesh.cas` | 非结构网格（`mode = uns` 必需） | 报错退出 |
| `unMesh.control` | 非结构求解器控制文件 | 报错退出 |

结果固定写 `unMesh.vtu`（以 `.cas` 的 basename 命名），所以 `run.sh` 跑完会 `mv` 改名归档。
**旧的 `<cas> <control>` 命令行接口（`uns_solver`）已删除**；仓库里残留的旧 README 提到
它时，以本节为准。串行二进制 `bin/mixnsolver` 与 MPI 版 `bin/mixnsolver_mpi` 都读同一组文件。

`mix.control` 里还可选写一行非结构并行路径开关（详见第 3 节）：

```
uns_parallel = auto | always | never      # 默认 auto；-np>1 时 auto 即 METIS 分块并行
```


---

## 2. 集群侧编译

```bash
git clone <repo> && cd MixNSSolver
make all -j1    # 必须带 all：裸 make 是空操作（第 2 节的坑）；初次构建要单线程，模块依赖靠 .d 落地
make all -j     # 之后随便并行增量编译
# → bin/mixnsolver          串行可执行（但注意：它是用 mpif90 链接的，见下）

make mpi -j1    # → bin/mixnsolver_mpi   （mpif90 + -DHAVE_MPI）
```

> **坑：光敲 `make` 什么都不会发生。** Makefile 的第一条规则是目录规则
> （`$(ALLDIRS): mkdir -p $@`，Makefile 第 184~186 行），所以 GNU make 的默认目标就是
> `build/ser/mod`，裸 `make` 只会 `mkdir -p build/ser/mod` 然后报
> `make: 'build/ser/mod' is up to date.` 就退出。**必须写 `make all`**（或 `make mpi`）。

* 依赖：`gfortran`（`-std=f2008 -cpp`）、`mpif90`、`make`、`python3`（可选，用于校验脚本）。
  **串行目标 `bin/mixnsolver` 也是用 `$(MPIFC)` 链接的**：`src/structured/` 与
  `src/coupling/`、`src/main_dispatch.f90` 无条件调用 MPI（跨组交换天生多进程），
  所以这两个目录在串行树里同样用 mpif90 编译；因此 `make all` 也要求 `mpif90` 存在
  （没有 MPI 的机器跑不了这个构建）。
* **ParMETIS 不参与链接**：`lib/parmetis` 里的预编译档案是 MPICH ABI
  （`MPI_Comm = int`），本机 OpenMPI 4.x（`MPI_Comm = void*`）链接后会在 `MPI_Comm_rank`
  里崩。分区统一走**串行 METIS**（仓库自带 `lib/metis/libmetis.a`），`make mpi` 会自动
  链接它。
* `DEBUG=1 make all` → `-O0 -g -fcheck=all -fbacktrace`；`make clean` 清掉 `build/` 与 `bin/`。
* 构建产物 `bin/` 与 `build/` 不进 git，集群上必须自己编一次（或让 `pack.sh --with-repo`
  把源码一起带过去：它只打包 `Makefile src lib docs cluster_cases`，实测 4.0 MB / 126 项，
  而仓库里的 `cases/` 有 1.2 GB —— 别整个 `scp -r`）。

---

## 3. 每个现场怎么跑

```bash
cd cluster_cases/<scene>
./run.sh                 # np=1，仓库根 = $MIXROOT 或 ../..
./run.sh 4               # 4 个 rank（走 bin/mixnsolver_mpi）
./run.sh 1 /path/to/repo # 显式指定持有 bin/ 的仓库根
```

`run.sh` 的行为：`cd` 到自己所在目录 → `mkdir -p logs out` → 用 `tee` 同时把日志打到
屏幕和 `logs/` → 把 `unMesh.vtu` 移到 `out/`。**现场目录除 `unMesh.control` 外不会被
改动**（BJ 现场会在三个 α 变体之间切换 `unMesh.control`，跑完复位到 `bj_a1.control`）。

`run.sh` 里有一行容易忽略但很关键的设置：

```bash
export GFORTRAN_UNBUFFERED_ALL=y
```

libgfortran 自己管理 I/O 缓冲，重定向到文件后**要攒满 4 KB 才落盘**，所以刚启动的几分钟
`logs/*.log` 是空的（用 `stdbuf -oL` 也**没用**，它管不到 Fortran 自己的缓冲）。
加这一行以后 `tail -f logs/*.log` 才能实时看到 `it=` 行。现场里第 1 次排查「算例好像卡住了」
的教训就是这一条（见第 6 节第 3 条）。

**np > 1 的两种走法**（`mix.control` 里可选键 `uns_parallel`，默认 `auto`）：

| `uns_parallel` | `-np 1` | `-np N > 1` | 结果 |
|---|---|---|---|
| `auto`（默认） | 全 rank 各持全网格 | **METIS 分块**（`<cas 基名>.part.map`） | np=1 与串行逐位一致 |
| `always` | METIS 分块（1 块） | METIS 分块 | 实测 np=1 时与整体解**逐位相同** |
| `never` | 全 rank 各持全网格 | 全 rank 各持全网格 | np=N 与 np=1 **逐位相同**（只做确定性检查） |

即：**`mode = uns` 且 `-np N > 1` 现在会真正把网格分块并行求解**（rank 0 算 METIS 分区并
广播，每个 rank 只解自己的子域 + 一层 halo，收尾 rank 0 汇聚写 `unMesh.vtu`）；有效核心数
上限是**单元数**。要用旧的「每 rank 一份全网格」冗余模式（当初的 MPI 确定性检查用法），
把 `uns_parallel` 设成 `never`。结构化侧的分区单位不是单元而是**网格块数**（`Mesh3d.x`
的 block，`partation.dat` 存块→rank 映射），所以 `mode = struct` 的 `-np` 上限是块数。

**核数声明与检查（`n_struct_ranks` / `n_uns_ranks`，推荐写）**：在 `mix.control` 里把两侧
核数写出来，二者之和就是总核数，**必须**和发作业的核数（`mpirun -np N` /
`#SBATCH --ntasks=N`）一致，否则 rank 0 打印诊断并**立即终止作业**（防「算例按 3 核准备、
却用 `-np 2` 提交」这类静默误用）：

```ini
mode = uns                # 或 struct / coupled
n_struct_ranks = 0        # 结构侧核数
n_uns_ranks    = 2        # 非结构侧核数（0+2 = 2 ⇒ 必须 mpirun -np 2 / --ntasks=2）
```

| 情形 | 行为 |
|---|---|
| 两个键都写、和 = 实际核数 | 打印 `rank split : n_struct_ranks = 0 + n_uns_ranks = 2 = 2 rank(s) == launched 2 rank(s) [checked]`，正常跑 |
| 只写一个键 | `incomplete rank-count declaration`，退出码 **6** |
| 和 ≠ 实际核数 | `rank-count mismatch`（打印声明值 / 和 / 实际核数），退出码 **7** |
| 与模式冲突（`uns` 给 `n_struct_ranks>0`、`struct` 给 `n_uns_ranks>0`、`coupled` 不是 1+1） | 退出码 **8** |
| 两个键都不写 | 保留旧隐式划分（耦合 1 结构 + 其余非结构；单求解器全部 rank），只打印一行 `rank split : undeclared`，**不检查** |

* 模式分工：`mode = uns` 要求 `n_struct_ranks = 0` 且 `n_uns_ranks =` 总核数；
  `mode = struct` 要求 `n_uns_ranks = 0` 且 `n_struct_ranks =` 总核数。
* **耦合只能 2 核（1 结构 + 1 非结构）**：跨组交换只在两个组根之间点对点、非结构组内部
  没有集合通信，多余 rank 会卡死在交换里（实测 `-np 3` **挂死**；结构侧界面提取也只在单
  rank 成立，`Interface_List` 是全局表而 `Mesh(1)%Block` 是本地块）。所以耦合算例请用
  `-np 2`：声明别的组合、或未声明却 `-np > 2`，都会在启动时**报错退出**而不是挂死。
* 键名与取值都大小写不敏感；`run.sh` 的 `np` 缺省值会读 `mix.control` 里的声明（没写则仍
  是 1），三个现场的 `mix.control` 里都留了注释示例。


---

## 4. 五个现场速查

| 现场 | 输入件 | 耗时 / 终态（本机实测） | 产物 | 校验脚本 |
|---|---|---|---|---|
| `channel_lowspeed` | `unMesh.cas`(7.9 MB) + `unMesh.control` + `variants/pressure_outlet.control` | ~4 min；**CONVERGED @ 773** | `logs/channel_lowspeed_np1.log`，`out/channel_lowspeed_np1.vtu`（**~50 MB** ASCII） | `plot_channel.py` |
| `cylinder_Re40` | `unMesh.cas` + `unMesh.control` + `reference_cylinder_cp.dat` | ~4 min；**不打 CONVERGED**（设计如此，8001 步上限收尾） | `logs/cylinder_Re40_np1.log`，`out/cylinder_Re40_np1.vtu`（4.2 MB） | `check_cylinder.py` |
| `bj_porous` | `unMesh.cas` + `unMesh.control` + `variants/bj_a{0,1,2}.control` | 每个 α ~2 min；**三个 α 都 CONVERGED**（~800 步） | `logs/bj_a*_np1.log`，`out/bj_a{0,1,2}_np1.vtu` | `plot_bj.py` |
| `cylinder_Re200` | `unMesh.cas` + `unMesh.control`（**瞬态 PIMPLE**） | 本机 np=2 **24.6 min**（12000 步，0.125 s/步；dt=0.005，t=60） | `logs/*.log`；末态 `out/*.vtu` + 60 个快照（**~250 MB**） | `check_cylinder_unsteady.py` |
| `cavity_Re1000` | `unMesh.cas`(1.6 MB) + `unMesh.control`（**瞬态 PIMPLE**） | 本机 np=4 **18.2 min**（6000 步，0.182 s/步；dt=0.01，t=60） | `logs/*.log`（尾部有 Ghia 表）；末态 `out/*.vtu` + 6 个快照（**~77 MB**） | `check_cavity.py` |

---

## 5. HEAD 上的验证结论

| 现场 | 关键量 | 实测 | 参考 | 偏差 |
|---|---|---|---|---|
| channel | u_max / U_bulk | 1.5006（max\|u\| = 2.354425e-4 = 1.5012·U_bulk） | 1.5 | +0.04% |
| channel | dp/dx（x ∈ [200,295] 拟合） | −3.472454e-10 | −3.474307e-10 | **−0.05%** |
| channel | 出口剖面误差 | max 0.484%·U / L1 0.248%·U | 0 | — |
| cylinder | Cp 前驻点 / 底部 / 吸力峰 | 1.5018 / −0.4925 / −0.9863 | 归档 1.5024 / −0.4920 / −0.9855 | ≤ 6e-4 |
| cylinder | 回流长度 L_w | 2.3420 D（与归档**完全相同**） | ~2.3 D | — |
| cylinder | 分离角 | 51.7° / 51.6° | ~53° | — |
| cylinder | 对 `reference_cylinder_cp.dat` | mean\|dCp\| 0.0082，max 0.1146 | 归档 0.0073 / 0.1152 | — |
| bj | α = 0 界面速度剖面 RMS | **0.22%** | 两层 ODE 解析解 | — |
| bj | α = 1 | **0.94%**（uF1 = 1.1768e-2 vs 1.2117e-2） | 归档 1.09%（uF1 1.173e-2） | < 0.3% |
| bj | α = 2 | **1.20%** | — | — |
| cylinder_Re200 | 尾迹主频（全精度快照 FFT） | **St = 0.156**（场量 max\|u\| 独立给出 0.157） | 0.196（Re = 200） | −20%（网格仅 80 个周向单元） |
| cylinder_Re200 | t = 60 瞬时 C_D = F_x/(½ρU²A_ref) | **1.285**（F_x = 3.2064，A_ref = 4.9923） | 平均 ~1.34 | −4% |
| cylinder_Re200 | 起振 / 饱和 | 对称基流 t = 5–10（max\|u\| ≡ 1.40000）→ t ≈ 10 起振 → t ≈ 40 饱和 | — | — |
| cavity_Re1000 | 竖中心线 u_x 相对 L2 | **3.01e-2**（t = 60 瞬态）vs **2.89e-2**（同网格稳态 SIMPLE） | Ghia 1982 表 I | 残差 4% |
| cavity_Re1000 | 网格收敛（稳态路径） | 40×40 **9.12e-2** → 80×80 **2.89e-2** | — | 加密一倍降 3.2 倍 |

**五个现场在 HEAD 上都通过**（cluster channel 的 it=1 行与 2025-10-05 归档日志逐位相同；
圆柱 Re=40 所有积分量与归档基线一致；BJ 三个 α 都能跑到显式 CONVERGED；两个**瞬态**现场
各自给出可判读的物理量：Re=200 的涡街频率 + C_D，Re=1000 的 Ghia L2）。
瞬态现场的**局限**（圆柱的 St 比文献低 20%，来自网格分辨率；方腔 3e-2 的 L2 由网格而非
时间推进决定）在各自 `README.md` 里写清楚了 —— 当作**路径/流程**的验证是可靠的，
当作定量结果前先看那些局限。详细数字与图在各现场自己的 `README.md` 里。

---

## 6. 坑与经验（按踩坑代价排序）

1. **低速算例的 `outer_tol` 得用绝对值思维。** `channel_lowspeed` 的 `du_max` 是绝对速度
   增量而 U ~ 1.6e-4，常规的 `1e-6` 等于 0.64%·U，会在进口信息还没填满 L = 300 时就
   "收敛"。现场里用的是 **1e-10（≈ 6e-7·U）**。
2. **收敛判据里除了 `du_max` 还有 `dT_max`。** `src/unstructured/mod_uns_simple.f90`
   把 `mass-imbal AND du_max AND dT_max < outer_tol` 当收敛条件，所以**等温算例也会被
   被动 T 场卡住**：`bj_porous` 归档 control 只给了一面墙 `tbc`，T 一路漂到 53 K；
   两面墙都钉 300 K 也不行（多孔床把场拉到 ~88 K 稳态）。现场改成**绝热壁 + T 初值 0**，
   被动场恒等于 0（`dT_max = 0.000E+00`），速度解逐位不变。
3. **看不到进度不等于卡住：libgfortran 有自己的缓冲。** 重定向到文件后 gfortran 要攒够
   4 KB 才落盘，`logs/*.log` 会几分钟是空的，**`stdbuf -oL` 管不到它**。正确开关是
   `export GFORTRAN_UNBUFFERED_ALL=y`（五个 `run.sh` 里都加了，并写了注释）。
   当初「channel 好像卡住了」其实是这一条。
4. **判据不总是日志里的 `CONVERGED` 行。** `cylinder_Re40` 的尾迹 `du_max` 有个 ~2e-4 的
   平台，永远压不到 1e-6，所以它**按设计跑满 8001 步收尾**（旧基线也这样）。判断它的收敛性
   要用 `check_cylinder.py` 的积分量。
5. **旧日志的 `max |u|` / `p_min` 别当回归目标。** 迁移前 `uns_solver` 的
   `--- Solution summary ---` 块口径不同（channel 旧日志的 `max |u| = 2.575e-4 = 1.64·U_bulk`
   比解析最大值 1.5·U_bulk 还大），而新二进制不打印该块。
6. **零参数接口**：只有 `mix.control` + `unMesh.cas` + `unMesh.control`，旧的
   `<cas> <control>` 命令行接口已删除；结果文件名固定是 `unMesh.vtu`，所以 `run.sh` 跑完
   会 `mv` 改名，否则下一个现场会覆盖它。
7. **裸 `make` 是空操作**（默认目标是目录规则 `build/ser/mod`），必须 `make all`；而且
   `make all` 也要 `mpif90`（structured/coupling/main_dispatch 无条件调 MPI）。见第 2 节。
8. **瞬态 PIMPLE 只存在于 MPI 驱动里，串行二进制会静默忽略。** `bin/mixnsolver` 读到
   `transient = true` 不报错、不警告，直接按稳态跑完——你会以为在跑涡街，其实拿到的是
   一个定常解。所以两个瞬态现场的 `mix.control` 写了 `uns_parallel = always`，`run.sh`
   **固定**调 `bin/mixnsolver_mpi`（np = 1 也走 MPI 路径）。判据也不同：**瞬态没有
   `CONVERGED` 行**，看 `du/dt` 平台 / St / L2（见各现场 README）。
9. **同一张网格上「稳态能用的 BC」，瞬态可能每步放大 2 倍。** `pressure-far-field` 同时给
   远场 Dirichlet 速度**和**压力；稳态 SIMPLE 容忍这个过约束（`cylinder_Re40` 就这么跑），
   冷启动的 PIMPLE 不接受：实测 max\|u\| 按 1.6 → 4.3 → 8.7 → 14 增长，`dt = 5e-3` 与
   `dt = 1e-2` **完全一样**（所以不是步长问题）。换成 `velocity-inlet 1 0 0` 就稳
   （CFL ~ 0.07）。**两个现场之间别照抄边界条件**。
10. **日志里的速度只有 3 位有效数字，不能拿来做频谱。** `max|u| = 1.44E+00` 这种量化
    会把圆柱的 St 从 0.157 歪成 0.144（假象）。频率类信息必须从**全精度快照**
    （`*.step<NNNN>.vtu`）拿，而且频率分辨率由**时间窗长度**决定：加密快照只减小峰值
    偏置、不提高分辨率，见第 11 条。
11. **快照是 ASCII VTU，比想象的占地方**：圆柱 4.1 MB/个、方腔 11 MB/个。圆柱默认
    `n_out_every = 200` → 60 个 = 250 MB；方腔 `1000` → 6 个 = 77 MB。`out/` 在集群上
    可能被清或超配额，先算一下再提交。
12. **日志尾部那两个块是「专用诊断」，换个物理就是噪声。** `ghia_compare` 是**方腔**
    专用（竖中心线 u_x），圆柱算例只会打印 `Ghia reference: no table for this Re, skipped.`；
    反过来 `--- Separation and wake ---` / `Cp distribution` 是**圆柱**专用，方腔日志里
    会打印 `Top separation: not detected`、并拿 Re = 40 的圆柱文献值来比 —— **忽略**。
    两个瞬态现场真正的末态判据分别是 `check_cylinder_unsteady.py` 与 `check_cavity.py`。
13. **两处「看不见的」口径要记住**：`max|u|` 是**单元中心值**（方腔报 0.9116 而不是 1.0，
    因为 1.0 是**面**上的盖速）；`Reference area` 是求解器按形心算出的 4.9923 而不是
    名义 5.0（差 0.15%，自己复算 C_D 会差同样的比例）。

---

## 7. 打包与传输

```bash
cd cluster_cases
./pack.sh                       # → cluster_cases_<日期>.tar.gz        (实测 2.2 MB, 46 项)
./pack.sh --with-repo           # 另外 → cluster_cases_<日期>_repo.tar.gz (实测 4.0 MB, 126 项)
./pack.sh /tmp/scenes.tar.gz    # 指定输出路径（--with-repo 可同时给）
```

* 两个包都**只含输入件**：`logs/`、`out/`、`*.vtu`、`__pycache__` 一律排除，所以打出来
  的永远是干净现场；`sha256sum` 会打印出来方便在集群侧核对。
* `cluster_cases_*.tar.gz` 解开后就是 `cluster_cases/…`，可以直接盖进一个已有仓库，
  也可以配合集群上自己编的 `bin/` 单独使用。
* `..._repo.tar.gz` 是**源码白名单** `Makefile src lib docs cluster_cases`（构建只读
  `src/` 和 `lib/metis/libmetis.a`，见 Makefile 的 `src/%.f90` 与 `METIS_LIB`），
  顶层带 `MixNSSolver/` 前缀，解开就能 `make all`。仓库里 `cases/` 有 1.2 GB，
  **不要整个 `scp -r`**。
* 传输后建议 `sha256sum` 两边对一遍再解包。

