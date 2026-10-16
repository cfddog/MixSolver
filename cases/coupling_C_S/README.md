# coupling_C_S 可压缩（struct）– 固体（uns）共轭传热（CHT）

> 日期：2026-10-16　求解器：`bin/mixnsolver_mpi -np 2`（rank0=struct，rank1=uns）
> 主题：**本仓库第一个"可压缩流体 ↔ 纯固体"共轭传热算例**——Ma=3 层流平板
> （struct，OpenCFD-EC Riemann 求解器）与 100 mm 不锈钢板（uns，phase-16
> `CZ_SOLID` 纯导热）在 y=0 面上弱耦合。界面类型
> `IFACE_COMP_FLUID_SOLID`（phase 17），交换配方 `T:dirichlet q:flux`。
> 本轮目标：**打通 CHT 通道 + 量化它的收敛代价**（不是跑出收敛的稳态——
> 见 §7、§8：这个算例的物理暂态需要 O(100) 个耦合迭代）。
> 网格与边界条件来自用户交付件 `cluster_cases/coupling_C_S/`，逐位保留，
> 只按交付 readme 第 (3) 条把 zone 4 从 mass-flow-inlet 改成 wall（§4）。

## 1. 布局、网格与界面

| 角色 | 文件 | 说明 |
|---|---|---|
| 结构侧（流体） | `Mesh3d.x`、`bc3d.inp`、`bc3d.fvbnd`、`control.ec`、`mix.control` | mm 单位；81×101×2 节点（80×100×1 单元）；x 0–1000、y 0–500、z 0–100；**首层 1 µm**；无普通壁面（唯一壁面＝耦合界面） |
| 非结构侧（固体） | `unMesh.cas`、`unMesh.control` | mm 单位；6642 节点 / 3200 六面体（80×40×1）；x 0–1000、**y −100–0**、z 0–100；**首层 1 µm**；cell zone 2 = `solid`（不锈钢） |
| 界面配对守卫 | `check_iface_alignment.py` | 由两套网格**独立重建**界面面序并逐面比对（§5） |
| 解析锚 | `estimate_cht_steady.py` | 由算例控制文件推导层流平板解析稳态（§6） |

界面：两套网格共享同一个 y=0 面离散（80 个四边形，每个 12.5 mm × 100 mm，
合计 0.1 m²），**节点坐标逐点相同**——守卫脚本实测逐面质心最大偏差
`0.000e+00 m`、两侧总面积都是 `0.1000000000 m²`。两侧首层厚度都是 1 µm
（几何增长：struct 前 5 层 1.000 / 1.014 / 1.041 / 1.083 / 1.139 µm；
slab 前 3 层 1.000 / 1.075 / 1.231 µm）。⇒ 界面无需插值，且
$k_f/d_1$ 极大（§8 的根因）。

## 2. 物理与参考量（`mix.control`）

| 量 | 值 | 来源 / 说明 |
|---|---|---|
| `rho_ref` | 0.44124 kg/m³ | 由 p=1 atm、T=800 K 反推（p_inf = ρR T_ref = 101329 Pa，误差 0.004%） |
| `T_ref` | 800.0 K | 来流静温（`control.ec` 的 `T_inf=800`） |
| `L_ref` | 1.0e-3 m | **1 个网格单位 = 1 mm**（沿用 `couple_porous` 的 mm 约定，网格不缩放） |
| `u_ref` | 1700.95 m/s | = Ma·a_ref = 3·√(1.4·287.058·800) = 1701.04（差 0.005%） |
| Ma | 3.0 | `control.ec` |
| Re | 1000.0 | **每个网格单位（1 mm）**：μ = ρu·L_ref/Re = **7.505e-4 Pa·s** |

OpenCFD-EC 无量纲化：`mu* = 1/Re`、`p* = p/(ρ_ref u_ref²)`（来流
`p* = 1/(γMa²)`）、`T* = T/T_ref`、能量通量尺度 `ρ_ref u_ref³`（见
`mod_struct_driver.f90:131-136` 的 q_SI 换算注释）。由 `γ`/`PrL` 反推出的
**模型流体**（不是真实空气）：

    cp = γR/(γ−1)   = 1004.7 J/kg/K
    k_f = μ cp/PrL  = 1.0772 W/m/K        α_f = 2.43e-3 m²/s
    ν = μ/ρ         = 1.701e-3 m²/s       Re_L(L=1 m) = 1.0e6
    δ(Blasius, x=1 m) = 5.0 mm  ⇒ 1 µm 首层 = δ 内 5000 层

即：Re=1000/格心单位 把黏性放大到真实空气的 ~20 倍（800 K 空气 μ≈3.7e-5），
是本算例（交付 readme 指定 Ma=3、Re=1000）的既定设定。

## 3. 边界条件

| 侧 | 面 | bc 码 | 含义 |
|---|---|---|---|
| struct | i=1（x=0） | 5 = `Inflow` | **超声速入口**：`boundary_Farfield(FLAG_INLET)` 强制 ghost=来流（四条特征线全部进入） |
| struct | i=81（x=1000） | 6 = `Outflow` | **超声速出口**：纯内点外推 `FLAG_OUTLET` |
| struct | k=1,2（z 面） | 3 = `Symmetry` | 对称 |
| struct | j=1（y=0） | 8 = `generic:8` | **耦合界面** → phase-17 CHT：等温无滑移壁面，温度＝固体表面温度，并回传壁面热流 q_f |
| struct | j=101（y=500） | 6 = `Outflow` | 顶部＝自由来流，超声速外推 |
| uns | zone 4（y=−100 mm） | `wall` + `tbc 300 K` | 固体底面等温 300 K（**交付件原为 mass-flow-inlet，本算例已改为 wall**） |
| uns | zone 5（x=0 / x=1000） | `wall`（无 tbc） | 绝热侧壁 |
| uns | zone 6（z=0 / z=100） | `symmetry` | 对称 |
| uns | zone 7（y=0） | 界面（无 bc/tbc 行） | 由耦合层每次迭代施加流体侧壁面热流 q_f（Neumann） |

> 交付 readme 要求"来流采用超声速入口，出口为超声速出口"：`bc3d.inp` 里
> code 5/6 在 `mod_struct_bc.f90:33-36` 直接走 `boundary_Farfield` 的
> 强制入口/出口分支（`:182-186`），**本身就是**超声速入口/出口，
> 无需新 BC 类型。注意：**struct 侧没有普通壁面**——唯一的壁面就是耦合界面，
> 这正是 CHT 平板所要的拓扑。

## 4. 交付件与本地修正（可追溯）

`cluster_cases/coupling_C_S/` 的 4 个输入件里 3 个**逐位照抄**，只有 `unMesh.cas`
按交付 readme 第 (3) 条（"cas 设置成了质量入口，请更正为固体温度壁面"）改了
**3 处记录**：

| 文件 | md5 | 说明 |
|---|---|---|
| `Mesh3d.x` | `f08c4e8a900dfd81720a61cf587049e0` | 与交付件逐位一致（392728 B） |
| `bc3d.inp` | `233df405266512e3dac484c7ba5da1f8` | 逐位一致 |
| `bc3d.fvbnd` | `d3e582f2e13c8cbfc5facecad971605f` | 逐位一致 |
| `unMesh.cas` | 交付 `ffe81ed4aa0d56ab89db4d4a66553a4f` → 本算例 `878be0a66c12ff3661e41e0fa7cf1d50` | 仅 3 处记录（下表）；保留 CRLF，唯一长度差 = +25 B |

| 行号 | 交付原文 | 本算例 |
|---|---|---|
| 12957 | `(0 "Zone 4 80 faces 6281..6360, BC: bc-2 mass-flow-inlet = 20")` | `(0 "Zone 4 80 faces 6281..6360, BC: bc-2 wall = 3 (edited from mass-flow-inlet = 20, coupling_C_S)")` |
| 12958 | `(13 (4 1889 18d8 14 4)(` | `(13 (4 1889 18d8 3 4)(`（Fluent 类型码 0x14=20 → 3=wall） |
| 13040 | `(45 (4 mass-flow-inlet bc-2)())` | `(45 (4 wall bc-2)())` |

真正生效的边界是 `unMesh.control` 里的 `bc = 4 wall` + `tbc = 4 1 300.0`
（BC 分组只由控制文件决定，见 `mod_uns_bc.f90:155-262`）；改 cas 是为了让交付件
本身不再自相矛盾，也避免后续工具按 zone 名误判。其余 3 个控制文件
（`control.ec`、`mix.control`、`unMesh.control`）为本算例新写；
`bc3d.inc` 是 struct 初始化时由 `bc3d.inp` 生成的**运行产物**，不入库。

## 5. 界面配对守卫 `check_iface_alignment.py`

**为什么必须查**：生产耦合路径按**索引**配对，不做几何匹配
（`mod_mix_driver.f90`：struct 的 `Interface_List` 第 c 个面 ↔ uns 交换数组第 c 行；
几何匹配器 `match_interfaces` 只走测试路径）。两侧面序由**各自的代码**独立建立：

* struct（`mod_struct_grid.f90:1046-1069`）：先按 `bc3d.inp` 的 generic:8 块，
  块内 `a`（第一平面内索引）**外**循环、`b`（第二索引）**内**循环 ⇒
  j 面上按 **i 递增**（本网格 k 只有 1 层单元）；
* uns（`mod_uns_bc.f90:311-318`）：按 `m%f(:)` 面号递增，即 **cas 文件里面序**。

两者若不一致，耦合照样能跑、结果看上去也"合理"，但每个面都在和**错位的伙伴**
交换温度/热流——静默的物理错误。脚本因此从两个网格文件**独立重建**两种面序，
逐面比对，并对本算例的已知陷阱设硬门槛。

**实测（PASS）**：

    --- structured side (Mesh3d.x + bc3d.inp) ---
      Mesh3d.x header word = 1, nodes = 81 x 101 x 2 = 16362
        stored block X: range [0, 1000] mm / Y: [0, 500] / Z: [0, 100]
        struct j=1 plane: y = 0 mm
      bc3d.inp bc blocks: 6, code-8 (generic) blocks: 1
      code-8 block: i 1..81, j 1..1, k 1..2 -> face 2 (j- plane)
      struct interface quads (registration order) = 80
    --- unstructured side (unMesh.cas) ---
      unMesh.cas integer tokens: hexadecimal
        nodes = 6642; face zones = [3, 4, 5, 6, 7]
      uns interface quads (cas face order) = 80 (faces 12841..12920)
    --- face-by-face pairing ---
      worst per-face centroid mismatch : 0.000e+00 m (face -1)
      sum area struct = 0.1000000000 m^2, uns = 0.1000000000 m^2
      min normal dot product (struct vs uns) = +1.000000
      face   1: struct (6.25, 0.0, 50.0) / uns (6.25, 0.0, 50.0) (mm)
      face  80: struct (993.75, 0.0, 50.0) / uns (993.75, 0.0, 50.0) (mm)
    RESULT: PASS (0 failures, 0 warning(s))

结论：**面序 1:1 吻合（最大偏差 0）**，可以信任本轮耦合结果。
脚本自带守卫（遇到就报错，不静默给数）：code-8 块必须唯一且是 j 面；
Mesh3d.x 三个坐标块跨度必须是 1000/500/100 mm（防换错网格）；cas 的
`(13 ...)` 十六进制区间必须与 `(0 "Zone ...")` 十进制注释逐条吻合（防进制误判）；
zone 4 必须叫 `wall`（即上面的交付件缺陷）、zone 7 必须是 `interface`；
两侧各 80 面、面积必须都是 0.1 m²。

## 6. 解析锚 `estimate_cht_steady.py`

流体膜 + 平板串联的稳态（$q=h(T_{aw}-T_w)=(k_s/L_s)(T_w-T_{base})$，
$T_{aw}=T_\infty(1+r\frac{\gamma-1}{2}Ma^2)$，$r=Pr^{1/3}$）：

| 取 h | h [W/m²K] | T_w [K] | q [W/m²] | Q=∫q dA [W] | 增益 G=h·R_s |
|---|---|---|---|---|---|
| $h_{avg}=0.664\sqrt{Re_L}Pr^{1/3}k/L$ | 635.1 | **1717.1** | **2.296e5** | **2.296e4** | 3.92 |
| $h_{loc}(x=L)$ | 317.5 | 1477.7 | 1.908e5 | 1.908e4 | 1.96 |
| $k_f/d_1$（首格电导） | 1.077e6 | 2078.3 | 2.881e5 | 2.881e4 | **6649.5** |

$R_s=L_s/k_s=6.1728e-3$ m²K/W，$T_{aw}=2078.6$ K（$r=\sqrt{Pr}$ 时 2004.8 K），
$T_{base}=300$ K。**CHT 的难点全在这张表里**：G 的物理上界只有 ~4，但只要界面的
散热/受热由**首格**的电导 $k_f/d_1$ 决定（1 µm 首层！），G 就放大到 6600。
脚本的 `--measure-resid R0 R1 --relax w` 还能由两次残差反推实测增益
（$G=(1-r_1/r_0)/w-1$），用于 §7 的判读。

## 7. 冒烟运行记录（3 段，全部串行单作业，每段 `timeout 600`）

机器：WSL2 笔记本（gfortran 13.3.0 + OpenMPI 4.1.6）。运行现场
`.dbg/coupling_C_S/`（运行产物不入库）。命令：

    cd .dbg/coupling_C_S && timeout 600 mpirun -np 2 --oversubscribe ../../bin/mixnsolver_mpi

一段 20 个耦合迭代的墙钟 ≈ **6.5 min**（≈20 s/迭代 = struct 1000 步 LU-SGS
+ uns 50 步 SIMPLE）。

### 7.1 初始化取证（各段日志一致）

    [struct] Registered 80 gridgen generic:8 coupling interface cell-face(s).
      [struct] block=1 face=2 nfaces=80
               centroid : 5.0000E+02 0.0000E+00 5.0000E+01   normal (0,-1,0)
               area     : 1.0000E+05                <- mm^2 = 0.1 m^2
    [uns] mesh coordinates scaled by 1.000E-03
      Registered interface zone 7 : bc-5   80 faces
              centroid (0.5, 0, 0.05) m   net normal (0,+1,0)   area 1.0000E-01
    --- Cell (volume) zones ---
      id condition  name  ncells  type   source
       2            solid  vc-2   3200  solid  control   <- cell zone 显式指定
    --- Interface dispatch table (phase 11) ---
      faces 80  comp-fluid<->solid
            recipe:T:dirichlet q:flux                     <- CHT 分派命中
      coupling interface: 80 faces -> BC_INTERFACE group
    --- Boundary conditions ---
      zone 4 wall faces: 80     fixed-T wall : 3.000E+02
      zone 5 wall faces: 80 / zone 6 symmetry faces: 6400
      zone 0 coupling-interface faces: 80
    [uns] CHT interface: pure-solid faces default to adiabatic (Neumann) mode

⇒ 80 面注册（两侧外法向相反，正常）、CHT 类型分派、固体 cell zone、底面 300 K
定温、界面走 Neumann，全部就位。

### 7.2 run 1：`iface_t_relax = 0.03`（cht_solid 的经验值）

| 耦合迭代 | 3 | 4 | 5 | 6 | 7 | 8 | 20 |
|---|---|---|---|---|---|---|---|
| `|r|max` [K] | 11.07 | 114.60 | 123.33 | 223.41 | 144.31 | 295.80 | **576.6** |

迭代 3→8 的逐迭代增长因子 ≈ **1.93**（按 $|1-w(1+G)|=1.93$ 反推 $G\approx97$
⇒ w = 0.03 会不稳，要求 w < 0.020）。迭代 20 时：struct 界面 T = 434.0 K、
p = 1.036e5 Pa、$\bar q=-1.15e8$ W/m²（Σ(qA) = −1.15e13 W）；uns 侧
Twall = 762.8–855.2 K、`q_iface_in` = −1.15e7 W、`q_base_in` = −8.2e3 W。
全程无 NaN、无 abort，20 迭代处正常存档
（`flow3d.dat` / `unMesh_restart.dat` / `couple_state.dat` 都落盘）。

### 7.3 run 2：`iface_t_relax = 0.005`（≈ 0.5/(1+G)）

| 耦合迭代 | 3 | 4 | 5 | 6 | 7 | 8 | 20 |
|---|---|---|---|---|---|---|---|
| `|r|max` [K] | 11.07 | 114.88 | 126.66 | 232.27 | 158.91 | 311.07 | **647.2** |

**与 run 1 的残差序列几乎相同（差 < 5%）——尽管 w 差了 6 倍。**
迭代 20 时：struct 界面 T = 326.4 K、p = 1.0285e5 Pa、$\bar q=-8.37e7$ W/m²；
uns 侧 Twall = 818.8–872.0 K、`q_iface_in` = −8.37e6 W、`q_base_in` = −8.6e3 W。

### 7.4 判读：这里的残差是**物理漂移**，不是松弛不稳定

w 从 0.03 降到 0.005 后轨迹几乎不变 ⇒ 残差 $|r|=|T_{recv}-T_{app}|$ 的演化由
**物理暂态**主导（固体表面被 ~1e8 W/m² 级的界面前期热流推着飞升、流体首格从
800 K 向壁温"淬火"），而不是由 T_app 的松弛反馈主导。因此：

* 本算例 **w 在 [0.005, 0.03] 区间内不可辨识**，两者都不是"错"值；
* 单靠调 w **无法**让它收敛——要收敛必须压下**物理暂态**（§8）。
  （cht_solid 的注记描述过同一现象："残差被物理漂移主导时 Aitken 估计会一直
  压到下限"；这里给出的是它的 w 不敏感实验证据。）

### 7.5 重启续跑（joint restart，phase 10）——发现一个 CHT 限制

现场 `.dbg/coupling_C_S/restart/`（拷 run 2 的 3 个存档 + `couple_restart=1`、
`n_couple=40`，其余不变）：

    [struct] joint restart: Iflag_init=1 forced
    [struct rank 0] restart from coupling iter 20
    Field dump read from: unMesh_restart.dat        (ncells=3200 ts_order=0 ver=1)
    [uns rank 1] restart from coupling iter 20

⇒ 两侧确实从 iter 21 续跑（没有重头、没有挂死）。**但 CHT 的界面状态不在存档里**：
重启后第一个交换的界面热流丢失（q = 0），uns 的稳态导热解立刻把平板拉回
"底面 300 K + 界面零通量"的均匀 300 K——第一个 SIMPLE 块的 `dT_max = 6.12E+02`
正说明读回来的场本来是热的（~900 K），随后 1–2 个外层迭代就被抹平。于是
迭代 25 时 `|r|` 掉回 22.6 K、Twall = 400.9–446.8 K：**这不是 run 2 的延续**。

**结论/待办**：CHT 长跑目前**不能靠分段重启**——要么一次跑完，要么先把界面
状态（q_f / 壁温 / T_app）纳入存档（`write_field_dump` 或 `couple_state.dat`）。
（对照：grid_BC 的低速–多孔重启是逐位可续的，见 worklog 阶段 10。）

### 7.6 重启段第二半（迭代 26–40）的走向

    iter 25: |r|max = 22.6 K, Twall = 400.9–446.8 K, q_iface_in = +1.02e7 W
    iter 40: |r|max = 426.8 K (rms 247.1), Twall = 452.1–648.7 K

（该段从"被重置的冷平板"重新走同一段物理暂态：Twall 400 → 649 K、|r| 22.6 →
426.8 K，struct 界面 p 稳定在 1.03–1.04e5 Pa，无 NaN。数值见 §7.7 的日志摘录。
**注意**：|r| 又把 ×~2/迭代 的"物理增长"重演了一遍——这正是 §7.4/§8 的结论
（不是松弛不稳，而是冷平板 soak 的物理反馈）在另一条轨迹上的复现。）

### 7.7 重启段终态（迭代 40）

    [struct rank 0] diag iter 40 struct SI: mean u_x, mean p(abs), mean T
                   = 0.0000E+00  1.0357E+05  3.1060E+02
    [struct rank 0] diag iter 40 CHT q_f: Sum(q*A) [W], mean, max
                   = -5.3977E+12  -5.3977E+07  7.9818E+07
    [uns rank 1]    diag iter 40 Twall min/max, q min/max
                   = 4.5205E+02  6.4870E+02  -7.9818E+07  -4.4345E+07
    [uns rank 1]    diag iter 40 CHT energy: q_iface_in, q_base_in, resid
                   = -5.3977E+06  -3.3697E+03  -5.4011E+06

⇒ 与 run 1/2 一样：界面通量量级 1e7–1e8 W/m²（物理稳态 2.3e5），
`q_base_in` 只有 −3.4e3 W（平板仍在储热）⇒ **20 个迭代只走完暂态的很小一段**，
与 §8 的机理分析一致。

## 8. 机理：为什么这个算例这么"硬"

1. **首层 1 µm ⇒ 界面增益被放大 ~1700×**。struct 的壁面热流按
   $q_{SI}=k_f(T_1-T_{wall})/d_1$ 计算（`mod_struct_driver.f90:131-136`），
   $d_1=1$ µm ⇒ $h_{eff}=k_f/d_1=1.077\times10^6$ W/m²K，而真实层流
   $h\approx635$ W/m²K。
2. **冷启动落差 500 K**：流体初值 800 K（`Iflag_init=0` 均匀来流）、固体
   `init_t=300` ⇒ 迭代 2 就出现 $\bar q=4.31\times10^8$ W/m²（物理稳态
   $2.296\times10^5$ 的 **1876 倍**）、$\Sigma(qA)=4.3\times10^{13}$ W。
   平板按稳态导热响应这个通量需要 $R_s\bar q\approx6\times10^5$ K 的表面温升
   ⇒ 两侧互相追逐、$|r|$ 逐迭代 ×~2、通量符号来回翻转（记录里 ±1e8 W/m²）。
3. **流体侧的热"淬火"极慢**：要把 1 µm 首格从 800 K 拉到近壁温，需要整条热
   边界层建立。实测首格温度：800 K（初值）→ 600 K（累计 ~1e3 步）→ ~400 K
   （~5e3 步）⇒ 时间常数 ~7e3 步；要落到 $q\approx h(T_{aw}-T_w)$ 的物理量级
   （即 $T_1-T_{wall}\lesssim1$ K）需要 **O(1e5) 步 ≈ O(10²) 个耦合迭代**。
4. 因此本轮 3 段各 20 迭代（≈2e4 步）**注定只走出暂态前段**——这不是代码缺陷，
   而是该网格/Re/初值的物理代价。同时这些量都"正常"：界面压力在
   0.77–1.04e5 Pa 之间摆动（无系统性漂移）、uns 侧 `|u|max = 0`（固体无滑移
   被强制）、`q_iface_in + q_base_in` 与 `q_iface_in` 同量级（平板在储热，
   未到稳态）、全程无 NaN、无 MPI 挂死。

## 9. 复现

    # 1) 界面配对与资产校验（必须 PASS，否则别信后面的数）
    python3 check_iface_alignment.py
    # 2) 解析锚（表 §6）
    python3 estimate_cht_steady.py
    python3 estimate_cht_steady.py --measure-resid 114.88 126.66 --relax 0.005
    # 3) 耦合冒烟：bin/mixnsolver_mpi 不带参数，按固定文件名从 cwd 读输入
    mkdir -p /tmp/cht && cd /tmp/cht
    cp <case>/Mesh3d.x <case>/bc3d.inp <case>/bc3d.fvbnd <case>/unMesh.cas \
       <case>/control.ec <case>/mix.control <case>/unMesh.control .
    timeout 600 mpirun -np 2 --oversubscribe <repo>/bin/mixnsolver_mpi
    #    20 迭代 ≈ 6.5 min；看日志里
    #      [struct rank 0] diag iter N  CHT q_f / CHT relax
    #      [uns rank 1]    diag iter N  Twall min/max / CHT energy
    # 4) 换 iface_t_relax 做敏感性：改 mix.control 一行即可（见 §7.3/§7.4）

## 10. 待办 / 建议

1. **拿到可与解析锚（$T_w\approx1717$ K、$q\approx2.3\times10^5$ W/m²）比较的
   稳态**（首要）：
   * (a) 集群上一次性跑 O(100) 个耦合迭代（≈2–3 h）。**注意 §7.5：目前不能
     分段重启**（界面热流不在存档里）；
   * (b) 或先把物理暂态压小：`n_struct_steps` 提到 5000–10000，让流体的壁面
     响应趋于物理值（实测增益 97 → 理论 ~4），此时可把 `iface_t_relax` 放大到
     0.1–0.2，所需迭代数大幅下降（每迭代更贵，但总数少）。
2. **界面处理改进**（可选，代码侧）：uns 侧现在放宽的是 T、而施加的是 q；
   可改为对 q 亚松弛；或把 Aitken 下限 0.02（`mod_mix_driver.f90:310`）
   降到 ~1e-3，让自适应模式在本算例也可用。
3. **重启补全**（§7.5）：把界面状态（q_f、壁温、T_app）纳入
   `write_field_dump` / `couple_state.dat`，CHT 长跑才能分段续跑。
4. struct 顶部 j=101 目前是"超声速出口外推"：Ma=3 自由来流下应无反射；
   若以后出现反射，把 `bc3d.inp` 第 6 行 code 6 改成 4（farfield）即可。
5. 网格侧（需要用户决定）：1 µm 首层对层流（δ=5 mm）其实过细（δ 内 5000 层），
   若允许改网格，首层放宽到 5–20 µm 可把界面增益降 5–20 倍。
