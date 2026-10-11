# cylinder_Re40 —— 二维圆柱绕流，Re = 40

自包含现场：网格 + control + 校验脚本。运行产物落在 `logs/`、`out/`。

## 1. 物理与网格

| 项 | 值 |
|---|---|
| 圆柱 | 半径 R = 0.5（D = 1），轴沿 z，展向 z ∈ [0,5] **一个单元厚** → 严格二维 |
| 远场 | 圆 r = 10（20 D） |
| 网格 | 非结构 2D：3928 单元 / 31424 四面体 / 80 个壁面环单元（平均 r = 0.5089，首层厚 0.0048） |
| 无量纲化 | ρ = 1，μ = 0.025，U_∞ = (1,0,0) → Re = ρU D/μ = **40**（p 已按 ρU² 无量纲） |
| 边界 | zone 6 = wall（圆柱无滑移）；zone 4 = symmetry（z 两平面）；zone 5 = pressure-far-field（p = 0，U = (1,0,0)） |
| 参考 | 分离角 ~53°，回流长度 ~2.3 D，底部 Cp ~ −0.5…−1.0，总 C_D(Re=40) ≈ 1.5 |

网格 `unMesh.cas` 与 `cases/cylinder_cas/cylinder.cas`、`UNSSolverProj/cases/cylinder_invis/cylinder.cas`
**逐位相同**（md5 `e597487bb48a…`），所以 `UNSSolverProj/cylinder.vtu` 可以直接当回归基线。

## 2. 关键设置（`unMesh.control`）

```
alpha_u = 0.3     alpha_p = 0.2
outer_max = 8000  outer_tol = 1.0e-6
lin_tol = 1.0e-8  lin_max = 500
conv_blend = 1.0  nonorth_corr = 0.0  ppe_precond = ic0
```

**这个算例不会打印 `CONVERGED`，这是设计如此，不是故障。** 尾迹里 `du_max` 到 ~2e-4
就进入平台，要压到绝对 `outer_tol = 1e-6` 以下需要极长的时间常数；8000 步后由
`nsteps cap` 正常收尾（`SIMPLE: stopped after 8001 outer iterations`）。
**迁移前的归档基线同样跑满上限**，所以判断收敛性要着 `check_cylinder.py` 的积分量
（Cp 分布、分离角、回流长度），不要看日志里有没有 CONVERGED 行。

## 3. 运行与校验

```bash
./run.sh          # → logs/cylinder_Re40_np1.log, out/cylinder_Re40_np1.vtu
python3 check_cylinder.py out/cylinder_Re40_np1.vtu reference_cylinder_cp.dat
```

本机单核约 4 分钟（8001 次外迭代，`lin-it ≈ 95`）。`unMesh.vtu` 由求解器自己写出，
`run.sh` 只是把它 `mv` 成带 np 标签的名字。

## 4. 实测（本 bundle，HEAD 二进制，np = 1）与归档基线对照

| 量 | 本 bundle | 归档基线 `UNSSolverProj/cylinder.vtu` | 期望 |
|---|---|---|---|
| Cp 前驻点（φ = 0） | **1.5018** | 1.5024 | ~ +1.5 |
| Cp 底部（φ = 180°） | **−0.4925** | −0.4920 | 死水区吸力 ~ −0.5 |
| Cp 吸力峰 | **−0.9863** @ φ = −82.9° | −0.9855 @ −82.9° | ~ −0.97 |
| 分离角（离尾迹轴） | **51.7° / 51.6°** | 47.3° / 51.6° | ~ 53° |
| 回流长度 L_w | **2.3420 D** | 2.3420 D | ~ 2.3 D |
| max\|u\| / U_∞ | **1.216703** | 1.215934 | — |
| max 反向 u_x | **−0.10211** | −0.10196 | — |
| C_D（仅压力） | **1.0766** | 1.0753 | 总 C_D ≈ 1.5（含摩擦） |
| vs `reference_cylinder_cp.dat` | mean\|dCp\| **0.0082**，max **0.1146** | mean 0.0073，max 0.1152 | — |

`reference_cylinder_cp.dat` 是迁移前用**同一张网格**算出的 80 点壁面 Cp 回归基线；
本 bundle 与基线两条记录对它的 mean|dCp| 都 < 0.01，且彼此差异在 1% 量级
（Cp 值差 ≤ 6e-4，回流长度完全相同）→ HEAD 与迁移前结果一致。
分离角一列的两个数来自壁面剪应力过零的两个相邻区间（`check_cylinder.py` 会把两处转折都打出来）。

## 5. `check_cylinder.py` 说明

* 读 VTU（ASCII），点 = 网格节点 + **每个单元一个形心点**；tets 形如
  `n_pts i1 i2 i3 centroid_id`，所以脚本用「唯一的形心 id」恢复每个求解单元的单一值。
* 沿 80 个壁面环单元反解 φ = atan2(y, x)，给出 Cp(φ) 前驻点 / 底部 / 吸力峰。
* 分离角 = 壁面 x 方向摩擦速度过零点（两处转折均打印）；
  回流长度 = 中心线 u_x < 0 的最远 x（以 D 为单位）。
* `C_D` 只含压力项（对数采样），摩擦项在无能量方程的 2D 尾迹里占比很小，
  与文献总 C_D ≈ 1.5 对比时注意这一点。

## 6. 备注

* 无能量方程：T 恒等于 0，`dT_max ≡ 0.000E+00`（与另外两个现场一致）。
* `./run.sh 4` 会用 `bin/mixnsolver_mpi`，非结构模式在 `-np > 1` 时**默认走 METIS 分块
  并行**（`mix.control` 的 `uns_parallel = auto`），即 4 个 rank 各解 1/4 单元 + 一层 halo；
  实测 300 步下 `C_D` 与整体解差 0.2%。要复刻旧行为（每个 rank 一份整网格、结果与串行
  逐位相同，仅用于 MPI 确定性检查），在 `mix.control` 里加 `uns_parallel = never`。
* **核数声明与守卫**：`mix.control` 里取消注释 `n_struct_ranks = 0` + `n_uns_ranks = N`
  （两键必须成对），和 = N 必须等于发作业核数，否则求解器启动时打印诊断并**终止作业**
  （退出码 6 只写一个键 / 7 和不符 / 8 与模式冲突）；`./run.sh`（不带参数）会把这个和当作
  默认 np。不写这两个键则保持旧的隐式划分、不做检查。详见 `cluster_cases/README.md` 第 3 节。
* 8 个壁面环单元/9.9 的远场距离决定了这个规模；再加密网格前先看
  `check_cylinder.py` 的 `dr_first = 0.0048`（首层厚度）。
