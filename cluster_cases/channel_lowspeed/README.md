# channel_lowspeed —— 低速平面槽道（Poiseuille），Re_H = 100

自包含现场：网格 + control + 校验脚本。运行产物落在 `logs/`、`out/`。

## 1. 物理与网格

| 项 | 值 |
|---|---|
| 网格 | 300(x) × 100(y) × 1(z) 六面体 = 30000 单元；x ∈ [0,300]，y ∈ [0,10]（壁面加密 dy_wall = 0.01、dy_c = 0.233），z ∈ [0,10] 一个单元厚 |
| 网格文件 | `unMesh.cas` = `cases/channel/channel_wall.cas`（md5 `cfcd790de8d2…`）。**注意不是** `cases/channel/channel.cas`（`5e94c34bfaff…`）——本现场用的是按 zone 拆过的版本 |
| 分区 | zone 2 = fluid；zone 4 = symmetry（z = 0 / z = 10）；zone 7 = wall（y = 0 / y = 10）；zone 5 = mass-flow-inlet（x = 0）；zone 6 = outflow（x = 300） |
| 流体 | 空气 300 K：ρ = 1.177 kg/m³，μ = 1.846e-5 Pa·s；**无能量方程** |
| 目标 | Re_H = ρ·U_bulk·H/μ = 100 → U_bulk = 100·1.846e-5/(1.177·10) = **1.568394e-4 m/s**，质量流密度 mdot = ρU_bulk = **1.8460e-4 kg/(m²·s)** |
| 解析参考 | u(y) = 6·U_bulk·η(1−η)（η = y/H），u_max = 1.5·U_bulk = 2.3526e-4；dp/dx = −12μU_bulk/H² = **−3.474307e-10**（每单位长度）；进口段 L_e ≈ 0.05·Re_H·H = 50 < 300 |

"准 SI 约定"：网格长度单位直接当米用，mdot 是单位面积质量流率。

## 2. 关键设置（`unMesh.control`）

```
alpha_u = 0.7    alpha_p = 0.3
outer_max = 8000
outer_tol = 1.0e-10      ← 绝对 du 阈值，故意取这么小，见下
lin_tol = 1.0e-8   lin_max = 300
conv_blend = 0.0   nonorth_corr = 1.0   ppe_precond = ic0
bc = 4 symmetry
bc = 7 wall
bc = 5 mass-flow-inlet 1.8460e-4
bc = 6 outflow
```

**为什么 `outer_tol` 是 1e-10 而不是常规的 1e-6**：`du_max` 是**绝对**速度增量，
而本算例 U ~ 1.57e-4，1e-6 就等于 0.64%·U —— 会在进口信息还没填满 L = 300 的域时
就宣布收敛。1e-10 ≈ 6e-7·U 才是这个低速算例该用的阈值（这是本现场最重要的一个坑）。

`variants/pressure_outlet.control` 是唯一的对照：**只把第 43 行**从
`bc = 6 outflow` 换成 `bc = 6 pressure-outlet 0.0`（其它 42 行完全相同），
用来对比两种出口处理（探索性，不是回归基线）。

## 3. 运行与校验

```bash
./run.sh          # → logs/channel_lowspeed_np1.log, out/channel_lowspeed_np1.vtu
python3 plot_channel.py out/channel_lowspeed_np1.vtu out/channel.png 100
```

`plot_channel.py <vtu> <png> [Re_H]` 会打印本节第 4 节那张表的所有数字并出 3 联图
（剖面 / 中线发展 / 压强）。单核到 CONVERGED 约 4 分钟（~773 步，每步线性 ~280 次）。
注意通道 VTU 是 **~50 MB ASCII**（3 万单元 + 6 万节点，网格大），搬运时留意磁盘。

## 4. 实测（本 bundle，HEAD 二进制，np = 1）

收敛历程（`logs/channel_lowspeed_np1.log`）：

* `it=1: mass-imbal = 5.113E-13, du_max = 2.840E-04, lin-it = 278`
  —— 与 2025-10-05 归档日志 `cases/channel/channel_Re100_out.log` 的 it=1 行**逐位相同**；
* 之后两者轨迹一致（`du_max` 每 50 步约 ×0.55），本 bundle **`CONVERGED after 773`**
  （`du_max = 9.950E-11`），归档基线在 790 步收敛 —— 差 2%，属残差抖动；
* `dT_max ≡ 0.000E+00`（无能量方程）。

`plot_channel.py`（`Re_H = 100`，出口 x = 299.5）：

| 量 | 实测 | 解析 | 偏差 |
|---|---|---|---|
| U_bulk（层平均）@ x = 50.5 | 1.567710e-04 | 1.568394e-04 | +0.04% |
| U_bulk @ x = 299.5（出口） | 1.565399e-04 | 1.568394e-04 | −0.19% |
| u_max/U_bulk（出口） | **1.5006** | 1.5 | +0.04% |
| max\|u\| 全域（VTU 单元值） | **2.354425e-04** @ (298.5, 5.116, 5.0) = 1.5012·U_bulk | 2.3526e-04 | +0.08% |
| 充分发展剖面误差（x = 299.5） | max 0.484%·U，L1 0.248%·U | 0 | — |
| dp/dx（拟合 [200,295]） | **−3.472454e-10** | −3.474307e-10 | **−0.05%** |
| 进口段长度 L_e | 55.5（5.55 H） | 估计 5.0 H | +11% |
| 进出口质量漂移 | +1.015%（U_bulk 沿程微降） | 0 | — |
| p 范围 | −1.317361e-07 … 0 | p_min = −L·\|dp/dx\| ≈ −1.04e-7 | 同量级 |
| 压强总升（进口→出口） | 1.1633e-07 | — | — |

结论：**出口剖面与解析 Poiseuille 解在 0.5%·U 以内，压强梯度在 0.05% 以内**，
进口段长度也和层流估计一致 —— 这个低速算例是用来验证「极低 Re + 质量流进口 +
外迭代绝对阈值」这条链路的。

**不要拿归档日志的 `max |u| : 2.57527E-04` 当回归目标**：那是迁移前旧求解器
`--- Solution summary ---` 块打印的另一个量（= 1.64·U_bulk，比解析最大值 1.5·U_bulk 还大），
新二进制的 `max|u|` 按上表是物理上正确的 1.5012·U_bulk。同理旧日志的
`p_min = −1.424427E-07` 与本 bundle 的 −1.3174e-07 属同一量级（旧跑未完全收敛）。

## 5. 备注

* 网格 `unMesh.cas` 用 `cases/channel/split_channel_zones.py` 从 Pointwise 原始
  `channel.cas` 拆出 zone 7 wall；这个拆分只影响边界分区，不影响网格点坐标。
* `./run.sh 4` 不再只是确定性检查：非结构模式 `-np > 1` 默认走 METIS 分块并行
  （`uns_parallel = auto`）；加 `uns_parallel = never` 可退回「每个 rank 一份整网格」的
  冗余模式（结果与串行逐位相同，仅用于 MPI 确定性检查）。
* **核数声明与守卫**：`mix.control` 里取消注释 `n_struct_ranks = 0` + `n_uns_ranks = N`
  （两键必须成对），和 = N 必须等于发作业核数，否则求解器启动时打印诊断并**终止作业**
  （退出码 6 只写一个键 / 7 和不符 / 8 与模式冲突）；`./run.sh`（不带参数）会把这个和当作
  默认 np。不写这两个键则保持旧的隐式划分、不做检查。详见 `cluster_cases/README.md` 第 3 节。
* `variants/pressure_outlet.control` 用法：
  `cp variants/pressure_outlet.control unMesh.control && MIXROOT=… ./run.sh`。
