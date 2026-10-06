# 不可压缩槽道流验证：mass-flow-inlet 启动稳定性 + outflow 出口边界

> 日期：2026-10-05　求解器：`bin/uns_solver` / `bin/uns_solver_mpi`（非结构 SIMPLE）
> 主题：① 验证 mass-flow-inlet 在纯流体槽道（无 Darcy 阻力）的启动稳定性；
> ② 诊断 pressure-outlet 出口段伪调整；③ 新增 `outflow`（充分发展出口）边界条件并验证。

---

## 1. 算例设置

- **几何**：二维平面槽道，x∈[0,300]、y∈[0,10]、z∈[0,10]（单层挤出，30000 hex）。
  网格 wall-clustered：dy_wall=0.01，dy_c=0.233，dx=1（dx/dy 高达 100:1）。
- **网格文件**：`channel.cas`（Pointwise 导出原件，/mnt/c/temp/channel.cas）；
  **`channel_wall.cas`**（实际使用）：Pointwise 把 z=0/z=10 对称面与 y=0/y=10 壁面
  合并导出一个 zone 4，用 `split_channel_zones.py` 按面号区间拆分为
  zone 4 symmetry + zone 7 wall（面号不变、cell 区段不动）。
- **边界**：zone 5 = mass-flow-inlet（x=0），zone 6 = pressure-outlet 或 outflow（x=300），
  zone 7 = wall，zone 4 = symmetry。
- **流体**：300 K 空气，rho=1.177 kg/m³，mu=1.846e-5 Pa·s（无能量方程）。
- **工况**：
  | 算例 | mdot [kg/(m²·s)] | U_bulk [m/s] | Re_H | 说明 |
  |---|---|---|---|---|
  | Re100  | 1.8460e-4 | 1.5684e-4 | 100  | 基层流验证 |
  | Re1000 | 1.8460e-3 | 1.5684e-3 | 1000 | 10× 启动冲量压力测试 |
- **数值**：SIMPLE 稳态，alpha_u=0.7、alpha_p=0.3，一阶迎风（conv_blend=0），
  nonorth_corr=1，ppe_precond=ic0，lin_tol=1e-8。
  **注意**：`outer_tol` 是绝对 du 阈值；U~1.6e-4 时 1e-6 = 0.64% U 会假收敛
  （首跑 43 步"收敛"但出口剖面 1.93U、质量漂移 5.4%）。本算例用 1e-10（~6e-7·U）。

## 2. 解析参考（plane Poiseuille，H=全高）

- u(y) = 6·U·η(1−η)，u_max/U = 1.5；dp/dx = −12μU/H²；
- 入口段估计 L_e ≈ 0.05·Re_H·H：Re100 → 5H，Re1000 → 50H（> L=30H，出口仍在发展）。

## 3. 结果

### 3.1 mass-flow-inlet 启动稳定性（主题①）

- Re100 从零场启动：**it=1 du_max = 2.84e-4 ≈ 1.8·U，mass-imbal ~5e-13（机器零）**，
  无发散、单调收敛（790–807 步）。
- Re1000（10× 冲量）：outflow 1033 步收敛，同样稳定。
- 结论：**上会话"纯流体 mdot-inlet 启动发散"确认为 LTNE 全无壁面网格的特异拓扑
  （发散源自出口零阻力+销钉 cell 1 假质量汇，已在 LTNE 任务修复），非 mdot-inlet 本身问题。**
  槽道算例 mdot-inlet 启动稳健。

### 3.2 内部流场精度（两出口条件一致到 7 位有效数字）

| 指标 | 数值 | 参考 | 偏差 |
|---|---|---|---|
| U(50.5)/U_bulk（Re100） | 1.567715e-4 | 1.568394e-4 | −0.04% |
| dp/dx（Re100, x∈[200,295] 拟合） | −3.4734e-10 | −3.4743e-10 | −0.02~−0.03% |
| 入口段 L_e | 65.5 = 6.55H | ~5H | 合理（0.5% 判据+网格） |
| x≤15H 剖面 | 落在解析抛物线上 | — | — |
| 串行 vs MPI(np2) | 同为 ~790 步收敛，全场速度差 ~5e-7 | — | 分区浮点噪声 |

### 3.3 出口段：pressure-outlet vs 新增 outflow（主题②③）

Re100 出口末列（x=299.5）对比：

| 指标 | pressure-outlet | **outflow** | 解析 |
|---|---|---|---|
| u_max/U | 1.7091 | **1.6420** | 1.5 |
| 横向速度 vmax/U | 0.1152 | **0.0402** | 0 |
| 剖面 L1 误差 | 11.45% | **5.63%** | 0 |
| 剖面最大误差 | 20.99% | **14.28%** | 0 |
| 末列体速度偏差 | −0.05% | +3.5%（=1/β，β=0.966） | 0 |
| 收敛步数 | 807 | 790（np2: 791，出口值一致） | — |

- pressure-outlet 异常（假"出口调整"：Ucl 1.5→1.71、横向速度向中心汇聚）经对照实验
  逐一排除入口 BC（velocity-inlet 同异常）、对流格式（conv_blend=1 同）、面法向、
  p′ 修正路径、并行因素，定性为**同位网格 + Rhie-Chow 在固定压力边界的离散不动点特性**。
- 新 `outflow` 边界：法向零梯度（u/p/T 外推）+ **全局质量缩放**（每次 Rhie-Chow 后
  把 outflow 面通量整体缩放使总出流=总入流），PPE 为纯 Neumann（触发 cell 1 销钉，
  缩放保证 Σrhs=0 相容，销钉不再成为假质量汇）。出口异常削弱约 3 倍。
- **已知限制**：末列仍有 ~9% 中心线凸起与 +3.5% 体速度外观偏差（守恒的面通量精确
  等于入流；偏差是末列 cell 中心速度与缩放通量的 1/β 不动点特性）。
  工程惯例规避：出口远离关注区 ~10H，不取末列数据。
  （尝试过动量方程改用滞后缩放通量：β 不收拢且略差，已回退。）

### 3.4 Re1000 对照（发展段流，Le≈50H>L=30H 全程发展中）

| 指标 | pressure-outlet | outflow |
|---|---|---|
| 收敛步数 | 1027 | 1033 |
| U(50.5)/U | 0.9997 | 0.9997 |
| u_max/U @exit | 1.691（末列额外抬高=伪调整） | 1.574（仍在发展，非异常） |
| 末列剖面 L1 误差 | 14.74% | 5.36% |

- 内部场两出口一致（U(50.5) 逐 7 位相同）；pressure-outlet 末列抬高机制与 Re100 相同。
- dp/dx 拟合偏差 ~28%（两者同）：出口段仍在发展，非误差。
- 注意脚本"mass drift in->out"比较的是入口第 0 列（近壁已减速、体平均偏低）
  与末列，不是质量不守恒指标。

## 4. 回归（确认 outflow 新增与 correct_fields p′ 修正无破坏）

| 算例 | 结果 |
|---|---|
| cavity Re100 | 518 步收敛；VTU 仅 z 分量 1e-20 机器噪声差异；Ghia 对比一致 |
| cylinder Re40 far-field | Cp 与基线平均差 0.0031（≈已记录残差平台噪声 0.004），Cp 范围一致 |
| porous_Ra10 | 748 步、Nu=1.07800 与归档完全一致；打印迭代历程逐字相同；末场 1e-4 相对漂移=非线性浮点放大（codegen 噪声，同 cavity np2 先例） |
| LTNE 1D | 3344 步（原 3343，差 1 步源自 p′ 修正路径）；ΔT_f(L)=995.02 K（−0.001%）等全部指标与归档一致 |

## 5. 复现命令

```bash
cd cases/channel
# 网格拆分（已执行，产物 channel_wall.cas）
python3 split_channel_zones.py
# Re100：pressure-outlet 对照 / outflow
../../bin/uns_solver channel_wall.cas channel_Re100.control     channel_Re100.vtu     > channel_Re100.log
../../bin/uns_solver channel_wall.cas channel_Re100_out.control channel_Re100_out.vtu > channel_Re100_out.log
# MPI 一致性（输出文件名由网格名派生，跑完改名）
mpirun -np 2 ../../bin/uns_solver_mpi channel_wall.cas channel_Re100_out.control > channel_Re100_out_np2.log
# Re1000 对照（pressure-outlet / outflow）
../../bin/uns_solver channel_wall.cas channel_Re1000.control     channel_Re1000.vtu     > channel_Re1000.log
../../bin/uns_solver channel_wall.cas channel_Re1000_out.control channel_Re1000_out.vtu > channel_Re1000_out.log
# 出图（第 3 参数为 Re_H，默认 100）
python3 plot_channel.py channel_Re100_out.vtu   images/channel_Re100_out.png   100
python3 plot_channel.py channel_Re1000.vtu      images/channel_Re1000.png      1000
python3 plot_channel.py channel_Re1000_out.vtu  images/channel_Re1000_out.png  1000
```

诊断对照算例（出口异常归因，保留备查）：`channel_Re100_vinlet.control`
（velocity-inlet 对照，与 mdot 版逐位相同）、`channel_Re100_ho.control`
（conv_blend=1 二阶对照，异常同）、`channel_Re100_rlx.control`
（αu=0.9/αp=0.5 强松弛，振荡未收敛，手动终止）。

## 6. 图片

- `images/channel_Re100.png` — pressure-outlet 基准（剖面/中心线发展/压力）
- `images/channel_Re100_out.png` — outflow 出口（本任务主结果）
- `images/channel_Re1000.png` — Re1000 发展段流 pressure-outlet 对照
- `images/channel_Re1000_out.png` — Re1000 发展段流 outflow

## 7. 结论

1. **mass-flow-inlet 启动稳定性验证通过**：纯流体槽道从零场无发散单调收敛，
   10× 冲量（Re1000）同样稳定；上会话发散为 LTNE 网格特异问题（已修）。
2. 内部流场精度优秀：充分发展段剖面贴解析解、dp/dx 偏差 0.02–0.03%、
   入口段长度合理、串行/MPI 一致。
3. 新增 `outflow`（充分发展出口）边界条件实现正确：纯 Neumann PPE + 全局质量
   缩放自洽，出口段伪调整较 pressure-outlet 削弱约 3 倍；末列 ~9% 凸起为
   一阶零梯度边界固有特性，记录为已知限制。
4. 全套回归无破坏。
