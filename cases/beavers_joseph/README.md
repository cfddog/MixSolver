# beavers-joseph：流体/多孔界面应力跳变（BJ）条件验证

阶段 9 收尾算例。开口槽道中上半为自由流体、下半为多孔床（ε=1，Brinkman 项
激活），界面用 `bj_alpha` 施加 Beavers-Joseph 型应力跳变通量；与**双层耦合
ODE 精确参考解**对拍（注意：不是经典 BJ 滑移公式，见 §2）。

## 1. 构型

- 二维开口槽道：x∈[0,0.1] m，y∈[0,0.01] m，网格 100×40×1（六面体）。
- 上半 y∈[5,10] mm 为流体（zone 2），下半 y∈[0,5] mm 为多孔（zone 3，
  `VC:porous` 标记），H=5 mm 每层。
- x=0 / x=L 均 pressure-outlet(0)，驱动力为全局体力
  `body_force = 0.05 0 0`（N/m³）。上下壁面 wall，z 两面 symmetry。
- 多孔参数：K=1e-6 m²、ε=1（λ=√K=1 mm，dy=0.25 mm → Brinkman 层 4 单元
  可分辨）；空气 ρ=1.177、μ=1.846e-5。
- 三工况：`bj_a0`（无 bj_alpha，标准内部面=速度+应力连续）、`bj_a1`
  （α=1）、`bj_a2`（α=2）。

## 2. 参考解（双层应力跳变模型，非经典 BJ）

离散 BJ 通量的连续极限是"应力连续 + 速度跳变"的串联阻力界面模型
（1/C = d_Pf/(μA) + λ/(μαA)）。由于多孔侧 Brinkman 项激活存在可分辨的
Brinkman 层，**经典 BJ 滑移公式（假设多孔侧纯 Darcy）不是正确参考**。

精确参考（y 自界面起算，流体 y∈[0,H]，多孔 y∈[−H,0]）：

```
fluid : μ u_f'' + fb = 0,                u_f(H)  = 0
porous: μ u_p'' − (μ/K) u_p + fb = 0,    u_p(−H) = 0
interface: μ u_f'(0) = μ u_p'(0) = τ
  α>0 : τ = μ (α/λ) (u_f(0) − u_p(0))    （应力跳变）
  α=0 : u_f(0) = u_p(0)                   （速度连续）
```

解形式 u_f = −fb y²/(2μ) + a1·y + a2，
u_p = u_D + E·e^{y/λ} + F·e^{−y/λ}（u_D=K·fb/μ，λ=√K），
系数由 4×4 线性系统解出（plot_bj.py `reference()`）。

## 3. 复现命令

```bash
python3 gen_bj.py bj.cas
../../bin/uns_solver bj.cas bj_a0.control > bj_a0.log 2>&1; mv bj.vtu bj_a0.vtu
../../bin/uns_solver bj.cas bj_a1.control > bj_a1.log 2>&1; mv bj.vtu bj_a1.vtu
../../bin/uns_solver bj.cas bj_a2.control > bj_a2.log 2>&1; mv bj.vtu bj_a2.vtu
python3 plot_bj.py
# MPI（与串行最大相对偏差 2.2e-7）:
mpirun -np 2 ../../bin/uns_solver_mpi bj.cas bj_a1.control
```

## 4. 结果（串行，outer_max=800）

x=L/2 剖面对拍（plot_bj.py 输出）：

| α | 首流体单元 u 数值 | 参考 u_f(dy/2) | 全剖面 RMS |
|---|---|---|---|
| 0 | 8.539e-3 | 8.498e-3 | 0.14% |
| 1 | 1.173e-2 | 1.212e-2 | 1.09% |
| 2 | 1.002e-2 | 1.045e-2 | 1.29% |

- α=0 精确复现连续界面（RMS 0.14%，纯离散误差量级）。
- α=1/2 复现界面速度跳变 + Brinkman 层，全剖面 RMS ~1%（λ/dy=4 分辨
  率下的预期截断误差）。
- 经典 BJ 滑移速度 u_B（α=1: 7.9e-3、α=2: 5.5e-3 m/s）与数值明显不符，
  反证正确参考必须是双层耦合系统；图中以竖点线标出对比。
- 图：`images/beavers_joseph.png`（左：全剖面；右：界面放大，实线为
  双层参考解）。

## 5. 建模细节与教训

- `bj_alpha` 挂在多孔 `cell_zone` 行上（如 `cell_zone = 3 perm=1.0e-6
  porosity=1.0 bj_alpha=1.0`），默认 0 = 禁用（位级回归不变）。界面通量
  等值反向作用于两侧单元，动量守恒。
- 闭盒体力会被压力梯度完全抵消，BJ 验证必须用**开口槽道**（两端
  pressure-outlet）。
- 温度在此算例中为纯被动标量（boussinesq=.false.），不影响速度场；
  当前 passive T 方程在长迭代下漂移使 dT_max 门槛不收敛，故取
  outer_max=800 直接截取已收敛速度场（du_max~1e-8）。
- VTU 输出名跟随 mesh 文件名（bj.cas → bj.vtu），并行/串行连续跑多工况
  需及时 `mv`。

## 6. 文件清单

- `gen_bj.py`：双层槽道网格（zone 2 fluid 上半 / 3 VC:porous 下半 /
  4 interior / 5 z-sym / 6 底壁 / 7 顶壁 / 8,9 pressure-outlet）
- `bj.cas`、`bj_a0/a1/a2.control`、`bj_a0/a1/a2.log`、`bj_a0/a1/a2.vtu`
- `bj_mpi.log`、`bj_mpi.vtu`（np=2 验证）
- `plot_bj.py`、`images/beavers_joseph.png`
