# porous-graetz：二维多孔介质强制对流（Graetz 热入口）验证

阶段 9 验证算例之一（需求 docs/程序功能说明.md L15）。Darcy 活塞流通过恒壁温
多孔填充平行板通道，热入口段温度场对拍 **Graetz 级数解析解**；同时用两个工况
（纯分子导热 / 叠加横向热弥散）检验热弥散模型在对流换热中的作用。

## 1. 构型

半平行板通道（中心线对称 + 单侧恒壁温）：

- y∈[0,a]，a=H/2=5 mm；x∈[0,L]，L=0.45 m；网格 180×40×1（wedge/tet）。
- y=0 中心线 symmetry（绝热）；y=a 壁面 wall，`tbc = 8 1 310.0`（Tw=310 K）；
  x=0 velocity-inlet（T0=300 K）；x=L pressure-outlet(0)；z 两面 symmetry。
- 多孔：ε=0.4，K=1e-10 m²（VC:porous 标记判型，control 只给系数）。
- 空气 ρ=1.177、μ=1.846e-5、cp=1005、k_cond=0.026；
  **k_s=k_f=0.026**，使 LTE 有效导热 k_eff=ε·k_f+(1−ε)·k_s=k_f，
  从而解析 κ_eff=k_eff/(ρcp) 干净无孔隙率因子。

## 2. 解析解（活塞流 Graetz）

θ=(T−Tw)/(T0−Tw)，η=y/a，x*=κ_eff·x/(U·a²)，λ_n=(n−1/2)π：

```
θ(x,η)  = Σ_n 2(−1)^(n+1)/λ_n · cos(λ_n η) · exp(−λ_n² x*)
θ_b(x)  = Σ_n 2/λ_n² · exp(−λ_n² x*)        （体均/混合杯温度）
Nu_a(x) = [−∂θ/∂η]_wall / θ_b
```

充分发展渐近值 Nu_a=λ_1²=π²/4=**2.4674**
（特征长度为半高 a；按水力直径 D_h=4a 为 9.87，按全高 H 为 4.93——
Nu 定义差异仅为常数倍，本文统一用 Nu_a）。

两工况：

| 工况 | U (m/s) | α_t (m) | κ_eff (m²/s) | Pe | 出口 x* |
|---|---|---|---|---|---|
| graetz_mol | 0.1 | 0 | 2.198e-5 | 22.7 | 3.96 |
| graetz_disp | 1.0 | 2.5e-4 | 2.720e-4 | 18.4 | 4.90 |

κ_eff=(k_eff+ρ·cp·α_t·U)/(ρcp)（Bear 横向弥散，α_l=0 不加流向弥散）。
两工况 Pe≈20，轴向导热影响限于入口约 1 列。

## 3. 复现命令

```bash
python3 gen_graetz.py graetz.cas
../../bin/uns_solver graetz.cas graetz_mol.control  > graetz_mol.log 2>&1;  mv graetz.vtu graetz_mol.vtu
../../bin/uns_solver graetz.cas graetz_disp.control > graetz_disp.log 2>&1; mv graetz.vtu graetz_disp.vtu
python3 plot_graetz.py
# MPI（np2 与串行同为 6238 步，du_max 位级一致）:
mpirun -np 2 ../../bin/uns_solver_mpi graetz.cas graetz_disp.control
```

## 4. 结果（串行，outer_tol=1e-9，5601/6238 步收敛）

剖面最大偏差 |θ_num−θ_series|、体均温度、局部 Nu_a：

| 截面 x* | mol 剖面 max | mol θ_b 误差 | disp 剖面 max | disp θ_b 误差 |
|---|---|---|---|---|
| 0.05/0.20 | 0.026 | 1.1% | 0.011 | 1.2% |
| 0.5 | 0.006 | 1.6% | 0.013 | 3.4% |
| 1.0 | 0.006 | 5.5%* | 0.008 | 7.7%* |
| 2.0 | 0.001 | 14%* | 0.002 | 17%* |

\* x*≥1 后 θ_b 本身已 <0.07，相对误差被小分母放大；绝对偏差始终 ≤7e-3。

- 充分发展 Nu_a：mol **2.457**（vs 2.4674，0.42%）；
  disp **2.434**（1.38%）。
- 两工况在各自 x* 坐标下体均温度曲线坍并到同一条级数曲线（见右图），
  定量证明 Bear 横向弥散在对流换热中恰好以 κ_eff=k_eff/(ρcp)+α_t·U
  增强横向输运。
- 图：`images/graetz.png`（两工况剖面 + 体均坍并）。

## 5. 建模细节与教训

- **近壁无弥散层**：Bear 弥散张量依赖格心速度，壁面无滑移 Brinkman 层
  （厚度 ~√K）使首排格速度亏缺 → 近壁弥散自动衰减，物理上真实存在，
  但解析 Graetz 假设弥散均匀到壁。首轮 K=1e-9 时 √K/dy=0.25、首排
  u=0.895U，渐近 Nu 偏高 16%；改 K=1e-10（√K/dy=0.08、u1=0.988U）
  后偏差降到 1.38%。要继续降低需更小 K（Δp∝1/K）或更细近壁网格。
- k_s 必须显式设为 k_f，否则 LTE 的 k_eff=ε·k_f 会改变解析 κ。
- VTU 每格输出 12 个四面体，分析脚本已剥离 connectivity 首元素 4。

## 6. 文件清单

- `gen_graetz.py`：半通道网格（zone3 内部/4 z-sym/5 中心线/6 入口/7 出口/8 热壁）
- `graetz.cas`、`graetz_mol.control`、`graetz_disp.control`
- `graetz_mol.log/.vtu`、`graetz_disp.log/.vtu`
- `plot_graetz.py`、`images/graetz.png`
