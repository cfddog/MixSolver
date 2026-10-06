# porous-disp：多孔介质横向热弥散验证（2D 温度混合层）

阶段 9 热弥散模型的验证算例。均匀 Darcy 流通过 `VC:porous` 多孔块，入口在
y=H/2 处分成冷热两股（310 K / 300 K），横向热弥散（Bear 分量式张量的 α_T 项）
使混合层展宽，对拍误差函数相似解。

## 1. 物理模型（本轮实现）

面扩散系数在静态有效导热之外叠加**弥散张量的法向投影** n·D·n（对角、轴对齐）：

```
D_dd = ρcp·( α_L·u_d² + α_T·(|u|²−u_d²)/|u| )   [W/m/K]
```

- LTE：加到混合有效导热 k_eff；LTNE：加到流体相导热 kf（固相不加）。
- α_L / α_T（纵向/横向弥散度，单位 m）按 cell_zone 逐区给定：`disp_l=.. disp_t=..`。
- 静止格（|u|≈0）与未给弥散度的格自动为 0。
- 非正交修正复用含弥散的 kf（正交网格下精确）。

同轮实现的**各向异性渗透率**（轴对齐对角张量）见 cases/porous_plug/README.md。

## 2. 算例设置

- 网格：`gen_mix.py` → 80×80×1 楔形单元，L=0.04 m，H=0.06 m（`mix.cas`）。
- 上下壁 symmetry，x=L pressure-outlet(0)；x=0 两个 velocity-inlet 区：
  下半 T=310 K、上半 T=300 K（zone 6 / zone 8）。
- 空气：ρ=1.177，μ=1.846e-5，cp=1005，k_cond=0.026。
- 多孔：ε=0.4，K=1e-8 m²（Darcy 强制活塞流，Re_p≈6.4），u=1 m/s。
- 控制文件（系数专用行，类型取自 VC 标记）：
  `cell_zone = 2 perm=1.0e-8 porosity=0.4 disp_t=1.0e-3`（disp_l=0，
  使流向弥散不扰动相似解形式）。
- 解析解：κ = Γy/(ρcp) = (k_cond + ρ·cp·α_t·u)/(ρcp) = 1.0220e-3 m²/s，

```
T(x,y) = 305 − 5·erf( (y−H/2)·sqrt( u/(4·κ·x) ) )
```

注意 erf 相似变量里是**热扩散率 κ**（m²/s），不是 Γy（W/m/K）。

## 3. 复现命令

```bash
python3 gen_mix.py mix.cas
../../bin/uns_solver mix.cas mix_disp.control  > mix_disp.log 2>&1;  mv mix.vtu mix_disp.vtu
../../bin/uns_solver mix.cas mix_nodisp.control > mix_nodisp.log 2>&1; mv mix.vtu mix_nodisp.vtu
python3 plot_mix.py
# MPI（结果与串行一致，188 步）:
mpirun -np 2 ../../bin/uns_solver_mpi mix.cas mix_disp.control
```

## 4. 结果

收敛：disp 工况 188 步（du_max=7e-15），nodisp 85 步；np=2 与串行同为 188 步。
u_x=0.9966（入口列恢复段拉低均值，内部列=1.0000），|u_y|max≈2e-4。

计算剖面 vs 相似解（偏差相对 5 K 全步长）：

| 截面 | x (mm) | max\|dT\| (K) | 偏差 | rms (K) |
|---|---|---|---|---|
| col16 | 8.2 | 0.052 | 1.03% | 0.018 |
| col40 | 20.2 | 0.030 | 0.60% | 0.013 |
| col60 | 30.2 | 0.025 | 0.50% | 0.013 |
| col76 | 38.2 | 0.023 | 0.46% | 0.012 |

- 弥散工况与 erf 解析解定量吻合（≤1%，首截面略大是台阶入口的离散化）；
- `disp_t=0` 对照：混合层宽 ~1.3 mm（约 2 个网格），欠分辨不可量化，
  仅作层宽对比（弥散使 Γy 增大 46 倍、层宽 ~9 mm）；
- 图：`images/mix_layer.png`。

## 5. 结论

Bear 分量式横向热弥散模型实现正确；α_T 驱动的混合层展宽与封闭解吻合。
（纵向 α_L 经同一张量机制生效，本算例设 0 以保持相似解形式。）

## 6. 文件清单

- `gen_mix.py`：网格生成器（VC:porous 标记 + 冷热双入口分区）
- `mix.cas`：网格
- `mix_disp.control` / `mix_nodisp.control`：带/不带弥散工况
- `mix_disp.log` / `mix_nodisp.log`、`mix_disp.vtu` / `mix_nodisp.vtu`
- `plot_mix.py`、`images/mix_layer.png`
