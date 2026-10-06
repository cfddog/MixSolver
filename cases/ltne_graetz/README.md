# LTNE Graetz 强制对流 + 横向热弥散组合验证

> 2026-10-05 | 阶段 9 收尾第二项：2D 局部热非平衡（LTNE）双温度 Graetz 热入口
> 问题，含/不含横向热弥散 α_t 两个工况，对拍**双相耦合特征值半解析参考**
> （numpy 离散 η 向二次特征值问题，一阶块系统）。参考解保留两相轴向导热，
> 入口取 T_f=T_in、固相入口绝热（与求解器一致）。

## 1. 算例设置

复用 `cases/porous_graetz` 的半平行板网格 `ltne_graetz.cas`
（由 gen 脚本复制：180×40 hex，半高 a=5 mm，L=0.45 m，180×40）：
zone 2=`VC:porous`、3 interior、4 z-sym、5 中心线 symmetry、
6 velocity-inlet（x=0）、7 pressure-outlet（x=L）、8 恒壁温 wall（y=a）。

| 项 | 值 |
|---|---|
| 流体（空气） | ρ=1.177，μ=1.846e-5，cp=1005，k_f=0.026，U=1 m/s，T_in=300 K |
| 多孔 | ε=0.4，K=1e-10 m²，陶瓷骨架 k_s=1.0/cp_s=800/ρ_s=2000 |
| 相间交换 | h_sf=2.0e5，a_sf=1（Bi=H a/(ρcp U)=0.845，交换长度 √(Ds/H)=1.73 mm≈14 dy） |
| 固相横向导热 | Ds=(1−ε)k_s=0.6 W/(m·K) |
| 壁面 | zone 8 恒壁温 Tw=310 K（两相均 Dirichlet） |
| 数值 | outer_max=20000，outer_tol=1e-9，lin_tol=1e-10，inlet_ramp=1 |

两工况（`ltne_grz_d.control` / `ltne_grz0.control`）：

| 工况 | α_t (m) | K_fy=εk_f+ρcp·α_t·U | κ_fy (m²/s) | Pe=Ua/κ |
|---|---|---|---|---|
| disp | 2.5e-4 | 0.3061 | 2.59e-4 | 19.3 |
| nodisp | 0 | 0.0104 | 8.79e-6 | 568.7 |

固相无量纲轴向导热 Λs=Ds/(ρcp U a)=0.1014。选 k_s=1.0（≈k_f 的 38 倍）
使固相与流体在全流场保持可观温差（真正的非平衡态），同时壁面热量经
高导骨架快速前传。

## 2. 半解析参考（`plot_ltne_graetz.py`）

θ=(T−Tw)/(T0−Tw)，η=y/a，ξ=x/a，稳态活塞流（表观 U）：

```
(1/Pex) θf,ξξ + (1/Pe) Lη θf = θf,ξ + Bi(θf−θs)
Λs (θs,ξξ + Lη θs)           = Bi(θs−θf)
```

横向 1/Pe=K_fy/(ρcp U a)；轴向 1/Pex=K_fx/(ρcp U a)，
K_fx=εk_f+ρcp·α_l·U（α_l=0 → 0.0104，横/轴向系数不同，必须分开）。
Lη=d²/dη²，η=0 绝热、η=1 两相均 θ=0。

模态 e^{μξ} 构成二次特征值问题，化为一阶块系统
y'=My，y=[θf,θs,χf,χs]（χ=∂/∂ξ），M 为 4m×4m；
取 Re μ<0 的 2m 个衰减模态（λ=−μ），入口两行条件
θf(0)=1、χs(0)=0（固相入口绝热，正是求解器处理）解模态系数。
η 向 m=180 节点有限差分（中心线 (−6,7,−1)/3 二阶单边、壁面 Dirichlet）。
特征值全实、入口残差 ~1e-13；Bi→0 极限 λ1·Pe=2.457≈π²/4=2.467（0.4%
FD 截断），退化为经典 Graetz。

## 3. 结果

两工况串行均 dT_max=0 正式收敛：disp **7018** 步、nodisp **6939** 步。

### 3.1 剖面对拍（x=20…350 mm 五站位，θ 满量程=1）

| 工况 | 网格 | max \|θf−ref\| | max \|θs−ref\| |
|---|---|---|---|
| disp α_t=2.5e-4 | 180×40 | **0.014** | 0.011 |
| disp | 360×40（x 加密） | **0.0073** | — |
| nodisp | 180×40 | **0.033** | 0.026 |
| nodisp | 360×40（x 加密） | **0.018** | 0.014 |

- 残差全部集中在**入口段 x=20 mm**（ξ≈4）：x≥50 mm 后 ≤0.005、
  x≥100 mm 后 ≤0.002、x≥200 mm 后 ≤1e-4。
- x 向加密（dx/a 0.5→0.25）入口误差近似减半，证实残差是入口段
  x 向离散误差（固相轴向导热 Λs=0.10 与相间弛豫长度 ~1.2a 在粗网格上
  分辨率有限），随加密单调趋向参考解。
- 流相体均温度 θ_f,b 沿 x* 与 ξ 两坐标与模态展开全程重合
  （图右上/右下：圆点=数值，线=参考）。
- 入口第一列固相剖面 vs χs=0 模态解偏差 0.05–0.06（同一入口离散效应）。

### 3.2 结果图

`images/ltne_graetz.png`：左列两工况 η 剖面（实心圆=流体、方框=固体、
空心三角=360×40 x 加密入口站位，实线=流体参考、虚线=固体参考）；
右上 θ_b(x*)，右下 θ_b(ξ)。

### 3.3 串行/MPI(np2) 对拍

`mpirun -np 2`（disp 工况，7016 步收敛 vs 串行 7018）：
T_f、T_s 全场最大绝对偏差均为 **2.0e-6 K（相对满量程 2e-7）**，
MPI VTU 含 temperature_solid 字段。

## 4. 复现

```bash
cp ../porous_graetz/graetz.cas ltne_graetz.cas
../../bin/uns_solver ltne_graetz.cas ltne_grz_d.control
mv ltne_graetz.vtu ltne_grz_d.vtu
../../bin/uns_solver ltne_graetz.cas ltne_grz0.control
mv ltne_graetz.vtu ltne_grz0.vtu
mpirun -np 2 ../../bin/uns_solver_mpi ltne_graetz.cas ltne_grz_d.control
mv ltne_graetz.vtu ltne_grz_d_mpi.vtu
# x-direction refinement (360x40):
python3 gen_graetz_xref.py ltne_graetz_xref.cas
../../bin/uns_solver ltne_graetz_xref.cas ltne_grz0_xref.control
mv ltne_graetz_xref.vtu ltne_grz0_xref.vtu
../../bin/uns_solver ltne_graetz_xref.cas ltne_grz_d_xref.control
mv ltne_graetz_xref.vtu ltne_grz_d_xref.vtu
python3 plot_ltne_graetz.py
```

## 5. 结论

LTNE 双温度 × 横向热弥散 2D 组合路径验证通过：5 站位两相剖面在
180×40 网格上最大偏差 1.4%（disp）/3.3%（nodisp，均仅入口段），
x 加密后减半至 0.7%/1.8%，体均温度与模态展开全程重合，
串并行两相场一致 2e-7。建模要点：弥散只进流体相横向有效导热、
两相各保留自己的轴向导热（α_l=0 时流体轴向仅 εk_f，横/纵 Pe 必须分开）、
恒壁温对两相 Dirichlet、固相入口绝热。
