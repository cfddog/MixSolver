# LTNE + 纵向热弥散 1D 发汗冷却组合验证

> 2026-10-05 | 阶段 9 收尾项：LTNE 双温度模型与 Bear 热弥散张量的**组合**验证。
> 1D 恒定壁面热流构型下解析解封闭，纵向弥散系数 α_l 直接进入流体有效导热，
> 可精确对拍 4 阶耦合 ODE 的半解析解。同算例完成串行/MPI(np2) 双温度对拍，
> 并在调试中修复两个真实 bug（纵向弥散公式、MPI 缺收 T_s）。

## 1. 算例设置

**物理场景**：多孔薄板 x∈[0,L]（L=0.05 m），空气以表观质量流速 G=2 kg/(m²·s)
正向通过，出口端面 x=L 施加恒定热流 q″=2.0e6 W/m²。热量沿固相骨架上溯、
经 h·a 交换给流体并被带出（transpiration cooling）。弥散只加流体相。

**网格**：`gen_ltne1d.py` 生成 `ltne_disp.cas`（NX×NY=100×2 hex 薄块，SI 米）：
606 节点 / 200 hex（VTU 输出为每 hex 12 tetra）/ 902 面；
zone 2=`VC:porous`、3 interior、4 z-sym、5 y-sym、
6 mass-flow-inlet（x=0）、7 pressure-outlet（x=L，兼热流面）。

**边界与参数**（两工况共用，差别仅 disp_l）：

| 项 | 值 |
|---|---|
| 流体（空气） | ρ=1.177，μ=1.846e-5，cp=1005，k_f=0.026 |
| 多孔 | ε=0.3，K=1e-8 m²，SS304 k_s=16.2/cp_s=500/ρ_s=7900 |
| 相间交换 | h_sf=2.0e5 W/(m³·K)，a_sf=1（交换长度 √(Ds/H)=7.54 mm≈15 dx） |
| 入口 | mass-flow-inlet G=2，T_in=300 K（u=G/ρ=1.699 m/s） |
| 出口 | p=0；`tbc = 7 2 2.0e6` 恒热流 |
| 数值 | outer_max=20000，outer_tol=1e-6，conv_blend=0，inlet_ramp=1 |

- `ltne_disp0.control`：无弥散，Kf=εk_f=0.0078 W/(m·K)
- `ltne_disp1.control`：disp_l=0.005 m，Kf=εk_f+ρcp·α_l·u=10.058 W/(m·K)

## 2. 半解析参考（`plot_ltne_disp.py`）

体积平均 1D 双温度方程（弥散作为流体轴向有效导热）：

```
G cp T_f' = (Kf T_f')' + H (T_s - T_f)
0 = Ds T_s'' + H (T_f - T_s)              Ds=(1-ε)k_s=11.34, H=h·a
```

设 T_f=c0 + Σ c_j exp(r_j x)，T_s=Σ m_j c_j exp(r_j x)，
m_j=1+(a r_j−Kf r_j²)/H，a=G cp；r 为三次多项式
`[Ds Kf, −Ds a, −H(Ds+Kf), H a]` 的三个根（一零根、一正一负）。

边界条件：T_f(0)=T_in、Ds T_s'(0)=0、Kf T_f'(L)=q_f、Ds T_s'(L)=q_s，
热流按微观导热分数 Nield 分配 q_f=q·Kf/(Kf+Ds)、q_s=q·Ds/(Kf+Ds)
（弥散增大 Kf → 流体承担份额在 disp1 工况显著上升）。
未知量取 [c0, d1,d2,d3]，d_j=c_j e^{r_j L}，指数写 exp(r(x−L))，
避免 disp_l=0 时刚性正根（r≈2.6e5）exp 上溢。

**总量守恒**：T_f(L)−T_in=q/(G cp)=995.0 K（与弥散无关；
弥散只改变热量回传分布与出口两相分配）。

## 3. 结果

两工况均正式收敛（dT_max=0）：disp0 **8095** 步、disp1 **11219** 步。

| 指标 | disp0 | disp1（α_l=5 mm） |
|---|---|---|
| T_f 剖面最大偏差（占 995 K） | 2.31%（出口陡段） | 0.78% |
| T_f 剖面中位偏差 | 0.42% | 0.41% |
| T_s 剖面最大偏差 | 0.83% | 0.30% |
| ΔT_f(L)，disp0 | 994.99 vs 995.01 | — |
| ΔT_f(L)，disp1 | 973.46 vs 974.41 | 回流传热（见下） |

\* 全局流体平衡 a·ΔT_f(L) = q − K_f T_f'(0)：出口热流按微观导热分数
Nield 分配（q_f=q·εk_f/(εk_f+Ds)、q_s 其余，**分配不含弥散贡献**，
与求解器实现一致）；disp1 工况弥散使流体有效轴向导热 Kf 很大，
入口又是 T_f=T_in 的 Dirichlet 面，部分热量被回流到入口端排出，
故出口流体温升小于 q/(Gcp)=995 K（参考解与求解器一致，差 0.1%）。
disp0 时 Kf 极小、回流通量可忽略，温升即 995 K。

图：`images/ltne_disp.png`（三联：T_f、T_s、T_s−T_f；点=求解器，线=半解析）。

### 3.1 串/并行对拍（np=2）

`mpirun -np 2 ../../bin/uns_solver_mpi ltne_disp.cas ltne_disp1.control`
11220 步收敛（串行 11219，末位判据噪声），gather 后全场：
**T_f 相对偏差 9.7e-8、T_s 6.4e-7 量级（绝对 ~1e-4 K）**，
与既有 B-J np2 一致性（2.2e-7）同量级。

## 4. 本算例修复的两个 bug

1. **纵向弥散缺 /|u|**（`src/unstructured/mod_uns_simple.f90`，温度装配）：
   Bear 张量纵向分量原为 ρcp·α_l·u_d²，量纲 m²/s²（多乘一个 |u|），
   横向分量本有 /umag。修正为
   `kdisp = ρcp/umag·(α_l u_d² + α_t(|u|²−u_d²))`。
   回归：纯横向弥散算例（porous_disp 混合层两工况、porous_graetz 两工况）
   修复前后 VTU **位级一致**（横向分量数学形式未变）。
2. **MPI gather 缺 T_s**（`mod_uns_gather.f90` + `main_uns_mpi.f90`）：
   `gather_fields_to_root` 原只 Gatherv u/p/T，MPI 的 VTU 无 temperature_solid；
   且 rank0 的 fld_g 未调 `setup_porous_fields`，输出 guard
   `maxval(h_sf*a_sf)>0` 不成立。已补 T_s pack/Gatherv/unpack 与
   rank0 多孔系数初始化。

## 5. 复现

```bash
python3 gen_ltne1d.py ltne_disp.cas
../../bin/uns_solver ltne_disp.cas ltne_disp0.control
mv ltne_disp.vtu ltne_disp0.vtu
../../bin/uns_solver ltne_disp.cas ltne_disp1.control   # 串行 VTU 已归档
mpirun -np 2 ../../bin/uns_solver_mpi ltne_disp.cas ltne_disp1.control
mv ltne_disp.vtu ltne_disp1_mpi.vtu
python3 plot_ltne_disp.py
```

## 6. 结论

LTNE 双温度 × Bear 纵向热弥散组合路径验证通过：两相剖面中位偏差 ≤0.5%、
能量守恒与 Nield 出口热流分配正确；弥散增强流体轴向回传的物理效应
（出口温升下降、两相分配改变）与 4 阶耦合 ODE 参考解一致；
串并行双温度场一致 ~1e-7。
