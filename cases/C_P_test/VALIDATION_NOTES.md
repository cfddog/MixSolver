# C_P_test 耦合验证记录（2026-10-09）

目标：可压主流（结构，3 块 `Mesh3d.x`）＋ 多孔发汗冷却板（非结构 `unstucutred.cas`）
弱耦合跑通，并验证结构/非结构**间隔轮流、逐步同频**的推进节奏。

## 1 运行方式

```bash
mpirun -np 2 bin/mixsolver_mpi mix.control Mesh3d.x control.ec unstucutred.cas unMesh.control
# 参数序：mix.control | 结构网格 | 结构 control | 非结构 cas | 非结构 control
# 需要 bc3d.inp（结构侧边界定义）与 Mesh3d.x 同目录
```

rank 0 = 结构，rank 1 = 非结构（`nproc>1` 时的默认切分）。

## 2 关键设置

| 项 | 值 | 说明 |
|---|---|---|
| `T_inf`（control.ec） | 309 K | 主流来流（readme：inlet 309 K） |
| `bc = 7 mass-flow-inlet` | `0.1226 300.0` | 冷却剂 **300 K 必须显式给**（见 §4.1） |
| `n_uns_steps / n_struct_steps` | 50 / 50 | 稳态步频（两侧相同 = 同频） |
| `n_struct_steps_start` | 1000 | 结构侧冷启动 1000 子步，每耦合步折半 |
| `iface_velocity` | `balance` | 界面按自身发展型出口＋总通量配平（见 §4.2） |
| `iface_p_anchor` | 1 | 表压基准对齐到结构侧界面压力（见 §4.3） |
| `n_couple` | 30 | |

结构/非结构步频实测：`1000 → 500 → 250 → 125 → 62 → 50 → 50 …`（同频后保持 50）。

## 3 结果（`iface_velocity=balance`＋`iface_p_anchor=1`，30 耦合步，全程无 NaN）

| 量 | 结构侧 | 多孔侧 |
|---|---|---|
| 界面平均 u_x | 102.4 m/s | 0.298 m/s |
| 界面压力 | 1.0148e5 Pa（≈远场） | 表压 ≈118–123 Pa |
| 界面平均温度 | 300.13 K | **300.12 K**（修复前 286.5 K，见 §4.4） |
| 板内 `T` 范围 | — | 299.63 … 303.46 K（修复前 286.08 … 300.00） |
| 板内 `|u|max` | — | 0.315 m/s（物理量级） |
| 每批 SIMPLE 收敛 | — | `du_max→1.4e-5`，`dT_max→1.1e-4` |
| 界面质量配平 | — | `ΣF_iface = +0.0625 m³/s` ＝ 冷却剂注入量，`du_n→7e-6` |

物理图像（修复后自洽）：主流 `T_inf = 309 K`，流经受发汗吹出的板面被冷却到 ≈300 K
（界面平均 300.13 K）；冷却剂 300 K 在板内吸收壁面热流后升温到最高 303.5 K。
修复前板内出现 286 K = 低于全部边界温度的反物理结果。

## 4 过程中发现并修复的四件事

1. **冷却剂入口温度缺省 0 K**：`mass-flow-inlet` 第 4 个 token 是入口静温，
   不写则取 0 K（`mod_uns_control.f90` 注释：能量方程必须给）。原写法
   `bc = 7 mass-flow-inlet 0.1226` 会把多孔板抽向 T→0（日志 `Tmin: 6.95e-16 …`）后 NaN。
2. **界面速度不能整矢量施加**：旧行为把主流 ~104 m/s 切向滑移压进
   `K = 5e-7 m²` 的多孔板，Darcy 阻力要求 `dp/dx = μu/K ≈ 5e5 Pa/m`，
   板内 `|u|max` 冲到 1000 m/s、表压 1.7e6 Pa → 发散。
   `iface_velocity = balance` 下板内 `|u|max = 0.32 m/s`。
3. **界面压力基准漂移**：封闭多孔板的 PPE 是纯 Neumann（除入口/界面无压力边界），
   表压基准是自由参数；而结构侧把 `p_uns(表压)+p_ref` 当**绝对**背压用。
   不锁定时单侧漂到 +1.2e5 Pa，结构侧"吸入"多孔板（`su·n → −0.65 m/s`）
   并把配平项顶到 33 m/s。`iface_p_anchor = 1` 每步整体平移表压场对齐（梯度不变）。
4. **（求解器缺陷，已修源码）能量方程"定温＋出流"边界面把流出焓扣了两次**：
   `temperature_assembly` 的边界 Dirichlet 分支原写
   `rhs += DT*T_face - FT*T_face`；出流（`FT = cp·ρu·S > 0`）时的迎风值本就是内部
   `T_c0`，已由 `ap += max(FT,0)` 隐式计入，多出的 `-FT·T_face` 是**与流出焓同量级的
   假热汇**（均匀 300 K 场残差 = `FT·T_face`；纯对流 1 胞模型给出
   `T_c = T_in − T_face`，应为 `T_in`）——这正是板内 286 K 的来源。
   **判据**：把 zone 4 改名绕开界面登记、改给普通 `velocity-inlet` 出流
   （同网格/同流量/同 300 K，完全不经耦合层）单跑，仍复现
   `it=1 dT_max=7.21`、`it=50 Tmin=286.0`（与耦合 run6 的 7.18/286.1 同量）。
   **修法**：`- FT*T_face` → `- min(FT,0)*T_face`（入流项不变）。
   **验证**：判据算例 `Tmin=Tmax=300.00`、`dT_max 7.21→2.84e-2`；
   本耦合算例界面两侧 `300.13/300.12 K`（此前 300.13/286.5），
   **速度/压力场与修复前逐位相同**（只去掉假热汇）；
   `cases/porous_plug` 与 `regress/m6wing` 回归通过。

## 5 遗留（已量化，未动）

- **动量方程同形项**：`mod_uns_simple.f90:1058`（定温出流面）与 `BC_FARFIELD` 入流支
  有同源的 `- F*uf` 重复计入，出流时应为 `- min(F,0)*uf`。本算例界面该假源项
  `ṁ|u_t| ≈ 4.3e-4 N`，仅为 Darcy 阻力（`≈3.8e-2 N`）的 ~1%，故本次未动，
  以免扰动已标定的界面动量配平（`du_n`）。
- **界面条件的模型层面升级**（文献：Betchen 2006 Eq.17/42/43/44、Zhang 2011
  Eq.27–31、Ochoa-Tapia & Whitaker 应力跳变）：
  ① 界面温度取孔隙率加权 `ε⟨T_f⟩+(1−ε)⟨T_s⟩` 且热流按面积比 ε 分配 ⇒ 只有把
  `thermal_model=ltne` 用到跨组界面时才需要（当前 LTE 的纯 Dirichlet 与 Betchen
  Eq.17 一致）；② 界面压力应由法向动量平衡＋`(1−ε)/ε` 面积突变动力压项给出，
  且 p–ṁ 需子迭代（当前 `iface_p_anchor` 只是整体基准平移）；③ 界面切向速度应按
  应力跳变条件（Beavers-Joseph 型）实现，而非整矢量 Dirichlet。

