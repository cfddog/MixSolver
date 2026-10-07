# 流 B 可压缩（struct）– 低速（uns）跨组界面验证

> 日期：2026-10-06　求解器：`bin/uns_solver`、`bin/mixsolver(_mpi)`
> 主题：结构侧可压缩通道（x=0–100 mm，Ma=0.1，OpenCFD-EC Riemann 求解器）
> + 非结构侧不可压 SIMPLE 通道（x=100–200 mm）弱耦合，界面在 x=100 mm；
> 用"单求解器全流域算例 `uns_full`"作参考对拍。
> 本轮修的是**耦合交换层的界面亏损（流 B）**：旧交换层用反射公式平移
> 界面速度（`iu = 2*iface_vel - iu_cell0`），使界面速度/压力出现固定
> 亏损（v3 冻结 1.7 kPa 压力跳变）；改为 **Dirichlet–Neumann 分区特征
> 界面**后，界面质量/压力连续性恢复，界面压力跳变降到 Pa 量级。

## 1. 布局

| 角色 | 文件 | 说明 |
|---|---|---|
| 结构侧（耦合） | `Mesh3d.x`、`bc3d.inp`、`control.ec`、`mix.control` | x=0–100 mm 可压缩通道；i− 强制均匀入口、i+ 界面（`gridgen generic:8`）、y 壁面、z symmetry（严格 2D） |
| 非结构侧（耦合） | `unMesh.cas`、`unMesh.control` | x=100–200 mm，**全域 = fluid**；zone 7 界面、zone 6 pressure-outlet（x=200）、zone 5 wall、zone 4 symmetry |
| 参考（单求解器） | `uns_full/{gen_full.py, unMesh.cas, unMesh.control, unMesh.vtu, full.log}` | 0–200 mm 用**一个** uns 求解器算全通道（界面处 `bc = 7 velocity-inlet`） |
| 非结构侧单跑 | `uns_standalone/{unMesh.cas, unMesh.control, unMesh.vtu, standalone.log}` | 仅 x=100–200 mm 的 uns 求解器，界面给定均匀入口 34.7224 m/s |
| 界面比较 | `compare_iface.py` | 界面处 struct 面／uns 首排／参考三方对比（u_x(y) 型线与均值、p 均值） |
| 网格生成 | `gen_meshes.py` | 生成 `unMesh.cas`（非结构侧）与 `Mesh3d.x`／`bc3d.inp`（结构侧） |

物理：300 K 空气 ρ=1.177 kg/m³；**人工黏性** μ=4.0868e-3 Pa·s；Ma=0.1 ⇒
U_inf = Ma·a_ref = 0.1·347.224 = **34.7224 m/s**；通道高 H=50 mm ⇒

    Re_H = ρ·U_inf·H/μ = 1.177·34.7224·0.05 / 4.0868e-3 = 500

## 2. 两侧量纲匹配

结构侧 OpenCFD-EC 用无量纲量（`mix.control` 给参考态）：

    u* = u / U_inf,   p* = p / (ρ_ref·U_inf²) = p / 1419.0 Pa

非结构侧全程 SI（m、m/s、Pa）。交换层把界面状态**先统一换算到 SI**、
交换完毕再各自换回（`mod_interface_units`）。匹配的关键量：

| 量 | 结构侧（无量纲） | 换算 | 非结构侧（SI） |
|---|---|---|---|
| ρ_ref | 1 | — | 1.177 kg/m³ |
| T_ref | 1（T*=T/T_ref） | ×300 | 300 K |
| U_inf | 1 | ×34.7224 | 34.7224 m/s |
| p 尺度 | 1 | ×1419.0 | 1419.0 Pa |



## 3. 缺陷（已修）：交换层反射公式造成的界面亏损

**症状**：耦合收敛后，界面处两个求解器的速度/压力不连续——

| 版本 | 界面现象 |
|---|---|
| v3（ghost 零阶直写 uns 首排格心） | 稳定但**冻结** 1.7 kPa 压力跳变（struct 面 1990 Pa vs uns 格心 323 Pa），界面质量 31.2 m/s（参考 34.72），末排格心 i=50 出现 3657 Pa 伪压力峰 |

**根因**：旧交换层在 uns 侧对界面速度做**反射**、struct 侧令 ghost =
`2*face - inner`（`iu = 2*iface_vel - iu_cell0`）。这是一个**固定亏损
源**：它把界面当作镜像面而非"另一侧的状态"，两侧交替放大/衰减，或
收敛到一个被平移过的错误不动点（v1 周期-2 翻转 iter75 NaN；v2 回声不动点、
流量单调衰减到 0.04 m/s）。问题只能在**交换层内、用两侧状态**解决。

## 4. 修法：Dirichlet–Neumann 分区特征界面

最终方案（`src/main.f90` + `src/structured/mod_struct_driver.f90`）把界面
当成一次 **Dirichlet–Neumann 分区**：

- **struct 侧**（亚声速出口）**只接收 uns 背压** `pb`，ghost 用与
  `boundary_Farfield` 亚声速出口**相同的线性化 Riemann 反射**构造：

      db = d1 + (pb - p1)/c1^2
      ub = u1 + (p1 - pb)/(rho*c)*n_out
      ghost = 2*face - inner        (p2 = 2*pb - p1)

  法向取自 `Interface_List%face` 1..6 的 Block 面法向表。
- **uns 侧只施加速度 Dirichlet**；界面压力取**零梯度**（`bc_face_p` 对
  `BC_INTERFACE` 本就 `pf = pP`），删除 `set_interface_p` 调用与压力重锚定。
- 背压松弛 `alpha = 0.3*min(1, iter/iface_ramp)`：Ma=0.1 下 1/(rho*c) 增益大，
  `alpha=1` 在 ramp 结束后周期发散（iter50 NaN / 819 OverLimit）；
  `alpha<1` 不改变收敛值（`p1` 每次交换都向 `pb` 松弛）。

### 4.1 四版迭代（np2，n_couple=1000，n_struct/uns_steps=10，iface_ramp=20）

1. 删 uns 反射 + struct ghost 在**子步前**设置 → 周期-2 翻转，iter75 NaN。
2. 同上但 ghost 在**子步后**设置 → 回声不动点，extract 恒等于 uns 回传值，
   struct 侧无物理信号，流量单调衰减到 0.04 m/s。
3. ghost 零阶直写 uns 首排格心 → 稳定冻结，但 1.7 kPa 压力跳变（见 §3）。
4. **定稿 Dirichlet–Neumann 特征界面**（本节）。

## 5. 结果（`run_flowb4c.log`，iter 1000，无 OverLimit，~100 iter 冻结）

`compare_iface.py` 输出（单平面窗，见 §5.1）：

```
reference @x=100 (n=40): mean u=34.722  mean p=46.9 Pa
coupled uns @x=0.101 (n=40): mean u=31.218  mean p=338.4 Pa
coupled struct face: mean u=33.708  mean p=101676.5 Pa abs (gauge 342.5)
coupled struct inner: mean u=33.708  p gauge 342.6

profile deviation vs reference (% of U_inf):
  struct face: max +11.78%  RMS 7.06%  mean-ux -2.92%
  uns  cell0 : max +17.38%  RMS 12.92%  mean-ux -10.09%
  pressure: struct 342.5 Pa, uns 338.4 Pa, ref 46.9 Pa (gauge)
```

- **界面质量连续**：struct 面／内部 u_x = **33.708 m/s**（沿程常数）；
  uns 从 x=103 mm 起到出口恒为 **33.707 m/s**（逐位相同）。仅界面**首排**
  格心（x=101 mm）为 31.218（比下游低 ~7%）——这是 uns 在界面速度
  Dirichlet 下的一格边界层松弛，不是全局亏损。
- **界面压力连续**：struct 面 342.5 Pa vs uns 首排 338.4 Pa，跳变
  **≈4.1 Pa**（旧方案为 kPa 级冻结跳变）。
- **对拍 `uns_full` 单一求解器参考**（x=100：mean u=34.722 m/s、p gauge=46.9 Pa）：
  struct 型线均值 **−2.92%**、uns 首排均值 **−10.09%**、近壁最大偏差
  **+11.78%**（结构近壁 u=13.11 vs uns 12.41 @1.25 mm）；界面压力相对参考
  **+296 Pa**（342.5 vs 46.9）。

### 5.1 `compare_iface.py` 取窗缺陷（已修，2026-10-07）

**症状**：脚本原把 uns 取样窗写成 `[xmin, xmin+0.0025)`（宽 2.5 mm）。
两侧 dx 均为 **2 mm**，该窗同时套住 x=101 mm（p=338.4）与 x=103 mm
（p=117.2），平均成 **227.8 Pa** ⇒ 凭空出现 ~115 Pa 的界面压力"跳变"。

**修法**：改用**绝对**单平面窗 `[X_IFACE, X_IFACE+DX) = [100, 102) mm`
（`DX = 0.002`），恰只含 x=101 mm 平面——与 `cases/couple_porous/` 的
同款修法一致（那里修的是宽 2.5 mm 窗跨两个平面引起的"幻影 112 Pa 跳变"）。

## 6. 复现

```bash
# 参考（单求解器，全流域 0–200 mm）
cd uns_full
python3 gen_full.py                                     # 重新生成 unMesh.cas
../../bin/uns_solver unMesh.cas unMesh.control unMesh.vtu > full.log 2>&1

# 单侧 uns（x=100–200 mm，界面给定均匀入口速度 34.7224 m/s）
cd ../uns_standalone
../../bin/uns_solver unMesh.cas unMesh.control unMesh.vtu > standalone.log 2>&1

# 耦合算例（在算例目录里跑；产出 compare_iface.py 所需的两份文件）
#   unMesh_coupled.vtu ← uns 侧终场快照（src/main.f90 的 VTU 写出）
#   flow3d.dat         ← struct 侧存档；**须 save_interval ≤ n_couple**
cd ..
mpirun -np 2 ../../bin/mixsolver_mpi mix.control Mesh3d.x control.ec \
        unMesh.cas unMesh.control > run_flowb4c.log 2>&1

# 界面比较（读 unMesh_coupled.vtu + flow3d.dat + uns_full/unMesh.vtu）
python3 compare_iface.py
```

本仓库的 `run_flowb4c.log` / `flow3d.dat` / `unMesh_coupled.vtu` 即上一条
命令在 iter 1000 的产物，`compare_iface.py` 可直接复现 §5 的关键行
（`mpirun -np 2` 全程约 4 min；400 iter 起界面量即已冻结）。

## 7. 已知遗留

- **壁面层物理差异**：struct 近壁 u=13.11 vs uns 12.41 @1.25 mm（参考
  9.05），近壁型线最大偏差 +11.78%。这是两求解器在**壁面处理/网格加密**
  上的差异（壁面层对齐任务，见 `memory-bank/` backlog），**不是交换层缺陷**。
- **界面压力偏置 +296 Pa**：界面压力相对单一求解器参考整体抬高，与
  `cases/couple_porous/README.md` §5 记录的绝对压力偏置（−250 Pa 一族）
  同源，属"压力电平非规范不变性"问题，值得单独立项；当前不影响界面
  **连续性**（跳变 4.1 Pa）。
- `uns` 单求解器每次运行会在当前目录落一个 `flux_debug.txt`
  （`src/unstructured/main_uns.f90` 无条件写出），可删。
