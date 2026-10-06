# 温度方程实现计划（SIMPLE + PISO）

## Context

当前求解器（阶段①~⑨）已完整实现低速不可压 SIMPLE（稳态）与 PISO（瞬态，含 BDF2），求解 u/p 并输出 .vtu + Ghia 基准对比。用户要求在 SIMPLE 和 PISO 中增加温度方程，同步求解流体温度。

设计决策（用户已确认）：
1. **耦合方式 = 先被动后加开关**：阶段 A 实现被动标量温度方程（T 被流场输运但不反作用于流场）；阶段 B 加 Boussinesq 浮力开关（`boussinesq=on/off`，ρ=ρ₀(1-β(T-T_ref)) 以重力体积力反馈到动量）。
2. **壁面热边界 = 两者都支持**：同时支持定温壁（Dirichlet, T_wall）和定热流壁（Neumann, q_wall）。

能量方程（不可压、低速、被动标量形式，不含黏性耗散/压力功）：
$$\rho c_p \left(\frac{\partial T}{\partial t} + \mathbf{u}\cdot\nabla T\right) = \nabla\cdot(k\nabla T)$$

离散与单个动量分量同构，复用 CSR/迎风对流/过松弛扩散分解/Barth-Jespersen 限量/BiCGSTAB+ILU0：
- 对流：`∑_f cp · fld%flux(i) · T_upwind`（`fld%flux` 已是 ρu·S 质量通量，**勿再乘 ρ**）
- 扩散：`k · |Sf|²/(Sf·d)·(T_c1-T_c0)` 隐式 + 切向交叉项滞后显式（复用 `nonorth_corr` 风格）
- 瞬态：`ρ·cp·V/dt·(T^{n+1}-T^n)` Euler；BDF2 `ap+=3·ρcp·V/(2dt)`，`rhs+=ρcp·V/(2dt)·(4·T_old-T_old_old)`，第一步退化 Euler，复用 `fld%ts_order`

关键约束：注释英文防乱码；CSR 0 偏移；对角贡献累加进 `ap()`；手册同步+xelatex 两遍；任务结束更新 memory-bank 三件套。

---

## 阶段 A：被动标量温度方程

### A1. mod_fields.f90（[src/mod_fields.f90](file:///home/sundong/Fortran_Project/UNSSolverProj/src/mod_fields.f90)）

`fields_t`（L27-38）新增：
- `T(ncells)` — 温度场
- `T_old(ncells)` / `T_old_old(ncells)` — 瞬态历史
- `gt(3,ncells)` — 温度梯度（Green-Gauss）

`init_fields`（L45-81）新增 allocate + 置 0。

`compute_gradients`（L86-128）末尾新增温度梯度段：面值用 `fld%lf` 线性插值 + 边界面调 `bc_face_T`（新增，见 A3），再调 `grad_scalar` 写入 `fld%gt`。

### A2. mod_control.f90（[src/mod_control.f90](file:///home/sundong/Fortran_Project/UNSSolverProj/src/mod_control.f90)）

`ctrl_t`（L43-77）新增字段：
- `cp = 1005.0_dp` — 比热容
- `k_cond = 0.026_dp` — 热导率
- `tref = 0.0_dp` / `beta = 0.0_dp` / `gravity(3) = 0.0_dp` / `boussinesq = .false.`（阶段 B 用，阶段 A 先占位）

`bc_spec_t`（L31-40）新增字段：
- `ttype = 0` — 热边界类型：0=绝热（缺省）、1=定温壁、2=定热流壁
- `tval = 0.0_dp` — 定温壁温度（Dirichlet）
- `qval = 0.0_dp` — 定热流壁热流（Neumann，正=向域内加热）
- `has_tbc_plane = .false.` / `tbc_dir / tbc_coord / tbc_ttype / tbc_tval / tbc_qval` — 平面过滤热边界（仿 `lid` 机制）

`read_control` select case（L120-177）新增 case：
- `cp` / `k_cond` / `tref` / `beta` / `gravity` / `boussinesq`（true/false）
- `tbc = <zone> <ttype> [tval] [qval]` — 全 zone 热边界（调 `parse_tbc`）
- `tbc_plane = <zone> <dir> <coord> <ttype> [tval] [qval]` — 平面过滤热边界（调 `parse_tbc_plane`，仿 `parse_lid`，按 zone 查找已有 wall/inlet bc 并挂 `has_tbc_plane` 等字段）

新增 `parse_tbc` / `parse_tbc_plane` 子程序（仿 `parse_bc`/`parse_lid`，L193-274）。

### A3. mod_bc.f90（[src/mod_bc.f90](file:///home/sundong/Fortran_Project/UNSSolverProj/src/mod_bc.f90)）

`bcgroup_t`（L18-30）新增字段：`ttype/tval/qval/has_tbc_plane/tbc_dir/tbc_coord/tbc_tol/tbc_ttype/tbc_tval/tbc_qval`（镜像 `bc_spec_t` 新字段）。

`build_bc`（L44-139）：复制 `sp%ttype/tval/qval/has_tbc_plane/tbc_dir/tbc_coord/tbc_ttype/tbc_tval/tbc_qval` 到 `gb`；`has_tbc_plane` 时设 `tbc_tol = 1e-4*maxval(coord_range)`（仿 `lid_tol`，L74）。

新增 `bc_face_T(bcs, i, T_cell, T_face, q_face, is_neumann)`：
- `BC_WALL`：若 `has_tbc_plane` 且面在 `tbc_dir/tbc_coord` 平面（容差 `tbc_tol`），用 `tbc_ttype/tbc_tval/tbc_qval`；否则用 zone 的 `ttype/tval/qval`
  - ttype=0（绝热）：`is_neumann=.true., q_face=0`（零梯度）
  - ttype=1（定温）：`T_face=tval, is_neumann=.false.`
  - ttype=2（定热流）：`is_neumann=.true., q_face=qval`
- `BC_SYMMETRY`：零梯度，`is_neumann=.true., q_face=0`
- `BC_VINLET`：定温（缺省 `tval`），`T_face=tval, is_neumann=.false.`
- `BC_POUTLET`：零梯度，`is_neumann=.true., q_face=0`

报告段（L122-137）补打印热边界信息。

### A4. mod_simple.f90（[src/mod_simple.f90](file:///home/sundong/Fortran_Project/UNSSolverProj/src/mod_simple.f90)）

**新增 `temperature_assembly(m, g, ctrl, bcs, fld, A, rhs, ap_T)`**（仿 [momentum_assembly](file:///home/sundong/Fortran_Project/UNSSolverProj/src/mod_simple.f90#L516-L683) 但去掉压力项、用 `cp/k_cond`）：
- `A%val=0; rhs=0; ap_T=0`
- 可选 Barth-Jespersen 限量（调 `temperature_limiter`，仅 `conv_blend>0`）
- 内部面（仿 L544-608）：
  - `F_T = ctrl%cp * fld%flux(i)`（**勿再乘 ρ**；fld%flux 已含 ρ）
  - `D_T = ctrl%k_cond * area²/(sf·d)`（过松弛；`nonorth_corr=0` 退化 `k*area/dn`）
  - 迎风对流：`ap_T(c0)+=max(F_T,0)+D_T`；`A(c0,c1)-=D_T+max(-F_T,0)`
  - 延迟二阶修正（`conv_blend>0`）：`rhs(c0)-=blend*F_T*dc`（`dc=psi_T(c0)*gt(:,c0)·(xf-xc)`）
  - 非正交修正（`nonorth_corr>0`）：`gf=lf*gt(:,c0)+(1-lf)*gt(:,c1)`；`noc=nonorth_corr*k_cond*gf·(sf-dvec*area²/sd)`；`rhs(c0)+=noc`；`rhs(c1)-=noc`
  - **无压力力项**（标量方程）
- 边界面（仿 L611-633）：
  - 调 `bc_face_T(bcs, i, T(c0), T_face, q_face, is_neumann)` 取热边界
  - `dw=|xf-xc0|`；`D_T=k_cond*area/dw`；`F_T=cp*flux(i)`
  - `if (.not. is_neumann)`（Dirichlet，定温壁/入口）：`ap_T(c0)+=D_T+max(F_T,0); rhs(c0)+=D_T*T_face - F_T*T_face`
  - `else`（Neumann，定热流壁/绝热/对称/出口）：
    - 若 `q_face != 0`（定热流壁）：`rhs(c0)+=q_face*area`（正=向域内加热）
    - 出口（`btype==BC_POUTLET`）：`ap_T(c0)+=max(F_T,0)`（迎风外流，零梯度扩散）
- 瞬态项（仿 L646-660，用 `fld%ts_order`）：
  - `ts_order=1` Euler：`at=rho*cp*V/dt; ap_T+=at; rhs+=at*T_old`
  - `ts_order=2` BDF2：`at=rho*cp*V/(2dt); ap_T+=3*at; rhs+=at*(4*T_old-T_old_old)`
- 对角+松弛（仿 L667-679）：
  - `ts_order>0`：`A(diag)=ap_T`（无松弛）
  - `ts_order==0`（稳态）：`A(diag)=ap_T/alpha_T`；`rhs+=(1-alpha_T)/alpha_T * ap_T * T`（稳态欠松弛）
  - **新增 `alpha_T` 参数**（缺省 0.9，标量场欠松弛，仿 `alpha_u`）—— 或复用 `alpha_u`？建议新增 `alpha_T` 避免 T 收敛慢被 u 拖累。**决策：复用 `alpha_u`（保持参数精简，T 与 u 同档松弛即可）**。
- 释放 `psi_T`

**新增 `temperature_limiter(m, g, bcs, fld, psi_T)`**（仿 [velocity_limiter](file:///home/sundong/Fortran_Project/UNSSolverProj/src/mod_simple.f90#L691-L753)）：把 `fld%u(comp,*)` 换成 `fld%T`，`fld%gu(comp,:,:)` 换成 `fld%gt`，边界面 T 用 `bc_face_T` 取。

**`simple_run`（L57-186）插入温度求解**：在第 5 步 `correct_fields` 之后（L153 后）、第 6 步收敛判定之前，新增：
```
call temperature_assembly(m, g, ctrl, bcs, fld, amat, rhs, ap)
call bicgstab_ilu0(amat, rhs, fld%T, ctrl%lin_tol, ctrl%lin_max, itl, resl, ierr)
if (ierr==2) then; ...; ier=22; return; end if
itl_max = max(itl_max, itl)
```
收敛判据仍只看 mass/du（T 是被动），但可加打印 `T_min/T_max`。

**`piso_run`（L201-389）插入温度求解**：
- 第 1 步推进历史（L280-291）补：`fld%T_old_old = fld%T_old; fld%T_old = fld%T`（与 u 同步，在覆写前 shift）
- `n_correct` 循环（L321-340）之后、第 5 步推进时间（L342）之前，新增温度求解（同 `temperature_assembly` + `bicgstab_ilu0`，含瞬态项按 `ts_order`）
- 诊断打印补 `T_min/T_max`

**注意 rho² bug**：[momentum_assembly L549](file:///home/sundong/Fortran_Project/UNSSolverProj/src/mod_simple.f90#L549) `F=ctrl%rho*fld%flux` 实为 ρ²u·S（fld%flux 已含 ρ），ρ=1 时不可见。**本任务不修动量此 bug**（超范围），只在 `temperature_assembly` 用 `F_T=cp*fld%flux`（正确单 ρ）。

### A5. mod_output.f90（[src/mod_output.f90](file:///home/sundong/Fortran_Project/UNSSolverProj/src/mod_output.f90)）

[vtk_write](file:///home/sundong/Fortran_Project/UNSSolverProj/src/mod_output.f90#L142-L168) CellData 段：在 `velocity` DataArray 之后、`</CellData>` 之前，新增 `temperature` 标量 DataArray（与 `pressure` 同结构，每个子四面体复制 `fld%T(c)`）。

### A6. 算例 cases/cavity_temp.control

被动标量 + 顶盖驱动：Re=100 稳态 SIMPLE，zone 5 壁面缺省绝热，用 `tbc_plane` 设底壁 y=0 定温 T=1、顶盖 y=1 定温 T=0：
```
rho = 1.0
mu = 0.01              # Re = 100
alpha_u = 0.7
alpha_p = 0.3
outer_max = 3000
outer_tol = 1.0e-6
lin_tol = 1.0e-6
lin_max = 300
conv_blend = 1.0
nonorth_corr = 1.0
ppe_precond = ic0
cp = 1.0               # 标度方便（被动标量，cp/k 只影响 T 量纲）
k_cond = 0.01          # Pr = mu*cp/k = 1.0
bc = 5 wall
bc = 4 symmetry
lid = 5 2 1.0 1.0 0.0 0.0
tbc_plane = 5 2 0.0 1 1.0    # bottom wall y=0, fixed T=1
tbc_plane = 5 2 1.0 1 0.0    # top/lid y=1, fixed T=0
```
验证：T 场对流输运（底壁热、顶壁冷、由顶盖涡输运），T_min≈0、T_max≈1，max|u|≈0.928 与原 Re=100 一致（被动标量不改流场）。

### A7. 文档 docs/程序使用手册.tex

- 在 `\section{低速求解器}` 下新增 `\subsection{能量方程（被动标量温度）}`：物理方程、离散（对流/扩散/瞬态）、边界处理（定温/定热流/绝热/入口/出口/对称）
- `\subsection{控制参数}` 表加行：`cp` / `k_cond` / `tbc` / `tbc_plane`
- `\section{算例}` 新增 `\subsection{温度场算例（被动标量，Re=100）}`
- `\section{更新记录}` 追加一条（日期/版本/功能/涉及源码/验证）
- `xelatex` 两遍

---

## 阶段 B：Boussinesq 浮力开关

### B1. mod_simple.f90 — 动量加浮力源

在 [momentum_assembly](file:///home/sundong/Fortran_Project/UNSSolverProj/src/mod_simple.f90#L516-L683) 边界面之后、瞬态项之前（L633 后、L646 前）新增 Boussinesq 浮力源（若 `ctrl%boussinesq`）：
```fortran
if ( ctrl%boussinesq ) then
   do kk = 1, m%ncells
      ! buoyancy force per volume = -rho0*beta*(T-Tref)*g_vec
      ! => rhs(comp) += -rho*beta*(T-Tref)*gravity(comp)*V
      rhs(kk) = rhs(kk) - ctrl%rho * ctrl%beta * (fld%T(kk) - ctrl%tref) &
                * ctrl%gravity(comp) * g%vol(kk)
   end do
end if
```
符号：`gravity` 为重力加速度向量（指向地心，如 y 向上则 `gravity=(0,-g,0)`）；T>Tref 时暖流体受正浮力（与重力反向），公式 `-ρβ(T-Tref)g_vec` 对 T>Tref、`g_vec=(0,-g,0)` 给出 +y 方向力（向上）。正确。

### B2. 算例 cases/cavity_boussinesq.control

差温方腔纯浮力驱动（无顶盖）：
```
rho = 1.0
mu = 0.01              # Re = 100
alpha_u = 0.7
alpha_p = 0.3
outer_max = 3000
outer_tol = 1.0e-6
lin_tol = 1.0e-6
lin_max = 300
conv_blend = 1.0
nonorth_corr = 1.0
ppe_precond = ic0
cp = 1.0
k_cond = 0.01          # Pr = 1
tref = 0.5
beta = 1.0
gravity = 0.0 -1.0 0.0  # g_vec pointing down (y up), dimensionless
boussinesq = true
bc = 5 wall
bc = 4 symmetry
# no lid (lid_vel = 0 by default)
tbc_plane = 5 1 0.0 1 1.0   # left wall x=0 hot T=1
tbc_plane = 5 1 1.0 1 0.0   # right wall x=1 cold T=0
```
Ra = g·β·ΔT·L³·ρ·cp/(μ·k) = 1·1·1·1·1·1/(0.01·0.01) = 1e4（温和浮力，稳态 SIMPLE 可收敛）。
验证：左侧热壁浮力上升、右侧冷壁下沉、形成顺时针大涡；T 场分层；max|u| 与 Ra 标度（约 1e-2 量级）。

### B3. 文档

- `\subsection{能量方程}` 下补 `\subsubsection{Boussinesq 浮力耦合}`
- 参数表加 `tref` / `beta` / `gravity` / `boussinesq`
- 算例加 `\subsection{差温方腔（Boussinesq 浮力，Ra=1e4）}`
- 更新记录追加
- `xelatex` 两遍

---

## 实现顺序

1. A1 mod_fields（T/T_old/T_old_old/gt + init + compute_gradients 温度段）
2. A2 mod_control（ctrl_t/bc_spec_t 新字段 + read_control 解析 + parse_tbc/parse_tbc_plane）
3. A3 mod_bc（bcgroup_t 新字段 + build_bc 复制 + bc_face_T + 报告）
4. A4 mod_simple（temperature_assembly + temperature_limiter + simple_run 插入 + piso_run 插入）
5. A5 mod_output（vtk_write 加 T CellData）
6. `make` 编译 + 跑 cavity_temp 验证（被动标量：max|u| 不变、T_min/T_max 合理）
7. B1 mod_simple（momentum_assembly 加 Boussinesq 源）
8. B2 cases/cavity_boussinesq.control + 跑验证（浮力大涡 + T 分层）
9. A7+B3 手册同步 + xelatex 两遍
10. memory-bank 三件套更新（activeContext/progress/worklog）

## 验证

- **编译**：`make clean && make`，无 warning/error
- **被动标量回归**：`./bin/unsolver Grid_BC/cavity.cas cases/cavity_temp.control` → max|u|≈0.928（与原 Re=100 一致，证明 T 不反作用于 u）、T_min≥0、T_max≤1、收敛 518 次左右
- **稳态 SIMPLE 原算例回归**：`./bin/unsolver Grid_BC/cavity.cas cases/cavity.control` → L2=0.0721 不变（无 tbc 时 T 全 0，不影响流场）
- **瞬态 PISO 回归**：`./bin/unsolver Grid_BC/cavity.cas cases/cavity_transient.control` → max|u|=0.914、L2=0.0823 不变
- **Boussinesq 算例**：`./bin/unsolver Grid_BC/cavity.cas cases/cavity_boussinesq.control` → 顺时针大涡、T 左热右冷分层、收敛
- **VTU 检查**：ParaView 打开 cavity.vtu 显示 temperature 标量场
- **手册**：`cd docs && xelatex 程序使用手册.tex` 两遍无 Error
