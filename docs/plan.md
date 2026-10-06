# MixNSSolver 混合网格 CFD 求解器 -- 编程计划

## 项目概述

将两个独立求解器集成为统一的混合求解器：
- **结构网格部分**：基于 OpenCFD-EC，采用 Riemann 求解器，适用于可压缩流
- **非结构网格部分**：基于 UNSSolverProj，采用 SIMPLE/PIMPLE，适用于低速流
- **耦合方式**：弱耦合（分步迭代），交界面传递边界条件
- **语言**：Fortran 2008, gfortran 8.3.0+

## 目录结构

```
MixNSSolver/
├── src/
│   ├── common/              -- 公共模块
│   │   ├── mod_precision.f90      -- 统一精度定义 (dp, sp)
│   │   ├── mod_constants.f90      -- 通用常数 (PI, gamma, R_gas 等)
│   │   └── mod_unit_convert.f90   -- 量纲转换工具 (无量纲 <-> SI)
│   │
│   ├── structured/          -- 结构网格求解器 (源自 OpenCFD-EC)
│   │   ├── mod_struct_types.f90     -- 结构网格数据类型 (Block_TYPE, BC_MSG_TYPE)
│   │   ├── mod_struct_grid.f90      -- 网格读取 (Plot3D / Gridgen .inp)
│   │   ├── mod_struct_geometry.f90  -- 几何量计算 (体积、面积、法向量)
│   │   ├── mod_struct_scheme.f90    -- 空间格式 (MUSCL/WENO/CD)
│   │   ├── mod_struct_flux.f90      -- Riemann 通量 (HLL/HLLC/Roe/AUSM)
│   │   ├── mod_struct_residual.f90  -- 残差计算
│   │   ├── mod_struct_bc.f90        -- 边界条件 (含 interface 类型)
│   │   ├── mod_struct_turbulence.f90-- 湍流模型 (SA/SST/BL)
│   │   ├── mod_struct_time.f90      -- 时间推进 (RK3/LU-SGS/Dual-time)
│   │   ├── mod_struct_init.f90      -- 流场初始化
│   │   ├── mod_struct_io.f90        -- I/O (Plot3D 流场输出)
│   │   ├── mod_struct_mpi.f90       -- MPI 通信 (halo 交换)
│   │   ├── mod_struct_partition.f90 -- 结构网格分块
│   │   └── mod_struct_solver.f90    -- 结构求解器主驱动
│   │
│   ├── unstructured/        -- 非结构网格求解器 (源自 UNSSolverProj)
│   │   ├── mod_uns_mesh.f90         -- 非结构网格数据类型与读取 (CAS 格式)
│   │   ├── mod_uns_connectivity.f90 -- 面-单元连接关系
│   │   ├── mod_uns_geometry.f90     -- 几何量计算
│   │   ├── mod_uns_linsolver.f90    -- 线性方程组求解器
│   │   ├── mod_uns_simple.f90       -- SIMPLE 算法
│   │   ├── mod_uns_pimple.f90       -- PIMPLE 算法 (瞬态)
│   │   ├── mod_uns_bc.f90           -- 边界条件 (含 interface 类型)
│   │   ├── mod_uns_fields.f90       -- 流场变量
│   │   ├── mod_uns_output.f90       -- 输出 (VTU/Tecplot)
│   │   ├── mod_uns_partition.f90    -- METIS/ParMETIS 分块
│   │   ├── mod_uns_halo.f90         -- Halo 交换
│   │   ├── mod_uns_local_mesh.f90   -- 局部网格提取
│   │   └── mod_uns_solver.f90       -- 非结构求解器主驱动
│   │
│   ├── coupling/            -- 耦合模块 (新建)
│   │   ├── mod_interface_def.f90    -- 交界面定义 (面片匹配关系)
│   │   ├── mod_interface_match.f90  -- 交界面网格匹配 (搜索最近点/投影)
│   │   ├── mod_interface_exchange.f90 -- 数据交换 (MPI 跨求解器通信)
│   │   └── mod_interface_units.f90  -- 量纲统一 (各自转 SI -> 交换 -> 转回)
│   │
│   └── main.f90             -- 主程序入口
│
├── cases/                   -- 测试算例
│   └── ...
├── lib/                     -- 预编译库 (metis, parmetis, tecplot)
├── docs/
└── Makefile
```

## 分阶段实施计划

### 阶段 1：项目骨架与公共模块

**目标**：建立项目目录结构、统一精度定义和构建系统。

**任务**：
1. 创建 `src/common/`, `src/structured/`, `src/unstructured/`, `src/coupling/` 目录
2. 编写 `src/common/mod_precision.f90`：
   - 统一定义 `dp` (双精度), `sp` (单精度)
   - 对应 MPI 数据类型常量
   - 注意：OpenCFD-EC 使用 `PRE_EC=8`，UNSSolverProj 使用 `dp`，需统一为 `dp`
3. 编写 `src/common/mod_constants.f90`：
   - 气体常数 R_gas, 比热比 gamma, PI 等
   - 参考量（rho_ref, T_ref, p_ref, u_ref）用于量纲转换
4. 编写顶层 `Makefile`：
   - 支持 `make structured` (仅结构求解器)
   - 支持 `make unstructured` (仅非结构求解器)
   - 支持 `make all` (混合求解器)
   - 支持 `make mpi` (MPI 并行版本)
   - 链接 lib/metis, lib/parmetis, lib/tecplot

---

### 阶段 2：结构网格求解器集成 (OpenCFD-EC)

**目标**：将 OpenCFD-EC 代码模块化并集成到 `src/structured/`。

**任务**：
1. **模块拆分**：将 OpenCFD-EC 的单体文件重构为独立模块
   - `sub_modules.f90` -> `mod_struct_types.f90` (类型定义) + `mod_constants.f90` (常数)
   - `sub_geometry.f90` -> `mod_struct_geometry.f90`
   - `sub_scheme.f90` -> `mod_struct_scheme.f90`
   - `sub_flux_split.f90` -> `mod_struct_flux.f90`
   - `sub_Residual.f90` -> `mod_struct_residual.f90`
   - `sub_boundary.f90` + `sub_boundary_user.f90` -> `mod_struct_bc.f90`
   - `sub_turbulence_*.f90` -> `mod_struct_turbulence.f90`
   - `sub_time_advance.f90` + `sub_LU_SGS.f90` -> `mod_struct_time.f90`
   - `sub_init.f90` -> `mod_struct_init.f90`
   - `sub_IO.f90` + `sub_Post.f90` -> `mod_struct_io.f90`
   - `sub_partation_mpi.f90` + `sub_update_buffer_mpi.f90` -> `mod_struct_mpi.f90`
   - `sub_read_parameter.f90` -> `mod_struct_grid.f90`
   - `opencfd_ec3d_v1.16a.f90` -> `mod_struct_solver.f90` (驱动逻辑)

2. **接口改造**：
   - 移除全局变量，改为模块变量或通过参数传递
   - 将 `include "mpif.h"` 替换为 `use mpi`
   - 统一精度为 `use mod_precision`
   - 保留原有边界类型定义，新增 `BC_Interface` 类型（用于标记交界面）

3. **验证**：确保结构求解器独立编译运行，结果与原 OpenCFD-EC 一致

---

### 阶段 3：非结构网格求解器集成 (UNSSolverProj)

**目标**：将 UNSSolverProj 代码迁移到 `src/unstructured/`。

**任务**：
1. **代码迁移**：直接复制并重命名模块
   - `mod_mesh.f90` -> `mod_uns_mesh.f90`
   - `mod_cas_reader.f90` -> 合并入 `mod_uns_mesh.f90`
   - `mod_connectivity.f90` -> `mod_uns_connectivity.f90`
   - `mod_geometry.f90` -> `mod_uns_geometry.f90`
   - `mod_linsolver.f90` -> `mod_uns_linsolver.f90`
   - `mod_simple.f90` -> `mod_uns_simple.f90`
   - `mod_simple_mpi.f90` -> 合并入 `mod_uns_simple.f90` (PIMPLE 同理)
   - `mod_bc.f90` -> `mod_uns_bc.f90`
   - `mod_fields.f90` -> `mod_uns_fields.f90`
   - `mod_output.f90` -> `mod_uns_output.f90`
   - `mod_partition.f90` -> `mod_uns_partition.f90`
   - `mod_halo.f90` -> `mod_uns_halo.f90`
   - `mod_local_mesh.f90` -> `mod_uns_local_mesh.f90`
   - `mod_mpi.f90` -> 与结构网格共享 MPI 初始化，或独立为 `mod_uns_mpi.f90`

2. **接口改造**：
   - 统一精度为 `use mod_precision` 中的 `dp`
   - 新增 `BC_INTERFACE` 边界类型（Fluent CAS 的 interface 边界）
   - 非结构网格块类型通过 CAS 文件中的 zone 名称或单独属性文件指定

3. **验证**：确保非结构求解器独立编译运行，结果与原 UNSSolverProj 一致

---

### 阶段 4：交界面定义与匹配

**目标**：实现结构/非结构网格交界面的识别和几何匹配。

**核心数据结构** (`mod_interface_def.f90`)：
```fortran
type :: interface_patch_t
   character(len=64) :: name           -- 交界面名称
   integer :: struct_block_id          -- 所属结构块编号
   integer :: struct_face              -- 结构面方向 (i/j/k min/max)
   integer :: struct_i_range(2)        -- 结构面 i 索引范围
   integer :: struct_j_range(2)        -- 结构面 j 索引范围
   integer :: struct_k_range(2)        -- 结构面 k 索引范围
   integer :: uns_zone_id             -- 非结构 zone 编号 (CAS interface)
   integer :: n_struct_points          -- 结构面节点数
   integer :: n_uns_points             -- 非结构面节点数
   real(dp), allocatable :: struct_xyz(:,:)  -- 结构面节点坐标 (3, n_struct)
   real(dp), allocatable :: uns_xyz(:,:)     -- 非结构面节点坐标 (3, n_uns)
   integer, allocatable :: match_map(:,:)    -- 匹配关系 (2, n_matches)
   real(dp), allocatable :: match_weights(:) -- 插值权重
end type
```

**匹配算法** (`mod_interface_match.f90`)：
1. 读取结构网格 `.inp` 中 `gridgen generic: 8` 类型的交界面定义
2. 读取非结构 CAS 文件中 `interface` 类型的边界
3. 对每对同名交界面：
   - 比较坐标范围，确认几何重叠
   - 对每个非结构面心点，在结构面片上搜索最近的单元/节点
   - 计算双线性插值权重
   - 建立双向映射关系

---

### 阶段 5：交界面数据交换与量纲统一

**目标**：实现跨求解器的交界面数据传递，确保量纲一致性。

**量纲转换** (`mod_interface_units.f90`)：
```
结构网格 OpenCFD-EC 无量纲方式：
  rho* = rho / rho_ref,  u* = u / a_ref,  T* = T / T_ref
  p* = p / (rho_ref * a_ref^2)

非结构网格 UNSSolverProj 使用 SI 单位制：
  rho [kg/m^3], u [m/s], T [K], p [Pa]

交换流程：
  1. 结构求解器：将无量纲量转换为 SI (rho_SI, u_SI, v_SI, w_SI, T_SI, p_SI)
  2. 非结构求解器：直接使用 SI 量
  3. 在交界面传递：
     结构 -> 非结构：rho, u, v, w, T (作为非结构的入口/远场边界)
     非结构 -> 结构：rho, u, v, w, p (作为结构的出口/远场边界)
  4. 各自转换回自己的无量纲形式
```

**数据交换** (`mod_interface_exchange.f90`)：
- 在同一个 MPI_COMM_WORLD 内，结构和非结构求解器共享进程组
- 通过 MPI_Send/Recv 在交界面传递数据
- 支持保守插值（保证质量/动量/能量通量守恒）

---

### 阶段 6：弱耦合迭代驱动

**目标**：实现主程序，交替调用两个求解器并在交界面传递数据。

**主程序流程** (`main.f90`)：
```
1. MPI 初始化
2. 读取全局配置（控制文件、交界面定义）
3. 读取结构网格 + 非结构网格
4. 建立交界面匹配关系
5. 初始化两个求解器的流场
6. 耦合迭代循环：
   do iter = 1, max_coupling_iter
     ! --- 结构网格求解步 ---
     ! 设置交界面边界条件（来自非结构的最新数据）
     call struct_set_interface_bc(interface_patches, fld_uns)
     ! 结构求解器推进 n_struct_steps 步
     call struct_solver_step(struct_grid, struct_fld, n_struct_steps)
     ! 提取交界面数据并转换为 SI
     call struct_extract_interface_data(struct_fld, interface_patches)
     call convert_struct_to_SI(interface_data)

     ! --- 交界面数据交换 ---
     call exchange_interface_data(interface_patches)

     ! --- 非结构网格求解步 ---
     ! 设置交界面边界条件（来自结构的最新数据）
     call uns_set_interface_bc(interface_patches, fld_struct_SI)
     ! 非结构求解器推进 n_uns_steps 步
     call uns_solver_step(uns_mesh, uns_fld, n_uns_steps)
     ! 提取交界面数据（已经是 SI）
     call uns_extract_interface_data(uns_fld, interface_patches)

     ! --- 交界面数据交换 ---
     call exchange_interface_data(interface_patches)

     ! --- 收敛检查 ---
     call check_interface_convergence(interface_patches, converged)
     if (converged) exit
   end do
7. 输出结果
   结构：Plot3D 格式 (rho, u, v, w, T)
   非结构：VTU 格式 (rho, u, v, w, p, tf, ts)
   定时保存重启文件
8. MPI finalize
```

---

### 阶段 7：构建系统与测试

**任务**：
1. 完善 Makefile，处理模块依赖顺序
2. 编写单元测试：
   - 量纲转换正确性测试
   - 交界面匹配精度测试（已知解析函数的插值误差）
   - 数据交换守恒性测试
3. 编写集成测试：
   - 简单算例（如结构矩形通道 + 非结构扩张段），验证交界面连续性

---

### 阶段 8：算例验证

**任务**：
1. 设计验证算例：
   - 亚声速通道流（结构) + 低速腔体流（非结构）
   - 验证交界面处质量/压力/速度连续性
2. 对比纯结构/纯非结构基准解
3. 检查量纲一致性（交界面两侧物理量偏差 < 1%）

---

## 后续计划（功能储备，源自 docs/程序功能说明.md 第 8–15 条，2026-10-05 录入）

> 标注规则：[已实现] / [部分实现] / [待办]。编号承接阶段 8。

### 阶段 9：多孔介质流体计算能力（需求第 8、12、15 条）

1. **块属性识别（VC 类型）** [已实现，2026-10-05]
   - `.control` 支持 `cell_zone = <id|name> [fluid|porous] [perm=.. inertial=.. porosity=.. k_s=.. cp_s=.. rho_s=.. h_sf=.. a_sf=..]`，
     CAS 体区域解析与 resolve_cell_zones 完成（默认 fluid）。
   - 已实现自动解析 CAS zone 名中的 VC 标记（`VC:porous` / `VC: fluid`，
     大小写不敏感；CAS 分词器不处理引号，网格中的标记须无空格）。
     优先级：显式 cell_zone 类型词 > VC 标记 > fluid 默认；类型词可省略
     （CZ_AUTO：只给系数，类型取自 VC 标记）。多孔区缺系数行时启动告警。
   - 验证：cases/porous_plug（source=VC tag）；porous_Ra10 回归 source=control。
2. **多孔源项与输运系数修正** [已实现，2026-10-05]
   - 已有：Darcy-Forchheimer 动量汇（Darcy 隐式对角 + Forchheimer 显式）、
     Brinkman mu_eff=mu/eps、孔隙率加权 k_eff/(rho cp)_eff、LTE/LTNE 双能量方程。
   - 已实现（2026-10-05）：**各向异性渗透率**（轴对齐对角张量 perm_xx/yy/zz，
     未给分量回退标量 perm；交叉项需块耦合动量装配，不支持）；
     **热弥散**（Bear 分量式张量 D_dd=ρcp/|u|·(α_L u_d²+α_T(|u|²−u_d²))，
     面扩散加 n·D·n 投影；LTE 加 k_eff、LTNE 加流体相 kf）。
     另修复 velocity-inlet 行可选静温 token 不被解析的 bug（独立运行温度排空根因）。
   - **2026-10-05 修复纵向弥散公式 bug**：纵向分量原缺 /|u|（量纲 m²/s²、
     结果放大 |u| 倍），横向本有 /umag；修正后纯横向算例（porous_disp、
     porous_graetz 四工况）VTU 位级一致。另修 MPI gather 缺 temperature_solid
     （mod_uns_gather 补 T_s Gatherv；rank0 fld_g 补 setup_porous_fields）。
   - [已实现，2026-10-06] **α/K 标定实验设计**：cases/calibration/，
     数值虚拟标定（合成观测→LM 反演→FIM 可辨识性）。六参数按物理特征
     解耦为四子问题：P1 多流速压降→K+C_F（1% 噪声 CRB 0.6%/5.5%，
     相关 0.65）；P2 混合层横向剖面→α_t（**必须 ΔT≥50 K**，TC 噪声按绝对
     温度缩放，0.2% 噪声 CRB 2.2%）；P3 LTNE 发汗冷却两相测温→α_l+h_sf
     联合反演（0.5% 噪声 CRB 1.5%/0.9%，相关 0.76）；P4 界面速度剖面→
     α_BJ（2% PIV 噪声 CRB 5%）。四子问题全部从 ~2× 离真值初值 4–6 步
     收敛、恢复误差在 CRB 界内。反演器与数据源解耦，可接求解器 VTU 采样。
3. **验证算例**
   - [已实现] 一维非平衡传热：cases/ltne/（LTNE 发汗冷却，vs 解析解，2026-10-04）。
   - [已实现] 二维自然热对流（含多孔方腔）：cases/natconv/（Darcy 侧壁加热
     Nu vs 文献，2026-10-04）。
   - [已实现，2026-10-05] **二维多孔介质强制对流（Graetz）**：
     cases/porous_graetz/，半平行板通道（中心线对称+恒壁温）Darcy 活塞流，
     热入口段 vs Graetz 级数解：充分发展 Nu_a=π²/4=2.4674，
     纯分子导热 0.42%、叠加横向热弥散 1.38%，剖面偏差≤0.026（入口）/≤0.013；
     两工况体均温度在 x* 坐标坍并，验证弥散以 κ+α_t·U 增强横向输运。
   - [已实现，2026-10-05] **porous-plug** 算例：cases/porous_plug/，
     一维均匀多孔塞，Darcy Δp=46.15 Pa 误差 0.000%、Darcy-Forchheimer
     Δp=520.35 Pa 误差 0.002%（内部压力斜率 vs 解析解）；
     同日追加各向异性渗透率验证（张量一致性 0.085%、Kxx 减半变参 0.122%）。
   - [已实现，2026-10-05] **热弥散混合层**算例：cases/porous_disp/，
     2D 多孔温度混合层，α_t=1e-3 m 时 T 剖面 vs erf 相似解偏差 ≤1%。
   - [已实现，2026-10-05] **Beavers-Joseph (B-J) 问题**：cases/beavers_joseph/，
     `cell_zone` 行新增 `bj_alpha=α`（默认 0=禁用，等值反向切向通量保证动量
     守恒），新增全局 `body_force = fx fy fz`（N/m³）驱动开口槽道。
     关键物理结论：离散 BJ 通量的连续极限是"应力连续+速度跳变"双层模型
     （多孔侧 Brinkman 层可分辨），经典 BJ 滑移公式不是正确参考；与 4×4 双层
     耦合 ODE 精确解对拍：α=0 RMS 0.14%，α=1/2 RMS ~1.1-1.3%（λ/dy=4）；
     MPI np=2 与串行最大相对偏差 2.2e-7。
     阶段 11 的低速-多孔界面复用同一 bj_alpha 机制。
   - [已实现，2026-10-05] **LTNE+纵向弥散 1D 组合验证**：cases/ltne_disp/，
     恒热流发汗冷却，disp_l=0/5mm 两工况对拍 4 阶耦合 ODE 半解析解
     （三次特征根+Nield 微观热流分配，exp(r(x−L)) 防刚性溢出）：
     Tf 中位偏差 0.4%、Ts ≤0.83%，弥散使热量回流入口、出口温升由
     995 K 降至 973 K，与参考一致（0.1%）；np2 双温度相对偏差 ~1e-7。
   - [已实现，2026-10-05] **LTNE Graetz 双温度+横向弥散 2D 组合验证**：
     cases/ltne_graetz/，半通道恒壁温（a=5mm，ε=0.4，k_s=1.0，Bi=0.845）
     α_t=0/2.5e-4 两工况，对拍耦合双相特征值参考（4m 一阶块系统，
     保留两相轴向导热、流体横/纵系数分开、入口固相绝热）：
     剖面最大偏差 3.3%/1.4%（仅 x=20mm 入口段），x 向加密（360×40）
     减半至 1.8%/0.7%，x≥50mm ≤0.5%，体均温度全程重合；np2 两相 2e-7。
     至此**阶段 9 全部完成**（含 α/K 标定实验设计，见上第 2 条末款）。

### 阶段 10：流场自动保存与耦合重启（需求第 9 条）

1. **[已实现，2026-10-06]** 自动保存 `save_interval`（mix.control 键，>0 生效）：
   结构侧 `output_flow` 写 flow3d.dat，非结构侧 `write_field_dump` 写
   unMesh_restart.dat；两侧在同一 coupling iter 末尾对齐保存；
   `couple_state.dat` 记录 iter 号，重启时同步恢复（struct uns）。
2. **[已实现，2026-10-06]** 结构侧节点插值 Plot3D 输出 `output_flow_nodes`：
   单元中心（含 ghost）插值到 Mesh3d.x 节点分布，格式自动跟随
   Mesh_File_Format（1=ascii，其它=unformatted）；
   由 `struct_solver_save` 在自动保存时一并调用。
3. **[已实现，2026-10-06]** 耦合驱动层联合重启入口：
   - `couple_restart=1` 时：struct 侧 `force_restart=.true.` 强制 `Iflag_init=1`
     读 flow3d.dat；uns 侧 `restart_file='unMesh_restart.dat'` 串行读 dump
     （`serial=.true.` 跳过 MPI_Bcast，避免耦合驱动中 COMM_WORLD 死锁）；
   - `couple_state.dat` 恢复 iter 号，循环从 `iter0+1` 继续，`iface_ramp`
     计数连续。
4. **验证（grid_BC，6 iter 连续 vs 3+3 重启）**：
   flow3d.dat 与 unMesh_restart.dat 均 **byte-identical**；ser/mpi 双树
   编译回归无新增 error。归档 `grid_BC/cont6/` 与 `grid_BC/restart33/`。

### 阶段 11：多类型交界面边界条件（需求第 14 条）

1. 三类界面处理，界面类型由**两侧块属性自动确定**（域归属 2026-10-06
   用户确认）：
   - [x] **C1 低速（纯流体）–多孔界面**：**只存在于非结构域内部**
     （fluid/porous cell zone 之间的内部面），不跨 solver、不经跨组
     交换层；直接复用阶段 9 的 BJ Robin 通量（bj_alpha），2026-10-06
     确认勾销，详见下方 C1 记录；
   - [x] **可压缩（struct）–低速（uns）界面**：即流 B（couple_channel，
     Dirichlet-Neumann 特征界面，2026-10-06 完成）；
   - [x] **可压缩（struct）–多孔（uns）界面**：跨组界面；2026-10-06 完成并验证，
     详见下方 C2 记录 ＋ `cases/couple_porous/README.md`。
2. [待办] 界面分派表：按两侧 cell_zone 类型（fluid/porous）与求解器类型
   （struct 可压 / uns 低速）在匹配阶段确定界面类型与交换量清单
   （状态量 Dirichlet / 通量型 / 跳跃条件）。
3. [待办] **`uns` 单求解器绝对压力水平 ≈−250 Pa 内部偏置的机理定位与修复**
   （C2 验收时暴露，2026-10-06 登记）。现象与四组对照取证：`cases/couple_porous/
   README.md` §5（复现配方 §6.1），工具 `uns_full/plane_profile.py`。
   特征：内部整段常数下移（同网格全流体对照内部严格平 −251.17 Pa）、紧邻入口边界
   那一层保持物理解水平、出口末列向上回收、凹陷深度严格 ∝u²；与多孔无关、与出口
   BC 类型无关（`outflow` 仅差 3 Pa）。
   候选根因：① `mass-flow-inlet` 在 collocated ＋ Rhie-Chow 下与投影步的相容性源项；
   ② 压力水平锚点（`cases/channel` §3.3 记录 `outflow` 的 PPE 为纯 Neumann 并 pin
   cell 1，而本网格 cell 1 正是入口邻格），但 pressure-outlet 下（PPE 有 Dirichlet）
   偏置仍在 ⇒ 锚点说不充分。
   影响面：**所有** uns 单求解器算例的**绝对**压力（既有验证只看梯度/型线，故一直
   未暴露）；入手处 `src/unstructured/` 投影步与 BC 层。不与阶段 12/13 冲突，可独立
   进行。约束（见 `.trae/rules/project_rules.md`）：修复前跨求解器一律**不得**比对
   绝对压力，只比界面连续性＋梯度。
4. 前置依赖：阶段 7/8 低速 fluid-fluid 耦合收敛验收；阶段 9 的 B-J 验证。

#### 2026-10-06 C1 低速-多孔界面：确认复用阶段 9 BJ 机制，勾销

- 用户明确三类界面的域归属：低速-多孔界面两侧均在**非结构域**内
  （低速流体 uns + 多孔 uns），属求解器内部面，不经过耦合跨组交换层；
  可压缩-低速界面为 struct↔uns（流 B，已完成）；可压缩-多孔界面为
  struct 可压缩 ↔ uns 多孔（**2026-10-06 完成并验证，见下方 C2 记录**）。
- 结论：C1 无新代码、无新算例——直接复用阶段 9 已实现并验证的内部
  BJ Robin 通量（mod_uns_simple momentum_assembly：切向
  C=μA(α/√K)/(1+α·d_Pf/√K)、法向两点扩散、交叉分量滞后等值反向），
  验证资产 cases/beavers_joseph/（对拍双层 4×4 ODE：α=0 RMS 0.14%，
  α=1/2 RMS 1.09%/1.29%；np2 vs 串行 2.2e-7）。
- 下一步：可压缩（struct）–多孔（uns）跨组界面。

#### 2026-10-06 C2 可压缩（struct）–多孔（uns）跨组界面：完成 ✅

- **算例**：`cases/couple_porous/`（struct 流体段 x=0–100 mm ↔ uns 多孔床
  x=100–200 mm 弱耦合，界面 x=100 mm）；对拍资产为单求解器全流域参考
  `uns_full/`（一个 uns：zone 2 fluid + zone 8 VC:porous）。**求解器本体零改动**
  ——两处缺陷都在网格生成/后处理：
  ① `uns_full/gen_full.py` 的 cell-id 排序把"x 串联"写成了"z 并联"
  （`1+i+j*NX+k*NX*NY`，k 最慢 ⇒ zone 2/8 各成一条 z 薄片）：床梯度只有解析值
  ~60%（1.78 vs 3.02 Pa/mm）、max|u|=42.87 > 入口、max p=398 Pa，且 zone 互换后
  逐位相同（z 镜像 ⇒ "换 zone 无影响"假象）；改为 `1+k+j*NZ+i*NY*NZ`（x 最慢）
  即修复，旧网格留档 `scratch/unMesh.cas.zsplit_bug`。
  ② `compare_iface.py` 取窗"中心 ±1.25 mm"跨了两个 dx=2 mm 平面，把 x=101 mm
  （295 Pa）与 x=103 mm（70 Pa）平均成 182.6 Pa ⇒ 幻影 112 Pa 跳变；改用绝对窗
  `[100,102) mm` 后界面压力跳变 **−0.08 Pa**。
- **验收**（`cases/couple_porous/README.md` §2–§4）：床梯度 −3021.3 vs 解析
  −3021.7（0.01%），`check_bed_gradient.py` PASS 且兼作回归守卫（旧网格上主动
  报错）；界面连续（struct 面 295.0 == uns 首排 295.0，跳变 −0.08 Pa；struct 面
  速度 −1.37%）；耦合侧床梯度 −2959.2 Pa/m 与"界面速度 34.233 代入
  μu/(εK)+ρβu² 得 2959.2 Pa/m"逐位吻合 ⇒ 剩余 −2.1% 全部来自 struct 侧界面
  速度的 −1.37%（既有壁面层残余），非界面通量缺陷。
- **新记录（跨算例共性，已列开放项）**：uns 单求解器的**绝对压力水平**存在
  ≈−250 Pa 内部偏置（内部整段常数下移 ＋ 紧邻入口边界那一层保持物理解水平
  ＋ 出口末列回收）；与多孔无关（全流体对照内部严格平 −251.17 Pa）、与出口
  BC 类型无关（`outflow` 只差 3 Pa）、严格 ∝u²。故 C2 判据取"界面连续性 ＋
  床梯度"，**不**与参考做绝对压力对齐；凡跨求解器比绝对压力都需先修正此项。
  取证见 README §5 与 §6.1，工具 `uns_full/plane_profile.py`（新建）。
- **涉及文件**：`cases/couple_porous/{README.md(新建), compare_iface.py,
  gen_meshes.py(注释), uns_full/{gen_full.py, check_bed_gradient.py,
  plane_profile.py(新建)}}`、`cases/couple_channel/{gen_meshes.py,
  uns_full/gen_full.py}`（同款 cid 排序注释；单 cell zone ⇒ 无害）。

#### 2026-10-06 方向修正：彻底弃用 IF_InnerFlow / IF_TurboMachinary 全套机制

- 背景：流 B（couple_channel 界面不连续）排查中一度尝试 struct 侧总压入口
  （IF_InnerFlow=1, P_In_Ratio=1.0123），np2 跑到 iter50 速度升至 78 m/s 发散，
  界面偏差依旧——证明**界面连续性问题与入口条件无关**，只能在耦合交换层内、
  用两侧状态解决。用户拍板：IF_InnerFlow 只保留外流部分；IF_TurboMachinary
  叶轮机模式全部弃用。
- 删除内容（src/structured/）：
  - mod_struct_global.f90：IF_TurboMachinary / IF_InnerFlow / P_In_Ratio /
    Turbo_Periodic_seta / Turbo_w 变量；
  - mod_struct_constants.f90：BC_Wall_Turbo=201；
  - mod_struct_bc.f90：Turbo/Inflow/Outflow 分派分支与 boundary_BC_Inflow_Turbo、
    boundary_BC_Outflow_Turbo、boundary_wall_Turbo 三个子程序，只保留
    wall/Farfield/Inflow/Outflow/Symmetry/Extrapolate/user/BC_INTERFACE 外流派生；
  - mod_struct_solver.f90：旋转惯性力（离心+科氏）源项块；
  - mod_struct_init.f90：Turbo_P0/T0/L0/w、Ref_medium_usrdef、namelist 键、
    默认值、bcast 槽位内容（槽位号保留不动）、叶轮机初始化 else 分支、
    默认介质（Ma=1/Re 反推）参数推导块；
  - mod_struct_mpi.f90：Umessage_Turbo_Periodic 子程序及调用点；
    Coordinate_Periodic / Mesh_Center_Periodic 删除旋转角分支，只保留平移周期。
  - 保留：AoS（侧滑角，原与 Turbo_* 同行声明，case-insensitive，曾误删后恢复）、
    Periodic_dX/dY/dZ 平移周期。
- 算例：couple_channel/control.ec 恢复强制均匀入口（Kstep_save=5000）；
  couple_channel 与 grid_BC 的 control.ec 删除 5 个废弃 namelist 键
  （键移除后 gfortran 不接受未知键）；删除 control.ec.forced、
  mix.control.ptot、run_ptot1.log 试验产物。
- 验证：ser/mpi 双树编译通过；M6-wing 对拍 external 原始基线
  （/tmp/m6reg，同 -O2 -std=legacy 标志，t_end=0.501 恰好 50 步），
  np1 与 np2 的 flow3d.dat / Step_mess.dat 均字节一致；couple_channel
  n_couple=2 冒烟通过（注：冒烟覆盖了 unMesh_coupled.vtu，流 B 需重跑）。
- 流 B 后续：界面连续性修法重新定位到交换层（main.f90 反射公式
  iu=2*iface_vel-iu_cell0 是固定亏损来源），参考 uns_full 单一求解器解。

#### 2026-10-06 流 B 交换层修法（Dirichlet-Neumann 特征界面，完成）

迭代过的方案（couple_channel，np2，n_couple=1000）：

1. 删除 uns 侧反射 + struct ghost=2*face-inner（子步前设置）：
   周期-2 翻转，iter75 前 NaN（uns/struct 交替 ±100 m/s）。
2. 同上但 ghost 在 10 个子步**之后**设置：回声不动点——extract 恒等于
   uns 回传值，struct 侧无物理信号，channel 流量单调衰减到 0.04 m/s。
3. ghost 直写 uns 首排格心（零阶保持）：稳定冻结，但界面压力跳变
   1.7 kPa（struct 面 1990 Pa vs uns 格心 323 Pa）、质量 31.2 m/s
   （参考 34.72），末排格心 i=50 出现 3657 Pa 伪压力峰。
4. **最终方案**：交换层改为 Dirichlet-Neumann 分区——
   - struct 侧（亚声速出口）只接收 uns 背压，ghost 用与
     boundary_Farfield 亚声速出口相同的线性化 Riemann 反射构造
     （db=d1+(pb-p1)/c1^2，ub=u1+(p1-pb)/(rho c)·n，
     ghost=2 face-inner）；法向量取自 Block 面法向表（face 1..6 分派）。
   - uns 侧只施加速度 Dirichlet；界面压力零梯度（bc_face_p 对
     BC_INTERFACE 本就 pf=pP），删除 set_interface_p 调用与压力重锚定。
   - 背压松弛 alpha=0.3·min(1,iter/iface_ramp)：Ma=0.1 下 1/(rho c)
     增益大，alpha=1 在 ramp 结束后周期发散；alpha<1 不改变收敛值
     （p1 每次交换都向 pb 松弛）。

最终结果（iter 1000，run_flowb4c.log，无 OverLimit，~100 iter 冻结）：
   - 质量连续：struct 面/内部 33.71 m/s == uns 下游 33.71（全域常数）；
   - 压力连续：struct 界面 342.5 Pa == uns 首排 342 Pa（跳变 <1 Pa，
     旧方案为 kPa 级冻结跳变）；
   - 对拍 uns_full（x=100：34.72 m/s，46.9 Pa）：均值 -2.9%，
     界面压力 +296 Pa；近壁型线最大偏差 11.8%。
剩余偏差来源是两求解器**壁面层物理差异**（struct 近壁 u=13.1 vs
uns 自然型线 9.05 @1.25mm，BL 发展程度不同），属物理一致性任务
（壁面处理/网格加密对齐），不是交换层缺陷。
另：compare_iface.py 修正 flow3d.dat 解析（output_flow 写
U(0:nx) 含 ghost，bNi=ni=51；i+ 内点索引 ni-1、ghost ni）。
已知遗留（既有，非本次引入）：`make all` serial 树在
mod_uns_restart.f90 失败（无条件 use mod_uns_mpi_core，需 #ifdef
HAVE_MPI 守护）；耦合只用 MPI 构建。

### 阶段 12：Gambit NEU 网格输入（需求第 13 条）

1. [待办] 增加 Gambit `.neu` 文件读取器作为非结构网格与边界条件的第二种输入格式
   （现有：Fluent CAS）。
2. [待办] NEU → 内部 mesh_t 的映射：体区域 fluid/porous 属性、面区域 BC、
   interface zone 登记，复用现有 build_bc / register_interface_zones 流程。

### 阶段 13：结构求解器演进（需求第 10、11 条）

1. [待办/暂缓] 结构 SST 湍流模型存在已知 bug，修复前**不启用**
   （当前算例一律 Iflag_turbulence_model=0 层流）；修复后补 SST 回归算例。
2. [待办，远期] 结构求解器改为 **Liao 的格心型有限差分方法**，
   要求兼容现有 Riemann 通量/求解器接口；切换需保持 M6-wing 等既有回归可对拍。

---

## 关键注意事项

1. **量纲一致性**（功能说明第6条重点强调）：
   - OpenCFD-EC 使用无量纲变量（基于参考量归一化）
   - UNSSolverProj 使用 SI 单位
   - 交界面交换时必须：各自转 SI -> 交换 -> 各自转回
   - 在 `mod_interface_units.f90` 中集中管理转换，便于检查和调试

2. **MPI 通信架构**：
   - 结构和非结构求解器共享 `MPI_COMM_WORLD`
   - 可通过 `MPI_Comm_split` 创建子通信域（结构组/非结构组）
   - 交界面交换使用跨子通信域的 MPI_Send/Recv

3. **非结构网格块类型**：
   - 通过 CAS 文件中的 zone 名称识别（如 `zone_name = "porous_zone"`）
   - 或通过单独的属性配置文件指定每个 zone 的求解器类型

4. **并行策略**：
   - 结构网格：沿用 OpenCFD-EC 原有的块级并行（每个 block 一个或一组进程）
   - 非结构网格：使用 METIS/ParMETIS 分块
   - 交界面交换在对应的进程间直接进行

5. **运行产物不入库**（2026-10-06 瘦身，硬规则）：`.vtu / *.log / *.dat / *.out /
   *.tmp / *.part.map / __pycache__` 等求解器写出文件一律不提交，无论大小
   （`.trae/rules/project_rules.md`「版本控制约定」）。算例目录只保留输入件
   （`*.cas/*.neu/*.cgns/*.x/*.control/mix.control/bc3d.*`）、工具脚本、`README.md`
   与 `images/*.png`；阶段验收的**定值与命令**写进 `cases/*/README.md`，
   靠重算复现，而非靠入库产物。历史剥离用了 `git filter-branch`（本机无
   `git-filter-repo`），安全网为**仓库外 bundle**（98 MB，含完整旧历史）：
   `/home/sundong/mixsolver_pre_slim_backup/repo_pre_slim.bundle`
   （另有 `cases`/`grid_BC` 硬链接快照）。实测：`.git` 120 MB → 46 MB、
   跟踪 324 文件、`main` 与 `origin/main` 0/0。

## 实施优先级

| 优先级 | 阶段 | 预计工作量 | 依赖 |
|--------|------|-----------|------|
| P0 | 阶段1: 项目骨架 | 1天 | 无 |
| P0 | 阶段2: 结构求解器集成 | 3-4天 | 阶段1 |
| P0 | 阶段3: 非结构求解器集成 | 2-3天 | 阶段1 |
| P1 | 阶段4: 交界面定义与匹配 | 3-4天 | 阶段2,3 |
| P1 | 阶段5: 数据交换与量纲 | 2-3天 | 阶段4 |
| P1 | 阶段6: 耦合驱动 | 2天 | 阶段4,5 |
| P2 | 阶段7: 构建与测试 | 2天 | 阶段6 |
| P2 | 阶段8: 算例验证 | 2-3天 | 阶段7 |
| P2 | 阶段9: 多孔介质能力完善（VC 自动识别/输运修正/porous-plug/B-J） | 视算例 | 阶段8 |
| P2 | 阶段10: 流场自动保存与耦合联合重启（Plot3D 节点插值） | 2-3天 | 阶段6 |
| P2 | 阶段11: 多类型界面（C1 低速-多孔已勾销；流 B 已完成；余可压缩-多孔 + 分派表） | 剩余2-3天 | 阶段8,9 |
| P3 | 阶段12: Gambit NEU 网格输入 | 2-3天 | 阶段3 |
| P3 | 阶段13: 结构求解器演进（SST 修复；Liao 格心型 FD，远期） | 远期 | — |
