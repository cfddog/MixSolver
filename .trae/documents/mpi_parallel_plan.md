# UNSSolver MPI 并行化方案（ParMETIS 域分解）

## 背景

当前 UNSSolver v0.1 为纯串行代码。用户要求通过 ParMETIS 将网格分区为 nzone 个计算域，每个域在一个 MPI 进程上运行，保留分区映射文件，输出时还原为完整原始网格，支持不同核数续算。

## 总体策略

- **最小侵入**：保留现有串行代码路径（`make` 不变），新增并行构建路径（`make mpi`）
- **本地网格抽取**：每个 rank 从全局网格抽取 owned + halo(1层) 单元，复用现有 `build_connectivity`/`compute_geometry`/`build_bc`，求解器主循环结构不变，仅在关键时机插入 halo exchange
- **输出聚合**：rank 0 持有全局 `mesh_t`，`MPI_Gatherv` 聚合字段后调用现有 `vtk_write`/`tecplot_write`（零修改）

## 库环境

- mpich: `/usr/bin/mpif90`, `/usr/bin/mpirun` (Ubuntu 13.3 gfortran)
- ParMETIS: `lib/parmetis/libparmetis_gnu_mpi.a`，导出 `ParMETIS_V3_PartKway`, `ParMETIS_V3_AdaptiveRepart`
- METIS: `lib/metis/libmetis.a`，导出 `METIS_PartGraphKway` 等
- 无头文件 → 手写 `iso_c_binding` 接口

## 7 个子步骤

### 步骤 1：MPI 引导 + ParMETIS 接口 + Makefile

**新建文件**：
- `src/mod_mpi.f90`：MPI_Init/Finalize/Comm_Rank/Comm_Size 薄封装，缓存 rank/nprocs
- `src/mod_parmetis_iface.f90`：iso_c_binding 接口绑定 `ParMETIS_V3_PartKway`、`ParMETIS_V3_AdaptiveRepart`（idx_t 先假设 `c_int32_t`，smoke test 确认）
- `tests/test_parmetis.f90`：最小 smoke test（4 顶点链图，2 proc 分区）

**Makefile 改动**：
- 新增 `mpi` 目标：`FC=mpif90`，链接 `lib/parmetis/libparmetis_gnu_mpi.a lib/metis/libmetis.a -lstdc++ -lm`
- 串行 `make` 路径不变
- 新增模块依赖图

**验证**：
- `make` 串行回归不变（518 次，L2=0.0721）
- `make mpi` 生成 `bin/unsolver_mpi`
- `mpirun -np 2 ./bin/test_parmetis` 对链图分区返回 edgecut=1, part=[0,0,1,1]

---

### 步骤 2：对偶图构造 + ParMETIS 分区 + 映射文件

**新建文件**：`src/mod_partition.f90`
- `build_dual_graph(conn, xadj, adjncy)`：从 `conn%c2c` 构造 0-based CSR dual graph
- `distribute_graph_naive(ncells, nprocs, vtxdist)`：naive block 分发
- `call_partkway(...)`：调 ParMETIS_V3_PartKway
- `gather_partition(...)`：MPI_Gatherv 拼回全局 part[]
- `write_partition_map(fname, part, nprocs, edgecut)`：ASCII 映射文件
- `read_partition_map(fname, part, nprocs_old)`：续算用

**验证**：
- `mpirun -np 2` 生成 `cavity.part`（1244 行，两分区均衡，edgecut 合理）
- `mpirun -np 4` 四分区均衡，edgecut 随 nprocs 递增
- 重复运行结果一致（确定种子）

---

### 步骤 3：本地网格抽取（owned + halo）

**新建文件**：`src/mod_local_mesh.f90`
- `extract_local_mesh(m_global, conn_global, part, my_rank, m_loc, g_to_l, l_to_g, ncells_loc, nghost)`
- 扫描内面：若 `part(c0)==my_rank` 或 `part(c1)==my_rank`，则该面/单元属于本 rank
- owned 单元 → 1..ncells_loc；halo 单元 → ncells_loc+1..ncells_loc+nghost
- 本地节点/面重新 1-based 编号，c0/c1 用本地索引
- **分区边界面** c1=halo（非零），被求解器当作普通内面处理
- 调用现有 `build_connectivity`/`compute_geometry`（零修改）
- `mod_bc.f90` 加一行：分区边界面（c1>0 的非物理边界）fgrp=0

**验证**：
- `sum(ncells_loc) == 1244`（无重无漏）
- `MPI_Allreduce(SUM(local owned vol)) == 1.0`
- 本地 GCL < 1e-12
- 分区边界面总数 == edgecut

---

### 步骤 4：Halo 通信 + 分布式线性求解器

**新建文件**：`src/mod_halo.f90`
- `halo_setup(m_loc, l_to_g, g_to_l, cell_rank, comm, halo_info)`：构造 send/recv 列表
- `halo_exchange_real1d/2d/3d(halo_info, vec, comm)`：非阻塞 Isend/Irecv + Waitall

**修改 `mod_linsolver.f90`**（新增，不破坏串行）：
- `csr_matvec_local(A, x, y, halo_info, comm)`：先 exchange 再算 owned 行
- `bicgstab_ilu0_local(...)` / `cg_ic0_local(...)`：matvec 用本地版，ILU/IC0 只分解 owned 行

**验证**：
- 单元测试：2 proc matvec 与串行 max diff < 1e-14
- 单元测试：2 proc BiCGSTAB/CG 收敛到同 tol
- `make test` 全 PASS

---

### 步骤 5：并行 SIMPLE/PISO

**修改 `mod_simple.f90`**（插入 halo exchange，逻辑不变）：
- `simple_run` 新增 `halo_info, comm` 参数（缺省退化串行）
- halo exchange 时机：① momentum_assembly 前 ② compute_gradients 后 ③ flux_rhiechow 前 ④ ppe_assembly 前 ⑤ correct_fields 后 ⑥ temperature_assembly 前
- `piso_run` 同理，corrector 循环内也要 exchange（含 `pp`）
- 收敛判据：MPI_Allreduce 全局残差
- rank 0 打印收敛历史

**验证**：
- `mpirun -np 1` 与串行 byte-identical
- `mpirun -np 2` Re=100：518±20 次收敛，L2 差 < 1%，max|u| 差 < 0.5%
- `mpirun -np 4` 同上
- PISO/Boussinesq/温度场全部通过

---

### 步骤 6：输出聚合 + 全局网格输出

**修改 `mod_output.f90`**：
- `gather_fields_to_root(m_global, fld_loc, l_to_g, ncells_loc, comm, fld_global)`：MPI_Gatherv
- rank 0 调现有 `vtk_write`/`tecplot_write`（零修改）
- PISO 快照同样聚合 + 输出
- `ghia_compare` 在 rank 0 的全局字段上跑

**验证**：
- 并行 `cavity.vtu` 与串行字段 max diff < 1e-12
- 并行 `cavity.plt` 结构正确
- PISO 快照与串行一致
- Ghia L2 一致

---

### 步骤 7：不同核数重启

**新建文件**：`src/mod_restart.f90`
- `write_field_dump(fname, m_global, fld_global)`：按全局单元 ID 写字段（u, p, T, 历史）
- `read_field_dump(fname, fld_global)`：rank 0 读回
- `restart_with_map(partfile, fieldfile, ...)`：
  - 若 nprocs_old == nprocs：直接用旧 part
  - 若 nprocs_old != nprocs：重新 PartKway 或 AdaptiveRepart

**文件格式**：分区映射文件 ASCII（人眼可读），字段 dump 二进制（效率优先）

**重分区策略**：AdaptiveRepart（利用旧分区，最小化数据迁移）

**控制参数新增**（`mod_control.f90`）：
- `restart`, `restart_file`, `partition_file`
- 变核数时自动调用 `ParMETIS_V3_AdaptiveRepart`（旧 part 作输入）

**验证**：
- 同核数重启：max|u - u_cont| < 1e-10
- 变核数重启（2→4）：max|u - u_4proc| < 1e-6
- 不重启时回归不变
- 手册同步更新（并行章节 + 参数表 + changelog，xelatex 两遍）

## 依赖链

```
Step 1 → Step 2 → Step 3 → Step 4 → Step 5 → Step 6 → Step 7
```

每步独立验证，用户确认后进入下一步。
