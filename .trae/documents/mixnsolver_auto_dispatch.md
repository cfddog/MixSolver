# mixnsolver：零参数、自动分发的统一可执行

## Context（背景）

用户不想在运行算例时敲一长串控制文件/网格文件名，希望**只有一个可执行文件名**，
运行程序自动根据当前目录里的输入文件判断该跑哪种模式，并直接执行：

- **可压缩 / 结构侧**（struct）
- **不可压缩 / 非结构侧**（uns）
- **弱耦合**（coupled = struct + uns）

现状是有多个入口：`struct_solver`（零参数，读 cwd 的 `control.ec`）、
`uns_solver(_mpi)`（需要 `<mesh.cas> <ctl>`）、`mixsolver(_mpi)`（需要 5 个文件名参数）。
目标交付物是一个新的 `bin/mixnsolver`（串）+ `bin/mixnsolver_mpi`（MPI），
**旧二进制全部保留、行为逐位不变**。

## 已确认决策（用户拍板）
1. 新增 `mixnsolver`，保留所有旧二进制。
2. 在 case 目录内运行，只探查/读取当前工作目录（cwd）的文件。
3. 非结构自动选文件：优先 `unMesh.cas`+`unMesh.control`；否则找 `a.cas` 与同名 `a.control`
   配对；多于一组配对则报错列出；没有则报错。

## 模式检测矩阵（cwd，rank 间零通信）

所有 rank 共享 cwd，`inquire(file=...,exist=)` 在每 rank 结果确定性一致，**无需广播**。
检测必须在任何 `MPI_Comm_split` 之前完成（struct/uns 隶属由模式决定）。按此优先级：

1. **COUPLED**：`inquire('mix.control', exist=...)` 为真。所有耦合算例与 `grid_BC/`
   都含 `mix.control`，而两个单求解器都不读它 → 干净触发。
2. **UNS**：无 `mix.control` 时，按下方非结构选文件算法在 cwd 找到唯一
   `<stem>.cas` + `<stem>.control` 配对。
3. **STRUCT**：无 UN 配对，且 `inquire('control.ec')` 与 `inquire('Mesh3d.x')` 均为真。
4. **错误/歧义**：无 `mix.control` 但 struct 输入与 UN 配对**同时存在** → 报歧义
   （建议补 `mix.control` 或清理目录）；全部匹配不到 → 报错并打印各候选文件的存在性。

## 非结构选文件算法（纯 Fortran；Linux 项目）

1. 规范对：`inquire('unMesh.cas')` 且 `inquire('unMesh.control')` → 命中。
2. 否则泛扫描：`execute_command_line('ls *.cas > .mixnsolver_caslist 2>/dev/null')`，
   读入列表，每行去掉扩展名得 `<stem>`，`inquire('<stem>.control')` 校验配对。
3. 判定：恰 1 个配对 → 用之；>1 个不同配对 → 报错列出全部；0 个 → 报错。
4. 用完删除临时列表文件（`execute_command_line('rm -f .mixnsolver_caslist')`）。
   `cases/betchen/`（多配对）即天然报错用例。

## 驱动复用（方案：模块化抽取，最低回归风险）

耦合驱动逻辑 `struct_group_driver` / `uns_group_driver` / `convert_to_si` /
`convert_uns_to_struct_nd` / `write_couple_state` / `read_couple_state` 目前是
`src/main.f90` 的 `program mixsolver` 内部 `contains` 子程序（无法被外部调用）。
它们已自包含（`use` 各自依赖模块、接收 `comm`/文件名参数），可直接迁移。

- **新建 `src/coupling/mod_mix_driver.f90`**（`module mod_mix_driver`）：把上述子程序
  **原样搬迁**为模块过程，导出 `struct_group_driver(comm,cas,ctl,rank,nproc)`、
  `uns_group_driver(comm,cas,ctl,rank,nproc)`（签名不变，`parallel public`）。
- **`src/main.f90` 变薄调用者**：保留 MPI_Init、参数解析（逐字节不变）、
  `read_mix_control`、`MPI_Comm_split`，然后 `use mod_mix_driver` 后改调
  `call struct_group_driver(...)` / `call uns_group_driver(...)`；清空 `contains`。
  → `bin/mixsolver(_mpi)` 行为不变（纯代码移动），以 `cases/couple_porous` 回归验证。

## 三种模式的执行分支（新建 `src/main_dispatch.f90`，`program mixnsolver`）

```
program mixnsolver
   use mpi; use mod_precision
   MPI_Init; MPI_Comm_rank/size(COMM_WORLD, rank, nproc)
   call detect_mode_i(nproc, mode, ucas, uctl)   ! 每 rank 独立（第 1/3 节）
   select case(mode)
   case(MODE_COUPLED)
      n_struct_ranks = merge(0,1,nproc==1)
      colour = STRUCT_GROUP if rank<n_struct else UNS_GROUP
      MPI_Comm_split(COMM_WORLD, colour, key, sub)
      if (rank<n_struct) call mod_mix_driver::struct_group_driver(sub,'Mesh3d.x','control.ec',rank,nproc)
      else               call mod_mix_driver::uns_group_driver  (sub,'unMesh.cas','unMesh.control',rank,nproc)
   case(MODE_STRUCT)
      call mod_struct_driver::struct_solver_run(COMM_WORLD,'control.ec')   ! 新增
   case(MODE_UNS)
      call uns_solver_init(ucas, uctl, m,c,g,ctrl,bcs,fld, ier)
      loop: call uns_solver_step(m,c,g,ctrl,bcs,fld, ctrl%outer_max, ier)   ! 驱动到收敛
      barrier; if(rank==0) 用 mod_uns_output 写 VTU + 报告（镜像 main_uns.f90）
   end select
   call MPI_Finalize(ierr)
end program
```
- **struct-only**：新增 `struct_solver_run(comm, ctlfile)` 到 `src/structured/mod_struct_driver.f90`
  （**additive**，不改既有导出），内部镜像 `src/structured/main.f90` 的 `program main`
  时间环（step / filter / comput_force / output_Res / output_flow / time_average / 收尾），
  供 standalone 与 dispatch 共用；`struct_solver_init` 已会写 `Struct_Comm`。
- **uns-only**：复用 `mod_uns_driver::uns_solver_init/step`——**全 rank 各持全网格**
  （已确认免分区，数值==串行；与耦合侧 uns_group_driver 一致），rank0 出结果。
  **取舍备注**：-np N 是"每 rank 冗余全网格"，非真并行分区；若要真分区后续接入
  `main_uns_mpi` 的 METIS 栈（本节点不做）。
- **coupled**：完全复用 `mod_mix_driver`（与 mixsolver 同一套逻辑）。

## Makefile 改动

- `MAIN_DISPATCH_F := $(wildcard src/main_dispatch.f90)`；并从 LIB 的 `filter-out` 排除。
- `src/coupling/mod_mix_driver.f90` 由既有 `COUPL_F` 自动纳入 `LIB_OBJ_S/M`
  （耦合层编译已用 `$(MPIFC)`，ser 树即 np1 MPI 构建，无需新编译规则）。
- 新增 pristine 边：`mod_mix_driver.o` 依赖 `mod_struct_driver` / `mod_uns_driver` /
  `mod_coupling_exchange`；`(S|M)/main_dispatch.o` 依赖耦合层（仿 `main.o` 规则，用 `$(MPIFC)`）。
- 新增链接规则（仿 `bin/mixsolver`，用 `$(MPIFC)`）：
  ```
  DISPATCH_OBJ_S := $(LIB_OBJ_S) $(S)/main_dispatch.o
  DISPATCH_OBJ_M := $(LIB_OBJ_M) $(M)/main_dispatch.o
  bin/mixnsolver:     $(DISPATCH_OBJ_S); $(MPIFC) $(FFLAGS_S) -o $@ $(DISPATCH_OBJ_S) $(LDLIBS_S)
  bin/mixnsolver_mpi: $(DISPATCH_OBJ_M); $(MPIFC) $(FFLAGS_M) -o $@ $(DISPATCH_OBJ_M) $(LDLIBS_M)
  ```
- `.PHONY` 加 `mixnsolver mixnsolver_mpi`；`all:` / `mpi:` 增对应目标；
  `ALL_MAIN_S/M` 的 boot 边增 `main_dispatch`。

## 验证方案

1. **干净构建**：新增模块后首建 `make -j1 all` 与 `make -j1 mpi` RC=0；再 `make -j4` 增量。
2. **旧二进制不回归**：`make mpi` 重建 `mixsolver*` 后，
   `cases/couple_porous` np2 重跑与基线对拍（床梯度 −3021 Pa/m、界面连续 −0.08 Pa）
   ；`make structured structured_mpi unstructured unstructured_mpi coupling_test
   iface_law_test match_test units_test` 全部照旧 RC=0。
3. **struct-only**：临时目录拷 `grid_BC/` 的 `Mesh3d.x/control.ec/bc3d.inp`
   （去 `mix.control` 与 `unMesh.cas`）→ `./bin/mixnsolver` 产 `flow3d.dat`
   md5 == `dc134a2d196422043ecad7c86ac8f898`；`mpirun -np 2 bin/mixnsolver_mpi`
   与 np1 byte-identical。
4. **uns-only**：`cases/porous_plug`（`plug.cas`+`plug_darcy.control`）与
   `cases/betchen/bj_dae2.*` → `./bin/mixnsolver` 产出 VTU md5 ==
   `bin/uns_solver <cas> <ctl>`；`mpirun -np 2 bin/mixnsolver_mpi` 数值==串行。
5. **coupled**：`cases/couple_porous` 下 `mpirun -np 2 bin/mixnsolver_mpi`
   复现 `bin/mixsolver_mpi` 关键量；`cases/couple_channel` 冒烟。
6. **错误路径**：空目录 → 报错；`cases/betchen` → 报"多配对"；缺 `.control` → 报错。

## 要点文件
- `src/main.f90`（把内部子程序抽到模块，变薄调用者）
- `src/coupling/mod_mix_driver.f90`（新建，驱动抽取）
- `src/main_dispatch.f90`（新建，模式检测 + 分发主程序）
- `src/structured/mod_struct_driver.f90`（新增 `struct_solver_run`，additive）
- `Makefile`（新主程序/对象/链接/pristine 边）

## 交付后文档维护（project_rules）
同步 `docs/程序使用手册.tex` §2（新增 mixnsolver 用法子节）、`docs/93_changelog.tex`
顶部追加记录，并重编译手册；更新 memory-bank 三文件后自动 git 提交并推送。