# 当前上下文 (activeContext)

> 最后更新：2026-10-06（**C2 可压缩–多孔跨组界面完成并验证；两处网格/后处理缺陷已修**）

> **自检缺口①②③④已全部修 2026-10-06**：
> ① `Makefile` 让 `src/coupling/*` 与 `src/main.f90` 在 ser 树也用 `$(MPIFC)` 编译
> （先例=structured），并补 `main.o -> 耦合层` 的 pristine 边；`src/main.f90` 的
> `mod_uns_restart`/`write_field_dump`（uns 自动保存）用 `#ifdef HAVE_MPI` 守护。
> ② `bin/units_test` 链接 `$(FC)`→`$(MPIFC)`（同 `bin/mixsolver`），并补
> `test_units_exchange.o -> 耦合层` 的 pristine 边（耦合层 MPI-intrinsic）。
> ③ `src/unstructured/mod_uns_driver.f90` 的 uns restart 读取由整段停用改为
> `#ifdef HAVE_MPI` 守护：MPI 树恢复 `read_field_dump(serial=.true.)`
> （`couple_restart=1` 逐位续跑），串行树走 `#else` 告警且不崩。
> ④ 重建 M6-wing 位级基线并**持久化**到 `regress/m6wing/`（`flow3d.dat` md5
> `dc134a2d196422043ecad7c86ac8f898`，与历史记录逐位一致）+ `run_regression.sh`。
> 验证：`make -j4 all`/干净 `make -j1 all` RC=0；回归 mpi/structured(_mpi)/
> unstructured(_mpi)/units_test/match_test/coupling_test 全 RC=0；units_test 3/0 PASS、
> coupling_test np2 6/0 PASS、match_test 匹配 0 误差；ser `bin/mixsolver` np2 跑通
> couple_channel；`regress/m6wing/run_regression.sh` PASS。

> **CHECKPOINT 2026-10-06**：用户要求保存当前节点，随时可回退。
> 已完成：阶段 1–6 全链路 + 阶段 9 多孔全套 + 阶段 10 自动保存/重启 +
> 阶段 11 流 B（Dirichlet-Neumann 特征界面）+ C1（BJ 复用勾销）+ Turbo 清理。
> 下一步：① 可压缩–多孔跨组界面 ② 界面分派表 ③ 远期 NEU/SST/Liao。
> （2026-10-06 更新：① 已完成并验证，见下节 C2；下一步 = ② 界面分派表 ＋
> ③′ 新登记待办：uns 绝对压力 ≈−250 Pa 偏置机理（plan 阶段 11 第 3 条）。
> **仓库瘦身（2026-10-06）**：212 个运行产物 / 761.8 MB 脱离跟踪并重写历史，
> 跟踪 536 → 324 文件、`.git` 120 MB → **46 MB**；远端 `main` 已 force-with-lease
> 覆盖（零告警）。安全网 = `/home/sundong/mixsolver_pre_slim_backup/`
> （`repo_pre_slim.bundle` 98 MB 含完整旧历史 + `cases`/`grid_BC` 硬链接快照）；
> 工作树产物留在磁盘且被 `.gitignore` 忽略。规则见
> `.trae/rules/project_rules.md`「运行产物不入库」。
> **清单收口（同日）**：勾销 9 条已完成但未勾销的旧条目 + 删 1 条陈旧残留；
> 新登记 `cases/couple_channel/README.md` 缺失。
> **本节点已保存**：`docs/plan.md`、`memory-bank/*.md`、`.trae/rules/project_rules.md`
> （新增「验证/诊断硬规则」5 条）、`cases/couple_porous/*`（README ＋ 5 个脚本）、
> `cases/couple_channel/{gen_meshes.py, uns_full/gen_full.py}`。
> **2026-10-06 起本仓库已纳入 git**：`git init -b main`，首个提交 `97f247f`
> （= 本节点快照；536 文件、`.git` ≈119 MB、分支 main、工作区干净）；
> 单文件回退 `git checkout 97f247f -- <path>`，未提交改动 `git stash`。
> **工作流（2026-10-06 起，已固化进规则）**：每完成一个**小节点**即自动
> `git add -A && git commit -F -`（**无需询问**），与 memory-bank 更新同批提交，
> 一个节点一个提交；验证未通过/半成品不提交。规则：
> `.trae/rules/project_rules.md`「版本控制约定」。
> **远端（2026-10-06）**：`origin = git@github.com:cfddog/MixSolver.git`（SSH）；
> 首次 `git push -u origin main` 成功（594 对象 / 93.5 MiB 包体、无大文件告警），
> HEAD == `origin/main` == `8f40a31`；此后每节点「提交后自动 push」。
> 复验命令见下节 C2 与 `cases/couple_porous/README.md` §6。）

## 阶段 11（2026-10-06）：C2 可压缩–多孔跨组界面 ✅（完成并验证）
- **算例** `cases/couple_porous/`：struct 流体段（x=0–100 mm）↔ uns 多孔床
  （x=100–200 mm）弱耦合，界面 x=100 mm；对拍资产 = 单求解器全流域参考
  `uns_full/`（一个 uns：zone 2 fluid + zone 8 VC:porous，床 dP/dx 解析
  3.0217 Pa/mm、压降 302.1 Pa）。
- **求解器本体零改动**：两处缺陷都在网格生成/后处理——
  ① `uns_full/gen_full.py` 的 cell-id 排序把"x 串联"写成"z 并联"
  （`1+i+j*NX+k*NX*NY`，k 最慢 ⇒ zone 2/8 各成一条 z 全长薄片）：床梯度只到
  解析值 ~60%（1.78 vs 3.02 Pa/mm）、max|u| 42.87 > 入口、max p 398 Pa，且
  zone 互换后**逐位相同**（z 镜像 ⇒ 假象）；改为 `1+k+j*NZ+i*NY*NZ`（x 最慢）
  修复（zone2=1..2000 为 x<100、zone8=2001..4000 为 x>100 连续区间），旧网格
  留档 `scratch/unMesh.cas.zsplit_bug`。
  ② `compare_iface.py` 取窗"中心 ±1.25 mm"跨了两个 dx=2 mm 平面（把 x=101 mm
  的 295 Pa 与 x=103 mm 的 70 Pa 平均成 182.6 Pa）⇒ 幻影 112 Pa 跳变；改用
  绝对窗 `[100,102) mm` 后界面跳变 **−0.08 Pa**。
- **验收**：床梯度 −3021.3 vs 解析 −3021.7（0.01%，`check_bed_gradient.py`
  PASS；该脚本兼作回归守卫，在旧网格上主动报错）；界面连续（struct 面 295.0
  == uns 首排 295.0；struct 面速度 −1.37%）；耦合侧床梯度 −2959.2 Pa/m 与
  μu/(εK)+ρβu² 代界面速度 34.233 得 2959.2 逐位吻合 ⇒ −2.1% 残差全部来自
  struct 侧界面速度的 −1.37%（既有壁面层残余），非界面通量缺陷。
- **跨算例共性（开放项）**：uns 单求解器的**绝对压力水平**存在 ≈−250 Pa 内部
  偏置（内部整段常数下移 ＋ 紧邻入口边界那一层保持物理解水平 ＋ 出口末列回收）；
  与多孔无关（同网格全流体对照内部严格平 −251.17 Pa）、与出口 BC 类型无关
  （`outflow` 内部 −248.1，仅差 3 Pa，且末列速度被抬到 36.79 m/s）、严格 ∝u²
  （凹陷 −285.76 / u2 −71.44 / 2u −1143.06 Pa）。故跨求解器比**绝对压力**前需
  先修此项；C2 判据用"界面连续 ＋ 床梯度"。
- **文件**：`cases/couple_porous/{README.md(新建，全文取证), compare_iface.py
  (取窗 + 输出 %% 修正), gen_meshes.py(注释), uns_full/{gen_full.py,
  check_bed_gradient.py, plane_profile.py(新建)}}`、
  `cases/couple_channel/{gen_meshes.py, uns_full/gen_full.py}`（同款 cid 注释，
  单 cell zone ⇒ 无害）、`docs/plan.md`（C2 记录 + 勾销）。
- **下一步**：① 界面分派表（按两侧 cell_zone + 求解器类型自动分派，plan 阶段 11
  第 2 条）；② **已登记待办**：`uns` 绝对压力 ≈−250 Pa 偏置机理定位与修复（plan
  阶段 11 第 3 条；入手处 `src/unstructured/` 投影步与 BC 层；修复前跨求解器不得
  比绝对压力，见 `.trae/rules/project_rules.md`）。

- **用户明确三类界面域归属**：
  - 低速-多孔：两侧均在**非结构域**（uns 流体 + uns 多孔），内部面，
    不跨 solver、不经跨组交换层；
  - 可压缩-低速：struct 可压 ↔ uns 低速 = 流 B（已完成）；
  - 可压缩-多孔：struct 可压 ↔ uns 多孔（跨组；**2026-10-06 完成并验证**，
    见上节 C2 与 `cases/couple_porous/README.md`）。
- **结论**：C1 无新代码/无新算例，直接复用阶段 9 已实现并验证的内部
  BJ Robin 通量（mod_uns_simple momentum_assembly：切向
  C=μA(α/√K)/(1+α·d_Pf/√K)、法向两点扩散、交叉滞后等值反向）。
  验证资产 cases/beavers_joseph/：对拍双层 4×4 ODE，α=0 RMS 0.14%、
  α=1/2 RMS 1.09%/1.29%；np2 vs 串行 2.2e-7。
- **下一步**：可压缩（struct）–多孔（uns）跨组界面；之后界面分派表
  （按 cell_zone + solver 类型自动分派，plan 阶段 11 第 2 条）。
- 文档：docs/plan.md 阶段 11 节更新（C1 勾销 + 域归属 + 优先级表）。

---

## 阶段 11（2026-10-06）：流 B 交换层修法 ✅（couple_channel 界面连续性）
- **结论**：界面连续性严格在交换层内用两侧状态解决，不借 struct 入口。
  最终方案为 **Dirichlet-Neumann 特征界面**（4 版迭代定稿）：
  - struct 侧（亚声速出口）只接收 uns 背压 pb，ghost 用 boundary_Farfield
    亚声速出口同款线性化 Riemann 反射：db=d1+(pb−p1)/c1²，
    ub=u1+(p1−pb)/(d1·c1)·n_out，ghost=2·face−inner；法向量按
    Interface_List%face（1..6）从 Block 面法向表构造外法向
    （未用注册顺序法 Interface_List%normal）。
  - uns 侧只施加速度 Dirichlet；BC_INTERFACE 压力本就零梯度
    （bc_face_p: pf=pP），删除 set_interface_p 与压力重锚定。
  - 背压松弛 alpha=0.3·min(1,iter/iface_ramp)（ALPHA_MAX=0.3 硬编码
    main.f90）；alpha=1 在 Ma=0.1 下因 1/(ρc) 增益过大，ramp 后周期
    发散（iter50 NaN）；alpha<1 不改变收敛值（p1 每轮向 pb 松弛）。
- **被否方案**：v1 uns 删反射+struct ghost=2f−inner（子步前设）→
  周期-2 翻转 iter75 NaN；v2 同式但 10 子步后设 → 回声不动点流量
  衰减到 0.04；v3 ghost 直写 uns 首排格心（零阶保持）→ 稳定但
  冻结 1.7 kPa 压力跳变、质量 31.2 m/s。
- **最终结果**（run_flowb4c.log，iter1000，无 OverLimit，~100 iter
  冻结）：质量精确连续（struct 面 33.71 == uns 下游 33.71 m/s）；
  压力连续（struct 界面 342.5 Pa ≈ uns 首排 342，<1 Pa）；
  对拍 uns_full 单一求解器解（x=100：34.72 m/s、46.9 Pa）均值
  −2.9%、界面压力 +296 Pa、近壁型线最大偏差 11.8%——剩余偏差判定
  为两解器壁面层物理差异（struct 近壁 u=13.1 vs uns 9.05 @1.25mm），
  非交换层缺陷。
- **改动文件**：src/main.f90（删 uns 反射块/sp_bc 重锚定/set_interface_p；
  ghost 设置移到 10 子步后并传 alpha）、
  src/structured/mod_struct_driver.f90（struct_set_iface_bc 改特征背压，
  optional alpha）、src/structured/mod_struct_bc.f90（清理为 915 行
  干净版，0 Turbo；boundary_user 分派器 + boundary_user_Inlet 从
  external 原始版重建）、cases/couple_channel/compare_iface.py
  （修正 flow3d.dat 解析：output_flow 写含 ghost 的 (ni+1) 数组）。
- **异常备注**：bc.f90 会话中曾被外部回退为 1263 行含 Turbo 脏版
  （mtime 2026-10-06 12:12:38，无 git 无法追溯），已手工修复并
  M6-wing 位级回归（md5 dc134a2d196422043ecad7c86ac8f898）。
- **下一步（用户既定 B→A→C→D）**：流 B 已完成；后续低速-多孔 BJ
  复用内部 Robin（C1）；可压缩-多孔 D 暂缓。
- 已知遗留（既有）：`make all` serial 树 mod_uns_restart.f90 失败
  （无条件 use mod_uns_mpi_core，需 #ifdef HAVE_MPI 守护）。

---

## 阶段 11（2026-10-06）：弃用内流/叶轮机机制 + 流 B 重新定位 ✅（清理部分）
- **决策**：总压入口试验（couple_channel，IF_InnerFlow=1/P_In_Ratio=1.0123）
  np2 iter50 速度 78 m/s 发散、界面偏差依旧，证明界面不连续与入口无关。
  用户拍板：IF_InnerFlow 只保留外流部分；IF_TurboMachinary 叶轮机模式全部弃用。
- **已删源码**（src/structured/，5 文件）：
  - global：IF_TurboMachinary/IF_InnerFlow/P_In_Ratio/Turbo_Periodic_seta/Turbo_w；
    constants：BC_Wall_Turbo=201；
  - bc：Turbo 分派分支 + boundary_BC_Inflow_Turbo/Outflow_Turbo/wall_Turbo
    三个子程序，分派只剩 wall/Farfield/Inflow/Outflow/Symmetry/
    Extrapolate/>=900 user/BC_INTERFACE；
  - solver：旋转惯性力（离心+科氏）源项块；
  - init：Turbo_P0/T0/L0/w、Ref_medium_usrdef、namelist 键、默认值、
    init 叶轮机 else 分支、Ma=1/Re 反推介质块、bcast 打包/解包赋值
    （rpara 34/35/39、Ipara 27/29 槽位号保留不动）；
  - mpi：Umessage_Turbo_Periodic 子程序+调用；Coordinate/Mesh_Center_
    Periodic 删旋转角分支只留平移。
  - **教训**：`Aos`（mod_struct_init）就是侧滑角 AoS（Fortran 大小写不敏感），
    误删后编译报错已恢复；删除前 grep 必须把大小写变体都查一遍。
- **算例**：couple_channel/control.ec 恢复强制均匀入口（Kstep_save=5000），
  与 grid_BC/control.ec 同删 5 个废键（IF_TurboMachinary/Turbo_Periodic_seta/
  Ref_medium_usrdef/Turbo_w/Turbo_L0；gfortran 遇未知 namelist 键报错）；
  删 control.ec.forced、mix.control.ptot、run_ptot1.log。
- **验证**：ser/mpi 双树编译 0 error；M6-wing 对拍 /tmp/m6reg 原始基线
  （同 -O2 -std=legacy，t_end=0.501 恰 50 步；local-dt 下 tt=Kstep*dt_global
  精确成立），np1、np2 的 flow3d.dat（13,600,032 B）与 Step_mess.dat 均
  **字节一致**；couple_channel n_couple=2 冒烟通过（注意冒烟覆盖了
  unMesh_coupled.vtu，流 B 需重跑）。
- **回归环境备注**：M6 输入只有 Mesh3d.dat（formatted Plot3D），struct 侧
  迁移后文件名为 Mesh3d.x——新算例需 cp 一份；运行参数：
  `mpirun -np N bin/mixsolver_mpi mix.control Mesh3d.x control.ec unMesh.cas unMesh.control`
  （5 参数均可覆盖默认 grid_BC/ 路径；standalone struct 用 bin/struct_solver）。
- **流 B 已完成**：见顶部「流 B 交换层修法」节（Dirichlet-Neumann
  特征界面，2026-10-06）。参考解 cases/couple_channel/uns_full/unMesh.vtu
  （x=100 面 mean u=34.72、p gauge≈46.9 Pa），工具 compare_iface.py。
- 日志：build_{mpi,ser}_phase11.log、build_base_m6.log、m6_{probe,base,new}*.log、
  couple_smoke_phase11.log（均在项目根）。手册 docs/程序使用手册.tex 不存在 → N/A。

---

## 阶段 10（2026-10-06）：自动保存 + 联合重启 ✅
- **mix.control 新键**：`save_interval`（>0 = 每 N 轮耦合迭代保存）、
  `couple_restart`（1 = 联合重启）。解析/存取在 mod_reference_state.f90。
- **自动保存（两侧 iter 末对齐）**：
  - struct：`struct_solver_save`（mod_struct_driver）→ `output_flow`
    （flow3d.dat+Step_mess.dat）+ 新增 `output_flow_nodes`
    （mod_struct_io，格心原始量 d/u/v/w/T 平均到 Mesh3d.x 节点，
    Plot3D 函数文件 flow3d_node.dat，格式随 Mesh_File_Format）。
  - uns：uns root（全局 rank = n_struct_ranks）`write_field_dump`
    → unMesh_restart.dat；standalone（nproc==1）时 uns root 兼写
    couple_state.dat。
  - struct rank0 写 couple_state.dat（iter 号）。
- **联合重启**（couple_restart=1）：两侧读 couple_state.dat 得 iter0；
  struct_solver_init 新增 `force_restart` 可选参（置 Iflag_init=1 读
  flow3d.dat）；uns_solver_init 新增 `restart_file` 可选参，经
  `read_field_dump(..., serial=.true.)` 读入。**关键修复**：
  read_field_dump 加 serial 选项跳过尾部 MPI_Bcast——耦合驱动不调
  mpi_bootstrap（mpi_comm=COMM_WORLD 默认），否则 struct rank 被卷入
  Bcast 挂死。循环从 iter0+1 继续，iface_ramp 计数连续。
- **验证（grid_BC np2，n_couple=6、save_interval=3）**：连续 6 轮 vs
  3 轮+重启 3 轮，flow3d.dat / unMesh_restart.dat 均 **byte-identical**
  （check_restart.py）；flow3d_node.dat 维度与 Mesh3d.x 一致。
  ser/mpi 双树编译无新增 error。归档 grid_BC/cont6/、restart33/。
- **已知限制**：dump v1 不含 T_s（LTNE 重启 T_s 回 init 值）；
  struct Kstep/tt 不恢复（仅影响输出文件名）。
- **涉及源码**：mod_reference_state.f90、mod_uns_restart.f90（serial）、
  mod_uns_driver.f90（restart_file）、mod_struct_driver.f90
  （force_restart + struct_solver_save）、mod_struct_io.f90
  （output_flow_nodes）、main.f90（两驱动 save/restart 分支 +
  read/write_couple_state）。
- 手册 docs/程序使用手册.tex 不存在 → N/A。
- **下一步候选**：阶段 11 多类型界面（低速-多孔 BJ / 低速可压缩 /
  可压缩-多孔，按块属性分派）；或挂起的 couple_channel 压力基准。

---

## 阶段 9 收尾（2026-10-06）：α/K 标定实验设计 ✅
- **交付**：cases/calibration/ 数值虚拟标定框架（4 Python 文件 + README
  设计报告 + images/calib_closure.png）。反演器与数据源解耦（forward 可
  换求解器 VTU 采样），正演构型复用 porous_plug/porous_disp/ltne_disp/
  beavers_joseph 四算例的半解析参考。
- **分层解耦设计**：P1 多流速压降→K+C_F；P2 混合层横向剖面→α_t；
  P3 LTNE 发汗冷却两相测温→α_l+h_sf 联合；P4 界面 u(y) 剖面→α_BJ。
  解耦依据：等温/非等温分离（P1 vs P2/P3）、流向 vs 横向 Pe 分离
  （α_l vs α_t）、界面通量独立（α_BJ）。
- **关键设计发现（仪器选型级）**：TC 噪声按绝对温度缩放（0.2%×305K≈0.6K），
  混合层 α_t 标定 ΔT=10 K 时 CRB 12.7% 不可用，**必须 ΔT≥50 K**（CRB 2.2%）。
- **闭环验证**：四子问题从 ~2× 离真值初值 4–6 步 LM 收敛，恢复误差
  K 0.19%/C_F 1.1%/α_t 3.8%/α_l 1.1%/h_sf 0.25%/α_BJ 1.8%，全部在
  CRB 界内；多参问题相关 0.65–0.76、条件数 ≤12.5，无病态。
- **CRB 选型表**：1% 压差→K 0.6%/C_F 5.5%；0.5% 测温→α_l 1.5%/h_sf 0.9%；
  2% PIV→α_BJ 5%。
- 无源码改动，无需回归。手册 docs/程序使用手册.tex 不存在 → N/A。
- **阶段 9 至此全部完成**。下一步候选：阶段 10（流场自动保存+耦合联合
  重启）或阶段 11 界面主线（含挂起的 couple_channel 压力基准、struct
  IF_InnerFlow=1 总压入口候选）。

---

## 阶段 9 收尾（2026-10-05）：LTNE + 热弥散组合验证 ✅
- **两 bug 修复**：①Bear 纵向弥散缺 /|u|（mod_uns_simple 温度装配，
  原结果放大 |u| 倍；横向形式不变，porous_disp/porous_graetz 四工况
  位级回归一致）；②MPI gather 缺 temperature_solid
  （mod_uns_gather 补 T_s Gatherv；main_uns_mpi rank0 fld_g 补
  setup_porous_fields 过输出 guard）。
- **1D cases/ltne_disp/**（恒热流发汗冷却，disp_l=0/5mm）：
  4 阶耦合 ODE 半解析（三次特征根，exp(r(x−L)) 防溢出，Nield 微观分配）；
  Tf 中位 0.4%/Ts≤0.83%，disp1 出口温升 973.46 vs 974.41
  （弥散回流到 Dirichlet 入口，aΔTf=q−Kf Tf'(0)）；8095/11219 步；
  np2 两相 rel 1e-7。
- **2D cases/ltne_graetz/**（恒壁温半通道，ε=0.4、k_s=1.0、Bi=0.845、
  Λs=0.101，disp_t=0/2.5e-4）：参考=双相二次特征值问题（4m 一阶块
  y'=My，2m 纯实衰减模态，入口 θf=1+χs=0；两相轴向导热保留、
  流体横/纵 Pe 分开；Bi→0 退经典 Graetz 2.457）；180×40 剖面 max
  1.4%/3.3%（仅 x=20mm 入口段），360×40 x 加密减半 0.7%/1.8%，
  x≥50mm ≤0.5%，bulk 全程重合；7018/6939 步；np2 两相 rel 2e-7。
  **关键**：固相入口不是 θs=1——固相无对流且入口绝热，θs(0) 是耦合
  模态给出的横向平衡剖面；忽略固相轴向导热（Λs=0.10）会使参考解偏差
  数%。两算例均已归档 README+图+脚本+日志。
- **阶段 9 剩余**：仅 α/K 标定实验设计（远期实验项）。
- 教训：hex→VTU 输出为每 hex 12 tetra（`4 n1..n4` 记录），解析必须剥
  首 token + ix 圆整分层（plot_bj.py 已正确，无需改）。
- 手册 docs/程序使用手册.tex 不存在 → N/A。
- 挂起主线：couple_channel 压力基准（struct IF_InnerFlow=1 总压入口候选）。

---

## 阶段 9 收尾（2026-10-05）：Beavers-Joseph 界面条件 ✅
- **实现**：`cell_zone` 行新增 `bj_alpha=α`（默认 0=禁用；Robin 型界面通量
  C=μA(α/√K)/(1+α·d_Pf/√K)，切向用 C、法向保 D、交叉分量滞后 rhs，
  两侧等值反向保动量守恒）；新增全局 `body_force = fx fy fz`（N/m³）。
  涉及 mod_uns_control（解析）、mod_uns_fields（bj_alpha 映射）、
  mod_uns_simple（通量+体力源项）；MPI 复用 momentum_assembly 零改动。
- **关键物理发现**：离散 BJ 通量的连续极限是"应力连续+速度跳变"双层模型
  （串联阻力 1/C=d_Pf/(μA)+λ/(μαA)），**不是**经典 BJ 滑移公式——多孔侧
  Brinkman 项激活（ε=1 时 μ_eff=μ）存在可分辨 Brinkman 层。正确参考为
  4×4 双层耦合 ODE 系统（plot_bj.py `reference()`）。
- **验证**（cases/beavers_joseph/，开口槽道 100×40、K=1e-6、λ/dy=4、
  体力驱动）：α=0 RMS 0.14%（纯离散误差）；α=1/2 RMS 1.09%/1.29%；
  经典 BJ 滑移值明显偏离（反证双层模型正确）。MPI np2 vs 串行最大相对
  偏差 2.2e-7。回归：cavity 697 步、porous_Ra10 748 步位级不变。
- **教训**：①闭盒体力被压力梯度抵消，BJ 验证必须开口槽道（两端
  pressure-outlet）；②本算例被动温度漂移使 dT_max 门槛不收敛，取
  outer_max=800 截取已收敛速度场（boussinesq=.false. 不影响 u）。
- **阶段 9 剩余**：α/K 标定实验设计（LTNE+弥散组合验证已于当日完成，见顶部）。B-J 机制供阶段 11
  低速-多孔界面复用。手册 docs/程序使用手册.tex 不存在 → N/A。
- 挂起主线：couple_channel 压力基准（struct IF_InnerFlow=1 总压入口候选）。

---

## 阶段 9 进展（2026-10-05）
- **第一批**：VC 标记自动识别；porous-plug Darcy/Forch 压降 0.000%/0.002%。
- **第二批**：各向异性渗透率 perm_xx/yy/zz（0.085%/0.122%）；热弥散 Bear 张量
  （混合层 erf ≤1%）；修复 velocity-inlet 静温 token bug。
- **第三批（本轮）**：二维多孔介质**强制**对流 cases/porous_graetz/（半平行板
  恒壁温 Graetz 热入口，K=1e-10、k_s=k_f）：渐近 Nu_a=π²/4，纯导热 0.42%、
  带横向弥散 1.38%；两工况 θ_b(x*) 坍并；np2=串行（6238 步）。
  自然对流多孔方腔早已在 cases/natconv/ 完成。
- 回归：Ra10 748 步/Nu=1.078 不变。
- **阶段 9 剩余**：Beavers-Joseph（依赖阶段 11 界面）；LTNE+弥散组合验证；
  α/K 参数标定实验设计。混合对流（穿流+浮升力）未做，需要时再排期。
- 关键细节：①Graetz 对均匀到壁解析解须 √K/dy≤0.08（近壁无弥散层效应）；
  ②LTE 对拍解析时 k_s 设 =k_f 使 k_eff=k_f；③Nu_a 用半高 a，渐近 π²/4；
  ④erf 相似变量用 κ=Γy/(ρcp)；⑤CAS (45 VC 标记须无空格。
- 挂起主线：couple_channel 压力基准（struct IF_InnerFlow=1 总压入口候选）。

---

## 文档维护（2026-10-05）
- 程序功能说明 L8–15 已录入 docs/plan.md「后续计划」= 阶段 9（多孔介质完善：
  VC 自动识别/porous-plug/B-J 等）、10（流场自动保存+耦合联合重启+Plot3D 节点插值）、
  11（低速-多孔/低速可压缩/可压缩-多孔界面按块属性分派）、12（Gambit NEU 输入）、
  13（SST 暂缓修复；Liao 格心型 FD 远期）。优先级表同步扩展。progress.md 已同步。
- 注：当前主线调试的 couple_channel（低速 fluid-fluid 界面）即阶段 11 低速可压缩界面的基础。

---

## 当前任务进展（2026-10-05：inlet_ramp + 界面压力/质量修正；grid_BC 判为物理不兼容）

**阶段 7 进行中。inlet_ramp 机制已实现（serial+MPI），耦合发散根因已查明并打两个补丁，但 grid_BC 物理本身不兼容弱耦合。**

本轮改动：

1. **inlet_ramp（代码就绪，耦合算例禁用）**：`ctrl_t%inlet_ramp`（默认 1）；
   mod_uns_bc 模块级 `g_inlet_ramp_factor` + `set_inlet_ramp_factor`；
   MASSINLET 面速乘 ramp；simple_run / simple_run_mpi 每外迭代设置。
2. **main.f90 uns_group_driver 两项耦合稳定性修正**：
   - 压力重锚定：`sp_bc -= mean(sp_bc) − mean(ip)`（struct 绝对压 ~1e5 Pa
     不可直接施加到 uns 不可压参考系；保留梯度、锚定均值）
   - 界面质量守恒：非界面边界用 bc_face_vel 求 F_other（含 ramp），
     界面面速加均匀法向修正 δu_n = −(F_other+F_iface)/A_iface，
     使总边界通量为零（全 Dirichlet 速度边界下纯 Neumann PPE 的相容条件）
3. **诊断结论**：
   - standalone uns 在 grid_BC 发散 = 预期（封闭域无出流，非 bug）
   - 修正后 it=1 mass-imbal 3.9e-10、du_max=151 有界，但 ~it=50 仍 ICC0
     非正主元 → NaN
   - **判定 grid_BC 物理不兼容**：struct Ma=3 可压超音速 vs uns 不可压
     1.7 m/s 腔体，Dirichlet 界面弱耦合无法收敛；inlet_ramp 与驱动级
     质量修正互斥（耦合算例须 inlet_ramp=1，unMesh.control 已注释）

**下一步候选**（待用户拍板）：
- ① 改双向物性兼容算例验证耦合（如两侧均不可压的腔体-腔体）
- ② 界面改通量型/混合 BC（而非 Dirichlet 状态）
- ③ 阶段 7 其余 TODO：uns→struct 面积加权（替全场平均占位）、struct→uns
  peer_w 双线性、extract 的 rho=ctrl%rho、界面收敛检查

---

**阶段6 核心收尾完成：struct 真实求解器接入驱动 + COMM_WORLD 分组隔离 + 双向界面交换全链路跑通。**

新增/修改文件（本轮）：

1. **`src/structured/mod_struct_driver.f90`**（新增）：struct 求解器封装
   - `struct_solver_init(comm, ctlfile)`：设 my_id/Total_proc（组内 rank），
     `Struct_Comm = comm`，附 Bsend buffer，symlink control.ec，
     read_parameter→Init→set_control_para→check_mesh_quality→Init_flow
   - `struct_solver_step()`：NS_Time_advance(1)
   - `struct_extract_iface(rho,u,v,w,T,p,nfaces)`：守恒量转原始量
     （T = p·γ·Ma²/ρ）
   - `struct_set_iface_bc(rho_nd,...,T_nd,nfaces)`：原始量转守恒量直写
     ghost cell（p=ρT/(γMa²)，E=p/(γ-1)+½ρ|u|²）

2. **`src/coupling/mod_coupling_exchange.f90`**（新增）：跨组 MPI 交换
   - count+payload 两段协议（支持 nfaces=0 空缓冲，两侧面数可不同）
   - exchange_struct_to_uns / exchange_uns_to_struct / recv_uns_iface_state
   - tags 110/111 (s→u)、210/211 (u→s)

3. **`src/common/mod_interface.f90`**：Interface_FACE_TYPE 增
   ic/jc/kc（内层 cell）+ ig/jg/kg（ghost cell）索引

4. **`src/structured/mod_struct_grid.f90`**：register_bc_interfaces 按
   Bc%face 分派填充 inner/ghost 索引

5. **`src/structured/mod_struct_global.f90`**：Global_Var 增
   `Struct_Comm = MPI_COMM_WORLD`（默认值保证 standalone 行为不变）

6. **6 个 struct 库文件批量替换** `MPI_COMM_WORLD → Struct_Comm`：
   mod_struct_fdm(2)/init(8)/solver(2)/io(37)/mpi(15) —— **关键修复**：
   替换前 np=2 时 struct 侧 collectives 在 COMM_WORLD 上等待 uns rank
   参与 → init 挂起；替换后各组只在自己的子通信域内同步。
   （structured/main.f90 的 Init_mpi 保留 COMM_WORLD，standalone 不受影响。）

7. **`src/main.f90`**：
   - struct 组：init(comm,ctl) → 循环{extract→exchange s→u→recv u→
     convert_uns_to_struct_nd（简单平均占位）→set_iface_bc→step}
   - uns 组：`has_struct = (nproc>1)` 包住交换+设BC（**nproc==1 修复**：
     此前 uns rank0 与自己 Recv 死锁）
   - nproc==1 时 n_struct_ranks=0（无 struct 组）

**验证（grid_BC，cwd=grid_BC 因 struct 硬编码文件名）**：
- `make mpi` / `make structured_mpi` 0 error
- np=2 耦合冒烟：3 轮耦合迭代 exit=0，双方向各交换 250 面，
  struct step 正常（`grid_BC/grid_BC_smoke_np2.log`）
- np=1（无 struct 组）：不再挂起，SIMPLE 正常启动
- standalone struct 回归：残差历程与改动前一致（Struct_Comm 默认值
  COMM_WORLD 生效）
- 已知遗留：uns 侧 `iface |u| NaN` 为 unMesh.cas mm 单位发散（旧问题，
  非本次耦合引入）

**遗留 TODO**：
- uns→struct 现为全场简单平均占位；应接 mod_interface_exchange 面积加权
- struct→uns 逐面直接设置，未用 peer_w 双线性插值
- uns_solver_extract_iface 的 rho 用 1.0 占位（应 ctrl%rho）
- n_couple=3 / n_uns_steps=50 硬编码（应从 mix.control 读）
- 界面收敛检查未实现
- unMesh.cas mm→m 单位缩放
- **grid_BC 真实耦合物理验证**（struct Ma=3 外流 + uns 低速内流）

---

## 历史：阶段6 弱耦合驱动框架（2026-10-05）

新增/修改文件：

1. **`src/main.f90`**（新增）：混合耦合驱动主程序
   - MPI_Init → 按 rank 拆 STRUCT_GROUP(0) / UNS_GROUP(1) 子通信域
     （`MPI_Comm_split`，默认 rank0=struct，其余=uns）
   - 读 `mix.control`（`read_mix_control`）
   - struct 组：占位驱动（循环 barrier 同步，真实求解器接入留后续）
   - uns 组：`uns_solver_init` → 设 interface BC（placeholder 均匀入流）
     → `uns_solver_step(n)` → `uns_solver_extract_iface` → barrier 同步
   - 固定 n_couple=3, n_uns_steps=50（TODO: 从 mix.control 读）

2. **`src/unstructured/mod_uns_driver.f90`**（新增）：非结构求解器驱动封装
   - `uns_solver_init(casfile, ctlfile, m,c,g,ctrl,bcs,fld, ier)`：镜像 main_uns 初始化
   - `uns_solver_step(..., nsteps, ier)`：包装 simple_run + nsteps 上限
   - `uns_solver_extract_iface(m,g,bcs,fld, faces, rho,u,T,p, ier)`：提取
     BC_INTERFACE 面的属主格 SI 状态（rho 用 1.0_dp 占位，TODO: 用 ctrl%rho）

3. **`src/unstructured/mod_uns_control.f90`**：新增 `BC_INTERFACE=9`，导出，
   bc_type_name → `'coupling-interface'`

4. **`src/unstructured/mod_uns_bc.f90`**：
   - `bc_t` 增 `iface_vel(3,nfaces)` / `iface_p(nfaces)` 数组
   - `build_bc`：interface zone 面归入 BC_INTERFACE 组（原 fgrp=0 跳过），
     分配 iface_vel/iface_p 零初始化
   - `bc_face_vel` / `bc_face_p`：BC_INTERFACE 分支读 iface_vel/iface_p
   - 新增 `set_interface_vel(bcs, faces, u(3,:))` / `set_interface_p(bcs, faces, p(:))`

5. **`src/unstructured/mod_uns_simple.f90`**：`simple_run` 增 optional `nsteps`
   参数，循环上限 `it_max = min(ctrl%outer_max, nsteps)`

6. **`Makefile`**：
   - `BOOTSTRAP_ORDER` 扩展为 mod_precision → mod_constants →
     mod_reference_state → mod_interface（保证 clean build 顺序）
   - 非结构层新增 `UNS_LAYER6_S = mod_uns_driver` / `UNS_LAYER9_M = mod_uns_driver`
     及对应 pristine 依赖边

**验证**：
- `make mpi` 编译通过，`bin/mixsolver_mpi` 生成
- cavity (rb_Ra1000) + 2 ranks：3 次耦合迭代完整完成，exit=0，
  struct/uns 同步正确，SIMPLE 收敛正常（mass-imbal 1e-13）
- unMesh.cas：interface zone 识别正常（250 面 → BC_INTERFACE 组），
  但 mm 单位导致数值发散（已知遗留问题，非框架问题）
- cavity 回归：697 步收敛，非结构求解器无破坏

**已知遗留 / TODO**：
- struct 侧真实求解器接入（init/step/extract/set BC）目前是 barrier 占位
- MPI 跨组数据交换（struct↔uns Send/Recv）未实现，当前用 barrier 占位
- `uns_solver_extract_iface` 的 rho 用 1.0_dp 占位，应改用 ctrl%rho
- n_couple / n_uns_steps 硬编码，应从 mix.control 读
- unMesh.cas 的 mm→m 单位缩放待 units 层统一

**下一步**：
- 接入结构侧真实求解器（提取界面顶点状态 + 设置界面 BC ghost cell）
- 实现 MPI 跨组数据交换包装（基于 phase-5 mod_interface_exchange）
- 界面收敛检查
- 修复 rho 占位 + 从 mix.control 读耦合参数

---

## 历史：槽道验证 + outflow 边界（2026-10-05）
**不可压缩槽道（plane Poiseuille）验证已归档：`cases/channel/`**。
- **主题① mdot-inlet 启动稳定性**：Re100 从零场 it=1 du_max≈1.8U、mass-imbal
  5e-13，无发散单调收敛；Re1000（10× 冲量）1033 步稳定。上会话"纯流体
  mdot-inlet 启动发散"确认为 LTNE 无壁面网格特异问题（已修），非 mdot-inlet 本身。
- **主题② pressure-outlet 出口段伪调整诊断**：末列 Ucl 1.5→1.71、vmax 0.115U
  向中心汇聚；对照实验排除入口 BC/对流格式/法向/p′ 路径/并行，定性为同位网格
  +Rhie-Chow 在固定压力边界的离散不动点。
- **主题③ 新增 `outflow` 边界（BC_OUTFLOW=8，用户拍板）**：法向零梯度（u/p/T
  外推）+ 全局质量缩放（flux_rhiechow 后缩放 outflow 面通量使总出流=总入流；
  POUTLET/FARFIELD 压力 Dirichlet 面自调故排除在 m_req 外）。PPE 为纯 Neumann
  走 cell 1 销钉，缩放保证 Σrhs=0 相容。实现：mod_uns_control（解析）+
  mod_uns_bc（bc_face_vel/T）+ mod_uns_simple（outflow_mass_sums/scale/rescale
  + 三驱动接入 + momentum/temperature 装配分支）+ mod_uns_simple_mpi
  （outflow_mass_rescale_mpi，3 元组 allreduce）。
- **结果（Re100 出口末列）**：Ucl 1.7091→1.6420、vmax 0.1152→0.0402、剖面
  L1 误差 11.45%→5.63%、807→790 步；内部场与 pressure-outlet 逐 7 位一致；
  np2 791 步出口值一致。**已知限制**：末列 +3.5% 体速度外观（=1/β 不动点，
  面通量精确守恒）与 ~9% 中心线凸起（一阶零梯度固有）；工程惯例出口远离
  关注区 ~10H。曾试动量用滞后缩放通量（更差，已回退）。
- **回归四算例 PASS**：cavity（518 步/Ghia）、cylinder far-field（Cp 差 0.0031）、
  porous_Ra10（748 步/Nu 一致）、LTNE（3344 步/全指标一致）。
- **Re1000 对照**：pressure-outlet 1027 步 vs outflow 1033 步；U(50.5)/U 同为
  0.9997；末列 u_max/U 1.691 vs 1.574、L1 14.74% vs 5.36%（伪调整机制同 Re100）。
  （旧后台任务随会话死亡、日志 0 字节，已重跑；README §3.4 已回填。）
- 手册 docs/程序使用手册.tex 不存在 → N/A（outflow/mass-flow-inlet 等新参数
  待手册创建后补录）。
- **下一步候选**：① 阶段 5（界面数据交换+量纲统一，主线 pending）；② 入口
  ramp 渐启（inlet_ramp=N，上会话建议，本任务证明槽道非必需）；③ LTNE 瞬态验证。

---

## 历史：1D LTNE 发汗冷却验证完成（2026-10-04）
**1D LTNE 发汗冷却算例已完成并归档：`<项目>/cases/ltne/`**（2026-10-04）。
- 算例：空气 2 kg/(m²·s)、300 K 强制通过多孔不锈钢板（ε=0.3, SS304），
  出口端面恒热流 q″=2e6 W/m²（压力出口 zone 兼作热流壁，`tbc = 6 2 2.0e6`）；
  网格 LNTE1D.cas（Pointwise, 12160 wedge，用户源文件 /mnt/c/temp/LNTE1D.cas）；
  准 SI 约定同 grid_BC。h·a=0.3 W/(m³·K) 使交换长度 6.15≈4.9dx 可分辨。
- **结果 vs 解析解**：能量守恒 ΔT_f(L)=995.02 K vs 995.02 K（−0.001%）；
  T_s 剖面偏差 max 0.074%；T_f 中位 0.198%（最大 9.7% 在出口陡峭段，dx 限制）；
  Darcy 压降 5.17e5 vs 理论 5.23e5 Pa；3343 次外迭代收敛，mass-imbal=1.8e-19。
- 归档：README.md（完整测试报告）、ltne1d.control、ltne1d.log、LNTE1D.vtu、
  plot_ltne.py、images/ltne1d_profiles.png（三联图 T_f/T_s/θ vs 解析解）。
- **本次修复 4 个真 bug**（调试中暴露，详见 worklog 2026-10-04 LTNE 条）：
  1. **PPE 钉死 cell 1 = 假质量汇**（核心）：有压力 Dirichlet 面（POUTLET/
     FARFIELD）时矩阵已非奇异仍钉死 cell 1 → 其连续性行被丢弃成永久源/汇。
     本算例 cell 1 恰在入口层：60% 流量"消失"、内部 u=0.4·u_in、imbal 卡死
     4.4e-2、该 cell T 爆炸。修复：仅全 Neumann 封闭域才钉死（串行 ppe_assembly
     + MPI ppe_assembly_mpi 带 has_pdir 的 MPI_LOR allreduce）。
  2. **POUTLET 倒流削弱动量对角**：ap += F 改为 max(F,0) 隐式 + min(F,0)·u_P
     显式（违反迎风、可致 ICC0 非正主元 NaN）。
  3. **flux_ref 归一化漏 massinlet uspeed**（6 处 uscale 补 gb%uspeed）。
  4. **MPI 本地网格丢 cell-zone**：mod_uns_local_mesh 未拷 czone/cztype/czt →
     **MPI 多孔/LTNE 物理整体缺失**（此前 porous np2"收敛"实为无 Darcy 的纯流体
     解，Ra10 np2 296 步 vs 串行 748 步即铁证）。修复 build_local_mesh 拷贝
     三数组后 np2 与串行逐迭代一致（du_max=1.699 等全程吻合）。
     另：build_bc 增 allow_empty 可选参（METIS 把小 zone 全分给单 rank 时，
     其他 rank 本地 build_bc 不再报 "no boundary faces"）。
- **回归**：porous_Ra10 串行位级一致（两次）；cavity np2 正常；LTNE np2 与串行
  迭代历程一致。**已知限制**：纯流体（无多孔阻力）mdot-inlet+出口在该网格
  启动期发散（it=1 速度修正 V/apc 放大 ~1e4，apc 仅粘性量级）——后续可做
  入口渐启/瞬态启动；多孔/含阻力算例受 Darcy 对角保护不受影响。
- natconv 归档（cases/natconv/）与归档约定不变；手册 docs/程序使用手册.tex
  仍不存在，N/A。
- **下一步候选**：① 阶段 5（界面数据交换+量纲统一，主线上一直 pending）；
  ② mass-flow-inlet 启动稳健性（BC 渐启）；③ LTNE 瞬态验证（init_t_s 非平衡
  初值 + 弛豫时间常数 τ 对比）。

---

## 历史：算例归档 + 测试报告（2026-10-04，cases/natconv）
- **用户约定（长期）**：以后每个完成的验证算例都归档到 `<项目>/cases/<主题>/` 并写
  README 测试报告+结果图；约定全文在 trae 项目记忆 project_memory.md。
- **沙箱注意**：cases 迁入项目后位于可写白名单内，求解器可就地读写、
  后台任务亦可；/tmp 会话间会被清空，仅作临时工作区。
- **Darcy 侧壁加热方腔基准**（左 301/右 300/上下绝热，eps=1，Da=1e-6）：
  Ra_K=10 Nu=1.078（文献 1.07）；Ra_K=100 Nu=3.096/3.097（文献 3.10）；
  Ra_K=1000 Nu=13.405(160²) 外推 13.45 vs Mahmud&Fraser 13.64（−1.4%）；
  热平衡 0.004%。
- **附加路径验证**：① 次临界 Ra_K=10<4π² → Nu=0.996 纯导热；② 超临界 Ra_K=100
  → Nu=2.13 成环；③ LTE k_eff →3.407 vs 3.4；④ 离散等价检验差 1e-5；
  ⑤ Forchheimer → Nu 3.096→2.76。

---

## 阶段 4 已完成（2026-10-04）：交界面几何匹配 + 插值权重
- 扩展 `Interface_FACE_TYPE`（mod_interface.f90）：增 `nv`、`verts(3,nv)`（面顶点坐标）、
  `peer_w(:)`（对侧顶点插值权重）。
- 结构侧 `register_bc_interfaces`（mod_struct_grid.f90）改为**逐单元面登记**：
  每个 interface 单元面 1 条条目、4 顶点 CCW、centroid/normal/area/bbox 逐面算；
  报告按 (block,face) 聚合。grid_BC block2 face2 由原 1 条→250 条。
- 非结构侧 `register_interface_zones`（mod_uns_geometry.f90）每条已有面补 `verts`（来自 m%f%nodes）。
- 新建 `src/coupling/mod_interface_match.f90`：`match_interfaces(tol=1e-3)`
  - 拆分 Interface_List 为 struct/uns 索引表
  - 对每个 uns 面质心：投影到各 struct quad 平面 → Newton 反演双线性 (u,v) → 取落在
    [0,1]² 内且法向距离最小者；无包含时回退最近 struct 质心
  - 存 4 顶点双线性权重 w=[(1-u)(1-v), u(1-v), uv, (1-u)v]（和恒为 1）
  - 双向 peer_id（uns→struct 带权重；struct→uns 反向链接，peer_w 留空待阶段 5 面积加权平均）
  - 报告：uns 匹配率、struct 被引用率、max/mean 投影距离、max|Σw−1|
- 新建 `src/coupling/test_match.f90` + Makefile `match_test` 目标：链接 structured(MPI)+
  unstructured(MPI)+coupling，读两侧真实网格→登记→匹配。
- **grid_BC 验证结果**：250/250 uns 面匹配(100%)、250/250 struct 面被引用(100%)、
  max 法向投影距离 1.03e-38（机器零，两侧共面）、max|Σw−1|=0。
- **回归**：cavity 串行位级一致；cylinder(vinlet) 位级一致；cavity np2 仅 epsilon 噪声
  （interface zone 为空时登记例程早退，新字段未触发）。

## 阶段 3 步骤 C 已完成（2026-10-03，真实网格 grid_BC 验证）：
- 用户提供 `grid_BC/`：Mesh3d.x（结构，unformatted Plot3D，3 块 43/26/34×101×11，单位 mm）、
  bc3d.inp（blk-1-split-2 的 j=1 面 generic:8）、unMesh.cas（Pointwise，5577 节点/9220 wedge/
  24422 面，单位 mm；体区域 zone 2 fluid；面区域 4 symmetry/5 wall/6 pressure-inlet/7 interface）。
- 体区域解析：`mesh_t` 新增 `czone(ncells)/nczone/czt/cztype`；cas_reader 的 sec_cells 存每单元体区域号，
  finalize_zones 把 m%zone 按"面/体引用"拆成面区域表与体区域表；新增 grow_iarray（顺带修复 ctype 多记录扩容）。
- `resolve_cell_zones(m,ctrl)`（mod_uns_control）：id/name 匹配，未列体区域默认 fluid，
  cztype 标 1/2，未知 zone 报错；主程序 read_control 后调用（MPI 各 rank 都算、仅 rank0 打印）。
- build_bc 豁免：interface zone 边界面不再强制要求 bc（fgrp=0，待阶段 5 交换），判定规则与登记一致。
- 新 BC：`mass-flow-inlet <mdot>`（BC_MASSINLET=7），uf=-(mdot/rho)·n_out，走定速度装配/PPE 零梯度压力；
  动量、温度、报告均已接入。
- 结构化源码 `Mesh3d.dat`→`Mesh3d.x`（4 文件 12 处）；control.ec 用 `Mesh_File_Format=0`。
- **界面几何两侧精确重合**：均为 y=0 平面、x∈[420,670]、z∈[0,50]、质心(545,0,25)、面积 12500 mm²、
  250 个四边形；法向相反（struct −y / uns +y），物理正确。阶段 4 可直接做点对匹配。
- 回归：cavity 串行位级一致；cylinder(vinlet) 与原始二进制位级一致；cavity np2 仅 epsilon 噪声。
- grid_BC/control.ec：Ma=3 / Re=1e4/mm / Tinf=108K / Tw=300K / AoA=0 / 层流（用户已确认）；
  grid_BC/unMesh.control：300K 空气 rho=1.177/mu=1.846e-5，mdot=2 kg/m²/s（用户已确认）。
- 已知现象：uns 独立运行在真实网格上发散（interface 面无出流→质量无出口→PPE 奇异 NaN），
  阶段 5 交换实现前不可独立跑物理场，属预期。
- 手册 N/A：MixNSSolver docs/ 下暂无 程序使用手册.tex（仅 plan.md、程序功能说明.md）。

## pressure-far-field BC 零解 bug 修复（2026-10-03，步骤 C 之前）
- 根因：`mod_uns_bc.f90` 的 `bc_face_vel` 中 BC_FARFIELD 用内场速度 `uP·n` 判定入流/出流；
  零初场 uP=0 → un=0 → 上游入流面被误判为出流（uf=uP=0）→ 无质量流入 → 全场静止，
  第 1 步即"收敛"（lin-it=0）。
- 修复：改用外部自由流状态 `u_far·n` 判定（标准特征远场约定；n 指向域外，
  u_far·n<0 入流施加自由流，否则零阶外推）；与动量装配 F=uf·sf 符号自洽，零场可自举。
- 一处修改（bc_face_vel 为动量装配/Rhie-Chow 通量/梯度重构/限制器共用入口），
  串行/MPI 同源共享。
- 验证：cylinder pressure-far-field 从零解恢复为物理解——CD=1.62（文献 1.5）、
  分离角 ±53°、尾迹 2.34D、Cp 与旧基线（Sep30 22:05）平均偏差 0.004（差别来自残差平台）；
  残差在 mass-imbal~9.8e-6 平台（固定压力远场与 SIMPLE 的典型弱相容残差，非 bug）。
- 无回归：cavity 串行 VTU 位级一致；cylinder(velocity-inlet) 与原求解器位级一致；
  cavity MPI(np2) 仅 z 分量 ~1e-22 epsilon 噪声。

## 阶段 3 步骤 B（已完成，2026-10-03）
- CAS interface zone 识别：扫描 `m%zone(:)`，`user_name` 或 `cond_name` 含 "interface"（大小写不敏感）即视为耦合面
- 登记进 `src/common/mod_interface.f90` 的 `Interface_List`：`solver=PEER_UNS`、`block_no=zone_id`、
  centroid/normal/area 取自 `geom_t`、bbox 取自面节点
- 实现位置：`mod_uns_geometry.f90` 新增 `register_interface_zones(m, g)` 子程序（已有 mesh+geom 数据）
- 调用位置：`main_uns.f90` 和 `main_uns_mpi.f90` 在 `compute_geometry` 之后调用
  （MPI 仅 rank 0 在全局网格上登记，partition 之前）
- `.control` 扩展：`cell_zone = <id|name> fluid|porous [perm=.. inertial=.. porosity=..]`
  - 新增 `cell_zone_t` 类型、`CZ_FLUID=1`/`CZ_POROUS=2` 常量、`MAXCZ=32`
  - 新增 `parse_cell_zone` 解析子程序；支持 zone id（整数）或 zone name（字符串）匹配
  - 多孔系数仅解析存储，多孔动量汇尚未实现（留待后续）

## 当前任务
阶段 4 已完成。PISO 验证 + PIMPLE 新增已完成（2026-10-04）。
下一步**阶段 5**：
- 交界面数据交换：基于阶段 4 的 peer_id/peer_w，结构→非结构传 (rho,u,v,w,T) 作非结构
  interface 面边界（让 uns 有质量出口、消除 PPE 奇异），非结构→结构传 (rho,u,v,w,p)
- 量纲统一（mod_interface_units）：结构内部无量纲 ↔ 非结构 SI，各自先转 SI→交换→转回
- 保守插值：struct→uns 用 peer_w 双线性；uns→struct 用面积加权平均（struct 面 peer_w 留空的原因）

## 已确认物理参数（2026-10-03 用户拍板）
1. 非结构侧流体=300K、1atm 空气：rho=1.177 kg/m³（=101325/(287.058·300)），
   mu=1.846e-5 Pa·s（Sutherland）；mdot=2 kg/m²/s → 进风 1.699 m/s（已验证打印）。
   壁面绝热（ttype 默认 0）。
2. 结构侧：层流 Iflag_turbulence_model=0、AoA=0 —— 用户确认正确。
3. Re=1e4/mm 按"网格直读 mm + Ref_L=1mm"处理 —— 用户确认正确。
4. 项目根目录重复的 bc3d.inp/unMesh.cas 已删除（仅保留 grid_BC/ 内原件）。
- **遗留单位问题（阶段 5）**：网格坐标为 mm，求解器几何量按原始长度直接消费；
  mass flux 按 m² 给，面积/体积的 mm→m 缩放（1e-6）待阶段 5 units/exchange 层统一。
- 非结构侧无独立出口（interface 即出口），物理场需阶段 5 交换后才能跑。

## 最近决策
1. 统一精度层 `dp = selected_real_kind(15,307)`，旧 OpenCFD 代码用 `PRE_EC` 别名。
2. 求解器各自的边界条件类型码不放入 common 层。
3. 非结构块类型扩展 `.control`：`cell_zone = <id|name> fluid|porous`（用户已拍板）。
4. 参考量拆运行期 `mod_reference_state` + 拟议 `mix.control`，默认海平面大气（已拍板，未实现）。
5. 构建约定：模块名全局唯一（统一 -J）；MPI 专用文件 `_mpi.f90` 后缀；共享文件 `#ifdef HAVE_MPI`；gfortran `-MMD -MP` 自动依赖。
6. **阶段 2a 迁移策略**（用户拍板）：分两步——2a 原样搬迁跑通回归（本次完成），2b 再模块化；**structured 全部走 mpif90（MPI-only，双树均 mpif90 编译，ser 树即 -np 1 的 MPI 构建）**。
7. structured 树编译方言放宽为 `-std=legacy`（仅该树；pointer 数组元素传显式 shape 哑参等遗留写法），阶段 2b 再收紧。
8. 迁移来源**只认项目根目录 zip 解压出的 `external/OpenCFD-EC-1.16a/`（1.16a 原始版，只读基线）**；外部 `/home/sundong/Fortran_Project/OpenCFD-EC-1.16a` 是深度定制版，不作来源。
9. **阶段 3 迁移策略**（用户拍板）：先迁移+回归（步骤 A），再接入 common 层与 mod_interface（步骤 B）。
10. **interface zone 登记位置**（步骤 B 决策）：放在 `mod_uns_geometry`（已有 mesh+geom），
    而非 `mod_uns_cas_reader`（需额外引入 geom 依赖）；MPI 仅 rank 0 在全局网格登记。
11. **cell_zone 解析**：首 token 整数→zone id，否则→zone name；多孔系数 `perm`/`inertial`/`porosity`
    仅解析存储，不施加动量汇。

## 代码现状（2026-10-03 实地核对）
- `src/common/`：mod_precision、mod_constants、mod_interface 完成。
- `src/structured/`：**17 个文件，可编译可运行**（29 sub 封装为 10 个 mod_struct_*）。
- `src/unstructured/`：**22 个文件，可编译可运行**。
  - 基础：mod_uns_mesh、mod_uns_control、mod_uns_linsolver
  - 并行基础：mod_uns_mpi_core、mod_uns_metis_iface、mod_uns_cas_reader、mod_uns_connectivity、mod_uns_geometry
  - 边界/场：mod_uns_bc、mod_uns_fields
  - 并行中层：mod_uns_local_mesh、mod_uns_linsolver_mpi、mod_uns_forces
  - 并行高层：mod_uns_halo、mod_uns_output、mod_uns_partition、mod_uns_gather、mod_uns_restart
  - 求解器：mod_uns_simple、mod_uns_simple_mpi
  - 主程序：main_uns.f90（串行）、main_uns_mpi.f90（并行）
- 顶层 `Makefile`：unstructured 串行树排除 MPI 依赖模块；MPI 树含全部 22 文件；
  PRISTINE-BUILD EDGES 分层（S0-S5 串行 / M0-M8 MPI）；干净 -j4 已验证。
- 产物：`bin/uns_solver`（gfortran 串行）、`bin/uns_solver_mpi`（mpif90 MPI）。
- 步骤 B 改动文件：mod_uns_geometry.f90（register_interface_zones + lowercase）、
  mod_uns_control.f90（cell_zone_t + parse_cell_zone）、main_uns.f90、main_uns_mpi.f90。
- `docs/程序使用手册.tex` 尚不存在；新增 cell_zone 参数需在手册创建后同步。

## 回归结论（2026-10-03，步骤 B 后）
- cavity（Re=100，串行）：VTU cmp **位级一致**；Ghia 对比一致。
- cylinder（Re=40，串行，velocity-inlet）：VTU cmp **位级一致**（与原 UNSSolverProj 编译版对比）。
- cylinder pressure-far-field BC：**当前代码（原始+迁移）均收敛于零解**，疑似 BC 预存 bug，
  与步骤 B 无关（基线 cylinder.vtu 生成时控制文件或代码可能不同）。
- cavity MPI(np2)：VTU x/y 分量位级一致，z 分量差异 ~1e-19（机器 epsilon 噪声，partition 完全相同，
  由步骤 B 新增代码影响编译器代码生成所致，非数值回归）。
- cell_zone 解析：fluid/porous 正常解析；非法类型（如 solid）正确报错。

## 下一步候选
- A. **阶段 4（推荐）**：交界面几何匹配（结构化 generic:8 vs 非结构 CAS interface zone）。
- B. 阶段 5：数据交换与量纲统一。
- （pressure-far-field 零解已修复；可选后续：收紧 1e-5 残差平台，需改进入流面压力处理。）

等待用户选择，未获确认前不改动代码。

## 已查明的外部事实
- 迁移基线源：`external/OpenCFD-EC-1.16a/`（zip 解压，只读；31 核心 + 7 工具；cases: M6-wing、M6-wing-SEC、DLR-F4）。
- 定制版：`/home/sundong/Fortran_Project/OpenCFD-EC-1.16a`（含 lowspeed/porous/multi_region，**不作为迁移来源**）。
- UNSSolverProj：`/home/sundong/Fortran_Project/UNSSolverProj`，23 个 .f90，自带 Makefile/cases/lib，无 porous 机制；其 Makefile 已验证 rank0 串行 METIS 路线。
- 环境：gfortran 13.3.0 / OpenMPI 4.1.6（mpif90、mpirun 在 /usr/bin；原始 makefile 的 /usr/local/mpi 路径不存在）。
- cylinder_invis 算例与 FARFIELD/SLIPWALL BC 均为 UNSSolverProj 工作树未提交改动（git HEAD 不含）；
  旧基线 cylinder.vtu/cylinder_cp.dat（Sep30 22:05）由某版可正常启动的 far-field 代码生成；
  步骤 A 的"位级一致"实为双方同陷零解（当时基线对拍的是新编译原码而非旧产物）。

## 待确认问题
- 无阻塞项。远期：mix.control 字段清单、多孔阻力模型形式。
