# 进度总览 (progress)

> 最后更新：2026-10-08（**Betchen 2006 三验证算例试跑 + BJ 参考解对比完成**）

## 已完成

- [x] **Betchen 2006 界面验证：三算例试跑 + BJ 参考解对拍（2026-10-08，无源码改动）**：
  - 输入/脚本在 `cases/betchen/`（`gen_{bj,plug,ht}.py`、`.cas`、`.control`、
    `compare_bj.py`/`plot_*.py`、`images/bj_compare_ref.png`）。
  - **算例1 BJ**（流上多孔下，Re_H=1，Da=1e-2/1e-3，默认连续界面=Betchen 式
    (13)/(14)）：对拍 `/mnt/c/temp/validate_case/BJ_1.csv`/`BJ_2.csv`（u 按
    **纯流体段均值 U0** 归一化，y 以 H 计）：峰值 1.408 vs 1.387（+1.5%）/
    1.459 vs 1.457（+0.1%）；整条 u(y)/U0 L2 残差 = 参考 RMS 的 4.09% / 2.88%，
    判定复现成功。
  - **算例2 PLUG**（流-多孔-流 3H/2H/3H，Re=1；+高Re 5H/5H/50H）：复现「三段
    线性压力 + 界面压力梯度不连续」；数字化 Fig7 缺失，定性判定。
  - **算例3 HT**（底部加热铝泡沫 LTNE，Calmidi-Mahajan：ε=0.9118、K=1.8e-7、
    h_sf=128.5 W/m²K、a_sf=139.3 /m）：情形A 完全收敛（dT→0 @1377），Tf/Ts
    300→310K；**情形B（+2mm 空气隙，多 cell zone）Ts 全场保持 300K**（疑固相
    能量方程底加热壁面热未接线，待查）；数字化 Fig9 缺失，定性判定。
  - **经验/教训**：CAS `(13` 区 id 令牌按**十六进制**解析（`10`→写 `a`）；
    参考 CSV 按 y 排序；参考 U0=流体段均值。

- [x] **文档节点：程序使用手册 ＋ 附录 A 参数总表 ＋ 更新记录 ＋ `control.ec.template`
  （2026-10-08）**：此前 4 个文件**均不存在**，本次补齐。
  `docs/程序使用手册.tex`（ctexart/xelatex，**23 页**）：§1 概述、§2 构建与运行、
  §3 输入文件、§4 结构化求解器、§5 低速求解器、§6 弱耦合（Dirichlet–Neumann 契约）、
  §7 算例集表、§8 验证与回归、§9 FAQ（10 条）、§10 更新记录、附录 A；
  `docs/91_appendix_params.tex`（A.1 `control.ec` / A.2 `*.control` / A.3 `mix.control`
  / A.4 `bc` 与边界类型码 / A.5 `cell_zone`（含 `h_sf`/`a_sf`）/ A.6 `tbc`/`tbc_plane`/`lid`，
  表 3–8）；`docs/93_changelog.tex`（逆序 2026-10-08 → 2026-10-03，阶段 1–11）；
  `docs/control.ec.template`。**验证**：`xelatex` 连跑 4 遍，末两遍零 error /
  零 undefined reference / 零 multiply defined；标签定义 26 / 引用 17 全命中。
  涉及文件：上述 4 个 docs 文件 + `.gitignore`（LaTeX 中间产物）。

- [x] **阶段 11：界面类型自动分派表（2026-10-07）**：`Interface_FACE_TYPE` 增
  `cz_type / loc_face / iface_type`；`classify_interface`（struct 侧恒 comp+fluid；
  跨求解器仅由 uns 侧 cztype 决定 fluid/porous；uns↔uns = `IFACE_UNS_FLUID_POROUS`）、
  `iface_exchange_recipe`（pure 子程序）、`iface_type_name/quantity_name/role_name/
  iface_recipe_string`、`report_interface_dispatch`。uns：`tag_interface_cell_zones`
  （读真实 `m%cztype`，须在 `resolve_cell_zones` 后）+ `register_interface_zones`
  记 `loc_face` + `mod_uns_driver` 生产接线（`PEER_STRUCT/IFACE_CZ_FLUID`）。
  coupling：`dispatch_interfaces`（匹配后分类/回填/打印）。struct：`mod_struct_grid`
  标 fluid。**仅分类+报告，不改 `main.f90` 数值行为**。
  验证：`bin/match_test`（grid_BC）250/250 全 `comp-fluid<->lowspeed-fluid`、
  `unclassified=0`、断言 PASS；`couple_channel` np2 → fluid、`couple_porous` np2 →
  porous（RC=0）；`regress/m6wing/run_regression.sh` 全 IDENTICAL（`flow3d.dat`
  md5 `dc134a2d196422043ecad7c86ac8f898` 不变）；全量 make RC=0。
  涉及文件：`src/common/mod_interface.f90`、`src/coupling/mod_interface_match.f90`、
  `src/coupling/test_match.f90`、`src/structured/mod_struct_grid.f90`、
  `src/unstructured/mod_uns_driver.f90`、`src/unstructured/mod_uns_geometry.f90`。

- [x] **仓库瘦身：算例运行产物脱离跟踪 + 历史重写 + 远端同步（2026-10-06）**：
  **212 个**运行产物（**761.8 MB**）从**全部历史**移除并同步到 GitHub：
  `cases/**`、`grid_BC/**` 下的 `*.vtu / *.log / *.dat / *.out / *.tmp /
  *.part.map / __pycache__`、根目录 `LNTE1D.vtu`、`cavity.vtu`、`build_*.log`、
  `m6_*.log`、`output_para.out`，以及 `checkpoint_2026-10-06_C2.tar.gz`。
  跟踪文件数 **536 → 324**。
  判定依据（逐个分类，非拍脑袋）：`grep` 源码确认 `Step_mess.dat`/`part_grid.dat`/
  `partation-auto.dat`/`mesh-quality.dat`/`*.part.map`/`flow3d*.dat`/`output_para.out`
  均为**求解器运行期写出**；各 `cases/*/README.md` 均把 `.vtu/.log` 记作"产出"并给出
  复现命令；删除前全文检索确认**没有任何** README/脚本把 `.vtu` 当输入或 oracle。
  手段：`git rm -r --cached --pathspec-from-file`（精确 212 条清单）→ `.gitignore`
  新增"运行产物不入库"规则 → `git filter-branch --index-filter`（本机**无**
  `git-filter-repo`）→ `git reflog expire --expire=now --all` +
  `git gc --prune=now` → `git push --force-with-lease origin main`。
  安全网：`/home/sundong/mixsolver_pre_slim_backup/`（`repo_pre_slim.bundle` 98 MB
  `git bundle verify` 通过 = 完整旧历史；`cases`/`grid_BC` 硬链接快照；checkpoint 副本）。
  ⚠️ 原拟用的本地分支 `backup/pre-slim-2026-10-06` 被 `filter-branch --all` 一并重写，
  已删除 ⇒ **恢复只能走 bundle**（详见下条校验记录）。
  实测结果：`.git` **120 MB → 46 MB**，跟踪 **324** 文件，历史中 `.vtu` 计数 0，
  `git push --force-with-lease` 零告警，`main` 与 `origin/main` 0/0。
  工作树里 761.8 MB 产物**保留在磁盘**（转为未跟踪 + 被 .gitignore 挡住），
  随时按各 README 复现小节重算。规则更新见 `.trae/rules/project_rules.md`
  「版本控制约定」（新增运行产物条款 + 一次性重写记录 + 体量提示刷新）。

- [x] **关闭自检缺口③④（2026-10-06）**：
  ③ `src/unstructured/mod_uns_driver.f90` 的 uns restart 读取由整段停用（TEMP
  WORKAROUND）改为 `#ifdef HAVE_MPI` 守护——MPI 树恢复 `read_field_dump(serial=.true.)`
  （couple_channel `couple_restart=1` 逐位续跑：dump → `restarted from dump` → iter3，
  RC=0），串行树走 `#else` 打印 WARNING 且不崩（RC=0）；`ier2/src_file` 声明并入
  `#ifdef`，串行树零 unused 警告。
  ④ 重建并**持久化** M6-wing 位级基线到 `regress/m6wing/`（原先仅存 /tmp，易失丢失）：
  按既定配方 iconv+legacy 重建 external 参考，np1 / t_end=0.501 / Kstep_save=50 →
  `flow3d.dat` 13,600,032 B，md5 **dc134a2d196422043ecad7c86ac8f898**（与历史记录逐位
  一致）；8 个产物对新 struct_solver 全 IDENTICAL；附 `run_regression.sh`
  （`M6_REBUILD_REF=1` 可从源重建参考）+ `README.md` + `baseline/md5sums.txt`。
  验证：`run_regression.sh` PASS；全量回归 `make -j4 all mpi structured(_mpi)
  unstructured(_mpi) units_test match_test coupling_test` RC=0，units_test 3/0、
  coupling_test np2 6/0。涉及文件：`src/unstructured/mod_uns_driver.f90` +
  新增 `regress/m6wing/*`。

- [x] **修复串行 `units_test`（2026-10-06，缺口②）**：`bin/units_test` 链接
  `$(FC)`→`$(MPIFC)`（同 `bin/mixsolver`），并补 `test_units_exchange.o -> 耦合层`
  的 pristine 边（耦合层 MPI-intrinsic）。验证：`make -j4 units_test` RC=0，
  模拟竞态 `-j4` RC=0，`./bin/units_test`/`mpirun -np1` 均 3/0 PASS。
  涉及文件：`Makefile`。仍未处理：mod_uns_driver restart 停用、/tmp/m6reg 基线丢失。

- [x] **修复串行 `make all`（2026-10-06，缺口①）**：耦合层（`src/coupling/*`）与
  混合驱动 `src/main.f90` 改为在 ser 树也用 `$(MPIFC)` 编译（先例=structured）；
  `main.f90` 的 `mod_uns_restart`/`write_field_dump` 用 `#ifdef HAVE_MPI` 守护。
  验证：`make -j4 all`/干净 `make -j1 all` RC=0，`bin/mixsolver` 可运行（np2 跑通
  couple_channel）；回归 mpi/structured(_mpi)/unstructured(_mpi)/coupling_test 全绿。
  涉及文件：`Makefile`、`src/main.f90`。仍未处理：串行 `units_test`、
  mod_uns_driver restart 停用、/tmp/m6reg 基线丢失。

- [x] **构建/回归自检（2026-10-06，无源码改动）**：`make -j4 mpi` +
  `structured(_mpi)` + `unstructured(_mpi)` 全 RC=0；`coupling_test` np2
  PASS（6/0）；couple_channel 拷贝 n_couple=3 冒烟 RC=0。**既有缺口（未修）**：
  串行 `make all` 与串行 `units_test` 因 main.f90 / mod_coupling_exchange.f90
  在 gfortran 串行树中无条件 `use mpi` 而失败；mod_uns_driver restart 读取
  被 TEMP WORKAROUND 停用；/tmp/m6reg 基线丢失。MPI 为唯一可用耦合构建。

- [x] **阶段11 C1：低速-多孔界面（2026-10-06 勾销，零新代码）**：
  用户明确三类界面域归属——低速-多孔界面**只存在于非结构域内部**
  （uns 流体 + uns 多孔 cell zone 间内部面），不跨 solver、不经跨组
  交换层。直接复用阶段 9 的内部 BJ Robin 通量（bj_alpha：
  C=μA(α/√K)/(1+α·d_Pf/√K)，等值反向），验证资产 cases/beavers_joseph/
  （双层 4×4 ODE：α=0 RMS 0.14%，α=1/2 RMS 1.09%/1.29%，np2 2.2e-7）。
  docs/plan.md 阶段 11 节标记 C1 完成。下一任务：可压缩（struct）–
  多孔（uns）跨组界面（已由下一条完成）；界面分派表仍待办。

- [x] **阶段11 C2：可压缩（struct）–多孔（uns）跨组界面（2026-10-06 完成并验证）**：
  算例 `cases/couple_porous/`（struct 流体段 x=0–100 mm ↔ uns 多孔床
  x=100–200 mm 弱耦合，界面 x=100 mm）；对拍资产 = 单求解器全流域参考
  `uns_full/`（一个 uns：zone 2 fluid + zone 8 VC:porous）。**求解器本体零改动**，
  两处缺陷都在网格生成/后处理：① `uns_full/gen_full.py` 的 `cid()` 把 x 串联写成
  z 并联（k 最慢 ⇒ zone 2/8 各成一条 z 全长薄片）：床梯度仅解析值 ~60%
  （1.78 vs 3.02 Pa/mm）、max|u| 42.87 > 入口、zone 互换后逐位相同（z 镜像假象）
  → 改为 `1+k+j*NZ+i*NY*NZ`（x 最慢）；② `compare_iface.py` 取窗 ±1.25 mm 跨了
  两个 dx=2 mm 平面（把 295 Pa 与 70 Pa 平均成 182.6 Pa ⇒ 幻影 112 Pa 跳变）
  → 改绝对窗 `[100,102) mm`。
  验收：床梯度 −3021.3 vs 解析 −3021.7（0.01%，`check_bed_gradient.py` PASS 兼
  回归守卫）、界面压力跳变 −0.08 Pa、struct 面速度 −1.37%、耦合侧床梯度
  −2959.2 Pa/m 与界面速度 34.233 代入 μu/(εK)+ρβu² 逐位吻合。
  新记录（开放项）：uns 单求解器**绝对压力水平**含 ≈−250 Pa 内部偏置——与多孔
  无关（全流体对照内部严格平 −251.17）、与出口 BC 类型无关（`outflow` 仅差 3 Pa）、
  严格 ∝u² ⇒ 跨求解器比对不得用绝对压力。详见 `cases/couple_porous/README.md`
  §3–§6.1，工具 `uns_full/plane_profile.py`。**（2026-10-06 已定位：整场电平/
  压力零模，入口印标签；见本文件 §「待办 → 阶段 11」两条 + README §5.1。）**
  **（2026-10-07 已修复：内部对流系数重复乘 ρ，`F = ctrl%rho*fld%flux(i)` →
  `F = fld%flux(i)`；平台 −251 →≈0、u₁/u_in 0.9168→1.000000、床组平移 +0.009 Pa，
  梯度不变；出口 ±10 另靠"出口压力与入口同速率 ramp"。见 README §5.3。）**

- [x] **阶段11 流 B：couple_channel 交换层界面连续性修法（2026-10-06 完成）**：
  界面连续性严格在耦合交换层内、用两侧状态解决（不借 struct 入口）。
  4 版迭代后定稿 **Dirichlet-Neumann 特征界面**：struct 亚声速出口只
  接收 uns 背压，ghost 用 boundary_Farfield 同款线性化 Riemann 反射
  （db=d1+(pb−p1)/c1²，ub=u1+(p1−pb)/(ρc)·n，ghost=2face−inner，
  法向按 face 1..6 从 Block 面法向表构造）；uns 侧只施加速度 Dirichlet，
  删 set_interface_p 与压力重锚定（BC_INTERFACE 压力本零梯度）；
  背压松弛 alpha=0.3·min(1,iter/iface_ramp)（alpha=1 在 Ma=0.1 下
  ramp 后周期发散，alpha<1 不改收敛值）。被否方案：v1 ghost=2f−inner
  子步前设（周期-2 翻转 NaN）、v2 子步后设（回声不动点流量衰减）、
  v3 零阶保持直写首排格心（冻结 1.7 kPa 跳变）。结果（iter1000，
  run_flowb4c.log，无 OverLimit）：质量精确连续 33.71==33.71 m/s、
  压力跳变 <1 Pa（342.5 vs 342）；对拍 uns_full 单一求解器解
  （34.72 m/s、46.9 Pa）均值 −2.9%、界面压力 +296 Pa、近壁型线
  最大偏差 11.8%，剩余偏差为两解器壁面层物理差异（非交换层缺陷）。
  改动：src/main.f90、src/structured/mod_struct_driver.f90、
  src/structured/mod_struct_bc.f90（会话中曾被外部回退为含 Turbo
  脏版，已重建为 915 行干净版并 M6 位级回归
  md5 dc134a2d196422043ecad7c86ac8f898）、
  cases/couple_channel/compare_iface.py（flow3d.dat 含 ghost 解析修正）。

- [x] **阶段11 清理：弃用 IF_InnerFlow / IF_TurboMachinary 全套机制（2026-10-06）**：
  总压入口试验（couple_channel，P_In_Ratio=1.0123）np2 iter50 发散（78 m/s）
  且界面偏差不变，证明界面连续性与入口无关；用户拍板两套机制全部弃用。
  删除 struct 侧 5 文件：global 变量（IF_TurboMachinary/IF_InnerFlow/
  P_In_Ratio/Turbo_Periodic_seta/Turbo_w）、constants（BC_Wall_Turbo）、
  bc（Turbo 分派+三个 Turbo BC 子程序）、solver（旋转惯性力源项）、
  init（Turbo_P0/T0/L0、Ref_medium_usrdef、namelist 键、默认值、叶轮机
  init 分支、介质推导块、bcast 赋值；槽位号保留）、mpi
  （Umessage_Turbo_Periodic + 两个周期子程序的旋转分支）。教训：
  Fortran 大小写不敏感，Aos 即侧滑角 AoS，误删后已恢复。
  couple_channel/control.ec 恢复强制均匀入口；两算例 control.ec 删 5 个
  废 namelist 键；删试验产物 control.ec.forced/mix.control.ptot/run_ptot1.log。
  验证：ser/mpi 双树 0 error；M6-wing 50 步对拍原始基线（/tmp/m6reg，
  统一 -O2 -std=legacy），np1/np2 flow3d.dat 与 Step_mess.dat 均
  byte-identical；couple_channel n_couple=2 冒烟 rc=0。流 B（交换层
  界面连续性修法）转待办。

- [x] **阶段10 流场自动保存与耦合联合重启（2026-10-06，全部完成）**：
  mix.control 新增 `save_interval`（>0 生效）与 `couple_restart`。
  两侧在同一 coupling iter 末尾对齐保存：struct 侧 `struct_solver_save`
  → 原生 `output_flow`（flow3d.dat+Step_mess.dat）+ 新增
  `output_flow_nodes`（flow3d_node.dat，格心→Mesh3d.x 节点插值 Plot3D
  函数文件，格式随 Mesh_File_Format）；uns 侧 uns root 写
  `unMesh_restart.dat`（es24.16 ASCII dump）。struct rank0 写
  `couple_state.dat`（iter 号）。重启：struct 侧 force_restart 强制
  `Iflag_init=1` 读 flow3d.dat；uns 侧 `uns_solver_init` 新增
  restart_file 可选参，经 `read_field_dump(serial=.true.)` 串行读入
  （新增 serial 选项跳过 MPI_Bcast——耦合驱动不调 mpi_bootstrap，
  Bcast 会把 struct rank 卷进挂死）；两侧从 couple_state 恢复 iter0，
  循环 iter0+1 起，iface_ramp/omega 计数连续。grid_BC np2 验证
  （n_couple=6、save_interval=3）：连续 6 轮 vs 3+重启 3 轮，
  flow3d.dat 与 unMesh_restart.dat 均 **byte-identical**；ser/mpi 双树
  编译无新增 error。归档 grid_BC/cont6/、grid_BC/restart33/、
  check_restart.py。已知限制：dump v1 不含 T_s（LTNE 算例重启 T_s
  回 init 值）；struct 的 Kstep/tt 不恢复（仅影响输出命名）。

- [x] **阶段9 α/K 标定实验设计（2026-10-06，阶段 9 全部收尾）**：
  cases/calibration/ 数值虚拟标定框架（calib_problems 问题定义 +
  calib_invert LM 反演/FIM + run_calib 闭环 + calib_design_scan 噪声选型）。
  六参数按物理特征解耦为四子问题：P1 多流速压降→K+C_F（CRB 0.6%/5.5%，
  相关 0.65，条件数 12.5）；P2 混合层横向剖面→α_t（**设计发现：TC 噪声按
  绝对温度缩放，ΔT 必须 ≥50 K**，0.2% 噪声 CRB 2.2%，10 K 步长时 12.7%
  不可用）；P3 LTNE 发汗冷却两相测温→α_l+h_sf 联合反演（CRB 1.5%/0.9%，
  相关 0.76）；P4 界面 u(y) 剖面→α_BJ（2% PIV CRB 5%）。四子问题从
  ~2× 离真值初值 4–6 步 LM 收敛，恢复误差全部在 CRB 界内
  （K 0.19%/C_F 1.1%/α_t 3.8%/α_l 1.1%/h_sf 0.25%/α_BJ 1.8%）。
  反演器与数据源解耦（forward 可换求解器 VTU 采样），正演构型复用
  porous_plug/porous_disp/ltne_disp/beavers_joseph 四算例。纯算例/脚本任务，
  无源码改动，无需回归。

- [x] **阶段9 LTNE+弥散组合验证（2026-10-05 收尾，1D+2D 全部完成）**：
  - **修复纵向弥散 bug**：Bear 张量纵向分量原写 ρcp·α_l·u_d²（缺 /|u|，
    量纲 m²/s²），改为 ρcp/umag·(α_l u_d²+α_t(umag²−u_d²))
    （mod_uns_simple.f90 温度装配）。纯横向四工况（porous_disp×2、
    porous_graetz×2）修复前后 VTU **位级一致**。
  - **修复 MPI gather 缺 T_s**：mod_uns_gather.f90 补 sbuf_Ts/rbuf_Ts
    pack/Gatherv/unpack；main_uns_mpi rank0 fld_g 补 setup_porous_fields
    （否则 VTU 输出 guard 跳过 temperature_solid 段）。
  - 1D cases/ltne_disp/：100×2 hex，恒热流发汗冷却，disp_l=0/5mm；
    4 阶耦合 ODE 半解析（三次特征根，未知量 [c0,d1..3]、exp(r(x−L))
    防溢出，Nield 微观热流分配）；Tf 中位 0.4%/Ts max 0.83%，
    disp1 出口温升 973.46 vs 974.41（热量经弥散回流至 Dirichlet 入口）；
    8095/11219 步收敛；np2 两相 rel 1e-7。
  - 2D cases/ltne_graetz/：复用 graetz 网格（ε=0.4, k_s=1.0, Bi=0.845,
    Λs=0.101），disp_t=0/2.5e-4；参考=双相耦合二次特征值问题
    （4m 一阶块 y'=My，保留两相轴向导热、流体横/纵 Pe 分开、
    入口 θf=1+χs=0，Bi→0 退经典 Graetz λ1Pe=2.457≈π²/4）；
    180×40 剖面 max 1.4%/3.3%（仅 x=20mm 入口段），360×40 x 加密
    减半 0.7%/1.8%，x≥50mm ≤0.5%，bulk 全程重合；7018/6939 步；
    np2 两相 abs 2e-6K（rel 2e-7）。两算例均归档 README+图+日志。
  - α/K 标定实验设计已于 2026-10-06 完成（见顶部条目），阶段 9 全部收尾。

- [x] **阶段9 Beavers-Joseph 界面条件（2026-10-05 收尾）**：
  `cell_zone ... bj_alpha=α`（默认 0=禁用）+ 全局 `body_force=fx fy fz`。
  Robin 型界面通量 C=μA(α/√K)/(1+α·d_Pf/√K)，切向用 C、交叉项滞后、
  等值反向保动量；MPI 复用 serial 装配零改动。验证 cases/beavers_joseph/
  （开口槽道+体力驱动，K=1e-6，λ/dy=4）：对拍 4×4 双层应力跳变 ODE 精确解
  （非经典 BJ 滑移公式——多孔侧 Brinkman 层可分辨），α=0 RMS 0.14%、
  α=1/2 RMS 1.09%/1.29%；np2 vs 串行最大相对偏差 2.2e-7；
  cavity/Ra10 回归位级不变。归档 README/plot_bj.py/images。

- [x] **阶段9 二维多孔介质强制对流 Graetz 验证（2026-10-05 第三批）**：
  cases/porous_graetz/，半平行板通道（a=5mm，L=0.45m，180×40），中心线
  symmetry + 恒壁温 Tw，Darcy 活塞流（K=1e-10, eps=0.4, k_s=k_f 使
  k_eff=k_f）。两工况（mol U=0.1 / disp U=1+disp_t=2.5e-4，Pe≈20）
  对拍 Graetz 级数：充分发展 Nu_a=π²/4=2.4674，误差 0.42%/1.38%；
  剖面偏差 ≤0.026（入口）/≤0.013；两工况体均温度在 x* 坐标坍并。
  5601/6238 步收敛，np2=串行。教训：Bear 弥散随速度在近壁 Brinkman 层
  自动衰减（无弥散层），√K/dy 须 ≤0.08（首排 u≥0.988U）才能对均匀到壁
  解析解；首轮 K=1e-9（√K/dy=0.25）渐近 Nu 偏 16%。

- [x] **阶段9 热弥散 + 各向异性渗透率（2026-10-05 第二批）**：
  - 各向异性渗透率：cell_zone 新增 perm_xx/perm_yy/perm_zz（轴对齐对角张量，
    未给分量回退标量 perm），fld%perm_dir(3,ncells)；Darcy 汇按动量分量取对角元。
    验证：张量一致性 0.085%、Kxx 减半→压降精确加倍（0.122%）、标量路径 0.000%。
  - 热弥散：Bear 分量式张量 D_dd=ρcp(α_L u_d²+α_T(|u|²−u_d²)/|u|)，
    fld%disp_l/disp_t；能量面扩散加 n·D·n 投影（LTE 的 k_eff、LTNE 流体相 kf）。
    验证：2D 温度混合层 α_t=1e-3 vs erf 相似解偏差 ≤1%（cases/porous_disp/）。
  - **修复 velocity-inlet 温度 token bug**：第 4 个 token（静温）此前不被解析，
    bc_face_T 锚定 tval=0 → 独立运行全场温度排空到 0；现解析可选 token。
  - Ra10 回归不变（748 步/Nu=1.078）；混合层 np2=串行（188 步）。
  - 注意：erf 相似变量用热扩散率 κ=Γy/(ρcp)，勿用 Γy 本身。

- [x] **阶段9 VC 自动识别 + porous-plug 验证（2026-10-05）**：CAS cell zone
  名 `VC:porous/fluid` 自动判定块类型（显式 cell_zone 优先，CZ_AUTO 系数行），
  zone 名字段 32→128；一维多孔塞 Darcy 压降误差 0.000%、Darcy-Forch 0.002%，
  归档 cases/porous_plug/；Ra10 回归 Nu=1.078、np2 与串行一致。
  顺手修复 setup_porous_fields 中 cz%h_s 旧笔误（应为 h_sf）。

- [x] **inlet_ramp + 耦合界面压力/质量修正（2026-10-05）**：`ctrl_t%inlet_ramp`
  （默认 1=禁用）+ mod_uns_bc `g_inlet_ramp_factor`/`set_inlet_ramp_factor` +
  simple_run(_mpi) 每外迭代设置；main.f90 uns 驱动加①界面压力重锚定
  （去 struct 绝对压偏置、保留梯度）②界面速度均匀法向修正保总边界通量为零。
  修正后耦合 it=1 mass-imbal 3.9e-10，但 grid_BC（Ma=3 可压 vs 不可压腔体）
  物理不兼容弱耦合，~it=50 仍 NaN；耦合算例须 inlet_ramp=1（与质量修正互斥）。
  下一步待拍板：换物性兼容算例 / 界面改通量型 BC。

- [x] **槽道验证 + outflow 出口边界（2026-10-05）**：
  - plane Poiseuille 槽道（300×100×1，L=30H）：mdot-inlet 从零场启动稳定
    （Re100 it=1 imbal 5e-13；Re1000 10× 冲量 1033 步）；上会话 LTNE 启动发散
    确认为该网格特异（已修），非 mdot-inlet 本身。
  - 诊断 pressure-outlet 出口段伪调整（末列 Ucl 1.5→1.71、vmax 0.115U），
    对照排除入口/格式/并行，定性为固定压力边界的离散不动点。
  - 新增 `bc = <zone> outflow`（BC_OUTFLOW=8）：法向零梯度 + 全局质量缩放
    （outflow_mass_sums/scale/rescale，MPI allreduce；POUTLET/FARFIELD 不计入
    m_req）；纯 Neumann PPE 销钉 cell 1 因缩放 Σrhs=0 而相容。
  - Re100 出口末列：Ucl 1.7091→1.6420、vmax 0.1152→0.0402、L1 11.45%→5.63%；
    内部场与 pressure-outlet 一致；np2 一致。已知限制：末列 1/β 外观 +3.5%、
    ~9% 凸起（一阶零梯度固有）。
  - 回归：cavity / cylinder(far-field) / porous_Ra10 / LTNE 全 PASS。
  - 归档 cases/channel/（README、双出口控制文件、log/vtu、plot_channel.py、
    split_channel_zones.py、诊断对照 vinlet/ho/rlx）。手册 N/A。

### 阶段 1：项目骨架与公共模块（基本完成）
- [x] 创建目录骨架：`src/common`、`src/structured`、`src/unstructured`、`src/coupling`
- [x] `src/common/mod_precision.f90`：dp/sp/ip/i8、PRE_EC 兼容别名、pi、eps/tiny
- [x] `src/common/mod_constants.f90`：气体常数、比热比、参考量、SI 转换系数、标准大气
- [x] 需求与计划文档就位（`docs/程序功能说明.md`、`docs/plan.md`）
- [x] memory-bank 五个文件初始化（2026-10-03）
- [x] 顶层 `Makefile`（2026-10-03）：all/mpi/structured(_mpi)/unstructured(_mpi)/common/clean/help；
      build/ser 与 build/mpi 分离；gfortran -MMD -MP 自动依赖；DEBUG=1；已验证串行/MPI、
      增量、-j4 干净构建（3 次）、缺失 main 友好报错；根目录游离 .mod 已清理
- [x] `mod_unit_convert.f90` + `mod_reference_state`（2026-10-06 勾销：阶段 5
      已落地，实际文件名 `mod_reference_state.f90` + `mod_interface_units.f90`；
      本条为「阶段 5 前做」的计划项，与阶段 1 收尾区那条同源）

### 环境准备
- [x] OpenCFD-EC-1.16a.zip 已放入项目根目录（另已有解压副本 `/home/sundong/Fortran_Project/OpenCFD-EC-1.16a`）
- [x] lib/ 库文件核对：metis/libmetis.a、parmetis/libparmetis(_gnu_mpi).a、tecplot/libtecio.a
- [x] 编译器核对：gfortran 13.3.0、OpenMPI 4.1.6（mpif90/mpirun 就位）
- [x] **版本控制初始化（2026-10-06 完成）**：`git init -b main` + `.gitignore`
      （`build/`、`bin/` 为原待办要求；另补 `*.o`、`*.mod`、`__pycache__/`、`*.pyc`、
      编辑器临时文件。**刻意不忽略** `lib/{metis,parmetis,tecplot}/*.a` —— vendored
      预编译第三方库，链接必需且不可由本仓库重建；也不忽略 `regress/m6wing/baseline/*`）。
      首个提交 `97f247f`「chore: 初始化版本库（首个提交 = C2 完成节点快照）」：
      536 文件、`.git` ≈119 MB（ASCII VTU 压缩 ~7×）、分支 `main`、工作区干净。
      说明：`cases/` 下生成产物（`*.vtu`、日志等）按原待办**未**忽略，一并入库；
      日后若要瘦身可 `git rm --cached` ＋ `commit --amend`（本地无远端，安全）。
      身份仅设仓库本地 `sundong <sundong@localhost>`（全局 user.name/email 为空）。
      **工作流已固化（2026-10-06）**：每完成一个小节点 → 跑验证 → 更新 memory-bank →
      自动 `git add -A && git commit -F -`（**无需询问**），一个节点一个提交；
      验证未通过/半成品不提交（保留脏工作区或 `git stash`）。规则：
      `.trae/rules/project_rules.md`「版本控制约定」。
      **远端同步（2026-10-06）**：`origin = git@github.com:cfddog/MixSolver.git`（SSH，
      密钥 `~/.ssh/id_ed25519`）；首次 `git push -u origin main` 成功 —— 594 对象 /
      93.5 MiB 包体 / ≈2.9 MiB/s，**无大文件告警**；HEAD == `origin/main` == `8f40a31`。
      工作流扩展为「验证 → 记账 → commit → **push**」，此后自动执行。

## 待办

### Betchen 2006 界面验证（2026-10-08 试跑）
- [x] **BJ 参考解对比（Da=1e-2/1e-3）**：峰值 1.408/1.459 vs 参考 1.387/1.457，
  L2 = 参考 RMS 的 4.09%/2.88% → 复现成功（`cases/betchen/compare_bj.py`）。
- [~] **PLUG 参考解对比（Da=1e-2=plug-1、Da=1e-3=plug-2）对比不佳，记为后续待办**：
  对拍 `/mnt/c/temp/validate_case/plug_{1,2}_{u,p}.csv`（中心线 `u/U`、`p/(ρU²)`
  沿 `x/H`，入口为抛物线 `6U(y/H)(1-y/H)`）。**问题**：
  ① 充分发展段 u/U 吻合（1.497 vs 1.503）但**多孔/界面段速度偏差大**
     （plug-1 ~4%、plug-2 界面处出现振荡尖峰 1.52→1.06，ref 平滑 1.24→1.14）；
  ② **压力绝对电平整体偏高 ~1.2–1.3×**（plug-1 入口 423 vs 345、plug-2 3035 vs
     2264），两个 Da 倍率相近 ⇒ 疑**渗透系数/有效黏度约定差异**（我方
     `mu_eff=mu/ε` 且 `inertial=0`，参考为外禀式 `ε·mu/K` 且带 `cE`），
     也可能是中心线 vs 面均值取样差异。
  **待下一步**：① 核对 `mod_uns_simple.f90` 多孔源项/有效黏度约定；② 支持抛物线
  入流；③ 界面加密网格复跑；④ 复核压力归一化（U 为抛物线均值而非均匀入口 U0）。
  相关脚本：`cases/betchen/compare_plug.py` + `images/plug_compare_ref.png`。

### 阶段 11（流 B、C1、C2、分派表、绝对压力偏置 **均已完成**）
- [x] 弃用并删除 IF_InnerFlow / IF_TurboMachinary 全套机制（2026-10-06 完成，位级回归通过）
- [x] **流 B：可压缩（struct）–低速（uns）界面（2026-10-06 完成，
  Dirichlet-Neumann 特征界面，详见「已完成」条目）**：
  iter1000 质量精确连续、压力跳变 <1 Pa；对拍 uns_full 剩余 −2.9%
  均值偏差判定为两解器壁面层物理差异。参考 uns_full/unMesh.vtu
  （x=100 面 mean u 34.72、p gauge≈46.9 Pa），工具 compare_iface.py。
- [x] **C1：低速-多孔界面（2026-10-06 勾销）**：纯 uns 内部界面，
  复用阶段 9 BJ Robin（详见「已完成」顶部条目）。
- [x] **可压缩（struct）–多孔（uns）跨组界面（2026-10-06 完成并验证）**：
  详见「已完成」顶部 C2 条目 ＋ `cases/couple_porous/README.md` §2–§6.1。
- [x] **（2026-10-06 新登记；诊断节点）`uns` 单求解器绝对压力水平 ≈−251 Pa 内部偏置**
  （plan 阶段 11 第 3 条）：**机理已定位（诊断节点完成）**——它是**整场电平
  （pressure null-mode）**而非梯度误差：内部 93 个平面**逐位相同**
  `p = −251.16752 Pa`、`u = 34.72218`（与入口通量 `mdot/ρ` 逐位相等）；
  电平由**入口边界单元**印上（入口邻格压力 = 物理解，但速度只有
  **0.91678·u_in**，即 −8.32%，与耦合"uns 首排 −9.67%"同族）后经零模内部
  **原样传播**，出口 Dirichlet 只钉住末列（p = −158.8）。系数
  `平台/ρu_in² = −0.176943`，且 **Δx（1/2/4 mm）、`inlet_ramp`（100 vs 10）、
  入口/出口 BC 类型全部无关**（逐位级），严格 ∝u²；床组整场平移 −251.4 Pa
  （std 0.71 Pa）⇒ 梯度精确不受影响。代码路径：
  `mod_uns_simple.f90:1035-1036`（入口动量行显式加**全量** ρu_in²A）
  ＋ `mod_uns_bc.f90:404-405`（入口面零梯度 `pf = pP`）
  ＋ `flux_rhiechow` `mod_uns_simple.f90:1239-1244`（**边界面跳过 Rhie-Chow
  压力梯度修正**）⇒ 压力 Dirichlet 补丁锚不住内部电平。
  **附带发现鲁棒性隐患**：出口 `pressure-outlet` 只有精确 `0.0` 收敛
  （`±10`/`100` 均 NaN/超时；`velocity-inlet`/far-field 入口亦 NaN）。
  取证：`cases/couple_porous/README.md` §5.1/§5.2 ＋ 新脚本
  `uns_full/bias_scan.sh`（含生成的 `sum.py`/`shift.py`）。C2/流 B 判据**不受
  影响**（只比梯度与界面连续性）。
  **（2026-10-07 更正＋修复）**：真正根因不是上面的入口路径，而是
  **内部对流系数重复乘 ρ**（`mod_uns_simple.f90` 原 `F = ctrl%rho*fld%flux(i)`；
  `fld%flux` 已是质量通量 kg/s）——内部对流被放大 ρ 倍、与边界面
  `F = ρ·(u_f·S_f)` 不一致，使开放边界盖章的偏移由 O(u²) 放大为 O((ρ−1)ρu²)，
  恰是量到的系数 `0.176943 = ρ−1`。改为 `F = fld%flux(i)` 后平台 −251→≈0、
  `u₁/u_in`→1.000000、床组平移 +0.009 Pa。修复与验证见下一条 ＋ README §5.3。
- [x] **（2026-10-06 登记；2026-10-07 修复并验证）修复上述偏置**：让绝对压力电平回到
  物理解，并顺带修边界鲁棒性。**根因**＝内部对流系数**重复乘 ρ**
  （`mod_uns_simple.f90` 原 `F = ctrl%rho*fld%flux(i)`；`fld%flux` 已是质量通量
  kg/s）⇒ 内部对流被放大 ρ 倍、与边界面 `F = ρ·(u_f·S_f)` 不一致，开放边界盖章的
  偏移由 O(u²) 放大为 O((ρ−1)ρu²)——正是量到的系数 `0.176943 = ρ−1`。**修法①**：
  `F = fld%flux(i)`。**修法②**：出口 `pressure-outlet` 冷启动失配（`inlet_ramp` 使首
  迭代 ap 极小、出口压力却满值 ⇒ `pval/(ρu²)`≈70 ⇒ 第二步爆）⇒ 出口/远场压力
  Dirichlet **与入口同速率 ramp**（`mod_uns_bc` 新增 `g_pval_ramp_factor`/
  `set_pval_ramp_factor`，`bc_face_p` 用之；`mod_uns_simple(_mpi)` 同步调用）。
  **验收全过**：baseline 平台 ≈0 Pa、`u₁/u_in`=1.000000、床组平移 +0.009 Pa
  （梯度 −3021.7 不变）、出口 `±10` 全收敛、nx50/nx200、m6wing 位级回归
  PASS、`units_test` 3/0、`coupling_test` np2 6/0；**耦合首排**（np2、400 iter）
  相对亏损 **−9.67% → −1.51%**、参考 @x=100 的 p 由 48.5 → **299.9 Pa**、界面跳变
  −0.49 Pa（§4）。详见 `cases/couple_porous/README.md` §5.3。
- [x] 界面分派表：按 cell_zone 类型 + solver 类型自动确定界面类型与
  交换量清单（plan 阶段 11 第 2 条，**2026-10-07 完成并验证**，见「已完成」顶部）。
- [ ] （可选，物理一致性）壁面层对齐：struct 近壁分辨率/壁面处理与 uns
  对齐以消除近壁型线 11.8% 偏差（非交换层任务）。

### 阶段 1 收尾（P0）
- [x] 顶层 Makefile（目标：structured / unstructured / all / mpi；自动模块依赖；链接 lib 库）
- [x] 确定 build 输出目录规范（build/ser、build/mpi），清理根目录游离的 .mod 文件
- [x] `mod_unit_convert.f90` + `mod_reference_state`（2026-10-06 勾销：阶段 5 已落地，
  实际实现文件名为 `mod_reference_state.f90` + `mod_interface_units.f90`，
  与条目原拟名不同；原条目移至"阶段 5 前与量纲重构合并"后**已完成但未勾销**）

### 阶段 2：结构求解器集成（P0，约 3-4 天）
- [x] **阶段 2a（2026-10-03 完成）原样搬迁 + 跑通回归**：
  - [x] 以 `external/OpenCFD-EC-1.16a/`（zip 解压的 1.16a 原始版，只读）为准，31 个核心文件迁入 `src/structured/`（7 个网格工具不迁）
  - [x] 29 个 sub_*.f90 GBK→UTF-8；sub_modules.f90 拆 5 个 mod_struct_*；主文件拆 mod_struct_flowvar + main.f90；real*8→PRE_EC
  - [x] precision_EC 薄壳（use mpi 保持 public；PRE_EC 自含定义，不泄漏 dp）
  - [x] Makefile 适配（双树均 mpif90；structured 树 -std=legacy -ffree-line-length-none；pristine 模块边）
  - [x] ser/mpi 双树 0 error，干净 -j4 并发构建通过；产物 bin/struct_solver(_mpi)
  - [x] M6-wing 回归（/tmp/ocfd_regression）：基线为原始源同选项编译；np1/np2 的 flow3d/SA3d/wall_dist/分区文件位级一致，残差/气动力/stdout 一致；唯一差异为原版 a0/d0 未初始化打印（不参与计算）
- [~] **阶段 2b（模块封装完成 2026-10-03，剩 s3/s4）**：按 plan 包成 mod_struct_* 模块，分 4 批自底向上、每批 np1 位级回归
  - [x] B1：mod_struct_scheme/flux/fdm/time，旧 6 文件已删，回归位级一致
  - [x] B2：mod_struct_bc/mpi/grid（Type_def1 内嵌模块保留在 grid 文件顶部），旧 8 文件已删，回归位级一致
  - [x] B3：mod_struct_solver（Residual/time_advance/turbulence×4/limitflow/filtering，SCC 循环同模块）
  - [x] B4：mod_struct_io + mod_struct_init + 删 sub_interfaces + Makefile 改 6 层模块链；干净 -j4 通过
  - [x] 29 个 sub_*.f90 全部封装为 10 个 mod_struct_*（src/structured 现 17 文件），四批 np1 回归均位级一致
- [x] s3：Global_Var 审计（84→71 成员；13 个单消费成员搬迁到 FDM_data/mod_struct_io/mod_struct_init；np1 位级一致；不深改签名）
- [x] s4a：structured 树收紧 -std=f2008（15 处机械修复：混合 kind 字面量/FORMAT 缺逗号/isnan→ieee_is_nan/access=append→position=append；干净双树 -j4 + np1/np2 回归位级一致）
- [x] s4b：新增 BC_Interface 耦合交界面类型（src/common/mod_interface.f90；bc=8 面登记+几何计算+BC 分派占位；M6 回归位级一致）

## 阶段 3：非结构求解器（UNSSolverProj）迁移 — 步骤 A/B/C 完成
- [x] 独立编译运行，与原 OpenCFD-EC 结果回归对比（2a/2b 均位级一致）

### 阶段 3：非结构求解器集成（P0，约 2-3 天）
- [x] 取得 UNSSolverProj 源码（位置：`/home/sundong/Fortran_Project/UNSSolverProj`）
- [x] **阶段 3 步骤 A（2026-10-03 完成）迁移 + 重命名 + 回归**：
  - 21 个 .f90 迁入 src/unstructured/，重命名为 mod_uns_*（mod_precision 删除，统一用 common 层）
  - main.f90→main_uns.f90，main_mpi.f90→main_uns_mpi.f90
  - Makefile：unstructured/unstructured_mpi 目标、分层 PRISTINE 边（S0-S5 / M0-M8）、
    MPI 依赖模块（gather/halo/restart/local_mesh/partition/mpi_core）串行排除
  - 干净 `make unstructured unstructured_mpi -j4` 0 error，双树通过
  - 回归：cavity/cylinder 串行 + MPI(np2) 与基线 cmp 位级一致
- [x] **阶段 3 步骤 B（2026-10-03 完成）接入 common 层与 mod_interface**：
  - [x] CAS interface zone 识别（user_name/cond_name 含 "interface"，大小写不敏感）
  - [x] interface zone 面登记进 Interface_List（mod_uns_geometry: register_interface_zones）
  - [x] `.control` 扩展 `cell_zone = <id|name> fluid|porous [perm=.. inertial=.. porosity=..]`
  - [x] 编译 0 error；cavity 串行位级一致；cylinder(velocity-inlet) 位级一致；
    cavity MPI(np2) x/y 位级一致（z 分量 ~1e-19 机器 epsilon 噪声）
- [x] **PISO 非定常验证（2026-10-04）**：
  cavity_transient (BDF2, dt=0.005, 600步) 串行 + MPI(np2) 与原版位级一致；
  MPI vs 串行差异 ~0.5% 为分区浮点累加固有现象，非回归。
- [x] **PIMPLE 非定常算法新增（2026-10-04）**：
  ctrl_t 增 pimple(逻辑) + n_outer_iter(默认1)；momentum/temperature_assembly 在
  PIMPLE 模式(ts_order>0 & pimple) 的时间项之上施加 alpha_u 欠松弛；pimple_run +
  pimple_run_mpi 每时间步 n_outer_iter 轮外迭代（动量+梯度+Rhie-Chow+n_correct 全量
  压力修正），历史在步首推进一次。验证：dt=0.1(CFL~1.4) PISO 发散 NaN、PIMPLE 稳定；
  稳态 cavity/cylinder 位级一致无回归。
- [x] **能量方程 + 多孔介质（2026-10-04，phase 12）**：
  - 多孔动量源（Darcy-Forchheimer）：S = -(mu_eff/perm)u - rho*inertial*|u|u，
    Darcy 线性项隐式进对角、Forchheimer 非线性项显式进 rhs；mu_eff = mu/porosity
    （Brinkman），内面 lf 插值、边界面取属主格。
  - 孔隙率修正输运系数：k_eff = eps*k_f + (1-eps)*k_s；(rho*cp)_eff =
    eps*(rho*cp)_f + (1-eps)*(rho*cp)_s；温度方程扩散/瞬态用有效值。
  - LTNE 双温度方程：流体 eps*rho*cp_f*(dT_f/dt+u·∇T_f) = ∇·(eps*k_f∇T_f) +
    h_sf*a_sf*(T_s-T_f)；固体 (1-eps)*rho_s*cp_s*dT_s/dt = ∇·((1-eps)*k_s∇T_s) +
    h_sf*a_sf*(T_f-T_s)。新建 solid_temperature_assembly（无对流，共享壁面热 BC）。
  - 控制参数：thermal_model=lte|ltne；cell_zone 扩展 k_s/cp_s/rho_s/h_sf/a_sf；
    setup_porous_fields 每单元格映射 zone 参数。
  - VTU 输出 temperature_solid（LTNE 时）。
  - 验证：cavity 稳态位级一致（max diff 1.8e-20）；LTE/LTNE/多孔串行+MPI 冒烟通过，
    T_f≠T_s 证实非平衡。
- [x] **Rayleigh-Bénard 自然对流验证（2026-10-08）**：
  - SIMPLE 稳态收敛判据补 dT_max（温度每外迭代变化），串行+MPI 同步；
    cavity 回归与原版完全一致（同迭代数、打印差 0）。
  - 次临界 Ra=1000/3000：Nu=0.996≈1、速度近零（纯导热解析解）；
    起对流离散区间 Ra∈(3000,5000)。
  - Ra=5e3/1e4/5e4：Nu=1.675/2.160/3.314 vs Hollands(1976) 1.90/2.39/3.44
    （偏差 −12%/−9.6%/−3.7%，侧壁约束所致）；80² 加密网格 Nu=2.145，网格收敛；
    冷热壁热平衡误差 <1.5%；单对流环结构正确。
  - 算例在 /tmp/natconv（含 gen_cas.py 结构化 N×N×1 CAS 生成器）。
- [x] **多孔介质自然对流验证（2026-10-04，无源码改动）**：
  - Darcy 侧壁加热方腔（eps=1，Da=1e-6）：Ra_K=10/100/1000 → Nu=1.078/3.096/
    13.405(160²)，外推 13.45；文献 1.07/3.10/13.64（Mahmud&Fraser），偏差 −1.4%；
    40→80→160 二阶网格收敛，热平衡误差 0.004%。
  - 次临界多孔 RB（Ra_K=10<4π²）Nu=0.996 无流动；超临界 Ra_K=100 Nu=2.13 成环。
  - LTE k_eff 解析解 3.407 vs 3.4；eps 权重严格等价检验差 1e-5；
    Forchheimer 强阻力使 Nu 3.096→2.76 单调抑制。
  - 待办：1D LTNE 非平衡导热验证（用户计划最后一项）。
  - 注意：/tmp 会被沙箱会话间清空；算例归档现位于项目内 cases/natconv
    （2026-10-04 由 ~/cases 迁入）。
- [x] **1D LTNE 发汗冷却验证（2026-10-04，LTNE 定量验证收尾）**：
  - 算例 cases/ltne/：空气 2 kg/m²/s、300K 穿多孔 SS304 板（ε=0.3），出口端面
    恒热流 2e6 W/m²（tbc 于压力出口 zone）；h·a=0.3 使交换长度 4.9dx 可分辨。
  - 结果 vs 闭式解析解：能量守恒 ΔT_f(L)=995.02K（−0.001%）、T_s max 0.074%、
    T_f 中位 0.198%、Darcy 压降 5.17e5 vs 5.23e5 Pa；3343 步收敛。
  - 归档 README/control/log/vtu/plot_ltne.py/images。
- [x] **求解器 4 个 bug 修复（2026-10-04，LTNE 调试暴露）**：
  - PPE 钉死 cell 1 假质量汇：有压力 Dirichlet 面时不再钉死（仅全 Neumann 钉死）；
    串行 ppe_assembly + MPI ppe_assembly_mpi（has_pdir allreduce MPI_LOR）。
  - POUTLET 倒流削弱动量对角：改 max(F,0) 隐式 + min(F,0)·u_P 显式。
  - flux_ref 漏 massinlet uspeed（6 处 uscale 补 gb%uspeed）。
  - **MPI 本地网格丢 cell-zone**（重大）：build_local_mesh 未拷 czone/cztype/czt
    → MPI 多孔/LTNE 物理整体缺失（此前 porous np2 实为纯流体解）；补拷后
    porous_Ra10 np2 = 748 步与串行一致（修复前 296 步无 Darcy）。
  - build_bc 增 allow_empty（MPI 本地空 bc zone 合法）。
  - 回归：porous_Ra10 串行位级一致×2；cavity np2 正常；LTNE np2 与串行一致。
    已知限制：纯流体 mdot-inlet+出口无阻力时启动发散（V/apc 放大），后续可做
    BC 渐启/瞬态启动。
- [x] **算例归档 + 测试报告（2026-10-04，用户约定的长期模式）**：
  - <项目>/cases/natconv 完整归档：README.md 测试报告（设置/结果表/热平衡/复现命令/结论）、
    fluid/ 6 算例（重建并重跑，结果逐位复现）、porous_darcy/、porous_extra/、
    images/ 3 张图（温度+流场云图、Nu 对比曲线）、gen_cas.py/plot_fields.py/plot_nu.py。
  - 归档约定已写入 trae 项目记忆，后续所有验证算例照此办理。
- [x] **pressure-far-field BC 零解 bug 修复（2026-10-03）**：
  bc_face_vel 的入/出流判据由内场 uP·n 改为外部自由流 u_far·n；
  cylinder 从零解恢复物理解（CD=1.62、分离角±53°、Cp 近基线）；cavity/vinlet 回归位级一致。
- [x] **阶段 3 步骤 C（2026-10-03 完成）真实网格 grid_BC 贯通 + 体区域落位**：
  - [x] 结构化源码 Mesh3d.dat→Mesh3d.x（4 文件 12 处）；Mesh3d.x 二进制用 Mesh_File_Format=0
  - [x] mesh_t 增 czone/nczone/czt/cztype；sec_cells 存体区域号；finalize 拆面/体区域表
  - [x] resolve_cell_zones（id/name 匹配 fluid/porous，默认 fluid，未知报错）；两主程序接入
  - [x] build_bc 豁免 interface zone（fgrp=0，阶段 5 交换）
  - [x] 新 BC mass-flow-inlet mdot（BC_MASSINLET=7，uf=-(mdot/rho)n_out）
  - [x] 真实网格两侧读取成功；界面几何精确重合（y=0,x420-670,z0-50,250 quad,法向相反）；
    结构 Ma=3 冒烟正常；cavity/vinlet/cavity-np2 回归通过
  - [x] 物理参数已确认（2026-10-03）：uns=300K 空气 rho=1.177/mu=1.846e-5/mdot=2（u_in=1.70m/s）；
    struct 层流/AoA=0/Ref_L=1mm；根目录重复文件已删
  - [ ] 遗留：mm→m 面积/体积缩放在阶段 5 units 层统一；多孔动量汇（后续）

### 阶段 4：交界面定义与匹配（P1，约 3-4 天）
- [x] **阶段 4（2026-10-04 完成）**：
  - [x] Interface_FACE_TYPE 增 nv/verts(3,nv)/peer_w（面顶点+对侧权重）
  - [x] 结构侧 register_bc_interfaces 改逐单元面登记+4 顶点 CCW（grid_BC: 1→250 条）
  - [x] 非结构侧 register_interface_zones 补每条面顶点坐标
  - [x] 新建 mod_interface_match.f90：投影+Newton 双线性反演+双向 peer_id+权重
  - [x] 新建 test_match.f90 + Makefile match_test 目标（struct+uns+coupling 联合构建）
  - [x] grid_BC 验证：250/250 双 100%、投影距 1e-38、|Σw−1|=0
  - [x] 回归：cavity/cylinder 串行位级一致，cavity np2 仅 epsilon

### 阶段 5：数据交换与量纲统一（P1，约 2-3 天）
- [x] **mod_reference_state.f90**（2026-10-05）：运行期可配参考量（rho_ref/T_ref/L_ref/u_ref + 派生 a_ref/p_ref/p_scale），默认海平面大气；`read_mix_control` 解析 mix.control（key=value，缺省不报错）。
- [x] **mod_interface_units.f90**（2026-10-05）：结构无量纲↔SI 互转（struct_to_SI/SI_to_struct + 速度向量版本），约定 ρ*=ρ/ρ_ref, u*=u/a_ref, T*=T/T_ref, p*=p/(ρ_ref·a_ref²)。
- [x] **mod_interface_exchange.f90**（2026-10-05）：iface_state_t 容器 + struct→uns 双线性插值（peer_w）+ uns→struct 面积加权平均（守恒）。MPI 通信包装留阶段6。
- [x] **单元测试 test_units_exchange.f90** + Makefile `units_test` 目标：3 测试全 PASS（量纲往返 1e-14、常量场插值精确、面积加权正确）。
- [x] MPI 子通信域 + 跨域 Send/Recv 包装（2026-10-06 勾销：阶段 6 已实现，
  `src/coupling/mod_coupling_exchange.f90`，count+payload 协议 tags 110/111/210/211，
  支持两侧面数不同与 nfaces=0）

### 阶段 6：弱耦合驱动（P1，约 2 天）
- [x] **main.f90 耦合驱动框架（2026-10-05）**：
  - [x] MPI_Init + `MPI_Comm_split` 拆 STRUCT_GROUP/UNS_GROUP 子通信域
  - [x] 读 mix.control（read_mix_control）
  - [x] 非结构侧：`mod_uns_driver.f90`（init/step/extract_iface）
    + `BC_INTERFACE=9` + `set_interface_vel/p` + simple_run nsteps 上限
  - [x] Makefile：BOOTSTRAP_ORDER 扩展 + UNS_LAYER6_S/9_M (mod_uns_driver)
  - [x] 冒烟测试通过：cavity 2 ranks × 3 耦合迭代 exit=0；
    cavity 回归 697 步收敛无破坏
- [x] **结构侧接入 + MPI 跨组交换（2026-10-05）**：
  - [x] `mod_struct_driver.f90`：init(comm,ctl)/step/extract_iface/
    set_iface_bc（守恒↔原始量互转，ghost cell 直写 B%U）
  - [x] `mod_coupling_exchange.f90`：count+payload 协议跨组 Send/Recv
    （tags 110/111/210/211，支持两侧面数不同与 nfaces=0）
  - [x] mod_interface 增 ic/jc/kc + ig/jg/kg；register_bc_interfaces 填充
  - [x] Global_Var 增 `Struct_Comm`（默认 COMM_WORLD）；6 个 struct 库文件
    64 处 `MPI_COMM_WORLD→Struct_Comm`（修复 np=2 struct init 挂起）
  - [x] main.f90：`has_struct` 保护修复 nproc=1 uns 自等 Recv 死锁
  - [x] 验证：grid_BC np=2 × 3 耦合迭代 exit=0、双向各 250 面、struct step
    正常；np=1 不挂起；standalone struct 回归一致
- [x] 界面插值接入驱动（2026-10-06 勾销：`src/coupling/mod_interface_exchange.f90`
  已实现 struct→uns `peer_w` 双线性、uns→struct 面面积加权平均（守恒），
  已在流 B / C2 耦合链路中实际使用；"简单平均占位"描述已过期）
- [ ] 界面收敛检查（**重新界定**：现状 = 固定 `n_couple` 迭代 + 每 25 步打印
  iface 均值（`src/main.f90`），缺"界面跳变 < tol 早停 / 界面残差历史"；
  动手前需先与用户定判据，故保留未勾销）
- [x] 修复 uns_solver_extract_iface 的 rho 占位（2026-10-06 勾销：
  `mod_uns_driver.f90` 增 optional `rho_in`，`src/main.f90` 调用处已传 `ctrl%rho`）
- [x] 从 mix.control 读 n_couple / n_uns_steps（2026-10-06 勾销：
  `src/main.f90` 用 `get_coupling_params(n_couple, n_uns_steps, iface_relax,
  iface_ramp, …)`）
- [x] unMesh.cas mm→m 单位缩放（2026-10-06 勾销：`mod_uns_control.read_mesh_scale`
  + `ctrl%mesh_scale`，默认 1.0；C2 算例用 `mesh_scale = 1.0e-3`）
- [x] grid_BC 真实耦合物理验证（2026-10-06 勾销：grid_BC np=2 双向各 250 面 +
  后续 `cases/couple_channel`（流 B）与 `cases/couple_porous`（C2）真实物理验收）

### 阶段 7：构建与测试（P2）
- [x] 完善 Makefile 依赖；单元测试（量纲转换、匹配精度、交换守恒性）（2026-10-06 勾销：
  `bin/units_test`（量纲往返/常量场插值/面积加权）、`bin/match_test`（匹配精度）、
  `bin/coupling_test`（np2 交换）三目标已就位且 rc=0；pristine/BOOTSTRAP 依赖链已建。
  "继续完善依赖"并入日常维护，不再单列）

### 阶段 8：算例验证（P2）
- [~] 亚声速通道（结构）+ 低速腔体（非结构）；界面物理量偏差 < 1%（2026-10-06 重界定：
  两侧算例均已在跑（`cases/couple_channel` = 流 B / `cases/natconv` 等），
  但实测界面速度偏差 **−2.9%**（流 B）/ 型线 11.8%（C2 §4），**尚未达到 <1% 判据**，
  已归因于两侧**壁面层物理差异**（非交换层误差）⇒ 与阶段 11 可选项"壁面层对齐"
  联动；需与用户重定判据或收口归因后方可勾销）

### 文档（随功能同步）
- [x] `docs/程序使用手册.tex` 及 `91_appendix_params.tex`、`93_changelog.tex`、`control.ec.template`
      **（2026-10-08 全部完成；23 页 PDF 编译通过，见「已完成」顶部文档节点）**
- [ ] **`cases/couple_channel/README.md` 缺失**（2026-10-06 新登记）：流 B 的复现命令、
  界面定值（质量精确连续、压力跳变 <1 Pa、对拍 uns_full 剩余 −2.9%）只散在
  `docs/plan.md` 流 B 节与 worklog 里；对照 `cases/couple_porous/README.md` 的
  体例补齐，作为算例文档一致性收口

### 阶段 9–13：后续功能储备（2026-10-05 录入 docs/plan.md，源自程序功能说明 L8–15）
- [~] **阶段9 多孔介质完善**：
  - [x] CAS zone 名 VC:porous/fluid 自动识别（2026-10-05）：优先级
    显式 cell_zone 类型词 > VC 标记 > fluid 默认；cell_zone 类型词可省略
    （CZ_AUTO 只给系数）；zone 名字段 32→128 字符；多孔区缺系数告警。
    涉及 mod_uns_control（CZ_AUTO/vc_tag_type/resolve/parse）、
    mod_uns_mesh（名字段加宽）、mod_uns_fields（AUTO 系数行生效+修 h_sf 旧笔误）。
  - [x] **porous-plug 验证（2026-10-05）cases/porous_plug/**：40×2×1 一维管，
    全 symmetry；Darcy Δp 46.15Pa 误差 0.000%、Darcy-Forch 520.35Pa 误差 0.002%；
    内部 u 偏差<0.02%；串行/np2 同为 61 步。Ra10 回归 Nu=1.078 一致。
  - [x] **热弥散混合层验证（2026-10-05）cases/porous_disp/**（见顶部第二批）。
  - [x] **二维多孔介质对流（2026-10-05）cases/porous_graetz/**：强制对流
    Graetz（见顶部第三批）；自然对流多孔方腔此前已在 cases/natconv/ 完成。
  - [x] **Beavers-Joseph（2026-10-05）cases/beavers_joseph/**（见顶部收尾条）。
  - [x] **α/K 标定实验设计（2026-10-06）cases/calibration/**（见顶部条目）。
  - 已完成历史：Darcy-Forchheimer、Brinkman、LTE/LTNE；1D LTNE ✅、多孔自然对流 ✅。
- [x] ~~**阶段10 自动保存与耦合联合重启**：两侧流场（含 halo）定时保存；struct 按节点
  插值写 Plot3D；main.f90 联合重启入口。~~（已完成并勾销 2026-10-06：见本文件
  「已完成」区阶段 10 条目——`write_field_dump`/`read_field_dump` +
  `couple_restart=1` 逐位续跑，本条为**陈旧残留**，故删除勾选）
- [x] **阶段11 多类型界面**：低速-多孔（BJ 滑移）✅ C1、低速可压缩 ✅ 流 B、
  可压缩-多孔 ✅ C2、**界面类型/交换量自动分派表 ✅（2026-10-07）**、
  **`uns` 绝对压力 ≈−250 Pa 偏置修复 ✅（2026-10-07）** **均已完成**（详见上文
  两条 ＋ README §5.3）。
- [ ] **阶段12 Gambit NEU 输入**：.neu 读取器（网格+BC+体区域属性），复用现有登记流程。
- [ ] **阶段13 结构求解器演进（远期）**：SST bug 修复前禁用；改用 Liao 格心型有限差分，兼容现有 Riemann。
