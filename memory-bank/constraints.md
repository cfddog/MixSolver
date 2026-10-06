# 关键约束 (constraints)

> 来源：`docs/程序功能说明.md`、`docs/plan.md`、项目规则 `.trae/rules/project_rules.md`、用户偏好。

## 1. 语言与编译
- Fortran 2008；目标编译器 gfortran 8.3.0（gcc 11 环境）。
- 统一使用 `mod_precision` 中的 `dp`；OpenCFD-EC 旧 `PRE_EC=8` 通过别名 `PRE_EC=dp` 兼容。
- MPI 代码用 `use mpi`，不再用 `include "mpif.h"`。

## 2. 量纲一致性（最高优先级，需求第 6 条）
- OpenCFD-EC 内部：无量纲（rho*=rho/rho_ref，u*=u/a_ref，T*=T/T_ref，p*=p/(rho_ref·a_ref²)）。
- UNSSolverProj 内部：SI（kg/m³, m/s, K, Pa）。
- 交界面交换铁律：**各自先转 SI → 交换 → 再各自转回**；转换集中在 `mod_interface_units.f90`。
- 界面传递方向（plan 阶段 5）：
  - 结构 → 非结构：rho, u, v, w, T（作非结构入口/远场边界）
  - 非结构 → 结构：rho, u, v, w, p（作结构出口/远场边界）
- 验证标准：界面两侧物理量偏差 < 1%。

## 3. 网格与交界面
- 结构与非结构网格分开输入。
- 结构侧交界面：Gridgen 格式中的 `gridgen generic: 8` 类型。
- 非结构侧交界面：Fluent CAS 文件的 interface 边界。
- 同名交界面配对；需建立面/单元对应关系（最近点/投影 + 双线性插值权重，双向映射）。
- 数据交换应支持保守插值，保证质量/动量/能量通量守恒。

## 4. 求解器适用范围
- 结构网格固定为可压缩 Riemann 求解器，不区分类型。
- 非结构网格按块区分（低速 SIMPLE / PIMPLE、多孔介质等），通过块名称或单独块属性控制文件指定。

## 5. 并行
- 结构网格并行沿用 OpenCFD-EC 原有方式，不更改。
- 非结构网格用 ParMETIS/METIS 分块。
- 两求解器共享 MPI_COMM_WORLD；可 MPI_Comm_split 建结构组/非结构组子域，跨组点对点交换。

## 6. 输出
- 一律国际单位制。
- 结构：Plot3D（rho, u, v, w, t）。
- 非结构：VTK（rho, u, v, w, p, tf, ts）。
- 定时保存重启文件。

## 7. 工作流规则（.trae/rules/project_rules.md）
- 每次新任务开始先读 memory-bank 五个文件（worklog 只读最近 5 条），先摘要汇报并询问是否继续，不直接改代码。
- 不凭猜测补信息，缺失就问用户；未被要求不通读整个代码库。
- 任务结束/关键决策后必须更新 activeContext.md、progress.md，并向 worklog.md 追加记录。
- 每新增/修改功能必须同步 `docs/程序使用手册.tex`：对应章节、附录 A 参数总表（`91_appendix_params.tex`）、`control.ec.template`、`93_changelog.tex` 追加变更记录，并重新编译。

## 8. 代码与注释风格（用户偏好）
- 所有新增注释一律用英文，防止乱码。
- 已有乱码注释直接删除；已有中文注释保留，待后续修改该代码时再转英文。

## 9. 已定稿的设计决策（2026-10-03 用户拍板）
- **非结构块类型**：不新建文件，扩展现有 `.control`（key=value 风格），新增
  `cell_zone = <id 或 name> fluid|porous [perm=.. inertial=.. porosity=..]`。
  阶段 3 需先给 CAS reader/mesh 补 cell zone（体区域）解析（现有代码只有 face zone）。
- **参考量**：mod_constants 只留普适常量；参考状态拆到运行期模块（mod_reference_state），
  默认海平面标准大气，启动时从耦合控制文件（拟名 `mix.control`，与 .control 同风格）读入覆盖。

## 10. 外部源码位置
- UNSSolverProj：`/home/sundong/Fortran_Project/UNSSolverProj`（独立 git 仓库，含 Makefile/cases/tests/lib，23 个 src 模块）。
- **结构求解器迁移唯一来源**：项目内 `external/OpenCFD-EC-1.16a/`（根目录 zip 解压的 1.16a 原始版，**只读基线，禁止改动**）。
- `/home/sundong/Fortran_Project/OpenCFD-EC-1.16a` 是深度定制版（含 lowspeed/porous/multi_region），**不作为迁移来源**，仅可作功能参考。

## 11. 已知未知（禁止臆造，需向用户确认）
- 耦合控制文件 mix.control 的完整字段清单（阶段 5/6 前确定即可）。
- 多孔介质具体阻力模型形式（达西/Brinkman-Forchheimer 系数约定），阶段 3 实现 porous 时确认。

## 12. 构建系统约定（2026-10-03 确立，顶层 Makefile）
- 目标：`make all/mpi`（混合）、`make structured(_mpi)`、`make unstructured(_mpi)`、
  `make common(_mpi)`、`make clean`、`make help`；`DEBUG=1` 开调试选项。
- 产物：`build/ser/...` 与 `build/mpi/...` 双树（.o/.d 镜像 src，.mod 集中 mod/），可执行文件在 `bin/`；
  串行/MPI 切换无需 clean；禁止在源码目录或根目录留 .mod/.o。
- 模块名**全局唯一**（共用一个 -J 目录）：迁移时一律用 mod_struct_*/mod_uns_* 前缀，
  删除两个源求解器自带的 mod_precision，统一 use common 层。
- MPI 专用源文件命名为 `*_mpi.f90`（串行构建自动排除）；共享文件中 MPI 代码用 `#ifdef HAVE_MPI` 保护。
- 依赖由 `gfortran -cpp -MMD -MP` 自动生成；干净树大批量新增模块后首次构建建议 `make -j1`，
  或在 Makefile "PRISTINE-BUILD EDGES" 段补显式边；基础模块加入 BOOTSTRAP_ORDER。
- 分块链接：默认串行 METIS 在 rank0 分块 + TECIO；libparmetis_gnu_mpi.a 为 MPICH ABI，
  与本机 OpenMPI 4.1.6 不兼容，重新编译前不要链接。
- 编译标准：`-cpp -std=f2008 -Wall`（本机 gfortran 13，需兼容用户的 gfortran 8.3.0）。

## 13. structured 树的例外约定（2026-10-03 阶段 2a 确立，2b 重新评估）
- structured 求解器 MPI-only：ser/mpi 两树**均用 mpif90** 编译链接；ser 树即"单进程 MPI 构建"
  （运行用 mpirun -np 1），不定义 HAVE_MPI。`%_mpi.f90` 串行排除规则仅适用于 unstructured。
- structured 树**单独**追加 `-std=legacy -ffree-line-length-none`（旧代码有 pointer 数组元素
  传显式 shape 哑参、REAL/INTEGER 哑参类型不匹配、>132 列长行）；common 与新代码仍守 f2008。
- precision_EC 薄壳（mod_struct_precision.f90）：`use mpi` 实体保持 public（等价原模块内
  include mpif.h 经 use 链传递 MPI 常量）；PRE_EC 用 selected_real_kind(15,307) 自含定义，
  **不得 import/重导出 common 的 dp**（旧例程存在名为 dp 的局部变量，会撞名）。
- 已知原版遗留问题（2a 忠实搬迁未修，留 2b）：程序结束不调 MPI_Finalize（mpirun RC=1 属正常）；
  sub_read_parameter.f90 的 a0/d0 仅涡轮机械分支赋值、普通算例未初始化即打印（不参与计算）。
- 回归基线做法：原始源拷贝到 /tmp（不在 external/ 内编译），iconv GBK→UTF-8、全角"！"转半角、
  mpif90 -O3 -std=legacy -ffree-line-length-none；回归比较用 cmp 二进制流场 + 文本 diff。
