# 项目简报 (projectbrief)

## 项目名称
MixNSSolver —— 结构网格 + 非结构网格混合 CFD 求解器

## 项目目标
将两个已有的独立求解器集成为一个统一的弱耦合混合求解器：
- **结构网格部分**：基于 OpenCFD-EC v1.16a，采用 Riemann 求解器，面向可压缩流体计算。
- **非结构网格部分**：基于已完成的 UNSSolverProj，采用 SIMPLE / PIMPLE 算法，面向低速流与多孔介质流动。
- **耦合方式**：弱耦合（分步迭代），在交界面互相传递边界条件。

## 来源依据
- 需求：`docs/程序功能说明.md`（用户编写的 6 条主体框架说明）
- 计划：`docs/plan.md`（8 个实施阶段）
- 专题规划：`.trae/documents/mpi_parallel_plan.md`、`.trae/documents/temperature_equation_plan.md`（尚未细读）

## 总体架构
```
src/
├── common/        公共模块（精度、常数、量纲转换）
├── structured/    结构求解器（源自 OpenCFD-EC，模块化拆分）
├── unstructured/  非结构求解器（源自 UNSSolverProj）
├── coupling/      交界面定义/匹配/交换/量纲统一
└── main.f90       耦合主驱动
```

## 耦合与网格匹配要点
- 结构与非结构网格**分开输入**，通过 interface 类型交界面匹配。
- 结构侧交界面：Gridgen `.inp` 中 `gridgen generic: 8` 类型。
- 非结构侧交界面：Fluent CAS 文件中的 interface 边界。
- 需建立相同交界面在两套网格中所属单元（面）的对应关系。
- 非结构网格不同块的类型（低速/多孔介质等）通过块名称或单独的块属性控制文件指定。

## 并行策略
- 结构网格：沿用 OpenCFD-EC 原有并行方式，不更改。
- 非结构网格：ParMETIS/METIS 分块（已在 UNSSolverProj 实现）。
- 非结构网格整合后，通过存储的界面关系为两个块设置界面边界条件。
- 两求解器共享 MPI_COMM_WORLD，可用 MPI_Comm_split 拆子通信域。

## 输出约定（国际单位制）
- 结构网格：Plot3D 流场文件（rho, u, v, w, t）。
- 非结构网格：VTK 格式（rho, u, v, w, p, tf, ts）。
- 定时保存重启文件。

## 头等风险
- **量纲一致性**：OpenCFD-EC 内部为无量纲变量，UNSSolverProj 为 SI；交界面必须先各自转 SI、交换、再各自转回。

## 技术环境
- 语言：Fortran 2008
- 编译器：gcc 11 / gfortran 8.3.0（用户说明）
- MPI + METIS/ParMETIS（预编译库位于 `lib/`）
