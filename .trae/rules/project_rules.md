# Cline 项目记忆规则

## 每次新任务开始时，必须执行
在分析代码、提出方案或修改文件之前，先按顺序读取：
1. memory-bank/projectbrief.md
2. memory-bank/activeContext.md
3. memory-bank/progress.md
4. memory-bank/worklog.md（只读最近 5 条）
5. memory-bank/constraints.md

读完后，先用不超过 10 行总结：
- 当前项目目标
- 已完成
- 当前任务
- 关键约束
- 建议下一步

然后问我是否继续。不要直接改代码。

## 禁止
- 未读完上述文件前，不要重新扫描整个项目。
- 不要凭猜测补充缺失信息，缺失就问我。
- 除非我明确要求，否则不要重新通读整个代码库。

## 每次任务结束或关键决策后，必须更新
- memory-bank/activeContext.md：当前任务、最近决策、下一步
- memory-bank/progress.md：已完成 / 待办
- memory-bank/worklog.md：追加一条简要记录

## 程序说明文档（LaTeX）维护规则
程序每新增/修改一个功能，必须同步更新 `docs/程序使用手册.tex：
1. 更新对应章节（物理模型、控制参数表、示例）；
2. 新增参数时同步附录 A 参数总表（`91_appendix_params.tex`）与 `control.ec.template`；
3. 在 `93_changelog.tex` 追加一条记录（日期/版本/功能/涉及源码/验证）；
4. 重新编译：新增功能以及控制参数含义，应更新到程序使用手册.tex。

## 验证/诊断硬规则（2026-10-06，源自 C2 可压缩–多孔跨组界面验证）
- **禁止跨求解器比对绝对压力**：uns 单求解器的绝对压力水平带 ≈−250 Pa 的内部
  常数偏置（内部整段下移、紧邻入口边界那一层保持物理解水平、出口末列回收；与
  多孔无关、与出口 BC 类型无关、严格 ∝u²）。跨组对比只比"界面连续性 ＋ 梯度/
  型线"；要用绝对压力必须先在 `src/unstructured/` 投影步与 BC 层修掉该项。
  取证：`cases/couple_porous/README.md` §5（四组对照 ＋ §6.1 复现配方）。
- **取窗/拟合窗宽度必须 ≤ 一个网格单元**（用绝对窗 `[x0, x0+dx)`，不要"中心
  ± 大于 dx/2"）：2.5 mm 窗跨了两个 dx=2 mm 平面，把 x=101 mm 的 295 Pa 与
  x=103 mm 的 70 Pa 平均成 182.6 Pa，造出幻影 112 Pa 界面跳变（C2 缺陷②）。
- **多 cell zone 网格的 cell-id 排序：分割方向索引必须最慢**，使每个 zone 对应
  一段**连续** id 区间（Fluent 的 cell-zone 记录就是单段连续区间）。`1+i+j*NX+
  k*NX*NY`（k 最慢）会把"x 串联"的 fluid/porous 写成了"z 并联"的两条全长薄片，
  症状：梯度只有解析值 ~60%、max|u| 超入口、zone 互换后结果**逐位相同**。
- **复现耦合算例必须设 `save_interval ≤ n_couple`**：`flow3d.dat` 只在周期存档
  处落盘，仓库默认 `n_couple=1000 / save_interval=1000` 时跑 400 迭代**什么都
  不写**。完整配方见 `cases/couple_porous/README.md` §6。
- **验证脚本必须自带回归守卫**：在错误（历史缺陷）资产上要能主动报错，而不是
  静默给出"看起来合理"的数，例如 `uns_full/check_bed_gradient.py` 的
  `zone 2 is not the x<100 mm half!`。

## 版本控制约定（2026-10-06 起，仓库已 `git init`）
- 分支 `main`；改完一个**可验证的节点**就提交一次，提交前先跑该节点对应的验证脚本
  （或位级回归 `regress/m6wing/run_regression.sh`）。
- `.gitignore` 只忽略构建/中间产物：`build/`、`bin/`、`*.o`、`*.mod`、
  `__pycache__/`、`*.pyc`、编辑器临时文件。**不得**忽略 `lib/*/*.a`
  （vendored METIS/ParMETIS/Tecplot 预编译库，链接必需且不可重建）与
  `regress/m6wing/baseline/*`（位级回归 oracle）。
- 提交信息格式：`<type>: <一句话>`（type ∈ feat / fix / docs / verify / chore），
  正文写清「验证方式 / 关键结果 / 遗留待办」，并与 `memory-bank/worklog.md`
  对应条目一致；本仓库文档为中文，提交信息用中文。
- 回退：单文件 `git checkout <sha> -- <path>`；未提交改动 `git stash`；
  无远端时可用 `commit --amend` 重写单个未共享提交。
