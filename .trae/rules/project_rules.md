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
- **git 提交（自动，无需询问）**：`git add -A && git commit -F -`，与上面三项同批提交，
  见下「版本控制约定」。

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
- 分支 `main`。
- **自动提交（默认行为，不必再问我）**：每完成一个**小节点**（一次可验证的改动 /
  缺陷修复 / 文档收尾）就走完「验证 → 更新 memory-bank → 提交」三步，直接执行
  `git add -A && git commit -F -`。**一个节点 = 一个提交**：不把互不相关的改动塞进
  同一提交，也不攒多个节点一次性提交。
  - 验证：跑该节点对应的验证脚本；动过解算器/耦合/网格代码时加跑位级回归
    `regress/m6wing/run_regression.sh`。
  - 记账：`memory-bank/{activeContext,progress,worklog}.md` 的更新与提交**同批**。
  - 只有「验证通过且状态可复现」才提交；验证未通过或半成品**不提交**，
    保持脏工作区或 `git stash` 并在 worklog 写明停在哪一步。
  - 算例的重算产物（`flow3d.dat`、`unMesh_coupled.vtu` 等）默认落在运行目录，
    提交前确认没有把一次性产物误加入仓库。
- `.gitignore`：忽略构建/中间产物（`build/`、`bin/`、`*.o`、`*.mod`、
  `__pycache__/`、`*.pyc`、编辑器临时文件）＋ **全部运行产物**
  （`*.vtu`、`*.log`、`*.part.map`、`*.tmp`、`cases/**/*.dat`、`grid_BC/**/*.dat`、
  `cases/**/output_para.out`、`checkpoint_*.tar.gz`）。
  **不得**忽略 `lib/*/*.a`（vendored METIS/ParMETIS/Tecplot 预编译库，链接必需
  且不可重建）、`regress/m6wing/baseline/*`（位级回归 oracle，含 `*.dat`！所以
  `.dat` 规则**只**写 `cases/` 与 `grid_BC/` 两棵树，禁写全局 `*.dat`）与
  `cases/*/**` 的输入件（`*.cas/*.neu/*.cgns/*.x/*.control/mix.control/bc3d.*`、
  工具 `*.py`、`README.md`、`images/*.png`、`scratch/unMesh.cas.zsplit_bug`）。
- **运行产物不入库**（2026-10-06 起，硬的）：`.vtu/*.log/*.dat/*.out/*.tmp/
  *.part.map/__pycache__` 等由求解器写出的文件**一律不提交**，无论大小。
  提交前 `git status` 里若出现这类文件即为误加，用 `git restore --staged` 退回，
  并确认 `.gitignore` 覆盖到位。需要留证据时把**关键数值/命令**写进
  `cases/*/README.md` 与 memory-bank，而不是把产物入库。
- 提交信息格式：`<type>: <一句话>`（type ∈ feat / fix / docs / verify / chore），
  正文写清「验证方式 / 关键结果 / 遗留待办」，并与 `memory-bank/worklog.md`
  对应条目一致；本仓库文档为中文，提交信息用中文。
- **远端同步**：`origin = git@github.com:cfddog/MixSolver.git`（SSH，密钥
  `~/.ssh/id_ed25519`）。**提交后自动 `git push origin main`**，与自动提交同一节点内
  完成，不单独询问。推送前 `git status -sb` 确认与 `origin/main` 无分叉（有分叉先
  `git pull --rebase`）；对已共享历史**禁止 `--force`**。网络不可达时跳过推送，并在
  worklog 注明"本节点未推送"。
- 仓库体量（2026-10-06 **瘦身后**）：跟踪文件 **324** 个（原 536），源码 + vendored
  库 + `external/` 原始包 + 算例**输入件**；212 个运行产物 / 761.8 MB 已用
  `git filter-branch` 从全部历史移除（`.git` 120 MB → 见 worklog 校验条目），
  工作树产物留在磁盘且被忽略。
- **一次性历史重写（已执行，2026-10-06）**：为剥离运行产物，在私有仓库（无他人
  clone，用户已确认）做过一次 `git push --force-with-lease=main:<远端旧值> origin main`
  （`3fa7b91 → f2ebf0b`，**零告警**，`.git` 120 MB → 46 MB）。
  安全网 = **仓库外 bundle**：
  `/home/sundong/mixsolver_pre_slim_backup/repo_pre_slim.bundle`（98 MB，
  `git bundle verify` = 完整旧历史）+ `cases_hardlinks/`、`grid_BC_hardlinks/` 硬链接快照。
  **本地分支不能当安全网**：`filter-branch --all` 会连备份分支与 `refs/remotes/*`
  一起重写（本次 `backup/pre-slim-2026-10-06` 就被一并改写，随后删除）。
  重写后推送前**必须** `git fetch origin` 复位远端跟踪引用，否则
  `--force-with-lease` 会拿被改写的本地缓存去比对而误判。
  **注意**：这次 `fetch` 会把旧远端对象拉回本地（实测 `.git` 46 → 140 MB），
  推送成功后要**再跑一次** `git reflog expire --expire=now --all &&
  git gc --prune=now` 回收（实测回到 46 MB / 1 pack / 434 对象）。
  **此后恢复铁律：对已共享历史禁止 `--force`/`--force-with-lease`**；
  再要重写必须先确认无他人 clone 并重新做 bundle 备份。
- 回退：单文件 `git checkout <sha> -- <path>`；未提交改动 `git stash`；
  历史重写后的旧历史只能从 bundle 找回：
  `git fetch /home/sundong/mixsolver_pre_slim_backup/repo_pre_slim.bundle
  'refs/heads/main:refs/heads/restored-pre-slim'`。
  远端已有对应提交后不得用 `commit --amend` 重写。
