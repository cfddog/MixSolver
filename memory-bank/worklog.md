# 工作日志 (worklog)

> 规则：每次任务结束或关键决策后追加一条简要记录（日期 | 内容 | 涉及文件）。
> 最近记录在最上方。

---

## 2026-10-07 | 阶段 11 补充：`pval`（出口压力）取值敏感性实测 + 文档纠错

- **背景**：追问 `pressure-outlet` 的 `pval`（表压）能否取 100/1000/10000。
- **实测**（全流体 baseline、默认 `inlet_ramp=100`）：ρ 为常数 ⇒ `pval` 只是**表压锚点**，
  收敛后整场压力均匀平移、速度/梯度不变（`+15`→p≡15.000、`−100`→p≡−100.000）。
  但 **`pval` 有收敛上限且正负不对称**：正值鲁棒区仅 ≈+15 Pa、负值到 ≈−200 Pa；
  **+20/+22/+100/+1000/+10000 一律发散**（`+100` 复现 5/5；`−100` 复现 5/5 收敛）；
  慢化 `inlet_ramp` 不救 `+100`（ramp=300/1000/3000 仍发散）。
- **机理线索**：失稳首现于**温度**（`pval=100` 首迭代 Tmin 300→≈114 K 后 NaN），而温度
  方程无显式 p 依赖 ⇒ 暴露**压力电平的非规范不变性**（与早期 −250 Pa 偏置同族），待立项。
- **文档**：README §5.3 表行纠错（原“±100 亦收敛”不成立）＋ 新增 §5.3.1。
- 涉及文件：`cases/couple_porous/README.md`、`memory-bank/{worklog,activeContext}.md`。

---

## 2026-10-07 | 阶段 11：uns 绝对压力 ≈−250 Pa 偏置修复（根因＝对流系数重复乘 ρ）＋ 出口压力 ramp

- **目标**：修掉 uns 全场绝对压力 ≈−251 Pa 偏置（机理已定位），并顺带修出口
  `pressure-outlet` 改值即发散。
- **根因**：`mod_uns_simple.f90` `momentum_assembly` 内部面对流系数
  `F = ctrl%rho*fld%flux(i)` **重复乘 ρ**——`fld%flux` 已是**质量通量 kg/s**
  （PPE `rhs -= flux`、温度装配 `FT = cp*flux`、边界面 `F = ρ*(u_f·S_f)` 三处佐证）。
  内部对流被放大 ρ 倍、与边界面不一致 ⇒ 开放边界盖章的偏移由 O(u²) 放大为
  O((ρ−1)ρu²)，正是量到的系数 `0.176943 = ρ−1`（ρ=1.177）。**修法①**：
  `F = fld%flux(i)`。
- **修法②（出口鲁棒性）**：出口/远场压力 Dirichlet 原先冷启动即满值（入口 `inlet_ramp`
  ⇒ 首迭代 ap 极小 ⇒ `pval/(ρu²)`~70 ⇒ 第二步 du_max 71.5→6.5e5→NaN）⇒
  `mod_uns_bc` 新增模块级 `g_pval_ramp_factor`/`set_pval_ramp_factor`，`bc_face_p`
  对 `BC_POUTLET/BC_FARFIELD` 用 `pf = pval*g_pval_ramp_factor`；`mod_uns_simple(_mpi)`
  在 `inlet_ramp>1` 时与 `set_inlet_ramp_factor` 同值调用（ramp 完成因子=1，收敛解不变）。
- **验证**：baseline 平台 −251.17→**≈0 Pa**、`u₁/u_in` 0.9168→**1.000000**、
  床组平移 −251→**+0.009 Pa**（梯度 −3021.7 Pa/m 不变）、出口 `±10` **全部收敛**、
  c2/c3/nx50/nx200 平台 ≈0、m6wing 结构化位级回归 **PASS**（`flow3d.dat` md5
  `dc134a2d196422043ecad7c86ac8f898` 不变）、`units_test` 3/0、`coupling_test` np2 6/0；
  耦合首排（np2/400iter）−9.67%→**−1.51%**、参考 @x=100 p 48.5→**299.9 Pa**、
  界面跳变 −0.49 Pa；`make all`/`mpi` RC=0（改动仅 3 文件）。→ **阶段 11 全部完成。**
- 涉及文件：`src/unstructured/mod_uns_simple.f90`、`src/unstructured/mod_uns_simple_mpi.f90`、
  `src/unstructured/mod_uns_bc.f90`、`cases/couple_porous/README.md`、
  `memory-bank/{activeContext,progress,worklog}.md`。

---

## 2026-10-07 | 阶段 11：界面类型自动分派表（分类 + 报告；不改数值行为）

- **目标**：把 `main.f90` 里硬编码的界面处理抽象为「按两侧 solver + cell-zone
  类别自动判定界面类型并报告交换量清单」。
- **common**（`src/common/mod_interface.f90`）：`Interface_FACE_TYPE` 增
  `cz_type / loc_face / iface_type` 三字段；枚举 `IFACE_UNKNOWN /
  IFACE_COMP_FLUID_FLUID / IFACE_COMP_FLUID_POROUS / IFACE_UNS_FLUID_POROUS`、
  cell-zone 类别 `IFACE_CZ_FLUID/POROUS`、交换量 `Q_U/Q_T/Q_P/Q_RHO` +
  角色 `XQ_DIRICHLET/XQ_CHARACTER/XQ_JUMP`；纯过程 `classify_interface`、
  `iface_exchange_recipe`（**用 pure 子程序**，因 pure function 不允许 intent(out)
  实参）、`iface_type_name/quantity_name/role_name/iface_recipe_string`、
  `report_interface_dispatch`。
- **踩坑**：① `pure function` 内 `intent(out)` 数组非法 → 改子程序；
  ② 局部逻辑名 `peer_struct/peer_uns` 与常量 `PEER_STRUCT/PEER_UNS` **大小写
  不敏感同名**互相遮蔽 → 去掉逻辑变量改写。
- **uns**：`mod_uns_geometry.tag_interface_cell_zones`（由
  `m%cztype(m%f(loc_face)%c0)` 打 `cz_type`，须在 `resolve_cell_zones` 后调用），
  `register_interface_zones` 记 `loc_face`；`mod_uns_driver` 接入生产路径
  （`PEER_STRUCT / IFACE_CZ_FLUID`，匹配前打标）。
- **struct**：`mod_struct_grid` 把界面面标 fluid。
- **coupling**：`mod_interface_match.dispatch_interfaces`（匹配后分类、两侧回填
  `iface_type`、打印分派表）。
- **判定规则**：struct 侧恒 comp+fluid；跨求解器仅由 uns 侧 `cz_type` 决定
  （porous → `IFACE_COMP_FLUID_POROUS`）；uns↔uns → `IFACE_UNS_FLUID_POROUS`。
- **验证**：`match_test`(grid_BC) 250/250 全 `comp-fluid<->lowspeed-fluid`、
  断言 PASS；`couple_channel` np2 → fluid、`couple_porous` np2 → porous；
  `regress/m6wing` PASS（`flow3d.dat` md5 不变）；全量 make RC=0。
- ~~**仍待办（同阶段 11 未勾销）**：uns 绝对压力 ≈−250 Pa 偏置修复。~~ **2026-10-07 已完成**（见上一条）。
- 涉及文件：`src/common/mod_interface.f90`、`src/coupling/mod_interface_match.f90`、
  `src/coupling/test_match.f90`、`src/structured/mod_struct_grid.f90`、
  `src/unstructured/mod_uns_driver.f90`、`src/unstructured/mod_uns_geometry.f90`、
  `docs/plan.md`、`memory-bank/{activeContext,progress,worklog}.md`。

## 2026-10-06 | 并行小节点：uns 绝对压力 −251 Pa 偏置**定位完成**（整场电平/零模，非梯度误差）

- **复现（逐位）**：新建 `cases/couple_porous/uns_full/bias_scan.sh`（工作目录参数
  + 自动拷 `gen_full.py`/`plane_profile.py` + 生成 `sum.py`/`shift.py`）跑 baseline
  （全流体 + `pressure-outlet 0.0` + `mass-flow-inlet 40.868`，4000 格，**3 s/跑**）
  得到 `max p = −0.051499`、`min p = −285.763`，与 README 组 2 逐位一致；床组
  （`unMesh.control` 原样）得 `max p = 301.8786`、`min p = −248.9305`，与组 1 逐位一致。
- **全 100 平面剖面**（关键新证据）：
  - x = 1 mm（入口邻格）：`p = −0.0532`、`u = 31.83194`（= **0.91678·u_in**）
  - x = 3–9 mm：凹陷到 `p = −285.763 @5 mm`
  - x = 11 … 195 mm：**93 个平面逐位相同** `p = −251.16752`、`u = 34.72218`
    （与入口通量 `mdot/ρ = 34.7222` **逐位相等**）
  - x = 197 mm：`p = −257.81`、`u = 33.831`；x = 199 mm（出口单元）：`p = −158.82`
  ⇒ 内部 = **精确的"压力零模"**（均匀流 + 平坦压力，内部面只贡献差值 ⇒ 任意常数
  都满足动量方程）；电平由**入口侧印上**后原样传播，出口 Dirichlet **只钉住末列**。
- **差分矩阵（全部落 /tmp，仓库无产物）**：

  | 变量 | 平台 p | u₁/u_in | 结论 |
  |---|---|---|---|
  | Δx = 2 / 1 / 4 mm | −251.16752 / −251.16755 / −251.16744 | 0.9168 | **Δx 无关**（±1e−4 Pa） |
  | `inlet_ramp` 100 → 10 | −251.16758 | 0.9168 | **与迭代历史无关** |
  | mdot/2（u/2） | −62.79188 = −251.16752/**4.0000** | 0.9168 | **严格 ∝ ρu²** |
  | `inlet_ramp` 100 → 1000 | NaN | — | 稳定性边界极窄 |
  | 出口 `pressure-outlet ±10` / `100` | NaN / 超时 | — | **只有精确 0.0 收敛** |
  | 入口 `velocity-inlet` / far-field | NaN | — | ramp 只作用于 `mass-flow-inlet` |
  | 出口 `outflow`（README 组 3） | −248.1（+3 Pa） | — | 与出口 BC 类型弱相关 |

  系数 `平台/ρu_in² = −0.176943`（baseline 与 u/2 **逐位相同**）。
- **床组整场平移**（`shift.py`，逐面对差）：流体半段 +50.76、床半段直到 x=195 mm
  都等于"物理解 − 251.4 Pa"（**std 0.71 Pa**）；入口邻格 +301.88（物理解 302）；
  出口末列 −155.9 ⇒ **梯度完全不受影响**（−3021.3 vs 解析 −3021.7）。
- **代码路径（定位到行）**：`mod_uns_simple.f90:1035-1036` 入口单元动量行显式加
  **全量** ρu_in²A 流入动量；`mod_uns_bc.f90:404-405` 入口面压力**零梯度**
  （`pf = pP`）；`flux_rhiechow`（`mod_uns_simple.f90:1239-1244`）对**边界面直接
  `flux = ρ·u_f·S_f`、跳过 Rhie-Chow 压力梯度修正** ⇒ 压力 Dirichlet 补丁无法把
  内部电平锚回物理值。三者叠加 = 入口把电平印成 −0.177ρu_in² 并被零模内部带到全场。
  首格速度亏损 −8.32% 与耦合算例"uns 首排 −9.67%"**同族（疑同源）**。
- **交付**：`cases/couple_porous/README.md` 新增 §5.1（机理定位 + 表格 + 代码路径）
  ＋ §5.2（复现命令）；新增 `uns_full/bias_scan.sh`（端到端校验通过：5 个收敛算例
  逐位复现 + 汇总表输出）；`memory-bank/progress.md` 阶段 11 该条改 `[~]` 已定位，
  并**新登记修复项**（验收 = baseline 平台 ≈ 0、床组保梯度、出口 ±10 能收敛、
  耦合首排复查，回归全 RC=0）；`activeContext.md` 下一步同步。
- **对既有判据影响**：**无**（C2/流 B 只比梯度与界面连续性）。但留下两条限制：
  ① 修复前禁止跨求解器比 uns **绝对**压力；② 该偏置伴随**鲁棒性隐患**——出口
  压力值一改就发散（原先未意识到），这本身值得列入修复验收。

---

## 2026-10-06 | 瘦身后校验：`.git` 120 MB → 46 MB，远端 `main` 已 force-with-lease 覆盖

- `git filter-branch --index-filter`（`--prune-empty -- --all`）**5 个提交全部重写**：
  新历史 `63af789`（init）→ `54b2571` → `7802068` → `25d86ac` → `f2ebf0b`（瘦身节点）。
- 清理链：`refs/original/*` 全删 → `git reflog expire --expire=now --all` →
  `git gc --prune=now`。
- **实测**：`.git` **120 MB → 46 MB**；跟踪文件 **324**；
  `git log --all --name-only --pretty=format: | grep -c '\.vtu$'` = **0**；
  `^cases/.*\.log$` = **0**；仍跟踪的 `.vtu/.log/.dat` 只剩
  `regress/m6wing/baseline/*`（7 个位级 oracle）+ `external/**/M6-wing/Mesh3d.dat`
  （上游原始输入）。工作区 `git status --porcelain` = **0 项**（761.8 MB 产物在磁盘、
  被 `.gitignore` 挡住）；`lib/*/*.a`、`src/**`、`cases/**` 输入件、`scratch/unMesh.cas.zsplit_bug`
  全部仍在跟踪（逐项 `git ls-files` 复查通过）。
- **推送**：`git fetch origin`（先把**被 filter-branch 一并改写**的
  `refs/remotes/origin/main` 复位为远端真值 `3fa7b91`）→
  `git push --force-with-lease=main:3fa7b917db557d38d451d543b2424eb87a0c7f73 origin main`
  → `+ 3fa7b91...f2ebf0b main -> main (forced update)`，**零 remote 告警**；
  之后 `main` 与 `origin/main` ahead/behind = **0/0**。
- **⚠️ 安全网更正（重要）**：原计划的本地分支 `backup/pre-slim-2026-10-06` 被
  `--all` **一并重写**（只剩同名骨架），已 `git branch -D` 删除。**唯一完整旧历史在
  `/home/sundong/mixsolver_pre_slim_backup/repo_pre_slim.bundle`（98 MB，
  `git bundle verify` = "records a complete history"）**；另有
  `cases_hardlinks/`、`grid_BC_hardlinks/`（硬链接快照）、`checkpoint_2026-10-06_C2.tar.gz`
  与 `LNTE1D.vtu` 副本。恢复：`git fetch /home/sundong/mixsolver_pre_slim_backup/repo_pre_slim.bundle 'refs/heads/main:refs/heads/restored-pre-slim'`。
- **两条经验（已写进规则）**：① `filter-branch --all` 会连**备份分支**与**远端跟踪引用**
  一起重写 ⇒ 备份必须用**仓库外的 bundle**，不能只靠分支；② 重写后推送前必须
  `git fetch` 复位 `origin/main`，否则 `--force-with-lease` 会拿被改写的本地缓存去比对
  而误报（或更糟：用 `--force` 绕过保护）。
- **操作坑（已实测）**：上面那条 `git fetch origin` 会把**旧的远端对象**一并拉回本地
  （`.git` 46 MB → 140 MB；`git count-objects -vH` = 2 packs / 138.53 MiB）。因为
  `--force-with-lease` 推送后**没有任何引用**指向旧历史，所以再跑一次
  `git reflog expire --expire=now --all && git gc --prune=now` 即回收（实测回到
  **1 pack / 434 对象 / 44.98 MiB = 46 MB**，`git fsck` 干净）。⇒ 规则补一条：
  **重写 + 推送后要再 gc 一次**。

---

## 2026-10-06 | 仓库瘦身：212 个运行产物脱离跟踪 + 历史重写（剥离 761.8 MB）

- **触发**：本次会话统计出跟踪总量 916 MB / `.git` 120 MB，其中 `cases/` + `grid_BC/`
  的运行产物占 **761.8 MB**，`cases/channel/channel_*.vtu` 单个 50.2 MB。
- **逐个判定（不按扩展名一刀切）**：
  - **清**：`cases/**`、`grid_BC/**` 的 `*.vtu/*.log/*.dat/*.out/*.tmp/*.part.map/
    __pycache__`、根目录 `LNTE1D.vtu`、`cavity.vtu`、`build_*.log`、`m6_*.log`、
    `output_para.out`、`checkpoint_2026-10-06_C2.tar.gz` → 共 **212 个**。
  - **保**：输入件 `*.cas/*.neu/*.cgns/*.x/*.control/mix.control/bc3d.*`、工具 `*.py`、
    `README.md`、`images/*.png`、`scratch/unMesh.cas.zsplit_bug`（plan 与 activeContext
    明确"留档"的失效网格）、`regress/m6wing/baseline/*`（位级 oracle）、`external/**`。
  - **依据**：源码 `grep` 证明全是运行期写出——`Step_mess.dat`
    (`mod_struct_io.f90:348`)、`part_grid.dat`/`partation-auto.dat`
    (`mod_struct_mpi.f90:160/169`)、`mesh-quality.dat` (`mod_struct_grid.f90:732`)、
    `*.part.map` (`main_uns_mpi.f90:83`)、`flow3d*.dat` (`mod_struct_io.f90:408`)、
    `output_para.out` (`mod_struct_init.f90:705`)；各 `cases/*/README.md` 均把
    `.vtu/.log` 记为"产出"并附复现命令；并全文检索确认**无**脚本/README 把 `.vtu`
    当输入或 oracle（唯一的"读"出现在求解器自报输出名与 `cases/ltne/README.md`
    的"先跑再画图"流程里）。
- **执行**：`git rm -r --cached --pathspec-from-file=<212 条 NUL 清单>`（35.1 M 行删除）
  → `.gitignore` 增"运行产物不入库"规则（`*.vtu/*.log/*.part.map/*.tmp` 全局，
  `*.dat` **只**限 `cases/**` 与 `grid_BC/**` 以保住 `regress/m6wing/baseline/*.dat`）
  → 提交节点 → `git filter-branch --index-filter` 重写全部历史
  （本机**无** `git-filter-repo`，故用 filter-branch）→ `rm -rf .git/refs/original`
  + `git reflog expire --expire=now --all` + `git gc --prune=now`
  → `git push --force-with-lease origin main`（**结果校验见下一条记录**）。
- **安全网**：`/home/sundong/mixsolver_pre_slim_backup/`：`repo_pre_slim.bundle`（98 MB，
  `git bundle verify` = "records a complete history"）、`cases_hardlinks/`、
  `grid_BC_hardlinks/`（硬链接快照，秒级、零额外空间）、`checkpoint_2026-10-06_C2.tar.gz`
  与 `LNTE1D.vtu` 副本。
  （原拟用作安全网的本地分支 `backup/pre-slim-2026-10-06` 在 `--all` 重写时被一并改写，
  已删除，详见上一条校验记录。）
- **工作树策略**：761.8 MB 产物**不删除**，留在原路径转为未跟踪（被 `.gitignore`
  挡住），需要时按各 README"复现"小节重算。
- **同批勾销**（清单与代码漂移的 9 条 + 1 条陈旧残留，逐条给了证据）：
  阶段 1 收尾 `mod_unit_convert/mod_reference_state`；阶段 4/6 `MPI 子通信域+跨域收发`、
  `界面插值接入驱动`、`rho 占位`、`n_couple/n_uns_steps 读取`、`mm→m 缩放`、
  `grid_BC 真实物理验证`；阶段 7 三单元测试目标；阶段 10（已完成区早已 `[x]`，
  待办区是陈旧残留）。
- **保留待办（未勾销，避免假完成）**：界面收敛检查（判据未定，已重新界定现状）、
  阶段 8 `<1%` 判据（实测 −2.9%，与"壁面层对齐"联动）。
- **新登记**：`cases/couple_channel/README.md` 缺失（流 B 复现命令与定值只散在
  `docs/plan.md` 与 worklog）。

---

## 2026-10-06 | 同步到 GitHub：origin = git@github.com:cfddog/MixSolver.git（首次推送成功）

- **预检**：`~/.ssh/id_ed25519` 已配 → `ssh -T git@github.com` 返回
  `Hi cfddog! You've successfully authenticated`；HTTPS 200；远端为空仓库
  （`git ls-remote` 无引用，故 `push` 为 fast-forward 建分支，无覆盖风险）。
- **`.gitignore` 无需改动**：最大单文件 8 个 `cases/channel/channel_*.vtu`
  = 50,206,535 B（50.2 MB = 47.9 MiB）< GitHub 的 50 MiB 告警阈值与 100 MB 硬限；
  跟踪总量 916 MB / `.git` 120 MB。
- **推送**：`git remote add origin …` ＋ `git push -u --progress origin main` →
  `* [new branch] main -> main`，`branch 'main' set up to track 'origin/main'`。
  统计：**594 对象 / 93.46 MiB pack / ≈2.9 MiB/s**；`remote: warning` **0 条**；
  推送后 `HEAD == origin/main == 8f40a31`，`git status -sb` = `## main...origin/main`
  （无 ahead/behind）。
- **规则扩展** `.trae/rules/project_rules.md`「版本控制约定」新增两条：
  ① **远端同步**：提交后自动 `git push origin main`（同一节点内、不单独询问）；
     推送前查 `git status -sb` 无分叉（有分叉先 `pull --rebase`）；已共享历史禁止
     `--force`；网络不可达时跳过并在 worklog 注明"本节点未推送"。
  ② **仓库体量提示**：记录 916 MB / 120 MB / 47.9 MiB 最大单文件与"≥100 MB 会被拒"。
     同时修正回退条目（已推送历史不再用 `--amend` 重写）。
- 备注：GitHub 仓库名 `MixSolver`、本地目录 `MixNSSolver`，路径不同不影响同步。

---


- **用户要求**：每次完成一个小节点（一次可验证的改动 / 缺陷修复 / 文档收尾）都
  **自动 commit**，并把该流程计入规则。
- **规则改动** `.trae/rules/project_rules.md`：
  ① 「每次任务结束或关键决策后，必须更新」新增一项 —— git 提交（自动、无需询问），
     与 memory-bank 三项同批提交；
  ② 「版本控制约定」首条改写为**自动提交（默认行为，不必再问我）**的完整流程：
     「验证 → 更新 memory-bank → `git add -A && git commit -F -`」，**一个节点 = 一个
     提交**（不混装、不攒批）；验证未通过/半成品**不提交**（保留脏工作区或 `git stash`，
     并在 worklog 写明停在哪一步）；动过解算器/耦合/网格时加跑位级回归。
- **本节点即该规则的首次执行**：规则文本 + memory-bank 更新（`activeContext` 的
  CHECKPOINT 工作流说明、`progress.md` 版本控制条目补充）**同批提交**，无需另行确认。
- 备注：提交信息统一用中文 `git commit -F -`（heredoc），避免 `-m` 多行转义问题。

---


- **`git init -b main`**：仓库原本**无** git（也无 `.git` 目录）；按 `progress.md`
  「环境准备」里的原待办补齐版本控制。全局 `user.name/email` 为空，故只设**仓库本地**
  身份 `sundong <sundong@localhost>`（如需改：`git config user.email …` ＋
  `git commit --amend --reset-author`）。
- **`.gitignore`**（`build/`、`bin/` 为原待办要求，另补通用垃圾）：`build/`、`bin/`、
  `*.o`、`*.mod`、`__pycache__/`、`*.pyc`、`*.swp`、`*~`、`.DS_Store`。
  **刻意不忽略** `lib/{metis,parmetis,tecplot}/*.a`（vendored 预编译第三方库，链接
  必需、不可由本仓库重建）与 `regress/m6wing/baseline/*`（位级回归 oracle）——
  该判断已用注释固化在 `.gitignore` 内。核实：源码树无散落 `*.o`/`*.mod`。
- **首个提交 `97f247f`**「chore: 初始化版本库（首个提交 = C2 完成节点快照）」：
  536 文件 = 排除 `build/`、`bin/` 与 3 个 `__pycache__/*.pyc` 后的**全部**文件
  （逐个核对无遗漏）；`.git` ≈119 MB（ASCII `*.vtu` 压缩 ~7×）；`git status` 干净。
  说明：`cases/` 下 757 MB 生成产物（`*.vtu`、日志、`Residual.dat` 等）按原待办
  **未**忽略，一并入库；日后若要瘦身：`git rm --cached <path>` ＋ `git commit --amend`
  （本地无远端，重写历史安全）。
- **回退方式**：单文件 `git checkout 97f247f -- <path>`；未提交改动 `git stash`。
- **遗留**：`checkpoint_2026-10-06_C2.tar.gz`（仓库根，上一轮打的文档快照）已随首个
  提交入库，有 git 后冗余 —— 可 `git rm checkpoint_2026-10-06_C2.tar.gz` 后删除磁盘
  文件（留着也不影响）。
- **新增规则**：`.trae/rules/project_rules.md` 增加「版本控制约定」小节（提交粒度、
  `.gitignore` 白/黑名单、提交信息格式、回退方式），避免后续会话误忽略 `lib/*/*.a`
  或 `regress/` 基线。

---


- **登记待办**：把 C2 暴露的唯一遗留问题 —— **uns 单求解器绝对压力水平
  ≈−250 Pa 内部偏置的机理定位与修复** —— 计入后续待办：
  `docs/plan.md` 阶段 11 **第 3 条**（原「前置依赖」顺延为第 4 条）、
  `memory-bank/progress.md`「待办 · 阶段 11」、`memory-bank/activeContext.md`
  的 CHECKPOINT 与 C2 节「下一步」。条目内含：现象特征（内部整段常数下移、
  入口邻格保物理解水平、出口末列回收、严格 ∝u²；与多孔/出口 BC 无关）、
  候选根因（mass-flow-inlet 与投影步相容性源项；压力锚点说不足）、影响面
  （所有 uns 单求解器算例的绝对压力）、入手处（`src/unstructured/` 投影步与
  BC 层）、取证指针（`cases/couple_porous/README.md` §5＋§6.1、
  `uns_full/plane_profile.py`）。
- **节点保存**：本仓库无 git，改动全部落盘 ——
  `cases/couple_porous/{README.md, compare_iface.py, gen_meshes.py,
  uns_full/{gen_full.py, check_bed_gradient.py, plane_profile.py}}`、
  `cases/couple_channel/{gen_meshes.py, uns_full/gen_full.py}`、`docs/plan.md`、
  `memory-bank/{activeContext,progress,worklog}.md`、`.trae/rules/project_rules.md`；
  另打快照 `checkpoint_2026-10-06_C2.tar.gz`（仓库根，含上述文档与脚本，约 0.1 MB）。
- **复验（保存后实跑）**：5 个 Python 文件 AST OK；`check_bed_gradient.py`
  PASS（302.13 vs 302.1 Pa）；`compare_iface.py` 界面跳变 −0.08 Pa；
  `plane_profile.py` 参考 −3021.3（0.105–0.195 窗）/ −3021.5（≥0.110 窗）、
  耦合 −2959.3、全流体平台 −251.17 Pa。

---

## 2026-10-06 | C2 可压缩（struct）–多孔（uns）跨组界面：完成验证（两处缺陷修复）

- **结论**：求解器本体零改动；`cases/couple_porous` 的两处缺陷都在网格生成/
  后处理，修复后 C2 验收通过（界面压力跳变 −0.08 Pa；床梯度 0.01%）。
- **缺陷①（网格）**：`uns_full/gen_full.py` 的 `cid()` 把 x 串联写成了 z 并联
  （`1+i+j*NX+k*NX*NY`，k 最慢 ⇒ zone 2/8 各为一条 z=2.5/7.5 mm 的 x 全长薄片，
  z 向**并联**）；症状 = 床梯度仅为解析值 ~60%（1.78 vs 3.02 Pa/mm）、
  max|u| 42.87 > 入口 34.72、max p 398 Pa、zone 互换逐位相同（z 镜像）。
  改为 `1+k+j*NZ+i*NY*NZ`（x 最慢）⇒ zone2=1..2000(x<100 mm)、
  zone8=2001..4000(x>100 mm) 各为连续区间，床梯度 −3021.3 Pa/m（解析 −3021.7）。
  旧网格留档 `scratch/unMesh.cas.zsplit_bug`。
- **缺陷②（后处理）**：`compare_iface.py` 取窗"中心 ±1.25 mm"（宽 2.5 mm）> dx=2 mm，
  把 x=101 mm（295 Pa）与 x=103 mm（70 Pa）平均成 182.6 Pa ⇒ 幻影 112 Pa 界面跳变；
  改为绝对窗 `[100,102) mm` 后 **−0.08 Pa**（struct 面 295.0 == uns 首排 295.0）。
- **新发现（写入 README §5 / plan.md）**：uns 单求解器**绝对压力水平**含
  ≈−250 Pa 内部偏置。四组同网格对照：全流体 + pressure-outlet 内部整段严格平
  −251.17 Pa；全流体 + `outflow` −248.1（仅差 3 Pa，且末列速度 34.72→36.79）；
  床 + pressure-outlet 满足"物理解 −251"；床 + outflow 最低点 −546.7
  （= 平台 −244 − 床压降 302）。⇒ 与多孔无关、与出口 BC 类型无关；凹陷深度
  严格 ∝u²（−285.76 / u2 −71.44 / 2u −1143.06 Pa）；结构 = 内部常数下移 ＋
  入口邻格保持物理解水平（+301.88≈302、−0.05≈0）＋ 出口末列回收。故梯度不受
  影响、绝对水平受影响；C2 判据取"界面连续 + 床梯度"，与参考绝对压力脱钩。
- **耦合侧自洽（新证据）**：耦合 uns 床梯度 −2959.2 Pa/m 与
  μu/(εK)+ρβu² 代界面速度 34.233 得 2959.2 Pa/m 逐位吻合 ⇒ −2.1% 残差全部
  来自 struct 侧界面速度的 −1.37%（既有壁面层残余）；uns 首排 −9.67% 是入口
  邻格固有特征（三算例比值 0.9168/0.9158/0.9168），非耦合缺陷。
- **新增工具**：`cases/couple_porous/uns_full/plane_profile.py`（逐 x 平面 p/u
  均值 + 拟合 dP/dx，兼容 12 子四面体/5 节点两种 VTU 布局）——§5 四组对照取证。
- **复现链闭环（本次新跑一遍验证 README §6）**：`mpirun -np 2 bin/mixsolver_mpi
  mix.control Mesh3d.x control.ec unMesh.cas unMesh.control`（`n_couple=400`、
  `save_interval=100`，约 4 min）→ 产出 `unMesh_coupled.vtu`（uns 侧终场快照，
  `src/main.f90:475`）与 `flow3d.dat`（struct 侧周期存档；**save_interval 必须
  ≤ n_couple**，仓库默认 1000/1000 跑满才落盘一次）；随后 `compare_iface.py`
  逐位复现 README §4：struct +1.38% / RMS 1.37% / mean-ux −1.37%、uns
  +9.68% / 9.67% / −9.67%、p 294.9 / 295.0 / 48.5 Pa、界面跳变 −0.08 Pa；
  `plane_profile.py` 给耦合床梯度 −2959.3 Pa/m。
- **规则固化**：`.trae/rules/project_rules.md` 新增「验证/诊断硬规则（2026-10-06，
  源自 C2）」5 条 —— ①禁止跨求解器比绝对压力（uns ≈−250 Pa 偏置）；②取窗/拟合
  窗宽 ≤ 一个网格单元；③多 cell zone 网格 cell-id 排序须让分割方向索引最慢；
  ④复现耦合算例须 `save_interval ≤ n_couple`；⑤验证脚本须自带回归守卫。
- **文件**：`cases/couple_porous/{README.md(新建), compare_iface.py, gen_meshes.py,
  uns_full/{gen_full.py, check_bed_gradient.py, plane_profile.py}}`、
  `cases/couple_channel/{gen_meshes.py, uns_full/gen_full.py}`（同款 cid 注释）、
  `docs/plan.md`、memory-bank/activeContext.md、memory-bank/worklog.md（本文件）。

---

## 2026-10-06 | 关闭自检缺口③（uns restart 停用）+ ④（M6 基线重建并持久化）

- **缺口③ 修**：`src/unstructured/mod_uns_driver.f90` 的 restart 读取此前被
  "TEMP WORKAROUND" 整段停用（`use mod_uns_restart` 被注释，串行/MPI 两树同时失效）。
  根因：`mod_uns_restart` 属 MPI-only 栈（-> mod_uns_mpi_core/partition/local_mesh），
  仅编入 MPI 树；把 `use` 一注释虽让串行树可编，却连 MPI 树的 restart 一起废掉。
  改法：改用 `#ifdef HAVE_MPI` 守护（与 main.f90 同款）——MPI 树恢复
  `read_field_dump(..., serial=.true.)`，串行树走 `#else` 打印不可用告警；
  `ier2`/`src_file` 声明并入 `#ifdef`，串行树零 unused 警告。
  行为验证（couple_channel 拷贝，save_interval=1，n_couple=3）：
  - 第 1 轮存盘 → `unMesh_restart.dat`(650 KB)/`flow3d.dat`/`couple_state.dat`；
  - 第 2 轮 `couple_restart=1`（`bin/mixsolver_mpi`）→ `Field dump read from:` +
    `[uns] restarted from dump:` + `restart from coupling iter 3`，RC=0；
  - 同场景用串行 `bin/mixsolver` RC=0，打印
    `[uns] WARNING: restart unavailable in the serial (non-HAVE_MPI) build`，
    struct 侧照常重启，无崩溃。
- **缺口④ 修**：`/tmp/m6reg` 基线已丢（/tmp 易失），M6 位级回归不可复现。
  按 constraints.md L91 / worklog L1108 的既定配方重建并**持久化进仓库**：
  - 参考源 `external/OpenCFD-EC-1.16a` 拷至 /tmp（不在 external/ 内编译），逐文件
    `iconv GBK→UTF-8`、全角"！"→半角；`mpif90 -O2 -std=legacy
    -ffree-line-length-none` 编译出 `opencfd-ec1.16a.out`（0 error）。
  - M6-wing 4 块算例 np1、t_end=0.501、Kstep_save=50（第 50 步存盘，共 51 步）。
  - 结果：`flow3d.dat` = 13,600,032 B，md5 **dc134a2d196422043ecad7c86ac8f898**
    ——与历史记录/旧基线**逐位一致**；SA3d/wall_dist/partation-auto/part_grid
    二进制 cmp 一致，Step_mess/bc3d.inc/mesh-quality 文本一致；唯一差异
    `output_para.out`（阶段11 删 Turbo/Ref_medium 键的 5 行良性差异）。
  - 新增持久资产 `regress/m6wing/`：`control.ec`、`baseline/`（8 个基线产物 +
    `md5sums.txt`）、`run_regression.sh`（建树→跑→cmp；`M6_REBUILD_REF=1` 可从源
    重建参考）、`README.md`（配方/结果/用法）。
  - 自检：`regress/m6wing/run_regression.sh` 端到端 **PASS**（8/8 IDENTICAL，md5 命中）。
- 回归：`make -j4 all mpi structured(_mpi) unstructured(_mpi) units_test match_test
  coupling_test` RC=0；units_test 3/0 PASS；coupling_test np2 6/0 PASS。
- 涉及文件：`src/unstructured/mod_uns_driver.f90`；新增 `regress/m6wing/{control.ec,
  run_regression.sh,README.md,baseline/*}`。
- 至此自检缺口 **①②③④ 全部关闭**。

---

## 2026-10-06 | 修复串行 `units_test` 构建（自检缺口②）

- 问题：`make units_test` 链接阶段 RC=2：`build/ser/coupling/mod_coupling_exchange.o`
  经 mpif90 编译后含 MPI 符号，但 `bin/units_test` 仍用 `$(FC)`=gfortran 链接
  → `undefined reference to 'mpi_recv_/mpi_send_'`。
- 根因：`UNITS_TEST_OBJ_S` 取 `$(filter $(S)/coupling/%.o,...)`，把整个耦合层
  （含 MPI-intrinsic 的 mod_coupling_exchange）都纳入；测试自身不用 MPI，但链接
  命令漏了 MPI 运行时。另 `test_units_exchange.o` 在 pristine 树无 .d 时可能先于
  耦合模块编译（顺序缺口）。
- 改动（`Makefile`，无源码改动）：
  - `bin/units_test` 链接 `$(FC)` → `$(MPIFC)`（与 `bin/mixsolver` 完全一致，
    ser 树 = 可单进程 MPI 构建，`-np 1` 运行）。
  - 新增 pristine 边 `$(S)/coupling/test_units_exchange.o -> 耦合层`。
  - 更新 L161 注释说明该测试为何用 MPIFC 链接。
- 验证：`make -j4 units_test` RC=0；模拟竞态（删耦合 .o/.d/.mod + 二进制后 `-j4`）
  RC=0；运行 `./bin/units_test` 与 `mpirun -np1` 均 3/0 PASS。
- 涉及文件：`Makefile`。
- 仍未处理：缺口③ mod_uns_driver restart 停用、④ /tmp/m6reg 基线丢失。

---


## 2026-10-06 | 修复串行 `make all`：耦合层 + 主程序按 MPI-intrinsic 构建（自检缺口①）

- 问题：`build/ser/coupling/mod_coupling_exchange.o` 与 `build/ser/main.o` 在
  gfortran 串行树报 "Cannot open module file 'mpi.mod'"（两文件无条件 `use mpi`）。
- 根因：`bin/mixsolver` 本就用 `$(MPIFC)` 链接，`src/structured/` 也早有
  "ser 树同样走 mpif90" 的规则；耦合层与 main.f90 同为 MPI-intrinsic，却仍落入
  通用 gfortran 规则；且 main.f90 是串行池里唯一引用 MPI-only 重启栈
  (`mod_uns_restart`) 的文件。
- 改动（2 个文件，零算法改动）：
  - `Makefile`：新增 `$(S)/coupling/%.o` 与 `$(S)/main.o` 用 `$(MPIFC)` 编译
    （先例 = structured）；PRISTINE-BUILD EDGES 增 `$(S)/main.o`、`$(M)/main.o
    -> 耦合层`（耦合层已传递依赖 common+structured+unstructured）。
  - `src/main.f90`：`use mod_uns_restart` 与 `write_field_dump` 调用（uns 自动
    保存）用 `#ifdef HAVE_MPI` 守护——MPI 树照常保存，串行"单进程 MPI"树省去
    MPI-only 重启栈。
- 验证：复现场景 `make -j4 all`（仅耦合/main 缺失）RC=0 ✓；干净树 `make -j1 all`
  RC=0（干净 `-j4` 仍受既有 common 引导竞态限制，header 已注明先 `-j1`）；回归
  `make -j4 mpi/structured(_mpi)/unstructured(_mpi)/coupling_test/match_test` 全
  RC=0，`coupling_test` np2 = 6/0 PASS；ser `bin/mixsolver` np2 跑通 couple_channel
  （3 轮，写 unMesh_coupled.vtu）；`save_interval=1` 对照——MPI 树写
  `unMesh_restart.dat`(650 KB)，ser 树跳过 uns dump、struct flow3d.dat 照写、无崩溃。
- 涉及文件：`Makefile`、`src/main.f90`。
- 仍未处理：缺口②串行 units_test、③mod_uns_driver restart 停用、④/tmp/m6reg 基线丢失。

---


## 2026-10-06 | 构建/回归自检：MPI 全绿，串行混合/units_test 既有缺口确认（无源码改动）

- 目的：不改源码，确认当前树可编译可运行。
- 环境：gfortran/mpif90 13.3.0，OpenMPI 4.1.6，GNU Make 4.3。
- 结果（均 `make -j4`）：
  - ✅ `mpi` → `bin/mixsolver_mpi`；`structured`、`structured_mpi`、
    `unstructured`、`unstructured_mpi` 全部 RC=0；无 error（仅 legacy
    unused-variable 警告）。
  - ✅ `coupling_test`/`match_test` 构建 RC=0；`mpirun -np2 bin/coupling_test`
    = "6 tests, 0 failures / PASS"。
  - ✅ 冒烟：/tmp 拷贝 couple_channel（n_couple=3）`mpirun -np2 bin/mixsolver_mpi`
    RC=0，两 group 各就位、3 轮耦合完成、写 unMesh_coupled.vtu，无 NaN。
    （match_test 需在算例目录/带参运行，从根目录直跑缺 control.ec 属调用问题。）
- **确认的既有缺口（未修，超出本次范围）**：
  1. 串行 `make all` RC=2：`build/ser/coupling/mod_coupling_exchange.o` 与
     `build/ser/main.o` 均 "Cannot open module file 'mpi.mod'"——两文件
     无条件 `use mpi`，却被 Makefile 放进 gfortran 串行树。根因比 memory
     旧述（mod_uns_restart）更靠前。
  2. 串行 `units_test` 不可构建：UNITS_TEST_OBJ_S 含 mod_coupling_exchange
     （use mpi）；且 `-j4` 下 test_units_exchange.o 缺 .d 先于
     mod_interface_units.o 编译（pristine 顺序缺口）。
  3. `src/unstructured/mod_uns_driver.f90` 的 restart 读取被 "TEMP
     WORKAROUND" 注释停用（serial restart not supported）——与阶段10
     "uns 侧可重启"记录不一致。
  4. `/tmp/m6reg` 基线已丢失，M6 位级回归暂不可复现（仅存 md5
     dc134a2d196422043ecad7c86ac8f898）。
- 涉及文件：无源码改动；仅 /tmp 自检日志与 scratch 拷贝。
- 结论：MPI 构建路径完全可用；串行混合/测试目标为既有缺口，未触碰。

---


## 2026-10-06 | 状态检查点：阶段 11 流 B + C1 完成，余可压缩-多孔 + 分派表

- 用户要求保存当前项目节点，随时可回退到此状态。
- **已完成全部工作的快照**：
  - 阶段 1–6：项目骨架 + 结构/非结构求解器迁移模块化 + 交界面匹配 +
    数据交换/量纲统一 + 弱耦合驱动 + MPI 跨组通信（全部位级回归通过）
  - 阶段 9：多孔介质全套（Darcy-Forchheimer/LTE/LTNE/热弥散/各向异性/
    BJ + 7 验证算例 + α/K 标定设计），全部完成
  - 阶段 10：流场自动保存 + 耦合联合重启（byte-identical），完成
  - 阶段 11 流 B：couple_channel Dirichlet-Neumann 特征界面（iter1000
    质量精确连续 33.71==33.71 m/s，压力跳变 <1 Pa），完成
  - 阶段 11 C1：低速-多孔界面确认复用 BJ，勾销
  - 阶段 11 清理：弃用 IF_InnerFlow/IF_TurboMachinary 全套（M6 位级回归）
- **当前待办（下一步）**：
  1. 可压缩（struct）–多孔（uns）跨组界面
  2. 界面分派表（按 cell_zone + solver 自动分派）
  3. 远期：阶段 12 Gambit NEU、阶段 13 SST/Liao FD
- **关键资产**：
  - 流 B 参考解 cases/couple_channel/uns_full/unMesh.vtu
  - 流 B 工具 cases/couple_channel/compare_iface.py
  - M6 回归基线 md5 dc134a2d196422043ecad7c86ac8f898
  - 运行命令 mpirun -np N bin/mixsolver_mpi mix.control Mesh3d.x control.ec unMesh.cas unMesh.control
- 涉及文件：memory-bank/{activeContext,progress,worklog}.md。

---

## 2026-10-06 | 阶段 11 C1 勾销：低速-多孔界面确认复用阶段 9 BJ Robin

- 用户澄清三类界面的域归属：低速-多孔界面**只存在于非结构域内部**
  （低速流体 uns + 多孔 uns，cell zone 间内部面，不跨 solver、不经
  跨组交换层）；可压缩-低速 = struct↔uns（流 B，已完成）；可压缩-
  多孔 = struct 可压 ↔ uns 多孔（下一任务）。
- 初始 C1 构型设想（struct 低速流体 + uns 多孔）被用户纠正。
- **结论（用户确认）**：C1 无新代码、无新算例，勾销并直接复用
  阶段 9 已验证的内部 BJ Robin 通量（mod_uns_simple
  momentum_assembly：切向 C=μA(α/√K)/(1+α·d_Pf/√K)、法向两点
  扩散、交叉分量滞后等值反向）。验证资产 cases/beavers_joseph/：
  双层 4×4 ODE 对拍，α=0 RMS 0.14%、α=1/2 RMS 1.09%/1.29%、
  np2 vs 串行 2.2e-7。
- **涉及文件**：docs/plan.md（阶段 11 节：C1 勾销、域归属、C1 记录、
  优先级表）、memory-bank/{activeContext,progress,worklog}.md。
- 手册 docs/程序使用手册.tex 不存在 → N/A（无功能/参数变更）。
- **下一步**：可压缩（struct）–多孔（uns）跨组界面。

---

## 2026-10-06 | 阶段 11 流 B 完成：交换层 Dirichlet-Neumann 特征界面修法

- **目标**：couple_channel 界面连续性只能在耦合交换层内、用两侧状态
  解决（前序总压入口试验已证明与 struct 入口无关）；以 uns_full
  单一求解器解（x=100：mean u=34.72 m/s、p gauge=46.9 Pa）为对拍基准。
- **4 版迭代**（np2，n_couple=1000，n_struct/uns_steps=10，iface_ramp=20）：
  v1 删 uns 反射 + struct ghost=2face−inner（子步前设）→ 周期-2 翻转
  iter75 NaN；v2 ghost 移到子步后设 → 回声不动点、流量衰减至 0.04 m/s；
  v3 ghost 零阶直写 uns 首排格心 → 稳定但冻结 1.7 kPa 压力跳变、
  质量 31.2 m/s；v4 定稿 Dirichlet-Neumann：struct 亚声速出口只收
  uns 背压，ghost 用 boundary_Farfield 同款线性化 Riemann 反射
  （db=d1+(pb−p1)/c1²，ub=u1+(p1−pb)/(ρc)·n_out，p2=2pb−p1），
  法向按 Interface_List%face 1..6 从 Block ni/nj/nk 面法向表构造；
  uns 侧只施速度 Dirichlet，删 set_interface_p 与压力重锚定；
  alpha=0.3·min(1,iter/iface_ramp)（alpha=1 因 Ma=0.1 的 1/(ρc)
  高增益在 ramp 后周期发散，iter50 NaN/819 OverLimit）。
- **结果**（run_flowb4c.log iter1000，无 OverLimit，~100 iter 冻结）：
  质量精确连续（struct 面 33.71 == uns 下游 33.71 m/s）；压力连续
  （342.5 vs 342 Pa，<1 Pa）；对拍 uns_full 均值 −2.9%、界面压力
  +296 Pa、近壁型线最大偏差 11.8%（struct 近壁 u=13.1 vs uns 9.05
  @1.25mm），判定为两解器壁面层物理差异，非交换层缺陷。
- **源码改动**：src/main.f90（删 uns 反射/压力重锚定/set_interface_p；
  ghost 设置移至子步后；ALPHA_MAX=0.3）；
  src/structured/mod_struct_driver.f90（struct_set_iface_bc 改特征
  背压 BC，optional alpha，0=纯外推/1=全背压）；
  src/structured/mod_struct_bc.f90（异常：会话中发现文件被外部回退为
  1263 行含 Turbo 脏版，mtime 12:12:38，无 git 不可追溯；手工重建
  boundary_user 分派器 + boundary_user_Inlet——取自 external/OpenCFD-
  EC-1.16a 原始 sub_boundary_user.f90，最终 915 行 0 Turbo）；
  cases/couple_channel/compare_iface.py（修正 flow3d.dat 解析：
  output_flow 写 U(0:nx) 含 ghost，形状 5*(nk+1)*(nj+1)*(ni+1)）；
  cases/couple_channel/mix.control（补 save_interval/couple_restart）。
- **验证**：make mpi / make structured 通过；M6-wing 50 步回归
  flow3d.dat md5 dc134a2d196422043ecad7c86ac8f898（13,600,032 B）
  与清理版基线一致。
- **教训**：分段删除 Fortran 代码后必须重新 grep 行号（Python 按旧
  行号删会误删边界并拼入残段）；并行 make 增量编译偶发不重编依赖
  （"Keyword argument not in procedure"假象，串行重 make 即解）。
- **涉及文件**：见上 + docs/plan.md（阶段 11 流 B 节）。
- 手册 docs/程序使用手册.tex 不存在 → N/A。

---

## 2026-10-06 | 阶段 11 方向修正：彻底删除 IF_InnerFlow / IF_TurboMachinary 全套机制

- **背景**：couple_channel 总压入口试验（IF_InnerFlow=1, P_In_Ratio=1.0123）
  np2 iter50 速度 78 m/s 发散、界面偏差不变；结论：界面连续性与入口无关，
  用户拍板弃用 IF_InnerFlow（只留外流）与 IF_TurboMachinary 叶轮机模式。
- **源码清理（src/structured/，5 文件）**：
  - mod_struct_global.f90：删 IF_TurboMachinary/IF_InnerFlow/P_In_Ratio/
    Turbo_Periodic_seta/Turbo_w；
  - mod_struct_constants.f90：删 BC_Wall_Turbo=201；
  - mod_struct_bc.f90：分派重写为外流单一路径，删三个 Turbo BC 子程序；
  - mod_struct_solver.f90：删旋转惯性力源项块（离心+科氏）；
  - mod_struct_init.f90：删 Turbo_P0/T0/L0/w、Ref_medium_usrdef、namelist
    键/默认值/叶轮机 init else 分支/介质推导块/output_para 行/bcast 赋值
    （rpara 34/35/39、Ipara 27/29 槽位保留）；中途误删 AoS（=侧滑角，
    Fortran 大小写不敏感，与 Aos 同名）导致编译失败，已恢复；
  - mod_struct_mpi.f90：删 Umessage_Turbo_Periodic 及调用，两个
    *_Periodic 坐标子程序只保留平移分支（另修删除范围漏掉的孤立 else）。
- **算例**：couple_channel/control.ec 用 control.ec.forced 恢复强制均匀入口
  并删 5 个废键；grid_BC/control.ec 同删 5 键；删除 control.ec.forced、
  mix.control.ptot、run_ptot1.log。
- **验证**：make mpi / make structured 双树 0 error；M6-wing 50 步对拍
  external 原始基线（/tmp/m6reg，统一 -O2 -std=legacy，t_end=0.501），
  np1 与 np2 flow3d.dat（13,600,032 B）、Step_mess.dat 均 byte-identical；
  couple_channel n_couple=2 冒烟 rc=0（unMesh_coupled.vtu 被 2 步结果覆盖，
  流 B 重跑即可）。
- **涉及文件**：上述 5 个 struct 源文件 + 2 个 control.ec + docs/plan.md。
- 手册 docs/程序使用手册.tex 不存在 → N/A。

---

## 2026-10-06 | 阶段 10 完成：流场自动保存 + 耦合联合重启

- **mix.control 新键**：`save_interval`（>0 生效）、`couple_restart`（0/1）。
  解析于 mod_reference_state.f90，get_coupling_params 扩展 optional 存取。
- **自动保存（两侧同一 coupling iter 末对齐）**：
  struct `struct_solver_save` → `output_flow`（flow3d.dat+Step_mess.dat）
  + 新增 `output_flow_nodes`（格心 d/u/v/w/T 插值到 Mesh3d.x 节点，
  flow3d_node.dat，ascii/unformatted 随 Mesh_File_Format）；
  uns root（全局 rank = n_struct_ranks）`write_field_dump` →
  unMesh_restart.dat。struct rank0 写 couple_state.dat（iter 号；
  standalone nproc==1 时由 uns root 写）。
- **联合重启**（couple_restart=1）：读 couple_state.dat 得 iter0，
  循环 iter0+1 起（iface_ramp 计数连续）；struct `force_restart` 强制
  Iflag_init=1 读 flow3d.dat；uns `restart_file` 分支经
  `read_field_dump(serial=.true.)` 串行读 dump。
- **关键坑（修复）**：read_field_dump 尾部 MPI_Bcast 走 mpi_comm——耦合
  驱动不调 mpi_bootstrap（默认 COMM_WORLD），会把 struct rank 卷入挂死；
  新增 serial 选项跳过全部 Bcast，各 rank 独立读同一文件。
- **验证（grid_BC np2，n_couple=6，save_interval=3）**：连续 6 轮 vs
  3+重启 3 轮，flow3d.dat 与 unMesh_restart.dat 均 **byte-identical**
  （check_restart.py）；flow3d_node.dat NB/dims 与 Mesh3d.x 一致；
  ser/mpi 双树编译无新增 error。
- **涉及文件**：src/common/mod_reference_state.f90、
  src/unstructured/mod_uns_restart.f90（serial）、
  src/unstructured/mod_uns_driver.f90（restart_file）、
  src/structured/mod_struct_driver.f90（force_restart+save 包装）、
  src/structured/mod_struct_io.f90（output_flow_nodes）、
  src/main.f90（save/restart 分支+couple_state 助手）；
  归档 grid_BC/cont6/、grid_BC/restart33/、grid_BC/check_restart.py。
- **已知限制**：dump v1 不含 T_s（LTNE 重启 T_s 回 init）；struct
  Kstep/tt 不恢复（仅影响输出命名）。
- 手册 docs/程序使用手册.tex 不存在 → N/A。

---

## 2026-10-06 | 阶段9 全部收尾：α/K 标定实验设计（数值虚拟标定）

- 纯算例/脚本任务（无源码改动，无需回归）：cases/calibration/。
- **框架**（4 个 Python 文件）：
  - calib_problems.py：四子问题定义（参数/真值/界/观测布置/解析正演），
    正演模型复用阶段 9 已验证的半解析参考（plug 二次压降、mix erf、
    ltne 4 阶耦合 ODE、BJ 双层 4×4 ODE 封闭解）；
  - calib_invert.py：LM 反演（θ=log10 空间、盒界、自适应 λ）+ FIM 分析
    （归一化灵敏度、相关矩阵、CRB 精度界、条件数）；
  - run_calib.py：闭环驱动（真值+噪声→离真值 +0.3 dex 初值→反演→图）；
  - calib_design_scan.py：CRB vs 噪声选型表。
- **设计核心结论**：
  - 六参数按物理特征解耦为四子问题：P1 多流速压降→K+C_F（Δp 低速线性
    钉 K、高速二次钉 C_F，u∈[0.05,2]，Re_p 跨 0.03–12.8）；
    P2 混合层横向剖面→α_t（3 截面×15 点）；P3 LTNE 发汗冷却两相测温→
    α_l+h_sf 联合（α_l 签名=轴向回流/出口温升，h_sf 签名=T_s−T_f 间隙，
    空间形态正交）；P4 界面 u(y) 剖面→α_BJ。
  - **仪器选型级发现**：TC 噪声按绝对温度缩放（0.2%×305K≈0.6K），混合层
    α_t 标定若 ΔT=10 K 则 CRB 12.7% 不可用；加大到 ΔT=50 K 后 CRB 2.2%。
  - 闭环（seed=1，1%/0.2%/0.5%/2% 噪声）：恢复误差 K 0.19%、C_F 1.1%、
    α_t 3.8%、α_l 1.1%、h_sf 0.25%、α_BJ 1.8%，全部在 CRB 界内；
    多参问题相关系数 0.65–0.76，无病态。
- **衔接**：反演器与数据源解耦，forward 换求解器 VTU 采样即可对真实数据
  标定；四构型控制文件在 case_ref 目录。
- 归档：cases/calibration/{README.md, calib_problems.py, calib_invert.py,
  run_calib.py, calib_design_scan.py, images/calib_closure.png}；
  plan.md 阶段 9 第 2 条末款勾销，标"阶段 9 全部完成"。
- 手册 docs/程序使用手册.tex 不存在 → N/A。

---

## 2026-10-05 | 阶段9 收尾：LTNE+热弥散组合验证（1D 纵向 + 2D Graetz 横向）

- **源码改动（两个真 bug）**：
  - mod_uns_simple.f90：Bear 纵向弥散分量缺 /|u|（u_d² 量纲 m²/s²），
    改为 ρcp/umag·(disp_l·u_d²+disp_t·(umag²−u_d²))，注释同步修正。
  - mod_uns_gather.f90：gather_fields_to_root 补 temperature_solid 的
    sbuf/rbuf、Gatherv、unpack；main_uns_mpi.f90 rank0 fld_g 调
    setup_porous_fields（否则输出 guard 跳过固相段，MPI VTU 无 T_s）。
  - 双树 make 通过。
- **回归（disp 公式改写）**：porous_disp 两工况、porous_graetz 两工况
  新旧 VTU 逐位差 0；cavity/Ra10 不涉及弥散不受影响。
- **1D cases/ltne_disp/**：gen_ltne1d.py（100×2 hex SI 米，mass-flow-inlet
  type20，出口 zone7 tbc 恒热流），ltne_disp0/1.control（disp_l=0/0.005），
  plot_ltne_disp.py 4 阶耦合 ODE 半解析（三次特征根，exp(r(x−L)) 防刚性根
  溢出，Nield 微观分配）。8095/11219 步 dT_max=0；Tf 中位 0.4%、
  Ts max 0.83%/0.30%；disp1 ΔTf(L)=973.46 vs 974.41
  （aΔTf=q−Kf Tf'(0)，弥散回流到 Dirichlet 入口）；np2 11220 步，
  两相 rel 1e-7。
- **2D cases/ltne_graetz/**：复用 porous_graetz 网格（180×40），
  ε=0.4/k_s=1.0/H=2e5（Bi=0.845、Λs=0.101），ltne_grz_d/0.control
  （disp_t=2.5e-4/0）。plot_ltne_graetz.py 参考解：双相二次特征值问题→
  4m×4m 一阶块 y'=My，2m 个纯实衰减模态，入口 θf=1、χs=0 解系数；
  流体横/纵 Pe 分开（α_l=0 轴向仅 εkf）；FD η 向 N=180；Bi→0 极限
  λ1Pe=2.457≈π²/4。7018/6939 步；剖面 max 1.4%/3.3%（仅 x=20mm），
  gen_graetz_xref.py 360×40 x 加密（7081/7029 步）降至 0.73%/1.76%；
  bulk θf,b 沿 x*/ξ 全程重合；np2（7016 步）两相 abs 2e-6 K。
- **教训**：①hex 输出为每单元 12 tetra（记录 `4 n1..n4`），解析必须剥首
  token 且用 ix 圆整分层（plot_bj.py 已正确）；②LTNE Graetz 固相在入口
  不是 θs=1：固相无对流、入口绝热，θs(0)=耦合模态给出的横向平衡剖面；
  ③1D SIMPLE 双温度耦合松弛慢（~8-11k 步，inlet_ramp=1）；
  ④28800 hex 全向加密外层线性求解过慢（lin-it 361），定向 x 加密更高效。
- 归档：两算例 README.md + images/*.png + 控制/日志/生成脚本；
  docs/plan.md 阶段 9 勾销 LTNE+弥散项。

---

## 2026-10-05 | 阶段9 收尾：Beavers-Joseph 界面条件实现与验证

- **源码改动**：
  - mod_uns_control.f90：cell_zone_t 增 `bj_alpha`（默认 0=禁用）+ 解析；
    ctrl_t 增 `body_force(3)` + 主循环 `body_force = fx fy fz` 分支。
  - mod_uns_fields.f90：fields_t 增 bj_alpha(:)，setup_porous_fields 填充。
  - mod_uns_simple.f90：momentum_assembly 内部面循环插 BJ 块——流体/多孔
    界面（cztype 相异且多孔侧 bj_alpha>0）Robin 通量
    C=μA(α/√Knn)/(1+α·d_Pf/√Knn)，切向分量用 C、法向保 D、斜交交叉项
    滞后 rhs，两侧等值反向保动量守恒；BJ 面禁 nonorth 修正；Darcy 源项前
    加 body_force 源项循环（any(body_force/=0) 守卫）。
  - MPI 零改动（simple_mpi 导入 serial momentum_assembly）。
- **回归**：cavity 697 步/Nu=0.996041、porous_Ra10 748 步/Nu=1.078 位级不变。
- **算例 cases/beavers_joseph/**：gen_bj.py 双层开口槽道（100×40，流体上半
  zone2 / VC:porous 下半 zone3，两端 pressure-outlet，body_force=0.05 N/m³
  驱动，K=1e-6 → λ/dy=4 Brinkman 层可分辨）。三工况 α=0/1/2。
- **关键物理发现**：离散实现的连续极限=应力连续+速度跳变双层模型
  （串联阻力 1/C=d_Pf/(μA)+λ/(μαA)），非经典 BJ 滑移公式（后者假设纯
  Darcy 多孔侧）。plot_bj.py 用 4×4 双层耦合 ODE 数值解作精确参考。
  结果：α=0 RMS 0.14%、α=1 RMS 1.09%、α=2 RMS 1.29%；经典 BJ u_B 明显
  偏离（反证）。MPI np2 vs 串行 max rel diff 2.2e-7。
- **教训**：①闭盒体力被压力梯度抵消，验证须开口槽道；②被动温度漂移使
  dT_max 门槛不收，outer_max=800 截取已收敛速度场；③perm 一度误改 1e-8
  已回退 1e-6（VTU 与 control 一致）；④VTU 名随 mesh 名，多工况须 mv。
- 涉及文件：上述 3 源码 + cases/beavers_joseph/{gen_bj.py, bj.cas,
  bj_a0/a1/a2.control/.log/.vtu, bj_mpi.log/.vtu, plot_bj.py, README.md,
  images/beavers_joseph.png}；docs/plan.md 阶段9 B-J 项勾销。
- 手册 docs/程序使用手册.tex 不存在 → N/A。

---

## 2026-10-05 | 阶段9 第三批：二维多孔介质强制对流 Graetz 验证算例

- 纯算例任务（无源码改动）：cases/porous_graetz/。
- 构型：半平行板通道 y∈[0,a]，a=5mm、L=0.45m、180×40×1 wedge；
  y=0 中心线 symmetry、y=a wall 固定 Tw=310（tbc=8 1 310.0）、
  x=0 velocity-inlet T0=300（第 4 token，用上轮修复的解析）、x=L p-outlet。
  K=1e-10、eps=0.4、**k_s=k_f=0.026 使 LTE k_eff=k_f**（否则解析 κ 带 ε 因子）。
- 两组：mol U=0.1/disp_t=0（κ=2.198e-5，Pe=22.7）；
  disp U=1/disp_t=2.5e-4（κ=2.720e-4，Pe=18.4，出口 x*=4.9）。
  解析：θ=Σ2(−1)^{n+1}/λ_n cos(λ_nη)e^{−λ_n²x*}，λ_n=(n−1/2)π；
  体均 θ_b=Σ2/λ_n²e；渐近 Nu_a=λ_1²=π²/4=2.4674（特征长度半高 a）。
- 结果：5601/6238 步（outer_tol=1e-9）；mol 渐近 Nu_a=2.457（0.42%）、
  disp 2.434（1.38%）；剖面 max|dθ| ≤0.026（入口 x*=0.05）/≤0.013；
  两工况 θ_b(x*) 坍并到同一条级数曲线；np2 与串行同为 6238 步。
- **关键教训**：Bear 弥散张量基于格心速度，壁面无滑移 Brinkman 层（厚 ~√K）使
  近壁弥散自动衰减（物理真实的无弥散层），但 Graetz 解析假设弥散均匀到壁。
  K=1e-9（√K/dy=0.25，首排 u=0.895U）渐近 Nu 偏高 16%；K=1e-10
  （√K/dy=0.08，u1=0.988U）降到 1.38%。判据：√K/dy≤0.08。
- 归档：gen_graetz.py/graetz.cas/两 control+log+vtu/plot_graetz.py/
  images/graetz.png/README.md。plan.md 阶段9 验证列表勾选（二维多孔介质对流）。
- 分析教训：①级数体均系数是 2/λ_n²（不是 4，x=0 须归一化到 1）；
  ②Nu 壁面梯度用非均匀三点公式 (3θ1−θ2/3)/h；③x*≥1 后 θ_b<0.07
  相对误差被小分母放大，应看绝对偏差。

---

## 2026-10-05 | 阶段9 第二批：热弥散 + 各向异性渗透率 + VINLET 温度 bug 修复

- **各向异性渗透率（轴对齐对角张量，用户确认方案）**：
  - cell_zone_t 新增 perm_xx/perm_yy/perm_zz（0=回退标量 perm），parse_cell_zone
    增同名 key；fld 新增 perm_dir(3,ncells)，setup_porous_fields 按 merge 填充。
  - 动量 Darcy 汇改 `ap += (mu/porosity)/fld%perm_dir(comp,kk)*vol`
    （momentum_assembly 按 comp 外层循环，天然逐方向取对角元）；
    交叉项需块耦合动量装配，注释明确不支持。
  - 验证（cases/porous_plug 追加）：plug_aniso1 仅 perm_xx=1e-8 → dpdx 误差
    0.085%（横向无阻力时入口列恢复段变长 col0=0.0175，需 outer_tol=1e-9 +
    拟合窗口避开前 10 列）；plug_aniso2 perm_xx=5e-9 → dpdx 精确加倍（0.122%）；
    标量路径 0.000%。
- **热弥散（Bear 分量式张量）**：
  - cell_zone_t 新增 disp_l/disp_t（m），fld%disp_l/disp_t；
  - temperature_assembly 预计算 kdisp(3,ncells)=ρcp(α_L u_d²+α_T(|u|²−u_d²)/|u|)，
    内面/边界面 kf += n·D·n（ndsf=(sf/area)²·kdisp 加权平均）；LTE 加 k_eff、
    LTNE 加流体相 kf（固相不加）；非正交修正复用含弥散 kf（正交网格精确）。
  - 验证 cases/porous_disp/（新算例）：80×80×1 混合层，双温速度入口（zone6 310K/
    zone8 300K，VC:porous），α_t=1e-3 → κ=1.022e-3 m²/s，T 剖面 vs erf 相似解
    max 偏差 0.46–1.03%（首截面 1% 为台阶入口离散化）；disp_t=0 对照层宽 45 倍差；
    188 步收敛，np2=串行。**教训：erf 相似变量填 κ=Γy/(ρcp)，不是 Γy**
    （首轮对比脚本用错量纲差点误判模型失效）。
- **修复 velocity-inlet 温度 token bug**（历史遗留，独立运行温度排空到 0 的根因）：
  parse_bc_line 的 velocity-inlet 分支只读 ux uy uz，bc_face_T 的 BC_VINLET
  Dirichlet 锚定 tval=0。现可选解析第 4 个 token 为静温（缺省仍 0）。
  副作用：plug darcy/forch 刷新后 T 恒 300，收敛加快（61→12、57→33 步，
  outer_tol 收紧到 1e-9），压降不变（0.0000%）。混合层算例依赖此修复。
- Ra10 回归：748 步、Nu=1.07800、du_max=1.902e-9，与改前一致。
- 归档：cases/porous_disp/（README/gen_mix.py/两工况 control+log+vtu/
  plot_mix.py/images/mix_layer.png）；cases/porous_plug/README 增第 4/5 节
  （各向异性验证表 + VINLET 温度 bug 说明），plug 工况刷新。
- 涉及源码：mod_uns_control.f90（cell_zone 字段/解析/VINLET 温度）、
  mod_uns_fields.f90（perm_dir/disp_l/disp_t）、mod_uns_simple.f90
  （Darcy 张量化 + 能量弥散投影）。

---

## 2026-10-05 | 阶段9：VC 标记自动识别 + porous-plug 一维压降验证

- **VC 自动识别（用户需求 docs/程序功能说明.md L8）**：
  - resolve_cell_zones 类型解析改为三级优先级：显式 cell_zone 类型词 >
    CAS zone 名内 `VC:porous/fluid` 标记 > fluid 默认；报告新增 source 列
    （control / VC tag / default）；VC 与显式冲突时打印 note；
    VC 判多孔但无系数行时 WARNING（perm=0 会静默关 Darcy 汇）。
  - 新增 vc_tag_type/count_cells_type（mod_uns_control）：词边界匹配 'vc'，
    跳 ':'/空格读后续 token，尾标点容忍；常量 CZ_AUTO=0 公开。
  - parse_cell_zone：类型词变为可选——首余 token 含 '=' 即 CZ_AUTO
    （系数专用行），cell_zone 可只写 `2 perm=.. porosity=..`。
  - setup_porous_fields：cztype=POROUS 且 ctrl 行 ztype=POROUS/AUTO 即给系数；
    顺带修复旧笔误 ctrl%cz(j)%h_s → h_sf（此前被陈旧 .o 掩盖）。
  - zone_t cond_name/user_name 32→128 字符（VC 标记可能在长名称尾部）。
  - 注意：CAS 分词器不处理引号，(45 记录里 VC 标记必须无空格（写 VC:porous）。
- **验证算例 cases/porous_plug/**：gen_plug.py 40×2×1 直管（0.1m，全
  symmetry，x=0 vel-inlet/x=L p-outlet，cell zone (45 带 VC:porous）。
  Darcy 工况 u=0.1：内部压力斜率 −461.500 Pa/m vs 理论 −461.500（0.000%），
  外推 Δp=46.15Pa；Forch 工况 u=1、C_F=500：−5203.492 vs −5203.500（0.002%），
  Δp=520.349 vs 520.350Pa；内部列 u 偏差 ≤0.016%；61/57 步收敛，
  np2=串行（61 步、du_max 位级一致）。入口/出口列速度外观偏差仍是 SIMPLE
  边界列固有特性，不影响内部压降。回归 porous_Ra10：748 步、Nu=1.078 与
  归档一致，source=control。
- 归档：README.md（设置/解析解/复现命令/结果表）、gen_plug.py、两工况
  control+log+vtu、plot_plug.py、images/plug_pressure.png。
- 涉及源码：src/unstructured/mod_uns_control.f90、mod_uns_mesh.f90、
  mod_uns_fields.f90。plan.md 阶段9 对应条目勾选。

---

## 2026-10-05 | plan.md 录入后续功能储备（程序功能说明 L8–15 → 阶段 9–13）

- 纯文档任务：docs/plan.md 在阶段 8 后新增"后续计划"章并扩展优先级表。
  - 阶段9 多孔介质（VC 自动识别；热弥散/各向异性待办；验证 LTNE ✅、
    多孔自然对流 ✅、二维多孔对流/porous-plug/B-J 待办）
  - 阶段10 两侧流场（含 halo）自动保存 + 耦合联合重启；struct Plot3D 节点插值、
    格式随 Mesh3d.x
  - 阶段11 低速-多孔/低速可压缩/可压缩-多孔界面，按两侧块属性自动分派
  - 阶段12 Gambit .neu 网格输入；阶段13 SST 修复（暂缓）与 Liao 格心型 FD（远期）
- 已实现项（cell_zone fluid/porous、Darcy-Forchheimer、LTE/LTNE、uns restart）
  标 [部分实现]/[已实现]，避免重复列为待办。
- 同步更新 progress.md（阶段 9–13 待办小节）。无源码改动。

---

## 2026-10-05 | inlet_ramp 实现 + 耦合发散根因诊断（压力重锚定/质量守恒修正）

- **inlet_ramp 实现（代码就绪）**：`ctrl_t%inlet_ramp`（默认 1=禁用）；
  mod_uns_bc 增模块级 `g_inlet_ramp_factor` + public `set_inlet_ramp_factor`；
  bc_face_vel 的 MASSINLET 分支乘 ramp 因子；simple_run 与 simple_run_mpi
  每外迭代设 `min(1, it/inlet_ramp)`。涉及：mod_uns_control / mod_uns_bc /
  mod_uns_simple / mod_uns_simple_mpi。
- **诊断链（重要）**：
  1. standalone uns 在 grid_BC 发散属**预期**——interface 面 standalone 下
     iface_vel=0 等效壁面 → 域封闭（有入流无出流）→ 无稳态解，ramp 与否都发散。
  2. 耦合发散首因 = **struct 绝对压力 ~1e5 Pa 砸到 uns 不可压 p~0 参考系**
     （omega=0.03 也有 3000 Pa 界面跳变，it=1 du_max=3e4）。
  3. 次因 = **全 Dirichlet 速度边界下净通量必须为零**，否则纯 Neumann PPE
     无解持续漂移（struct 施加的界面通量一般不与 uns 入流平衡）。
- **main.f90 uns_group_driver 两项修正**：
  ① 压力重锚定 `sp_bc -= mean(sp_bc) - mean(ip)`（保留界面压力梯度、锚定
     到 uns 自身水平）；② 界面速度均匀法向修正：先对非界面边界用
     bc_face_vel 求 F_other（含 ramp 因子），再加 δu_n = −(F_other+F_iface)/A_iface
     使总边界通量为零。
- **结果**：it=1 mass-imbal 3.9e-10（机器精度，修正有效），du_max=151 有界，
  但 ~it=50 仍 ICC0 非正主元 → NaN。判定为 **grid_BC 物理不兼容**：
  struct Ma=3 可压超音速 vs uns 不可压低速腔体，弱耦合+Dirichlet 界面
  无法收敛（omega=0.03 也压不住 10 m/s 界面速度与 1.7 m/s 内部的阶跃）。
- **inlet_ramp 与界面质量修正互斥**：驱动级修正按满流计算，simple_run 内
  逐步 ramp 会造成持续质量汇；耦合算例 inlet_ramp 应保持 1（unMesh.control
  已注释说明）。
- **下一步建议**：改双向物性兼容算例（如腔体-腔体耦合，两侧均不可压），
  或在界面上改用通量型/混合 BC。
- 涉及文件：src/main.f90、src/unstructured/mod_uns_{control,bc,simple,simple_mpi}.f90、
  grid_BC/unMesh.control（inlet_ramp=1 + 注释，outer_max 维持 10）。
- 日志：grid_BC/uns_ramp_np2.log、couple_fix_np2.log、couple_fix2_np2.log。

---

## 2026-10-05 | unMesh.cas mm→m 缩放 + control.ec 事故恢复

- **mm→m 缩放**：`ctrl_t` 新增 `mesh_scale`（默认 1.0）；`mod_uns_control` 新增
  `read_mesh_scale` 预扫描（read_control 在 compute_geometry 之后，缩放必须在
  read_cas 后立即做）。三个入口均接入：mod_uns_driver / main_uns / main_uns_mpi。
  `grid_BC/unMesh.control` 设 `mesh_scale = 1.0e-3`。验证：坐标缩放生效，
  界面面积 1.25e-2 m² 与 0.25m×0.05m 一致，耦合协议不受影响。
- **事故**：`struct_solver_init` 的 `ln -sf <ctlfile> control.ec` 在 ctlfile 与
  cwd 下 control.ec 同名时生成自指 symlink，原 `grid_BC/control.ec` 被毁。
  **修复**：`read_parameter`/`read_parameter_ec` 增加 optional 文件名参数，
  彻底移除 symlink hack。control.ec 已按 `output_para.out` 回显重建。
- **遗留**：unMesh SIMPLE 仍发散（du_max→NaN）。standalone 发散属预期
  （interface 面无出流→封闭域）；耦合态发散因 struct Ma=3 初值直接作为
  uns 界面 BC（~340 m/s 量级的速度砸到 1.7 m/s 内流）→ 需界面状态
  ramp/松弛，列入阶段 7 TODO。

---

## 2026-10-05 | 阶段6 验证：常量场 MPI 端到端测试

- **任务**：验证 mod_coupling_exchange 跨组交换协议的正确性（阶段5/6 验证①）。
- **新增**：`src/coupling/test_coupling_exchange.f90`（np=2，rank0=struct root /
  rank1=uns root，不加载求解器，仅 common+coupling）；4 项测试：
  ① 常量场 s→u（协议+链内 struct_to_SI 转换对独立手算 SI 期望值）；
  ② 斜坡场 u→s（逐面值保序，查索引错位）；
  ③ struct 侧 0 面（count=0 不死锁、对端计数正确）；
  ④ uns 侧 0 面（③的镜像）。
- **Makefile**：新增 `coupling_test` 目标（common+coupling，MPI 树）。
- **结果**：`mpirun -np 2 bin/coupling_test` → 4 tests, 0 failures, PASS。

---

## 2026-10-05 | 阶段6：结构侧接入 + MPI 跨组交换

- **任务**：接入真实 struct 求解器驱动并实现跨 MPI 进程组数据交换。
- **新增**：
  - `src/structured/mod_struct_driver.f90`：struct_solver_init(comm,ctlfile)/
    step/extract_iface/set_iface_bc；守恒量↔原始量互转
    （T=p·γ·Ma²/ρ；ghost 写 U=(ρ,ρu,ρv,ρw,p/(γ-1)+½ρ|u|²)）；
    init 内设 my_id/Total_proc、`Struct_Comm=comm`、附 Bsend buffer、
    symlink control.ec。
  - `src/coupling/mod_coupling_exchange.f90`：count+payload 两段协议
    （tags 110/111 s→u、210/211 u→s，支持 nfaces=0 空缓冲）。
- **修改**：
  - `mod_interface.f90`：Interface_FACE_TYPE 增 ic/jc/kc + ig/jg/kg。
  - `mod_struct_grid.f90`：register_bc_interfaces 按 face 分派填 inner/ghost。
  - `mod_struct_global.f90`：Global_Var 增 `Struct_Comm=MPI_COMM_WORLD`。
  - **6 个 struct 库文件 64 处 `MPI_COMM_WORLD→Struct_Comm`**
    （fdm/init/solver/io/mpi）——修复 np=2 时 struct collectives 等待
    uns rank 参与 COMM_WORLD 导致 init 挂起；structured/main.f90 保留
    COMM_WORLD，standalone 行为不变（默认值）。
  - `main.f90`：uns 组 `has_struct=(nproc>1)` 包住交换+设BC
    （修复 nproc=1 uns 自等 Recv 死锁）；nproc==1 时 n_struct_ranks=0。
- **验证**（grid_BC，cwd=grid_BC）：np=2 × 3 耦合迭代 exit=0、双向各
  250 面、struct step 正常（grid_BC_smoke_np2.log）；np=1 不挂起；
  standalone struct 回归残差一致。uns 侧 iface|u| NaN 为 unMesh.cas
  mm 单位遗留问题，非耦合引入。
- **遗留**：界面插值未接入驱动（s→u 逐面直设、u→s 简单平均占位）；
  rho=1.0 占位；n_couple/n_uns_steps 硬编码；界面收敛检查；mm→m 缩放；
  grid_BC 真实耦合物理验证。

---

## 2026-10-05 | 阶段6：弱耦合驱动框架 main.f90

- **任务**：实现混合结构/非结构弱耦合驱动主程序（plan 阶段6）。
- **新增**：
  - `src/main.f90`：MPI_Init → `MPI_Comm_split` 拆
    STRUCT_GROUP(rank0) / UNS_GROUP(其余)；读 mix.control；
    struct 组占位驱动（循环 barrier）；uns 组
    init→set_interface_vel/p→step(n)→extract_iface→barrier。
  - `src/unstructured/mod_uns_driver.f90`：uns_solver_init/step/extract_iface
    封装，镜像 main_uns 初始化流程。
- **修改**：
  - `mod_uns_control.f90`：BC_INTERFACE=9。
  - `mod_uns_bc.f90`：bc_t 增 iface_vel(3,nfaces)/iface_p(nfaces)；
    build_bc 把 interface zone 面归入 BC_INTERFACE 组；bc_face_vel/p 读
    iface 数组；新增 set_interface_vel(bcs,faces,u(3,:))/set_interface_p。
  - `mod_uns_simple.f90`：simple_run 增 optional nsteps 上限。
  - `Makefile`：BOOTSTRAP_ORDER 扩展（precision→constants→reference_state
    →interface）；UNS_LAYER6_S / UNS_LAYER9_M = mod_uns_driver 及依赖边。
- **验证**：`make mpi` 通过；cavity(rb_Ra1000) 2 ranks × 3 耦合迭代
  exit=0，SIMPLE mass-imbal 1e-13；unMesh.cas interface 识别 250 面正常
  （mm 单位致数值发散，已知遗留）；cavity 回归 697 步收敛无破坏。
- **遗留**：struct 侧真实求解器占位；MPI 跨组交换用 barrier 占位；
  extract_iface 的 rho=1.0 占位；n_couple/n_uns_steps 硬编码。

---

## 2026-10-05 | 阶段5：数据交换与量纲统一

- **任务**：实现结构↔非结构耦合界面的数据交换与量纲统一（plan 阶段5，
  constraints §2 量纲铁律：各自转 SI → 交换 → 转回）。
- **新增模块**：
  1. `src/common/mod_reference_state.f90`：`reference_state_t`（rho_ref/
     T_ref/L_ref/u_ref + 派生 a_ref/p_ref/p_scale），默认海平面大气；
     `read_mix_control(filename)` 解析 mix.control（key=value，`#` 注释；
     文件不存在不报错，缺省键用默认值）。mod_constants 的编译期常量保留
     作默认值来源。
  2. `src/coupling/mod_interface_units.f90`：基于 reference_state 的结构
     无量纲↔SI 互转——`struct_to_SI`/`SI_to_struct`（ρ,u,v,w,T,p 全套）、
     `struct_vel_to_SI`/`SI_vel_to_struct`（速度向量）。约定：
     ρ*=ρ/ρ_ref, u*=u/a_ref, T*=T/T_ref, p*=p/(ρ_ref·a_ref²)。
  3. `src/coupling/mod_interface_exchange.f90`：基于阶段4 peer_id/peer_w
     的保守插值。`iface_state_t`（与 Interface_List 平行的 SI 状态数组 +
     struct 顶点状态）；`struct_to_uns_interpolate`（双线性 Σpeer_w·顶点
     状态，Σw=1）；`uns_to_struct_interpolate`（面积加权平均 Σ(area·q)/
     Σarea，守恒）。纯数组函数，无 MPI 依赖，可单元测试。
  4. `src/coupling/test_units_exchange.f90` + Makefile `units_test` 目标：
     3 测试全 PASS——量纲往返 <1e-14、常量场插值精确、面积加权正确。
- **MPI 交换**：阶段5仅实现插值核心；MPI_Comm_split 子域 + 跨域 Send/Recv
  包装留阶段6（那时才有真实子通信域与数据流）。
- **构建/回归**：struct_solver/uns_solver/uns_solver_mpi/units_test 全部
  编译通过；cavity Re100 518 步收敛（与之前一致），VTU 仅 z 分量 1e-20
  codegen 噪声差异。
- 涉及文件：`src/common/mod_reference_state.f90`、
  `src/coupling/{mod_interface_units,mod_interface_exchange,test_units_exchange}.f90`、
  `Makefile`（MAIN_UNITS_TEST_F / UNITS_TEST_OBJ_S / bin/units_test）。
  手册 N/A（docs/程序使用手册.tex 不存在）。

---

## 2026-10-05 | 槽道验证：mdot-inlet 启动稳定性 + 新增 outflow 出口边界（cases/channel）

- **任务**：① 验证 mass-flow-inlet 在纯流体槽道的启动稳定性（上轮遗留：LTNE
  纯流体启动发散曾归因存疑）；② 诊断 pressure-outlet 出口段伪调整；③ 用户
  拍板实现 Fluent 风格 `outflow`（充分发展出口）边界。
- **启动稳定性**：Re100 从零场 it=1 du_max=2.84e-4≈1.8U、imbal 5e-13，无发散
  单调收敛（807/790 步）；Re1000（10× 冲量）1033 步。LTNE 发散为无壁面网格
  特异拓扑（销钉假质量汇，已修），mdot-inlet 本身稳健。
- **出口诊断**：pressure-outlet 末列 Ucl 1.5→1.71、vmax 0.115U 向心汇聚；
  vinlet/conv_blend=1/强松弛/并行对照逐一排除，定性为同位网格+Rhie-Chow 在
  固定压力边界的离散不动点。correct_fields 修复（p′ Dirichlet 面+出口通量
  修正）使 it=1 imbal 3.5e-4→3.6e-13 但不改终场。
- **outflow 实现**（BC_OUTFLOW=8）：`bc = <zone> outflow`；bc_face_vel 零阶
  外推、bc_face_T 同 POUTLET 合并、p 默认零梯度；新增 outflow_mass_sums/
  outflow_mass_scale/outflow_mass_rescale（串行）与 outflow_mass_rescale_mpi
  （局部 sums→MPI_SUM allreduce 3 元组→scale），在 simple/piso/pimple 三驱动
  flux_rhiechow 后各调一次；momentum_assembly 与 POUTLET 合并走迎风分裂；
  temperature_assembly 出流隐式项补 BC_OUTFLOW。m_req 排除 POUTLET/FARFIELD
  （压力 Dirichlet 面自调）。启动期无出流时按面积均匀分配（β 返回 0 标志）。
- **结果**：Re100 末列 Ucl 1.7091→1.6420、vmax 0.1152→0.0402、L1 11.45%→
  5.63%；内部场（x=0.5/50.5）与 pressure-outlet 逐 7 位相同；np2 791 步一致。
  已知限制：末列体速度 +3.5%（1/β 不动点，面通量精确守恒）与 ~9% 中心线凸起
  （一阶零梯度固有）；动量改用滞后缩放通量实验更差已回退。
  Re1000 对照：pressure-outlet 1027 步 vs outflow 1033 步；内部场一致
  （U(50.5)/U=0.9997）；末列 u_max/U 1.691 vs 1.574、L1 14.74% vs 5.36%。
- **回归**：cavity 518 步 Ghia 一致（z 分量 1e-20 噪声）；cylinder far-field
  Cp 均差 0.0031≈平台噪声；porous_Ra10 748 步 Nu=1.07800 打印历程逐字相同；
  LTNE 3344 步 ΔT_f=995.02K 全指标一致。
- 涉及文件：`src/unstructured/mod_uns_control.f90`、`mod_uns_bc.f90`、
  `mod_uns_simple.f90`、`mod_uns_simple_mpi.f90`；
  `cases/channel/{README.md,channel_wall.cas,channel_Re100*.control/log/vtu,
  plot_channel.py,split_channel_zones.py,images/}`。手册 N/A（不存在）。

---

## 2026-10-04 | 1D LTNE 发汗冷却验证 + 修复 4 个求解器 bug（cases/ltne）

- **任务**：用户要求将 /mnt/c/temp/LNTE1D.cas 归档至 cases/ltne 并验证：入口
  2 kg/m²/s、300K 空气，压力出口端面恒热流 2e6 W/m²，多孔 ε=0.3 不锈钢，
  与解析解对比。
- **算例设计**：无 wall zone → 出口端面兼作热流壁（tbc=6 2 2e6，复用前会话
  bc_face_T POUTLET 热流扩展）；h·a=0.3 使相间交换长度 6.15=4.9dx 可分辨；
  解析解 = 双温度指数解（T_f' 特征根 λ±），总量守恒 ΔT_f=q″/(G·cp)=995.0K。
- **结果**：3343 步收敛；能量守恒 −0.001%；T_s 剖面 max 0.074%；T_f 中位
  0.198%；Darcy 压降 5.17e5 vs 理论 5.23e5 Pa。归档 README/control/log/vtu/
  plot_ltne.py/images（三联图全对数程吻合）。
- **修复 4 个 bug**（调试过程暴露）：
  1. `mod_uns_simple.f90:ppe_assembly` / `mod_uns_simple_mpi.f90:ppe_assembly_mpi`：
     有压力 Dirichlet 面（POUTLET/FARFIELD）时不再钉死 cell 1（原实现无条件钉死，
     被钉 cell 连续性行被丢弃 → 永久假质量汇；本算例 cell 1 在入口层，60% 流量
     "消失"、内部 u=0.4u_in、imbal 卡死 4.4e-2、T 爆炸 1e8）。MPI 侧 has_pdir 用
     MPI_ALLREDUCE(MPI_LOR) 全局归约。
  2. `momentum_assembly` POUTLET 分支：ap+=F 改为 max(F,0) 隐式 +
     min(F,0)·u_P 显式（倒流时原实现违反迎风、削弱对角、可致 ICC0 NaN）。
  3. uscale 补 gb%uspeed（串行/MPI ×3 驱动），mass-imbal 归一化恢复正确基准。
  4. **`mod_uns_local_mesh.f90:build_local_mesh` 补拷 czone/cztype/czt**（重大）：
     MPI 本地网格此前完全丢失 cell-zone → 多孔/LTNE 物理在 MPI 下整体缺失；
     铁证：porous_Ra10 np2 修复前 296 步（纯流体解）→ 修复后 748 步=串行。
     另 `mod_uns_bc.f90:build_bc` 增 allow_empty 可选参（MPI 本地空 bc zone 合法，
     主并行入口 main_uns_mpi 传 .true.；串行/全局保持严格报错）。
- **回归**：porous_Ra10 串行位级一致（修复后两次验证）；cavity np2 正常；
  LTNE np2 与串行迭代历程一致（du_max=1.699@it1 等逐步吻合，末位分区噪声）。
- **已知限制**：纯流体（无多孔阻力）mdot-inlet+出口算例启动期发散（it=1 速度
  修正 V/apc·∇p′ 放大 ~1e4，apc 仅粘性量级 → NaN）——候选后续：BC 渐启/
  瞬态启动。多孔/含阻力算例受 Darcy 对角保护。
- 涉及文件：`src/unstructured/mod_uns_simple.f90`、`mod_uns_simple_mpi.f90`、
  `mod_uns_bc.f90`、`mod_uns_local_mesh.f90`、`main_uns_mpi.f90`、
  `cases/ltne/{README.md,ltne1d.control,ltne1d.log,LNTE1D.cas,LNTE1D.vtu,plot_ltne.py,
  images/ltne1d_profiles.png}`、memory-bank 3 文件。
- 手册：docs/程序使用手册.tex 仍不存在，N/A（mass-flow-inlet 可选 Tin、tbc 于
  出口 zone 等新参数待手册创建后补录）。

---

## 2026-10-04 | 算例归档 + 测试报告（cases/natconv，用户约定的长期模式）

- **路径变更**：初建在 ~/cases/natconv，当日用户要求迁入项目内
  **MixNSSolver/cases/natconv**（97MB）；迁后处于沙箱可写区，求解器可就地运行，
  原 /tmp 工作区方案作废。
- **归档结构**：cases/natconv/{README.md 测试报告, gen_cas.py, plot_fields.py,
  plot_nu.py, 三套网格, fluid/ 6 算例, porous_darcy/ 3 算例+网格收敛,
  porous_extra/ 6 附加验证, images/ 3 张图}。
- **纯流体 RB 算例重建**：/tmp 清空导致 fluid 控制/log/VTU 全失，从 memory-bank
  参数重建 6 个控制文件并全部重跑，结果逐位复现（Ra=1e3/3e3/5e3/1e4/5e4 →
  Nu=0.9960/0.9960/1.6746/2.1601/3.3139；80² Nu=2.1448/8264 步）。
- **结果图**（matplotlib，plot_fields.py/plot_nu.py，已随归档）：
  fluid_Ra1e4_fields.png、porous_RaK1000_fields.png、nu_comparison.png。
  VTU 解析要点：连接记录="4 n1 n2 n3 质心id"（首列是节点数），六面体拆 12 个四面体、
  单元数据逐六面体复制 12 份；取第 5 列 id 即母六面体质心。绘图用双线性溅射分箱
  （最近邻分箱在粗网格上出棋盘格）。
- **测试报告**：cases/natconv/README.md——设置、两张结果表+热平衡、附加验证表、
  图片、复现命令、结论。
- **用户约定（长期有效）**：以后每个完成的验证算例都按此模式归档到
  <项目>/cases/<主题>/ 并写 README 报告+结果图；已写入 trae 项目记忆
  project_memory.md 的 "Validation Case Archiving Convention"。
- **环境**：迁移前 ~/cases 在沙箱白名单外，需 /tmp 工作区+拷回；迁移后 cases/
  在项目内即可就地运行（后台任务亦可写）。
- 涉及源码：无。报告/图/脚本在 <项目>/cases/natconv。

---

## 2026-10-04 | 多孔介质自然对流验证（Darcy/Forchheimer + LTE，未改源码）

- **标准 Darcy 侧壁加热方腔基准**（/tmp/natconv，cavity.cas 40×40×1 hex 挤出=2D；
  左壁 301K 热/右壁 300K 冷/上下绝热；eps=1 使实现退化为标准 Darcy：阻力 μ/K、
  Brinkman μ_eff=μ、k_eff=k_f；K=1e-6 即 Da=1e-6，Brinkman 可忽略；
  rho=cp=1，μ=k_f=√(gβΔT·K/Ra_K)，gβΔT=3.2667e-2）：
  - Ra_K=10：Nu=**1.078**（文献 Walker&Homsy/Bejan ≈1.07），umax=1.9e-4。
  - Ra_K=100：Nu=**3.096**（40²）/3.097（80²）（文献 3.10–3.16）。
  - Ra_K=1000：40² 12.481→80² 13.286→**160² 13.405**（增量比 0.148≈1/4，二阶收敛），
    Richardson 外推 **13.45**；vs Mahmud&Fraser 13.64（−1.4%）、Baytaş&Pop 14.06（−4.3%，
    文献间本身差 ~3%）；2252 步收敛；冷热壁热平衡误差 0.004%。
  - VTU 确认边界层环流：热壁上升/冷壁下降、核心稳定分层、薄壁面温度层。
- **附加验证（全部通过，无代码改动）**：
  - 底部加热多孔层次临界 Ra_K=10 < 4π²=39.48：Nu=0.996、umax=1.3e-4（解析纯导热 Nu=1）。
  - 超临界 Ra_K=100（t_pert=0.15）：Nu=2.13，稳态对流环形成。
  - LTE 有效导热解析：eps=0.4、k_s=5k_f → Nu=3.407 vs 解析 3.4（0.2%）。
  - **eps 权重严格等价检验**：eps=0.4 配 K=K0/eps=2.5e-8、k_s=k_f（K0=1e-8 使
    Brinkman 占比 K/dx²=1.6e-5），Nu/umax 与 eps=1 基线相差仅 1e-5/5e-5 →
    Darcy 阻力(μ/eps/K)、LTE 混合导热、壁面热流三处 eps 权重全部正确。
  - **Forchheimer 路径**：Ra_K=100 加 inertial=1e4 → Nu 3.096→2.76、umax −25%，
    二次阻力单调抑制符合物理。
- **重要环境事件**：会话间 /tmp 被沙箱清空过一次（算例全丢，后台任务被杀）；
  从 /mnt/c/temp 源 cas + 仓库内记录重建；gen_cas.py 与全部控制文件已重建且
  40/80 网格结果逐位复现。算例已于 2026-10-04 迁至 ~/cases/natconv 持久保留
  （/tmp/natconv 副本仍在但随时可能被清）。
- 涉及源码：无（本阶段仅验证）。算例/控制/日志/VTU 在 ~/cases/natconv（含
  gen_cas.py 结构化 N×N×1 CAS 生成器、porous_Ra*.control 共 8 个）。手册 N/A。

---

## 2026-10-08 | Rayleigh-Bénard 自然对流验证（能量方程+浮力耦合）

- **稳态 SIMPLE 收敛判据补温度项**：原判据仅 mass-imbal + du_max，浮力耦合算例
  速度先停滞而温度场仍在演化导致提前退出（Ra=1000 曾误报 Nu=1.115）。
  simple_run/simple_run_mpi 增 Tprev/dT_max（每外迭代 T 最大变化量，MPI Allreduce MAX），
  三条件同时满足才判收敛，打印行同步加 dT_max；PISO/PIMPLE（时间推进）不动。
- **算例**（/tmp/natconv，cavity_natural_convection_porous.cas，1m³，40×40×1 hex
  挤出+z 向对称面=2D；热底 y=0 T=301，冷顶 y=1 T=300，侧壁绝热，Pr=0.71，
  Boussinesq，β=3.33e-3，cp=1，μ,k 按 Ra 配比）：
  - Ra=1000、3000：Nu=0.996，max|u|<5e-4，纯导热（次临界，解析 Nu=1）✓
  - Ra=5000：强扰动(t_pert=0.15)→Nu=1.675 稳定对流；弱扰动(0.02)→陷导热支。
    离散起对流 Ra∈(3000,5000)（连续理论刚性方形腔 ~2585/无限层 1708）。
  - Ra=1e4：Nu=2.160（max|u|=0.048）；80×80 加密网格（自编结构化 CAS 生成器
    /tmp/natconv/gen_cas.py）Nu=2.145（8264 步收敛），网格变化仅 0.7%。
  - Ra=5e4：Nu=3.314。
  - Hollands(1976) 无限水平空气层关联式：Ra=5e3/1e4/5e4 → 1.90/2.39/3.44；
    本算例偏差 −12%/−9.6%/−3.7%，系统偏低源于宽高比1侧壁无滑移约束（关联式适用大宽高比）。
  - 冷热壁热平衡误差 0.4–1.4%；VTU 显示单一对流环（右侧热流体上升、左侧冷流体下降，
    温度场倾斜），物理结构正确。
- **注意**：稳态求解器+对称初值 sin(πx) 对单环反对称模无直接投影，临界附近需有限幅值
  扰动才能到对流支，属分岔/初值盆现象非代码缺陷。
- **回归**：cavity 稳态与原求解器完全一致（同 518 步、Ghia 表打印差 0）；
  MPI(np2) Ra=1000 收敛 700 步正常。
- 涉及文件：src/unstructured/mod_uns_simple.f90、mod_uns_simple_mpi.f90。
  算例/控制文件/日志在 /tmp/natconv（仓库外）。手册 N/A。

---

## 2026-10-04 | 能量方程 + 多孔介质（Darcy-Forchheimer + LTE/LTNE）

- **多孔动量源**：momentum_assembly 加 Darcy-Forchheimer
  S = -(mu_eff/perm)u - rho*inertial*|u|u；Darcy 线性项隐式进对角 ap，
  Forchheimer 非线性项显式进 rhs（lagged 速度）。mu_eff = mu/porosity（Brinkman），
  内面 lf 插值、边界面取属主格。
- **孔隙率修正输运系数**：温度方程扩散用 k_eff=eps*k_f+(1-eps)*k_s，
  瞬态用 (rho*cp)_eff=eps*(rho*cp)_f+(1-eps)*(rho*cp)_s；对流始终用流体 cp。
- **LTNE 双温度方程**：新建 solid_temperature_assembly
  （(1-eps)*rho_s*cp_s*dT_s/dt = ∇·((1-eps)*k_s∇T_s) + h_sf*a_sf*(T_f-T_s)，
  无对流，共享壁面热 BC）；流体温度方程在 LTNE 模式用 eps*k_f 扩散、
  eps*rho*cp 瞬态、加界面源 +h_sf*a_sf*(T_s-T_f)。
- **数据结构**：cell_zone_t 加 k_s/cp_s/rho_s/h_sf/a_sf；ctrl_t 加 thermal_model
  (lte/ltne)；fields_t 加 T_s/T_s_old/T_s_old_old/gts + 每单元格 perm/inertial/
  porosity/k_s/cp_s/rho_s/h_sf/a_sf；setup_porous_fields 映射 zone→cell。
- **驱动/输出**：simple/piso/pimple(_mpi) 在 LTNE 时解 T_f 后解 T_s，历史同步推进；
  VTU 输出 temperature_solid（h_sf*a_sf>0 时）。
- **验证**：cavity 稳态位级一致（max 1.8e-20）；LTE/LTNE/多孔 串行+MPI 冒烟通过，
  T_f≠T_s 证实非平衡；cylinder far-field 无 NaN（差异为前期 BC 修复，非本次回归）。
- 涉及文件：mod_uns_control.f90、mod_uns_fields.f90、mod_uns_simple.f90、
  mod_uns_simple_mpi.f90、mod_uns_output.f90、main_uns.f90、main_uns_mpi.f90。
  手册 N/A（docs/程序使用手册.tex 不存在）。

---

## 2026-10-04 | PISO 非定常验证 + PIMPLE 算法新增

- **PISO 验证**：cavity_transient (BDF2, dt=0.005, 600步) 串行 7 个 VTU 与原版
  位级一致；MPI(np2) 与原版 np2 位级一致。MPI vs 串行差异 max 4.6e-3（分区累加
  固有现象），非回归。
- **PIMPLE 新增**（phase 11）：
  - ctrl_t 增 `pimple`(逻辑, 默认.false.) + `n_outer_iter`(默认 1)；解析 `pimple`、
    `n_outer_iter`。
  - momentum_assembly / temperature_assembly 松弛分支改为
    `if (ts_order>0 .and. .not. pimple)` → PISO 不松弛；否则（SIMPLE/PIMPLE）松弛。
    PIMPLE = 时间项 + alpha_u 欠松弛。
  - 新建 pimple_run（mod_uns_simple.f90）+ pimple_run_mpi（mod_uns_simple_mpi.f90）：
    每时间步 n_outer_iter 轮外迭代，每轮=动量(松弛, apc=ap/alpha_u)+梯度+Rhie-Chow+
    n_correct 次全量压力修正(alpha_p=1)。历史(u_old/u_old_old) 在步首推进一次，
    外迭代只收敛步内非线性。温度方程每步解一次。
  - 主程序分派：transient & pimple → pimple_run(_mpi)；transient & !pimple → piso_run。
- **验证**：
  - dt=0.1 (CFL~1.4)：PISO step 21 发散 NaN；PIMPLE(alpha_u=0.5) 稳定出物理解。
  - n_outer_iter=2 正常运行，与 n_outer=1 结果差 0.05（非线性收敛改善）。
  - 回归：cavity/cylinder 稳态位级一致；cavity PISO 位级一致；PIMPLE 串行+MPI 物理合理。
- 涉及文件：mod_uns_control.f90、mod_uns_simple.f90、mod_uns_simple_mpi.f90、
  main_uns.f90、main_uns_mpi.f90。手册 N/A。

---

## 2026-10-04 | 阶段 4 完成：交界面几何匹配 + 双线性插值权重

- Interface_FACE_TYPE 扩展：增 `nv`、`verts(3,nv)`、`peer_w(:)`（mod_interface.f90）。
- 结构侧 register_bc_interfaces 改逐单元面登记：每个 interface 单元面 1 条条目、4 顶点 CCW、
  逐面算 centroid/normal/area/bbox，报告按 (block,face) 聚合（grid_BC 1→250 条）。
- 非结构侧 register_interface_zones 每条面补 verts（来自 m%f%nodes）。
- 新建 src/coupling/mod_interface_match.f90（match_interfaces）：投影到 struct quad 平面 +
  Newton 反演双线性 (u,v) + 取 [0,1]² 内法向距最小者，存 4 顶点权重 w=[(1-u)(1-v),u(1-v),uv,(1-u)v]，
  双向 peer_id（struct→uns peer_w 留空，阶段 5 面积加权平均）。
- 新建 src/coupling/test_match.f90 + Makefile match_test 目标（struct MPI + uns MPI + coupling 联合）。
- grid_BC 验证：250/250 双 100% 覆盖、max 法向投影距 1.03e-38（机器零）、max|Σw−1|=0。
- 回归：cavity/cylinder 串行位级一致；cavity np2 仅 epsilon 噪声（无 interface zone 时登记早退）。
- 涉及文件：mod_interface.f90、mod_struct_grid.f90、mod_uns_geometry.f90、
  mod_interface_match.f90（新）、test_match.f90（新）、Makefile。
- 手册 N/A（docs/ 无 程序使用手册.tex）。

---

## 2026-10-03 | 步骤 C 物理参数确认（用户拍板）

- 非结构侧=300K、1atm 空气：unMesh.control 改 rho=1.177 kg/m³、mu=1.846e-5 Pa·s；
  重跑验证进风速度 mdot/rho=1.699 m/s 打印正确。
- 结构侧层流/AoA=0、Re=1e4/mm 按 Ref_L=1mm 直读，均确认无需改动。
- 删除项目根目录重复文件 bc3d.inp、unMesh.cas（仅留 grid_BC/ 原件）。
- 记录遗留：mm→m 面积/体积缩放（×1e-6）待阶段 5 units/exchange 层统一。

---

## 2026-10-03 | 阶段 3 步骤 C：真实耦合网格贯通 + CAS 体区域落位 + mass-flow-inlet

- 输入：grid_BC/{Mesh3d.x(二进制 Plot3D 3 块), bc3d.inp(generic:8 在 block2 j=1), unMesh.cas
  (Pointwise wedge 网格；体 zone2 fluid；面 zone4 sym/5 wall/6 pressure-inlet/7 interface)}；
  单位均 mm；结构 Ma=3, Re=1e4/mm, Tinf=108K, Tw=300K；非结构 mdot=2 kg/m²/s、绝热壁。
- 网格文件名：结构化源码 Mesh3d.dat→Mesh3d.x（mod_struct_init/grid/io/mpi，12 处）。
- 体区域：mesh_t 增 czone/nczone/czt/cztype；sec_cells 记录体区域号（grow_iarray 扩容）；
  finalize_zones 按面/体引用拆分 zone 表；新增 resolve_cell_zones（id/name 匹配 fluid/porous，
  默认 fluid，未知报错），main_uns/main_uns_mpi 在 read_control 后调用。
- build_bc 豁免 interface zone 边界面（zone_is_interface，与登记规则一致），不再强制 bc。
- 新增 BC_MASSINLET=7（mass-flow-inlet mdot）：bc_spec/bcgroup 增 mdot/uspeed=mdot/rho，
  bc_face_vel uf=-uspeed·n_out，动量定速度装配、温度、PPE/报告接入。
- 验证：两侧独立读取真实网格成功；界面几何完全重合（y=0, x420-670, z0-50, 质心(545,0,25),
  面积12500, 250 quad，法向相反 struct-y/uns+y）；结构 Ma=3 冒烟推进无 NaN；
  cavity 串行、cylinder(vinlet) 对原始二进制位级一致；cavity np2 仅 epsilon 噪声。
- 已知：uns 独立运行无质量出口→PPE 奇异 NaN（阶段 5 交换解决，预期）。
- 文件：mod_uns_mesh/mod_uns_cas_reader/mod_uns_control/mod_uns_bc/mod_uns_simple/
  mod_uns_geometry/main_uns/main_uns_mpi + 结构化 4 文件；grid_BC/{control.ec,unMesh.control}。
- 手册：MixNSSolver docs/ 无 程序使用手册.tex，暂 N/A。

---

## 2026-10-03 | pressure-far-field BC 零解 bug 修复

- 现象：cylinder_invis 算例 it=1 即"收敛"，lin-it=0，max|u|=0（原始 UNSSolverProj 与
  迁移代码同陷此 bug）。
- 根因：mod_uns_bc.f90 的 bc_face_vel 中 BC_FARFIELD 以**内场** uP·n 判入/出流；
  零初场 un=0，上游入流半球被误判为出流（uf=uP=0），边界零通量 → 零解静止不动。
- 修复：改为以外部自由流 u_far·n 判定（u_far·n<0 施加自由流，否则零阶外推）；
  bc_face_vel 是动量装配/Rhie-Chow 通量/梯度/限制器的共用入口，一处修复全局生效。
- 验证：cylinder(far-field) 恢复物理解，CD=1.62（文献1.5）、分离角±53°（55-60°）、
  尾迹2.34D、Cp 与旧基线平均偏差0.004；残差平台~9.8e-6（固定p远场的弱相容残差，非bug）。
- 无回归：cavity 串行位级一致；cylinder(vinlet) 与原码位级一致；cavity np2 仅 z 向 epsilon 噪声。
- 涉及文件：src/unstructured/mod_uns_bc.f90（1 处分支）。手册 N/A（bug 修复，无新参数）。
- 备注：仅修 MixNSSolver 树；UNSSolverProj 工作树未改动。

---

## 2026-10-03 | 阶段 3 步骤 B：CAS interface zone 登记 + .control cell_zone 扩展

- CAS interface zone 识别：`mod_uns_geometry.f90` 新增 `register_interface_zones(m, g)`，
  扫描 zone 表中 `user_name`/`cond_name` 含 "interface"（大小写不敏感）的 zone，
  将其所有面登记进 `Interface_List`（solver=PEER_UNS, block_no=zone_id,
  centroid/normal/area 取自 geom_t, bbox 取自面节点）。
- 调用位置：`main_uns.f90` 和 `main_uns_mpi.f90` 在 `compute_geometry` 之后；
  MPI 仅 rank 0 在全局网格上登记（partition 之前）。
- 设计决策：登记逻辑放在 mod_uns_geometry（已有 mesh+geom），而非 mod_uns_cas_reader
  （需额外引入 geom 依赖，破坏分层）。
- `.control` 扩展：`cell_zone = <id|name> fluid|porous [perm=.. inertial=.. porosity=..]`；
  新增 `cell_zone_t` 类型、`CZ_FLUID=1`/`CZ_POROUS=2`、`MAXCZ=32`、`parse_cell_zone` 子程序；
  首 token 整数→zone id，否则→zone name；多孔系数仅解析存储（动量汇未实现）。
- 编译：`make unstructured unstructured_mpi -j4` 0 error。
- 回归：cavity 串行 VTU 位级一致；cylinder(velocity-inlet) VTU 与原求解器位级一致；
  cavity MPI(np2) x/y 位级一致（z 分量 ~1e-19 机器 epsilon 噪声，partition 相同，
  由新增代码影响编译器代码生成所致，非数值回归）。
- 发现：cylinder pressure-far-field BC 在当前原始+迁移代码中均收敛于零解（lin-it=0），
  疑似预存 BC bug，与步骤 B 无关。
- 涉及文件：src/unstructured/mod_uns_geometry.f90、mod_uns_control.f90、
  main_uns.f90、main_uns_mpi.f90。
- 手册：`docs/程序使用手册.tex` 尚不存在；cell_zone 参数待手册创建后同步。
- 下一步：阶段 4（交界面几何匹配）。

---

## 2026-10-03 | 阶段 3 步骤 A：非结构求解器迁移完成，回归位级一致

- 勘察：UNSSolverProj 23 个 .f90（7936 行）已是模块化代码，与阶段 2 的"封装自由子程序"性质不同。
- 迁移：21 个 .f90 迁入 src/unstructured/（mod_precision 删除，统一用 common 层）；
  main.f90→main_uns.f90，main_mpi.f90→main_uns_mpi.f90。
- 重命名：所有 mod_*→mod_uns_*；mod_precision use 语句改为 `use mod_precision, only: dp, ip, pi`；
  mod_mpi→mod_uns_mpi_core（避免与文件名后缀冲突）。
- Makefile：unstructured/unstructured_mpi 目标；分层 PRISTINE 边（S0-S5 串行 / M0-M8 MPI）；
  MPI 依赖模块（gather/halo/restart/local_mesh/partition/mpi_core）串行排除；
  干净 `make unstructured unstructured_mpi -j4` 0 error，双树通过。
- 排错：mod_uns_fields 依赖 mod_uns_bc 需分层；mod_uns_linsolver_mpi 依赖 mod_uns_halo 需分层；
  MPI 依赖模块补 `use mpi`（原 UNSSolverProj 中这些模块隐式依赖 use mpi 传递）。
- 回归：cavity/cylinder 串行 + MPI(np2) 与基线 cmp 位级一致；Ghia 对比、Cp、分离角 stdout 一致。
- 涉及文件：src/unstructured/ 22 个文件、Makefile。
- 手册：纯迁移、无新增物理模型/控制参数，`docs/程序使用手册.tex` N/A。
- 下一步：阶段 3 步骤 B（CAS interface zone 识别 + 登记进 mod_interface + .control 扩展 cell_zone）。

---

## 2026-10-03 | 阶段 2b-s4b：新增 BC_Interface 耦合交界面类型，回归位级一致

- 设计输入（用户拍板）：纯标记/注册类型；与 BC_MSG_TYPE 并存登记（不动旧类型）；
  结构面与非结构 zone 用几何自动匹配；SI 转换留给 Phase 3 独立耦合层；
  类型放 common 层新文件；几何信息 2b 就算好。
- 关键发现：原始 OpenCFD-EC 无 bc=8（generic）常量与分派分支，bc=8 会落 else 报错；
  交界面是 MixNS 新增类型。
- 落地：
  1. 新建 src/common/mod_interface.f90：BC_INTERFACE=8、PEER_STRUCT/UNS、
     Interface_FACE_TYPE（登记身份 + centroid/normal/area/bbox 几何 + match_state/peer_id
     配对状态）、Interface_List/Num_Interface。
  2. mod_struct_grid 加 register_bc_interfaces（含 iface_node/iface_bbox/iface_tri
     三个内部子程序）：扫描 Mesh(1) 的 bc_msg，bc==8 的面登记索引并按 face 方向
     用节点坐标做面积加权几何计算（quad 拆两三角形）。
  3. mod_struct_init 的 init 里 update_Mesh_Center(1) 后调用 register_bc_interfaces。
  4. mod_struct_bc 边界分派（非叶轮机 + 叶轮机两分支）加 bc==BC_INTERFACE 占位分支：
     2b 阶段 cycle 跳过（Phase 3 耦合前不施边界条件），避免 bc=8 落 else 报错。
  5. Makefile 加 common→structured 依赖边（mod_interface 先于 structured 编译）。
- 调试：修复 a2/b2 循环变量漏声明 + 内部子程序参数 b 与 B（Block_TYPE）在
  Fortran 大小写不敏感下重复（b→bv）。
- 验证：rm -rf build bin 干净 `make structured structured_mpi -j4` 0 error；
  M6 np1 回归 5 二进制 + 6 文本位级一致；run.log 无新增输出（无交界面时静默）。
- 遗留：几何计算（centroid/normal/area/bbox）在 M6 上未执行（无 bc=8），其正确性
  需 Phase 3 有交界面算例时验证。
- **至此阶段 2b 全部完成**（B1-B4 封装 + s3 收窄 + s4a f2008 + s4b BC_Interface）。

---

## 2026-10-03 | 阶段 2b-s4a：structured 树收紧 -std=f2008，15 处机械修复，双回归位级一致

- 评估（仅试编译不改源码）：f2008 全量违规 15 处、5 文件、4 类，无算法改写。
- 修复（全部机械、位级零风险）：
  - A 混合 kind 字面量 ×6：mod_struct_time L29/141/170/199、mod_struct_solver L2569/2952，
    `1.0`→`1.0_PRE_EC`、`4./3.`→`4._PRE_EC/3._PRE_EC`、`10.0`→`10.0_PRE_EC`。
  - B FORMAT 缺逗号 ×6：mod_struct_bc L867/977、mod_struct_grid L742/890、
    mod_struct_io L408/1183，补 `,`（格式串等价）。
  - C isnan ×2 调用点（3 处）：mod_struct_grid L698、mod_struct_solver L1854
    `Isnan`→`ieee_is_nan`，所在子程序（check_mesh_quality_onemesh、comput_max_Res_onemesh）
    加 `use, intrinsic :: ieee_arithmetic`。
  - D OPEN access=append ×1：mod_struct_grid L725 → `position="append"`。
- Makefile：删除 structured 两树 pattern 规则的 `-std=legacy`，与项目其余部分统一
  `-std=f2008`（FFLAGS 行 74 不变）。
- 验证：rm -rf build bin 后干净 `make structured structured_mpi -j4` 一次通过 0 error；
  M6-wing np1（ser）+ np2（mpi）双回归：flow3d/SA3d/wall_dist/partation-auto/part_grid
  二进制 cmp 一致；文本产物一致；仅 output_para.out a0/d0 已知良性差异。
- 遗留（不阻塞，另开任务）：~2700 条 `-Wtabs`（源码 tab 缩进，纯外观）、22 条
  unused-variable。
- 至此 2b 仅余 s4b（BC_Interface 类型设计，gridgen generic:8 需与用户讨论）。

---

## 2026-10-03 | 阶段 2b-s3：Global_Var 审计 + 13 个单消费成员搬迁，回归位级一致

- 审计（脚本逐成员统计跨模块使用）：Global_Var 84 成员，9 模块 use。
  关键结论：Fortran 中 private 实体禁止其他模块 use 关联，而 Global_Var 自身无过程，
  故不能原地标 private；收窄唯一途径是把"单一消费模块"成员搬到拥有模块。
- 搬迁 13 个（用户批准）：
  - FD_Flux/FD_scheme → FDM_data（mod_struct_fdm_data.f90）；仅 read_para_FDM 与
    Residual_FDM 包装子程序引用全局副本（二者均已 use FDM_data）；
    Residual_FDM_local/fp/fm 内同名项是哑参，不受影响。
  - Istep_average=0 → mod_struct_io 模块级（init_average/Time_average/output_flow_average
    全在该模块内）。
  - Iflag_init/Kstep_init_smooth/NUM_THREADS/IF_Walldist/Ref_medium_usrdef/Pre_Step_Mesh(3)
    （integer）与 Aos/Turbo_P0/Turbo_T0/Turbo_L0（real，模块头 use precision_EC）
    → mod_struct_init 模块级。
- 安全性预检：接收模块内无同名局部/哑参遮蔽（io/init）；每个被含子程序均有 implicit none，
  杜绝隐式名误绑宿主变量。
- 审计记录的遗留点（本轮不动）：71 个成员跨 2+ 模块须保持 public；全局标量与
  Mesh_TYPE 组件双重表示（NVAR/Iflag_Scheme/IFlag_flux/IFlag_Reconstruction/
  Iflag_turbulence_model/Bound_Scheme/Num_block），依赖 init 时全局→per-Mesh 拷贝约定；
  非 mod_struct_ 前缀模块名（Global_Var/const_var/Flow_Var/FDM_data/Type_def1/
  Wall_dist/filting_Var）改名属纯 churn，留待后续。
- 回归（M6-wing np1，101 步，RC=0）：5 二进制 cmp 一致 + 6 文本一致，
  output_para.out 除已知 a0/d0 行外一致。

---

## 2026-10-03 | 阶段 2b-B4：mod_struct_io + mod_struct_init，29 个 sub 全封装完毕

- 新建 mod_struct_io.f90（1525 行）：sub_IO + sub_comput_dw + sub_Post +
  sub_Post_timeAverage；Wall_dist 小模块原样保留在文件顶部。
- 新建 mod_struct_init.f90（1035 行）：sub_init + sub_read_parameter。
- 补丁：sub_init 三处 use（init→io: read_main_Mesh/read_inc/comput_dist_wall；
  Init_flow→io: read_flow_data；init_flow_zero→io: smoothing_oneMesh）；
  main.f90 program use io（comput_force/output_flow/output_vt/Time_average/
  output_flow_average）+ init（read_parameter/init/Init_flow）。
- sub_interfaces.f90 已无任何 use 点（get_U_conner/get_xyz_conner 在 mod_struct_mpi、
  prolongation 在 mod_struct_solver、allocate_mem_Blocks 在 mod_struct_init 均有实现），
  与最后 6 个 sub 一起删除（用户批准 7 文件）。
- Makefile：ST_BASE_NAMES 旧链替换为 6 层真实 use 边链（L0 precision →
  L1 constants/types/fdm_data/flowvar → L2 global →
  L3 scheme/flux/fdm/time/bc/mpi/grid → L4 solver/io → L5 init → main），
  ser/mpi 两树各一组，st_objs/st_objm foreach 生成。
- 验证：rm -rf build bin 后干净 `make structured structured_mpi -j4` 一次通过
  （0 circular、0 error）；增量 no-op 正常；M6-wing np1 回归（101 步，RC=0）
  5 个二进制 cmp 一致 + 6 个文本 diff 一致，仅两处已知良性差异。
- 已知外观问题（2a 起，非新增）：增量构建时 gfortran -MMD 对同文件定义并 use 的
  辅助模块（wall_dist/filting_var/type_def1）生成自引用边，make 报
  "Circular <mod>.mod dependency dropped" 并丢弃自边；干净树首建无 .d 不触发，
  不影响正确性，暂不处理。
- 至此 src/structured 仅 17 文件（main + 16 个 mod_struct_*），29 个自由子程序
  全部模块封装完成。2b 剩余：s3 Global_Var 审计（只加 private/public 不深改）；
  s4 收紧 -std=legacy 评估 + BC_Interface 类型设计（gridgen generic:8 需与用户讨论）。

---

## 2026-10-03 | 阶段 2b-B3：mod_struct_solver 封装，np1 回归位级一致

- 新建 mod_struct_solver.f90（4142 行），封装 8 个文件 46 个子程序：sub_Residual +
  sub_time_advance + sub_turbulence_{SA,NewSA,SST,BL} + sub_limitflow + sub_filtering。
  SCC(6) 循环依赖（Residual↔time_advance↔turbulence）要求必须同模块。
- filting_Var 小模块原样保留在 mod_struct_solver 文件顶部（名字不改，f/f0 模块变量）。
- 删除 sub_time_advance 内 `use interface_defines`（prolongation 已在同模块）；8 文件内
  无同名子程序冲突、无局部变量遮蔽兄弟过程（脚本预检）。
- 调用方补丁：main.f90 program use solver（NS_Time_advance/NS_2stge/NS_3stge/
  Filtering_oneMesh/output_Res）；sub_init 的 init_flow_zero use solver
  （NS_Time_advance/output_Res/prolong_U）。
- 排错：-j4 竞态（main/sub_init 先于 mod_struct_solver.mod 编译，Fatal Cannot open module），
  先 make 该 .o 再全量；旧 8 文件未删时 filting_Var 的 f/f0 链接重复定义，删源后消失。
- 回归（M6-wing np1，101 步，RC=0）：5 个二进制 cmp 一致 + 6 个文本 diff 一致，
  仅 a0/d0 与 run.log 两处已知良性差异。
- 删除旧文件（用户批准）：sub_Residual/sub_time_advance/sub_turbulence_{SA,NewSA,SST,BL}/
  sub_limitflow/sub_filtering 共 8 个 .f90 及两树 .o/.d。
- 下一步 B4：mod_struct_io（sub_IO+sub_comput_dw+sub_Post+sub_Post_timeAverage，
  Wall_dist 小模块原样搬入）+ mod_struct_init（sub_init+sub_read_parameter）；
  删 sub_interfaces；更新 Makefile pristine edges；干净 -j4；回归。

---

## 2026-10-03 | 阶段 2b-B1/B2：7 个 mod_struct_* 模块封装，两批回归均位级一致

- 2b 手法（用户拍板）：最小侵入封装——新文件 = module 头（无模块级 use/implicit none）+ contains
  + 旧 sub 逐字 + end module；子程序体内 use/implicit none 不动；每批回归通过后再删旧文件。
  分 4 批 B1→B4 自底向上，每批跑 M6-wing np1 回归。
- B1（完成）：mod_struct_scheme（sub_scheme）/mod_struct_flux（sub_flux_split）/
  mod_struct_fdm（sub_Finite_Difference1+2）/mod_struct_time（sub_time_acceleraction+sub_LU_SGS）。
  修复：局部 real 变量 minmod/minmod2/minmod4 遮蔽同名模块函数（删局部声明）；显式接口下
  pointer 数组元素传显式 shape 哑参变硬错误（改传整数组，哑参界核实完全一致，未用
  -fallow-argument-mismatch）；旧文件无尾换行致拼接粘连（加空行分隔）。旧 6 文件已删。
- B2（本次完成）：mod_struct_bc（sub_boundary+sub_boundary_user）/mod_struct_mpi
  （sub_partation_mpi+sub_update_buffer_mpi+sub_udate_Meshlink_mpi）/mod_struct_grid
  （sub_convert_inp+sub_geometry+sub_debug）。Type_def1 内嵌模块原样保留在
  mod_struct_grid 之外（文件顶部），convert_inp_inc/Convert_bc 入 mod_struct_grid。
- B2 补丁（外部子程序获得显式接口必须 use）：sub_time_advance 6 个 ns_* 子程序→bc+mpi；
  sub_init 的 init（partation/Update_coordinate_buffer_onemesh/update_Mesh_Center→mpi，
  Comput_Goemetric_var/Output_mesh_debug→grid）、Creat_Mesh（同）、init_flow_zero→bc+mpi；
  sub_IO 的 read_inc→grid（convert_inp_inc）；sub_Post 的 smoothing_oneMesh→bc+mpi；
  main.f90 program→grid（check_mesh_multigrid/set_control_para/check_mesh_quality）。
- B2 排错 4 处：
  1. mod_struct_grid 首次生成时 contains 误置于 Type_def1 之后（convert_inp_inc 漏在模块外），
     修正为 Type_def1 结束 → module mod_struct_grid → contains。
  2. partation 内残留 part 的 interface 块（遮蔽同模块 part 过程，链接报 undefined `part_`），删除。
  3. init/Creat_Mesh 漏 use Update_coordinate_buffer_onemesh（链接 undefined），补入。
  4. 删旧源后旧 .d 文件仍把 type_def1.mod 绑到 sub_convert_inp.f90（make "No rule"），
     删除两树中 8 个已删单元的陈旧 .d；旧 Type_def1 与新模块重复定义链接错误随旧源删除消失。
- B2 回归（M6-wing np1，101 步，t_end=1，RC=0）：flow3d/SA3d/wall_dist/partation-auto/
  part_grid 二进制 cmp 全一致；Residual/Step_mess/bc3d.inc/mesh-quality/force*.log 文本一致；
  仅两处已知良性差异（output_para.out 的 a0/d0；run.log 末尾 mpirun 警告）。
- 删除旧文件（用户批准）：sub_boundary/sub_boundary_user/sub_partation_mpi/sub_update_buffer_mpi/
  sub_udate_Meshlink_mpi/sub_convert_inp/sub_geometry/sub_debug 共 8 个 .f90。
- 涉及文件：src/structured/mod_struct_{bc,mpi,grid,scheme,flux,fdm,time}.f90、
  sub_{Residual,init,IO,Post,time_advance}.f90、main.f90；build/{ser,mpi} 陈旧 .d 清理。
- 下一步：B3 封装 mod_struct_solver（sub_Residual+sub_time_advance+
  sub_turbulence_{SA,NewSA,SST,BL}+sub_limitflow+sub_filtering，SCC(6) 循环依赖必须同模块；
  注意删 sub_time_advance 内 `use interface_defines`，prolongation 同文件）。

---



## 2026-10-03 | 阶段 2a 完成：结构求解器原样搬迁 + 位级回归通过

- 来源：仅以根目录 zip 解压的 `external/OpenCFD-EC-1.16a/`（1.16a 原始版，只读基线）为准；
  外部 `…/OpenCFD-EC-1.16a` 为深度定制版，明确不作来源。
- 搬迁：31 个核心文件 → `src/structured/`（7 个独立网格工具不迁）；29 个 sub_*.f90 经
  GBK→UTF-8 原样拷贝；sub_modules.f90 拆 5 个 mod_struct_*（precision/constants/types/
  global/fdm_data）；opencfd 主文件拆 mod_struct_flowvar.f90 + main.f90；real*8→PRE_EC。
- 关键排错 3 处：
  1. f2008 拒绝遗留写法（pointer 数组元素传显式 shape 哑参等）→ structured 树单独加 `-std=legacy`；
     另加 `-ffree-line-length-none`（50 行超长）。
  2. 薄壳重导出 dp 与旧例程局部变量 dp 撞名；且原模块内 `include "mpif.h"` 使 MPI 常量经
     use 链全局可见 → 薄壳改为 `use mpi`（实体 public）+ PRE_EC 用 selected_real_kind(15,307)
     自含定义、不 import dp。
  3. 基线编译：原始 GBK 源需 iconv UTF-8、全角"！"转半角（列 1 不被认作注释符）、-std=legacy。
- Makefile 改动：structured 在 ser/mpi 两树均用 mpif90（MPI-only；ser 树即 -np 1 构建）；
  %_mpi.f90 过滤仅对 unstructured；ser 树 structured 专用编译/链接规则；
  "PRISTINE-BUILD EDGES" 增 9 个基础模块链；干净 `make -j4 structured_mpi` 通过。
- 回归（/tmp/ocfd_regression，基线=原始源 mpif90 -O3 -std=legacy）：
  - M6-wing（4 块，约 40 万点）np1（t_end=1，106 步）与 np2（t_end=0.3，31 步）；
  - flow3d.dat/SA3d.dat/wall_dist.dat/partation-auto.dat/part_grid.dat 二进制 cmp 位级一致；
    Residual.dat/force*.log/Step_mess.dat/bc3d.inc/mesh-quality.dat 文本一致；stdout 仅计时行
    与多进程打印交错差异。
  - 唯一差异 output_para.out 的 a0/d0：sub_read_parameter 局部变量，仅涡轮分支赋值，
    普通算例未初始化即打印（原版 latent UB，不参与计算）；2a 不修，留 2b。
  - 原版结束不调 MPI_Finalize（mpirun RC=1）为基线固有行为。
- 产物：bin/struct_solver、bin/struct_solver_mpi。
- 涉及文件：src/structured/ 下 36 文件、顶层 Makefile；external/ 只读未动。
- 手册：纯迁移、无新增物理模型/控制参数，`docs/程序使用手册.tex` N/A。
- 下一步：等用户确认进入阶段 2b（模块化封装、去 Global_Var、BC_Interface）。

---

## 2026-10-03 | 阶段 1：顶层 Makefile 完成

- 用户拍板两项设计：非结构块类型扩展 `.control`（`cell_zone = <id|name> fluid|porous`）；
  参考量拆为运行期 mod_reference_state + 拟议 mix.control（均记入 constraints）。
- 调研：参考 UNSSolverProj/Makefile；核对 lib（libmetis.a、libparmetis(_gnu_mpi).a、
  libtecio.a）；环境为 gfortran 13.3.0 + OpenMPI 4.1.6；OpenCFD src 为 37 个 .f90。
- 新建顶层 `Makefile`：
  - 目标 all/mpi/structured(_mpi)/unstructured(_mpi)/common(_mpi)/clean/help，DEBUG=1。
  - build/ser 与 build/mpi 双树（切换串行/MPI 无需 clean）；.o/.d 镜像 src 目录，.mod 集中到 mod/。
  - gfortran `-cpp -MMD -MP` 自动依赖；BOOTSTRAP_ORDER 给 mod_precision 真实依赖边保证干净 -j 构建。
  - 约定：模块名全局唯一；`*_mpi.f90` 仅 MPI 构建；共享文件 `#ifdef HAVE_MPI`。
  - MPI 默认链串行 METIS（rank0 分块，规避 ParMETIS MPICH/OpenMPI ABI 问题）+ TECIO。
- 验证：串行/MPI 编译通过；增量 no-op 正确；接口变更级联重建正确（gfortran 接口未变时
  不重写 .mod，属预期）；干净 `make -j4` 连续 3 次成功；缺失 main 给出友好报错。
- 删除根目录游离的 mod_constants.mod、mod_precision.mod（会干扰模块解析）。
- 涉及文件：`Makefile`（新建）、根目录 2 个 .mod（删除）、memory-bank 4 个文件更新。
- 手册：本次为构建配置、无物理模型/控制参数变化，`docs/程序使用手册.tex` N/A（该手册尚不存在）。

---

## 2026-10-03 | 查明三个待确认问题

- UNSSolverProj 位于 `/home/sundong/Fortran_Project/UNSSolverProj`（git 仓库，23 个 src 模块、Makefile、cases、lib、tests 齐备）；OpenCFD-EC 已解压于 `/home/sundong/Fortran_Project/OpenCFD-EC-1.16a`。
- 现状：非结构侧已有 `.control`（key=value）输入体系与按 zone 号指定 BC 的机制；CAS reader 已解析 face zone 的 cond_name/user_name；尚无 cell zone（体区域/多孔块）表、无 porous 机制。
- 给出建议：块类型扩展 `.control`（`cell_zone = <id|name> fluid|porous`）而非新建文件；参考量由编译期 parameter 改为运行期读耦合控制文件。均待用户拍板。
- 未改动 MixNSSolver 源码。

---

## 2026-10-03 | 初始化 memory-bank

- 通读 `docs/程序功能说明.md`（6 条框架需求）与 `docs/plan.md`（8 阶段计划），熟悉工程。
- 实地核对目录：common 下两个模块已存在且可编译（根目录留有 .mod 产物）；structured / unstructured / coupling 为空；无 Makefile、无 main.f90。
- 创建 memory-bank 五个文件：projectbrief、activeContext、progress、worklog、constraints。
- 涉及文件：
  - `memory-bank/projectbrief.md`（新建）
  - `memory-bank/activeContext.md`（新建）
  - `memory-bank/progress.md`（新建）
  - `memory-bank/worklog.md`（本文件）
  - `memory-bank/constraints.md`（新建）
- 未改动任何源码。

### 待确认
- UNSSolverProj 源码位置未知。
- 下一步方向（Makefile / 解压 OpenCFD 启动阶段 2 / 细读专题规划）等待用户选择。
