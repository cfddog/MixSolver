# bj_porous —— Beavers-Joseph 流体/多孔界面（α = 0 / 1 / 2）

自包含现场：网格 + 4 个 control（主控制 + 三个滑移系数变体）+ 运行/校验脚本。
运行产物落在 `logs/`、`out/`。

## 1. 物理与网格

| 项 | 值 |
|---|---|
| 几何 | 二维开槽道 x ∈ [0, 0.1] m，y ∈ [0, 0.01] m；100 × 40 × 1 六面体（z 一个单元厚） |
| 分区 | 上半 y ∈ [5,10] mm = 流体（zone 2，`upper`）；下半 = 多孔床（zone 3 `bed`，CAS 里带 `VC:porous` 标记） |
| 层厚 | H = 5 mm（每层）；K = 1e-6 m²，ε = 1（→ μ_eff = μ） |
| 驱动 | 均匀体积力 **f = 0.05 N/m³**（x 向）；两端 zone 8 / 9 都是 pressure-outlet p = 0，所以驱动压强梯度恰为 f |
| 流体 | 空气 ρ = 1.177 kg/m³，μ = 1.846e-5 Pa·s（常物性、无浮力 → 无能量耦合） |
| 边界 | zone 5 = symmetry（z 两平面）；zone 6 = wall（y = 10 mm，流体侧）；zone 7 = wall（y = 0，床底）；zone 8/9 = pressure-outlet |

## 2. 参考解（不是经典 BJ 滑移公式）

本算例的参考量是**两层耦合 ODE（流体层 + Brinkman 多孔层）**的解析解，
推导与背景见 `cases/beavers_joseph/README.md` 第 2 节：
经典 BJ 滑移公式在这里**不适用**（它会给出 3.4e-2 ~ 5.5e-3 量级的界面速度，
比数值解高/低数倍，见校验脚本的 `uB_class` 列）。

| 量 | 解析值 |
|---|---|
| Darcy 塞流速度 u_D = K·G/μ | 2.710e-3 m/s |
| 界面速度 u_B（α = 0） | 3.385e-2 m/s |
| 界面速度 u_B（α = 1） | 7.893e-3 m/s |
| 界面速度 u_B（α = 2） | 5.536e-3 m/s |

## 3. 运行

```bash
./run.sh          # 依次跑 α = 0/1/2，np = 1
```

`run.sh` 逐个把 `variants/bj_a{0,1,2}.control` 拷成 `unMesh.control`（工作副本永远是
`unMesh.control`，跑完复位到 α = 1，与归档 `bj_a1` 对齐），日志进 `logs/bj_a*_np1.log`，
VTU 进 `out/bj_a{0,1,2}_np1.vtu`。单核每个 α 约 2 分钟（~810 次外迭代到 CONVERGED）。

校验：

```bash
BJ_PLOT_OUT=out/bj_check_np1.png python3 plot_bj.py \
    0.0=out/bj_a0_np1.vtu 1.0=out/bj_a1_np1.vtu 2.0=out/bj_a2_np1.vtu
```

## 4. 实测（本 bundle，HEAD 二进制，np = 1，均跑到显式 CONVERGED）

```
alpha     uF1_num     uF1_ref      uD_num      uD_ref    uB_class   rms_err
   0.0  8.5512e-03  8.4977e-03  2.9946e-03  2.7086e-03  3.3857e-02     0.22%
   1.0  1.1768e-02  1.2117e-02  2.9180e-03  2.7086e-03  7.9000e-03     0.94%
   2.0  1.0038e-02  1.0447e-02  2.9592e-03  2.7086e-03  5.5402e-03     1.20%
```

* `uF1` = 界面第一个流体单元中心的 u_x，`uF1_ref` = 两层 ODE 在同一位置的取值，
  `rms_err` = 流体层整条剖面的 RMS 相对误差（`plot_bj.py` 计算）。
* 归档值（`cases/beavers_joseph/README.md` 第 4 节，运行被 800 步上限截断）为
  α = 1：uF1 = 1.173e-2，RMS 1.09% —— 与本 bundle 的 1.1768e-2 / 0.94% 相差 < 0.3%；
  差异来自归档跑在 800 步上限处（`du_max` 还有 3.1e-9），而本 bundle 的 α=1 多跑到 882 步
  才停（`du_max < 1e-10`，更收敛）✓
* `uD_num ≈ 2.92 … 3.00e-3` 略高于 Darcy 值 2.7086e-3 属预期（床内 Brinkman 黏性修正，
  见归档 README）。
* 三个 α 的外迭代步数（`run.sh` 顺序 0 → 1 → 2）：**α=0 @754、α=1 @882、α=2 @809**。
  收敛快慢与精度无关（α=0 的 RMS 最小却收敛最快）。
* **tarball → 编译 → 跑现场**这条链复核过：用 `--with-repo` 打出的源码在干净目录里
  `make all -j1` 重新编出 `bin/mixnsolver`，跑 α=1 得到与本表**完全相同**的 it=1 行
  （`mass-imbal = 3.110E-18, du_max = 2.120E-04, lin-res = 1.36E-09`）、
  `CONVERGED after 882` 以及 `uF1 = 1.1768e-02 / RMS 0.94%` —— 说明现场对二进制版本不敏感
  （同一份源码同一编译器）。

## 5. 与归档 `cases/beavers_joseph/*.control` 的两处差异（重要）

1. **热边界被关掉**（没有 `tbc`，也没有 `init_t`）。本算例是等温的（常物性、无浮力），
   但求解器仍然会输运一个被动 T 场，而收敛判据要求
   `mass-imbal AND du_max AND dT_max` 同时小于 `outer_tol`
   （`src/unstructured/mod_uns_simple.f90` 第 235 行）。归档 control 只把 zone 7
   （床底）钉在 300 K，流体侧壁面默认 0 K，于是 T 一路漂到 53 K；把**两面墙都钉 300 K
   也不行**——多孔床会把场拉到 ~88 K 的稳态（实测 dT_max 在 1200 步后仍有 1.5e-4），
   永远收敛不了。改成「壁面绝热（`ttype` 默认 0）+ T 初值 0」后被动场**恒等于 0**
   （`dT_max = 0.000E+00`，从 it = 1 就是），速度单独收敛。
   **速度解不受影响**：`dT_max ≡ 0` 的那次运行与 T 漂移的那次运行逐行相同
   （对照用 α = 0 控制 `variants/bj_a0.control`：两次的 `it=1` 行只有 `dT_max` 不同
   —— `1.533E-01` vs `0.000E+00` —— 而 `du_max = 1.706E-04`、`lin-res = 1.01E-09`
   完全相同，`it=600 du_max = 1.229E-08` 也完全相同），
   第 4 节的 uF1 / RMS 也与归档值一致。
2. `outer_max` 800 → 2000：归档运行在 800 步上限处 du_max 还有 3.1e-9（判据是 1e-9），
   即「reached outer_max」收尾；2000 步让现场能打出显式 `CONVERGED` 行。

## 6. 备注

* `gen_bj.py` 是网格生成脚本（记录几何怎么来的），`plot_bj.py` 是校验/绘图脚本；
  两者都随现场一起带走，集群上只需要 python3 + numpy + matplotlib。
* `./run.sh 4` 不再只是确定性检查：非结构模式 `-np > 1` 默认走 METIS 分块并行
  （`uns_parallel = auto`）；加 `uns_parallel = never` 可退回「每个 rank 一份整网格」的
  冗余模式（结果与串行逐位相同，仅用于 MPI 确定性检查）。
* **核数声明与守卫**：`mix.control` 里取消注释 `n_struct_ranks = 0` + `n_uns_ranks = N`
  （两键必须成对），和 = N 必须等于发作业核数，否则求解器启动时打印诊断并**终止作业**
  （退出码 6 只写一个键 / 7 和不符 / 8 与模式冲突）；`./run.sh`（不带参数）会把这个和当作
  默认 np。不写这两个键则保持旧的隐式划分、不做检查。详见 `cluster_cases/README.md` 第 3 节。
* 想改 α 又不想动主控制：直接 `cp variants/bj_a2.control unMesh.control` 再跑零参数二进制。
