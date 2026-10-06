# porous-plug：一维多孔塞 Darcy / Darcy-Forchheimer 压降验证

阶段 9 首个验证算例。验证两件事：

1. **VC 标记自动识别**：CAS cell zone 名称带 `VC:porous` 标记时，求解器自动把该块判为
   多孔介质，`.control` 只需给系数（无 `fluid|porous` 类型词，`CZ_AUTO`）；
2. **多孔动量源项定量正确**：一维均匀多孔塞的稳态压降对拍 Darcy / Darcy-Forchheimer
   解析解。

2026-10-05 追加：**各向异性渗透率（轴对齐对角张量）验证**，见第 4 节。

## 1. 算例设置

- 网格：`gen_plug.py` 生成结构化 hex 管道，40×2×1，L=0.1 m，H=0.01 m，
  z 厚 0.002 m（`plug.cas`，80 cells，坐标已是米，`mesh_scale=1`）。
- 四个侧面全部 symmetry（无壁面边界层），x=0 velocity-inlet，x=L pressure-outlet(0)，
  因此流动严格一维。
- 流体：300 K 空气，ρ=1.177 kg/m³，μ=1.846e-5 Pa·s。
- 多孔介质：孔隙率 ε=0.4，渗透率 K=1e-8 m²（Brinkman 有效粘度 μ_eff=μ/ε）。
- cell zone 的 (45 记录写为 `(45 (2 VC:porous plug)())`（CAS 分词器不处理引号，
  标记必须无空格；`VC:` 后允许冒号/空格分隔，大小写不敏感）。
- 控制文件系数行不含类型词：
  `cell_zone = 2 perm=1.0e-8 porosity=0.4`（Forch 工况另加 `inertial=500.0`）。

解析解（稳态一维，μ_eff=μ/ε）：

```
dp/dx = -(μ/ε/K) u                    [Darcy]
dp/dx = -(μ/ε/K) u - ρ C_F |u| u      [Darcy-Forchheimer]
```

| 工况 | u_in (m/s) | C_F (1/m) | Re_p=ρu√K/μ | 理论 Δp (Pa) |
|---|---|---|---|---|
| plug_darcy | 0.1 | 0 | 0.64 | 46.150 |
| plug_forch | 1.0 | 500 | 6.4 | 520.350（461.5 Darcy + 58.85 Forch） |

## 2. 复现命令

```bash
python3 gen_plug.py plug.cas
../../bin/uns_solver plug.cas plug_darcy.control > plug_darcy.log 2>&1
mv plug.vtu plug_darcy.vtu
../../bin/uns_solver plug.cas plug_forch.control > plug_forch.log 2>&1
mv plug.vtu plug_forch.vtu
python3 plot_plug.py
# MPI:
mpirun -np 2 ../../bin/uns_solver_mpi plug.cas plug_darcy.control
```

## 3. 结果

启动日志确认类型来自 VC 标记（无需显式声明）：

```
--- Cell (volume) zones ---
    id  condition          name                ncells  type    source
     2           VC:porous              plug       80  porous  VC tag
```

内部列（col 3..36，避开入口/出口边界列）线性拟合压力梯度：

| 工况 | 拟合 dp/dx (Pa/m) | 理论 (Pa/m) | 误差 | 外推 Δp (Pa) | Δp 误差 | 内部 u 最大偏差 |
|---|---|---|---|---|---|---|
| Darcy | −461.500 | −461.500 | **0.000%** | 46.150 | **0.000%** | 0.007% |
| Darcy-Forch | −5203.492 | −5203.500 | **0.002%** | 520.349 | **0.002%** | 0.016% |

- 压力沿 x 严格线性，斜率与解析解一致到机器精度量级；
- 内部速度均匀（=u_in，偏差 <0.02%）；首列/末列的速度外观偏差
  （SIMPLE 入口列固有欠修正 + pressure-outlet 边界列不动点，
  与 cases/channel、couple_channel 中记录的现象一致）不影响内部压降；
- Darcy 工况 61 外迭代收敛（du_max=1.8e-11，mass-imbal≈1e-18），
  Forch 工况 57 步；np=2 与串行同为 61 步且 du_max 位级一致。

图：`images/plug_pressure.png`（两工况 p(x) 对解析直线 + u/u_in 剖面）。

## 4. 结论

- VC 标记自动识别 + `CZ_AUTO` 系数行路径正确；显式 `cell_zone ... porous|fluid`
  优先级高于 VC 标记（回归 porous_Ra10 报告 source=control，Nu=1.078 与归档一致）。
- 多孔介质 Darcy 线性项与 Forchheimer 二次项的 1D 压降定量正确。

## 4. 各向异性渗透率验证（2026-10-05 追加）

`cell_zone` 新增可选对角张量分量 `perm_xx=` `perm_yy=` `perm_zz=`（轴对齐各向异性；
未给分量回退标量 `perm=`，交叉项需要块耦合动量装配、不支持）。Darcy 汇按动量分量
取对应对角元：`ap += (μ/ε)/K_d·V`。

| 工况 | 设置 | 拟合 dp/dx (Pa/m) | 理论 (Pa/m) | 误差 |
|---|---|---|---|---|
| plug_aniso1(b) | 仅 `perm_xx=1.0e-8`（yy/zz 回退 perm=0） | −461.109 | −461.500 | **0.085%** |
| plug_aniso2 | `perm=1.0e-8 perm_xx=5.0e-9` | −921.876 | −923.000（=2×Darcy） | **0.122%** |
| plug_darcy（标量路径对照） | `perm=1.0e-8` | −461.500 | −461.500 | 0.000% |

- 张量路径与标量路径一致；Kxx 减半 → 压降精确加倍（对角分量逐方向选取正确）；
- aniso1 横向无 Darcy 阻力，入口列恢复段变长（col0=0.0175、col5≈0.1，仍收敛到
  正确固定点 u=0.09997@col20），拟合窗口须避开前 10 列（`outer_tol=1e-9`）；
- 热弥散（disp_l/disp_t）的验证另见 cases/porous_disp/。

新增文件：`plug_aniso1.control`（及高容差版 `plug_aniso1b`）、`plug_aniso2.control`、
对应 log/vtu。

## 5. 已修复的入口温度 bug（影响所有独立 velocity-inlet 算例）

`bc = <zone> velocity-inlet ux uy uz [T]` 的第 4 个 token（静温）此前**不被解析**，
`bc_face_T` 将 BC_VINLET 的 Dirichlet 温度锚定到缺省 tval=0，独立运行时全场温度
被排空到 0（耦合运行不受影响，界面温度走 set_interface_T）。现已解析可选温度
token（缺省仍 0）。本目录 control 均写明 T=300；刷新后 darcy/forch 收敛加快
（12/33 步，`outer_tol=1e-9`），压降不变。

## 6. 文件清单

- `gen_plug.py`：网格生成器（cell zone 名称内嵌 VC:porous）
- `plug.cas`：网格
- `plug_darcy.control` / `plug_forch.control`：两工况控制文件
- `plug_darcy.log` / `plug_forch.log`：运行日志
- `plug_darcy.vtu` / `plug_forch.vtu`：收敛流场
- `plot_plug.py`、`images/plug_pressure.png`：对比图
