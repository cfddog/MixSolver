#!/usr/bin/env python3
# ============================================================================
# check_heat_conduction.py -- cases/heat_conduction 5 个纯导热算例 vs 解析解
#
# 读取 .dbg/heat_conduction/<case>/unMesh.vtu（由 run_all.sh 跑出），
# 取 CellData 的 "temperature_solid"（固体导热温度场）与解析解对拍。
#
# .vtu 读取约定（与 cases/porous_plug/plot_plug.py 一致）：
#   * vtk_write 把每个六面体胞拆成 12 个四面体（"face fan + 胞心"）；
#   * connectivity 每行就是该四面体的节点号 "id0 id1 id2 cen"（VTK 规范写法，
#     无前导顶点数 token），offsets 为累积值 4,8,12,...；2026-10-15 之前写的是
#     legacy VTK CELLS 风格 "4 id0 id1 id2 cen"（offsets 5,10,15,...），
#     本脚本两种都能读（见 load_vtu 的 token 兼容处理）；
#   * points 末尾追加 m%ncells 个胞心坐标 g%xc，其 0-based id = nnodes+c-1，
#     即每个四面体的**最大节点号**就是它所属原胞的胞心 id；
#   * CellData 每个数组按原胞顺序写 nsub=12 个相同值（块状）。
#   本脚本用 offsets 切块读单元节点号，再用
#   "块内恰好 4 个节点号 / 胞心 id = 块编号" 双重自校验。
#
# 用法：
#   python3 check_heat_conduction.py               # 检查全部已运行算例
#   python3 check_heat_conduction.py c1 c5         # 只检查指定算例
#   python3 check_heat_conduction.py --tables      # 附打印沿 x / r 剖面表
#   python3 check_heat_conduction.py -v c1         # 附打印前 4 个胞明细
#   python3 check_heat_conduction.py --refine      # 追加 c3f（径向加密）并报收敛比
#   python3 check_heat_conduction.py --selftest    # c5 回归守卫自测（坏资产必须报错）
#
# 退出码：0 = 全部通过；1 = 有算例超差 / 缺 .vtu / 缺 CellData / 守卫失效。
# ============================================================================
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DBG = os.path.join(ROOT, '.dbg', 'heat_conduction')

# ---- 与 gen_1d.py / gen_annulus.py / *.control 保持一致的几何与物性 ---------
L_1D = 0.1                 # m，1D 棒长
T_COLD, T_HOT = 300.0, 500.0
K_S = 15.0                 # W/(m*K)
Q_SRC = 1.0e6              # W/m^3，c5 热源区
X_A, X_B = 0.03, 0.07      # m，c5 热源区 [X_A, X_B)
RI, RO = 0.01, 0.03        # m，圆环内/外半径


def ref_1d_linear(x, y):
    """c1：两端定温 300 / 500 K -> 线性剖面。"""
    return T_COLD + (T_HOT - T_COLD) * x / L_1D


def ref_400(x, y):
    """c2 / c4：一处定温 400 K + 其余绝热 -> 处处 400 K。"""
    return np.full_like(x, 400.0)


def ref_annulus(x, y):
    """c3：圆环内外壁定温 -> 对数剖面 T(r)。"""
    r = np.hypot(x, y)
    return T_COLD + (T_HOT - T_COLD) * np.log(r / RI) / np.log(RO / RI)


def ref_1d_source(x, y):
    """c5：三段解析解。

    q/k = Q_SRC/K_S；a = X_A，b = X_B。左/中/右段斜率 Al/Am/Ar 与截距
    Bl/Bm/Br 由 (i) 分段斜率连续性 Al = Am - (q/k)a、Ar = Am - (q/k)b，
    (ii) 斜率连续处的温度连续性，(iii) 两端定温 T(0)=T_COLD、
    T(L)=T_HOT 唯一确定，闭式结果如下（与 c5_1d_source.control 注释一致）：
        Bl = T_COLD
        Bm = Bl - (q/2k) a^2
        Br = Bm + (q/2k) b^2
        Ar = (T_HOT - Br) / L
        Am = Ar + (q/k) b
        Al = Am - (q/k) a
    """
    qk = Q_SRC / K_S
    a, b = X_A, X_B
    Bl = T_COLD
    Bm = Bl - 0.5 * qk * a * a
    Br = Bm + 0.5 * qk * b * b
    Ar = (T_HOT - Br) / L_1D
    Am = Ar + qk * b
    Al = Am - qk * a
    left = Al * x + Bl
    mid = -0.5 * qk * x * x + Am * x + Bm
    right = Ar * x + Br
    return np.where(x < a, left, np.where(x < b, mid, right))

def note_c5(ctr, Ts, ref, err, cd):
    """c5 专用报告 ＋ 回归守卫。

    ① 回归守卫（对应 .trae/rules/project_rules.md「多 cell zone 网格的 cell-id 排序」）：
       q_src>0 的胞必须**恰好**是 x in [X_A, X_B) 那些胞。gen_1d.py 早期版本在 NY=2 下
       按 `nz2 = ix0*NY` 切体区，结果 zone 8 的 id 段实际落在 j=0 行的 x in [0.06, 0.1)
       （体区记录只能是单段连续 id），q_src 漏到区间外，最大误差 26.6 K。此守卫能让这种
       资产立刻报错，而不是静默给出"看起来合理"的温度场。
    ② 离散化说明：源区边界正好落在网格面上，解析解在 x=0.03/0.07 处二阶导跳变
       （"拐点"），格心 FV 格式在这两个界面胞上留下截断误差，其闭式解为

           delta = q * h^2 / (8 k)      （表现为源区内部的一个常数偏移）

       本算例 q=1e6 W/m^3, k=15 W/(m*K), h=0.0025 m -> delta = 0.0520833 K。
       实测偏差与之相符 ⇒ 剩余偏差是**二阶截断误差**而非实现错误（同一离散系统的独立
       numpy 直接求解与求解器结果只差 ~3e-7 K）。
    """
    h, qk = 0.0025, Q_SRC / K_S
    qs = cd.get('heat_source_solid')
    if qs is not None:
        on = np.flatnonzero(qs > 0.0)
        allex = np.unique(np.round(ctr[:, 0], 9))
        # 体区 zone 8 的胞心应恰好是落在 (X_A, X_B) 内的那些 x 面
        exp_x = allex[(allex > X_A + 1e-9) & (allex < X_B - 1e-9)]
        got_x = np.unique(np.round(ctr[on, 0], 9))
        print('      回归守卫 : q_src>0 的胞 %d 个，覆盖 %d 个 x 面（应 %d），'
              'x in [%.6g, %.6g]，应为 x in [%.6g, %.6g]'
              % (len(on), len(got_x), len(exp_x),
                 got_x.min() if len(got_x) else float('nan'),
                 got_x.max() if len(got_x) else float('nan'),
                 exp_x.min(), exp_x.max()))
        if not np.array_equal(got_x, exp_x):
            raise AssertionError(
                'q_src 落点与体区 id 段不匹配：q>0 覆盖 %d 个 x 面 [%.6g, %.6g]，'
                '应为 (%.6g, %.6g) 内的 %d 个面；见 .trae/rules/project_rules.md'
                '「多 cell zone 网格的 cell-id 排序」（分割方向索引必须最慢 ⇒ '
                '--src 必须 --ny 1）'
                % (len(got_x),
                   got_x.min() if len(got_x) else float('nan'),
                   got_x.max() if len(got_x) else float('nan'),
                   X_A, X_B, len(exp_x)))
    d = qk * h * h / 8.0
    insrc = (ctr[:, 0] >= X_A) & (ctr[:, 0] < X_B)
    print('      理论"拐点"截断误差 q*h^2/(8k) = %+.8f K' % d)
    print('      实测偏差：源区 mean = %+.8f K（≈上式），区外 max|误差| = %.2e K'
          % (err[insrc].mean(), np.abs(err[~insrc]).max()))


def guard_selftest():
    """回归守卫自测（.trae/rules/project_rules.md「验证脚本必须自带回归守卫」）。

    在**历史缺陷资产**上跑守卫，必须报错：旧 `gen_1d.py --src`（未强制 --ny 1）
    把 zone 8 的 id 段错当成 x in [0.06, 0.1) 的 j=0 行胞（16 个），q_src 漏到
    区间外、最大误差 26.6 K。这里用当前 c5 的网格 + 人工把 q_src 挪到那 16 个胞，
    复现坏资产；守卫必须抛 AssertionError。再在正确资产上确认守卫不误报。

    用法 : python3 check_heat_conduction.py --selftest
    """
    vtu = os.path.join(DBG, 'c5', 'unMesh.vtu')
    if not os.path.isfile(vtu):
        print('  SKIP  %s 不存在（先 bash run_all.sh c5）' % vtu)
        return 1
    ncell, nsub, ctr, cd = load_vtu(vtu)
    Ts = cd['temperature_solid']
    ref = ref_1d_source(ctr[:, 0], ctr[:, 1])
    err = Ts - ref
    print('  [1/2] 正确资产（q_src 在 x in [%.2f, %.2f)）：守卫应通过' % (X_A, X_B))
    note_c5(ctr, Ts, ref, err, cd)

    print('  [2/2] 坏资产（复现旧 --src 未强制 --ny 1 的落点 x in [0.06, 0.1) j=0 行）：'
          '守卫应报错')
    bad = dict(cd)
    q = np.zeros(ncell)
    on = np.flatnonzero((ctr[:, 0] >= 0.06) & (ctr[:, 0] < L_1D))
    q[on] = Q_SRC
    bad['heat_source_solid'] = q
    try:
        note_c5(ctr, Ts, ref, err, bad)
    except AssertionError as e:
        print('  OK  : 守卫按预期报错 -> %s' % e)
        return 0
    print('  FAIL: 守卫**没有**在坏资产上报错（回归守卫失效）')
    return 1


# 算例表：每个条目给出控制文件、剖面坐标轴、解析解函数、判定容差（K）
CASES = {
    'c1': dict(ctl='c1_1d_fixedT.control', axis='x', ref=ref_1d_linear,
               tol=1.0e-6,
               desc='1D，两端定温 300/500 K，侧面绝热（解析：线性）'),
    'c2': dict(ctl='c2_1d_adiabatic.control', axis='x', ref=ref_400,
               tol=1.0e-6,
               desc='1D，x=0 定温 400 K + 其余绝热（解析：均匀 400 K）'),
    'c3': dict(ctl='c3_annulus_fixedT.control', axis='r', ref=ref_annulus,
               tol=1.0,
               desc='2D 圆环，内/外壁定温 300/500 K，z 绝热（解析：对数）'),
    'c4': dict(ctl='c4_annulus_adiabatic.control', axis='r', ref=ref_400,
               tol=1.0e-6,
               desc='2D 圆环，内壁定温 400 K + 外壁绝热（解析：均匀 400 K）'),
    'c5': dict(ctl='c5_1d_source.control', axis='x', ref=ref_1d_source,
               tol=0.1, extra=note_c5,
               desc='1D，两端定温 300/500 K + 中段 q=1e6 W/m^3（解析：分段）'),
}

# 附加（非必须）算例：c3 的径向加密网格，用于网格收敛性检查。
# 网格在 .dbg/heat_conduction/c3f/ 里现场生成（不入库）：
#   bash run_all.sh c3f        # Nr=40, Nt=72（C3F_NR= 可改 Nr）
# 控制文件与 c3 完全相同（zone id 一致）。
REFINE = {
    'c3f': dict(ctl='c3_annulus_fixedT.control', axis='r', ref=ref_annulus,
                tol=1.0,
                desc='c3 的径向加密网格 Nr=40（网格收敛性检查，附加算例）'),
}
ALL = dict(CASES)
ALL.update(REFINE)




def _array(name, text):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii">(.*?)</DataArray>'
                  % re.escape(name), text, re.S)
    if m is None:
        return None
    return np.fromstring(m.group(1), sep=' ')


def load_vtu(fn):
    """返回 (ncell, nsub, ctr[ncell,3], celldata{name: vec[ncell]})。"""
    s = open(fn).read()
    pts = _array('Coordinates', re.search(r'<Points>(.*?)</Points>',
                                          s, re.S).group(1)).reshape(-1, 3)
    cells = re.search(r'<Cells>(.*?)</Cells>', s, re.S).group(1)
    conn = _array('connectivity', cells).astype(int)
    off = _array('offsets', cells).astype(int)
    cdat = re.search(r'<CellData>(.*?)</CellData>', s, re.S)
    assert cdat is not None, 'no <CellData> in ' + fn
    cdat = cdat.group(1)

    # 用 offsets 切块 -> 单元节点号。兼容两种写法：
    #   * 规范（2026-10-15 起）：块内就是 4 个节点号；
    #   * 旧文件：块内 5 个 token "4 n0 n1 n2 cen"（首 token 是顶点数）。
    st = np.concatenate(([0], off[:-1]))
    ntets = len(off)
    cen_id = np.empty(ntets, dtype=int)
    for i in range(ntets):
        ids = conn[st[i]:off[i]]
        if ids[0] == 4 and len(ids) == 5:   # 旧文件：legacy "4 n0 n1 n2 cen"
            ids = ids[1:]                   # （不能用 len(ids)==ids[0]+1 判：
                                            #  新格式里首个节点号恰为 3 时会误判！）
        assert len(ids) == 4, 'bad CELLS block at tet %d: %s' % (i, ids)
        cen_id[i] = ids.max()          # 胞心点 id = 该四面体最大节点号
    nnodes0 = int(cen_id.min())        # 0-based 首个胞心 id = 网格节点数
    ncell = len(pts) - nnodes0
    assert ncell > 0 and ntets % ncell == 0, (ntets, ncell)
    nsub = ntets // ncell
    # 自校验：前 nsub 个四面体属于胞 0，依次类推
    assert np.array_equal(cen_id - nnodes0, np.repeat(np.arange(ncell), nsub)), \
        'tet -> cell mapping is not block-ordered'
    ctr = pts[nnodes0:]

    celldata = {}
    for nm in ('pressure', 'velocity', 'temperature', 'temperature_solid',
               'heat_source_solid'):
        a = _array(nm, cdat)
        if a is None:
            continue
        if nm == 'velocity':
            a = a.reshape(-1, 3)[:, 0]      # 只关心 x 分量（全为 0）
        assert len(a) == ntets, (nm, len(a), ntets)
        b = a.reshape(ncell, nsub)
        assert np.allclose(b, b[:, :1], rtol=0.0, atol=1e-12 * max(1.0, np.abs(b).max())), \
            'CellData "%s" is not replicated per cell' % nm
        celldata[nm] = b[:, 0]
    return ncell, nsub, ctr, celldata


def log_summary(case):
    """从 run.log 里提取最后一步的外迭代号与 CONVERGED 行（尽力而为）。"""
    path = os.path.join(DBG, case, 'run.log')
    if not os.path.isfile(path):
        return None
    it, conv = None, None
    for ln in open(path, 'r', errors='replace'):
        m = re.match(r'\s*it=\s*(\d+)', ln)
        if m:
            it = int(m.group(1))
        if 'CONVERGED' in ln or 'TIMEOUT' in ln:
            conv = ln.strip()
    return it, conv



def main():
    flags = [a for a in sys.argv[1:] if a.startswith('-')]
    named = [a for a in sys.argv[1:] if not a.startswith('-')]
    if named and named[0] == '--selftest':
        named = named[1:]
    cases = named or sorted(CASES)
    verbose, tables = '-v' in flags, '--tables' in flags
    if '--refine' in flags and not named:
        cases = sorted(CASES) + ['c3f']

    print('=' * 76)
    print('cases/heat_conduction -- 纯导热算例 vs 解析解')
    print('取量 : <case>/unMesh.vtu 的 CellData "temperature_solid"（固体温度）')
    print('结果目录 : %s' % DBG)
    print('=' * 76)

    if '--selftest' in flags:
        print('回归守卫自测（--selftest）')
        return guard_selftest()

    nfail = 0
    emax_of = {}
    for c in cases:
        if c not in ALL:
            print('[%s] unknown case (known: %s)' % (c, ' '.join(sorted(ALL))))
            nfail += 1
            continue
        info = ALL[c]
        vtu = os.path.join(DBG, c, 'unMesh.vtu')
        print('')
        print('-' * 76)
        print('[%s] %s' % (c, info['desc']))
        print('      控制文件 : %s' % info['ctl'])

        lg = log_summary(c)
        if lg:
            print('      run.log  : 最后 it=%s , %s' % (lg[0], lg[1] or '(无 CONVERGED 行)'))

        if not os.path.isfile(vtu):
            print('  MISSING %s' % vtu)
            print('  -> 先运行 : bash %s/run_all.sh %s' % (HERE, c))
            nfail += 1
            continue

        ncell, nsub, ctr, cd = load_vtu(vtu)
        assert 'temperature_solid' in cd, 'no temperature_solid in ' + vtu
        Ts = cd['temperature_solid']
        ref = info['ref'](ctr[:, 0], ctr[:, 1])
        err = Ts - ref
        ax = np.hypot(ctr[:, 0], ctr[:, 1]) if info['axis'] == 'r' else ctr[:, 0]
        i = int(np.argmax(np.abs(err)))

        print('      胞数 %d（每个六面体胞在 .vtu 里拆成 %d 个四面体）' % (ncell, nsub))
        print('      T_solid      : %.6f .. %.6f K' % (Ts.min(), Ts.max()))
        print('      解析解       : %.6f .. %.6f K' % (ref.min(), ref.max()))
        print('      最大绝对误差 : %.3e K  (%.4f %% of 200 K)  @ %s = %.6g m'
              % (abs(err[i]), abs(err[i]) / (T_HOT - T_COLD) * 100.0,
                 info['axis'], ax[i]))
        print('      RMS 误差     : %.3e K' % math.sqrt(float((err ** 2).mean())))

        if 'temperature' in cd:
            dlt = np.abs(cd['temperature'] - Ts).max()
            print('      |T_lte - T_solid|_max = %.3e K' % dlt)
            print('        （固体胞里 T_lte 被钉成"上一轮"的 T_s；本类纯固体算例第 1 次')
            print('          外迭代即满足收敛判据，故 T_lte 落后一轮，物理量请看 '
                  'temperature_solid）')
        if 'velocity' in cd:
            print('      max |u_x| = %.3e m/s（纯固体：应为 0）'
                  % np.abs(cd['velocity']).max())
        if 'heat_source_solid' in cd:
            qs = cd['heat_source_solid']
            on = np.flatnonzero(qs > 0.0)
            if len(on):
                print('      heat_source_solid : %d 个胞 > 0，q = %.6g W/m^3，'
                      'x in [%.6g, %.6g] m'
                      % (len(on), qs[on[0]], ctr[on, 0].min(), ctr[on, 0].max()))
            else:
                print('      heat_source_solid : 全为 0')

        if info.get('extra') is not None:
            info['extra'](ctr, Ts, ref, err, cd)

        ok = abs(err).max() <= info['tol']
        emax_of[c] = float(abs(err).max())
        print('      判定 : %s   (tol = %.1e K)' % ('PASS' if ok else 'FAIL', info['tol']))
        if not ok:
            nfail += 1

        if tables:
            order = np.argsort(ax)
            step = max(1, ncell // 10)
            print('      %-14s %-14s %-14s %-12s'
                  % (info['axis'] + ' [m]', 'T_solid [K]', 'analytic [K]', 'diff [K]'))
            for k in order[::step]:
                print('      %-14.6g %-14.6f %-14.6f %+.3e'
                      % (ax[k], Ts[k], ref[k], err[k]))
        if verbose:
            for k in range(min(4, ncell)):
                print('      cell %d : T_s=%.9f  ref=%.9f  (x=%.6g, y=%.6g)'
                      % (k, Ts[k], ref[k], ctr[k, 0], ctr[k, 1]))

    if 'c3' in emax_of and 'c3f' in emax_of:
        e1, e2 = emax_of['c3'], emax_of['c3f']
        print('')
        print('-' * 76)
        print('网格收敛性（c3 圆环，径向加密 Nr: 20 -> 40，周向 Nt=72 不变）')
        print('  最大绝对误差 : %.3e K  ->  %.3e K   比值 %.2f'
              % (e1, e2, (e1 / e2) if e2 > 0 else float('inf')))
        print('  （C3F_NR=20/40/80 可反复运行 run_all.sh c3f 做更细的序列；误差来源：'
              '内/外壁')
        print('    Dirichlet 的“半格心”壁面距离对圆柱几何为 1 阶，环内径向扩散为 2 阶，'
              '故比值介于 1~2。）')

    print('')
    print('=' * 76)
    print('结果 : %s' % ('全部 PASS' if nfail == 0 else '%d 项未通过' % nfail))
    print('=' * 76)
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
