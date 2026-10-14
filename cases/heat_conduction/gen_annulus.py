#!/usr/bin/env python3
# ============================================================================
# gen_annulus.py -- 2D 圆环（极坐标）纯导热验证算例的六面体 .cas 网格生成器
#
# 几何：r in [RI, RO]，theta in [0, 2*pi)（闭合整环，theta 方向周期），
#       z in [0, LZ]（z 方向只 1 层胞）。
#       网格 Nr（径向）x Nt（周向）x 1（轴向），坐标单位 m。
#
# theta = 0 与 theta = 2*pi 处为**重合节点**（seam），两侧胞通过 it = Nt 的
# 内部面共享，从而形成拓扑自洽的周期环（不需要 special periodic 处理）。
#
# 面区：
#   3  Interior      内部面（r = const 面 + theta = const 面，含 seam 面）
#   4  symmetry (z)  z = 0 / z = LZ 两个面 -> 绝热
#   5  inner wall    r = RI（.control 里给 bc + tbc）
#   6  outer wall    r = RO（.control 里给 bc + tbc）
#
# 体区：
#   zone 2  整个圆环（单一固体区）
#
# 用法：  python3 gen_annulus.py <out.cas> [Nr] [Nt]
#         默认 RI=0.01 m, RO=0.03 m, LZ=0.002 m, Nr=20, Nt=72
#
# 所有 Fluent 索引按十六进制输出。
# ============================================================================
import math
import sys

RI, RO, LZ = 0.01, 0.03, 0.002     # m
NR, NT = 20, 72

argv = sys.argv[1:]
if len(argv) < 1:
    raise SystemExit('usage: gen_annulus.py <out.cas> [Nr] [Nt]')
out = argv[0]
if len(argv) >= 2:
    NR = int(argv[1])
if len(argv) >= 3:
    NT = int(argv[2])

hx = lambda n: '%x' % n

NODE = (NR + 1) * (NT + 1) * 2
NCELL = NR * NT


def nid(ir, k, iz):
    """节点编号：ir 最慢、iz 最快。k = 0..NT（k = NT 与 k = 0 位置重合）。"""
    return 1 + ir + k * (NR + 1) + iz * (NR + 1) * (NT + 1)


def cell(ir, it):
    """单元编号：it 最快（单区，编号连续性只需整体连续）。"""
    return 1 + it + ir * NT


def r_of(ir):
    return RI + (RO - RI) * ir / NR


def th_of(k):
    """周向角；k = 0..NT，k = NT 与 k = 0 位置重合（seam）。"""
    return 2.0 * math.pi * k / NT


# ---- 面数统计 --------------------------------------------------------------
NF_RINT = (NR - 1) * NT            # r = const 内部面
NF_TINT = NR * NT                  # theta = const 内部面（含 seam）
NF_INT = NF_RINT + NF_TINT
NF_ZSYM = 2 * NR * NT              # z = 0 / z = LZ（zone 4）
NF_MIN = NT                        # r = RI（zone 5）
NF_MAX = NT                        # r = RO（zone 6）
NF = NF_INT + NF_ZSYM + NF_MIN + NF_MAX

L = []
L.append('(0 "MixNSSolver 2D annulus solid-conduction case, gen_annulus.py")')
L.append('')
L.append('(2 3)')                                  # 3D
L.append('(10 (0 1 %s 0 3))' % hx(NODE))
L.append('(12 (0 1 %s 0))' % hx(NCELL))
L.append('(13 (0 1 %s 0))' % hx(NF))
L.append('')

# ---- 节点 ------------------------------------------------------------------
L.append('(0 "Zone 1  Number of Nodes : %d")' % NODE)
L.append('(10 (1 1 %s 1 3)(' % hx(NODE))
for iz in range(2):
    for k in range(NT + 1):
        tt = th_of(k)
        for ir in range(NR + 1):
            r = r_of(ir)
            L.append('  %.15e   %.15e   %.15e'
                     % (r * math.cos(tt), r * math.sin(tt), iz * LZ))
L.append('))')
L.append('')

# ---- 体区 ------------------------------------------------------------------
L.append('(0 "Zone 2 %d cells 1..%d, solid annulus")' % (NCELL, NCELL))
L.append('(12 (2 1 %s 1 4))' % hx(NCELL))
L.append('(45 (2 solid unspecified)())')
L.append('')

# ---- 面 --------------------------------------------------------------------
face_no = 0

L.append('(0 "Zone 3 %d faces %d..%d, Interior")'
         % (NF_INT, face_no + 1, face_no + NF_INT))
L.append('(13 (3 %s %s 2 4)(' % (hx(face_no + 1), hx(face_no + NF_INT)))
for ir in range(1, NR):                            # r = const 内部面
    for it in range(NT):
        L.append(' '.join(hx(v) for v in
                          (nid(ir, it, 0), nid(ir, it + 1, 0),
                           nid(ir, it + 1, 1), nid(ir, it, 1)))
                 + ' ' + hx(cell(ir - 1, it)) + ' ' + hx(cell(ir, it)))
face_no += NF_RINT
for ir in range(NR):                               # theta = const 内部面
    for it in range(1, NT + 1):
        it0 = it - 1
        it1 = 0 if it == NT else it                # seam: it = NT 接回 it = 0
        L.append(' '.join(hx(v) for v in
                          (nid(ir, it, 0), nid(ir + 1, it, 0),
                           nid(ir + 1, it, 1), nid(ir, it, 1)))
                 + ' ' + hx(cell(ir, it0)) + ' ' + hx(cell(ir, it1)))
face_no += NF_TINT
L.append('))')
L.append('')

L.append('(0 "Zone 4 %d faces %d..%d, Symmetry (z)")'
         % (NF_ZSYM, face_no + 1, face_no + NF_ZSYM))
L.append('(13 (4 %s %s 7 4)(' % (hx(face_no + 1), hx(face_no + NF_ZSYM)))
for iz in range(2):
    for ir in range(NR):
        for it in range(NT):
            L.append(' '.join(hx(v) for v in
                              (nid(ir, it, iz), nid(ir, it + 1, iz),
                               nid(ir + 1, it + 1, iz), nid(ir + 1, it, iz)))
                     + ' ' + hx(cell(ir, it)) + ' 0')
face_no += NF_ZSYM
L.append('))')
L.append('')

L.append('(0 "Zone 5 %d faces %d..%d, inner wall r = RI")'
         % (NF_MIN, face_no + 1, face_no + NF_MIN))
L.append('(13 (5 %s %s 3 4)(' % (hx(face_no + 1), hx(face_no + NF_MIN)))
for it in range(NT):
    L.append(' '.join(hx(v) for v in
                      (nid(0, it, 0), nid(0, it + 1, 0),
                       nid(0, it + 1, 1), nid(0, it, 1)))
             + ' ' + hx(cell(0, it)) + ' 0')
face_no += NF_MIN
L.append('))')
L.append('')

L.append('(0 "Zone 6 %d faces %d..%d, outer wall r = RO")'
         % (NF_MAX, face_no + 1, face_no + NF_MAX))
L.append('(13 (6 %s %s 3 4)(' % (hx(face_no + 1), hx(face_no + NF_MAX)))
for it in range(NT):
    L.append(' '.join(hx(v) for v in
                      (nid(NR, it, 0), nid(NR, it + 1, 0),
                       nid(NR, it + 1, 1), nid(NR, it, 1)))
             + ' ' + hx(cell(NR - 1, it)) + ' 0')
face_no += NF_MAX
L.append('))')
L.append('')

assert face_no == NF, (face_no, NF)

with open(out, 'w') as f:
    f.write('\n'.join(L) + '\n')

print('wrote %s : Nr=%d Nt=%d -> %d cells, %d faces'
      % (out, NR, NT, NCELL, NF))
print('  annulus RI=%.6g RO=%.6g LZ=%.6g m, dr=%.6g m, dtheta=%.6g rad (%.4g deg)'
      % (RI, RO, LZ, (RO - RI) / NR, 2 * math.pi / NT, 360.0 / NT))
print('  zone 5 = r = RI (inner wall), zone 6 = r = RO (outer wall)')

def th_of(k):
    return 2.0 * math.pi * k / NT
