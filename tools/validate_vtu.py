#!/usr/bin/env python3
# ============================================================================
# validate_vtu.py -- MixNSSolver 写出的 VTK XML (.vtu) 格式守卫
#
# 目的：把「.vtu 必须严格符合 VTK XML 规范」变成可自动执行的检查，防止
# 2026-10-15 修掉的老毛病回归 —— 当时 <connectivity> 的每个单元前面多写了一个
# legacy VTK "CELLS" 顶点数 token（四面体写成 "4 a b c cen"、offsets 为
# 5,10,15,...）。ParaView 靠 offsets-顶点数 反推容错，所以看不出问题；Tecplot
# 360 (2020) 严格解析，把这个 token 当节点号 -> 网格错乱/打不开。
#
# 检查项（对每个 <Piece>，inline ascii 数据）：
#   1. 根元素 <VTKFile type="UnstructuredGrid">；
#   2. <Points> 的 3 分量数组长度 == 3*NumberOfPoints；
#   3. <connectivity> 只含节点号：长度 == offsets[-1]，且每个 id ∈ [0,npts)；
#   4. <offsets> 严格递增，且每单元节点数 offsets[i]-offsets[i-1] 必须等于该
#      单元 types[i] 的规范节点数（四面体=4、六面体=8、...）；
#      **若等于规范节点数 + 1，直接判定为 legacy 前导 token 回归**；
#   5. <types> 长度 == NumberOfCells，且类型码已知；
#   6. <PointData>/<CellData> 每个数组长度 == NumberOfComponents*(npts|ncells)；
#      标签**带属性**时给出告警：仓库内 cases/**/*.py 用严正则 r'<CellData>(.*?)'
#      取数据块，mod_vtu_writer 故意写无属性裸标签，改回带属性会让它们全部失效；
#   7. --vtk：额外用 vtk 模块（与 ParaView 同源的读取器）复查单元类型/节点数，
#      未安装 vtk 时给出提示并跳过。
#
# 用法：
#   python3 tools/validate_vtu.py FILE.vtu [FILE2.vtu ...]   # 校验，退出码 0/1
#   python3 tools/validate_vtu.py -q FILE.vtu                # 只打印错误
#   python3 tools/validate_vtu.py --vtk FILE.vtu             # 附 vtk 读取器复查
#   python3 tools/validate_vtu.py --selftest                 # 守卫自测（必须 0）
#
# 退出码：0 = 全部通过；1 = 有文件不合规 / 自测失败；2 = 用法错误。
# ============================================================================
import os
import re
import sys
import xml.etree.ElementTree as ET

# ---- VTK XML 规范：固定节点数的单元类型码 --------------------------------
CELL_NNODE = {1: 1, 3: 2, 5: 3, 9: 4, 10: 4, 12: 8, 13: 6, 14: 5}
CELL_NAME = {1: 'VERTEX', 3: 'LINE', 5: 'TRIANGLE', 9: 'QUAD',
             10: 'TETRA', 12: 'HEXAHEDRON', 13: 'WEDGE', 14: 'PYRAMID'}


def _nums(text):
    """DataArray 文本 -> float 列表（空白/换行/逗号都能切分）。"""
    return [float(t) for t in re.split(r'[\s,]+', text or '') if t]


def _da(elem, name=None):
    """在 elem 下找第一个 <DataArray>（可选按 Name 过滤）；找不到返回 None。"""
    if elem is None:
        return None
    for da in elem.findall('DataArray'):
        if name is None or (da.get('Name') or '') == name:
            return da
    return None


def _int_attr(elem, attr, tag, err):
    v = elem.get(attr)
    try:
        return int(v)
    except (TypeError, ValueError):
        err.append('%s: 缺少/非法属性 %s="%s"' % (tag, attr, v))
        return None


def _da_ints(da, tag, err, warn):
    """取整数型 DataArray 的值；缺失/非 ascii 时记录并返回 None。"""
    if da is None:
        err.append('%s: 缺少 DataArray' % tag)
        return None
    fmt = da.get('format') or 'ascii'
    if fmt != 'ascii':
        warn.append('%s: format="%s"，本脚本只解析 inline ascii，跳过长度校验'
                    % (tag, fmt))
        return None
    return [int(v) for v in _nums(da.text)]


def _da_len(da, tag, err, warn):
    """取 DataArray 的元素个数（= ncomp*nent），非 ascii 时返回 None。"""
    if da is None:
        err.append('%s: 缺少 DataArray' % tag)
        return None
    fmt = da.get('format') or 'ascii'
    if fmt != 'ascii':
        warn.append('%s: format="%s"，本脚本只解析 inline ascii，跳过长度校验'
                    % (tag, fmt))
        return None
    return len(_nums(da.text))


def _check_data_section(piece, sec, nent, tag, err, warn):
    """PointData / CellData：每个数组长度必须 = ncomp*nent。"""
    e = piece.find(sec)
    if e is None:
        return
    # 仓库内 cases/**/*.py 的读取器用严格的 r'<CellData>(.*?)</CellData>' 取数据块
    # （历史上 MixNSSolver 写的就是无属性裸标签），所以一旦给标签加上属性，
    # 那 20 多个验证脚本会全部报 "no <CellData>"。这里给告警而不是报错：
    # 属性本身合规范，但会破坏本仓库的读取约定。
    if e.attrib:
        warn.append('%s: <%s> 带属性 %s —— 仓库内 cases/**/*.py 的严正则读取器'
                    '（r\'<%s>(.*?)</%s>\'）会取不到该数据块，'
                    'mod_vtu_writer 故意写成裸标签'
                    % (tag, sec, sorted(e.attrib), sec, sec))
    for da in e.findall('DataArray'):
        nm = da.get('Name') or '?'
        full = '%s.%s[%s]' % (tag, sec, nm)
        try:
            nc = int(da.get('NumberOfComponents') or 1)
        except ValueError:
            err.append('%s: 非法 NumberOfComponents="%s"'
                       % (full, da.get('NumberOfComponents')))
            continue
        k = _da_len(da, full, err, warn)
        if k is not None and nc > 0 and k != nc * nent:
            err.append('%s: 数据长度 %d != NumberOfComponents(%d)*%d'
                       % (full, k, nc, nent))


def _check_piece(piece, tag):
    """校验一个 <Piece>，返回 (errors, warnings)。"""
    err, warn = [], []
    npts = _int_attr(piece, 'NumberOfPoints', tag, err)
    ncell = _int_attr(piece, 'NumberOfCells', tag, err)
    if npts is None or ncell is None:
        return err, warn

    # ---- points ----
    n = _da_len(_da(piece.find('Points')), tag + '.Points', err, warn)
    if n is not None and n != 3 * npts:
        err.append('%s: Points 数据 %d 个 != 3*NumberOfPoints(%d)'
                   % (tag, n, npts))

    # ---- cells ----
    cells = piece.find('Cells')
    if cells is None:
        err.append('%s: 缺少 <Cells>' % tag)
        return err, warn
    conn = _da_ints(_da(cells, 'connectivity'), tag + '.connectivity', err, warn)
    offs = _da_ints(_da(cells, 'offsets'), tag + '.offsets', err, warn)
    typs = _da_ints(_da(cells, 'types'), tag + '.types', err, warn)

    if typs is not None and len(typs) != ncell:
        err.append('%s: types 长度 %d != NumberOfCells(%d)'
                   % (tag, len(typs), ncell))
    if offs is not None:
        if len(offs) != ncell:
            err.append('%s: offsets 长度 %d != NumberOfCells(%d)'
                       % (tag, len(offs), ncell))
        if any(offs[i] <= offs[i - 1] for i in range(1, len(offs))):
            err.append('%s: offsets 不是严格递增' % tag)
    if conn is not None and offs is not None and len(offs) > 0:
        if offs[-1] != len(conn):
            err.append('%s: offsets 末值 %d != connectivity 长度 %d'
                       % (tag, offs[-1], len(conn)))
        bad = [v for v in conn if v < 0 or v >= npts]
        if bad:
            err.append('%s: connectivity 含越界节点号 %s (合法范围 [0,%d))'
                       % (tag, bad[:4], npts))

    # ---- 每单元节点数（核心：防 legacy 前导 token 回归） ----
    if offs is not None and typs is not None and len(offs) == len(typs):
        prev = 0
        legacy, wrong, unknown = 0, 0, 0
        exleg, exwrong = [], []
        for i, (o, t) in enumerate(zip(offs, typs)):
            nv = o - prev
            prev = o
            want = CELL_NNODE.get(t)
            if want is None:
                unknown += 1
                continue
            if nv == want + 1:
                legacy += 1
                if len(exleg) < 3:
                    exleg.append('#%d(type %s)' % (i, CELL_NAME[t]))
            elif nv != want:
                wrong += 1
                if len(exwrong) < 3:
                    exwrong.append('#%d(type %s: 差 %d != %d)'
                                   % (i, CELL_NAME[t], nv, want))
        if legacy:
            err.append('%s: %d/%d 个单元的 offsets 差 == 规范节点数+1 -> '
                       'connectivity 里混入了 legacy VTK 顶点数前导 token！'
                       '（例 %s）' % (tag, legacy, len(offs), ', '.join(exleg)))
        if wrong:
            err.append('%s: %d 个单元的 offsets 差 != 规范节点数（例 %s）'
                       % (tag, wrong, ', '.join(exwrong)))
        if unknown:
            warn.append('%s: %d 个单元的类型码未知（可变节点数），已跳过节点数校验'
                        % (tag, unknown))

    # ---- point / cell data ----
    _check_data_section(piece, 'PointData', npts, tag, err, warn)
    _check_data_section(piece, 'CellData', ncell, tag, err, warn)
    return err, warn


def validate_text(text, name='<memory>'):
    """校验一段 .vtu 文本，返回 (errors, warnings) 两个字符串列表。"""
    err, warn = [], []
    try:
        root = ET.fromstring(text)
    except ET.ParseError as e:
        return ['XML 解析失败: %s' % e], warn
    if root.tag != 'VTKFile':
        return ['根元素是 <%s>，应为 <VTKFile>' % root.tag], warn
    if root.get('type') != 'UnstructuredGrid':
        err.append('<VTKFile type="%s">，应为 "UnstructuredGrid"' % root.get('type'))
    if root.get('version') != '1.0':
        warn.append('<VTKFile version="%s">（期望 "1.0"）' % root.get('version'))
    grid = root.find('UnstructuredGrid')
    if grid is None:
        err.append('缺少 <UnstructuredGrid>')
        return err, warn
    pieces = grid.findall('Piece')
    if not pieces:
        err.append('缺少 <Piece>')
    for i, p in enumerate(pieces):
        e, w = _check_piece(p, '%s Piece[%d]' % (name, i))
        err += e
        warn += w
    return err, warn


def validate_file(path):
    """读文件并校验，返回 (errors, warnings)。"""
    with open(path, 'r', errors='replace') as f:
        return validate_text(f.read(), name=os.path.basename(path))


def vtk_check(path, err, warn):
    """可选：用 vtk 模块（与 ParaView 同源的读取器）复查单元类型/节点数/体积。"""
    try:
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy
    except ImportError:
        warn.append('%s: 未安装 vtk 模块，跳过 --vtk 复查'
                    % os.path.basename(path))
        return
    rd = vtk.vtkXMLUnstructuredGridReader()
    rd.SetFileName(path)
    rd.Update()
    g = rd.GetOutput()
    npts, ncell = g.GetNumberOfPoints(), g.GetNumberOfCells()
    if npts <= 0 or ncell <= 0:
        err.append('%s: vtk 读取到 %d 点 / %d 单元' % (path, npts, ncell))
        return
    pts = vtk_to_numpy(g.GetPoints().GetData())
    badtype, badnn, negvol, first = 0, 0, 0, []
    for i in range(ncell):
        if g.GetCellType(i) != vtk.VTK_TETRA:
            badtype += 1
            continue
        c = g.GetCell(i)
        ids = [c.GetPointId(k) for k in range(c.GetNumberOfPoints())]
        if len(ids) != 4:
            badnn += 1
            continue
        a, b, cc, d = (pts[ids[0]], pts[ids[1]], pts[ids[2]], pts[ids[3]])
        v = (np_dot(np_cross(b - a, cc - a), d - a)) / 6.0
        if v <= 0.0:
            negvol += 1
            if len(first) < 3:
                first.append((i, v))
    print('  vtk: npts=%d ncells=%d (全部 %s)' % (npts, ncell, 'TETRA'
          if badtype == 0 else '含非 TETRA'))
    if badtype:
        err.append('%s: vtk 读到 %d 个非 TETRA 单元' % (path, badtype))
    if badnn:
        err.append('%s: vtk 读到 %d 个节点数 != 4 的 TETRA 单元' % (path, badnn))
    if negvol:
        err.append('%s: vtk 读到 %d 个非正体积四面体，例 %s' % (path, negvol, first))


def np_cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def np_dot(a, b):
    return sum(x * y for x, y in zip(a, b))


# ---------------- 守卫自测：坏样例必须被识别为不合规 ------------------------
_TPL = """<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="5" NumberOfCells="2">
      <Points><DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="ascii">
0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.2 0.2 0.2</DataArray></Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">%s</DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">%s</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">10 10</DataArray>
      </Cells>
      <CellData>
        <DataArray type="Float64" Name="pressure" format="ascii">1.0 2.0</DataArray>
      </CellData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
"""
GOOD = _TPL % ('0 1 2 4 0 2 3 4', '4 8')
BAD_LEGACY = _TPL % ('4 0 1 2 4 4 0 2 3 4', '5 10')     # 2026-10-15 前的缺陷
BAD_SHORT = _TPL % ('0 1 2 4 0 2 3', '4 7')             # offsets/连通性不一致
BAD_RANGE = _TPL % ('0 1 2 9 0 2 3 4', '4 8')           # 节点号越界
# 带属性的 <CellData>：格式合法，但仓库内严正则读取器取不到 -> 必须出告警
ATTR_CD = GOOD.replace('<CellData>', '<CellData Scalars="pressure">')


def selftest():
    """守卫自测：规范样例 PASS，三类坏样例必须被识别（退出码 0/1）。"""
    ok = True
    e, w = validate_text(GOOD, 'good')
    if e:
        print('SELFTEST FAIL: 规范样例被误判: %s' % e)
        ok = False
    else:
        print('selftest: 规范样例 PASS（期望）')
    if w:
        print('SELFTEST FAIL: 规范样例出现告警: %s' % w)
        ok = False
    for nm, txt, key in (('legacy token', BAD_LEGACY, '前导 token'),
                         ('offsets 错', BAD_SHORT, 'offsets'),
                         ('节点号越界', BAD_RANGE, '越界')):
        e, _ = validate_text(txt, nm)
        if any(key in s for s in e):
            print('selftest: %-12s 被正确判为不合规（期望）' % nm)
        else:
            print('selftest: %-12s 未被识别，检查失效！%s' % (nm, e))
            ok = False
    e, w = validate_text(ATTR_CD, 'attr')
    if not e and any('裸标签' in s for s in w):
        print('selftest: 带属性数据块   触发读取兼容性告警（期望）')
    else:
        print('selftest: 带属性数据块   未触发兼容性告警！errors=%s warn=%s'
              % (e, w))
        ok = False
    print('selftest %s' % ('OK: 格式守卫有效' if ok else 'FAILED'))
    return 0 if ok else 1


def main(argv):
    quiet = '-q' in argv or '--quiet' in argv
    use_vtk = '--vtk' in argv
    files = [a for a in argv if not a.startswith('-')]
    if '--selftest' in argv:
        return selftest()
    if not files:
        print('用法: python3 tools/validate_vtu.py [--vtk] [-q] FILE.vtu [...]')
        print('      python3 tools/validate_vtu.py --selftest')
        return 2
    rc = 0
    for f in files:
        if not os.path.isfile(f):
            print('  ERROR %s: 文件不存在' % f)
            rc = 1
            continue
        err, warn = validate_file(f)
        if use_vtk:
            vtk_check(f, err, warn)
        if not quiet:
            for w in warn:
                print('  WARN  %s' % w)
        for e in err:
            print('  ERROR %s' % e)
        print('%s  %s  (%d errors, %d warnings)'
              % ('PASS' if not err else 'FAIL', f, len(err), len(warn)))
        if err:
            rc = 1
    return rc


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
