#!/usr/bin/env python3
"""Phase-10 bitwise restart consistency check.

Compares a continuous 6-iter run against a 3-iter + restart + 3-iter run:
  - struct side: flow3d.dat  (unformatted; compared byte-wise after the run)
  - uns side:    unMesh_restart.dat (ASCII es24.16 dump)

Byte-identical files => bit-level consistency.
"""
import hashlib, os, sys

def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def strip_gen_line(path):
    """Remove the '# generated:' timestamp line from the uns dump header."""
    with open(path, 'r') as f:
        lines = f.readlines()
    out = [l for l in lines if not l.startswith('# generated:')]
    return ''.join(out)

base = os.path.dirname(os.path.abspath(__file__))
pairs = [
    ('flow3d.dat',        'run A (cont)', 'flow3d.dat',        'run B (restart)', True),
    ('unMesh_restart.dat','run A (cont)', 'unMesh_restart.dat','run B (restart)', False),
]

dirA = os.path.join(base, 'cont6')
dirB = os.path.join(base, 'restart33')

ok = True
for fa, la, fb, lb, raw in pairs:
    pa, pb = os.path.join(dirA, fa), os.path.join(dirB, fb)
    if not (os.path.exists(pa) and os.path.exists(pb)):
        print(f'MISSING: {pa} or {pb}'); ok = False; continue
    if raw:
        ha, hb = md5(pa), md5(pb)
    else:
        ha = hashlib.md5(strip_gen_line(pa).encode()).hexdigest()
        hb = hashlib.md5(strip_gen_line(pb).encode()).hexdigest()
    same = (ha == hb)
    ok = ok and same
    print(f'{"IDENTICAL" if same else "DIFFER"}: {fa}  ({la} vs {lb})')
    print(f'   {ha}  vs  {hb}')

sys.exit(0 if ok else 1)
