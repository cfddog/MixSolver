# M6-wing bit-level regression baseline

Persistent replacement for the lost `/tmp/m6reg` baseline (issue ④): the
ephemeral `/tmp` copy was wiped, so M6-wing bit-level regression could no longer
be reproduced.  The baseline is now regenerated and stored **inside the repo**.

## What this checks

The in-tree structured solver (`bin/struct_solver`, ser tree, `mpif90` compiled)
must reproduce the pristine OpenCFD-EC 1.16a reference **byte-for-byte** on the
4-block M6-wing case, `np1`, `t_end=0.501`, `Kstep_save=50`.

Reference binary: `external/OpenCFD-EC-1.16a` (pristine 1.16a), built with
`mpif90 -O2 -std=legacy -ffree-line-length-none` after a GBK→UTF-8 conversion.

## Provenance / result (2026-10-06)

| item | value |
|------|-------|
| reference build flags | `mpif90 -O2 -std=legacy -ffree-line-length-none` |
| run | `mpirun -np 1`, 51 steps (save at Kstep=50), `t_end=0.501` |
| `flow3d.dat`        | 13,600,032 B, md5 **dc134a2d196422043ecad7c86ac8f898** |
| in-tree vs reference | flow3d/SA3d/wall_dist/partation-auto/part_grid, Step_mess/bc3d.inc/mesh-quality **all identical** |

The `flow3d.dat` md5 is **identical to the historical `/tmp/m6reg` baseline**
(the only md5 ever recorded for it), i.e. the lost baseline is faithfully
restored.

Only known benign difference: `output_para.out` (not in the compared set) — the
pristine binary prints 5 extra lines for the phase-11 removed keys
(`IF_TurboMachinary`, `Turbo_*`, `Ref_medium_usrdef`, `a0`).  `a0`/`d0` are
uninitialised prints in the pristine code and never feed the computation.

## Layout

```
regress/m6wing/
  control.ec            regression control (t_end=0.501, Kstep_save=50)
  run_regression.sh     self-contained rebuild + run + compare script
  README.md             this file
  baseline/             pristine reference outputs (the baseline)
    flow3d.dat          <- primary bit-level baseline (13.6 MB)
    SA3d.dat wall_dist.dat partation-auto.dat part_grid.dat
    Step_mess.dat bc3d.inc mesh-quality.dat
    md5sums.txt         md5 manifest for all of the above
```

## Running it

```bash
# verify the current solver against the stored baseline (fast, no ref rebuild)
regress/m6wing/run_regression.sh

# also rebuild + rerun the pristine reference from source, then compare
M6_REBUILD_REF=1 regress/m6wing/run_regression.sh
```

Knobs: `M6_WORK` (scratch dir, default `/tmp/m6wing_regress`), `M6_NP` (ranks,
default 1), `M6_OPT` (reference flags), `M6_J` (make -j).

## Notes

* The in-tree structured solver reads `Mesh3d.x`; the pristine reference reads
  `Mesh3d.dat`.  The script copies the case mesh under both names.
* The pristine binary never calls `MPI_Finalize`, so its `mpirun` exit code is
  non-zero (pre-existing; see `memory-bank/constraints.md`).  The script
  tolerates this and judges success from the file comparison only.
* M6 input mesh is only shipped as `Mesh3d.dat` in
  `external/OpenCFD-EC-1.16a/cases/M6-wing/` (13.7 MB) and is intentionally not
  duplicated here.
