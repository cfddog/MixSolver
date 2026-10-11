# Static libraries required to compile and run UCNS3D

Precompiled **metis** and **parmetis** static libraries. They can be used to link
during compilation.

- `metis/libmetis.a` : serial METIS, used by the MPI build on rank 0 (always linked).
- `parmetis/`        : prebuilt for the MPICH ABI; NOT linked by default (see the
  note in the top-level `Makefile`).

The formerly vendored `tecplot/libtecio.a` (TecIO) has been **removed**: the
binary Tecplot `.plt` writer was deleted from the solver, which now writes ASCII
VTU output only, so the archive is neither needed nor linked.   
