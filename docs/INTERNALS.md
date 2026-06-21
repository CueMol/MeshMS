# MeshMS — indexing & data-layout conventions

MeshMS implements the **MolSurfComp** SES algorithm (Quan & Stamm), whose data
structures are deeply index-interdependent. To reproduce them faithfully the C++
code keeps the algorithm's original **1-based indexing** rather than renumbering;
this document is the reference for that convention and the record layouts the
internal stages pass around. (The public API in `meshms/meshms.hpp` is fully
0-based and dense — these conventions are internal to the library.)

Pipeline order:

```
read_xyzr -> interstructure -> data_I_Cir -> data_Seg_Pat -> data_ext
          -> convex patches + concave patches + toroidal patches
          -> accumulate mesh -> orient_faces -> (V, F, N)
```

## Indexing convention (CRITICAL)

To mirror the index-interdependent algorithm we use **1-based indexing** internally:

- **Entity arrays** carry a dummy slot at index 0; real entities are `1..N`.
- **Coordinate arrays** (centers, `I` points): dummy ROW 0, real columns 0,1,2 =
  x,y,z. Access a whole point as `I[s]` (a length-3 vector); `C[i]` is atom `i`'s
  center.
- **Record matrices** (circle, segment, ncrasegment, segment0, direction, Iijk,
  loops, patches, satom, Ii, circleindex, …) carry a dummy ROW 0 **and** a dummy
  COLUMN 0, so a record's fields live at columns `1..k`.
- Convert to 0-based ONLY for the final PLY face indices (handled by `MeshState`),
  which is also what the public `MeshResult.faces` exposes.
- A circular-next index `mod(m, N) + 1` (1-based) stays `m % N + 1` in this scheme
  — do NOT shift it. Be careful when porting modular index arithmetic.

## 1-based record layouts

The internal structs mirror these field layouts (columns are 1-based; column 0 is
the dummy):

- `Geom` — `M`, `centers[(M+1,3)]`, `R[(M+1,)]`; built by `read_xyzr`.
- `Para` — `radius_probe`, `mesh_size`.
- `Neighbors` (`interstructure`) — `M_int[(M+1,kmax+1)]`, `num_int[(M+1,)]`;
  `M_int[i, 1:num_int[i]+1]` = sorted neighbour atoms of atom `i` (CSR-backed).
- `DataI` (`data_I_Cir`):
  - `I[(nI+1,3)]`, `nI`, `direction[(nI+1,4)]` cols 1..3, `Iijk[(nI+1,4)]`
    cols 1..3 = i,j,k, `Ii[(M+1,*)]`, `In[(M+1,)]`,
    `I_circle[(M+1,kmax+1,*)]`, `I_circle_num[(M+1,kmax+1)]`,
    `high_I[(nI+1,)]` (1 = non-singular triangle, 0 = singularity candidate),
    `hightvalue[(nI+1,)]`.
- `DataCir` (`data_I_Cir`):
  - `circle[(ncircle+1,10)]`: cols `[_, i, j, A(3)=3:6, n(3)=6:9, r=9]`, `ncircle`,
    `circleindex[(M+1,11)]`, `ncircleindex[(M+1,)]`.
- `DataSeg` (`data_Seg_Pat`):
  - `segment[(nseg+1,6)]` cols `[_, i, j, p1, p2, direct]` (p1,p2 index into `I`;
    direct ∈ {1,-1}); `ncrasegment[(nseg+1,9)]` cols `[_, n(3)=1:4, A(3)=4:7, r=7,
    radian=8]`; `satom[(M+1,*)]`, `nsatom[(M+1,)]`.
- `DataLoop` (`data_Seg_Pat`):
  - `loops[(nloops+1,smax+1)]` (entries = global segment indices),
    `loopsize[(nloops+1,)]`, `nloops`, `loops_index[(M+1,3)]` cols `[_,start,end]`.
- `DataPat` (`data_Seg_Pat`):
  - `patches[(npatch+1,pmax+1)]` (signed entries: `+k` = loop index local to atom,
    `-k` = circle index local to atom), `patchesize[(npatch+1,)]`, `npatches`,
    `patches_index[(M+1,3)]` cols `[_,start,end]`, `patch_atom[(npatch+1,)]`.

`MeshState` is the single mesh sink for ALL patch families:
`add_patch(P, T, face_normals)` where `P[1:]` are real `(3,)` points (dummy
`P[0]`), `T` is an iterable of 1-based `(a,b,c)` index triples into `P`, and
`face_normals` is aligned with `T` (or empty).

## Spherical-patch mesher contract

The convex and concave drivers both mesh through one spherical-patch mesher. Its
inputs are a sphere center/radius, the local loops (lists of segments), the
segment table `segment0` (cols `[_, c(3)=1:4, n(3)=4:7, r=7, spoint(3)=8:11,
angle=11]`, `n` pointing OUTWARD), an optional circle table `circle0`, and a
signed `patches` list (`+k` → `loops[k]`, `-k` → `circle0[k]`). Key behaviour:

- `tolerance = 0.8 * min(d, r_sphere)`.
- For **convex** patches (`r_sphere = R_i + Rp > Rp`) the boundary and sampled
  points are contracted toward the sphere center by `(r_sphere - Rp) / r_sphere`
  (the SAS→SES contraction). For **concave** patches (`r_sphere = Rp`) there is no
  contraction.
- Normals: convex `arg_NV = +1` (outward = away from atom center), concave
  `arg_NV = -1` (outward = toward probe center), then `add_patch(P, T, NV)`.

## Numerical faithfulness (do not "optimize" away)

These are deliberate invariants that keep the mesh bit-stable under the DEFAULT
`MESHMS_FP=strict` policy; see `meshms/vec3.hpp`. A `MESHMS_FP=fast` deploy build
intentionally drops the two *faithfulness* invariants (`pysq`, the contraction ban)
and is verified by the equivalence gate `tests/test_fp_gate.cpp` instead of the
golden suite. The *robustness* invariants below (the `acos` clip, the `sqrt` guard)
hold in BOTH policies and must never be removed.

When adding code, keep the faithful expression on the strict side of the `#if`:
never write a bare `x*x` where the reference squared a scalar with `x**2`. The
pattern to follow is `angle_sphere` / `angle_vectors` in `src/meshing.cpp` and
the coverage tests in `src/sas.cpp`: the `#else` branch is the reference
implementation verbatim, the `MESHMS_FP_FAST` branch is the rewrite, and the
comment states the algebraic identity that justifies it.

- Clip every `acos` argument to `[-1, 1]` (`acos_clamped`).
- Guard `sqrt(max(x, 0))` where the original silently took the real part.
- `pysq(x)` (= `pow(x, 2.0)`) is used wherever the reference squared a SCALAR with
  `x**2`; it differs from `x*x` by up to 1 ULP, which is enough to flip a `<`/`<=`
  boundary test and change the SAS triple-point set. Keep the float term/evaluation
  order identical to the reference in dot/cross/circlecenter.
- `-ffp-contract=off` is mandatory in the strict policy (see
  [`OPTIMIZATION.md`](OPTIMIZATION.md)): an FMA contraction would change the
  last-bit rounding and break the golden gate. `MESHMS_FP=fast` enables it
  deliberately. `-ffast-math`/`-Ofast` are rejected in BOTH policies -- they imply
  `-ffinite-math-only`, which folds `isfinite()` away and so disables the NaN
  tripwire the deploy gate rests on.
