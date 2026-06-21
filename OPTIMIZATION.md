# C++ optimization — 28-core parallelism (output-preserving)

Branch `feature/cpp-optimize`. Profiling-driven, **no CGAL**. Every optimization is
**output-preserving**: all 14 exact-match ctest stay green under both the parallel
(`OMP_NUM_THREADS=28`) and serial (`=1`) build, and clean under ASan/UBSan + 10× repeat
runs. The mesh is bit-for-bit identical to the pre-optimization faithful port.

## Result

`build_surface(fuse=false)`, native `-O3 -march=native -ffp-contract=off`, best-of-7.

| molecule | atoms | serial (1T) | parallel (28T) | **speedup** |
|---|--:|--:|--:|--:|
| ArgArg  |   38 |   2.85 ms |  1.74 ms | 1.6× |
| 1YJO    |   67 |   6.88 ms |  2.70 ms | 2.6× |
| 1ETN    |  160 |  14.08 ms |  3.92 ms | 3.6× |
| 1crn    |  327 |  25.84 ms |  6.74 ms | 3.8× |
| 1B17    |  483 |  44.39 ms | 11.33 ms | 3.9× |
| 101M    | 1413 | 109.28 ms | 20.56 ms | **5.3×** |
| barstar | 1426 |  91.21 ms | 18.20 ms | **5.0×** |

Larger molecules clear ~5×; small ones are Amdahl-bound by the fixed per-call overhead
and the remaining serial stages. Combined with the ~300× the faithful C++ port already
had over Python, **101M is ~1500× faster than the numpy reference** (30.76 s → 20.6 ms).

## How (output-preserving parallelism)

The pattern throughout: do the expensive **independent per-entity** work in parallel into
**per-entity buffers**, then **merge serially in ascending entity-index order** so the
global numbering / append order is identical to serial. No float reductions, no reordered
accumulation, no changed float term order or `pysq` usage.

| stage | % serial | technique | stage speedup @101M |
|---|--:|---|--:|
| `data_I_Cir` | 55% | two-phase: `#pragma omp parallel for` per-atom event discovery (read-only) → **serial replay** (zero float math) reproducing the exact `s`-numbering + `Ii`/`I_circle` order | ~12–13× |
| `interstructure` | 8% | parallel per-atom neighbour-row build (const grid reads) + serial CSR concat | ~6× |
| convex mesh | 15% | per-atom parallel; **thread-local `segment0`/`Rj`** (fixes a shared-buffer race) + ordered `add_patch` | ~5× |
| concave mesh | 13% | per-probe + per-simple-triangle parallel (call-local scratch) + ordered merge | ~3× |
| toroidal mesh | 7% | per-segment / per-circle parallel + ordered merge | ~1.7× |

Enabling refactor: the meshers (`mesh_sphpat`, `mesh_toroide`, `mesh_cusp`) now **return a
`LocalMesh`** instead of writing to `MeshState`, so a driver meshes patches in parallel and
calls `add_patch` in order (skipping empty meshes exactly as before).

Single-thread micro-opts: `col_of` map → `std::lower_bound` over the sorted CSR row;
per-frame `Dist` scratch reuse in the advancing front; `reserve()` hints.

## Build flags (critical)

`-ffp-contract=off` is **mandatory and always on**: with `-march=native`, g++ otherwise
fuses `a*b+c` into an FMA, changing the last-bit rounding and flipping the ULP-sensitive SAS
boundary tests (the bit-for-bit golden gate). Options: `MESHMS_NATIVE` (`-march=native -O3`),
`MESHMS_OPENMP` (default ON; `#pragma omp` is `#if defined(_OPENMP)`-guarded so an
OpenMP-off build runs serially and identically).

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DMESHMS_NATIVE=ON -DMESHMS_OPENMP=ON
cmake --build build && ctest --test-dir build          # 17/17, OMP_NUM_THREADS=28
OMP_NUM_THREADS=1 ctest --test-dir build               # 17/17, serial (same output)
# bench: g++ -std=c++20 -O3 -march=native -fopenmp -ffp-contract=off -I include \
#   tools/bench.cpp build/libMeshMS.a -o /tmp/bench ; OMP_NUM_THREADS=28 /tmp/bench tests/data 101M:0.6
```

## Remaining serial tail (future, smaller payoff)

At 28T@101M the serial-ish remainder is `data_Seg_Pat` (~1.6 ms), `data_ext` (flood-fill,
~0.9 ms), `orient_faces` (~0.6 ms) and the `data_I_Cir` serial replay. `data_Seg_Pat` could
take the same two-phase parallel treatment (per-atom segment/loop/patch discovery + serial
renumber) for ~7% more; `data_ext` is inherently serial. The advancing-front mesher is the
least-scaling parallel stage (load-imbalanced large patches) — a cell-list (`CPP_PORT_DESIGN.md`
§6) would speed it but is byte-fragile and was deliberately **not** attempted (it changes the
`index_np` tie-break → output). The cell-list and `data_Seg_Pat` parallel are the next levers
if more is needed.

## Comparison vs MSMS

MSMS 2.6.1 (Sanner 1994 — the established reduced-surface SES mesher) is single-threaded.
Comparison at **matched mesh resolution** (MSMS `-density` tuned so its vertex count matches
molsurf at `mesh-size 0.6`), **compute-time only** (MSMS internal "Total Time real" with
`-no_area`; molsurf `build_surface`; file I/O excluded on both), best-of-N:

| molecule | atoms | molsurf nV | molsurf serial | molsurf 28T | MSMS nV | MSMS | **vs serial** | **vs 28T** |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| ArgArg  |   38 |  2464 |   2.5 ms | 1.5 ms |  2444 |  20 ms | 8.1× | 13.6× |
| 1ETN    |  160 |  9895 |  13.0 ms | 4.4 ms |  9874 |  60 ms | 4.6× | 13.6× |
| 1crn    |  327 | 19973 |  24.8 ms | 7.4 ms | 21325 | 110 ms | 4.4× | 14.9× |
| 1B17    |  483 | 28113 |  36.5 ms | 8.4 ms | 30504 | 160 ms | 4.4× | 19.1× |
| 101M    | 1413 | 63198 | 111.5 ms | 23.5 ms| 70973 | 370 ms | 3.3× | 15.8× |
| barstar | 1426 | 41847 |  85.8 ms | 15.4 ms| 43180 | 250 ms | 2.9× | 16.3× |

- **molsurf single-threaded is ~3× faster than MSMS** at protein scale (101M 3.3×, barstar
  2.9×; more on small molecules where MSMS's fixed overhead dominates, ArgArg 8×).
- **molsurf with 28 cores is ~13–19× faster than MSMS.**
- Both compute the same SES (cross-check: ArgArg SES area MSMS 261.7 vs molsurf 260.2 Å²).

Fairness notes: vertex counts matched within ~5–12% (MSMS slightly higher on the large
molecules, i.e. slightly more work — molsurf still wins); MSMS's internal timer is 10 ms
resolution so the small-molecule MSMS values are coarse; MSMS is a single-thread 1994 code,
so molsurf wins on BOTH a faster serial baseline and 28-core parallelism. (MSMS additionally
computes analytic area/volume by default — disabled here via `-no_area` for a meshing-only
comparison; molsurf defers analytic area/volume entirely in this build.)
