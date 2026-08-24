# C++ optimization — 28-core parallelism (output-preserving)

Profiling-driven, **no CGAL**. Every optimization is
**output-preserving**: all exact-match ctest stay green under both the parallel
(oneTBB) and serial (`-DMESHMS_TBB=OFF`) build, and clean under ASan/UBSan + 10× repeat
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
had over Python, **101M is ~1500× faster than the Python reference** (30.76 s → 20.6 ms).

## How (output-preserving parallelism)

The pattern throughout: do the expensive **independent per-entity** work in parallel into
**per-entity buffers**, then **merge serially in ascending entity-index order** so the
global numbering / append order is identical to serial. No float reductions, no reordered
accumulation, no changed float term order or `pysq` usage.

| stage | % serial | technique | stage speedup @101M |
|---|--:|---|--:|
| `data_I_Cir` | 55% | two-phase: `meshms::parallel_for` per-atom event discovery (read-only) → **serial replay** (zero float math) reproducing the exact `s`-numbering + `Ii`/`I_circle` order | ~12–13× |
| `interstructure` | 8% | parallel per-atom neighbour-row build (const grid reads) + serial CSR concat | ~6× |
| convex mesh | 15% | per-atom parallel; **thread-local `segment0`/`Rj`** (fixes a shared-buffer race) + ordered `add_patch` | ~5× |
| concave mesh | 13% | per-probe + per-simple-triangle parallel (call-local scratch) + ordered merge | ~3× |
| toroidal mesh | 7% | per-segment / per-circle parallel + ordered merge | ~1.7× |

Enabling refactor: the meshers (`mesh_sphpat`, `mesh_toroide`, `mesh_cusp`) now **return a
`LocalMesh`** instead of writing to `MeshState`, so a driver meshes patches in parallel and
calls `add_patch` in order (skipping empty meshes exactly as before).

Single-thread micro-opts: `col_of` map → `std::lower_bound` over the sorted CSR row;
per-frame `Dist` scratch reuse in the advancing front; `reserve()` hints.

## Build flags: the FP policy (critical)

`MESHMS_FP` selects the floating-point semantics the library is compiled with. It is the
only knob that changes the numbers; `-O3`, `-march=native` and oneTBB never do.

| | `strict` (default) | `fast` (deploy) |
|---|---|---|
| GNU / Clang | `-ffp-contract=off` | `-ffp-contract=fast -fno-math-errno -fno-trapping-math -fno-signed-zeros -freciprocal-math` |
| MSVC | `/fp:precise /fp:contract:off` | `/fp:precise /fp:contract` |
| `pysq(x)` | `pow(x, 2.0)` | `x * x` |
| gated by | the golden suite (bit-for-bit) | `test_fp_gate` (equivalence) |

**Why strict has to be strict.** With FMA contraction, `a*b+c` fuses, the last-bit rounding
changes, and the ULP-sensitive SAS boundary tests (`disc <= 0`, `pysq(q+rij)+pp-Rext2 < 0`,
`c > 0` in `sas.cpp`) flip. That changes the triple-point set and with it the whole integer
topology. `pysq` is the same story from the source side: `pow(x, 2.0)` and `x*x` differ by up
to 1 ULP on ~0.08% of inputs.

**What fast is allowed to break, and what it is not.** `fast` gives up the *faithfulness*
constructs -- `pysq` and the contraction ban. It does NOT touch the *robustness* constructs:
`acos_clamped`, `sqrt(max(x, 0))` and the divide-by-zero guards stay in both policies.

**Flags deliberately excluded from `fast`, and rejected outright at configure time**
(`-ffast-math`, `-Ofast`, `-ffinite-math-only`, `-funsafe-math-optimizations`,
`-fapprox-func`, `/fp:fast`; override with `MESHMS_ALLOW_UNSAFE_FP_FLAGS=ON`):

- **`-ffinite-math-only`** folds `isfinite()`/`isnan()` to a constant. The deploy gate's
  primary tripwire is `MeshReport::nonfinite_vertices`, so this flag would disable exactly
  the check that tells a consumer the mesh is unusable rather than merely different.
- **`-ffast-math` / `-Ofast` / `-funsafe-math-optimizations`** all match GCC's crtfastmath.o
  link spec (`%{Ofast|ffast-math|funsafe-math-optimizations:%{!shared:crtfastmath.o%s}}`),
  which sets FTZ/DAZ for the **whole host process**. `libMeshMS.a` is linked into cuemol2;
  silently changing the host application's FP mode is not ours to do. The individual flags
  `fast` does use match no such spec.
- **`-fassociative-math`** (implied by `-funsafe-math-optimizations`) buys almost nothing
  here -- the parallel design deliberately has no float reductions -- while it *would*
  reassociate the sequential sums in `mesh_area`/`signed_volume`, the very quantities the
  gate measures.

`meshms/vec3.hpp` carries an `#error` on `__FAST_MATH__` and `mesh_check.cpp` one on
`__FINITE_MATH_ONLY__`, so the ban also holds when the flag arrives from outside CMake.

### Measured (Apple M2, AppleClang 15, Release + `-march=native` + oneTBB, best-of-15)

| mol | d | strict | fast | gain |
|---|--:|--:|--:|--:|
| 1crn | 0.5 | 3.95 ms | 3.89 ms | ~0% |
| barstar | 0.5 | 10.40 ms | 10.06 ms | 3% |
| 101M | 0.5 | 13.81 ms | 13.03 ms | 4% |
| 101M | 0.25 | 30.33 ms | 28.92 ms | 4% |
| barstar | 0.25 | 19.54 ms | 18.87 ms | 3% |

**On macOS the gain is small, and the reason is specific**: `nm` shows a *strict* build of
`libMeshMS.a` contains no `pow` reference at all -- AppleClang already folds `pow(x, 2.0)` at
all 56 `pysq` sites. So `-fno-math-errno` has nothing left to win, and the source-level
`pysq` switch is a no-op there. Full `-ffast-math` was measured too and reached only 4-5%,
so the excluded flags are not where the time is either.

**On GCC and MSVC the picture should differ**: neither folds `pow(x, 2.0)` under the default
FP settings, so `fast` removes 56 real libm calls -- 11 of them in `data_I_Cir` (55% of the
serial time) and 16 in the concave mesher (13%). On MSVC that source-level switch is the
*only* lever, because `/fp:precise` never folds `pow` and MSVC has no `-fno-math-errno`.

Contraction needs the instruction to exist. On Apple Silicon FMA is always there; on x86-64
the default baseline has none, so `/fp:contract` and `-ffp-contract=fast` emit nothing until
the target says otherwise -- pass `MESHMS_ARCH=avx2` for an Intel/AMD deploy build.

None of this has been measured here; measure before assuming, on the target platform, with
`meshms_bench` (its banner prints the policy).

Other options: `MESHMS_NATIVE` (`-march=native -O3`), `MESHMS_ARCH` (ISA baseline for
redistributable builds; `avx2` expands to `-march=x86-64-v3` / `/arch:AVX2` on x86-64 and is
ignored elsewhere, any other value passes through verbatim; `native` bakes in the build
machine's ISA and SIGILLs elsewhere), `MESHMS_LTO` (off by default; a static library built
with LTO carries bitcode and requires a matching consumer toolchain),
`MESHMS_TBB` (default ON; `meshms::parallel_for` falls back to a serial loop when
`MESHMS_WITH_TBB` is undefined, so a TBB-off build runs serially and identically).

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DMESHMS_NATIVE=ON \
      -DMESHMS_TBB=ON -DCMAKE_PREFIX_PATH=<deplibs root>   # oneTBB from cuemol2 deplibs
cmake --build build && ctest --test-dir build          # 21/21, parallel (oneTBB)
cmake -S . -B build-serial -G Ninja -DMESHMS_TBB=OFF
ctest --test-dir build-serial                          # 21/21, serial (same output)

# deploy build, and the A/B against it
cmake -S . -B build-fast -G Ninja -DCMAKE_BUILD_TYPE=Release -DMESHMS_FP=fast
cmake --build build-fast && ctest --test-dir build-fast   # 9 skipped, 12 green
BENCH_REPS=15 ./build/meshms_bench      tests/data 101M:0.5 101M:0.25 barstar:0.25
BENCH_REPS=15 ./build-fast/meshms_bench tests/data 101M:0.5 101M:0.25 barstar:0.25
```

### Safety valves in a relaxed-FP build

The frame loop of the advancing front is already bounded (`while (k <= N)`, `k` always
increments). The unbounded path is the mutual recursion
`advancing_front_approach` <-> `collapse_nonneighbor{1,2}_sphere`, where
`collapse_nonneighbor2_sphere` can hand the callee a *larger* front. `Ctx::depth` +
`kMaxFrontDepth` (4096) caps it: a frame that hits the cap returns without closing its
front, which leaves a hole the gate sees as extra boundary edges -- a detectable degradation
rather than a stack overflow. A strict build never reaches it (the golden suite is the
proof), and the guard costs one compare per frame, never per triangle.

`geom_from_array` rejects non-finite coordinates and negative radii at the facade boundary
(radius 0 stays valid -- barstar.xyzr has 503 such atoms), and `manifold_report` reports
`nonfinite_vertices` so a consumer can check a deploy-built mesh with one `analyze_mesh()`
call.

## Round 3a (2026-08): allocation + serial-merge pass (bit-exact)

Gated like every round: golden ctest green on the TBB and serial builds, unchanged
bytes. Three commits (Apple M2, best-of-15, 101M d=0.25; cumulative TBB 30.1 -> 23.1 ms
(-23%), serial 116 -> 96.7 ms (-17%)):

- **Allocator traffic** (was ~30% of serial self time; tbbmalloc_proxy measured
  SLOWER than the system allocator on macOS, so the fix is fewer allocations, not a
  different allocator): ae_from_rows reserves; the sweep-2 Dist scratch moved into
  Ctx and stopped being zero-filled (the assign()'s O(Nae) memset per iteration was
  pure waste); arc_division tests its early return before generating and pushes
  straight into P/Ae; add_patch appends by bulk insert; the concave simple-triangle
  pass compacts its sparse iteration space.
- **Patch-parallel scatter merge** (was ~14 of 31 ms at 8T): merge_local_meshes
  keeps the exact merge order via a serial prefix sum of the base indices, then
  copies patches into disjoint slices in a parallel_for, freeing each dead
  LocalMesh on the worker that copied it. torus went from 1.7x to ~2.5x scaling.
- **Advancing front as a ring buffer**: the per-triangle front rewrites are O(1)
  head/tail moves instead of full O(Nae) rebuilds; splits use per-recursion-depth
  spare rings (a deque -- parents hold references into earlier slots).

Also fixed: tools/bench.cpp timed a ~10 MB V/F deep copy inside the orient span
(the real pipeline moves them), inflating that column ~6x.

## Remaining serial tail (future, smaller payoff)

At 28T@101M the serial-ish remainder is `data_Seg_Pat` (~1.6 ms), `data_ext` (flood-fill,
~0.9 ms), `orient_faces` (~0.6 ms) and the `data_I_Cir` serial replay. `data_Seg_Pat` could
take the same two-phase parallel treatment (per-atom segment/loop/patch discovery + serial
renumber) for ~7% more; `data_ext` is inherently serial. The advancing-front mesher is the
least-scaling parallel stage (load-imbalanced large patches) — a cell-list would speed it
but is byte-fragile and was deliberately **not** attempted (it changes the
`index_np` tie-break → output). The cell-list and `data_Seg_Pat` parallel are the next levers
if more is needed.

## Comparison vs MSMS

MSMS 2.6.1 (Sanner 1994 — the established reduced-surface SES mesher) is single-threaded.
Comparison at **matched mesh resolution** (MSMS `-density` tuned so its vertex count matches
MeshMS at `mesh-size 0.6`), **compute-time only** (MSMS internal "Total Time real" with
`-no_area`; MeshMS `build_surface`; file I/O excluded on both), best-of-N:

| molecule | atoms | MeshMS nV | MeshMS serial | MeshMS 28T | MSMS nV | MSMS | **vs serial** | **vs 28T** |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| ArgArg  |   38 |  2464 |   2.5 ms | 1.5 ms |  2444 |  20 ms | 8.1× | 13.6× |
| 1ETN    |  160 |  9895 |  13.0 ms | 4.4 ms |  9874 |  60 ms | 4.6× | 13.6× |
| 1crn    |  327 | 19973 |  24.8 ms | 7.4 ms | 21325 | 110 ms | 4.4× | 14.9× |
| 1B17    |  483 | 28113 |  36.5 ms | 8.4 ms | 30504 | 160 ms | 4.4× | 19.1× |
| 101M    | 1413 | 63198 | 111.5 ms | 23.5 ms| 70973 | 370 ms | 3.3× | 15.8× |
| barstar | 1426 | 41847 |  85.8 ms | 15.4 ms| 43180 | 250 ms | 2.9× | 16.3× |

- **MeshMS single-threaded is ~3× faster than MSMS** at protein scale (101M 3.3×, barstar
  2.9×; more on small molecules where MSMS's fixed overhead dominates, ArgArg 8×).
- **MeshMS with 28 cores is ~13–19× faster than MSMS.**
- Both compute the same SES (cross-check: ArgArg SES area MSMS 261.7 vs MeshMS 260.2 Å²).

Fairness notes: vertex counts matched within ~5–12% (MSMS slightly higher on the large
molecules, i.e. slightly more work — MeshMS still wins); MSMS's internal timer is 10 ms
resolution so the small-molecule MSMS values are coarse; MSMS is a single-thread 1994 code,
so MeshMS wins on BOTH a faster serial baseline and 28-core parallelism. (MSMS additionally
computes analytic area/volume by default — disabled here via `-no_area` for a meshing-only
comparison; MeshMS defers analytic area/volume entirely in this build.)

## Round 2 (2026-08): output-preserving micro/structural pass

A second optimization pass, gated exactly like the first: all golden ctest green on
BOTH the TBB and the serial (`-DMESHMS_TBB=OFF`) build, plus FNV-hash snapshots of
V/F/N/NV/atom_id/ftype on 8 large molecules (carbo 6609 / ubch5b 2364 /
glutaredoxin 1276 / 101M / barstar / 1B17 / 1crn / 1GZI; d=0.5 and 0.25, fuse
on/off) bit-identical before and after every change, on both builds. Measured on an
8-core Apple M-series (best-of-N; the MolSurfComp UFF-radius set, which is
neighbour-denser than the 101M-style xyzr set, so `data_I_Cir` dominates more).

What landed (all output-preserving):

- **pipeline**: MeshState arrays moved (not copied) into Surface; `orient_faces`
  parallel per face.
- **mesh_check**: ordered maps -> hash maps (uint64 edge keys); core templated so
  `manifold_report`/`boundary_loops` gained facade-layout overloads -- the capi
  `analyze_mesh`/`boundary_diagnostics` no longer convert the whole mesh.
- **weld**: `fill_small_holes`/`remove_nonmanifold_flaps` take V by value and move
  it through (V was copied unchanged).
- **MeshState**: merge-phase `reserve_extra` (sizes pre-summed from the LocalMesh
  slots) and an `add_patch` overload that MOVES the per-vertex TagLists.
- **toroidal**: exact `T` reserves; the per-face analytic normal fused into the
  face-orientation loop (the intermediate `normal` array is gone).
- **data_I_Cir**: `Rext/Rext2 = R+Rp / pysq(R+Rp)` tables (pysq is a real libm
  `pow` on GCC -- Apple clang folds it, so this line item pays on the Linux/GCC
  targets); the coverage pre-loop's `sqrt(1-ctheta^2)*dist_a` / `pysq(ctheta*dist_a)`
  cached per row and reused by the whole-circle test; integer guards before the
  discarded norms in the Y1/Y2 loops.
- **data_ext**: owner map-of-maps -> flat 2-slot arrays; dense neighbor rectangle
  -> CSR. Stage -70%.
- **data_Seg_Pat**: `interiorloop` per-atom precompute (unit vectors, second acos
  term, alpha2 -- all point-independent); tree sets moved + reserved; the atom loop
  split into a serial segment pass and a **parallel per-atom loop/patch pass**
  (satom[a] is complete once atom a's rows are done) with a serial ascending-i
  merge. Stage -30..-38% (TBB).
- **convex**: the per-atom full-size `segment0`/`Rj` alloc+memset (O(M x nsegment)
  per build) replaced by persistent thread_local scratch -- no clearing needed
  because every read row is rewritten per atom. Stage -50%.
- **concave**: `build_neighbors` parallel per probe (fixed slots, no merge).
- **RS cache (feature)**: the density-independent concave probe decomposition is
  now precomputed into `RSComponents` (`precompute_concave` -> `ConcaveDecomp`;
  `SESconcavepat_mesh` re-meshes it), so multi-density `build_mesh` calls skip it.
- **fuse_by_id**: two-phase -- complete flat (tag,cell) grid, then PARALLEL
  per-vertex partner discovery (j<i pairs: the identical pair set and float
  expression), then a serial union replay; the union-find partition is
  order-independent (component root == minimum member), so the output is
  bit-identical.

Tried and REVERTED (serial A/B showed a regression): caching sweep-1 distances /
flags for sweep 2 of the advancing front (+8..19% on convex/concave -- the
per-iteration scratch fills outweigh the rarely-taken recomputes).

End-to-end effect (TBB, 8 cores, best-of-7, d=0.5): total -17..-28% (101M
17.6 -> 13.5 ms, barstar 13.3 -> 10.3 ms, carbo 141 -> 117 ms); d=0.25 -13..-21%.
Serial totals -5..-8% (dominated by `data_I_Cir`, whose float work is already
minimal under clang's pow folding). Per-density re-mesh from a cached RS
(`build_mesh` only): **fuse=true -48..-60%** (101M d=0.5: 42.5 -> 18.0 ms, carbo:
103.5 -> 41.0 ms), fuse=false -5..-11%.

Remaining levers: the `data_I_Cir` coverage/verification loops (a margin-guarded
triangle-inequality prune is outcome-safe but was not attempted); an advancing-front
ring buffer for the O(Nt x Nae) front copies in `addnewpoint_sphere`; hoisting the
(small, post-fix) convex `mod_seg_loop_cir` slice work into the RS cache.
