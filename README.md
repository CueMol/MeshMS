# MeshMS / libMeshMS

Faithful C++20 port of the `molsurf` SES (Solvent Excluded Surface) pipeline,
packaged as the static library `libMeshMS.a` (CMake target `MeshMS`, namespace
`meshms::`). It reads an `xyzr` atom file and produces a watertight, triangulated
PLY of the molecular SES — no third-party C++ dependencies by default.

Reference docs in this repo:

- [`PORTING_CONTRACT.md`](PORTING_CONTRACT.md) — algorithm semantics + the
  1-based-with-dummy-row-0 indexing convention used throughout the code.
- [`OPTIMIZATION.md`](OPTIMIZATION.md) — the 28-core, output-preserving parallel design.

**Faithful-first: zero CGAL.** Acceptance is *equivalence* (watertight, counts
exact, area/volume within ~1e-6 of the Python golden), not byte-identity.

## Build & test

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

The CMake uses `file(GLOB CONFIGURE_DEPENDS)`: every `src/*.cpp` joins the
`MeshMS` library, every `tests/test_*.cpp` becomes its own CTest case — add files,
no `CMakeLists.txt` edit needed. Options: `MESHMS_WITH_CGAL` (default OFF, zero
third-party deps), `MESHMS_SANITIZE`, `MESHMS_NATIVE`, `MESHMS_OPENMP` (default ON).

## Use as a library

`find_package(MeshMS REQUIRED)` then link `MeshMS::MeshMS`. The only header an
external consumer needs is the C++17-clean facade `<meshms/capi.hpp>`
(`compute_rs_from_array`, `build_mesh_from_cache`, `build_surface_from_array`).
The library is compiled with C++20 but the facade links cleanly into a C++17 program.

## Status — faithful CGAL-free port COMPLETE

Every Python module is ported and cross-checked against a frozen Python golden
(per-stage integer topology EXACT + floats within 1e-9, or exact-match meshes):

| C++ unit | ports | verified against Python |
|---|---|---|
| `vec3.hpp` (+ `pysq`) | `mathutil.py` | `test_vec3` (term-order faithful) |
| `geom`,`params` | `geom.py`,`params.py` | `test_geom` |
| `intersection` (→ CSR `Neighbors`) | `intersection.py` | `test_intersection` — exact `M_int` |
| `mesh_check`,`ply` | `check_mesh.py`,`ply.py` | `test_mesh_check` — gate on the golden mesh |
| `sas` (CSR/append) | `sas.py` `data_I_Cir` | `test_sas` — DataI/DataCir exact |
| `sas_patches` | `sas_patches.py` `data_Seg_Pat` | `test_segpat` — seg/loop/patch exact |
| `exterior` | `exterior.py` `data_ext` | `test_ext` — eSAS flags exact |
| `mesh_state`,`meshing` | `mesh_state.py`,`meshing.py` | `test_meshing` — mesh_sphpat exact replay |
| `concave` | `concave.py` `SESconcavepat` | `test_concave` — V/F/N exact |
| `convex` (+`mod_seg_loop_cir`) | `convex.py` | `test_convex` — V/F/N exact |
| `toroidal` | `toroidal.py` (toroide+cusp) | `test_toroidal` — V/F/N exact (incl. cusp2/3) |
| `pipeline` | `pipeline.py` | `test_surface` — **end-to-end vs GOLDEN_MESHES** |
| `fusion` | `fusion.py` `fuse_by_id` | `test_fusion` — fused mesh exact (closes C1 seams) |
| `weld` | `weld.py` | `test_weld` — weld/loops/fill/flaps exact |
| `capi` (public facade) | — | `test_capi`, `test_rscache`, `test_atomid` |
| `meshms_cli` | `__main__.py` | PLY **byte-identical** to `python -m molsurf` |

**End-to-end:** `build_surface` reproduces the Python SES — topology EXACT, geometry
within 1e-9 (bit-exact on 4/6 golden molecules; 3spheres/1YJO differ only at the last
ULP from libm/numpy transcendental rounding), and the CLI's `%.6f` PLY is byte-identical.
The optional fuse path closes all C1 seams to a watertight 2-manifold. 17/17 ctest green.

Out of scope (deferred): the analytic area/volume (`DataAV`), the numpy-RNG
symmetry-jitter fallback, the large-scale memory-ledger sweep, and all CGAL
acceleration (power-diagram).

## Golden reference fixtures

`tests/ref/` holds the frozen Python oracle the C++ tests cross-check against
(per stage: `*.interstructure.txt`, `*.sas.txt`, `*.segpat.txt`, `*.ext.txt`,
`*.meshcalls.txt`, `*.concave.txt`, `*.convex.txt`, `*.toroidal.txt`, `*.surface.txt`,
`*.fused.txt`, `*.weld.txt`, plus `*.golden.ply` + `*.gate.txt`). `tests/data/`
holds the `xyzr` input molecules. The fixtures are pre-generated and committed; the
`dump_*.py` scripts that regenerate them depend on the original Python `molsurf`
package (the PoC repo) and are not part of this C++-only library.
