# MeshMS / libMeshMS

[![CI](https://github.com/CueMol/MeshMS/actions/workflows/ci.yml/badge.svg)](https://github.com/CueMol/MeshMS/actions/workflows/ci.yml)

A C++ library that computes the **Solvent Excluded Surface (SES)** of a molecule
analytically and returns it as a watertight, triangulated mesh. It is packaged as
the static library `libMeshMS.a` (CMake target `MeshMS`, namespace `meshms::`).
Parallelism uses **Intel oneTBB** (from the cuemol2 deplibs bundle); a serial
build (`-DMESHMS_TBB=OFF`) has **no third-party C++ dependencies**.

The surface algorithm is a faithful implementation of **MolSurfComp** (Quan &
Stamm, 2016/2017): the SES is decomposed into convex spherical, toroidal, and
concave spherical patches, each meshed analytically, then assembled into one mesh.

The **primary consumer is cuemol2/cuemol3**; the bundled `meshms_cli` is a thin
reference front-end that drives the library through its public API only.

## Build & test

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

The CMake uses `file(GLOB CONFIGURE_DEPENDS)`: every `src/*.cpp` joins the
`MeshMS` library, every `tests/test_*.cpp` becomes its own CTest case — add files,
no `CMakeLists.txt` edit needed.

| option | default | effect |
|---|---|---|
| `MESHMS_FP` | `strict` | FP policy: `strict` (bit-exact, golden-gated) or `fast` (deploy) — see **Numerical modes** |
| `MESHMS_TBB` | `ON` | oneTBB from the cuemol2 deplibs — point `CMAKE_PREFIX_PATH` / `TBB_DIR` at it, or `OFF` for a serial, dependency-free build |
| `MESHMS_NATIVE` | `OFF` | `-march=native -O3` |
| `MESHMS_ARCH` | *(empty)* | ISA baseline for redistributable builds. `avx2` is portable — it expands to `-march=x86-64-v3` / `/arch:AVX2` on x86-64 and is ignored elsewhere, so one build script can pass it everywhere. Any other value goes through verbatim as `-march=<v>` / `/arch:<v>`. Use this, not `native` |
| `MESHMS_LTO` | `OFF` | interprocedural optimization. A static library built with LTO carries bitcode, so the consumer must use a matching toolchain |
| `MESHMS_SANITIZE` | `OFF` | ASan/UBSan for the fragile mesher |
| `MESHMS_WITH_CGAL` | `OFF` | reserved; the CGAL gate cross-checker is not built |

## Numerical modes

| mode | for | golden bit gate | what it guarantees |
|---|---|---|---|
| `strict` (default) | development, bug reproduction, reference comparison | on (21/21) | reproduces the reference implementation to the last bit |
| `fast` | deploy (cuemol2/cuemol3) | off (9 tests skip) | never throws, no NaN vertices, `close_cusps` still watertight, area/volume within 0.1% of strict |

`fast` enables FMA contraction and a small unsafe-math subset, and switches
`pysq(x)` from `pow(x, 2.0)` to `x*x`. It gives up bit-for-bit reproducibility, so
the nine golden bit-regression tests skip themselves there and `tests/test_fp_gate.cpp`
— an equivalence gate against a frozen strict baseline — takes over. `-ffast-math`
and `-Ofast` are rejected in **both** modes (they would disable the NaN tripwire and
set FTZ/DAZ for the whole host process); see [`OPTIMIZATION.md`](OPTIMIZATION.md).

```sh
# deploy build for linking into cuemol2/cuemol3
cmake -S . -B build-deploy -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DMESHMS_FP=fast -DMESHMS_ARCH=avx2 \
      -DMESHMS_TBB=ON -DCMAKE_PREFIX_PATH=<deplibs root> \
      -DCMAKE_INSTALL_PREFIX=<install root>
cmake --build build-deploy && cmake --install build-deploy
```

Notes for a deploy build: `MESHMS_ARCH=avx2` is the right x86-64 baseline —
without FMA in the target ISA, neither `/fp:contract` nor `-ffp-contract=fast` can
emit a single FMA, so `fast` loses half its point on the default x86-64 baseline;
it is ignored on Apple Silicon, which always has FMA. Never use `native` for a
redistributable build (it bakes in the build machine's ISA and SIGILLs elsewhere).
Enable `MESHMS_LTO` only when MeshMS and the final binary use the same toolchain.
Reproduce any bug on a `strict` build, which is what the golden suite covers.
`build_info()` reports which policy a binary was built with, and the CLI stamps it
into the `.vert`/`.face` headers.

## Using the library

`find_package(MeshMS REQUIRED)` then link `MeshMS::MeshMS`. The ONLY header a
consumer includes is the C++17-clean facade `<meshms/meshms.hpp>` — it trades only
in `std` types, so no MeshMS-internal type crosses the boundary, and a C++17
program links the C++20-built library directly. See [`docs/API.md`](docs/API.md)
for the complete public API reference.

```cpp
#include <meshms/meshms.hpp>
using namespace meshms;

std::vector<std::array<double,4>> atoms = /* {x,y,z,radius} per atom */;

// One-shot: build the SES mesh (fuse=true => welded, watertight C1 surface).
MeshResult m = build_surface_from_array(atoms, /*probe=*/1.4, /*mesh_size=*/0.5, /*fuse=*/true);
// m.verts, m.vnormals, m.faces (0-based triangles), m.atom_id (per-vertex owner),
// m.face_type (per-face SES component: 3=convex, 2=concave, 1=toroidal; MSMS codes).

// Multi-density: precompute the geometry once, mesh cheaply at many densities.
auto rs = compute_rs_from_array(atoms, 1.4);
MeshResult coarse = build_mesh_from_cache(rs, 0.5);
MeshResult fine   = build_mesh_from_cache(rs, 0.25);
```

The facade also exposes the post-processing and inspection a consumer needs to
finish and validate a surface:

| API | purpose |
|---|---|
| `compute_rs_from_array`, `build_mesh_from_cache` | density-independent precompute, then mesh per density |
| `build_surface_from_array` | one-shot mesh from an atom array |
| `close_cusps(mesh, weld_tol)` | weld + fan-fill cusp/singular seams into a fully closed manifold (drops `atom_id`) |
| `remove_flaps(mesh)` | drop spurious doubled-"flap" triangles (no-op on a clean mesh; keeps `atom_id`) |
| `vertex_normals(verts, faces)` | area-weighted per-vertex normals (to refresh after a consumer-side edit) |
| `analyze_mesh(mesh)` → `MeshReport` | vertices/faces, area, signed volume, watertightness, boundary/non-manifold edge counts, non-finite vertex count |
| `boundary_diagnostics(mesh)` → loops + non-manifold edges | localize open holes/seams (centroid, edge count, closed?) |
| `version()`, `build_info()` | semver string; and the build configuration (FP policy, parallel backend) for logs and bug reports |

`MeshResult` always carries `vnormals` (per-vertex normals) and, on the build and
`remove_flaps` paths, `atom_id` (1-based per-vertex owning atom; value `i` refers
to `atoms[i-1]`, `0` = unknown) so a consumer can colour the surface by atom, plus
`face_type` (per-face SES component, aligned with `faces`: `3` convex/contact,
`2` concave/spherical-reentrant, `1` toroidal — the MSMS face-type codes; carried
through `remove_flaps`, empty after `close_cusps`).

**Multi-component inputs** are handled by the facade: the atoms are split into
connected components of the SAS-intersection graph, the (unchanged) pipeline runs
per component, and the meshes are concatenated. Isolated atoms — single atoms
whose SAS sphere touches nothing, which the faithful exterior extraction cannot
represent — are meshed directly as full vdW spheres (icosphere at the requested
`mesh_size`). So a protein plus stray waters/ions, or several detached chains,
all appear in one mesh, with `atom_id` still referring to the input array. A
single-component input takes the exact pre-existing pipeline path bit-for-bit.

## CLI

```sh
meshms_cli INPUT.xyzr -o OUT [--probe 1.4] [--mesh-size 0.5]
           [--no-fuse] [--fuse-cusps] [--weld-tol 1e-4] [--ply] [--vertex-normals] [-v]
```

`meshms_cli` reads an `xyzr` file, builds the surface, and writes it. The **default
output is the MSMS `.vert`/`.face` pair** (`<OUT>.vert` + `<OUT>.face`; any
`.vert`/`.face`/`.ply` suffix on `OUT` is stripped to form the base name) — each
vertex line is `x y z  nx ny nz  0  closest_sphere  vertex_type` and each face line
`v1 v2 v3 (1-based)  face_type  0`, with the SES component in `closest_sphere`
(= `atom_id`) and the type fields. `--ply` (or an `OUT` ending in `.ply`) writes a
single ASCII PLY instead. It is implemented entirely on top of `<meshms/meshms.hpp>`
(argument/PLY/MSMS/xyzr I/O and reporting are the only things it does itself) — the
worked example of consuming libMeshMS the way cuemol2/3 does.

## Status & validation

The faithful, CGAL-free port is complete. Every pipeline stage is cross-checked
against a frozen reference oracle (per-stage integer topology EXACT and floats
within 1e-9, or exact-match meshes); 21/21 ctest green.

| C++ unit | computes | golden test |
|---|---|---|
| `vec3` | fused length-3 vector math | `test_vec3` |
| `geom`, `params` | xyzr input, run parameters | `test_geom` |
| `intersection` | SAS sphere-intersection graph (CSR `Neighbors`) | `test_intersection` (exact `M_int`) |
| `sas` | SAS arrangement (`data_I_Cir`) | `test_sas` (DataI/DataCir exact) |
| `sas_patches` | segment/loop/patch decomposition (`data_Seg_Pat`) | `test_segpat` |
| `exterior` | exterior-SAS extraction (`data_ext`) | `test_ext` |
| `mesh_state`, `meshing` | advancing-front spherical-patch mesher | `test_meshing` |
| `concave` | concave SES patches (`SESconcavepat`) | `test_concave` (V/F/N exact) |
| `convex` | convex SES patches | `test_convex` (V/F/N exact) |
| `toroidal` | toroidal patches + cusps | `test_toroidal` (incl. cusp2/3) |
| `pipeline` | end-to-end SES assembly | `test_surface` (vs golden meshes) |
| `fusion` | ID-based boundary fusion (`fuse_by_id`) | `test_fusion` |
| `weld`, `mesh_check` | weld/fill/flap repair + manifold report | `test_weld`, `test_mesh_check` |
| `meshms` (public facade) | the API above | `test_meshms`, `test_meshms_post`, `test_rscache`, `test_atomid` |
| *(end-to-end, both FP modes)* | no throw / no NaN / bounded vertices / aggregates vs a frozen baseline / `close_cusps` watertight | `test_fp_gate` (11 molecule×density cases) |

End-to-end the mesh reproduces the reference SES: topology EXACT, geometry within
1e-9 (bit-exact on 4/6 golden molecules; the rest differ only at the last ULP from
libm transcendental rounding). The optional close-cusps path produces a watertight
2-manifold.

A `MESHMS_FP=fast` build is held to `test_fp_gate` instead: the nine golden
bit-regression tests skip themselves, the other eleven plus the gate stay green.
The gate compares aggregates (vertex/face counts within 2%, area and signed volume
within 0.1%, no increase in duplicate faces or non-manifold edges) against a
strict-built baseline, and checks the properties that need no baseline at all — the
build does not throw, no vertex is non-finite, every vertex stays inside the atom
bounding box, and `close_cusps` still closes what it closed before. A consumer of a
fast build should call `analyze_mesh()` once and confirm `nonfinite_vertices == 0`.

Strongly symmetric molecules (e.g. fullerene) whose faithful coordinates collapse
a toroidal patch into a NaN-poisoned mesh are handled by an automatic
symmetry-jitter fallback: `build_surface_from_array` defaults to `Jitter::Auto`,
which meshes faithfully first and only re-meshes with small deterministic center
perturbations when the result is non-finite (clean molecules are unaffected and
stay bit-exact). See [`docs/API.md`](docs/API.md) for the `Jitter` enum. The
jittered mesh is std-lib-RNG dependent and thus not byte-stable across platforms,
which is why fullerene stays excluded from `test_fp_gate`.

Out of scope (deferred): analytic area/volume and all CGAL acceleration
(power-diagram).

## Reference docs & fixtures

- [`docs/API.md`](docs/API.md) — the complete public API reference (types,
  functions, data contracts, output formats).
- [`docs/INTERNALS.md`](docs/INTERNALS.md) — the indexing & data-layout
  conventions the code follows (1-based-with-dummy-row-0).
- [`docs/OPTIMIZATION.md`](docs/OPTIMIZATION.md) — the 28-core, output-preserving parallel design.

`tests/data/` holds the `xyzr` input molecules and `tests/ref/` the frozen golden
oracle the tests cross-check against. The fixtures are pre-generated and committed;
the scripts that regenerate them live with the original reference implementation
and are not part of this C++ library.
