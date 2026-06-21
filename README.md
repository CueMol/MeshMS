# MeshMS / libMeshMS

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
no `CMakeLists.txt` edit needed. Options: `MESHMS_WITH_CGAL` (default OFF, zero
third-party deps), `MESHMS_SANITIZE`, `MESHMS_NATIVE`, `MESHMS_TBB` (default ON;
oneTBB from the cuemol2 deplibs — point `CMAKE_PREFIX_PATH` / `TBB_DIR` at it, or
`-DMESHMS_TBB=OFF` for a serial, dependency-free build).

## Using the library

`find_package(MeshMS REQUIRED)` then link `MeshMS::MeshMS`. The ONLY header a
consumer includes is the C++17-clean facade `<meshms/capi.hpp>` — it trades only
in `std` types, so no MeshMS-internal type crosses the boundary, and a C++17
program links the C++20-built library directly.

```cpp
#include <meshms/capi.hpp>
using namespace meshms;

std::vector<std::array<double,4>> atoms = /* {x,y,z,radius} per atom */;

// One-shot: build the SES mesh (fuse=true => welded, watertight C1 surface).
MeshResult m = build_surface_from_array(atoms, /*probe=*/1.4, /*mesh_size=*/0.5, /*fuse=*/true);
// m.verts, m.vnormals, m.faces (0-based triangles), m.atom_id (per-vertex owner).

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
| `analyze_mesh(mesh)` → `MeshReport` | vertices/faces, area, signed volume, watertightness, boundary/non-manifold edge counts |
| `boundary_diagnostics(mesh)` → loops + non-manifold edges | localize open holes/seams (centroid, edge count, closed?) |
| `version()` | library version string |

`MeshResult` always carries `vnormals` (per-vertex normals) and, on the build and
`remove_flaps` paths, `atom_id` (1-based per-vertex owning atom; value `i` refers
to `atoms[i-1]`, `0` = unknown) so a consumer can colour the surface by atom.

## CLI

```sh
meshms_cli INPUT.xyzr -o OUT.ply [--probe 1.4] [--mesh-size 0.5]
           [--no-fuse] [--fuse-cusps] [--weld-tol 1e-4] [--vertex-normals] [-v]
```

`meshms_cli` reads an `xyzr` file, builds the surface, and writes an ASCII PLY.
It is implemented entirely on top of `<meshms/capi.hpp>` (argument/PLY/xyzr I/O
and reporting are the only things it does itself) — the worked example of
consuming libMeshMS the way cuemol2/3 does.

## Status & validation

The faithful, CGAL-free port is complete. Every pipeline stage is cross-checked
against a frozen reference oracle (per-stage integer topology EXACT and floats
within 1e-9, or exact-match meshes); 18/18 ctest green.

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
| `capi` (public facade) | the API above | `test_capi`, `test_capi_post`, `test_rscache`, `test_atomid` |

End-to-end the mesh reproduces the reference SES: topology EXACT, geometry within
1e-9 (bit-exact on 4/6 golden molecules; the rest differ only at the last ULP from
libm transcendental rounding). The optional close-cusps path produces a watertight
2-manifold.

Out of scope (deferred): analytic area/volume, the symmetry-jitter degenerate
fallback, and all CGAL acceleration (power-diagram).

## Reference docs & fixtures

- [`PORTING_CONTRACT.md`](PORTING_CONTRACT.md) — the indexing & data-layout
  conventions the code follows (1-based-with-dummy-row-0).
- [`OPTIMIZATION.md`](OPTIMIZATION.md) — the 28-core, output-preserving parallel design.

`tests/data/` holds the `xyzr` input molecules and `tests/ref/` the frozen golden
oracle the tests cross-check against. The fixtures are pre-generated and committed;
the scripts that regenerate them live with the original reference implementation
and are not part of this C++ library.
