# MeshMS — public API reference

The complete reference for the public C++ facade, `<meshms/meshms.hpp>`. This is
the ONLY header an external consumer needs to include. It is **C++17-clean** (it
includes none of the internal C++20 MeshMS headers) and trades only in `std`
types (`std::vector` / `std::array` / `std::uint32_t`), so no MeshMS-internal
type ever crosses the boundary: a C++17 program can include this header and link
the C++20-built `libMeshMS.a` directly. Everything lives in `namespace meshms`.

For the library's internal indexing & data-layout conventions (the 1-based-with-
dummy-row-0 record matrices the pipeline passes around), see
[`INTERNALS.md`](INTERNALS.md). They are internal — the public API below is fully
0-based and dense.

## 1. Overview & usage

Building a surface, then optionally post-processing it, follows this order:

```
build_surface_from_array   ─┐
            or              ├─►  [ close_cusps ]  ─►  [ remove_flaps ]
build_mesh_from_cache      ─┘     (watertight)         (repair flaps)
```

There are two ways to build. Use the **one-shot** path for a single mesh, or the
**cached** path when you mesh the same molecule at several densities — it
precomputes the density-independent RS components once and re-runs only the
(cheap) mesher per density.

```cpp
#include <meshms/meshms.hpp>
using namespace meshms;

std::vector<std::array<double,4>> atoms = /* {x,y,z,radius} per atom */;

// One-shot: build the SES mesh (fuse=true => welded, watertight C1 surface).
MeshResult m = build_surface_from_array(atoms, /*probe=*/1.4, /*mesh_size=*/0.5, /*fuse=*/true);

// Multi-density: precompute the geometry once, mesh cheaply at many densities.
auto rs = compute_rs_from_array(atoms, 1.4);
MeshResult coarse = build_mesh_from_cache(rs, 0.5);
MeshResult fine   = build_mesh_from_cache(rs, 0.25);
```

The one-shot `build_surface_from_array` is bit-for-bit identical (verts/faces) to
building from the equivalent `.xyzr` file, and to the cached path at the same
density.

## 2. Types

### `MeshResult`

The triangle mesh returned to the caller. All arrays are **0-based and dense**.

```cpp
struct MeshResult {
  std::vector<std::array<double, 3>> verts;     // vertex positions
  std::vector<std::array<double, 3>> vnormals;  // per-vertex outward normals (aligned with verts)
  std::vector<std::array<std::uint32_t, 3>> faces;  // triangles, 0-based indices into verts
  std::vector<std::uint32_t> atom_id;           // per-vertex owning atom (1-based; 0 = unknown)
  std::vector<std::uint8_t>  face_type;         // per-face SES component code (3/2/1)
};
```

- `verts` / `vnormals` are aligned: `vnormals[i]` is the area-weighted outward
  normal at `verts[i]`. `MeshResult` always carries `vnormals`.
- `faces` are triangles, 0-based indices into `verts`.
- `atom_id` is the per-vertex owning atom, aligned with `verts`. Value `i` is
  **1-based** (`0` = unknown) and refers to `atoms[i-1]` — the `(i-1)`-th entry
  of the input array the caller passed. The caller maps this index back to its
  own atom id. Carried on the build and `remove_flaps` paths; **empty** after
  `close_cusps` (see the carry rules in §3).
- `face_type` is the per-face SES component, aligned with `faces`, using the
  MSMS face-type codes: `3` = convex (contact), `2` = concave (spherical
  reentrant), `1` = toroidal (toric reentrant). Populated by the build; carried
  (filtered) through `remove_flaps`; **empty** after `close_cusps`.

### `RSCache`

```cpp
struct RSCache;   // opaque; consumers hold only a std::shared_ptr<RSCache>
```

Opaque handle to the precomputed, density-independent RS components. The
definition lives in `meshms.cpp`; consumers only ever hold a `shared_ptr` to it,
produced by `compute_rs_from_array` and consumed by `build_mesh_from_cache`.

### `Jitter`

```cpp
enum class Jitter { None, Auto };
```

Symmetry-jitter policy for the one-shot `build_surface_from_array` path
(default `Auto`):

- **`None`** — mesh the faithful coordinates exactly (bit-for-bit reproducible).
- **`Auto`** — mesh faithfully first and, only if that produces a degenerate,
  NaN-poisoned mesh (a strongly symmetric molecule such as fullerene), retry with
  small deterministic Gaussian perturbations of the atom centers until the mesh
  is finite.

The jittered retry is std-lib-RNG based and therefore **not reproducible across
platforms**, and it applies **only to the one-shot path** — the cache
(`compute_rs_from_array` / `build_mesh_from_cache`) is always faithful. Clean
molecules never trigger a retry, so `Auto` is identical to `None` for them.

### `MeshReport`

Whole-mesh quality report returned by `analyze_mesh`.

```cpp
struct MeshReport {
  std::uint32_t n_vertices = 0;
  std::uint32_t n_faces = 0;
  std::uint32_t degenerate_faces = 0;   // repeated index or |cross| < 1e-12
  std::uint32_t duplicate_faces = 0;    // same vertex triple seen more than once
  std::uint32_t boundary_edges = 0;     // edges used by exactly one face
  std::uint32_t nonmanifold_edges = 0;  // edges used by >= 3 faces
  bool watertight = false;              // boundary == 0 && nonmanifold == 0 && nF > 0
  double area = 0.0;                    // sum of triangle areas
  double signed_volume = 0.0;          // divergence-theorem signed volume
};
```

### `BoundaryLoopInfo`, `NonmanifoldEdgeInfo`, `BoundaryDiagnostics`

The localized boundary/non-manifold diagnostics returned by
`boundary_diagnostics` (these back the CLI's `-v` output).

```cpp
// One connected component of the boundary-edge subgraph (a hole/seam).
struct BoundaryLoopInfo {
  std::uint32_t n_verts = 0;
  std::uint32_t n_edges = 0;
  bool closed = false;               // every member has exactly 2 boundary neighbours
  std::array<double, 3> centroid{};  // mean of the component's vertex coordinates
};

// A non-manifold edge (u, v) shared by `count` >= 3 faces.
struct NonmanifoldEdgeInfo {
  std::uint32_t u = 0;
  std::uint32_t v = 0;
  std::uint32_t count = 0;
};

struct BoundaryDiagnostics {
  std::vector<BoundaryLoopInfo>     loops;       // sorted largest-first by n_edges
  std::vector<NonmanifoldEdgeInfo>  nonmanifold; // first-seen order
};
```

## 3. Functions

Nine functions make up the API. Each post-processing function takes a
`MeshResult` and returns a new one; `vnormals` in the result are recomputed
(area-weighted) from the returned geometry, so they stay valid.

| function | signature (defaults) | purpose |
|---|---|---|
| `compute_rs_from_array` | `(const std::vector<std::array<double,4>>& xyzr, double radius_probe) -> std::shared_ptr<RSCache>` | density-independent precompute from an atom array |
| `build_mesh_from_cache` | `(const std::shared_ptr<RSCache>& rs, double mesh_size, bool fuse = false) -> MeshResult` | mesh the SES from a cache at the given density |
| `build_surface_from_array` | `(const std::vector<std::array<double,4>>& xyzr, double radius_probe, double mesh_size, bool fuse = false, Jitter jitter = Jitter::Auto) -> MeshResult` | one-shot: `compute_rs_from_array` + `build_mesh_from_cache`, with an optional symmetry-jitter NaN fallback |
| `version` | `() -> const char*` | library version string, e.g. `"0.1.0"` |
| `close_cusps` | `(const MeshResult& mesh, double weld_tol = 1e-4) -> MeshResult` | weld + fan-fill cusp/singular seams into a watertight 2-manifold |
| `remove_flaps` | `(const MeshResult& mesh, int passes = 4) -> MeshResult` | drop spurious doubled-"flap" (non-manifold) triangles |
| `vertex_normals` | `(const std::vector<std::array<double,3>>& verts, const std::vector<std::array<std::uint32_t,3>>& faces) -> std::vector<std::array<double,3>>` | area-weighted per-vertex outward normals |
| `analyze_mesh` | `(const MeshResult& mesh) -> MeshReport` | whole-mesh quality report |
| `boundary_diagnostics` | `(const MeshResult& mesh) -> BoundaryDiagnostics` | localize open holes/seams and non-manifold edges |

Details:

- **`compute_rs_from_array`** — `xyzr[a]` is `{x, y, z, radius}` for atom `a`
  (0-based input; an output `atom_id` value `i` refers to `xyzr[i-1]`).
  `radius_probe` is the solvent probe radius Rp.
- **`build_mesh_from_cache`** — `fuse=true` welds tagged boundary vertices and
  the result then carries per-vertex normals recomputed from the fused geometry.
- **`build_surface_from_array`** — one-shot convenience. With `jitter ==
  Jitter::None` (or for any molecule that meshes cleanly) the result is bit-for-
  bit identical (verts/faces) to building from the equivalent `.xyzr` file. The
  default `jitter == Jitter::Auto` adds a degenerate-symmetry fallback: it meshes
  the faithful coordinates first and, only if that mesh is NaN-poisoned (a
  strongly symmetric molecule such as fullerene), silently re-meshes with small
  deterministic Gaussian center perturbations until the mesh is finite. The retry
  is std-lib-RNG based and **not cross-platform reproducible**, and it applies
  **only to this one-shot path** — the cache (`compute_rs_from_array` /
  `build_mesh_from_cache`) is always faithful. Clean molecules never trigger a
  retry, so `Auto` is identical to `None` for them.
- **`vertex_normals`** — for each face accumulate `cross(b - a, c - a)` onto its
  three vertices, then normalize. This is the same normal the build and
  post-processing fill in; exposed for callers that edit the mesh themselves and
  need to refresh the normals afterwards.

### Carry rules for `atom_id` and `face_type`

| stage | `atom_id` | `face_type` |
|---|---|---|
| build (`build_surface_from_array` / `build_mesh_from_cache`) | **set** (1-based owner, `0` = unknown) | **set** (3/2/1) |
| `remove_flaps` | **carried** unchanged — vertices are left untouched, so ids stay aligned | **filtered** to the surviving faces (stays aligned with `faces`) |
| `close_cusps` | **empty** — the weld merges vertices of possibly different owners, so per-vertex ownership is not preserved | **empty** — the fan-fill rebuilds faces with no SES component |

`close_cusps` welds near-coincident boundary vertices (round(`v / weld_tol`)
integer-key buckets), then fan-fills the small open boundary cycles — the heavy
"fully closed" option. `remove_flaps` is a no-op on a clean mesh.

## 4. Data contract

The public API is fully **0-based and dense**: `faces` index into `verts` from
`0`, and `MeshReport`/diagnostics counts are plain sizes. The sole 1-based value
that crosses the boundary is `MeshResult.atom_id` — value `i` refers to the
caller's `xyzr[i-1]`, with `0` reserved for "unknown" (see §2).

Internally the library mirrors the MolSurfComp algorithm with **1-based indexing
and dummy padding rows/columns**; that convention never leaks through this
facade. See [`INTERNALS.md`](INTERNALS.md) for the internal record layouts.

## 5. Output formats (CLI-only, not part of the library API)

The MSMS `.vert`/`.face` pair, the ASCII PLY mesh, and the `.xyzr` atom input are
**file formats the `meshms_cli` reference front-end reads/writes — they are NOT
part of the library API.** The facade trades only in `std`-typed in-memory arrays
(`MeshResult` etc.); serialization to these on-disk formats lives entirely in the
CLI (`src/cli/`). They are documented here once so the formats are described in a
single place:

- **MSMS `.vert`/`.face`** (the CLI default; `<OUT>.vert` + `<OUT>.face`). Each
  vertex line is `x y z  nx ny nz  0  closest_sphere  vertex_type`; each face line
  is `v1 v2 v3 (1-based)  face_type  0`. The SES component lands in
  `closest_sphere` (= `atom_id`) and the type fields (= `face_type`).
- **ASCII PLY** (`--ply`, or an `OUT` ending in `.ply`). A single ASCII PLY mesh;
  `--vertex-normals` adds per-vertex normals.
- **`.xyzr`** input — one atom per line, `x y z radius` (blank/`#` rows skipped,
  first four columns taken), parsed like `numpy.loadtxt`.

## 6. Implementation note: the facade boundary copy

The facade deliberately performs a plain element-wise copy at the boundary, in
both directions, rather than handing out internal types:

- **out** — the internal `Surface` (`Vec3`/`Tri`) is copied field-by-field into a
  `MeshResult` of `std::array` rows.
- **in** — `MeshResult` arrays are copied back into internal `Vec3`/`Tri` for the
  post-processing entry points, the inverse of the out-copy.

This copy is what keeps the boundary C++17-clean: the internal types pull in
C++20 (`vec3.hpp` uses `<numbers>`), so they cannot appear in this header. The
copy is cold relative to meshing, and it is the price of the C++17/C++20 split —
no internal type ever crosses the API boundary.
