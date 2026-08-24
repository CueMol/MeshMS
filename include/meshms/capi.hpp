#pragma once
// Public, C++17-safe facade for the MeshMS (libMeshMS) SES surface mesher.
//
// This is the ONLY header an external consumer needs to include; cuemol2/cuemol3
// is the primary consumer. It deliberately includes NONE of the internal meshms
// headers (which require C++20: vec3.hpp uses <numbers>, etc.) -- only standard-
// library headers that are valid in C++17. A C++17 program can therefore include
// this header and link the (C++20-built) libMeshMS static library directly.
//
// The whole API lives in namespace meshms and trades only in std types
// (std::vector / std::array / std::uint32_t), so no MeshMS-internal type ever
// crosses the boundary. Beyond meshing, it exposes the post-processing
// (close_cusps, remove_flaps), the diagnostics (analyze_mesh,
// boundary_diagnostics) and the normal recompute (vertex_normals) a consumer
// needs to finish and validate a surface.
//
// Usage:
//   * One-shot:  build_surface_from_array(xyzr, Rp, mesh_size, fuse)
//   * Cached:    auto rs = compute_rs_from_array(xyzr, Rp);          // once
//                auto m0 = build_mesh_from_cache(rs, 0.5);           // density 1
//                auto m1 = build_mesh_from_cache(rs, 0.25);          // density 2
//     The cache holds the density-independent RS components; only the (cheap)
//     mesher re-runs per density.
//
// Multi-component inputs are supported at this layer: the atoms are split into
// connected components of the SAS-intersection graph, the pipeline runs per
// component, and the meshes are concatenated (components in ascending order of
// their smallest atom index). An isolated atom (no SAS neighbour) is meshed
// directly as its full vdW sphere (icosphere at the requested mesh_size,
// face_type 3). atom_id always refers to the caller's input array. An empty
// input yields an empty mesh.
#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace meshms {

// Opaque handle to the precomputed, density-independent RS components. The
// definition lives in capi.cpp; consumers only ever hold a shared_ptr to it.
struct RSCache;

// Triangle mesh returned to the caller. All arrays are 0-based and dense.
struct MeshResult {
  std::vector<std::array<double, 3>> verts;     // vertex positions
  std::vector<std::array<double, 3>> vnormals;  // per-vertex outward normals (aligned with verts)
  std::vector<std::array<std::uint32_t, 3>> faces;  // triangles, 0-based indices into verts
  // Per-vertex owning atom, aligned with verts. Value i (1-based, 0 = unknown)
  // refers to xyzr[i-1] -- the (i-1)-th entry of the input array the caller
  // passed. The caller maps this index back to its own atom id.
  std::vector<std::uint32_t> atom_id;
  // Per-face SES component type, aligned with faces (MSMS face-type codes):
  //   3 = convex (contact), 2 = concave (spherical reentrant),
  //   1 = toroidal (toric reentrant).
  // Populated by build_*; carried (filtered) through remove_flaps. Empty after
  // close_cusps, which rebuilds faces (fan-fill), so the types no longer apply.
  std::vector<std::uint8_t> face_type;
};

// Compute the density-independent RS components from an atom array. xyzr[a] is
// {x, y, z, radius} for atom a (0-based input; output atom_id value i refers to
// xyzr[i-1]). `radius_probe` is the solvent probe radius Rp.
std::shared_ptr<RSCache> compute_rs_from_array(
    const std::vector<std::array<double, 4>>& xyzr, double radius_probe);

// Mesh the SES from a precomputed cache at the given mesh density. fuse=true
// welds tagged boundary vertices and the result then carries per-vertex normals
// recomputed from the fused geometry.
MeshResult build_mesh_from_cache(const std::shared_ptr<RSCache>& rs,
                                 double mesh_size, bool fuse = false);

// One-shot convenience: compute_rs_from_array + build_mesh_from_cache. Bit-for-
// bit identical (verts/faces) to building from the equivalent .xyzr file.
MeshResult build_surface_from_array(
    const std::vector<std::array<double, 4>>& xyzr, double radius_probe,
    double mesh_size, bool fuse = false);

// Library version string, e.g. "0.1.0". Safe to call at any time.
const char* version();

// Human-readable build configuration: version, FP policy, parallel backend.
// Unlike version() -- which is a semver string embedded in the MSMS .vert/.face
// headers -- this is free-form and meant for logs and bug reports. A surface
// produced by a relaxed-FP deploy build is identifiable from it.
const char* build_info();

// FP policy the LIBRARY was compiled with: 0 = strict (bit-exact, golden-gated),
// 1 = fast (deploy). Compare against MESHMS_FP_FAST to detect a header/library
// mismatch, which would mean the inline math in the internal headers was
// compiled under a different policy than the library.
int fp_mode();

// ===== Mesh post-processing =================================================
// Each takes a MeshResult and returns a new one. vnormals in the result are
// recomputed (area-weighted) from the returned geometry, so they stay valid.

// Close cusp/singular seams into a watertight 2-manifold: weld near-coincident
// boundary vertices (round(v / weld_tol) integer-key buckets), then fan-fill the
// small open boundary cycles. This is the heavy "fully closed" option.
// NOTE: welding merges vertices of possibly different owners, so the returned
// atom_id is EMPTY -- per-vertex atom ownership is not preserved across a weld.
// face_type is likewise EMPTY: the fan-fill adds faces with no SES component.
MeshResult close_cusps(const MeshResult& mesh, double weld_tol = 1e-4);

// Drop spurious doubled-"flap" (non-manifold) triangles; a no-op on a clean
// mesh. Vertices are left untouched, so atom_id is carried through unchanged;
// face_type is filtered to the surviving faces (stays aligned with faces).
MeshResult remove_flaps(const MeshResult& mesh, int passes = 4);

// Area-weighted per-vertex outward normals from (verts, faces): for each face
// accumulate cross(b - a, c - a) onto its three vertices, then normalize. This
// is the same normal the post-processing above fills in; exposed for callers
// that edit the mesh themselves and need to refresh the normals afterwards.
std::vector<std::array<double, 3>> vertex_normals(
    const std::vector<std::array<double, 3>>& verts,
    const std::vector<std::array<std::uint32_t, 3>>& faces);

// ===== Mesh diagnostics =====================================================

// Whole-mesh quality report: edge incidence, watertightness, area and volume.
struct MeshReport {
  std::uint32_t n_vertices = 0;
  std::uint32_t n_faces = 0;
  // Vertices with a NaN or Inf component. A build with a relaxed FP policy
  // (MESHMS_FP=fast) has no bit-exact reference to diff against, so a consumer
  // should check this once after building a surface: 0 means usable.
  std::uint32_t nonfinite_vertices = 0;
  std::uint32_t degenerate_faces = 0;   // repeated index or |cross| < 1e-12
  std::uint32_t duplicate_faces = 0;    // same vertex triple seen more than once
  std::uint32_t boundary_edges = 0;     // edges used by exactly one face
  std::uint32_t nonmanifold_edges = 0;  // edges used by >= 3 faces
  bool watertight = false;              // boundary == 0 && nonmanifold == 0 && nF > 0
  double area = 0.0;                    // sum of triangle areas
  double signed_volume = 0.0;           // divergence-theorem signed volume
};
MeshReport analyze_mesh(const MeshResult& mesh);

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

// Localize open boundaries/holes and non-manifold edges (the CLI -v diagnostics).
// loops are sorted largest-first by n_edges; nonmanifold keeps first-seen order.
struct BoundaryDiagnostics {
  std::vector<BoundaryLoopInfo> loops;
  std::vector<NonmanifoldEdgeInfo> nonmanifold;
};
BoundaryDiagnostics boundary_diagnostics(const MeshResult& mesh);

}  // namespace meshms
