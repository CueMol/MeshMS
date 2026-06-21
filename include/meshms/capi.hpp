#pragma once
// Public, C++17-safe facade for the meshms_core SES mesher.
//
// This is the ONLY header an external consumer (e.g. cuemol2, BALL) needs to
// include. It deliberately includes NONE of the internal meshms headers (which
// require C++20: vec3.hpp uses <numbers>, etc.) -- only standard-library headers
// that are valid in C++17. A C++17 program can therefore include this header and
// link the (C++20-built) static library directly.
//
// Usage:
//   * One-shot:  build_surface_from_array(xyzr, Rp, mesh_size, fuse)
//   * Cached:    auto rs = compute_rs_from_array(xyzr, Rp);          // once
//                auto m0 = build_mesh_from_cache(rs, 0.5);           // density 1
//                auto m1 = build_mesh_from_cache(rs, 0.25);          // density 2
//     The cache holds the density-independent RS components; only the (cheap)
//     mesher re-runs per density.
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

}  // namespace meshms
