#pragma once
// Vertex welding + boundary/flap repair --- faithful port of the weld module.
//
// The patch meshers emit each spherical/toroidal patch independently, so
// adjacent patches duplicate vertices along shared boundary arcs. These
// post-processing routines make the accumulated mesh a clean, weldable surface:
//   weld                   coordinate weld (round(V/tol) integer-key buckets)
//   boundary_loops         localize holes (boundary-edge components) + nonmanifold
//   fill_small_holes       fan-triangulate tiny open boundary cycles
//   remove_nonmanifold_flaps  drop doubled "flap" face fans hanging off an edge
//
// manifold_report is ALREADY ported in mesh_check.hpp (different dict keys); it is
// NOT duplicated here. Indexing here is the 0-based PLY/mesh convention.
#include <cstdint>
#include <vector>

#include "meshms/mesh.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

// weld(V, F, tol): merge near-coincident vertices via a round(V/tol) integer-key
// spatial hash (numpy np.round == banker's rounding, round-half-to-even), remap
// the faces onto the welded vertices, and drop degenerate triangles (any two
// welded indices coincide). Returns the welded (V2, F2).
//
// The first vertex hashing to a bucket is its representative coordinate, and the
// new-index order follows first-seen order (range(n)), exactly like the Python.
struct WeldResult {
  std::vector<Vec3> V;
  std::vector<Tri> F;
};
WeldResult weld(const std::vector<Vec3>& V, const std::vector<Tri>& F, double tol = 1e-6);

// One connected component of the boundary-edge subgraph = one hole/seam.
struct BoundaryLoop {
  int n_verts = 0;          // #vertices in the component
  int n_edges = 0;          // #boundary edges in the component
  Vec3 centroid{};          // mean of the component's vertex coordinates
  bool closed = false;      // every member has exactly 2 boundary neighbours
  std::vector<int> vids;    // member vertex ids, sorted ascending
};

// A non-manifold edge (u, v) used by `count` >= 3 faces.
struct NonmanifoldEdge {
  int u = 0;
  int v = 0;
  int count = 0;
};

// boundary_loops(V, F): group the boundary edges (used by exactly one face) into
// connected components (each one hole/seam) and collect non-manifold edges (used
// by 3+ faces). loops are sorted largest-first by n_edges (stable, matching the
// Python list.sort); nonmanifold preserves first-seen edge order.
struct BoundaryLoopsResult {
  std::vector<BoundaryLoop> loops;
  std::vector<NonmanifoldEdge> nonmanifold;
};
BoundaryLoopsResult boundary_loops(const std::vector<Vec3>& V, const std::vector<Tri>& F);

// fill_small_holes(V, F, max_loop): follow boundary edges as directed half-edges
// into closed cycles and fan-triangulate every cycle of at most max_loop edges,
// orienting each new triangle to agree with the adjacent-face normal average.
// Larger genuine boundaries are left untouched. V is returned unchanged; only F
// gains the fill triangles (appended after the original faces).
struct FillResult {
  std::vector<Vec3> V;
  std::vector<Tri> F;
};
FillResult fill_small_holes(const std::vector<Vec3>& V, const std::vector<Tri>& F,
                            int max_loop = 12);

// remove_nonmanifold_flaps(V, F, passes): topological flap removal. On a
// non-manifold edge (>=3 incident faces) an incident triangle whose apex hangs
// ONLY off that edge (every face at the apex contains the edge) is a flap; the
// apex's whole face fan is dropped. Repeats until no non-manifold edge or no
// flap remains. Not a coordinate weld; removes nothing on a clean mesh.
struct FlapResult {
  std::vector<Vec3> V;
  std::vector<Tri> F;
};
FlapResult remove_nonmanifold_flaps(const std::vector<Vec3>& V, const std::vector<Tri>& F,
                                    int passes = 4);

}  // namespace meshms
