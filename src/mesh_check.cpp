// Faithful port of tests/check_mesh.py mesh oracle (mesh_area, signed_volume,
// manifold_report). Per-triangle float term/eval order matches the numpy helpers
// in vec3.hpp. The OUTER reduction is a sequential left-to-right sum rather than
// numpy's pairwise summation, so the totals drift from the Python by ~1e-16
// relative on the golden meshes -- far inside the equivalence gate (counts exact;
// area/volume within ~1e-6 relative, CPP_PORT_PLAN §5). The C++ port is gated on
// equivalence, NOT byte-identity, so this reduction-order difference is expected.
#include "meshms/mesh_check.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>

namespace meshms {

namespace {

// Vertex i of a triangle (Tri stores 0-based int32 indices).
inline const Vec3& vert(const std::vector<Vec3>& V, int32_t i) {
  return V[static_cast<std::size_t>(i)];
}

}  // namespace

double mesh_area(const std::vector<Vec3>& V, const std::vector<Tri>& F) {
  if (F.empty()) return 0.0;
  double total = 0.0;
  for (const Tri& f : F) {
    const Vec3& a = vert(V, f[0]);
    const Vec3& b = vert(V, f[1]);
    const Vec3& c = vert(V, f[2]);
    Vec3 cr = cross(b - a, c - a);
    // np.sqrt((cr*cr).sum(axis=1)) summed: per-triangle magnitude then 0.5*sum.
    total += norm(cr);
  }
  return 0.5 * total;
}

double signed_volume(const std::vector<Vec3>& V, const std::vector<Tri>& F) {
  if (F.empty()) return 0.0;
  double acc = 0.0;
  for (const Tri& f : F) {
    const Vec3& a = vert(V, f[0]);
    const Vec3& b = vert(V, f[1]);
    const Vec3& c = vert(V, f[2]);
    // (a * np.cross(b, c)).sum() == a . (b x c) == triple(a, b, c).
    acc += triple(a, b, c);
  }
  return acc / 6.0;
}

ManifoldReport manifold_report(const std::vector<Vec3>& V, const std::vector<Tri>& F) {
  ManifoldReport rep{};
  const int nV = static_cast<int>(V.size());
  const int nF = static_cast<int>(F.size());
  rep.n_vertices = nV;
  rep.n_faces = nF;

  // Degenerate triangles: repeated vertex index OR |cross| < 1e-12.
  int degen = 0;
  for (const Tri& f : F) {
    bool repeated = (f[0] == f[1]) || (f[1] == f[2]) || (f[0] == f[2]);
    const Vec3& a = vert(V, f[0]);
    const Vec3& b = vert(V, f[1]);
    const Vec3& c = vert(V, f[2]);
    Vec3 cr = cross(b - a, c - a);
    bool zero = norm(cr) < 1e-12;
    if (repeated || zero) ++degen;
  }
  rep.degenerate_faces = degen;

  // Undirected edge incidence (key: sorted (u,v)) and sorted-triple face counts.
  std::map<std::pair<int32_t, int32_t>, int> edge_count;
  std::map<std::array<int32_t, 3>, int> face_seen;
  for (const Tri& f : F) {
    int32_t i = f[0], j = f[1], k = f[2];
    const std::array<std::pair<int32_t, int32_t>, 3> edges = {
        std::make_pair(i, j), std::make_pair(j, k), std::make_pair(k, i)};
    for (const auto& e : edges) {
      int32_t u = e.first, v = e.second;
      std::pair<int32_t, int32_t> key = (u < v) ? std::make_pair(u, v) : std::make_pair(v, u);
      ++edge_count[key];
    }
    std::array<int32_t, 3> tri = {i, j, k};
    std::sort(tri.begin(), tri.end());
    ++face_seen[tri];
  }

  int boundary = 0;
  int nonmanifold = 0;
  for (const auto& kv : edge_count) {
    int n = kv.second;
    if (n == 1) ++boundary;
    else if (n >= 3) ++nonmanifold;
  }
  int dup = 0;
  for (const auto& kv : face_seen) {
    if (kv.second > 1) ++dup;
  }

  rep.boundary_edges = boundary;
  rep.nonmanifold_edges = nonmanifold;
  rep.duplicate_faces = dup;
  rep.watertight = (boundary == 0) && (nonmanifold == 0) && (nF > 0);
  rep.area = mesh_area(V, F);
  rep.signed_volume = signed_volume(V, F);
  return rep;
}

}  // namespace meshms
