// Faithful port of the mesh-check oracle (mesh_area, signed_volume,
// manifold_report). Per-triangle float term/eval order matches the numpy helpers
// in vec3.hpp. The OUTER reduction is a sequential left-to-right sum rather than
// numpy's pairwise summation, so the totals drift from the Python by ~1e-16
// relative on the golden meshes -- far inside the equivalence gate (counts exact;
// area/volume within ~1e-6 relative). The C++ port is gated on
// equivalence, NOT byte-identity, so this reduction-order difference is expected.
//
// The core is templated over the vertex/face container element types so the
// facade-layout overloads (std::array<double,3> / std::array<uint32_t,3>, the
// meshms::capi MeshResult arrays) run the IDENTICAL math without converting the
// whole mesh to Vec3/Tri first. The array<double,3> -> Vec3 load is an
// element-wise copy of the same three doubles, so results are bit-identical.
#include "meshms/mesh_check.hpp"

// This TU reports non-finite vertices, so it needs a working isfinite().
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__
#error "mesh_check.cpp requires a working isfinite(): -ffinite-math-only is forbidden"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>

namespace meshms {

namespace {

// Element adapters: position of vertex i and 0-based index k of a face, for both
// the internal (Vec3/Tri) and the facade (array<double,3>/array<uint32,3>) types.
inline Vec3 vget(const std::vector<Vec3>& V, std::size_t i) { return V[i]; }
inline Vec3 vget(const std::vector<std::array<double, 3>>& V, std::size_t i) {
  return Vec3{V[i][0], V[i][1], V[i][2]};
}
inline int32_t iget(const Tri& f, int k) { return f[static_cast<std::size_t>(k)]; }
inline int32_t iget(const std::array<std::uint32_t, 3>& f, int k) {
  return static_cast<int32_t>(f[static_cast<std::size_t>(k)]);
}

// FNV-1a over the raw bytes of a sorted index triple (hash-map key hash only;
// key equality is exact, so collisions cost time, never correctness).
struct TriKeyHash {
  std::size_t operator()(const std::array<int32_t, 3>& t) const {
    std::uint64_t h = 1469598103934665603ULL;
    const unsigned char* p = reinterpret_cast<const unsigned char*>(t.data());
    for (std::size_t i = 0; i < sizeof(t); ++i) {
      h ^= p[i];
      h *= 1099511628211ULL;
    }
    return static_cast<std::size_t>(h);
  }
};

template <class VT, class FT>
double mesh_area_impl(const std::vector<VT>& V, const std::vector<FT>& F) {
  if (F.empty()) return 0.0;
  double total = 0.0;
  for (const FT& f : F) {
    const Vec3 a = vget(V, static_cast<std::size_t>(iget(f, 0)));
    const Vec3 b = vget(V, static_cast<std::size_t>(iget(f, 1)));
    const Vec3 c = vget(V, static_cast<std::size_t>(iget(f, 2)));
    Vec3 cr = cross(b - a, c - a);
    // np.sqrt((cr*cr).sum(axis=1)) summed: per-triangle magnitude then 0.5*sum.
    total += norm(cr);
  }
  return 0.5 * total;
}

template <class VT, class FT>
double signed_volume_impl(const std::vector<VT>& V, const std::vector<FT>& F) {
  if (F.empty()) return 0.0;
  double acc = 0.0;
  for (const FT& f : F) {
    const Vec3 a = vget(V, static_cast<std::size_t>(iget(f, 0)));
    const Vec3 b = vget(V, static_cast<std::size_t>(iget(f, 1)));
    const Vec3 c = vget(V, static_cast<std::size_t>(iget(f, 2)));
    // (a * np.cross(b, c)).sum() == a . (b x c) == triple(a, b, c).
    acc += triple(a, b, c);
  }
  return acc / 6.0;
}

template <class VT, class FT>
ManifoldReport manifold_report_impl(const std::vector<VT>& V,
                                    const std::vector<FT>& F) {
  ManifoldReport rep{};
  const int nV = static_cast<int>(V.size());
  const int nF = static_cast<int>(F.size());
  rep.n_vertices = nV;
  rep.n_faces = nF;

  // NaN/Inf tripwire. Cheap (one O(nV) pass) and always compiled: a relaxed-FP
  // deploy build has no golden to diff against, so this is what tells a consumer
  // the mesh is unusable rather than merely different.
  int nonfinite = 0;
  for (std::size_t i = 0; i < V.size(); ++i) {
    const Vec3 v = vget(V, i);
    if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z)) ++nonfinite;
  }
  rep.nonfinite_vertices = nonfinite;

  // Degenerate triangles: repeated vertex index OR |cross| < 1e-12.
  int degen = 0;
  for (const FT& f : F) {
    const int32_t i = iget(f, 0), j = iget(f, 1), k = iget(f, 2);
    bool repeated = (i == j) || (j == k) || (i == k);
    const Vec3 a = vget(V, static_cast<std::size_t>(i));
    const Vec3 b = vget(V, static_cast<std::size_t>(j));
    const Vec3 c = vget(V, static_cast<std::size_t>(k));
    Vec3 cr = cross(b - a, c - a);
    bool zero = norm(cr) < 1e-12;
    if (repeated || zero) ++degen;
  }
  rep.degenerate_faces = degen;

  // Undirected edge incidence (key: sorted (u,v) packed into a uint64) and
  // sorted-triple face counts. Hash maps instead of ordered maps: the consumers
  // below only COUNT entries by incidence, so iteration order is irrelevant and
  // the O(log n) tree walks / per-node allocations are pure overhead.
  std::unordered_map<std::uint64_t, int> edge_count;
  std::unordered_map<std::array<int32_t, 3>, int, TriKeyHash> face_seen;
  edge_count.reserve(F.size() * 3);
  face_seen.reserve(F.size());
  for (const FT& f : F) {
    const int32_t i = iget(f, 0), j = iget(f, 1), k = iget(f, 2);
    const std::array<std::pair<int32_t, int32_t>, 3> edges = {
        std::make_pair(i, j), std::make_pair(j, k), std::make_pair(k, i)};
    for (const auto& e : edges) {
      int32_t u = e.first, v = e.second;
      if (u > v) std::swap(u, v);
      const std::uint64_t key =
          (static_cast<std::uint64_t>(static_cast<std::uint32_t>(u)) << 32) |
          static_cast<std::uint32_t>(v);
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
  rep.area = mesh_area_impl(V, F);
  rep.signed_volume = signed_volume_impl(V, F);
  return rep;
}

}  // namespace

double mesh_area(const std::vector<Vec3>& V, const std::vector<Tri>& F) {
  return mesh_area_impl(V, F);
}

double signed_volume(const std::vector<Vec3>& V, const std::vector<Tri>& F) {
  return signed_volume_impl(V, F);
}

ManifoldReport manifold_report(const std::vector<Vec3>& V, const std::vector<Tri>& F) {
  return manifold_report_impl(V, F);
}

ManifoldReport manifold_report(const std::vector<std::array<double, 3>>& V,
                               const std::vector<std::array<std::uint32_t, 3>>& F) {
  return manifold_report_impl(V, F);
}

}  // namespace meshms
