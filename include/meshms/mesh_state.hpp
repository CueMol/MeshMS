#pragma once
// Global mesh accumulator (faithful port of the mesh_state module): the running
// vertex / triangle / per-face-normal / per-vertex-tag lists that every patch
// mesher (mesh_sphpat, visutorpat, ...) pushes into via add_patch.
//
// Each patch mesher builds a local 1-based point list P (P[0] dummy, real points
// at P[1..Np]) and a 1-based triangle list T of index triples into P, plus an
// aligned per-face normal list. add_patch rebases the triangles onto the growing
// global vertex list (0-based output, like MeshState.to_arrays()).
#include <array>
#include <cstdint>
#include <vector>

#include "meshms/mesh.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

// Pick the nearest of `n` candidate atoms to point `p`. `a[k]` is the (1-based)
// atom id of candidate k and `c[k]` is its centre. Returns a[k] minimising the
// SQUARED distance; ties resolve to the smallest atom id (a[k]) so the result is
// fully deterministic and thread-count independent. Shared by the meshers that
// attribute a vertex to a neighbouring atom (concave). n must be >= 1.
inline int32_t nearest_atom(const Vec3& p, const int32_t* a, const Vec3* c, int n) {
  int best = 0;
  double best_d = dot(p - c[0], p - c[0]);
  for (int k = 1; k < n; ++k) {
    double dk = dot(p - c[static_cast<std::size_t>(k)], p - c[static_cast<std::size_t>(k)]);
    if (dk < best_d || (dk == best_d && a[k] < a[best])) {
      best = k;
      best_d = dk;
    }
  }
  return a[best];
}

// Per-vertex boundary tag for ID-based fusion (meshms.fusion.fuse_by_id). The
// Python uses None / ('atom',a) / ('probe',p) / ('tseam',k,row) tuples; here a
// flat struct: kind 0=none, 1=atom, 2=probe, 3=tseam. For atom/probe only i is
// used; for tseam i=k and j=row.
struct Tag {
  int kind = 0;
  int i = 0;
  int j = 0;
};

inline bool operator==(const Tag& a, const Tag& b) {
  return a.kind == b.kind && a.i == b.i && a.j == b.j;
}
inline bool operator!=(const Tag& a, const Tag& b) { return !(a == b); }

// Per-vertex tag LIST: a single vertex may carry several tags at once (e.g. a
// toroidal corner lying on BOTH an atom sphere and a probe). An empty TagList
// means interior (never fused), mirroring the Python None / [] entries.
using TagList = std::vector<Tag>;

struct MeshState {
  std::vector<Vec3> V;        // global vertices, 0-based, raw coords
  std::vector<Tri> F;         // global triangles, 0-based indices
  std::vector<Vec3> N;        // per-face outward normals
  std::vector<TagList> tags;  // per-vertex boundary tag LIST (empty => interior)
  std::vector<int32_t> vatom; // per-vertex owning atom (1-based, 0 => unknown)

  // Append one patch mesh.
  //   P            : P[0] dummy, P[1..Np] the real points (size Np+1).
  //   T            : 1-based (a,b,c) index triples into P.
  //   face_normals : per-face normals aligned with T (empty => none, like the
  //                  Python face_normals is None: no normals are appended).
  //   vids         : per-vertex tag lists aligned with P[1..Np] (empty whole =>
  //                  all interior; a per-vertex empty TagList is also interior),
  //                  like the Python vids=None default / None entries.
  //   patch_vatom  : per-vertex owning atom aligned with P[1..Np] (empty whole =>
  //                  all 0/unknown). Parallel to vids; never affects V/F/N.
  void add_patch(const std::vector<Vec3>& P,
                 const std::vector<std::array<int, 3>>& T,
                 const std::vector<Vec3>& face_normals,
                 const std::vector<TagList>& vids = {},
                 const std::vector<int32_t>& patch_vatom = {}) {
    const std::size_t base = V.size();
    // P[1:] are the real points; P[0] is the dummy slot.
    const std::size_t Np = P.empty() ? 0 : P.size() - 1;
    for (std::size_t k = 1; k <= Np; ++k) {
      V.push_back(P[k]);
      tags.push_back(vids.empty() ? TagList{} : vids[k - 1]);
      vatom.push_back(patch_vatom.empty() ? 0 : patch_vatom[k - 1]);
    }
    const bool have_nrm = !face_normals.empty();
    for (std::size_t idx = 0; idx < T.size(); ++idx) {
      const auto& t = T[idx];
      F.push_back(Tri{static_cast<int32_t>(t[0] - 1 + static_cast<int>(base)),
                      static_cast<int32_t>(t[1] - 1 + static_cast<int>(base)),
                      static_cast<int32_t>(t[2] - 1 + static_cast<int>(base))});
      if (have_nrm) N.push_back(face_normals[idx]);
    }
  }

  std::size_t nverts() const { return V.size(); }
  std::size_t ntris() const { return F.size(); }
};

}  // namespace meshms
