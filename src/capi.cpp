// Implementation of the C++17-safe public facade (capi.hpp). Compiled as part of
// meshms_core (C++20); only the HEADER is constrained to C++17.
#include "meshms/capi.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>

#include "meshms/geom.hpp"
#include "meshms/mesh.hpp"
#include "meshms/pipeline.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

// Opaque cache definition: just wraps the internal RSComponents.
struct RSCache {
  RSComponents rs;
};

namespace {

// Build a Geom from the xyzr array, matching read_xyzr's row order EXACTLY
// (geom.cpp: 1-based atoms 1..M, dummy row 0, centers={x,y,z}, R=radius) so that
// from-array == from-file bit-for-bit.
Geom geom_from_array(const std::vector<std::array<double, 4>>& xyzr) {
  const int m = static_cast<int>(xyzr.size());
  Geom g;
  g.M = m;
  g.centers.assign(static_cast<std::size_t>(m) + 1, Vec3{0.0, 0.0, 0.0});  // row 0 dummy
  g.R.assign(static_cast<std::size_t>(m) + 1, 0.0);                        // entry 0 dummy
  for (int a = 1; a <= m; ++a) {
    const std::array<double, 4>& c = xyzr[static_cast<std::size_t>(a - 1)];
    g.centers[static_cast<std::size_t>(a)] = Vec3{c[0], c[1], c[2]};
    g.R[static_cast<std::size_t>(a)] = c[3];
  }
  return g;
}

// Surface (internal) -> MeshResult (facade): plain element-wise copy.
MeshResult to_result(const Surface& s) {
  MeshResult r;
  r.verts.reserve(s.V.size());
  for (const Vec3& v : s.V) r.verts.push_back({v.x, v.y, v.z});
  r.vnormals.reserve(s.NV.size());
  for (const Vec3& n : s.NV) r.vnormals.push_back({n.x, n.y, n.z});
  r.faces.reserve(s.F.size());
  for (const Tri& f : s.F) {
    r.faces.push_back({static_cast<std::uint32_t>(f[0]),
                       static_cast<std::uint32_t>(f[1]),
                       static_cast<std::uint32_t>(f[2])});
  }
  r.atom_id.reserve(s.atom_id.size());
  for (std::int32_t a : s.atom_id) r.atom_id.push_back(static_cast<std::uint32_t>(a));
  return r;
}

}  // namespace

std::shared_ptr<RSCache> compute_rs_from_array(
    const std::vector<std::array<double, 4>>& xyzr, double radius_probe) {
  auto cache = std::make_shared<RSCache>();
  cache->rs = compute_rs(geom_from_array(xyzr), radius_probe);
  return cache;
}

MeshResult build_mesh_from_cache(const std::shared_ptr<RSCache>& rs,
                                 double mesh_size, bool fuse) {
  return to_result(build_mesh(rs->rs, mesh_size, fuse));
}

MeshResult build_surface_from_array(
    const std::vector<std::array<double, 4>>& xyzr, double radius_probe,
    double mesh_size, bool fuse) {
  RSComponents rs = compute_rs(geom_from_array(xyzr), radius_probe);
  return to_result(build_mesh(rs, mesh_size, fuse));
}

}  // namespace meshms
