// Implementation of the C++17-safe public facade (capi.hpp). Compiled as part of
// meshms_core (C++20); only the HEADER is constrained to C++17.
#include "meshms/capi.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "meshms/geom.hpp"
#include "meshms/mesh.hpp"
#include "meshms/mesh_check.hpp"
#include "meshms/pipeline.hpp"
#include "meshms/vec3.hpp"
#include "meshms/weld.hpp"

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
  r.face_type = s.ftype;  // per-face SES type, aligned with faces
  return r;
}

// MeshResult arrays -> internal Vec3/Tri (the inverse of to_result's copy).
std::vector<Vec3> to_vec3(const std::vector<std::array<double, 3>>& verts) {
  std::vector<Vec3> V;
  V.reserve(verts.size());
  for (const std::array<double, 3>& v : verts) V.push_back(Vec3{v[0], v[1], v[2]});
  return V;
}
std::vector<Tri> to_tri(const std::vector<std::array<std::uint32_t, 3>>& faces) {
  std::vector<Tri> F;
  F.reserve(faces.size());
  for (const std::array<std::uint32_t, 3>& f : faces)
    F.push_back(Tri{static_cast<std::int32_t>(f[0]), static_cast<std::int32_t>(f[1]),
                    static_cast<std::int32_t>(f[2])});
  return F;
}
std::vector<std::array<double, 3>> from_vec3(const std::vector<Vec3>& V) {
  std::vector<std::array<double, 3>> out;
  out.reserve(V.size());
  for (const Vec3& v : V) out.push_back({v.x, v.y, v.z});
  return out;
}
std::vector<std::array<std::uint32_t, 3>> from_tri(const std::vector<Tri>& F) {
  std::vector<std::array<std::uint32_t, 3>> out;
  out.reserve(F.size());
  for (const Tri& f : F)
    out.push_back({static_cast<std::uint32_t>(f[0]), static_cast<std::uint32_t>(f[1]),
                   static_cast<std::uint32_t>(f[2])});
  return out;
}

// Area-weighted per-vertex normals: cross(b-a, c-a) accumulated per vertex, then
// normalized. Term-for-term identical to ply.vertex_normals_from_faces so the
// CLI's --vertex-normals output is unchanged.
std::vector<Vec3> vnormals_from_faces(const std::vector<Vec3>& V,
                                      const std::vector<Tri>& F) {
  std::vector<Vec3> Nv(V.size(), Vec3{});
  for (const Tri& f : F) {
    const Vec3& a = V[static_cast<std::size_t>(f[0])];
    const Vec3& b = V[static_cast<std::size_t>(f[1])];
    const Vec3& c = V[static_cast<std::size_t>(f[2])];
    Vec3 fn = cross(b - a, c - a);  // magnitude == 2 * triangle area
    Nv[static_cast<std::size_t>(f[0])] = Nv[static_cast<std::size_t>(f[0])] + fn;
    Nv[static_cast<std::size_t>(f[1])] = Nv[static_cast<std::size_t>(f[1])] + fn;
    Nv[static_cast<std::size_t>(f[2])] = Nv[static_cast<std::size_t>(f[2])] + fn;
  }
  for (Vec3& nv : Nv) {
    double len = std::sqrt(pysq(nv.x) + pysq(nv.y) + pysq(nv.z));
    if (len == 0.0) len = 1.0;
    nv = nv / len;
  }
  return Nv;
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

// version() is defined in version.cpp (declared in capi.hpp).

MeshResult close_cusps(const MeshResult& mesh, double weld_tol) {
  WeldResult w = weld(to_vec3(mesh.verts), to_tri(mesh.faces), weld_tol);
  FillResult fh = fill_small_holes(std::move(w.V), w.F);  // w.V dead after this
  MeshResult r;
  r.verts = from_vec3(fh.V);
  r.faces = from_tri(fh.F);
  r.vnormals = from_vec3(vnormals_from_faces(fh.V, fh.F));
  // atom_id intentionally left empty: the weld merges vertices of different owners.
  return r;
}

MeshResult remove_flaps(const MeshResult& mesh, int passes) {
  std::vector<std::uint32_t> kept;
  FlapResult fl =
      remove_nonmanifold_flaps(to_vec3(mesh.verts), to_tri(mesh.faces), passes, &kept);
  MeshResult r;
  r.verts = from_vec3(fl.V);
  r.faces = from_tri(fl.F);
  r.vnormals = from_vec3(vnormals_from_faces(fl.V, fl.F));
  r.atom_id = mesh.atom_id;  // remove_nonmanifold_flaps keeps V intact -> ids align
  // Filter face_type to the surviving faces (kept holds original indices into faces).
  if (!mesh.face_type.empty()) {
    r.face_type.reserve(kept.size());
    for (std::uint32_t k : kept) r.face_type.push_back(mesh.face_type[k]);
  }
  return r;
}

std::vector<std::array<double, 3>> vertex_normals(
    const std::vector<std::array<double, 3>>& verts,
    const std::vector<std::array<std::uint32_t, 3>>& faces) {
  return from_vec3(vnormals_from_faces(to_vec3(verts), to_tri(faces)));
}

MeshReport analyze_mesh(const MeshResult& mesh) {
  // Facade-layout overload: no whole-mesh Vec3/Tri conversion copies.
  ManifoldReport rep = manifold_report(mesh.verts, mesh.faces);
  MeshReport out;
  out.n_vertices = static_cast<std::uint32_t>(rep.n_vertices);
  out.n_faces = static_cast<std::uint32_t>(rep.n_faces);
  out.degenerate_faces = static_cast<std::uint32_t>(rep.degenerate_faces);
  out.duplicate_faces = static_cast<std::uint32_t>(rep.duplicate_faces);
  out.boundary_edges = static_cast<std::uint32_t>(rep.boundary_edges);
  out.nonmanifold_edges = static_cast<std::uint32_t>(rep.nonmanifold_edges);
  out.watertight = rep.watertight;
  out.area = rep.area;
  out.signed_volume = rep.signed_volume;
  return out;
}

BoundaryDiagnostics boundary_diagnostics(const MeshResult& mesh) {
  // Facade-layout overload: no whole-mesh Vec3/Tri conversion copies.
  BoundaryLoopsResult res = boundary_loops(mesh.verts, mesh.faces);
  BoundaryDiagnostics out;
  out.loops.reserve(res.loops.size());
  for (const BoundaryLoop& L : res.loops) {
    BoundaryLoopInfo info;
    info.n_verts = static_cast<std::uint32_t>(L.n_verts);
    info.n_edges = static_cast<std::uint32_t>(L.n_edges);
    info.closed = L.closed;
    info.centroid = {L.centroid.x, L.centroid.y, L.centroid.z};
    out.loops.push_back(info);
  }
  out.nonmanifold.reserve(res.nonmanifold.size());
  for (const NonmanifoldEdge& e : res.nonmanifold) {
    NonmanifoldEdgeInfo info;
    info.u = static_cast<std::uint32_t>(e.u);
    info.v = static_cast<std::uint32_t>(e.v);
    info.count = static_cast<std::uint32_t>(e.count);
    out.nonmanifold.push_back(info);
  }
  return out;
}

}  // namespace meshms
