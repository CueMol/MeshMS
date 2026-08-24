// Implementation of the C++17-safe public facade (capi.hpp). Compiled as part of
// meshms_core (C++20); only the HEADER is constrained to C++17.
//
// Multi-component support lives HERE, above the faithful pipeline: the input is
// split into connected components of the SAS-intersection graph (the same
// dist < (R_a+Rp)+(R_b+Rp) predicate interstructure() uses), the unchanged
// pipeline runs per component, and the meshes are concatenated. Isolated atoms
// (no SAS neighbour), which the faithful exterior extraction cannot represent
// (it throws "isolated SAS-ball" or silently drops unreachable components), are
// meshed directly as full vdW spheres. A single-component input takes the exact
// pre-existing path, so the golden gates are unaffected.
#include "meshms/capi.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

#include "meshms/geom.hpp"
#include "meshms/intersection.hpp"
#include "meshms/mesh.hpp"
#include "meshms/mesh_check.hpp"
#include "meshms/pipeline.hpp"
#include "meshms/vec3.hpp"
#include "meshms/weld.hpp"

namespace meshms {

// Opaque cache definition: the density-independent RS components of every
// connected component (>= 2 atoms each), plus the isolated atoms meshed as
// full spheres. `orig` maps a component-local 1-based atom index back to the
// 0-based index into the caller's xyzr array, so the merged mesh's atom_id
// keeps referring to the caller's array.
struct RSCache {
  struct Component {
    RSComponents rs;
    std::vector<std::uint32_t> orig;  // comp-local atom a -> orig[a-1] (0-based)
  };
  struct Isolated {
    std::uint32_t orig;          // 0-based index into the caller's xyzr
    std::array<double, 4> xyzr;  // {x, y, z, radius}
  };
  std::vector<Component> comps;
  std::vector<Isolated> isolated;
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
    // Reject non-finite or negative input here, at the facade boundary: a NaN
    // coordinate from a malformed structure would otherwise propagate silently
    // through every stage and surface as a NaN mesh. read_xyzr stays untouched
    // -- it is the faithful reader for the golden fixtures.
    // Radius 0 IS valid and stays accepted: the xyzr convention uses it for
    // atoms that do not contribute to the surface (barstar.xyzr has 503).
    if (!std::isfinite(c[0]) || !std::isfinite(c[1]) || !std::isfinite(c[2]) ||
        !std::isfinite(c[3]) || c[3] < 0.0) {
      throw std::invalid_argument(
          "meshms: atom " + std::to_string(a - 1) +
          " has a non-finite coordinate or a negative radius");
    }
    g.centers[static_cast<std::size_t>(a)] = Vec3{c[0], c[1], c[2]};
    g.R[static_cast<std::size_t>(a)] = c[3];
  }
  return g;
}

// MeshResult arrays -> internal Vec3/Tri (the inverse of append_surface's copy).
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

// Connected components of the SAS-intersection graph, via interstructure()
// (the exact predicate the pipeline itself uses) + union-find. Returns one
// ascending 0-based index list per component, components ordered by their
// smallest member -- fully deterministic.
std::vector<std::vector<std::uint32_t>> sas_components(
    const std::vector<std::array<double, 4>>& xyzr, double radius_probe) {
  const int m = static_cast<int>(xyzr.size());
  std::vector<int> parent(static_cast<std::size_t>(m) + 1);
  for (int a = 0; a <= m; ++a) parent[static_cast<std::size_t>(a)] = a;
  auto find = [&parent](int a) {
    while (parent[static_cast<std::size_t>(a)] != a) {
      parent[static_cast<std::size_t>(a)] =
          parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(a)])];
      a = parent[static_cast<std::size_t>(a)];
    }
    return a;
  };

  const Neighbors nb = interstructure(geom_from_array(xyzr), radius_probe);
  for (int a = 1; a <= m; ++a) {
    for (std::int32_t b : nb.of(a)) {
      const int ra = find(a), rb = find(static_cast<int>(b));
      if (ra != rb) parent[static_cast<std::size_t>(std::max(ra, rb))] = std::min(ra, rb);
    }
  }

  // Group 1-based atoms by root; map roots to components in first-seen
  // (= ascending smallest-member) order.
  std::vector<std::vector<std::uint32_t>> comps;
  std::vector<int> comp_of_root(static_cast<std::size_t>(m) + 1, -1);
  for (int a = 1; a <= m; ++a) {
    const int r = find(a);
    if (comp_of_root[static_cast<std::size_t>(r)] < 0) {
      comp_of_root[static_cast<std::size_t>(r)] = static_cast<int>(comps.size());
      comps.emplace_back();
    }
    comps[static_cast<std::size_t>(comp_of_root[static_cast<std::size_t>(r)])]
        .push_back(static_cast<std::uint32_t>(a - 1));
  }
  return comps;
}

// Append one component's Surface to the merged MeshResult: vertex indices are
// offset by the current vertex count and atom_id is remapped through `orig`
// back to the caller's array. With base == 0 and an identity `orig` (the
// single-component case) every emitted value equals to_result's -- bit-for-bit.
void append_surface(MeshResult& out, const Surface& s,
                    const std::vector<std::uint32_t>& orig) {
  const std::uint32_t base = static_cast<std::uint32_t>(out.verts.size());
  out.verts.reserve(out.verts.size() + s.V.size());
  for (const Vec3& v : s.V) out.verts.push_back({v.x, v.y, v.z});
  out.vnormals.reserve(out.vnormals.size() + s.NV.size());
  for (const Vec3& n : s.NV) out.vnormals.push_back({n.x, n.y, n.z});
  out.faces.reserve(out.faces.size() + s.F.size());
  for (const Tri& f : s.F) {
    out.faces.push_back({base + static_cast<std::uint32_t>(f[0]),
                         base + static_cast<std::uint32_t>(f[1]),
                         base + static_cast<std::uint32_t>(f[2])});
  }
  out.atom_id.reserve(out.atom_id.size() + s.atom_id.size());
  for (std::int32_t a : s.atom_id) {
    out.atom_id.push_back(
        a == 0 ? 0u : orig[static_cast<std::size_t>(a - 1)] + 1u);
  }
  out.face_type.insert(out.face_type.end(), s.ftype.begin(), s.ftype.end());
}

// Append a full sphere mesh (icosphere) for an isolated atom: its SES is its
// vdW sphere. The subdivision level is chosen so the edge length is <=
// mesh_size; normals are the exact radial directions; faces wind outward
// (positive signed volume); face_type is 3 (convex) for every face.
void append_isolated_sphere(MeshResult& out, const RSCache::Isolated& iso,
                            double mesh_size) {
  const double r = iso.xyzr[3];
  if (!(r > 0.0)) return;  // degenerate radius: nothing to mesh

  // Unit icosahedron (outward-wound), t = golden ratio.
  const double t = (1.0 + std::sqrt(5.0)) / 2.0;
  std::vector<Vec3> dirs = {
      {-1, t, 0}, {1, t, 0},   {-1, -t, 0}, {1, -t, 0},
      {0, -1, t}, {0, 1, t},   {0, -1, -t}, {0, 1, -t},
      {t, 0, -1}, {t, 0, 1},   {-t, 0, -1}, {-t, 0, 1}};
  for (Vec3& d : dirs) d = d / std::sqrt(pysq(d.x) + pysq(d.y) + pysq(d.z));
  std::vector<Tri> faces = {
      {0, 11, 5}, {0, 5, 1},  {0, 1, 7},   {0, 7, 10}, {0, 10, 11},
      {1, 5, 9},  {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
      {3, 9, 4},  {3, 4, 2},  {3, 2, 6},   {3, 6, 8},  {3, 8, 9},
      {4, 9, 5},  {2, 4, 11}, {6, 2, 10},  {8, 6, 7},  {9, 8, 1}};

  // Unit icosahedron edge length is 4/sqrt(10+2*sqrt(5)) ~= 1.05146; each
  // subdivision halves it. Cap the level so a tiny mesh_size cannot explode.
  const double edge0 = 4.0 / std::sqrt(10.0 + 2.0 * std::sqrt(5.0));
  int level = 0;
  while (level < 7 && edge0 * r / static_cast<double>(1 << level) > mesh_size)
    ++level;

  for (int s = 0; s < level; ++s) {
    std::map<std::pair<std::int32_t, std::int32_t>, std::int32_t> midpoint;
    auto mid = [&](std::int32_t a, std::int32_t b) {
      const std::pair<std::int32_t, std::int32_t> key = std::minmax(a, b);
      auto it = midpoint.find(key);
      if (it != midpoint.end()) return it->second;
      Vec3 mvec = dirs[static_cast<std::size_t>(a)] + dirs[static_cast<std::size_t>(b)];
      mvec = mvec / std::sqrt(pysq(mvec.x) + pysq(mvec.y) + pysq(mvec.z));
      const std::int32_t idx = static_cast<std::int32_t>(dirs.size());
      dirs.push_back(mvec);
      midpoint.emplace(key, idx);
      return idx;
    };
    std::vector<Tri> next;
    next.reserve(faces.size() * 4);
    for (const Tri& f : faces) {
      const std::int32_t ab = mid(f[0], f[1]);
      const std::int32_t bc = mid(f[1], f[2]);
      const std::int32_t ca = mid(f[2], f[0]);
      next.push_back({f[0], ab, ca});
      next.push_back({f[1], bc, ab});
      next.push_back({f[2], ca, bc});
      next.push_back({ab, bc, ca});
    }
    faces = std::move(next);
  }

  const std::uint32_t base = static_cast<std::uint32_t>(out.verts.size());
  const Vec3 c{iso.xyzr[0], iso.xyzr[1], iso.xyzr[2]};
  for (const Vec3& d : dirs) {
    out.verts.push_back({c.x + d.x * r, c.y + d.y * r, c.z + d.z * r});
    out.vnormals.push_back({d.x, d.y, d.z});
    out.atom_id.push_back(iso.orig + 1u);
  }
  for (const Tri& f : faces) {
    out.faces.push_back({base + static_cast<std::uint32_t>(f[0]),
                         base + static_cast<std::uint32_t>(f[1]),
                         base + static_cast<std::uint32_t>(f[2])});
    out.face_type.push_back(3);  // convex (contact), MSMS code
  }
}

}  // namespace

std::shared_ptr<RSCache> compute_rs_from_array(
    const std::vector<std::array<double, 4>>& xyzr, double radius_probe) {
  auto cache = std::make_shared<RSCache>();
  for (const std::vector<std::uint32_t>& comp : sas_components(xyzr, radius_probe)) {
    if (comp.size() == 1) {
      cache->isolated.push_back(
          RSCache::Isolated{comp[0], xyzr[static_cast<std::size_t>(comp[0])]});
      continue;
    }
    // Sub-array in ascending original order, so a single-component input
    // reproduces geom_from_array(xyzr) exactly (identity `orig`).
    std::vector<std::array<double, 4>> sub;
    sub.reserve(comp.size());
    for (std::uint32_t i : comp) sub.push_back(xyzr[static_cast<std::size_t>(i)]);
    RSCache::Component c;
    c.rs = compute_rs(geom_from_array(sub), radius_probe);
    c.orig = comp;
    cache->comps.push_back(std::move(c));
  }
  return cache;
}

MeshResult build_mesh_from_cache(const std::shared_ptr<RSCache>& rs,
                                 double mesh_size, bool fuse) {
  MeshResult out;
  for (const RSCache::Component& comp : rs->comps)
    append_surface(out, build_mesh(comp.rs, mesh_size, fuse), comp.orig);
  for (const RSCache::Isolated& iso : rs->isolated)
    append_isolated_sphere(out, iso, mesh_size);
  return out;
}

MeshResult build_surface_from_array(
    const std::vector<std::array<double, 4>>& xyzr, double radius_probe,
    double mesh_size, bool fuse) {
  return build_mesh_from_cache(compute_rs_from_array(xyzr, radius_probe),
                               mesh_size, fuse);
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
  out.nonfinite_vertices = static_cast<std::uint32_t>(rep.nonfinite_vertices);
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
