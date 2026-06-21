// Verify the multi-component support in the public facade (meshms.cpp):
//   * an isolated atom is meshed as its full vdW sphere (watertight, area ~
//     4*pi*r^2, radial normals, atom_id set, face_type == 3);
//   * a multi-component input produces every component, each bit-for-bit equal
//     to that component built alone (vertex offset + atom_id remap only);
//   * the cache path equals the one-shot path on multi-component input;
//   * an empty input yields an empty mesh (no throw).
// Everything is driven through the public facade only, with in-code inputs.
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <vector>

#include "meshms/meshms.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

using Xyzr = std::vector<std::array<double, 4>>;

// Two overlapping spheres (the 2spheres golden molecule's layout): one
// ordinary connected component for the faithful pipeline.
Xyzr two_spheres(double dx) {
  return {{dx + 0.0, 0.0, 0.0, 1.2}, {dx + 0.0, -2.0, 0.0, 1.2}};
}

// Bit-exact comparison of a component slice of `all` against `sub` built
// alone: vertices/normals equal, faces equal after the vertex offset, atom_id
// equal after the index offset.
void check_component_slice(const MeshResult& all, std::size_t v0, std::size_t f0,
                           const MeshResult& sub, std::uint32_t id_off) {
  CHECK(all.verts.size() >= v0 + sub.verts.size());
  CHECK(all.faces.size() >= f0 + sub.faces.size());
  std::size_t vdiff = 0, ndiff = 0, fdiff = 0, adiff = 0, tdiff = 0;
  for (std::size_t k = 0; k < sub.verts.size(); ++k) {
    if (all.verts[v0 + k] != sub.verts[k]) ++vdiff;
    if (all.vnormals[v0 + k] != sub.vnormals[k]) ++ndiff;
    if (all.atom_id[v0 + k] != sub.atom_id[k] + id_off) ++adiff;
  }
  for (std::size_t k = 0; k < sub.faces.size(); ++k) {
    for (int j = 0; j < 3; ++j) {
      if (all.faces[f0 + k][static_cast<std::size_t>(j)] !=
          sub.faces[k][static_cast<std::size_t>(j)] + static_cast<std::uint32_t>(v0)) {
        ++fdiff;
        break;
      }
    }
    if (all.face_type[f0 + k] != sub.face_type[k]) ++tdiff;
  }
  CHECK(vdiff == 0);
  CHECK(ndiff == 0);
  CHECK(fdiff == 0);
  CHECK(adiff == 0);
  CHECK(tdiff == 0);
}

void test_isolated_atom() {
  const double r = 1.5;
  const Xyzr one = {{2.0, -3.0, 5.0, r}};
  const MeshResult m = build_surface_from_array(one, 1.4, 0.5, /*fuse=*/true);

  CHECK(!m.verts.empty());
  CHECK(!m.faces.empty());
  CHECK(m.vnormals.size() == m.verts.size());
  CHECK(m.atom_id.size() == m.verts.size());
  CHECK(m.face_type.size() == m.faces.size());

  // Every vertex sits on the vdW sphere; the normal is the radial direction.
  std::size_t off_sphere = 0, off_normal = 0, bad_id = 0, bad_type = 0;
  for (std::size_t k = 0; k < m.verts.size(); ++k) {
    const double dx = m.verts[k][0] - 2.0;
    const double dy = m.verts[k][1] + 3.0;
    const double dz = m.verts[k][2] - 5.0;
    const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (std::fabs(d - r) > 1e-12) ++off_sphere;
    const double dot =
        (dx * m.vnormals[k][0] + dy * m.vnormals[k][1] + dz * m.vnormals[k][2]) / d;
    if (std::fabs(dot - 1.0) > 1e-12) ++off_normal;
    if (m.atom_id[k] != 1u) ++bad_id;
  }
  for (std::uint8_t t : m.face_type)
    if (t != 3) ++bad_type;
  CHECK(off_sphere == 0);
  CHECK(off_normal == 0);
  CHECK(bad_id == 0);
  CHECK(bad_type == 0);

  // Watertight closed sphere with the analytic area/volume (icosphere at
  // mesh_size 0.5 underestimates both by a few percent -- inscribed polyhedron).
  const MeshReport rep = analyze_mesh(m);
  CHECK(rep.watertight);
  const double area = 4.0 * std::numbers::pi * r * r;
  const double vol = area * r / 3.0;
  CHECK_NEAR(rep.area, area, 0.05 * area);
  CHECK(rep.signed_volume > 0.0);
  CHECK_NEAR(rep.signed_volume, vol, 0.05 * vol);
}

void test_component_plus_isolated() {
  const Xyzr base = two_spheres(0.0);
  Xyzr all = base;
  all.push_back({100.0, 0.0, 0.0, 1.5});  // far away: isolated

  const MeshResult sub = build_surface_from_array(base, 1.4, 0.5, true);
  const MeshResult m = build_surface_from_array(all, 1.4, 0.5, true);

  // The connected component is emitted first, bit-for-bit as built alone.
  check_component_slice(m, 0, 0, sub, 0);

  // The isolated sphere follows: its vertices all carry atom_id 3.
  CHECK(m.verts.size() > sub.verts.size());
  std::size_t bad_id = 0;
  for (std::size_t k = sub.verts.size(); k < m.verts.size(); ++k)
    if (m.atom_id[k] != 3u) ++bad_id;
  CHECK(bad_id == 0);

  const MeshReport rep = analyze_mesh(m);
  CHECK(rep.watertight);
}

void test_two_components() {
  const Xyzr a = two_spheres(0.0);
  const Xyzr b = two_spheres(50.0);
  Xyzr all = a;
  all.insert(all.end(), b.begin(), b.end());

  const MeshResult ma = build_surface_from_array(a, 1.4, 0.5, true);
  const MeshResult mb = build_surface_from_array(b, 1.4, 0.5, true);
  const MeshResult m = build_surface_from_array(all, 1.4, 0.5, true);

  CHECK(m.verts.size() == ma.verts.size() + mb.verts.size());
  CHECK(m.faces.size() == ma.faces.size() + mb.faces.size());
  check_component_slice(m, 0, 0, ma, 0);
  check_component_slice(m, ma.verts.size(), ma.faces.size(), mb, 2);
}

void test_cache_equals_oneshot_multicomponent() {
  Xyzr all = two_spheres(0.0);
  all.push_back({100.0, 0.0, 0.0, 1.5});

  const std::shared_ptr<RSCache> rs = compute_rs_from_array(all, 1.4);
  for (double d : {0.5, 0.25}) {
    const MeshResult cached = build_mesh_from_cache(rs, d, true);
    const MeshResult oneshot = build_surface_from_array(all, 1.4, d, true);
    CHECK(cached.verts == oneshot.verts);
    CHECK(cached.vnormals == oneshot.vnormals);
    CHECK(cached.faces == oneshot.faces);
    CHECK(cached.atom_id == oneshot.atom_id);
    CHECK(cached.face_type == oneshot.face_type);
  }
}

void test_empty_input() {
  const MeshResult m = build_surface_from_array({}, 1.4, 0.5, true);
  CHECK(m.verts.empty());
  CHECK(m.faces.empty());
}

}  // namespace

int main() {
  test_isolated_atom();
  test_component_plus_isolated();
  test_two_components();
  test_cache_equals_oneshot_multicomponent();
  test_empty_input();
  TEST_MAIN_RETURN();
}
