// Verify Feature 3 (public C++ facade, meshms.hpp):
//   * build_surface_from_array(xyzr,...) is bit-for-bit identical (verts/faces/
//     normals/atom_id) to the internal build_surface(file,...) -- guarantees the
//     from-array Geom construction matches read_xyzr exactly;
//   * the cache path compute_rs_from_array + build_mesh_from_cache reproduces the
//     one-shot build_surface_from_array across two densities.
#include <array>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "meshms/meshms.hpp"
#include "meshms/geom.hpp"
#include "meshms/pipeline.hpp"
#include "meshms/vec3.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

// Read the xyzr file via the internal reader, then re-express it as the facade's
// {x,y,z,r} array (0-based). Using read_xyzr here keeps the array == file by
// construction, so the test isolates the facade's geom_from_array path.
std::vector<std::array<double, 4>> load_array(const std::string& xyzr) {
  Geom g = read_xyzr(xyzr);
  std::vector<std::array<double, 4>> out;
  out.reserve(static_cast<std::size_t>(g.M));
  for (int a = 1; a <= g.M; ++a) {
    out.push_back({g.centers[static_cast<std::size_t>(a)].x,
                   g.centers[static_cast<std::size_t>(a)].y,
                   g.centers[static_cast<std::size_t>(a)].z,
                   g.R[static_cast<std::size_t>(a)]});
  }
  return out;
}

// Bit-exact: facade MeshResult vs internal Surface.
void check_vs_surface(const std::string& mol, const MeshResult& m,
                      const Surface& s) {
  CHECK(m.verts.size() == s.V.size());
  CHECK(m.faces.size() == s.F.size());
  CHECK(m.vnormals.size() == s.NV.size());
  CHECK(m.atom_id.size() == s.atom_id.size());

  std::size_t vdiff = 0, fdiff = 0, ndiff = 0, adiff = 0;
  const std::size_t nV = std::min(m.verts.size(), s.V.size());
  for (std::size_t k = 0; k < nV; ++k) {
    if (m.verts[k][0] != s.V[k].x || m.verts[k][1] != s.V[k].y ||
        m.verts[k][2] != s.V[k].z)
      ++vdiff;
    if (m.vnormals[k][0] != s.NV[k].x || m.vnormals[k][1] != s.NV[k].y ||
        m.vnormals[k][2] != s.NV[k].z)
      ++ndiff;
    if (m.atom_id[k] != static_cast<std::uint32_t>(s.atom_id[k])) ++adiff;
  }
  const std::size_t nF = std::min(m.faces.size(), s.F.size());
  for (std::size_t k = 0; k < nF; ++k) {
    if (m.faces[k][0] != static_cast<std::uint32_t>(s.F[k][0]) ||
        m.faces[k][1] != static_cast<std::uint32_t>(s.F[k][1]) ||
        m.faces[k][2] != static_cast<std::uint32_t>(s.F[k][2]))
      ++fdiff;
  }
  CHECK(vdiff == 0);
  CHECK(fdiff == 0);
  CHECK(ndiff == 0);
  CHECK(adiff == 0);
  std::printf("  %-9s facade==internal V=%zu F=%zu (vdiff=%zu fdiff=%zu)\n",
              mol.c_str(), m.verts.size(), m.faces.size(), vdiff, fdiff);
}

// Bit-exact: two facade MeshResults (verts/faces/atom_id).
void check_results_equal(const MeshResult& a, const MeshResult& b) {
  CHECK(a.verts.size() == b.verts.size());
  CHECK(a.faces.size() == b.faces.size());
  CHECK(a.atom_id.size() == b.atom_id.size());
  std::size_t vdiff = 0, fdiff = 0, adiff = 0;
  const std::size_t nV = std::min(a.verts.size(), b.verts.size());
  for (std::size_t k = 0; k < nV; ++k) {
    if (a.verts[k] != b.verts[k]) ++vdiff;
    if (a.atom_id[k] != b.atom_id[k]) ++adiff;
  }
  const std::size_t nF = std::min(a.faces.size(), b.faces.size());
  for (std::size_t k = 0; k < nF; ++k)
    if (a.faces[k] != b.faces[k]) ++fdiff;
  CHECK(vdiff == 0);
  CHECK(fdiff == 0);
  CHECK(adiff == 0);
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  std::vector<std::array<double, 4>> arr = load_array(xyzr);

  // (1) from-array one-shot == from-file internal build.
  MeshResult m = build_surface_from_array(arr, 1.4, 0.5, /*fuse=*/false);
  Surface s = build_surface(xyzr, 1.4, 0.5, /*fuse=*/false);
  check_vs_surface(mol, m, s);

  // (2) cache path reproduces the one-shot result at two densities.
  std::shared_ptr<RSCache> rs = compute_rs_from_array(arr, 1.4);
  for (double d : {0.5, 0.25}) {
    MeshResult cached = build_mesh_from_cache(rs, d, /*fuse=*/false);
    MeshResult oneshot = build_surface_from_array(arr, 1.4, d, /*fuse=*/false);
    check_results_equal(cached, oneshot);
  }
}

}  // namespace

int main() {
  check_mol("2spheres");
  check_mol("3spheres");
  check_mol("tetra");
  TEST_MAIN_RETURN();
}
