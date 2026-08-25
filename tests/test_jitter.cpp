// Verify Part D (jitter="auto" NaN-retry fallback) on the degenerate-symmetry
// fullerene fixture:
//   * Jitter::None reproduces the documented degenerate symptom -- the faithful
//     mesh contains non-finite (NaN) vertices;
//   * Jitter::Auto recovers a usable mesh -- no NaN vertices, V and F are non-
//     empty, and analyze_mesh(...).area is finite and strictly positive.
// PROPERTY-BASED ONLY: jitter_centers uses std::normal_distribution, whose output
// is std-lib-implementation dependent, so the jittered mesh is NOT byte-stable
// across platforms. We therefore assert properties, never byte-compare.
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "meshms/geom.hpp"
#include "meshms/meshms.hpp"
#include "meshms/vec3.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

// Read the xyzr file via the internal reader, then re-express it as the facade's
// {x,y,z,r} array (0-based) -- same load pattern as test_meshms.cpp.
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

// True if any vertex component of the result is non-finite (NaN or inf).
bool result_has_nan(const MeshResult& m) {
  for (const std::array<double, 3>& v : m.verts) {
    if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
      return true;
  }
  return false;
}

}  // namespace

int main() {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/fullerene.xyzr";
  std::vector<std::array<double, 4>> arr = load_array(xyzr);
  CHECK(!arr.empty());

  // (a) Jitter::None -- the strongly symmetric fullerene meshes faithfully into a
  // degenerate, NaN-poisoned surface. This pins the failure mode the fallback
  // exists to fix; if a future change makes the faithful mesh clean, this CHECK
  // flags that the premise (and thus the need for Auto) has shifted. STRICT
  // BUILDS ONLY: the relaxed-FP (MESHMS_FP=fast) angle comparison resolves the
  // degenerate symmetry differently and happens to mesh fullerene cleanly, so
  // the premise pin does not hold there (the Auto properties below still do).
#ifndef MESHMS_FP_FAST
  MeshResult none =
      build_surface_from_array(arr, 1.4, 0.5, /*fuse=*/false, Jitter::None);
  CHECK(result_has_nan(none));
  std::printf("  None  : V=%zu F=%zu has_nan=%d\n", none.verts.size(),
              none.faces.size(), static_cast<int>(result_has_nan(none)));
#endif

  // (b) Jitter::Auto -- the NaN-retry fallback recovers a finite, usable mesh.
  MeshResult various =
      build_surface_from_array(arr, 1.4, 0.5, /*fuse=*/false, Jitter::Auto);
  CHECK(!result_has_nan(various));
  CHECK(!various.verts.empty());
  CHECK(!various.faces.empty());

  MeshReport rep = analyze_mesh(various);
  CHECK(std::isfinite(rep.area));
  CHECK(rep.area > 0.0);
  std::printf("  Auto  : V=%zu F=%zu has_nan=%d area=%.6g\n",
              various.verts.size(), various.faces.size(),
              static_cast<int>(result_has_nan(various)), rep.area);

  TEST_MAIN_RETURN();
}
