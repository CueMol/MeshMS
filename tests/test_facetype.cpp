// Verify the per-face SES type (MeshResult.face_type) -- the MSMS face-type codes
// 3=convex(contact) / 2=concave(spheric reentrant) / 1=toroidal(toric reentrant)
// that back the CLI's MSMS .vert/.face output:
//   * aligned with faces and valued in {1,2,3} on the build path (fuse off & on);
//   * a 2-atom blob has only convex+toroidal (no probe seats on 3 atoms => no
//     concave); a real fragment (ArgArg) exercises all three codes;
//   * stays aligned through remove_flaps; empty after close_cusps (faces rebuilt);
//   * deterministic: identical with 1 and 8 TBB threads (output-preserving).
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#if defined(MESHMS_WITH_TBB)
#include <tbb/global_control.h>
#endif

#include "meshms/meshms.hpp"
#include "meshms/geom.hpp"
#include "meshms/vec3.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

std::vector<std::array<double, 4>> load_array(const std::string& mol) {
  Geom g = read_xyzr(std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr");
  std::vector<std::array<double, 4>> out;
  out.reserve(static_cast<std::size_t>(g.M));
  for (int a = 1; a <= g.M; ++a)
    out.push_back({g.centers[static_cast<std::size_t>(a)].x,
                   g.centers[static_cast<std::size_t>(a)].y,
                   g.centers[static_cast<std::size_t>(a)].z,
                   g.R[static_cast<std::size_t>(a)]});
  return out;
}

// counts[t] = #faces with type t (t in 1..3); returns the #out-of-range codes.
std::size_t tally(const MeshResult& m, std::array<std::size_t, 4>& counts) {
  counts = {0, 0, 0, 0};
  std::size_t bad = 0;
  for (std::uint8_t t : m.face_type) {
    if (t >= 1 && t <= 3) ++counts[t];
    else ++bad;
  }
  return bad;
}

// Build (fuse=false), returning the per-type counts; checks alignment + range.
std::array<std::size_t, 4> check_build(const std::string& mol) {
  const std::vector<std::array<double, 4>> atoms = load_array(mol);

  MeshResult mf = build_surface_from_array(atoms, 1.4, 0.5, /*fuse=*/false);
  std::array<std::size_t, 4> c{};
  CHECK(mf.face_type.size() == mf.faces.size());
  CHECK(tally(mf, c) == 0);  // every face carries a valid SES code
  CHECK(c[3] > 0);           // convex/contact always present
  CHECK(c[1] > 0);           // toroidal always present (>= 2 overlapping atoms)

  // Fused build: fuse drops degenerate faces, face_type must follow the survivors.
  MeshResult mfu = build_surface_from_array(atoms, 1.4, 0.5, /*fuse=*/true);
  std::array<std::size_t, 4> cu{};
  CHECK(mfu.face_type.size() == mfu.faces.size());
  CHECK(tally(mfu, cu) == 0);

  // remove_flaps drops flap faces; alignment must survive the filter.
  MeshResult mrf = remove_flaps(mfu);
  std::array<std::size_t, 4> cr{};
  CHECK(mrf.face_type.size() == mrf.faces.size());
  CHECK(tally(mrf, cr) == 0);

  // close_cusps rebuilds faces (fan-fill) -> face_type intentionally empty.
  MeshResult mc = close_cusps(mfu);
  CHECK(mc.face_type.empty());

  std::printf("  %-8s nF=%zu  convex=%zu concave=%zu toroidal=%zu  fused nF=%zu\n",
              mol.c_str(), mf.faces.size(), c[3], c[2], c[1], mfu.faces.size());
  return c;
}

#if defined(MESHMS_WITH_TBB)
MeshResult build_threads(const std::vector<std::array<double, 4>>& atoms, int nthreads) {
  tbb::global_control gc(tbb::global_control::max_allowed_parallelism,
                         static_cast<std::size_t>(nthreads));
  return build_surface_from_array(atoms, 1.4, 0.5, /*fuse=*/false);
}
void check_threads(const std::string& mol) {
  const std::vector<std::array<double, 4>> atoms = load_array(mol);
  MeshResult a = build_threads(atoms, 1);
  MeshResult b = build_threads(atoms, 8);
  CHECK(a.face_type.size() == b.face_type.size());
  std::size_t diff = 0;
  const std::size_t n = std::min(a.face_type.size(), b.face_type.size());
  for (std::size_t i = 0; i < n; ++i)
    if (a.face_type[i] != b.face_type[i]) ++diff;
  CHECK(diff == 0);
}
#endif

}  // namespace

int main() {
  // A 2-atom blob: convex caps + the toroidal neck, but no concave (a probe never
  // seats on 3 atoms at once).
  std::array<std::size_t, 4> two = check_build("2spheres");
  CHECK(two[2] == 0);

  // A real fragment exercises all three SES primitive kinds.
  std::array<std::size_t, 4> arg = check_build("ArgArg");
  CHECK(arg[2] > 0);  // concave present

  check_build("tetra");

#if defined(MESHMS_WITH_TBB)
  check_threads("ArgArg");
#endif
  TEST_MAIN_RETURN();
}
