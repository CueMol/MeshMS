// Verify Feature 2 (per-vertex atom_id output):
//   * every vertex has 1 <= atom_id <= M (no unknowns escape);
//   * a vertex lying on exactly ONE atom's VdW sphere (a convex / rim vertex) is
//     attributed to THAT atom (catches mis-attribution);
//   * atom_id is deterministic: identical with 1 and 8 OpenMP threads;
//   * under fuse=true, atom_id stays aligned with V (size == V.size()) and in
//     range [1, M].
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "meshms/geom.hpp"
#include "meshms/params.hpp"
#include "meshms/pipeline.hpp"
#include "meshms/vec3.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

Surface build_threads(const Geom& geom, double d, bool fuse, int nthreads) {
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#else
  (void)nthreads;
#endif
  Para para;
  para.radius_probe = 1.4;
  para.mesh_size = d;
  return run(geom, para, fuse);
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  Geom geom = read_xyzr(xyzr);
  const int M = geom.M;

  Surface s = build_threads(geom, 0.5, /*fuse=*/false, 1);

  // atom_id is aligned with V and 1-based in range.
  CHECK(s.atom_id.size() == s.V.size());
  std::size_t out_of_range = 0;
  for (std::int32_t a : s.atom_id)
    if (a < 1 || a > M) ++out_of_range;
  CHECK(out_of_range == 0);

  // A vertex on exactly one atom's VdW sphere must be owned by that atom.
  const double tol = 1e-6;
  std::size_t on_sphere = 0, mis = 0;
  for (std::size_t k = 0; k < s.V.size(); ++k) {
    const Vec3& v = s.V[k];
    int hit = 0, which = 0;
    for (int a = 1; a <= M; ++a) {
      double dist = norm(v - geom.centers[static_cast<std::size_t>(a)]);
      if (std::fabs(dist - geom.R[static_cast<std::size_t>(a)]) < tol) {
        ++hit;
        which = a;
      }
    }
    if (hit == 1) {
      ++on_sphere;
      if (s.atom_id[k] != which) ++mis;
    }
  }
  CHECK(on_sphere > 0);   // convex/rim vertices must exist
  CHECK(mis == 0);        // and be attributed to their sphere's atom

  // Determinism across thread counts (only meaningful when OpenMP is enabled).
  Surface s8 = build_threads(geom, 0.5, /*fuse=*/false, 8);
  CHECK(s8.atom_id.size() == s.atom_id.size());
  std::size_t tdiff = 0;
  const std::size_t n = std::min(s8.atom_id.size(), s.atom_id.size());
  for (std::size_t k = 0; k < n; ++k)
    if (s8.atom_id[k] != s.atom_id[k]) ++tdiff;
  CHECK(tdiff == 0);

  // fuse=true: atom_id remapped, still aligned and in range.
  Surface sf = build_threads(geom, 0.5, /*fuse=*/true, 1);
  CHECK(sf.atom_id.size() == sf.V.size());
  std::size_t f_out = 0;
  for (std::int32_t a : sf.atom_id)
    if (a < 1 || a > M) ++f_out;
  CHECK(f_out == 0);

  std::printf("  %-9s M=%d nV=%zu on_sphere=%zu mis=%zu thr_diff=%zu  fused nV=%zu\n",
              mol.c_str(), M, s.V.size(), on_sphere, mis, tdiff, sf.V.size());
}

}  // namespace

int main() {
  check_mol("2spheres");
  check_mol("3spheres");
  check_mol("tetra");
  TEST_MAIN_RETURN();
}
