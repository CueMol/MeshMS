// Verify the density-independent RS cache: compute_rs() once, then build_mesh()
// at several densities must produce the SAME mesh (bit-for-bit V/F/N) as the
// monolithic run(geom, {Rp, d}). This guards Feature 1 (RS separation): the
// split must be output-preserving and the cache reusable across densities.
#include <array>
#include <cstdio>
#include <string>
#include <vector>

#include "meshms/geom.hpp"
#include "meshms/params.hpp"
#include "meshms/pipeline.hpp"
#include "meshms/vec3.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

// Bit-exact equality of two Surfaces' V/F/N (the golden-checked fields).
// `what` labels the comparison in the printed line.
void check_same(const std::string& mol, double d, const char* what,
                const Surface& a, const Surface& b) {
  CHECK(a.V.size() == b.V.size());
  CHECK(a.F.size() == b.F.size());
  CHECK(a.N.size() == b.N.size());
  const std::size_t nV = std::min(a.V.size(), b.V.size());
  std::size_t vdiff = 0;
  for (std::size_t k = 0; k < nV; ++k) {
    if (a.V[k].x != b.V[k].x || a.V[k].y != b.V[k].y || a.V[k].z != b.V[k].z)
      ++vdiff;
  }
  CHECK(vdiff == 0);
  const std::size_t nF = std::min(a.F.size(), b.F.size());
  std::size_t fdiff = 0;
  for (std::size_t k = 0; k < nF; ++k) {
    if (a.F[k][0] != b.F[k][0] || a.F[k][1] != b.F[k][1] ||
        a.F[k][2] != b.F[k][2])
      ++fdiff;
  }
  CHECK(fdiff == 0);
  const std::size_t nN = std::min(a.N.size(), b.N.size());
  std::size_t ndiff = 0;
  for (std::size_t k = 0; k < nN; ++k) {
    if (a.N[k].x != b.N[k].x || a.N[k].y != b.N[k].y || a.N[k].z != b.N[k].z)
      ++ndiff;
  }
  CHECK(ndiff == 0);
  std::printf("  %-9s d=%.3g nV=%zu nF=%zu  %s: %s\n", mol.c_str(), d,
              a.V.size(), a.F.size(), what,
              (vdiff == 0 && fdiff == 0 && ndiff == 0) ? "yes" : "NO");
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  Geom geom = read_xyzr(xyzr);

  // Compute the RS components ONCE; reuse for every density.
  RSComponents rs = compute_rs(geom, 1.4);

  for (double d : {0.5, 0.25}) {
    Surface cached = build_mesh(rs, d, /*fuse=*/false);
    Para para;
    para.radius_probe = 1.4;
    para.mesh_size = d;
    Surface direct = run(geom, para, /*fuse=*/false);
    check_same(mol, d, "cached==direct", cached, direct);
  }

  // Idempotency / call-order independence: rebuilding at a density already meshed
  // from the SAME cache -- after meshing OTHER densities (and a fuse pass) in
  // between -- must reproduce the first result bit-for-bit. Guards against any
  // state leak in build_mesh or mutation of the (const) cache across calls.
  Surface first = build_mesh(rs, 0.5, /*fuse=*/false);
  Surface other = build_mesh(rs, 0.25, /*fuse=*/false);  // interleave a 2nd density
  Surface fpass = build_mesh(rs, 0.5, /*fuse=*/true);     // interleave a fuse pass
  Surface again = build_mesh(rs, 0.5, /*fuse=*/false);    // rebuild the 1st density
  check_same(mol, 0.5, "rebuild==first", first, again);
  (void)other;
  (void)fpass;

  // The cache must also be reusable after a fuse=true mesh (no mutation of rs).
  Surface fused_cached = build_mesh(rs, 0.5, /*fuse=*/true);
  Para para;
  para.radius_probe = 1.4;
  para.mesh_size = 0.5;
  Surface fused_direct = run(geom, para, /*fuse=*/true);
  CHECK(fused_cached.V.size() == fused_direct.V.size());
  CHECK(fused_cached.F.size() == fused_direct.F.size());
  CHECK(fused_cached.N.empty() && fused_direct.N.empty());
}

}  // namespace

int main() {
  check_mol("2spheres");
  check_mol("3spheres");
  check_mol("tetra");
  TEST_MAIN_RETURN();
}
