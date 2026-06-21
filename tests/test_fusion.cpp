// Cross-check ID-based boundary fusion (build_surface, fuse=True) against the
// Python golden <mol>.fused.txt.
//
// The golden is the full build_surface mesh AFTER orient_faces AND fuse_by_id
// (warn=False): fused vertices V2 and remapped 0-based faces F2 (no normals --
// N is dropped when fused). The golden floats are %.17g (full IEEE-double
// round-trip), so each parsed golden double equals the Python double EXACTLY.
//
// PASS criterion: nV/nF EXACT, F2 exact 0-based topology, V2 within 1e-9.
// ADDITIONALLY (reported): manifold_report boundary_edges BEFORE (fuse=False)
// and AFTER (fuse=True) fusion -- the fused mesh must have FEWER boundary edges
// (the C1 seams close), checked here as a regression guard.
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "meshms/mesh_check.hpp"
#include "meshms/pipeline.hpp"
#include "meshms/vec3.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

struct Golden {
  int nV = 0, nF = 0;
  std::vector<Vec3> V;                // 0-based store (V2[k-1])
  std::vector<std::array<int, 3>> F;  // 0-based vertex indices (post-fusion)
};

Golden read_golden(const std::string& path) {
  std::ifstream in(path);
  CHECK(static_cast<bool>(in));
  Golden g;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream ls(line);
    std::string tag;
    ls >> tag;
    if (tag == "nV") {
      ls >> g.nV;
    } else if (tag == "nF") {
      ls >> g.nF;
    } else if (tag == "V") {
      int k = 0;
      Vec3 v{};
      ls >> k >> v.x >> v.y >> v.z;
      g.V.push_back(v);
    } else if (tag == "F") {
      int k = 0;
      std::array<int, 3> t{};
      ls >> k >> t[0] >> t[1] >> t[2];
      g.F.push_back(t);
    }
  }
  return g;
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  const std::string refpath =
      std::string(MESHMS_REF_DIR) + "/" + mol + ".fused.txt";

  // Unfused mesh (for the boundary-edge before/after regression guard).
  Surface raw = build_surface(xyzr, 1.4, 0.5, /*fuse=*/false);
  ManifoldReport mr_before = manifold_report(raw.V, raw.F);

  // Fused mesh -- the subject under test.
  Surface s = build_surface(xyzr, 1.4, 0.5, /*fuse=*/true);
  Golden gold = read_golden(refpath);

  const int gotNV = static_cast<int>(s.V.size());
  const int gotNF = static_cast<int>(s.F.size());
  CHECK(gotNV == gold.nV);
  CHECK(gotNF == gold.nF);
  if (gotNV != gold.nV)
    std::fprintf(stderr, "  %s nV got=%d want=%d\n", mol.c_str(), gotNV, gold.nV);
  if (gotNF != gold.nF)
    std::fprintf(stderr, "  %s nF got=%d want=%d\n", mol.c_str(), gotNF, gold.nF);

  // fuse=True drops normals (Python returns N=None).
  CHECK(s.N.empty());

  // --- vertices: 1e-9 --------------------------------------------------------
  const int nV = std::min(gotNV, gold.nV);
  bool vfirst = true;
  for (int k = 0; k < nV; ++k) {
    const Vec3& got = s.V[static_cast<std::size_t>(k)];
    const Vec3& want = gold.V[static_cast<std::size_t>(k)];
    const bool near = std::fabs(got.x - want.x) <= 1e-9 &&
                      std::fabs(got.y - want.y) <= 1e-9 &&
                      std::fabs(got.z - want.z) <= 1e-9;
    CHECK(near);
    if (!near && vfirst) {
      std::fprintf(stderr,
                   "  %s V[%d] got=(%.17g,%.17g,%.17g) want=(%.17g,%.17g,%.17g)\n",
                   mol.c_str(), k + 1, got.x, got.y, got.z, want.x, want.y, want.z);
      vfirst = false;
    }
  }

  // --- faces: EXACT 0-based topology (post-fusion) ---------------------------
  const int nF = std::min(gotNF, gold.nF);
  bool ffirst = true;
  for (int t = 0; t < nF; ++t) {
    const Tri& f = s.F[static_cast<std::size_t>(t)];
    const std::array<int, 3>& w = gold.F[static_cast<std::size_t>(t)];
    const bool ok = (f[0] == w[0]) && (f[1] == w[1]) && (f[2] == w[2]);
    CHECK(ok);
    if (!ok && ffirst) {
      std::fprintf(stderr, "  %s F[%d] got=(%d,%d,%d) want=(%d,%d,%d)\n",
                   mol.c_str(), t + 1, f[0], f[1], f[2], w[0], w[1], w[2]);
      ffirst = false;
    }
  }

  // --- boundary-edge regression guard: fusion must CLOSE seams ---------------
  ManifoldReport mr_after = manifold_report(s.V, s.F);
  std::printf(
      "  %-9s nV %d->%d nF %d->%d  boundary_edges %d->%d  area=%.6f\n",
      mol.c_str(), static_cast<int>(raw.V.size()), gotNV,
      static_cast<int>(raw.F.size()), gotNF, mr_before.boundary_edges,
      mr_after.boundary_edges, mr_after.area);
  CHECK(mr_after.boundary_edges < mr_before.boundary_edges);
}

}  // namespace

int main() {
  check_mol("2spheres");
  check_mol("3spheres");
  check_mol("tetra");
  check_mol("ArgArg");
  TEST_MAIN_RETURN();
}
