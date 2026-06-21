// Cross-check the full SES pipeline (build_surface, fuse=False) against the
// Python golden <mol>.surface.txt for the golden molecule set.
//
// The golden is the full build_surface mesh AFTER orient_faces (NO fusion):
// vertices V, 0-based faces F (AFTER the winding fix), per-face normals N, with
// nV/nF. The golden floats are %.17g (full IEEE-double round-trip), so each
// parsed golden double equals the Python double EXACTLY -- we can check BOTH
// bit-equality (==) and 1e-9 closeness.
//
// PASS criterion: nV/nF EXACT, F exact 0-based topology, V and N within 1e-9.
// ADDITIONALLY (reported, not required): how many V/N doubles are bit-equal (==)
// to the golden -- a 'bit-exact: yes/no (k of n differ)' line per molecule.
// Watertightness is NOT checked (fuse=False meshes have boundary edges by
// design), but manifold_report area is printed for sanity.
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
  std::vector<Vec3> V;                // 0-based store (V[k-1])
  std::vector<std::array<int, 3>> F;  // 0-based vertex indices (after orient)
  std::vector<Vec3> N;                // per-face normal
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
    } else if (tag == "N") {
      int k = 0;
      Vec3 v{};
      ls >> k >> v.x >> v.y >> v.z;
      g.N.push_back(v);
    }
  }
  return g;
}

// Count how many of the 3 components differ within 1e-9 (CHECK) and how many are
// NOT bit-equal (==). Returns the number of components that are not bit-equal.
int check_vec(const std::string& mol, const char* what, int idx, const Vec3& got,
              const Vec3& want, double tol, bool report_first_diff) {
  const bool near = std::fabs(got.x - want.x) <= tol &&
                    std::fabs(got.y - want.y) <= tol &&
                    std::fabs(got.z - want.z) <= tol;
  if (!near) {
    ++g_fail;
    if (report_first_diff) {
      std::fprintf(stderr,
                   "FAIL %s %s[%d] got=(%.17g,%.17g,%.17g) "
                   "want=(%.17g,%.17g,%.17g)\n",
                   mol.c_str(), what, idx, got.x, got.y, got.z, want.x, want.y,
                   want.z);
    }
  }
  int ndiff = 0;
  ndiff += (got.x != want.x);
  ndiff += (got.y != want.y);
  ndiff += (got.z != want.z);
  return ndiff;
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  const std::string refpath =
      std::string(MESHMS_REF_DIR) + "/" + mol + ".surface.txt";

  Surface s = build_surface(xyzr, 1.4, 0.5, /*fuse=*/false);
  Golden gold = read_golden(refpath);

  const int gotNV = static_cast<int>(s.V.size());
  const int gotNF = static_cast<int>(s.F.size());
  CHECK(gotNV == gold.nV);
  CHECK(gotNF == gold.nF);
  if (gotNV != gold.nV)
    std::fprintf(stderr, "  %s nV got=%d want=%d\n", mol.c_str(), gotNV, gold.nV);
  if (gotNF != gold.nF)
    std::fprintf(stderr, "  %s nF got=%d want=%d\n", mol.c_str(), gotNF, gold.nF);

  // --- vertices: 1e-9 (pass) + bit-equal counting (report only) -------------
  const int nV = std::min(gotNV, gold.nV);
  long vbits = 0;        // doubles NOT bit-equal
  bool vfirst = true;    // report only the first V divergence
  for (int k = 0; k < nV; ++k) {
    const Vec3& got = s.V[static_cast<std::size_t>(k)];
    const Vec3& want = gold.V[static_cast<std::size_t>(k)];
    const long before = g_fail;
    const int nd = check_vec(mol, "V", k + 1, got, want, 1e-9, vfirst);
    if (g_fail != before) vfirst = false;
    vbits += nd;
  }

  // --- faces: EXACT 0-based topology (after orient_faces) -------------------
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

  // --- per-face normals: 1e-9 (pass) + bit-equal counting (report only) -----
  const int gotN = static_cast<int>(s.N.size());
  const int nN = std::min(gotN, static_cast<int>(gold.N.size()));
  CHECK(gotN == gold.nF);  // one normal per face
  long nbits = 0;
  bool nfirst = true;
  for (int t = 0; t < nN; ++t) {
    const Vec3& got = s.N[static_cast<std::size_t>(t)];
    const Vec3& want = gold.N[static_cast<std::size_t>(t)];
    const long before = g_fail;
    const int nd = check_vec(mol, "N", t + 1, got, want, 1e-9, nfirst);
    if (g_fail != before) nfirst = false;
    nbits += nd;
  }

  // --- bit-exactness report (not a pass/fail criterion) ---------------------
  const long vtotal = static_cast<long>(nV) * 3;
  const long ntotal = static_cast<long>(nN) * 3;
  const long differ = vbits + nbits;
  std::printf("  %-9s nV=%d nF=%d  bit-exact: %s (%ld of %ld differ)\n",
              mol.c_str(), gotNV, gotNF, (differ == 0 ? "yes" : "no"), differ,
              vtotal + ntotal);

  // --- manifold sanity (area printed; watertightness NOT required) ----------
  ManifoldReport mr = manifold_report(s.V, s.F);
  std::printf(
      "  %-9s area=%.6f  boundary_edges=%d nonmanifold_edges=%d degenerate=%d\n",
      mol.c_str(), mr.area, mr.boundary_edges, mr.nonmanifold_edges,
      mr.degenerate_faces);
}

}  // namespace

int main() {
  check_mol("2spheres");
  check_mol("3spheres");
  check_mol("tetra");
  check_mol("ArgArg");
  check_mol("1YJO");
  check_mol("1ETN");
  TEST_MAIN_RETURN();
}
