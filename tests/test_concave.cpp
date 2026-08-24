// Cross-check SESconcavepat (concave SES mesh, av=None) against the Python golden
// <mol>.concave.txt for {2spheres, 3spheres, tetra, ArgArg}.
//
// The golden is the RAW MeshState after SESconcavepat (NO orient_faces, NO
// fusion): vertices V (1e-9), faces F (EXACT 0-based topology), per-face normals
// N (1e-9), nV/nF EXACT. mesh_sphpat is already exact-verified, so a divergence
// means the concave GEOMETRY (the probe arrangement) differs -- the failure print
// reports the first differing vertex/face.
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "meshms/concave.hpp"
#include "meshms/exterior.hpp"
#include "meshms/geom.hpp"
#include "meshms/intersection.hpp"
#include "meshms/mesh_state.hpp"
#include "meshms/sas.hpp"
#include "meshms/sas_patches.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

struct Golden {
  int nV = 0, nF = 0;
  std::vector<Vec3> V;                  // 0-based store (V[k-1])
  std::vector<std::array<int, 3>> F;    // 0-based vertex indices
  std::vector<Vec3> N;                  // per-face normal
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

void check_vec(const std::string& mol, const char* what, int idx, const Vec3& got,
               const Vec3& want, double tol) {
  if (!(std::fabs(got.x - want.x) <= tol && std::fabs(got.y - want.y) <= tol &&
        std::fabs(got.z - want.z) <= tol)) {
    std::fprintf(stderr,
                 "FAIL %s %s[%d] got=(%.17g,%.17g,%.17g) "
                 "want=(%.17g,%.17g,%.17g)\n",
                 mol.c_str(), what, idx, got.x, got.y, got.z, want.x, want.y,
                 want.z);
    ++g_fail;
  }
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  const std::string refpath =
      std::string(MESHMS_REF_DIR) + "/" + mol + ".concave.txt";

  Geom g = read_xyzr(xyzr);
  Neighbors nb = interstructure(g, 1.4);
  auto [di, dc] = data_I_Cir(g, nb, 1.4);
  auto [ds, dl, dp] = data_Seg_Pat(g, nb, di, dc, 1.4);
  Ext ext = data_ext(g, nb, di, dc, ds, dl, dp, 1.4);

  MeshState st;
  SESconcavepat(st, g, di, ext, 1.4, 0.5, nb);

  Golden gold = read_golden(refpath);

  const int gotNV = static_cast<int>(st.V.size());
  const int gotNF = static_cast<int>(st.F.size());
  CHECK(gotNV == gold.nV);
  CHECK(gotNF == gold.nF);
  if (gotNV != gold.nV)
    std::fprintf(stderr, "  %s nV got=%d want=%d\n", mol.c_str(), gotNV, gold.nV);
  if (gotNF != gold.nF)
    std::fprintf(stderr, "  %s nF got=%d want=%d\n", mol.c_str(), gotNF, gold.nF);

  const int nV = std::min(gotNV, gold.nV);
  for (int k = 0; k < nV; ++k) {
    check_vec(mol, "V", k + 1, st.V[static_cast<std::size_t>(k)],
              gold.V[static_cast<std::size_t>(k)], 1e-9);
  }

  const int nF = std::min(gotNF, gold.nF);
  for (int t = 0; t < nF; ++t) {
    const Tri& f = st.F[static_cast<std::size_t>(t)];
    const std::array<int, 3>& w = gold.F[static_cast<std::size_t>(t)];
    bool ok = (f[0] == w[0]) && (f[1] == w[1]) && (f[2] == w[2]);
    CHECK(ok);
    if (!ok) {
      std::fprintf(stderr, "  %s F[%d] got=(%d,%d,%d) want=(%d,%d,%d)\n",
                   mol.c_str(), t + 1, f[0], f[1], f[2], w[0], w[1], w[2]);
    }
  }

  const int gotN = static_cast<int>(st.N.size());
  const int nN = std::min(gotN, static_cast<int>(gold.N.size()));
  for (int t = 0; t < nN; ++t) {
    check_vec(mol, "N", t + 1, st.N[static_cast<std::size_t>(t)],
              gold.N[static_cast<std::size_t>(t)], 1e-9);
  }
}

}  // namespace

int main() {
  MESHMS_SKIP_IF_FAST();
  check_mol("2spheres");
  check_mol("3spheres");
  check_mol("tetra");
  check_mol("ArgArg");
  TEST_MAIN_RETURN();
}
