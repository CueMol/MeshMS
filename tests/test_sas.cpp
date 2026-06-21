// Cross-check data_I_Cir() against the Python golden .sas.txt for each of
// {2spheres, 3spheres, tetra, ArgArg}. Integer topology (nI, Iijk, direction,
// high_I, the Ii id sequences, every IC row, ncircle, the circle atom pair, the
// circleindex per atom) MUST match EXACTLY; coordinates / hightvalue / circle
// A,n,r compare within 1e-9 (the vec3 helpers reproduce the Python term order).
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "meshms/geom.hpp"
#include "meshms/intersection.hpp"
#include "meshms/sas.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

// Relative-or-absolute near check (rel/abs 1e-9), matching the contract.
void check_close(const std::string& mol, const char* what, int idx, double got,
                 double want) {
  const double tol = 1e-9 * (1.0 + std::fabs(want));
  if (!(std::fabs(got - want) <= tol)) {
    std::fprintf(stderr, "FAIL %s %s[%d] got=%.17g want=%.17g (d=%.3g)\n",
                 mol.c_str(), what, idx, got, want, std::fabs(got - want));
    ++g_fail;
  }
}

struct PRow {
  double Ix, Iy, Iz;
  int i, j, k;
  int dij, djk, dki;
  int high_I;
  double hightvalue;
};

struct CRow {
  int ci, cj;
  double Ax, Ay, Az, nx, ny, nz, r;
};

// Parsed golden .sas.txt.
struct Ref {
  int nI{0};
  int M{0};
  std::vector<PRow> P;                                 // 1-based: P[s], P[0] dummy
  std::vector<std::vector<int>> Ii;                    // Ii[a] id list
  // IC keyed by (i,row) -> id list, plus the count token for cross-check.
  std::map<std::pair<int, int>, std::vector<int>> IC;
  int ncircle{0};
  std::vector<CRow> C;                                 // 1-based: C[k], C[0] dummy
  std::vector<std::vector<int>> CI;                    // CI[a] circle id list
};

Ref read_ref(const std::string& path) {
  std::ifstream in(path);
  CHECK(static_cast<bool>(in));
  Ref ref;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream ls(line);
    std::string tag;
    ls >> tag;
    if (tag == "nI") {
      ls >> ref.nI;
      ref.P.assign(static_cast<std::size_t>(ref.nI) + 1, PRow{});
    } else if (tag == "M") {
      ls >> ref.M;
      ref.Ii.assign(static_cast<std::size_t>(ref.M) + 1, {});
      ref.CI.assign(static_cast<std::size_t>(ref.M) + 1, {});
    } else if (tag == "P") {
      int s = 0;
      PRow p{};
      ls >> s >> p.Ix >> p.Iy >> p.Iz >> p.i >> p.j >> p.k >> p.dij >> p.djk >>
          p.dki >> p.high_I >> p.hightvalue;
      ref.P[static_cast<std::size_t>(s)] = p;
    } else if (tag == "Ii") {
      int a = 0, In = 0;
      ls >> a >> In;
      std::vector<int> ids;
      int v;
      while (ls >> v) ids.push_back(v);
      CHECK(static_cast<int>(ids.size()) == In);
      ref.Ii[static_cast<std::size_t>(a)] = std::move(ids);
    } else if (tag == "IC") {
      int i = 0, row = 0, count = 0;
      ls >> i >> row >> count;
      std::vector<int> ids;
      int v;
      while (ls >> v) ids.push_back(v);
      CHECK(static_cast<int>(ids.size()) == count);
      ref.IC[{i, row}] = std::move(ids);
    } else if (tag == "ncircle") {
      ls >> ref.ncircle;
      ref.C.assign(static_cast<std::size_t>(ref.ncircle) + 1, CRow{});
    } else if (tag == "C") {
      int k = 0;
      CRow c{};
      ls >> k >> c.ci >> c.cj >> c.Ax >> c.Ay >> c.Az >> c.nx >> c.ny >> c.nz >>
          c.r;
      ref.C[static_cast<std::size_t>(k)] = c;
    } else if (tag == "CI") {
      int a = 0, n = 0;
      ls >> a >> n;
      std::vector<int> ids;
      int v;
      while (ls >> v) ids.push_back(v);
      CHECK(static_cast<int>(ids.size()) == n);
      ref.CI[static_cast<std::size_t>(a)] = std::move(ids);
    }
  }
  return ref;
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  const std::string refpath =
      std::string(MESHMS_REF_DIR) + "/" + mol + ".sas.txt";

  Geom g = read_xyzr(xyzr);
  Neighbors nb = interstructure(g, 1.4);
  auto [di, dc] = data_I_Cir(g, nb, 1.4);
  Ref ref = read_ref(refpath);

  CHECK(g.M == ref.M);

  // --- nI EXACT; per-point Iijk/direction/high_I EXACT, coords + hightvalue near.
  CHECK(di.nI == ref.nI);
  if (di.nI == ref.nI) {
    for (int s = 1; s <= ref.nI; ++s) {
      const PRow& p = ref.P[static_cast<std::size_t>(s)];
      const auto& ijk = di.Iijk[static_cast<std::size_t>(s)];
      const auto& dir = di.direction[static_cast<std::size_t>(s)];
      CHECK(ijk[0] == p.i && ijk[1] == p.j && ijk[2] == p.k);
      CHECK(dir[0] == p.dij && dir[1] == p.djk && dir[2] == p.dki);
      CHECK(di.high_I[static_cast<std::size_t>(s)] == p.high_I);
      if (ijk[0] != p.i || ijk[1] != p.j || ijk[2] != p.k ||
          dir[0] != p.dij || dir[1] != p.djk || dir[2] != p.dki ||
          di.high_I[static_cast<std::size_t>(s)] != p.high_I) {
        std::fprintf(stderr,
                     "  %s P[%d] got ijk=(%d,%d,%d) dir=(%d,%d,%d) high=%d "
                     "want ijk=(%d,%d,%d) dir=(%d,%d,%d) high=%d\n",
                     mol.c_str(), s, ijk[0], ijk[1], ijk[2], dir[0], dir[1],
                     dir[2], di.high_I[static_cast<std::size_t>(s)], p.i, p.j,
                     p.k, p.dij, p.djk, p.dki, p.high_I);
      }
      const Vec3& I = di.I[static_cast<std::size_t>(s)];
      check_close(mol, "Ix", s, I.x, p.Ix);
      check_close(mol, "Iy", s, I.y, p.Iy);
      check_close(mol, "Iz", s, I.z, p.Iz);
      check_close(mol, "hightvalue", s, di.hightvalue[static_cast<std::size_t>(s)],
                  p.hightvalue);
    }
  }

  // --- In[a] and the Ii[a] id sequence EXACT for every atom.
  for (int a = 1; a <= ref.M; ++a) {
    const std::vector<int>& want = ref.Ii[static_cast<std::size_t>(a)];
    const std::vector<int32_t>& got = di.Ii[static_cast<std::size_t>(a)];
    CHECK(static_cast<int>(got.size()) == static_cast<int>(want.size()));
    if (got.size() == want.size()) {
      for (std::size_t t = 0; t < want.size(); ++t) {
        CHECK(got[t] == want[t]);
      }
    } else {
      std::fprintf(stderr, "  %s Ii[%d] size got=%zu want=%zu\n", mol.c_str(), a,
                   got.size(), want.size());
    }
  }

  // --- every IC (i,row,ids) EXACT (both directions: ref<->got).
  for (const auto& kv : ref.IC) {
    const int i = kv.first.first, row = kv.first.second;
    const std::vector<int>& want = kv.second;
    CHECK(i >= 1 && i <= ref.M);
    CHECK(row >= 1 && row < static_cast<int>(di.I_circle[static_cast<std::size_t>(i)].size()));
    const std::vector<int32_t>& got =
        di.I_circle[static_cast<std::size_t>(i)][static_cast<std::size_t>(row)];
    CHECK(static_cast<int>(got.size()) == static_cast<int>(want.size()));
    bool ok = got.size() == want.size();
    if (ok) {
      for (std::size_t t = 0; t < want.size(); ++t) {
        CHECK(got[t] == want[t]);
        if (got[t] != want[t]) ok = false;
      }
    }
    if (!ok) {
      std::fprintf(stderr, "  %s IC[%d,%d] mismatch (got size=%zu want=%zu)\n",
                   mol.c_str(), i, row, got.size(), want.size());
    }
  }
  // Every non-empty got IC must be present in the ref (no spurious rows).
  for (int i = 1; i <= ref.M; ++i) {
    const auto& rows = di.I_circle[static_cast<std::size_t>(i)];
    for (int row = 1; row < static_cast<int>(rows.size()); ++row) {
      if (!rows[static_cast<std::size_t>(row)].empty()) {
        CHECK(ref.IC.find({i, row}) != ref.IC.end());
      }
    }
  }

  // --- ncircle EXACT; circle ci/cj EXACT, A/n/r near; circleindex per atom EXACT.
  CHECK(dc.ncircle == ref.ncircle);
  if (dc.ncircle == ref.ncircle) {
    for (int k = 1; k <= ref.ncircle; ++k) {
      const CRow& c = ref.C[static_cast<std::size_t>(k)];
      const auto& g9 = dc.circle[static_cast<std::size_t>(k)];
      CHECK(static_cast<int>(g9[1]) == c.ci);
      CHECK(static_cast<int>(g9[2]) == c.cj);
      check_close(mol, "Ax", k, g9[3], c.Ax);
      check_close(mol, "Ay", k, g9[4], c.Ay);
      check_close(mol, "Az", k, g9[5], c.Az);
      check_close(mol, "nx", k, g9[6], c.nx);
      check_close(mol, "ny", k, g9[7], c.ny);
      check_close(mol, "nz", k, g9[8], c.nz);
      check_close(mol, "r", k, g9[9], c.r);
    }
  }
  for (int a = 1; a <= ref.M; ++a) {
    const std::vector<int>& want = ref.CI[static_cast<std::size_t>(a)];
    const std::vector<int32_t>& got = dc.circleindex[static_cast<std::size_t>(a)];
    CHECK(static_cast<int>(got.size()) == static_cast<int>(want.size()));
    if (got.size() == want.size()) {
      for (std::size_t t = 0; t < want.size(); ++t) CHECK(got[t] == want[t]);
    }
  }
}

}  // namespace

int main() {
  check_mol("2spheres");
  check_mol("3spheres");
  check_mol("tetra");
  check_mol("ArgArg");
  TEST_MAIN_RETURN();
}
