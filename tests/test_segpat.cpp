// Cross-check data_Seg_Pat() (mesh path, want_area=false) against the Python
// golden <mol>.segpat.txt for {2spheres, 3spheres, tetra, ArgArg}.
//
// Integer topology MUST match EXACTLY: nsegment; each segment i/j/p1/p2/direct;
// nsatom/satom per atom; nloops, each loop size+entries; loops_index; npatches,
// each patch size+patch_atom+signed entries; patches_index. The ncrasegment
// (n,A,r,radian) and Rj compare within 1e-9. A topology mismatch means a float
// term-order / pysq / stable-sort deviation flipped a decision.
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "meshms/geom.hpp"
#include "meshms/intersection.hpp"
#include "meshms/sas.hpp"
#include "meshms/sas_patches.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

void check_close(const std::string& mol, const char* what, int idx, double got,
                 double want) {
  const double tol = 1e-9 * (1.0 + std::fabs(want));
  if (!(std::fabs(got - want) <= tol)) {
    std::fprintf(stderr, "FAIL %s %s[%d] got=%.17g want=%.17g (d=%.3g)\n",
                 mol.c_str(), what, idx, got, want, std::fabs(got - want));
    ++g_fail;
  }
}

struct SegRow {
  int i, j, p1, p2, direct;
  double nx, ny, nz, Ax, Ay, Az, r, radian, Rj;
};

struct Ref {
  int nsegment{0};
  int M{0};
  std::vector<SegRow> seg;                       // 1-based: seg[k], seg[0] dummy
  std::vector<std::vector<int>> SA;              // SA[a] = seg ids on atom a
  int nloops{0};
  std::vector<std::vector<int>> LP;              // LP[k] = global seg ids
  std::vector<std::array<int, 2>> LPI;           // LPI[a] = {start,end}
  int npatches{0};
  std::vector<std::vector<int>> PT;              // PT[k] = signed entries
  std::vector<int> PTatom;                       // PTatom[k] = patch atom
  std::vector<std::array<int, 2>> PTI;           // PTI[a] = {start,end}
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
    if (tag == "nsegment") {
      ls >> ref.nsegment;
      ref.seg.assign(static_cast<std::size_t>(ref.nsegment) + 1, SegRow{});
      // LP/PT vectors are sized lazily once nloops/npatches are read.
    } else if (tag == "M") {
      ls >> ref.M;
      ref.SA.assign(static_cast<std::size_t>(ref.M) + 1, {});
      ref.LPI.assign(static_cast<std::size_t>(ref.M) + 1, {0, 0});
      ref.PTI.assign(static_cast<std::size_t>(ref.M) + 1, {0, 0});
    } else if (tag == "SEG") {
      int k = 0;
      SegRow s{};
      ls >> k >> s.i >> s.j >> s.p1 >> s.p2 >> s.direct >> s.nx >> s.ny >> s.nz >>
          s.Ax >> s.Ay >> s.Az >> s.r >> s.radian >> s.Rj;
      ref.seg[static_cast<std::size_t>(k)] = s;
    } else if (tag == "SA") {
      int a = 0, nsa = 0;
      ls >> a >> nsa;
      std::vector<int> ids;
      int v;
      while (ls >> v) ids.push_back(v);
      CHECK(static_cast<int>(ids.size()) == nsa);
      ref.SA[static_cast<std::size_t>(a)] = std::move(ids);
    } else if (tag == "nloops") {
      ls >> ref.nloops;
      ref.LP.assign(static_cast<std::size_t>(ref.nloops) + 1, {});
    } else if (tag == "LP") {
      int k = 0, lsz = 0;
      ls >> k >> lsz;
      std::vector<int> ids;
      int v;
      while (ls >> v) ids.push_back(v);
      CHECK(static_cast<int>(ids.size()) == lsz);
      ref.LP[static_cast<std::size_t>(k)] = std::move(ids);
    } else if (tag == "LPI") {
      int a = 0, start = 0, end = 0;
      ls >> a >> start >> end;
      ref.LPI[static_cast<std::size_t>(a)] = {start, end};
    } else if (tag == "npatches") {
      ls >> ref.npatches;
      ref.PT.assign(static_cast<std::size_t>(ref.npatches) + 1, {});
      ref.PTatom.assign(static_cast<std::size_t>(ref.npatches) + 1, 0);
    } else if (tag == "PT") {
      int k = 0, psz = 0, pa = 0;
      ls >> k >> psz >> pa;
      std::vector<int> ids;
      int v;
      while (ls >> v) ids.push_back(v);
      CHECK(static_cast<int>(ids.size()) == psz);
      ref.PT[static_cast<std::size_t>(k)] = std::move(ids);
      ref.PTatom[static_cast<std::size_t>(k)] = pa;
    } else if (tag == "PTI") {
      int a = 0, start = 0, end = 0;
      ls >> a >> start >> end;
      ref.PTI[static_cast<std::size_t>(a)] = {start, end};
    }
  }
  return ref;
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  const std::string refpath =
      std::string(MESHMS_REF_DIR) + "/" + mol + ".segpat.txt";

  Geom g = read_xyzr(xyzr);
  Neighbors nb = interstructure(g, 1.4);
  auto [di, dc] = data_I_Cir(g, nb, 1.4);
  auto [ds, dl, dp] = data_Seg_Pat(g, nb, di, dc, 1.4);
  Ref ref = read_ref(refpath);

  CHECK(g.M == ref.M);

  // --- segments: count EXACT; i/j/p1/p2/direct EXACT; ncra+Rj within 1e-9. ---
  CHECK(ds.nsegment == ref.nsegment);
  if (ds.nsegment == ref.nsegment) {
    for (int k = 1; k <= ref.nsegment; ++k) {
      const SegRow& s = ref.seg[static_cast<std::size_t>(k)];
      const auto& seg = ds.segment[static_cast<std::size_t>(k)];
      const auto& ncra = ds.ncrasegment[static_cast<std::size_t>(k)];
      bool topo_ok = seg[0] == s.i && seg[1] == s.j && seg[2] == s.p1 &&
                     seg[3] == s.p2 && seg[4] == s.direct;
      CHECK(topo_ok);
      if (!topo_ok) {
        std::fprintf(stderr,
                     "  %s SEG[%d] got (i,j,p1,p2,d)=(%d,%d,%d,%d,%d) "
                     "want (%d,%d,%d,%d,%d)\n",
                     mol.c_str(), k, seg[0], seg[1], seg[2], seg[3], seg[4], s.i,
                     s.j, s.p1, s.p2, s.direct);
      }
      check_close(mol, "nx", k, ncra[0], s.nx);
      check_close(mol, "ny", k, ncra[1], s.ny);
      check_close(mol, "nz", k, ncra[2], s.nz);
      check_close(mol, "Ax", k, ncra[3], s.Ax);
      check_close(mol, "Ay", k, ncra[4], s.Ay);
      check_close(mol, "Az", k, ncra[5], s.Az);
      check_close(mol, "r", k, ncra[6], s.r);
      check_close(mol, "radian", k, ncra[7], s.radian);
      check_close(mol, "Rj", k, ds.Rj[static_cast<std::size_t>(k)], s.Rj);
    }
  }

  // --- nsatom / satom per atom EXACT. ---
  for (int a = 1; a <= ref.M; ++a) {
    const std::vector<int>& want = ref.SA[static_cast<std::size_t>(a)];
    const std::vector<int32_t>& got = ds.satom[static_cast<std::size_t>(a)];
    CHECK(static_cast<int>(got.size()) == static_cast<int>(want.size()));
    if (got.size() == want.size()) {
      for (std::size_t t = 0; t < want.size(); ++t) CHECK(got[t] == want[t]);
    } else {
      std::fprintf(stderr, "  %s SA[%d] size got=%zu want=%zu\n", mol.c_str(), a,
                   got.size(), want.size());
    }
  }

  // --- loops: count EXACT; each loop size+entries EXACT; loops_index EXACT. ---
  CHECK(dl.nloops == ref.nloops);
  if (dl.nloops == ref.nloops) {
    for (int k = 1; k <= ref.nloops; ++k) {
      const std::vector<int>& want = ref.LP[static_cast<std::size_t>(k)];
      // dl.loops[k] is 1-based with [0] dummy.
      const std::vector<int32_t>& got = dl.loops[static_cast<std::size_t>(k)];
      const int gsz = static_cast<int>(got.size()) - 1;  // drop [0] dummy
      CHECK(gsz == static_cast<int>(want.size()));
      bool ok = gsz == static_cast<int>(want.size());
      if (ok) {
        for (std::size_t t = 0; t < want.size(); ++t) {
          if (got[t + 1] != want[t]) ok = false;
          CHECK(got[t + 1] == want[t]);
        }
      }
      if (!ok) {
        std::fprintf(stderr, "  %s LP[%d] size got=%d want=%zu\n", mol.c_str(), k,
                     gsz, want.size());
      }
    }
  }
  for (int a = 1; a <= ref.M; ++a) {
    CHECK(dl.loops_index[static_cast<std::size_t>(a)][0] ==
          ref.LPI[static_cast<std::size_t>(a)][0]);
    CHECK(dl.loops_index[static_cast<std::size_t>(a)][1] ==
          ref.LPI[static_cast<std::size_t>(a)][1]);
  }

  // --- patches: count EXACT; each patch size+atom+signed entries EXACT;
  //     patches_index EXACT. ---
  CHECK(dp.npatches == ref.npatches);
  if (dp.npatches == ref.npatches) {
    for (int k = 1; k <= ref.npatches; ++k) {
      const std::vector<int>& want = ref.PT[static_cast<std::size_t>(k)];
      const std::vector<int32_t>& got = dp.patches[static_cast<std::size_t>(k)];
      const int gsz = static_cast<int>(got.size()) - 1;  // drop [0] dummy
      CHECK(gsz == static_cast<int>(want.size()));
      CHECK(dp.patch_atom[static_cast<std::size_t>(k)] ==
            ref.PTatom[static_cast<std::size_t>(k)]);
      bool ok = gsz == static_cast<int>(want.size());
      if (ok) {
        for (std::size_t t = 0; t < want.size(); ++t) {
          if (got[t + 1] != want[t]) ok = false;
          CHECK(got[t + 1] == want[t]);
        }
      }
      if (!ok || dp.patch_atom[static_cast<std::size_t>(k)] !=
                     ref.PTatom[static_cast<std::size_t>(k)]) {
        std::fprintf(stderr,
                     "  %s PT[%d] got size=%d atom=%d want size=%zu atom=%d\n",
                     mol.c_str(), k, gsz,
                     dp.patch_atom[static_cast<std::size_t>(k)], want.size(),
                     ref.PTatom[static_cast<std::size_t>(k)]);
      }
    }
  }
  for (int a = 1; a <= ref.M; ++a) {
    CHECK(dp.patches_index[static_cast<std::size_t>(a)][0] ==
          ref.PTI[static_cast<std::size_t>(a)][0]);
    CHECK(dp.patches_index[static_cast<std::size_t>(a)][1] ==
          ref.PTI[static_cast<std::size_t>(a)][1]);
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
