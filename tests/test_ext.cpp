// Cross-check data_ext() (exterior / eSAS extraction) against the Python golden
// <mol>.ext.txt for {2spheres, 3spheres, tetra, ArgArg}.
//
// All integer flags (0/1) MUST match EXACTLY: ext.I[1..nI], ext.circle[1..ncircle],
// ext.segment[1..nsegment], ext.patch[1..npatches]. A mismatch means the
// neighbour graph, leftmost-patch seed, flood-fill or marking diverged.
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "meshms/exterior.hpp"
#include "meshms/geom.hpp"
#include "meshms/intersection.hpp"
#include "meshms/sas.hpp"
#include "meshms/sas_patches.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

struct Ref {
  int nI{0}, ncircle{0}, nsegment{0}, npatches{0};
  std::vector<int> extI, extC, extS, extP;  // 0-based lists of flags (count long)
};

// Parse "<tag> <count>" header lines and "ext<X> <f1> ... <fn>" flag lines.
Ref read_ref(const std::string& path) {
  std::ifstream in(path);
  CHECK(static_cast<bool>(in));
  Ref ref;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::string tag;
    if (!(ls >> tag)) continue;  // skip blank lines
    if (tag == "nI") {
      ls >> ref.nI;
    } else if (tag == "ncircle") {
      ls >> ref.ncircle;
    } else if (tag == "nsegment") {
      ls >> ref.nsegment;
    } else if (tag == "npatches") {
      ls >> ref.npatches;
    } else if (tag == "extI") {
      int v;
      while (ls >> v) ref.extI.push_back(v);
    } else if (tag == "extC") {
      int v;
      while (ls >> v) ref.extC.push_back(v);
    } else if (tag == "extS") {
      int v;
      while (ls >> v) ref.extS.push_back(v);
    } else if (tag == "extP") {
      int v;
      while (ls >> v) ref.extP.push_back(v);
    }
  }
  return ref;
}

// Compare a 1-based-with-dummy flag vector (got) against the 0-based golden list
// (want) of length 'count'.
void check_flags(const std::string& mol, const char* what, int count,
                 const std::vector<int32_t>& got, const std::vector<int>& want) {
  CHECK(static_cast<int>(want.size()) == count);
  // got is sized count+1 ([0] dummy, 1..count real).
  CHECK(static_cast<int>(got.size()) == count + 1);
  if (static_cast<int>(want.size()) != count ||
      static_cast<int>(got.size()) != count + 1) {
    std::fprintf(stderr, "  %s ext%s size mismatch: count=%d got=%zu want=%zu\n",
                 mol.c_str(), what, count, got.size(), want.size());
    return;
  }
  for (int k = 1; k <= count; ++k) {
    if (got[static_cast<std::size_t>(k)] != want[static_cast<std::size_t>(k - 1)]) {
      std::fprintf(stderr, "  %s ext%s[%d] got=%d want=%d\n", mol.c_str(), what,
                   k, got[static_cast<std::size_t>(k)],
                   want[static_cast<std::size_t>(k - 1)]);
      ++g_fail;
    }
  }
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  const std::string refpath =
      std::string(MESHMS_REF_DIR) + "/" + mol + ".ext.txt";

  Geom g = read_xyzr(xyzr);
  Neighbors nb = interstructure(g, 1.4);
  auto [di, dc] = data_I_Cir(g, nb, 1.4);
  auto [ds, dl, dp] = data_Seg_Pat(g, nb, di, dc, 1.4);
  Ext ext = data_ext(g, nb, di, dc, ds, dl, dp, 1.4);

  Ref ref = read_ref(refpath);

  // Counts the port produced must match the golden header counts.
  CHECK(di.nI == ref.nI);
  CHECK(dc.ncircle == ref.ncircle);
  CHECK(ds.nsegment == ref.nsegment);
  CHECK(dp.npatches == ref.npatches);

  check_flags(mol, "I", ref.nI, ext.I, ref.extI);
  check_flags(mol, "C", ref.ncircle, ext.circle, ref.extC);
  check_flags(mol, "S", ref.nsegment, ext.segment, ref.extS);
  check_flags(mol, "P", ref.npatches, ext.patch, ref.extP);
}

}  // namespace

int main() {
  check_mol("2spheres");
  check_mol("3spheres");
  check_mol("tetra");
  check_mol("ArgArg");
  TEST_MAIN_RETURN();
}
