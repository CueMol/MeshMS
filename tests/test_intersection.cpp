// Cross-check interstructure() against the Python reference CSR for each of
// {2spheres, 3spheres, tetra, ArgArg}. The ref file (.interstructure.txt) is:
//   line 1: "M <int>"
//   then M lines, atom a = 1..M: "<a> <count> <nbr1> <nbr2> ..." (ascending).
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "meshms/geom.hpp"
#include "meshms/intersection.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

// Parsed reference: M and per-atom (1..M) ascending neighbour lists.
struct Ref {
  int M{0};
  std::vector<std::vector<int>> nbr;  // size M+1, [0] dummy/empty
};

Ref read_ref(const std::string& path) {
  std::ifstream in(path);
  CHECK(static_cast<bool>(in));
  Ref ref;
  std::string line;

  // line 1: "M <int>"
  std::getline(in, line);
  {
    std::istringstream ls(line);
    std::string tag;
    ls >> tag >> ref.M;
    CHECK(tag == "M");
  }
  ref.nbr.assign(static_cast<std::size_t>(ref.M) + 1, {});

  // M lines: "<a> <count> <nbr...>"
  for (int i = 0; i < ref.M; ++i) {
    if (!std::getline(in, line)) break;
    std::istringstream ls(line);
    int a = 0, count = 0;
    ls >> a >> count;
    std::vector<int> nbrs;
    nbrs.reserve(static_cast<std::size_t>(count));
    int v;
    while (ls >> v) nbrs.push_back(v);
    CHECK(static_cast<int>(nbrs.size()) == count);
    ref.nbr[static_cast<std::size_t>(a)] = std::move(nbrs);
  }
  return ref;
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  const std::string refpath =
      std::string(MESHMS_REF_DIR) + "/" + mol + ".interstructure.txt";

  Geom g = read_xyzr(xyzr);
  Neighbors nb = interstructure(g, 1.4);
  Ref ref = read_ref(refpath);

  CHECK(nb.M == ref.M);
  CHECK(g.M == ref.M);

  for (int a = 1; a <= ref.M; ++a) {
    const std::vector<int>& want = ref.nbr[static_cast<std::size_t>(a)];
    auto got = nb.of(a);
    // Count must match.
    CHECK(nb.count(a) == static_cast<int>(want.size()));
    if (nb.count(a) != static_cast<int>(want.size())) {
      std::fprintf(stderr, "  mol=%s atom=%d count got=%d want=%zu\n",
                   mol.c_str(), a, nb.count(a), want.size());
      continue;
    }
    // Exact ascending values must match.
    for (std::size_t t = 0; t < want.size(); ++t) {
      CHECK(got[t] == want[t]);
      if (got[t] != want[t]) {
        std::fprintf(stderr,
                     "  mol=%s atom=%d pos=%zu got=%d want=%d\n",
                     mol.c_str(), a, t, got[t], want[t]);
      }
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
