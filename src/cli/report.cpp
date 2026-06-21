#include "report.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

namespace meshms_cli {

bool ends_with(const std::string& s, const std::string& suffix) {
  return s.size() >= suffix.size() &&
         s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string nearest_atoms(const std::vector<Atom>& atoms,
                          const std::array<double, 3>& ctr) {
  std::vector<std::pair<double, int>> d;
  d.reserve(atoms.size());
  for (std::size_t a = 0; a < atoms.size(); ++a) {
    const double dx = atoms[a][0] - ctr[0];
    const double dy = atoms[a][1] - ctr[1];
    const double dz = atoms[a][2] - ctr[2];
    d.emplace_back(std::sqrt(dx * dx + dy * dy + dz * dz), static_cast<int>(a) + 1);
  }
  const int k = std::min<int>(3, static_cast<int>(d.size()));
  std::partial_sort(d.begin(), d.begin() + k, d.end());
  std::string out;
  for (int i = 0; i < k; ++i) {
    if (i) out += ", ";
    out += "a" + std::to_string(d[static_cast<std::size_t>(i)].second);
  }
  return out;
}

void usage() {
  std::fprintf(stderr,
               "usage: meshms_cli INPUT.xyzr -o OUT [--probe 1.4] [--mesh-size 0.5]\n"
               "                   [--no-fuse] [--fuse-cusps] [--weld-tol 1e-4] [--ply]\n"
               "                   [--vertex-normals] [--no-jitter] [-v|--verbose]\n"
               "  default output: MSMS <OUT>.vert + <OUT>.face   (--ply -> ASCII PLY)\n"
               "  --no-jitter: mesh faithful coords only (skip the auto symmetry\n"
               "               jitter that de-NaNs degenerate molecules, e.g. fullerene)\n");
}

}  // namespace meshms_cli
