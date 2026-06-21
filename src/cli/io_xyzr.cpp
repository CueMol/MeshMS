#include "io_xyzr.hpp"

#include <cctype>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace meshms_cli {

std::vector<Atom> read_xyzr_file(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    std::fprintf(stderr, "error: cannot open input file: %s\n", path.c_str());
    std::exit(2);
  }
  std::vector<Atom> atoms;
  std::string line;
  int row = 0;
  while (std::getline(in, line)) {
    std::size_t i = 0;
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    if (i == line.size() || line[i] == '#') continue;  // blank / comment line
    std::istringstream ls(line);
    std::vector<double> cols;
    double v;
    while (ls >> v) cols.push_back(v);
    if (cols.empty()) continue;
    ++row;
    if (cols.size() < 4) {
      std::fprintf(stderr, "error: xyzr row %d has <4 columns\n", row);
      std::exit(2);
    }
    atoms.push_back({cols[0], cols[1], cols[2], cols[3]});
  }
  if (atoms.empty()) {
    std::fprintf(stderr, "error: no data rows in %s\n", path.c_str());
    std::exit(2);
  }
  return atoms;
}

}  // namespace meshms_cli
