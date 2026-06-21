// Faithful port of meshms/geom.py (read_xyzr, jitter_centers).
#include "meshms/geom.hpp"

#include <cctype>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace meshms {

Geom read_xyzr(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("read_xyzr: cannot open file: " + path);
  }

  // Parse like numpy.loadtxt: skip blank lines and '#' comment lines, split
  // each remaining line on whitespace into float columns.
  std::vector<std::vector<double>> rows;
  std::string line;
  while (std::getline(in, line)) {
    // Find first non-space char to detect blank / comment lines.
    std::size_t i = 0;
    while (i < line.size() &&
           std::isspace(static_cast<unsigned char>(line[i]))) {
      ++i;
    }
    if (i == line.size()) continue;     // blank line
    if (line[i] == '#') continue;       // comment line

    std::istringstream ls(line);
    std::vector<double> cols;
    double v;
    while (ls >> v) cols.push_back(v);
    if (!cols.empty()) rows.push_back(std::move(cols));
  }

  if (rows.empty()) {
    throw std::runtime_error("read_xyzr: no data rows in file: " + path);
  }
  if (rows[0].size() < 4) {
    throw std::runtime_error(
        "xyzr file must have >=4 columns, got " +
        std::to_string(rows[0].size()));
  }

  const int m = static_cast<int>(rows.size());
  Geom g;
  g.M = m;
  g.centers.assign(m + 1, Vec3{0.0, 0.0, 0.0});  // row 0 dummy
  g.R.assign(m + 1, 0.0);                          // entry 0 dummy
  for (int a = 1; a <= m; ++a) {
    const std::vector<double>& c = rows[a - 1];
    // numpy.loadtxt raises on a ragged row; mirror that for every row, not just
    // the first (defensive -- the well-formed golden .xyzr files never hit this).
    if (c.size() < 4) {
      throw std::runtime_error(
          "xyzr file row " + std::to_string(a) + " has <4 columns");
    }
    g.centers[a] = Vec3{c[0], c[1], c[2]};
    g.R[a] = c[3];
  }
  return g;
}

Geom jitter_centers(const Geom& geom, double magnitude, int seed) {
  // NOTE: deterministic but NOT bit-identical to numpy's default_rng(seed)
  // PCG64 + ziggurat normal. Only the degenerate-symmetry fallback; never used
  // for the golden molecules. The radii and the dummy row 0 are left untouched.
  std::mt19937_64 rng(static_cast<std::uint64_t>(seed));
  std::normal_distribution<double> gauss(0.0, magnitude);

  Geom out;
  out.M = geom.M;
  out.centers = geom.centers;  // copy (row 0 dummy preserved)
  out.R = geom.R;              // copy
  for (int a = 1; a <= geom.M; ++a) {
    out.centers[a].x += gauss(rng);
    out.centers[a].y += gauss(rng);
    out.centers[a].z += gauss(rng);
  }
  return out;
}

}  // namespace meshms
