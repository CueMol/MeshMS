// Convert a .xyz or .pdb molecule file to the .xyzr format read_xyzr expects.
//
// Faithful reimplementation of the MolSurfComp reference readers
// (Read/read_xyz.m and Read/read_PDB.m): UFF VdW radii scaled by 1.1.
//   H 1.443, C 1.9255, O 1.75, N 1.83, P 2.0735, S 2.0175, Mg 1.5105 (xyz only)
// xyz: line 1 is a header (atom count and/or name) and is skipped; every other
//      line whose first token is a known element yields an atom; unknown-element
//      lines are skipped (read_xyz.m `otherwise, continue`).
// pdb: ATOM/HETATM records; the element is the single character at column 14
//      (1-based, i.e. the second character of the atom-name field, exactly as
//      read_PDB.m's `switch(line(14))`); unknown elements fall back to 1.75;
//      coordinates come from columns 31-54.
//
// Usage: meshms_mol2xyzr [--format=xyz|pdb] <input> <output.xyzr>
// Without --format the input extension decides (.pdb -> pdb, anything else xyz).
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {

constexpr double kScale = 1.1;  // MolSurfComp scalingfactor

// UFF VdW radius for a one/two-letter element symbol; 0.0 = unknown.
double uff_radius(const std::string& elem) {
  if (elem == "H") return 1.443;
  if (elem == "C") return 1.9255;
  if (elem == "O") return 1.75;
  if (elem == "N") return 1.83;
  if (elem == "P") return 2.0735;
  if (elem == "S") return 2.0175;
  if (elem == "Mg") return 1.5105;
  return 0.0;
}

int convert_xyz(std::istream& in, std::ostream& out) {
  std::string line;
  int n = 0, lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    if (lineno == 1) continue;  // header (count and/or name), as in read_xyz.m
    std::istringstream ss(line);
    std::string tok;
    double x, y, z;
    if (!(ss >> tok)) continue;
    const double r = uff_radius(tok);
    if (r == 0.0) continue;  // unknown element / stray count line: skip
    if (!(ss >> x >> y >> z)) continue;
    char buf[128];
    std::snprintf(buf, sizeof(buf), "%.6f %.6f %.6f %.6f\n", x, y, z, kScale * r);
    out << buf;
    ++n;
  }
  return n;
}

int convert_pdb(std::istream& in, std::ostream& out) {
  std::string line;
  int n = 0;
  while (std::getline(in, line)) {
    if (line.compare(0, 6, "HETATM") != 0 && line.compare(0, 4, "ATOM") != 0)
      continue;
    if (line.size() < 54) continue;  // needs the coordinate columns
    // read_PDB.m: element = line(14); unknown -> 1.75.
    double r = uff_radius(std::string(1, line[13]));
    if (r == 0.0) r = 1.75;
    // read_PDB.m: coordinates from columns 31-54 (1-based).
    double x, y, z;
    if (std::sscanf(line.substr(30, 24).c_str(), "%lf %lf %lf", &x, &y, &z) != 3)
      continue;
    char buf[128];
    std::snprintf(buf, sizeof(buf), "%.6f %.6f %.6f %.6f\n", x, y, z, kScale * r);
    out << buf;
    ++n;
  }
  return n;
}

}  // namespace

int main(int argc, char** argv) {
  std::string format;
  int ai = 1;
  if (ai < argc && std::strncmp(argv[ai], "--format=", 9) == 0) {
    format = argv[ai] + 9;
    ++ai;
  }
  if (argc - ai != 2) {
    std::fprintf(stderr,
                 "usage: meshms_mol2xyzr [--format=xyz|pdb] <input> <output.xyzr>\n");
    return 2;
  }
  const std::string in_path = argv[ai], out_path = argv[ai + 1];
  if (format.empty()) {
    const auto dot = in_path.rfind('.');
    const std::string ext = dot == std::string::npos ? "" : in_path.substr(dot);
    format = (ext == ".pdb" || ext == ".PDB") ? "pdb" : "xyz";
  }
  if (format != "xyz" && format != "pdb") {
    std::fprintf(stderr, "error: unknown format '%s'\n", format.c_str());
    return 2;
  }

  std::ifstream in(in_path);
  if (!in) {
    std::fprintf(stderr, "error: cannot open %s\n", in_path.c_str());
    return 1;
  }
  std::ofstream out(out_path);
  if (!out) {
    std::fprintf(stderr, "error: cannot open %s\n", out_path.c_str());
    return 1;
  }
  const int n = format == "pdb" ? convert_pdb(in, out) : convert_xyz(in, out);
  if (n == 0) {
    std::fprintf(stderr, "error: no atoms converted from %s\n", in_path.c_str());
    return 1;
  }
  std::printf("%s: %d atoms -> %s\n", in_path.c_str(), n, out_path.c_str());
  return 0;
}
