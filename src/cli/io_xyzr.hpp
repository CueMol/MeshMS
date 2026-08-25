#pragma once
// CLI xyzr reader. Atom is the CLI's local {x, y, z, radius} row type and the
// single inter-module dependency (report.hpp also includes this header for it).
//
// INTENTIONAL DUPLICATION: read_xyzr_file mirrors the library's internal
// read_xyzr; it is kept CLI-local on purpose so the front-end stays on the
// C++17-clean facade boundary (the library's reader lives behind C++20 headers).
#include <array>
#include <string>
#include <vector>

namespace meshms_cli {

using Atom = std::array<double, 4>;  // x, y, z, radius

// Parse an xyzr file like numpy.loadtxt (mirrors the library's read_xyzr): skip
// blank and '#'-comment lines, split each row on whitespace, keep the first 4
// columns. Atoms are returned in file order, so a MeshResult.atom_id value i
// refers to the (i-1)-th atom here.
std::vector<Atom> read_xyzr_file(const std::string& path);

}  // namespace meshms_cli
