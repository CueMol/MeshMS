#pragma once
// CLI helper bundle: usage text, suffix test, and the nearest-atom labelling
// used by the boundary diagnostics. Includes io_xyzr.hpp purely for Atom (the
// only edge between CLI modules).
#include <array>
#include <string>
#include <vector>

#include "io_xyzr.hpp"

namespace meshms_cli {

bool ends_with(const std::string& s, const std::string& suffix);

// Nearest-3 atom labels (1-based "a<k>") to a point, by Euclidean distance --
// the same selection the diagnostic used (argsort of |atom - ctr|, first 3).
std::string nearest_atoms(const std::vector<Atom>& atoms,
                          const std::array<double, 3>& ctr);

void usage();

}  // namespace meshms_cli
