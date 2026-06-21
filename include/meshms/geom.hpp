#pragma once
// Molecule geometry and xyzr input --- faithful port of meshms/geom.py.
//
// 1-based indexing with a dummy padding row 0: centers[i]/R[i] refer to atom i
// for i in 1..M (PORTING_CONTRACT.md). Row 0 is left zero.
#include <string>
#include <vector>

#include "meshms/vec3.hpp"

namespace meshms {

struct Geom {
  int M{0};
  std::vector<Vec3> centers;  // size M+1, [0] dummy
  std::vector<double> R;      // size M+1, [0] dummy
};

// Read an xyzr file (>=4 whitespace-separated columns: x y z radius ...).
// Mirrors numpy.loadtxt: blank lines and lines starting with '#' are skipped.
Geom read_xyzr(const std::string& path);

// Return a copy of geom with atom centers displaced by a Gaussian of std
// `magnitude` (symmetry-breaking perturbation; see geom.py jitter_centers).
//
// NOTE: this CANNOT bit-match numpy's PCG64 + ziggurat normal sampler. It is
// ONLY the degenerate-symmetry fallback (e.g. fullerene) and is never used for
// the golden molecules (which build with no jitter). A deterministic
// std::mt19937_64(seed) gaussian is therefore acceptable here.
Geom jitter_centers(const Geom& geom, double magnitude, int seed = 0);

}  // namespace meshms
