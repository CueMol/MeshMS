#pragma once
// Molecule geometry and xyzr input --- faithful port of the geom module.
//
// 1-based indexing with a dummy padding row 0: centers[i]/R[i] refer to atom i
// for i in 1..M (docs/INTERNALS.md). Row 0 is left zero.
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
// `magnitude` (symmetry-breaking perturbation; see geom's jitter_centers).
//
// This is the wired degenerate-symmetry fallback used by pipeline::run_auto
// (the build_surface_from_array jitter="auto" path): a strongly symmetric
// molecule whose faithful mesh is NaN-poisoned is re-meshed from perturbed
// centers. NOTE: it CANNOT bit-match numpy's PCG64 + ziggurat normal sampler --
// a deterministic std::mt19937_64(seed) gaussian is used, so the perturbed
// result is platform-implementation dependent (property-tested, not byte-
// compared). It never fires for the golden molecules, which build with no
// jitter, so their results stay bit-exact.
Geom jitter_centers(const Geom& geom, double magnitude, int seed = 0);

}  // namespace meshms
