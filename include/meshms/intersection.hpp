#pragma once
// Intersection structure between SAS balls --- faithful port of
// the intersection module's interstructure().
//
// Two SAS balls a, b (radii R+Rp) intersect iff
//     dist(c_a, c_b) < (R_a + Rp) + (R_b + Rp).
// Reproduces the +/-2 spatial-hash search of the Python and emits the result
// directly as a CSR instead of a dense
// (M+1, kmax+1) rectangle. Each row is strictly ascending-unique.
#include "meshms/csr.hpp"
#include "meshms/geom.hpp"

namespace meshms {

// Build the per-atom intersecting-ball CSR. 1-based-with-dummy-row-0 indexing
// is preserved: nb.of(0) is empty, real atoms are nb.of(1..M).
Neighbors interstructure(const Geom& geom, double Rp);

}  // namespace meshms
