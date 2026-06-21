#pragma once
// Exterior (eSAS) extraction --- faithful port of meshms/exterior.py data_ext().
//
// Classifies which spherical patches, segments, circles and intersection points
// lie on the *exterior* solvent-accessible surface:
//   1. Build the patch-neighbour graph (per boundary segment/circle, the patch
//      on the OTHER sphere). Uses the seg_owner/cir_owner reverse-index version
//      of the Python source (the _owner lookup), which is the ACTUAL code path.
//   2. Pick an initial exterior patch via the leftmost sphere / leftmost boundary
//      point (smallest x cannot be buried inside the molecule).
//   3. Flood-fill the neighbour graph from that patch.
//   4. Mark segments/circles bounding exterior patches, then the intersection
//      points on exterior segments.
//
// NOTE: the module-level neighbor_patch() helper in exterior.py is DEAD (defined
// but never called by data_ext --- data_ext uses the _owner reverse index). It is
// therefore NOT ported here.
//
// 1-based-with-dummy-row-0 indexing is preserved: every flag array has a [0]
// dummy slot and real entries 1..count; flags are 0/1.
#include <cstdint>
#include <vector>

#include "meshms/csr.hpp"
#include "meshms/geom.hpp"
#include "meshms/sas.hpp"
#include "meshms/sas_patches.hpp"

namespace meshms {

// Exterior-surface membership flags (1 = on the eSAS). All 1-based-with-dummy.
struct Ext {
  std::vector<int32_t> I;        // (nI+1)      ext_I[i] = 1 if point i is exterior
  std::vector<int32_t> circle;   // (ncircle+1) ext_circle[k]
  std::vector<int32_t> segment;  // (nsegment+1) ext_segment[k]
  std::vector<int32_t> patch;    // (npatches+1) ext_patch[k]
};

// Port of data_ext.m / exterior.py data_ext.
Ext data_ext(const Geom& geom, const Neighbors& nb, const DataI& data_i,
             const DataCir& data_c, const DataSeg& data_seg,
             const DataLoop& data_loop, const DataPat& data_pat, double Rp);

}  // namespace meshms
