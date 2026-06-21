#pragma once
// SAS segments, loops and spherical patches --- faithful port of
// the sas_patches module's data_Seg_Pat() (MESH PATH ONLY, want_area=False).
//
// The want_area block (mod_seg_loop_cir + area_spherical + DataAV + segment0 +
// Area_sphpat) is DEAD in the mesh path and is DEFERRED here (see the want_area
// parameter; passing true is not yet implemented).
//
// Float term / eval order, the np.argsort(kind="stable") -> std::stable_sort and
// the pysq() (Python scalar `x ** 2`) calls are kept identical to the numpy
// source so the per-decision integer topology (segments, loops, patches) matches
// the Python golden EXACTLY. 1-based-with-dummy-row-0 indexing is preserved
// throughout (every record vector has a [0] dummy entry).
#include <array>
#include <cstdint>
#include <tuple>
#include <vector>

#include "meshms/csr.hpp"
#include "meshms/geom.hpp"
#include "meshms/sas.hpp"

namespace meshms {

// SAS circular segments (DataSeg). All record vectors carry a [0] dummy entry.
struct DataSeg {
  std::vector<std::array<int32_t, 5>> segment;     // [k] = {i, j, p1, p2, direct}; [0] dummy
  int nsegment{0};
  std::vector<std::array<double, 8>> ncrasegment;  // [k] = {nx,ny,nz, Ax,Ay,Az, r, radian}; [0] dummy
  std::vector<std::vector<int32_t>> satom;         // satom[a] = seg ids on atom a (nsatom == .size())
  std::vector<double> Rj;                          // [k] per segment; [0] dummy (mesh path: all 0)
};

// SAS loops (DataLoop). loops[k] = the global segment ids on the k-th loop.
struct DataLoop {
  std::vector<std::vector<int32_t>> loops;         // loops[k] = global seg ids (loopsize == .size()); [0] dummy
  int nloops{0};
  std::vector<std::array<int32_t, 2>> loops_index; // loops_index[a] = {start, end}, a in 0..M
};

// SAS spherical patches (DataPat). patches[k] entries are SIGNED: +loop / -circle
// (local indices on the patch atom).
struct DataPat {
  std::vector<std::vector<int32_t>> patches;       // patches[k] = signed entries (patchesize == .size()); [0] dummy
  int npatches{0};
  std::vector<std::array<int32_t, 2>> patches_index; // patches_index[a] = {start, end}, a in 0..M
  std::vector<int32_t> patch_atom;                 // patch_atom[k] -> atom; [0] dummy
};

// Compute Segments, Loops & Patches on the SAS (mesh path, want_area=false).
//
// nb is the intersection CSR (interstructure(geom, Rp)): Row[i] == nb.count(i),
// inter.M_int[i,row] == nb.of(i)[row-1] (1-based row, ascending).
//
// want_area is reserved for the deferred Gauss-Bonnet area/volume report and is
// NOT yet implemented; passing true is unsupported (asserts/throws).
std::tuple<DataSeg, DataLoop, DataPat> data_Seg_Pat(const Geom& geom,
                                                    const Neighbors& nb,
                                                    const DataI& data_i,
                                                    const DataCir& data_c,
                                                    double Rp,
                                                    bool want_area = false);

}  // namespace meshms
