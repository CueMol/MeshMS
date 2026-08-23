#pragma once
// Concave SES patch construction --- faithful port of the concave module's
// SESconcavepat (MESH PATH ONLY, av=None).
//
// Meshes the concave (probe-reentrant) SES patches into a MeshState. The mesh is
// produced by two passes only:
//   * the cSAS pass (arg_eSAS=0): for every NON-simple singular triple point
//     (high_I==0), build the local probe arrangement (_data_concavepat) and mesh
//     each resulting concave patch via mesh_sphpat with the FULL neighbour trim;
//   * the simple-triangle pass (high_I==1): mesh the 3-segment concave triangle.
//
// The av-only area/volume accounting (V_cSAS/V_eSAS, _area_spherical,
// _comp_area_vol_concave, the eSAS interior-trim pass) is DEFERRED: it produces no
// mesh and is dead when av==None.
//
// FAITHFULNESS: every scalar `x ** 2` in the Python source becomes pysq(x); the
// angle sorts use std::stable_sort (ties keep input order); np.sign(triple(...))
// in _test1/_test2; EPSILON=1e-10; 1-based-with-dummy-row-0 indexing throughout.
#include <array>
#include <cstdint>
#include <vector>

#include "meshms/csr.hpp"
#include "meshms/exterior.hpp"
#include "meshms/geom.hpp"
#include "meshms/mesh_state.hpp"
#include "meshms/sas.hpp"

namespace meshms {

// Port of SESconcavepat.m (mesh path, av=None). `inter` is the SAS-ball
// intersection CSR (interstructure(geom, Rp)): inter.num_int[i] == inter.count(i),
// inter.M_int[i,row] == inter.of(i)[row-1].
// Equivalent to precompute_concave(...) followed by SESconcavepat_mesh(...).
void SESconcavepat(MeshState& state, const Geom& geom, const DataI& di,
                   const Ext& ext, double Rp, double d, const Neighbors& inter);

// ===== Density-independent concave decomposition ============================
// The expensive half of SESconcavepat -- the per-probe patch decomposition
// (neighbour pruning, boundary/interior intersections, loop + patch trees) --
// depends only on (geom, Rp), NOT on the mesh density d. ProbePatchSet holds
// one cSAS probe's decomposition: exactly the inputs its mesh_sphpat calls
// need, plus the vatom attribution atoms. Compute once via
// precompute_concave(), then SESconcavepat_mesh() re-meshes at any density.
struct ProbePatchSet {
  int k1 = 0;                                    // probe point id (center I[k1])
  std::vector<std::vector<int>> loops;           // [k] 1-based-with-dummy
  std::vector<std::array<double, 12>> segment0;  // record matrix, [0] dummy
  std::vector<std::array<double, 9>> circle0;    // record matrix, [0] dummy
  std::vector<std::vector<int>> patches;         // [k] 1-based-with-dummy
  std::vector<int> patchesize;                   // 1-based-with-dummy
  int npatches = 0;
  // fill_vatom_nearest3 inputs (per-vertex owning-atom attribution):
  std::int32_t a_i = 0, a_j = 0, a_k = 0;
  Vec3 ci{}, cj{}, ck{};
};

struct ConcaveDecomp {
  int nhight = 0;                    // number of cSAS probes
  std::vector<ProbePatchSet> probes; // [1..nhight]; [0] dummy
};

// The density-independent half of SESconcavepat (cSAS probe decomposition).
ConcaveDecomp precompute_concave(const Geom& geom, const DataI& di, double Rp,
                                 const Neighbors& inter);

// Mesh the concave SES patches from a precomputed decomposition at density d.
// Bit-identical to the corresponding part of SESconcavepat.
void SESconcavepat_mesh(MeshState& state, const Geom& geom, const DataI& di,
                        const ConcaveDecomp& decomp, double Rp, double d);

}  // namespace meshms
