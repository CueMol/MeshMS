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
#include "meshms/csr.hpp"
#include "meshms/exterior.hpp"
#include "meshms/geom.hpp"
#include "meshms/mesh_state.hpp"
#include "meshms/sas.hpp"

namespace meshms {

// Port of SESconcavepat.m (mesh path, av=None). `inter` is the SAS-ball
// intersection CSR (interstructure(geom, Rp)): inter.num_int[i] == inter.count(i),
// inter.M_int[i,row] == inter.of(i)[row-1].
void SESconcavepat(MeshState& state, const Geom& geom, const DataI& di,
                   const Ext& ext, double Rp, double d, const Neighbors& inter);

}  // namespace meshms
