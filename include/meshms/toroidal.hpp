#pragma once
// Toroidal SES patches: regular thin tori (mesh_toroide) and singular cusp tori
// (mesh_cusp) --- faithful port of the toroidal module's data_SEStorpat().
//
// Every patch is accumulated UNCONDITIONALLY into the shared MeshState (the
// MATLAB ext/int figure routing only chooses a figure; both branches mesh, so
// the accumulated mesh is identical). `ext` is therefore unused (pass nullptr).
//
// FAITHFULNESS (see docs/INTERNALS.md / SHARED CONTRACT):
//   * Python scalar `x ** 2` -> meshms::pysq(x); _sqrt_pos(x) = (x>0?sqrt(x):0).
//   * math.cos/sin/acos/tan/copysign -> std::cos/sin/acos/tan/copysign (libm).
//   * 1-based-with-dummy P (P[0] dummy) and 1-based T triples.
//   * the structured (N_probe+1)x(N_arc+1) grid (mesh_toroide / index_P_toroide)
//     and the VARIABLE-ROW grid (mesh_cusp; t[] array, row_off prefix sums,
//     index_P, and the intricate t[i+1]==t[i]+1 / t[i+1]==t[i+2]-1 connectivity
//     branches) are reproduced EXACTLY.
//   * orient_face_toroidal orients each face by the analytic torus normal.
//   * vids ID-fusion tags use Tag{}; a multi-tag corner stores the FIRST listed
//     tag (the mesh is identical; full multi-tag fusion is a later step).
#include "meshms/exterior.hpp"
#include "meshms/geom.hpp"
#include "meshms/mesh_state.hpp"
#include "meshms/sas.hpp"
#include "meshms/sas_patches.hpp"

namespace meshms {

// Mesh all toroidal SES patches (arc segments + full circles) into `state`.
// `ext` is unused (both ext/int patches are meshed); pass nullptr. `av` deferred.
void data_SEStorpat(MeshState& state, const Geom& geom, const DataI& di,
                    const DataSeg& ds, const DataCir& dc, const Ext* ext,
                    double Rp, double d);

}  // namespace meshms
