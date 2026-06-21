#pragma once
// Convex (SAS->SES contracted) spherical SES patch construction --- faithful port
// of the convex module's data_SESsphpat_convex (MESH PATH ONLY, ext unused).
//
// Per atom i it slices the loops / patches built on the SAS (DataLoop/DataPat),
// rebuilds the direction-corrected local loops_i0, the global segment0 record
// matrix and the local circle0 record matrix via mod_seg_loop_cir (ported here as
// a FREE FUNCTION -- the want_area block where it lived in sas_patches was deferred
// there), writes the per-segment neighbour-radius Rj, and meshes every spherical
// patch on atom i with mesh_sphpat (the SAS->SES contraction is INSIDE mesh_sphpat).
//
// FAITHFULNESS: 1-based-with-dummy-row-0 indexing throughout; segment0 is the
// GLOBAL (nsegment+1)x12 buffer reused across atoms (allocated once, overwritten
// per atom by mod_seg_loop_cir); Rj is a mutable (nsegment+1) vector written by
// mod_seg_loop_cir and read by mesh_sphpat. boundary_tag = ("atom", i).
#include "meshms/exterior.hpp"
#include "meshms/geom.hpp"
#include "meshms/mesh_state.hpp"
#include "meshms/sas.hpp"
#include "meshms/sas_patches.hpp"

namespace meshms {

// Mesh all convex (SAS->SES contracted) spherical patches into `state`.
//
// `ext` is unused (both ext/int patches are meshed -- the MATLAB ext/int branch
// only chooses the figure); pass nullptr.
void data_SESsphpat_convex(MeshState& state, const Geom& geom, const DataI& di,
                           const DataCir& dc, const DataSeg& ds,
                           const DataLoop& dl, const DataPat& dp, const Ext* ext,
                           double Rp, double d);

}  // namespace meshms
