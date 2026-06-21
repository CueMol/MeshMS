#pragma once
// CLI-local MSMS .vert/.face writer. This is deliberately NOT part of the
// library/facade: MSMS is a legacy format-specific output whose only consumer is
// this CLI (promoting it would pin the format to the ABI). The format is
// documented once in docs/API.md.
#include <string>

#include "meshms/meshms.hpp"

namespace meshms_cli {

// Write the mesh as an MSMS vertex/face pair (<base>.vert + <base>.face), the
// reduced-surface mesher's native format (Sanner). Each file has three header
// lines (two comments + "count nsph density probe") then one record per line:
//   .vert: x y z  nx ny nz  face_number  closest_sphere  vertex_type
//   .face: v1 v2 v3 (1-based)  face_type  analytic_face_number
// closest_sphere is the 1-based owning atom (MeshResult.atom_id; 0 = unknown).
// face_type / vertex_type are the SES component codes (3 contact, 2 spheric
// reentrant, 1 toric reentrant); vertex_type is the smallest incident face_type
// (toric < reentrant < contact), so seam/edge vertices read as the lower-numbered
// component, matching MSMS's convention. The analytic face_number fields are 0
// (libMeshMS does not number analytical faces). When face_type is unavailable
// (e.g. after --fuse-cusps) the type fields are written as 0.
void write_msms(const std::string& base, const meshms::MeshResult& m, int natom,
                double density, double probe);

}  // namespace meshms_cli
