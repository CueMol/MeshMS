#pragma once
// ID-based boundary vertex fusion --- faithful port of meshms/fusion.py
// fuse_by_id (the warn=False pipeline path only).
//
// Each patch mesher tags every boundary vertex with the SES primitive(s) it lies
// on (Tag: ('atom',a)/('probe',p)/('tseam',k,row)). Fusion merges vertices that
// share a tag AND coincide within eps, via a union-find so a dual-tagged toroidal
// corner bridges the convex (atom) and concave (probe) groups at a junction.
// Interior vertices (empty TagList == Python None) are never fused.
//
// The eps test is only a coincidence confirmation scoped to a shared-tag group:
// per-tag spatial hash grids (cell = eps) find the eps-coincident partners in
// O(1) -- if dd < eps^2 then every axis-difference is < eps, so the two cells
// differ by at most 1 per axis, and probing the 27 neighbour cells catches every
// coincident pair exactly (identical V2/F2 to a linear scan).
//
// warn=False: the displaced-partner near-miss diagnostic (the coarse grid_warn
// grid) is dropped -- the pipeline always calls warn=False, and the near-miss
// path never affects the fusion result.
#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

#include "meshms/mesh.hpp"
#include "meshms/mesh_state.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

// (V2, F2) = fuse_by_id(V, F, tags, eps): union-find fusion of vertices sharing a
// Tag AND coinciding within eps. F2 is remapped onto V2 with degenerate triangles
// (two equal indices) dropped. V2 is built from union-find roots in first-seen
// order. Faithful to fuse_by_id(..., warn=False).
std::pair<std::vector<Vec3>, std::vector<Tri>> fuse_by_id(
    const std::vector<Vec3>& V, const std::vector<Tri>& F,
    const std::vector<TagList>& tags, double eps = 1e-6);

// (V2, F2, atom_id2): same fusion, additionally remapping per-vertex atom ids.
// The surviving vertex of each merged group is the union-find root (first-seen
// order, the same vertex whose coordinate is kept in V2), so atom_id2 takes that
// root's atom id -- a neighbouring atom, sufficient for the colouring use case.
// If `atom_id` is empty, atom_id2 is returned empty (the 2-output overload above
// delegates here and drops it). V2/F2 are bit-identical to the 2-output version.
std::tuple<std::vector<Vec3>, std::vector<Tri>, std::vector<int32_t>> fuse_by_id(
    const std::vector<Vec3>& V, const std::vector<Tri>& F,
    const std::vector<TagList>& tags, const std::vector<int32_t>& atom_id,
    double eps = 1e-6);

}  // namespace meshms
