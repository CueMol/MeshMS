#pragma once
// Advancing-front spherical patch mesher --- faithful port of the meshing module
// (Mesh_Generation/mesh_sphpat.m, advancing_front_approach.m, compute_NV.m,
// activefront.m).
//
// BYTE-FRAGILE: every scalar `x ** 2` in the Python source becomes pysq(x); every
// float term and evaluation order is reproduced with the vec3 helpers (cross/dot/
// norm/sign). The advancing front recurses (collapse_nonneighbor1/2 recurse into
// advancing_front_approach) -- a plain recursive port (Python raises the limit).
//
// Indexing (1-based-with-dummy, see PORTING_CONTRACT.md):
//   * P : std::vector<Vec3>, P[0] dummy, real points at P[1..Np] (size Np+1).
//   * T : std::vector<std::array<int,3>> of 1-based [a,b,c] triples; NO dummy slot,
//         Nt == T.size().
//   * Ae: std::vector<std::array<int,2>> of 1-based [tail,head] pairs with a dummy
//         Ae[0]; MATLAB Ae(a:b,:) becomes Ae[a..b] re-wrapped with a fresh dummy.
#include <array>
#include <cstdint>
#include <optional>
#include <vector>

#include "meshms/mesh_state.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

// loop[1..loopsize] = local segment0 indices; loop[0] dummy.
using Loop = std::vector<int>;

// One inactive front loop (other loops of the same patch) -- port of activefront.m.
struct ActiveFront {
  int meshed = 0;
  std::vector<std::array<int, 2>> Ae;  // 1-based [tail,head], Ae[0] dummy
  int Nae = 0;
};

// A locally-built patch mesh returned by a mesher so the driver can later mesh
// patches in parallel and merge them in order. P is 1-based (P[0] dummy); T is
// 1-based [a,b,c] triples into P. NV is per-face aligned with T; vids are
// per-vertex tag lists aligned with P[1..Np]. emit=false means "no add_patch"
// (the old r_sphere==0 / Nt==0 no-add paths): the driver must skip it so vertex
// bases do not shift.
struct LocalMesh {
  std::vector<Vec3> P;                       // 1-based, P[0] dummy
  std::vector<std::array<int, 3>> T;         // 1-based [a,b,c] triples
  std::vector<Vec3> NV;                      // per-face normals aligned with T
  std::vector<TagList> vids;                 // per-vertex tags aligned with P[1..]
  std::vector<int32_t> vatom;                // per-vertex owning atom (0 => unknown)
  bool emit = false;                         // false => driver must NOT add_patch
};

// Per-face outward normal for each triangle of a spherical patch (compute_NV.m).
//   T : 1-based [a,b,c] triples into P; P : 1-based (P[0] dummy).
//   arg_NV = +1 convex (r_sphere>Rp), -1 concave.
std::vector<Vec3> compute_NV(const std::vector<std::array<int, 3>>& T,
                             const std::vector<Vec3>& P, const Vec3& c_sphere,
                             double arg_NV);

// Mesh a (convex or concave) spherical patch with the advancing front and RETURN
// the result as a LocalMesh (the driver calls state.add_patch(lm.P, lm.T, lm.NV,
// lm.vids) only when lm.emit). emit is true exactly on the old add_patch path
// (r_sphere != 0 AND Nt > 0); otherwise emit is false and nothing is added.
//   c_sphere, r_sphere : the patch sphere.
//   loops[k][1..loops[k].size()-1] : local segment0 indices for loop k (loops[0]/[..][0] dummy).
//   segment0 : (*,12) record matrix ([s][0] dummy; [1..3]=c, [4..6]=n, [7]=r,
//              [8..10]=spoint, [11]=angle), [0] dummy row.
//   circle0  : (*,9) record matrix ([k][1..3]=c, [4..6]=n, [7]=r, [8]=torus radius
//              (0 if none)), [0] dummy row, [k][0] dummy col.
//   patches[1..patchesize] : signed (+k -> loops[k], -k -> circle0[k]); patches[0] dummy.
//   Rp, d : probe radius / mesh size.
//   Rj : nullptr or maps GLOBAL segment index -> neighbour VdW radius (Rj[0] dummy);
//        when null, ctx.Rj is none.
//   boundary_tag : tag for boundary vertices (kind 0 => none -> interior tags).
LocalMesh mesh_sphpat(const Vec3& c_sphere, double r_sphere,
                      const std::vector<Loop>& loops,
                      const std::vector<std::array<double, 12>>& segment0,
                      const std::vector<std::array<double, 9>>& circle0,
                      const std::vector<int>& patches, int patchesize, double Rp,
                      double d, const std::vector<double>* Rj, Tag boundary_tag);

}  // namespace meshms
