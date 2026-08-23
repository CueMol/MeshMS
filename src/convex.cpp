// Convex (SAS->SES contracted) spherical SES patch construction --- faithful port
// of the convex module's data_SESsphpat_convex, INCLUDING mod_seg_loop_cir (from
// the sas_patches module lines 260-353, deferred there in the want_area block) ported as a
// free function. See convex.hpp for scope. MESH PATH ONLY (ext unused).
#include "meshms/convex.hpp"
#include "meshms/parallel.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "meshms/meshing.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

namespace {

// ----- record-matrix helpers (1-based-with-dummy rows AND columns) ----------
template <std::size_t N>
inline Vec3 rget3(const std::array<double, N>& row, int c0) {
  return Vec3{row[static_cast<std::size_t>(c0)], row[static_cast<std::size_t>(c0 + 1)],
              row[static_cast<std::size_t>(c0 + 2)]};
}
template <std::size_t N>
inline void rset3(std::array<double, N>& row, int c0, const Vec3& v) {
  row[static_cast<std::size_t>(c0)] = v.x;
  row[static_cast<std::size_t>(c0 + 1)] = v.y;
  row[static_cast<std::size_t>(c0 + 2)] = v.z;
}

// ----- mod_seg_loop_cir (port of Core/mod_seg_loop_cir.m) -------------------
// Build the direction-corrected local loops_i0, overwrite the GLOBAL segment0
// rows for atom i's segments, build the local circle0, and WRITE Rj per segment.
//
//   loops_i  : 1-based-with-dummy rows; loops_i[k] is a 1-based-with-dummy vector
//              of GLOBAL segment ids; loopsize_i[k] = loops_i[k].size()-1.
//   segment0 : GLOBAL (nsegment+1)x12 record matrix (reused across atoms); cols
//              [1..3]=c, [4..6]=n, [7]=r, [8..10]=spoint, [11]=angle.
//   Rj       : GLOBAL (nsegment+1) neighbour-atom VdW radius per segment (WRITTEN).
// Returns loops_i0 (segment order reversed per loop) and circle0 (width-9 rows,
// local to atom i). segment0/Rj are mutated in place.
struct ModOut {
  std::vector<std::vector<int32_t>> loops_i0;       // 1-based-with-dummy rows
  std::vector<std::array<double, 9>> circle0;       // [k][1..7]=c,n,r; [k][8]=torus R
};
ModOut mod_seg_loop_cir(int i, int nloops_i,
                        const std::vector<std::vector<int32_t>>& loops_i,
                        const std::vector<int>& loopsize_i, const Geom& geom,
                        const DataI& di, const DataCir& dc, const DataSeg& ds,
                        std::vector<std::array<double, 12>>& segment0,
                        std::vector<double>& Rj) {
  const std::vector<Vec3>& I = di.I;
  const std::vector<std::array<double, 10>>& circle = dc.circle;
  const std::vector<std::array<int32_t, 5>>& segment = ds.segment;
  const std::vector<std::array<double, 8>>& ncrasegment = ds.ncrasegment;
  const std::vector<int32_t>& satom_i = ds.satom[static_cast<std::size_t>(i)];
  const std::vector<int32_t>& circleindex_i = dc.circleindex[static_cast<std::size_t>(i)];
  const std::vector<double>& R = geom.R;

  const int nsatom_i = static_cast<int>(satom_i.size());
  const int ncircleindex_i = static_cast<int>(circleindex_i.size());

  ModOut out;

  // loops_i0 starts as a copy of loops_i (entries = global segment indices).
  out.loops_i0 = loops_i;

  // circle0: dummy row 0; data columns 1..7 = c(3),n(3),r; col 8 = torus partner R.
  out.circle0.assign(static_cast<std::size_t>(ncircleindex_i) + 1, std::array<double, 9>{});

  // --- modify the segment's direction (build into segvect, then scatter) ----
  std::vector<std::array<double, 12>> segvect(static_cast<std::size_t>(nsatom_i) + 1,
                                              std::array<double, 12>{});
  for (int k = 1; k <= nsatom_i; ++k) {
    const int sa = satom_i[static_cast<std::size_t>(k - 1)];
    const auto& seg = segment[static_cast<std::size_t>(sa)];   // {i,j,p1,p2,direct}
    const auto& ncra = ncrasegment[static_cast<std::size_t>(sa)];  // {n(3),A(3),r,radian}

    const int seg_i = seg[0];   // segment[sa,1]
    const int seg_j = seg[1];   // segment[sa,2]
    const int seg_p1 = seg[2];  // segment[sa,3]
    const int seg_p2 = seg[3];  // segment[sa,4]
    const int seg_direct = seg[4];  // segment[sa,5]

    // direct_s == -1 means the right hand of the segment is interior.
    const int direct_s = seg_direct * (2 * (seg_i == i ? 1 : 0) - 1);
    // direct_n == -1 means the normal points inside.
    const double direct_n = static_cast<double>(2 * (seg_i == i ? 1 : 0) - 1);

    const Vec3 n{ncra[0], ncra[1], ncra[2]};  // ncra[1:4]
    const Vec3 A{ncra[3], ncra[4], ncra[5]};  // ncra[4:7]
    const double r = ncra[6];                 // ncra[7]
    const double radian = ncra[7];            // ncra[8]

    auto& sv = segvect[static_cast<std::size_t>(k)];
    if (direct_s == -1) {
      // use endpoint p2 as start point: P1 = I[segment[sa,4]]
      const Vec3 P1 = I[static_cast<std::size_t>(seg_p2)];
      rset3(sv, 1, A);
      rset3(sv, 4, direct_n * n);
      sv[7] = r;
      rset3(sv, 8, P1);
      sv[11] = radian;
    } else {
      const Vec3 P2 = I[static_cast<std::size_t>(seg_p1)];  // I[segment[sa,3]]
      rset3(sv, 1, A);
      rset3(sv, 4, direct_n * n);
      sv[7] = r;
      rset3(sv, 8, P2);
      sv[11] = radian;
    }

    if (seg_i == i) {
      Rj[static_cast<std::size_t>(sa)] = R[static_cast<std::size_t>(seg_j)];
    } else {
      Rj[static_cast<std::size_t>(sa)] = R[static_cast<std::size_t>(seg_i)];
    }
  }

  // scatter the modified segments into the global segment0 at their indices.
  for (int k = 1; k <= nsatom_i; ++k) {
    const int sa = satom_i[static_cast<std::size_t>(k - 1)];
    segment0[static_cast<std::size_t>(sa)] = segvect[static_cast<std::size_t>(k)];
  }

  // --- modify the loop's direction (reverse each loop's segment order) -------
  for (int k = 1; k <= nloops_i; ++k) {
    const int ls = loopsize_i[static_cast<std::size_t>(k)];
    std::vector<int32_t>& dst = out.loops_i0[static_cast<std::size_t>(k)];
    const std::vector<int32_t>& src = loops_i[static_cast<std::size_t>(k)];
    // loops_i0[k, 1:ls+1] = loops_i[k, 1:ls+1][::-1]
    for (int t = 1; t <= ls; ++t) {
      dst[static_cast<std::size_t>(t)] = src[static_cast<std::size_t>(ls - t + 1)];
    }
  }

  // --- modify the circle's direction ----------------------------------------
  for (int k = 1; k <= ncircleindex_i; ++k) {
    const int cidx = circleindex_i[static_cast<std::size_t>(k - 1)];
    const auto& crow = circle[static_cast<std::size_t>(cidx)];  // {_, i,j, c(3), n(3), r}
    const int c_i = static_cast<int>(crow[1]);  // circle[cidx,1]
    const int c_j = static_cast<int>(crow[2]);  // circle[cidx,2]
    const Vec3 c{crow[3], crow[4], crow[5]};    // circle[cidx,3:6]
    const Vec3 cn{crow[6], crow[7], crow[8]};   // circle[cidx,6:9]
    const double cr = crow[9];                  // circle[cidx,9]

    auto& c0 = out.circle0[static_cast<std::size_t>(k)];
    if (c_i == i) {
      rset3(c0, 1, c);
      rset3(c0, 4, cn);
      c0[7] = cr;
      c0[8] = R[static_cast<std::size_t>(c_j)];  // torus partner VdW radius
    } else {
      rset3(c0, 1, c);
      rset3(c0, 4, -1.0 * cn);
      c0[7] = cr;
      c0[8] = R[static_cast<std::size_t>(c_i)];  // torus partner VdW radius
    }
  }

  return out;
}

}  // namespace

// ----- data_SESsphpat_convex (driver) ---------------------------------------
void data_SESsphpat_convex(MeshState& state, const Geom& geom, const DataI& di,
                           const DataCir& dc, const DataSeg& ds,
                           const DataLoop& dl, const DataPat& dp, const Ext* ext,
                           double Rp, double d) {
  (void)ext;  // ext only routes ext/int figures in MATLAB; both are meshed here.

  const int M = geom.M;
  const std::vector<Vec3>& C = geom.centers;
  const std::vector<double>& R = geom.R;

  const int nsegment = ds.nsegment;

  // PARALLEL S8: each atom meshes independently into its own vector<LocalMesh>
  // (one slot per patch). CRITICAL: segment0 and Rj were GLOBAL reused buffers;
  // Rj[sa] and segment0[sa] differ by which atom is meshing (a segment is shared
  // between two atoms but mod_seg_loop_cir scatters direction/spoint relative to
  // the CURRENT atom), so a shared buffer would be a data race. They are now
  // THREAD-LOCAL: allocated fresh per atom inside the parallel body. A SERIAL
  // ordered merge then add_patch each patch in atom-then-patch order.
  std::vector<std::vector<LocalMesh>> atom_lm(static_cast<std::size_t>(M) + 1);

  meshms::parallel_for(1, M + 1, [&](int i) {
    const int nsatom_i = static_cast<int>(ds.satom[static_cast<std::size_t>(i)].size());
    const int ncircleindex_i =
        static_cast<int>(dc.circleindex[static_cast<std::size_t>(i)].size());

    // --- slice loops_i / loopsize_i ----------------------------------------
    std::vector<std::vector<int32_t>> loops_i;  // 1-based-with-dummy rows; [0] dummy
    std::vector<int> loopsize_i;                // 1-based-with-dummy; [0] dummy
    int nloops_i = 0;

    if (nsatom_i > 0) {
      const int start = dl.loops_index[static_cast<std::size_t>(i)][0];
      const int end = dl.loops_index[static_cast<std::size_t>(i)][1];
      nloops_i = end - start + 1;
      loops_i.emplace_back();   // loops_i[0] dummy
      loopsize_i.push_back(0);  // loopsize_i[0] dummy
      for (int k = 1; k <= nloops_i; ++k) {
        const std::vector<int32_t>& glb = dl.loops[static_cast<std::size_t>(start + k - 1)];
        loops_i.push_back(glb);
        loopsize_i.push_back(static_cast<int>(glb.size()) - 1);  // drop [0] dummy
      }
    } else if (ncircleindex_i > 0) {
      loops_i.emplace_back();   // single dummy row ([] -> dummy)
      loopsize_i.push_back(0);
      nloops_i = 0;
    } else {
      // no patches on this atom; skip (neither MATLAB branch runs).
      return;
    }

    // nsatom_i > 0 || ncircleindex_i > 0 is guaranteed here.

    // --- slice patches_i / patchesize_i ------------------------------------
    const int pstart = dp.patches_index[static_cast<std::size_t>(i)][0];
    const int pend = dp.patches_index[static_cast<std::size_t>(i)][1];
    const int npatches_i = pend - pstart + 1;

    // --- THREAD-LOCAL segment0 / Rj (one per atom; no shared write) --------
    // segment0 : (nsegment+1)x12 record matrix, rebuilt for atom i's segments.
    std::vector<std::array<double, 12>> segment0(static_cast<std::size_t>(nsegment) + 1,
                                                 std::array<double, 12>{});
    // Rj (size nsegment+1, init 0): written by mod_seg_loop_cir, read by
    // mesh_sphpat for the convex near-cusp arc refinement.
    std::vector<double> Rj(static_cast<std::size_t>(nsegment) + 1, 0.0);

    // --- rebuild loops_i0 / segment0 / circle0 and write Rj ----------------
    ModOut mo = mod_seg_loop_cir(i, nloops_i, loops_i, loopsize_i, geom, di, dc,
                                 ds, segment0, Rj);

    // mesh_sphpat derives loopsize from loops[k].size()-1; trim each loop to
    // loopsize_i[k]+1 so phantom dummy (segment0[0]) entries are not meshed.
    std::vector<Loop> loops_mesh(static_cast<std::size_t>(nloops_i) + 1);
    loops_mesh[0] = Loop(1, 0);  // dummy row 0
    for (int k = 1; k <= nloops_i; ++k) {
      const int ls = loopsize_i[static_cast<std::size_t>(k)];
      Loop lk(static_cast<std::size_t>(ls) + 1, 0);  // [0] dummy + ls entries
      for (int t = 1; t <= ls; ++t) {
        lk[static_cast<std::size_t>(t)] =
            mo.loops_i0[static_cast<std::size_t>(k)][static_cast<std::size_t>(t)];
      }
      loops_mesh[static_cast<std::size_t>(k)] = std::move(lk);
    }

    // --- mesh every spherical patch on atom i ------------------------------
    std::vector<LocalMesh>& out_i = atom_lm[static_cast<std::size_t>(i)];
    out_i.reserve(static_cast<std::size_t>(npatches_i));
    const Tag boundary_tag{1, i, 0};  // ("atom", i)
    for (int j = 1; j <= npatches_i; ++j) {
      const std::vector<int32_t>& prow = dp.patches[static_cast<std::size_t>(pstart + j - 1)];
      const int patchesize_i = static_cast<int>(prow.size()) - 1;  // drop [0] dummy
      // mesh_sphpat wants patches as std::vector<int> (1-based-with-dummy).
      std::vector<int> patches_j(prow.begin(), prow.end());
      // Convex SAS-patch on atom i: sphere centre C[i], radius R[i]+Rp.
      // mesh_sphpat applies the SAS->SES contraction and uses Rj for the convex
      // near-cusp arc refinement. Add only when emit (old r_sphere!=0 && Nt>0).
      out_i.push_back(mesh_sphpat(
          C[static_cast<std::size_t>(i)], R[static_cast<std::size_t>(i)] + Rp,
          loops_mesh, segment0, mo.circle0, patches_j, patchesize_i, Rp, d, &Rj,
          boundary_tag));
      // The whole convex patch lies on atom i: every vertex is owned by i
      // (same atom id as boundary_tag = ("atom", i)). vatom is aligned with
      // P[1..Np] and never affects V/F/N.
      LocalMesh& lm = out_i.back();
      const std::size_t Np = lm.P.empty() ? 0 : lm.P.size() - 1;
      lm.vatom.assign(Np, static_cast<int32_t>(i));
    }
  });

  // --- SERIAL ordered merge (atom ascending, then patch ascending) -------
  // Pre-sum the emitted sizes so the accumulator reserves once (capacity-only),
  // and move the per-vertex tag lists out of the dead LocalMesh.
  std::size_t add_v = 0, add_f = 0;
  for (int i = 1; i <= M; ++i) {
    for (const LocalMesh& lm : atom_lm[static_cast<std::size_t>(i)]) {
      if (lm.emit) {
        add_v += lm.P.empty() ? 0 : lm.P.size() - 1;
        add_f += lm.T.size();
      }
    }
  }
  state.reserve_extra(add_v, add_f);
  for (int i = 1; i <= M; ++i) {
    for (LocalMesh& lm : atom_lm[static_cast<std::size_t>(i)]) {
      if (lm.emit) state.add_patch(lm.P, lm.T, lm.NV, std::move(lm.vids), lm.vatom);
    }
  }
}

}  // namespace meshms
