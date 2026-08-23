// Faithful port of the exterior module's data_ext() (the _owner reverse-index code
// path). See exterior.hpp for the algorithm overview.
//
// The module-level neighbor_patch() helper in the exterior module is DEAD (defined but
// never called by data_ext) and is intentionally NOT ported.
#include "meshms/exterior.hpp"

#include <array>
#include <stdexcept>

namespace meshms {

Ext data_ext(const Geom& geom, const Neighbors& /*nb*/, const DataI& data_i,
             const DataCir& data_c, const DataSeg& data_seg,
             const DataLoop& data_loop, const DataPat& data_pat, double /*Rp*/) {
  // ---- unpack (mirrors DataGlob.m bindings) ----
  const int M = geom.M;
  const std::vector<Vec3>& C = geom.centers;
  const std::vector<double>& R = geom.R;

  // Row = inter.num_int is only used in Python to size the heuristic neighbor
  // matrix width; the C++ port computes the true max_boundary below, so nb is
  // unused here.

  const std::vector<Vec3>& I = data_i.I;
  const int s = data_i.nI;

  const std::vector<std::array<double, 10>>& circle = data_c.circle;
  // circleindex[a] is 0-based per-atom: Python circleindex[i, c] (1-based c) ->
  // circleindex[i][c-1].
  const std::vector<std::vector<int32_t>>& circleindex = data_c.circleindex;
  const int ncircle = data_c.ncircle;

  const std::vector<std::array<int32_t, 5>>& seg = data_seg.segment;
  const int nsegment = data_seg.nsegment;

  // loops[ln] is 1-based with [0] dummy: loopsize[ln] == loops[ln].size()-1,
  // loops[ln, t] == loops[ln][t].
  const std::vector<std::vector<int32_t>>& loops = data_loop.loops;
  const std::vector<std::array<int32_t, 2>>& loops_index = data_loop.loops_index;

  // patches[j] is 1-based with [0] dummy: patchesize[j] == patches[j].size()-1,
  // patches[j, k] == patches[j][k].
  const std::vector<std::vector<int32_t>>& patches = data_pat.patches;
  const int npatches = data_pat.npatches;
  const std::vector<std::array<int32_t, 2>>& patches_index =
      data_pat.patches_index;
  const std::vector<int32_t>& patch_atom = data_pat.patch_atom;

  // Convenience accessors mirroring the MATLAB *size globals (jagged in C++).
  auto loopsize = [&](int ln) {
    return static_cast<int>(loops[static_cast<std::size_t>(ln)].size()) - 1;
  };
  auto patchesize = [&](int j) {
    return static_cast<int>(patches[static_cast<std::size_t>(j)].size()) - 1;
  };
  // circleindex[i, c] for 1-based c (Python uses circleindex[i, -patches[j,k]]).
  auto cidx = [&](int i, int c) {
    return circleindex[static_cast<std::size_t>(i)]
                      [static_cast<std::size_t>(c - 1)];
  };

  // ================================================================
  // neighbourship between patches
  // ================================================================
  // neighbor[j] records all neighbour patches to the j-th spherical patch.
  // Size the width to the true maximum per-patch boundary element count
  // (sum of loop sizes + number of circles) --- matches the exterior module.
  // Exact per-patch boundary element count (sum of loop sizes + #circles); the
  // counts drive a CSR layout below instead of the old dense
  // npatches x max_boundary rectangle (npatches allocations of the global
  // maximum width). Pure integer bookkeeping: the appended values and their
  // per-patch order are unchanged.
  std::vector<int64_t> pcnt(static_cast<std::size_t>(npatches) + 1, 0);
  for (int i = 1; i <= M; ++i) {
    if (patches_index[static_cast<std::size_t>(i)][0] > 0) {
      for (int j = patches_index[static_cast<std::size_t>(i)][0];
           j <= patches_index[static_cast<std::size_t>(i)][1]; ++j) {
        int cnt = 0;
        for (int k = 1; k <= patchesize(j); ++k) {
          if (patches[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)] >
              0) {
            const int ln =
                loops_index[static_cast<std::size_t>(i)][0] +
                patches[static_cast<std::size_t>(j)][static_cast<std::size_t>(
                    k)] -
                1;
            cnt += loopsize(ln);
          } else {
            cnt += 1;
          }
        }
        pcnt[static_cast<std::size_t>(j)] = cnt;
      }
    }
  }
  // CSR: patch j's neighbour slots live at nval[noff[j] .. noff[j]+Nneighbor[j]-1].
  std::vector<int64_t> noff(static_cast<std::size_t>(npatches) + 2, 0);
  for (int j = 1; j <= npatches; ++j) {
    noff[static_cast<std::size_t>(j) + 1] =
        noff[static_cast<std::size_t>(j)] + pcnt[static_cast<std::size_t>(j)];
  }
  std::vector<int32_t> nval(static_cast<std::size_t>(noff[static_cast<std::size_t>(npatches) + 1]), 0);
  std::vector<int64_t> Nneighbor(static_cast<std::size_t>(npatches) + 1, 0);

  // Reverse index: per (id, atom) -> the single patch on that atom owning the
  // given segment/circle id. A segment lies on exactly the two atoms seg[sn][0]
  // and seg[sn][1], a circle on circle[cid][1] and circle[cid][2], so the old
  // map-of-maps (one heap hash map per id, <=2 entries each) flattens to one
  // 2-slot row per id. _owner only ever queries those two atoms, so a write by
  // any other atom (never queried, dead in the map too) is dropped. Same
  // last-writer-wins overwrite semantics, zero hashing.
  std::vector<std::array<int32_t, 2>> seg_owner(
      static_cast<std::size_t>(nsegment) + 1, std::array<int32_t, 2>{0, 0});
  std::vector<std::array<int32_t, 2>> cir_owner(
      static_cast<std::size_t>(ncircle) + 1, std::array<int32_t, 2>{0, 0});
  auto seg_slot = [&](int sn, int atom) -> int {
    if (atom == seg[static_cast<std::size_t>(sn)][0]) return 0;
    if (atom == seg[static_cast<std::size_t>(sn)][1]) return 1;
    return -1;
  };
  auto cir_slot = [&](int cid, int atom) -> int {
    if (atom == static_cast<int>(circle[static_cast<std::size_t>(cid)][1])) return 0;
    if (atom == static_cast<int>(circle[static_cast<std::size_t>(cid)][2])) return 1;
    return -1;
  };
  for (int i = 1; i <= M; ++i) {
    if (patches_index[static_cast<std::size_t>(i)][0] > 0) {
      for (int j = patches_index[static_cast<std::size_t>(i)][0];
           j <= patches_index[static_cast<std::size_t>(i)][1]; ++j) {
        for (int k = 1; k <= patchesize(j); ++k) {
          const int pjk =
              patches[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)];
          if (pjk > 0) {
            const int ln =
                loops_index[static_cast<std::size_t>(i)][0] + pjk - 1;
            for (int s0 = 1; s0 <= loopsize(ln); ++s0) {
              const int sn =
                  loops[static_cast<std::size_t>(ln)][static_cast<std::size_t>(
                      s0)];
              const int sl = seg_slot(sn, i);
              if (sl >= 0) seg_owner[static_cast<std::size_t>(sn)][static_cast<std::size_t>(sl)] = j;
            }
          } else {
            const int cid = cidx(i, -pjk);
            const int sl = cir_slot(cid, i);
            if (sl >= 0) cir_owner[static_cast<std::size_t>(cid)][static_cast<std::size_t>(sl)] = j;
          }
        }
      }
    }
  }

  // Patch on the OTHER sphere (i0) owning signed id n (n>0 segment, n<0 circle).
  // Falls back to patches_index[i0,2] (last patch on i0) when not found, matching
  // neighbor_patch's fall-through (which leaves j == patches_index[i0,2]).
  auto _owner = [&](int n, int sp0, int sp1, int i) -> int {
    int i0;
    if (sp0 == i) {
      i0 = sp1;
    } else if (sp1 == i) {
      i0 = sp0;
    } else {
      throw std::runtime_error("error");
    }
    const int key = (n > 0) ? n : -n;
    const int sl = (n > 0) ? seg_slot(key, i0) : cir_slot(key, i0);
    if (sl >= 0) {
      const int32_t owner = (n > 0)
          ? seg_owner[static_cast<std::size_t>(key)][static_cast<std::size_t>(sl)]
          : cir_owner[static_cast<std::size_t>(key)][static_cast<std::size_t>(sl)];
      if (owner != 0) return owner;  // 0 == never written (map-miss fallback)
    }
    return patches_index[static_cast<std::size_t>(i0)][1];
  };

  for (int i = 1; i <= M; ++i) {
    if (patches_index[static_cast<std::size_t>(i)][0] > 0) {  // i-th sphere
      for (int j = patches_index[static_cast<std::size_t>(i)][0];
           j <= patches_index[static_cast<std::size_t>(i)][1]; ++j) {
        // j denotes the j-th spherical patch
        for (int k = 1; k <= patchesize(j); ++k) {
          const int pjk =
              patches[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)];
          if (pjk > 0) {
            // patches(j,k)>0 -> index of loop on the i-th sphere
            const int ln =
                loops_index[static_cast<std::size_t>(i)][0] + pjk - 1;
            for (int s0 = 1; s0 <= loopsize(ln); ++s0) {
              const int sn =
                  loops[static_cast<std::size_t>(ln)][static_cast<std::size_t>(
                      s0)];
              // seg[sn, 1:3] == (seg(sn,1), seg(sn,2)) == segment[sn][0..1]
              const int j0 =
                  _owner(sn, seg[static_cast<std::size_t>(sn)][0],
                         seg[static_cast<std::size_t>(sn)][1], i);
              nval[static_cast<std::size_t>(
                  noff[static_cast<std::size_t>(j)] +
                  Nneighbor[static_cast<std::size_t>(j)])] = j0;
              Nneighbor[static_cast<std::size_t>(j)] += 1;
            }
          } else {
            // patches(j,k)<0 -> index of circle on the i-th sphere
            const int cn = -cidx(i, -pjk);  // circle number (negative)
            // circle[-cn, 1:3] == (circle(-cn,1), circle(-cn,2))
            const int sp0 = static_cast<int>(
                circle[static_cast<std::size_t>(-cn)][1]);
            const int sp1 = static_cast<int>(
                circle[static_cast<std::size_t>(-cn)][2]);
            const int j0 = _owner(cn, sp0, sp1, i);
            nval[static_cast<std::size_t>(
                noff[static_cast<std::size_t>(j)] +
                Nneighbor[static_cast<std::size_t>(j)])] = j0;
            Nneighbor[static_cast<std::size_t>(j)] += 1;
          }
        }
      }
    }
  }

  // ================================================================
  // ext_patch, ext_segment, ext_I
  // ================================================================
  // obtain an initial patch outside: the leftmost sphere.
  // [~,i_left] = min(C(1:M,1)-R); -> argmin over atoms 1..M of x - R, numpy
  // first-min tie-break (strict '<', keep first).
  int i_left = 1;
  {
    double best = C[1].x - R[1];
    for (int i = 2; i <= M; ++i) {
      const double v = C[static_cast<std::size_t>(i)].x - R[static_cast<std::size_t>(i)];
      if (v < best) {
        best = v;
        i_left = i;
      }
    }
  }
  int t = 1;
  if (patches_index[static_cast<std::size_t>(i_left)][0] == 0) {
    // secondary argmin on .y, same first-min tie-break.
    i_left = 1;
    double best = C[1].y - R[1];
    for (int i = 2; i <= M; ++i) {
      const double v = C[static_cast<std::size_t>(i)].y - R[static_cast<std::size_t>(i)];
      if (v < best) {
        best = v;
        i_left = i;
      }
    }
  }

  double x_left = 0.0;  // set on first encountered boundary point (t toggles)
  int j_left = 0;
  for (int j = patches_index[static_cast<std::size_t>(i_left)][0];
       j <= patches_index[static_cast<std::size_t>(i_left)][1]; ++j) {
    if (j == patches_index[static_cast<std::size_t>(i_left)][0]) {
      j_left = j;

      if (j == 0) {
        throw std::runtime_error("Error: there exists isolated SAS-ball!");
      }

      for (int k = 1; k <= patchesize(j); ++k) {
        const int pjk =
            patches[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)];
        if (pjk > 0) {
          const int ln =
              loops_index[static_cast<std::size_t>(i_left)][0] + pjk - 1;
          for (int s0 = 1; s0 <= loopsize(ln); ++s0) {
            const int sn =
                loops[static_cast<std::size_t>(ln)][static_cast<std::size_t>(
                    s0)];
            // I[seg[sn,3], 0] == I[ segment[sn][2] ].x
            const double xv =
                I[static_cast<std::size_t>(seg[static_cast<std::size_t>(sn)][2])]
                    .x;
            if (t == 1) {
              x_left = xv;
              t = 0;
            } else if (x_left > xv) {
              x_left = xv;
            }
          }
        } else {
          const int cn = -cidx(i_left, -pjk);  // circle number (negative)
          // circle[-cn, 3] == Ax (center x)
          const double xv = circle[static_cast<std::size_t>(-cn)][3];
          if (t == 1) {
            x_left = xv;
            t = 0;
          } else if (x_left > xv) {
            x_left = xv;
          }
        }
      }
    } else {
      for (int k = 1; k <= patchesize(j); ++k) {
        const int pjk =
            patches[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)];
        if (pjk > 0) {
          const int ln =
              loops_index[static_cast<std::size_t>(i_left)][0] + pjk - 1;
          for (int s0 = 1; s0 <= loopsize(ln); ++s0) {
            const int sn =
                loops[static_cast<std::size_t>(ln)][static_cast<std::size_t>(
                    s0)];
            const double xv =
                I[static_cast<std::size_t>(seg[static_cast<std::size_t>(sn)][2])]
                    .x;
            if (x_left > xv) {
              j_left = j;
              x_left = xv;
            }
          }
        } else {
          const int cn = -cidx(i_left, -pjk);  // circle number (negative)
          const double xv = circle[static_cast<std::size_t>(-cn)][3];
          if (x_left > xv) {
            j_left = j;
            x_left = xv;
          }
        }
      }
    }
  }

  // ext_patch(i) = 1 if the i-th spherical patch is on the eSAS
  std::vector<int32_t> ext_patch(static_cast<std::size_t>(npatches) + 1, 0);
  ext_patch[static_cast<std::size_t>(j_left)] = 1;

  // ext_patchset = the set of spherical patches on the eSAS
  std::vector<int32_t> ext_patchset(static_cast<std::size_t>(npatches) + 1, 0);
  ext_patchset[1] = j_left;
  int nset = 1;

  std::vector<int32_t> ext_segment(static_cast<std::size_t>(nsegment) + 1, 0);
  std::vector<int32_t> ext_circle(static_cast<std::size_t>(ncircle) + 1, 0);

  // flood-fill: walk the discovered exterior patches and add new neighbours.
  for (int i = 1; i <= npatches; ++i) {
    if (i <= nset) {
      const int j = ext_patchset[static_cast<std::size_t>(i)];
      for (int k = 1; k <= Nneighbor[static_cast<std::size_t>(j)]; ++k) {
        const int n = nval[static_cast<std::size_t>(
            noff[static_cast<std::size_t>(j)] + k - 1)];
        if (ext_patch[static_cast<std::size_t>(n)] == 0) {
          nset += 1;
          ext_patchset[static_cast<std::size_t>(nset)] = n;
          ext_patch[static_cast<std::size_t>(n)] = 1;
        }
      }
    } else {
      break;
    }
  }

  // mark exterior segments / circles
  for (int i = 1; i <= nset; ++i) {
    const int j = ext_patchset[static_cast<std::size_t>(i)];  // j-th patch
    const int pa = patch_atom[static_cast<std::size_t>(j)];
    for (int k = 1; k <= patchesize(j); ++k) {
      const int pjk =
          patches[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)];
      if (pjk > 0) {
        const int ln = loops_index[static_cast<std::size_t>(pa)][0] + pjk - 1;
        for (int s0 = 1; s0 <= loopsize(ln); ++s0) {
          const int sn =
              loops[static_cast<std::size_t>(ln)][static_cast<std::size_t>(s0)];
          if (ext_segment[static_cast<std::size_t>(sn)] == 0) {
            ext_segment[static_cast<std::size_t>(sn)] = 1;
          }
        }
      } else {
        const int cn = -cidx(pa, -pjk);  // circle number (negative)
        if (ext_circle[static_cast<std::size_t>(-cn)] == 0) {
          ext_circle[static_cast<std::size_t>(-cn)] = 1;
        }
      }
    }
  }

  // ext_I(i) = 1 if the i-th intersection point is on the eSAS
  std::vector<int32_t> ext_I(static_cast<std::size_t>(s) + 1, 0);
  for (int i = 1; i <= nsegment; ++i) {
    if (ext_segment[static_cast<std::size_t>(i)] == 1) {
      int j = seg[static_cast<std::size_t>(i)][2];  // seg[i,3] == p1
      if (ext_I[static_cast<std::size_t>(j)] == 0) ext_I[static_cast<std::size_t>(j)] = 1;
      j = seg[static_cast<std::size_t>(i)][3];  // seg[i,4] == p2
      if (ext_I[static_cast<std::size_t>(j)] == 0) ext_I[static_cast<std::size_t>(j)] = 1;
    }
  }

  Ext out;
  out.I = std::move(ext_I);
  out.circle = std::move(ext_circle);
  out.segment = std::move(ext_segment);
  out.patch = std::move(ext_patch);
  return out;
}

}  // namespace meshms
