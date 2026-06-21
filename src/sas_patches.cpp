// Faithful port of the sas_patches module's data_Seg_Pat() (MESH PATH, want_area
// =False) and its helpers interiorloop / loopconstruct / patchesconstruct /
// TreeNode. The float term/eval order, np.argsort(kind="stable") ->
// std::stable_sort, np.clip-then-arccos and pysq() (Python scalar `x ** 2`) are
// reproduced exactly so the per-decision integer topology (segments, loops,
// patches) matches the Python golden to the branch.
//
// The want_area block (mod_seg_loop_cir + area_spherical + DataAV + segment0 +
// Area_sphpat) is DEAD in the mesh path and is DEFERRED (see the TODO below).
//
// Indexing: 1-based with a dummy entry at index 0 throughout, matching the
// Python record matrices (dummy row 0). Coordinate vectors I[] (data_i.I) and
// centers C (geom.centers) also use a [0] dummy. The CSR I_circle stores its
// per-(i,row) point ids 0-based without a dummy (di.I_circle[i][row][t] for
// t = 0..count-1), so Python I_circle[i,row,1] (start point) == [i][row][0].
#include "meshms/sas_patches.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>

#include "meshms/vec3.hpp"

namespace meshms {

namespace {

// ----------------------------------------------------------------------------
// TreeNode (port of Classes/treenode.m). `set` is a 1-based-by-Python-index but
// here a plain 0-based std::vector<int32_t> of signed loop/circle entries; the
// Python S[s1-1] access maps to set[s1-1].
// ----------------------------------------------------------------------------
struct TreeNode {
  int activenode{0};
  std::vector<int32_t> set;
  int activeelement{0};
  int n1{0};
  int n2{0};
  int leftnode{0};
  int rightnode{0};
};

// Local alpha(direct,u,v,n) == arc_angle(u,v,n,direct).
inline double alpha_local(int direct, const Vec3& u, const Vec3& v, const Vec3& n) {
  return arc_angle(u, v, n, direct);
}

// ----------------------------------------------------------------------------
// interiorloop (port of Core/interiorloop.m): is `point` inside the loop?
//
// segment / ncrasegment use the C++ DataSeg layout: segment[k] = {i,j,p1,p2,
// direct} (Python cols 1..5), ncrasegment[k] = {nx,ny,nz, Ax,Ay,Az, r, radian}
// (Python cols 1..8). loop is the 1-based loop vector (loop[1..loopsize] = global
// segment ids; loop[0] dummy). center == geom.centers, I == data_i.I.
// ----------------------------------------------------------------------------
int interiorloop(const Vec3& point, const Vec3& ci, const std::vector<Vec3>& center,
                 double ri, int i, const std::vector<int32_t>& loop, int loopsize,
                 const std::vector<Vec3>& I,
                 const std::vector<std::array<int32_t, 5>>& segment,
                 const std::vector<std::array<double, 8>>& ncrasegment) {
  int nearest = 0;
  double theta0 = 0.0;
  for (int k = 1; k <= loopsize; ++k) {
    const auto& seg = segment[static_cast<std::size_t>(loop[static_cast<std::size_t>(k)])];
    int j;
    if (seg[0] == i) {  // segment[loop[k], 1] == i
      j = seg[1];
    } else {
      j = seg[0];
    }

    Vec3 v1 = (center[static_cast<std::size_t>(j)] - ci) /
                  norm(center[static_cast<std::size_t>(j)] - ci) * ri;
    Vec3 v2 = I[static_cast<std::size_t>(seg[2])] - ci;  // I[segment[loop[k], 3]] - ci
    Vec3 v3 = point - ci;
    // theta = arccos(clip(dot(v3,v1)/ri**2)) - arccos(clip(dot(v2,v1)/ri**2))
    double theta = acos_clamped(dot(v3, v1) / pysq(ri)) -
                   acos_clamped(dot(v2, v1) / pysq(ri));

    if (k == 1 || theta < theta0) {
      theta0 = theta;
      nearest = k;
    }
  }

  // K[] collects the loop entries whose "other atom" j == j_nearest.
  std::vector<int> K(static_cast<std::size_t>(loopsize) + 1, 0);  // [0] dummy
  int nK = 0;
  int j_nearest;
  {
    const auto& segn = segment[static_cast<std::size_t>(loop[static_cast<std::size_t>(nearest)])];
    if (segn[0] == i) {
      j_nearest = segn[1];
    } else {
      j_nearest = segn[0];
    }
  }

  for (int k = 1; k <= loopsize; ++k) {
    const auto& seg = segment[static_cast<std::size_t>(loop[static_cast<std::size_t>(k)])];
    int j;
    if (seg[0] == i) {
      j = seg[1];
    } else {
      j = seg[0];
    }
    if (j == j_nearest) {
      ++nK;
      K[static_cast<std::size_t>(nK)] = k;
    }
  }

  int j = j_nearest;

  Vec3 v1 = (center[static_cast<std::size_t>(j)] - ci) /
                norm(center[static_cast<std::size_t>(j)] - ci);
  Vec3 v3 = point - ci;
  Vec3 v = v3 - dot(v1, v3) * v1;

  for (int k0 = 1; k0 <= nK; ++k0) {
    int k = K[static_cast<std::size_t>(k0)];
    const auto& seg = segment[static_cast<std::size_t>(loop[static_cast<std::size_t>(k)])];
    const auto& ncra = ncrasegment[static_cast<std::size_t>(loop[static_cast<std::size_t>(k)])];

    int direct = seg[4];                         // segment[loop[k], 5]
    Vec3 n{ncra[0], ncra[1], ncra[2]};           // ncrasegment[loop[k], 1:4]
    Vec3 A{ncra[3], ncra[4], ncra[5]};           // ncrasegment[loop[k], 4:7]
    Vec3 p3 = I[static_cast<std::size_t>(seg[2])];   // I[segment[loop[k], 3]]
    Vec3 p4 = I[static_cast<std::size_t>(seg[3])];   // I[segment[loop[k], 4]]
    double alpha1 = alpha_local(direct, p3 - A, v, n);
    double alpha2 = alpha_local(direct, p3 - A, p4 - A, n);

    if (alpha1 < alpha2) {
      return 1;
    }
  }

  return 0;
}

// ----------------------------------------------------------------------------
// loopconstruct (local helper of data_Seg_Pat.m): build all loops on sphere i.
//
// sa is the 1-based segment id list on atom i (sa[1..ns], sa[0] dummy). Returns
// loops (loops[k] = global segment ids; loops[0] dummy), loopsize implicit in
// loops[k].size(), and nloops.
// ----------------------------------------------------------------------------
struct LoopResult {
  std::vector<std::vector<int32_t>> loops;  // loops[k] (1-based; [0] dummy)
  int nloops{0};
};

LoopResult loopconstruct(int i, const std::vector<int32_t>& sa, int ns,
                         const std::vector<std::array<int32_t, 5>>& segment) {
  // s[k] = directed endpoints (p_start, p_end) of segment sa[k]; [0] dummy in k.
  std::vector<std::array<int, 2>> s(static_cast<std::size_t>(ns) + 1, {0, 0});
  for (int k = 1; k <= ns; ++k) {
    int sak = sa[static_cast<std::size_t>(k)];
    const auto& seg = segment[static_cast<std::size_t>(sak)];
    // direction so the left hand is interior:
    // segment[sak,5] * (2*(segment[sak,1]==i) - 1) == 1
    if (seg[4] * (2 * (seg[0] == i ? 1 : 0) - 1) == 1) {
      s[static_cast<std::size_t>(k)][0] = seg[3];  // segment[sak, 4]
      s[static_cast<std::size_t>(k)][1] = seg[2];  // segment[sak, 3]
    } else {
      s[static_cast<std::size_t>(k)][0] = seg[2];
      s[static_cast<std::size_t>(k)][1] = seg[3];
    }
  }

  std::vector<int> used(static_cast<std::size_t>(ns) + 1, 0);  // Python `true[]`
  int number = 0;
  int nloops = 0;

  LoopResult res;
  res.loops.emplace_back();  // loops[0] dummy

  for (int k1 = 1; k1 <= ns; ++k1) {
    std::vector<int32_t> loop;  // 0-based collected ids; we prepend a dummy at end
    if (used[static_cast<std::size_t>(k1)] == 0 && number < ns) {
      ++nloops;
      int pointer = k1;
      ++number;
      loop.push_back(sa[static_cast<std::size_t>(pointer)]);
      used[static_cast<std::size_t>(k1)] = 1;

      while (s[static_cast<std::size_t>(k1)][0] != s[static_cast<std::size_t>(pointer)][1]) {
        bool advanced = false;
        for (int k2 = k1 + 1; k2 <= ns; ++k2) {
          if (used[static_cast<std::size_t>(k2)] == 0 &&
              s[static_cast<std::size_t>(k2)][0] == s[static_cast<std::size_t>(pointer)][1]) {
            pointer = k2;
            loop.push_back(sa[static_cast<std::size_t>(pointer)]);
            ++number;
            used[static_cast<std::size_t>(pointer)] = 1;
            advanced = true;
          }
        }
        if (!advanced) {
          // No segment continues the chain (degenerate tangencies): break to
          // avoid an infinite loop (MATLAB would spin forever here).
          break;
        }
      }

      // Store as 1-based with [0] dummy: loops[nloops][0] dummy, [1..] = loop.
      std::vector<int32_t> stored;
      stored.push_back(0);  // dummy
      for (int32_t id : loop) stored.push_back(id);
      res.loops.push_back(std::move(stored));
    }
  }

  res.nloops = nloops;
  return res;
}

// ----------------------------------------------------------------------------
// patchesconstruct (local helper of data_Seg_Pat.m): binary-tree classification
// of loops (and whole circles) into patch boundaries via interiorloop.
//
// loops_i[k] is the 1-based loop vector (loops_i[k][0] dummy); circleindex_i is
// the 1-based circleindex row on atom i (circleindex_i[1..ncircleindex_i],
// circleindex_i[0] dummy). Returns patches (patches[k] signed entries; [0]
// dummy) and npatches.
// ----------------------------------------------------------------------------
struct PatchResult {
  std::vector<std::vector<int32_t>> patches;  // patches[k] (1-based; [0] dummy)
  int npatches{0};
};

PatchResult patchesconstruct(int i, const std::vector<Vec3>& C,
                             const std::vector<double>& R,
                             const std::vector<std::array<int32_t, 5>>& segment,
                             const std::vector<std::array<double, 8>>& ncrasegment,
                             const std::vector<Vec3>& I,
                             const std::vector<std::array<double, 10>>& circle,
                             const std::vector<int32_t>& circleindex_i,
                             int ncircleindex_i,
                             const std::vector<std::vector<int32_t>>& loops_i,
                             int nloops_i, double Rp) {
  PatchResult res;
  res.patches.emplace_back();  // patches[0] dummy
  int npatches = 0;

  // S0 = [1, 2, ..., nloops_i, -1, -2, ..., -ncircleindex_i]
  std::vector<int32_t> S0;
  for (int t = 1; t <= nloops_i; ++t) S0.push_back(t);
  for (int t = 1; t <= ncircleindex_i; ++t) S0.push_back(-t);

  // Tree workspace (1-based; index 0 dummy). Grows as needed.
  std::vector<TreeNode> tree;
  tree.emplace_back();  // dummy slot at index 0
  tree.emplace_back();  // tree[1]
  tree[1].activenode = 1;
  tree[1].set = S0;
  tree[1].activeelement = S0.empty() ? 0 : S0[0];
  tree[1].n1 = nloops_i;
  tree[1].n2 = ncircleindex_i;

  int ntree = 1;
  int j = 1;

  // Python: for _s in range(1, 2*(nloops_i+ncircleindex_i) + 1 + 1) -> _s = 1..2N+1
  // (range upper is exclusive), i.e. exactly 2N+1 iterations.
  const int iters = 2 * (nloops_i + ncircleindex_i) + 1;
  for (int _s = 1; _s <= iters; ++_s) {
    if (j > ntree) break;  // the set has been completely divided

    if (j <= ntree && tree[static_cast<std::size_t>(j)].activenode == 1) {
      std::vector<int32_t> S1, S2;
      int k1 = 0, k2 = 0;
      int t1 = 1, t2 = 1;
      int left_n1 = 0, left_n2 = 0, right_n1 = 0, right_n2 = 0;

      const std::vector<int32_t>& S = tree[static_cast<std::size_t>(j)].set;
      int k = tree[static_cast<std::size_t>(j)].activeelement;  // the fixed loop

      if (static_cast<int>(S.size()) == 1) {
        tree[static_cast<std::size_t>(j)].activenode = 0;
        ++j;
        continue;
      }

      if (k > 0) {
        const std::vector<int32_t>& loopk = loops_i[static_cast<std::size_t>(k)];
        const int loopsize_k = static_cast<int>(loopk.size()) - 1;  // drop [0] dummy
        const Vec3 ci = C[static_cast<std::size_t>(i)];
        const double ri = R[static_cast<std::size_t>(i)] + Rp;

        // loop entries (1..n1)
        for (int s1 = 1; s1 <= tree[static_cast<std::size_t>(j)].n1; ++s1) {
          int Ss1 = S[static_cast<std::size_t>(s1 - 1)];
          // point = I[segment[loops_i[Ss1, 1], 3]]
          int seg0 = loops_i[static_cast<std::size_t>(Ss1)][1];  // loops_i[Ss1, 1]
          Vec3 point = I[static_cast<std::size_t>(segment[static_cast<std::size_t>(seg0)][2])];

          if (k == Ss1 ||
              interiorloop(point, ci, C, ri, i, loopk, loopsize_k, I, segment,
                           ncrasegment)) {
            S1.push_back(Ss1);
            ++left_n1;
            if (t1 == 1 && Ss1 > k) {
              t1 = 0;
              k1 = Ss1;
            }
          } else {
            S2.push_back(Ss1);
            ++right_n1;
            if (t2 == 1 && Ss1 > k) {
              t2 = 0;
              k2 = Ss1;
            }
          }
        }

        // circle entries (n1+1 .. n1+n2)
        for (int s2 = tree[static_cast<std::size_t>(j)].n1 + 1;
             s2 <= tree[static_cast<std::size_t>(j)].n1 + tree[static_cast<std::size_t>(j)].n2;
             ++s2) {
          int Ss2 = S[static_cast<std::size_t>(s2 - 1)];
          int cidx = circleindex_i[static_cast<std::size_t>(-Ss2)];  // circleindex_i[-Ss2]
          const auto& crow = circle[static_cast<std::size_t>(cidx)];
          Vec3 cnorm{crow[6], crow[7], crow[8]};   // circle[cidx, 6:9]
          auto [vector1, v2unused] = orthogonalvectors(cnorm);
          (void)v2unused;
          Vec3 ccenter{crow[3], crow[4], crow[5]};  // circle[cidx, 3:6]
          double cr = crow[9];                      // circle[cidx, 9]
          Vec3 point = ccenter + cr * vector1;

          if (interiorloop(point, ci, C, ri, i, loopk, loopsize_k, I, segment,
                           ncrasegment)) {
            S1.push_back(Ss2);
            ++left_n2;
            if (t1 == 1) {
              t1 = 0;
              k1 = Ss2;
            }
          } else {
            S2.push_back(Ss2);
            ++right_n2;
            if (t2 == 1) {
              t2 = 0;
              k2 = Ss2;
            }
          }
        }

        if (S2.empty()) {
          if (t1 == 1) {
            tree[static_cast<std::size_t>(j)].activenode = 0;
            ++j;
          } else {
            tree[static_cast<std::size_t>(j)].activeelement = k1;
          }
        } else {
          tree[static_cast<std::size_t>(j)].activenode = 0;
          tree[static_cast<std::size_t>(j)].leftnode = ntree + 1;
          tree[static_cast<std::size_t>(j)].rightnode = ntree + 2;

          ++ntree;
          tree.emplace_back();
          tree[static_cast<std::size_t>(ntree)].activenode = 1;
          tree[static_cast<std::size_t>(ntree)].set = S1;
          tree[static_cast<std::size_t>(ntree)].activeelement = k1;
          tree[static_cast<std::size_t>(ntree)].n1 = left_n1;
          tree[static_cast<std::size_t>(ntree)].n2 = left_n2;
          if (t1 == 1) tree[static_cast<std::size_t>(ntree)].activenode = 0;

          ++ntree;
          tree.emplace_back();
          tree[static_cast<std::size_t>(ntree)].activenode = 1;
          tree[static_cast<std::size_t>(ntree)].set = S2;
          tree[static_cast<std::size_t>(ntree)].activeelement = k2;
          tree[static_cast<std::size_t>(ntree)].n1 = right_n1;
          tree[static_cast<std::size_t>(ntree)].n2 = right_n2;
          if (t2 == 1) tree[static_cast<std::size_t>(ntree)].activenode = 0;

          ++j;
        }
      } else {
        tree[static_cast<std::size_t>(j)].activenode = 0;
        ++j;
      }
    } else {
      ++j;
    }
  }

  // Leaf collection: leaves (activenode==0 && leftnode==0) become patches.
  for (int jj = 1; jj <= ntree; ++jj) {
    const TreeNode& node = tree[static_cast<std::size_t>(jj)];
    if (node.activenode == 0 && node.leftnode == 0) {
      ++npatches;
      std::vector<int32_t> entry;
      entry.push_back(0);  // [0] dummy
      for (int32_t v : node.set) entry.push_back(v);
      res.patches.push_back(std::move(entry));
    }
  }

  res.npatches = npatches;
  return res;
}

}  // namespace

// ----------------------------------------------------------------------------
// data_Seg_Pat (port of Core/data_Seg_Pat.m), mesh path (want_area=false).
// ----------------------------------------------------------------------------
std::tuple<DataSeg, DataLoop, DataPat> data_Seg_Pat(const Geom& geom,
                                                    const Neighbors& nb,
                                                    const DataI& data_i,
                                                    const DataCir& data_c,
                                                    double Rp,
                                                    bool want_area) {
  if (want_area) {
    // TODO(want_area): the Gauss-Bonnet area/volume report (mod_seg_loop_cir +
    // area_spherical + DataAV + segment0 + Area_sphpat) is DEAD in the mesh path
    // and is not yet ported. Implement when the area/volume report is needed.
    throw std::runtime_error("data_Seg_Pat: want_area=true is not implemented");
  }

  const int M = geom.M;
  const std::vector<Vec3>& C = geom.centers;
  const std::vector<double>& R = geom.R;

  // Row[i] -> nb.count(i); inter.M_int[i,row] -> nb.of(i)[row-1] (1-based row).
  auto Row = [&](int a) { return nb.count(a); };
  auto Mint = [&](int a, int row) {
    return static_cast<int>(nb.of(a)[static_cast<std::size_t>(row - 1)]);
  };

  const std::vector<Vec3>& I = data_i.I;
  const auto& Iijk = data_i.Iijk;
  const auto& direction = data_i.direction;
  const auto& I_circle = data_i.I_circle;

  const auto& circle = data_c.circle;
  const auto& circleindex = data_c.circleindex;
  auto ncircleindex = [&](int a) {
    return static_cast<int>(circleindex[static_cast<std::size_t>(a)].size());
  };

  // --- Outputs (all 1-based with [0] dummy) ---------------------------------
  DataSeg ds;
  ds.segment.emplace_back(std::array<int32_t, 5>{0, 0, 0, 0, 0});      // [0] dummy
  ds.ncrasegment.emplace_back(std::array<double, 8>{});                // [0] dummy
  ds.satom.assign(static_cast<std::size_t>(M) + 1, {});                // satom[a]
  ds.Rj.assign(1, 0.0);                                                // Rj[0] dummy
  int nsegment = 0;

  DataLoop dl;
  dl.loops.emplace_back();                                             // loops[0] dummy
  dl.loops_index.assign(static_cast<std::size_t>(M) + 1, {0, 0});
  int nloops = 0;

  DataPat dp;
  dp.patches.emplace_back();                                           // patches[0] dummy
  dp.patches_index.assign(static_cast<std::size_t>(M) + 1, {0, 0});
  dp.patch_atom.assign(1, 0);                                          // patch_atom[0] dummy
  int npatches = 0;

  // nsatom counter per atom (== ds.satom[a].size()); maintained alongside.
  for (int i = 1; i <= M; ++i) {
    const int Ri = Row(i);
    for (int row = 1; row <= Ri; ++row) {
      int j = Mint(i, row);
      const std::vector<int32_t>& icirc =
          I_circle[static_cast<std::size_t>(i)][static_cast<std::size_t>(row)];
      const int npt = static_cast<int>(icirc.size());  // I_circle_num[i,row]

      if (j > i && npt > 0) {
        Vec3 A = circlecenter(C[static_cast<std::size_t>(i)], C[static_cast<std::size_t>(j)],
                              R[static_cast<std::size_t>(i)] + Rp,
                              R[static_cast<std::size_t>(j)] + Rp);
        // rij = sqrt(max((R[i]+Rp)**2 - norm(C[i]-A)**2, 0))  (scalar ** -> pysq)
        double rij = std::sqrt(std::max(
            pysq(R[static_cast<std::size_t>(i)] + Rp) -
                pysq(norm(C[static_cast<std::size_t>(i)] - A)),
            0.0));
        Vec3 nij = (C[static_cast<std::size_t>(j)] - C[static_cast<std::size_t>(i)]) /
                   norm(C[static_cast<std::size_t>(j)] - C[static_cast<std::size_t>(i)]);

        // point[] = I_circle[i,row,1:npt+1]; spoint = I_circle[i,row,1].
        // icirc is 0-based: point(Python k, 1-based) == icirc[k-1]; spoint == icirc[0].
        int spoint = icirc[0];

        // direct lookup via Iijk[spoint] / direction[spoint].
        const auto& index = Iijk[static_cast<std::size_t>(spoint)];  // {i,j,k}
        const auto& dirsp = direction[static_cast<std::size_t>(spoint)];  // {dij,djk,dki}
        int direct = 0;
        if (i == index[0] && j == index[1]) {
          direct = dirsp[0];  // direction[spoint, 1]
        } else if (i == index[1] && j == index[2]) {
          direct = dirsp[1];  // direction[spoint, 2]
        } else if (i == index[0] && j == index[2]) {
          direct = dirsp[2];  // direction[spoint, 3]
        }

        // angles relative to start point (1-based; alpha1[1] == 0 implicit).
        std::vector<double> alpha1(static_cast<std::size_t>(npt) + 1, 0.0);
        Vec3 Isp = I[static_cast<std::size_t>(spoint)];
        for (int k = 2; k <= npt; ++k) {
          // alpha(direct, I[spoint]-A, I[point[k-1]]-A, nij). Python point[] is
          // the 0-based I_circle[i,row,1:npt+1]; point[m]==icirc[m], so
          // point[k-1]==icirc[k-1].
          int pk = icirc[static_cast<std::size_t>(k - 1)];
          alpha1[static_cast<std::size_t>(k)] =
              alpha_local(direct, Isp - A, I[static_cast<std::size_t>(pk)] - A, nij);
        }

        // STABLE argsort of alpha1[1..npt]; order is 1-based (+1).
        std::vector<int> order(static_cast<std::size_t>(npt));  // 0-based indices into [1..npt]
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
          return alpha1[static_cast<std::size_t>(a + 1)] < alpha1[static_cast<std::size_t>(b + 1)];
        });

        // alpha2[1..npt] = alpha1[order]; point reordered (icirc[order]).
        std::vector<double> alpha2(static_cast<std::size_t>(npt) + 1, 0.0);
        std::vector<int> pointord(static_cast<std::size_t>(npt));  // reordered point ids (0-based)
        for (int t = 0; t < npt; ++t) {
          alpha2[static_cast<std::size_t>(t + 1)] = alpha1[static_cast<std::size_t>(order[static_cast<std::size_t>(t)] + 1)];
          pointord[static_cast<std::size_t>(t)] = icirc[static_cast<std::size_t>(order[static_cast<std::size_t>(t)])];
        }

        // pair points (2k-1, 2k) into segments.
        for (int k = 1; k <= npt / 2; ++k) {
          ++nsegment;
          int p1 = pointord[static_cast<std::size_t>(2 * k - 2)];  // point(2k-1)
          int p2 = pointord[static_cast<std::size_t>(2 * k - 1)];  // point(2k)
          ds.segment.push_back({i, j, p1, p2, direct});
          double radian = alpha2[static_cast<std::size_t>(2 * k)] -
                          alpha2[static_cast<std::size_t>(2 * k - 1)];
          ds.ncrasegment.push_back({nij.x, nij.y, nij.z, A.x, A.y, A.z, rij, radian});
          ds.Rj.push_back(0.0);  // mesh path: Rj filled only in want_area block

          ds.satom[static_cast<std::size_t>(i)].push_back(nsegment);
          ds.satom[static_cast<std::size_t>(j)].push_back(nsegment);
        }
      }
    }

    // --- loop construction on the i-th SAS-ball -------------------------------
    std::vector<std::vector<int32_t>> loops_i;  // [k] 1-based; [0] dummy
    int nloops_i = 0;
    const int nsatom_i = static_cast<int>(ds.satom[static_cast<std::size_t>(i)].size());

    if (nsatom_i > 0) {
      // satom[i] is 0-based here; loopconstruct wants 1-based sa with [0] dummy.
      std::vector<int32_t> sa_i;
      sa_i.push_back(0);  // [0] dummy
      for (int32_t id : ds.satom[static_cast<std::size_t>(i)]) sa_i.push_back(id);

      LoopResult lr = loopconstruct(i, sa_i, nsatom_i, ds.segment);
      loops_i = std::move(lr.loops);
      nloops_i = lr.nloops;

      // accumulate into global loops; loops_index[i] = [start, end].
      for (int k = 1; k <= nloops_i; ++k) {
        dl.loops.push_back(loops_i[static_cast<std::size_t>(k)]);
      }
      dl.loops_index[static_cast<std::size_t>(i)][0] = nloops + 1;
      dl.loops_index[static_cast<std::size_t>(i)][1] = nloops + nloops_i;
      nloops += nloops_i;
    } else if (ncircleindex(i) > 0) {
      // loops_i = [] -> a single dummy row.
      loops_i.emplace_back();  // loops_i[0] dummy
      nloops_i = 0;
    }

    // --- patch construction ---------------------------------------------------
    if (nsatom_i > 0 || ncircleindex(i) > 0) {
      // circleindex[i] is 0-based; patchesconstruct wants 1-based with [0] dummy.
      std::vector<int32_t> circleindex_i;
      circleindex_i.push_back(0);  // [0] dummy
      for (int32_t cid : circleindex[static_cast<std::size_t>(i)]) circleindex_i.push_back(cid);

      PatchResult pr = patchesconstruct(i, C, R, ds.segment, ds.ncrasegment, I,
                                        circle, circleindex_i, ncircleindex(i),
                                        loops_i, nloops_i, Rp);

      for (int k = 1; k <= pr.npatches; ++k) {
        dp.patches.push_back(pr.patches[static_cast<std::size_t>(k)]);
        dp.patch_atom.push_back(i);
      }
      dp.patches_index[static_cast<std::size_t>(i)][0] = npatches + 1;
      dp.patches_index[static_cast<std::size_t>(i)][1] = npatches + pr.npatches;
      npatches += pr.npatches;

      // TODO(want_area): the per-patch Gauss-Bonnet SAS area (mod_seg_loop_cir +
      // area_spherical -> data_av) is omitted in the mesh path.
    }
  }

  ds.nsegment = nsegment;
  dl.nloops = nloops;
  dp.npatches = npatches;

  return {std::move(ds), std::move(dl), std::move(dp)};
}

}  // namespace meshms
