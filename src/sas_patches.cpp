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
#include <utility>

#include "meshms/parallel.hpp"

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
// Per-loop precompute for interiorloop: every quantity that depends only on the
// atom and the loop ELEMENT (not on `point`) -- the unit vector u1 towards the
// other atom, the other-atom id, the second acos term B, alpha2, and the
// point-independent alpha_local inputs. interiorloop was re-evaluating all of
// these on every call (once per set element per tree node); precomputing them
// once per atom substitutes bit-identical values (same expressions, same
// inputs) and changes nothing else.
struct LoopPre {
  std::vector<Vec3> u1;        // (center[j]-ci)/norm(center[j]-ci), [0] dummy
  std::vector<int32_t> jat;    // other atom of element k
  std::vector<double> B;       // acos_clamped(dot(v2, u1*ri)/pysq(ri))
  std::vector<double> alpha2;  // alpha_local(direct, p3-A, p4-A, n)
  std::vector<int32_t> direct; // segment[loop[k], 5]
  std::vector<Vec3> p3A;       // p3 - A
  std::vector<Vec3> n;         // ncrasegment[loop[k], 1:4]
};

LoopPre build_loop_pre(int i, const Vec3& ci, double ri,
                       const std::vector<int32_t>& loop, int loopsize,
                       const std::vector<Vec3>& center,
                       const std::vector<Vec3>& I,
                       const std::vector<std::array<int32_t, 5>>& segment,
                       const std::vector<std::array<double, 8>>& ncrasegment) {
  LoopPre pk;
  pk.u1.assign(static_cast<std::size_t>(loopsize) + 1, Vec3{});
  pk.jat.assign(static_cast<std::size_t>(loopsize) + 1, 0);
  pk.B.assign(static_cast<std::size_t>(loopsize) + 1, 0.0);
  pk.alpha2.assign(static_cast<std::size_t>(loopsize) + 1, 0.0);
  pk.direct.assign(static_cast<std::size_t>(loopsize) + 1, 0);
  pk.p3A.assign(static_cast<std::size_t>(loopsize) + 1, Vec3{});
  pk.n.assign(static_cast<std::size_t>(loopsize) + 1, Vec3{});
  for (int k = 1; k <= loopsize; ++k) {
    const auto& seg = segment[static_cast<std::size_t>(loop[static_cast<std::size_t>(k)])];
    const auto& ncra = ncrasegment[static_cast<std::size_t>(loop[static_cast<std::size_t>(k)])];
    int j = (seg[0] == i) ? seg[1] : seg[0];
    Vec3 u1 = (center[static_cast<std::size_t>(j)] - ci) /
              norm(center[static_cast<std::size_t>(j)] - ci);
    Vec3 v1 = u1 * ri;  // == (center[j]-ci)/norm(center[j]-ci) * ri exactly
    Vec3 v2 = I[static_cast<std::size_t>(seg[2])] - ci;
    Vec3 n{ncra[0], ncra[1], ncra[2]};
    Vec3 A{ncra[3], ncra[4], ncra[5]};
    Vec3 p3 = I[static_cast<std::size_t>(seg[2])];
    Vec3 p4 = I[static_cast<std::size_t>(seg[3])];
    pk.u1[static_cast<std::size_t>(k)] = u1;
    pk.jat[static_cast<std::size_t>(k)] = j;
    pk.B[static_cast<std::size_t>(k)] = acos_clamped(dot(v2, v1) / pysq(ri));
    pk.alpha2[static_cast<std::size_t>(k)] = alpha_local(seg[4], p3 - A, p4 - A, n);
    pk.direct[static_cast<std::size_t>(k)] = seg[4];
    pk.p3A[static_cast<std::size_t>(k)] = p3 - A;
    pk.n[static_cast<std::size_t>(k)] = n;
  }
  return pk;
}

int interiorloop(const Vec3& point, const Vec3& ci, double ri, int loopsize,
                 const LoopPre& pk, std::vector<int>& K) {
  int nearest = 0;
  double theta0 = 0.0;
  // v3 and pysq(ri) are iteration-invariant (identical value every k).
  const Vec3 v3 = point - ci;
  const double ri2 = pysq(ri);
  for (int k = 1; k <= loopsize; ++k) {
    const Vec3 v1 = pk.u1[static_cast<std::size_t>(k)] * ri;
    // theta = arccos(clip(dot(v3,v1)/ri**2)) - arccos(clip(dot(v2,v1)/ri**2))
    double theta = acos_clamped(dot(v3, v1) / ri2) - pk.B[static_cast<std::size_t>(k)];
    if (k == 1 || theta < theta0) {
      theta0 = theta;
      nearest = k;
    }
  }

  // K[] collects the loop entries whose "other atom" j == j_nearest.
  K.assign(static_cast<std::size_t>(loopsize) + 1, 0);  // [0] dummy
  int nK = 0;
  const int j_nearest = pk.jat[static_cast<std::size_t>(nearest)];
  for (int k = 1; k <= loopsize; ++k) {
    if (pk.jat[static_cast<std::size_t>(k)] == j_nearest) {
      ++nK;
      K[static_cast<std::size_t>(nK)] = k;
    }
  }

  // pk.u1[nearest] IS (center[j_nearest]-ci)/norm(center[j_nearest]-ci).
  const Vec3& v1 = pk.u1[static_cast<std::size_t>(nearest)];
  Vec3 v = v3 - dot(v1, v3) * v1;

  for (int k0 = 1; k0 <= nK; ++k0) {
    int k = K[static_cast<std::size_t>(k0)];
    double alpha1 = alpha_local(pk.direct[static_cast<std::size_t>(k)],
                                pk.p3A[static_cast<std::size_t>(k)], v,
                                pk.n[static_cast<std::size_t>(k)]);
    if (alpha1 < pk.alpha2[static_cast<std::size_t>(k)]) {
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

  const Vec3 ci = C[static_cast<std::size_t>(i)];
  const double ri = R[static_cast<std::size_t>(i)] + Rp;

  // Per-atom precompute: everything the tree walk below re-evaluated per node
  // that depends only on the atom and a loop/circle entry. Identical
  // expressions, evaluated once.
  std::vector<LoopPre> pre(static_cast<std::size_t>(nloops_i) + 1);
  std::vector<Vec3> loop_point(static_cast<std::size_t>(nloops_i) + 1, Vec3{});
  for (int k = 1; k <= nloops_i; ++k) {
    const std::vector<int32_t>& loopk = loops_i[static_cast<std::size_t>(k)];
    const int lsz = static_cast<int>(loopk.size()) - 1;
    pre[static_cast<std::size_t>(k)] =
        build_loop_pre(i, ci, ri, loopk, lsz, C, I, segment, ncrasegment);
    // point = I[segment[loops_i[k, 1], 3]]
    const int seg0 = loopk[1];
    loop_point[static_cast<std::size_t>(k)] =
        I[static_cast<std::size_t>(segment[static_cast<std::size_t>(seg0)][2])];
  }
  std::vector<Vec3> circle_point(static_cast<std::size_t>(ncircleindex_i) + 1, Vec3{});
  for (int t = 1; t <= ncircleindex_i; ++t) {
    const int cidx = circleindex_i[static_cast<std::size_t>(t)];
    const auto& crow = circle[static_cast<std::size_t>(cidx)];
    Vec3 cnorm{crow[6], crow[7], crow[8]};   // circle[cidx, 6:9]
    auto [vector1, v2unused] = orthogonalvectors(cnorm);
    (void)v2unused;
    Vec3 ccenter{crow[3], crow[4], crow[5]};  // circle[cidx, 3:6]
    double cr = crow[9];                      // circle[cidx, 9]
    circle_point[static_cast<std::size_t>(t)] = ccenter + cr * vector1;
  }
  std::vector<int> Kscratch;  // interiorloop K[] buffer, reused across calls

  // S0 = [1, 2, ..., nloops_i, -1, -2, ..., -ncircleindex_i]
  std::vector<int32_t> S0;
  S0.reserve(static_cast<std::size_t>(nloops_i + ncircleindex_i));
  for (int t = 1; t <= nloops_i; ++t) S0.push_back(t);
  for (int t = 1; t <= ncircleindex_i; ++t) S0.push_back(-t);

  // Tree workspace (1-based; index 0 dummy). Reserved to its maximum size (the
  // walk adds at most 2 nodes per iteration of the 2N+1 loop) so growth never
  // reallocates; sets are MOVED into nodes (they are dead at the source after).
  std::vector<TreeNode> tree;
  tree.reserve(static_cast<std::size_t>(2 * (nloops_i + ncircleindex_i) + 4));
  tree.emplace_back();  // dummy slot at index 0
  tree.emplace_back();  // tree[1]
  tree[1].activenode = 1;
  tree[1].activeelement = S0.empty() ? 0 : S0[0];
  tree[1].set = std::move(S0);
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
        const LoopPre& pk = pre[static_cast<std::size_t>(k)];

        // loop entries (1..n1)
        for (int s1 = 1; s1 <= tree[static_cast<std::size_t>(j)].n1; ++s1) {
          int Ss1 = S[static_cast<std::size_t>(s1 - 1)];
          // point = I[segment[loops_i[Ss1, 1], 3]] (precomputed per loop)
          const Vec3& point = loop_point[static_cast<std::size_t>(Ss1)];

          if (k == Ss1 ||
              interiorloop(point, ci, ri, loopsize_k, pk, Kscratch)) {
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
          // point precomputed per circle entry (circleindex_i[-Ss2])
          const Vec3& point = circle_point[static_cast<std::size_t>(-Ss2)];

          if (interiorloop(point, ci, ri, loopsize_k, pk, Kscratch)) {
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
          tree[static_cast<std::size_t>(ntree)].set = std::move(S1);  // S1 dead after
          tree[static_cast<std::size_t>(ntree)].activeelement = k1;
          tree[static_cast<std::size_t>(ntree)].n1 = left_n1;
          tree[static_cast<std::size_t>(ntree)].n2 = left_n2;
          if (t1 == 1) tree[static_cast<std::size_t>(ntree)].activenode = 0;

          ++ntree;
          tree.emplace_back();
          tree[static_cast<std::size_t>(ntree)].activenode = 1;
          tree[static_cast<std::size_t>(ntree)].set = std::move(S2);  // S2 dead after
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

  // ---- PASS 1 (serial): segment creation for every atom -------------------
  // Split from the loop/patch construction below: every segment touching atom a
  // is created at atom min(i,j) <= a, so satom[a] is complete once atom a's own
  // row scan finishes; the loop/patch code for atom a reads only ds rows listed
  // in satom[a], all of which the old interleaved order had already created.
  // The split is therefore byte-identical, and it makes PASS 2 per-atom
  // independent (parallelizable with a serial ascending-i merge).
  // Scratch buffers reused across (i,row): every element is (re)assigned before
  // use, so contents match the old per-iteration fresh vectors exactly.
  std::vector<double> alpha1, alpha2;
  std::vector<int> order, pointord;
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
        alpha1.assign(static_cast<std::size_t>(npt) + 1, 0.0);
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
        order.resize(static_cast<std::size_t>(npt));  // 0-based indices into [1..npt]
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
          return alpha1[static_cast<std::size_t>(a + 1)] < alpha1[static_cast<std::size_t>(b + 1)];
        });

        // alpha2[1..npt] = alpha1[order]; point reordered (icirc[order]).
        alpha2.assign(static_cast<std::size_t>(npt) + 1, 0.0);
        pointord.resize(static_cast<std::size_t>(npt));  // reordered point ids (0-based)
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
  }

  // ---- PASS 2 (parallel per atom): loop + patch construction --------------
  // Reads only PASS-1 outputs (all const from here) plus the const inputs, and
  // writes only its own per_atom slot; the serial ascending-i merge below then
  // replays the exact accumulation (numbering included) of the old interleaved
  // loop. Same LocalMesh-style pattern as the mesher drivers.
  struct AtomLP {
    std::vector<std::vector<int32_t>> loops_i;  // [k] 1-based; [0] dummy
    int nloops_i = 0;
    bool has_loops = false;    // nsatom_i > 0: contributes to dl
    bool has_patches = false;  // nsatom_i > 0 || ncircleindex(i) > 0
    PatchResult pr;
  };
  std::vector<AtomLP> per_atom(static_cast<std::size_t>(M) + 1);

  meshms::parallel_for(1, M + 1, [&](int i) {
    AtomLP& out = per_atom[static_cast<std::size_t>(i)];
    const int nsatom_i = static_cast<int>(ds.satom[static_cast<std::size_t>(i)].size());

    // --- loop construction on the i-th SAS-ball ----------------------------
    if (nsatom_i > 0) {
      // satom[i] is 0-based here; loopconstruct wants 1-based sa with [0] dummy.
      std::vector<int32_t> sa_i;
      sa_i.reserve(static_cast<std::size_t>(nsatom_i) + 1);
      sa_i.push_back(0);  // [0] dummy
      for (int32_t id : ds.satom[static_cast<std::size_t>(i)]) sa_i.push_back(id);

      LoopResult lr = loopconstruct(i, sa_i, nsatom_i, ds.segment);
      out.loops_i = std::move(lr.loops);
      out.nloops_i = lr.nloops;
      out.has_loops = true;
    } else if (ncircleindex(i) > 0) {
      // loops_i = [] -> a single dummy row.
      out.loops_i.emplace_back();  // loops_i[0] dummy
      out.nloops_i = 0;
    }

    // --- patch construction ------------------------------------------------
    if (nsatom_i > 0 || ncircleindex(i) > 0) {
      // circleindex[i] is 0-based; patchesconstruct wants 1-based with [0] dummy.
      std::vector<int32_t> circleindex_i;
      circleindex_i.reserve(circleindex[static_cast<std::size_t>(i)].size() + 1);
      circleindex_i.push_back(0);  // [0] dummy
      for (int32_t cid : circleindex[static_cast<std::size_t>(i)]) circleindex_i.push_back(cid);

      out.pr = patchesconstruct(i, C, R, ds.segment, ds.ncrasegment, I,
                                circle, circleindex_i, ncircleindex(i),
                                out.loops_i, out.nloops_i, Rp);
      out.has_patches = true;

      // TODO(want_area): the per-patch Gauss-Bonnet SAS area (mod_seg_loop_cir +
      // area_spherical -> data_av) is omitted in the mesh path.
    }
  });

  // ---- SERIAL merge (ascending i): exact old accumulation order -----------
  for (int i = 1; i <= M; ++i) {
    AtomLP& a = per_atom[static_cast<std::size_t>(i)];
    if (a.has_loops) {
      // accumulate into global loops; loops_index[i] = [start, end].
      for (int k = 1; k <= a.nloops_i; ++k) {
        dl.loops.push_back(std::move(a.loops_i[static_cast<std::size_t>(k)]));
      }
      dl.loops_index[static_cast<std::size_t>(i)][0] = nloops + 1;
      dl.loops_index[static_cast<std::size_t>(i)][1] = nloops + a.nloops_i;
      nloops += a.nloops_i;
    }
    if (a.has_patches) {
      for (int k = 1; k <= a.pr.npatches; ++k) {
        dp.patches.push_back(std::move(a.pr.patches[static_cast<std::size_t>(k)]));
        dp.patch_atom.push_back(i);
      }
      dp.patches_index[static_cast<std::size_t>(i)][0] = npatches + 1;
      dp.patches_index[static_cast<std::size_t>(i)][1] = npatches + a.pr.npatches;
      npatches += a.pr.npatches;
    }
  }

  ds.nsegment = nsegment;
  dl.nloops = nloops;
  dp.npatches = npatches;

  return {std::move(ds), std::move(dl), std::move(dp)};
}

}  // namespace meshms
