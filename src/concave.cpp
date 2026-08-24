// Concave SES patch construction --- faithful port of the concave module's
// SESconcavepat (MESH PATH ONLY, av=None). See concave.hpp for scope.
#include "meshms/concave.hpp"
#include "meshms/parallel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <numeric>
#include <vector>

#include "meshms/meshing.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

namespace {

constexpr double EPSILON = 1.0e-10;

// ----- record-matrix helpers (1-based-with-dummy rows AND columns) ----------
// A "record matrix" row stores 1-based columns; column 0 is a dummy slot.

// Read a length-3 coordinate from columns [c0, c0+1, c0+2] of a record row.
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

// ----- small numerical helpers (mirror the MATLAB local functions) ----------
inline double det3(const Vec3& r0, const Vec3& r1, const Vec3& r2) {
  return dot(r0, cross(r1, r2));
}

// Local alpha(u,v,n) == arc_angle(u,v,n,1).
inline double alpha_(const Vec3& u, const Vec3& v, const Vec3& n) {
  return arc_angle(u, v, n, 1);
}

// Port of test2: is x1 inside the tetrahedron (x, xi, xj, xk)?
inline bool test2(const Vec3& x1, const Vec3& x, const Vec3& xi, const Vec3& xj,
                  const Vec3& xk) {
  double D0 = sign(det3(xi - x, xj - x, xk - x));
  double D1 = sign(det3(x1 - x, xj - x, xk - x));
  double D2 = sign(det3(xi - x, x1 - x, xk - x));
  double D3 = sign(det3(xi - x, xj - x, x1 - x));
  return (D0 == D1) && (D0 == D2) && (D0 == D3);
}

// Port of test1. `circle` is the 1-based record matrix circle_ij/jk/ki:
// columns [_, c(1:3)=1:3? -> cols 1..3, r=4, j=5]; rows 1..ncircle.
inline bool test1(const Vec3& x1, const Vec3& x, const Vec3& xi, const Vec3& xj,
                  const Vec3& xk, int k,
                  const std::vector<std::array<double, 6>>& circle, int ncircle) {
  for (int i = 1; i <= ncircle; ++i) {
    if (i != k) {
      if (norm(x1 - rget3(circle[static_cast<std::size_t>(i)], 1)) <=
          circle[static_cast<std::size_t>(i)][4] - EPSILON) {
        return false;
      }
    }
  }
  double D0 = sign(det3(xi - x, xj - x, xk - x));
  double D1 = sign(det3(x1 - x, xj - x, xk - x));
  double D2 = sign(det3(xi - x, x1 - x, xk - x));
  if (D0 != D1 || D0 != D2) return false;
  return true;
}

// Port of the local ``intersection`` (three probes of radius r).
// Returns true if the three spheres meet, filling x1/x2.
inline bool intersection_three(const Vec3& c1, const Vec3& c2, const Vec3& c3,
                               double r, Vec3& x1, Vec3& x2) {
  Vec3 A = 0.5 * (c1 + c2);
  Vec3 B = 0.5 * (c2 + c3);
  Vec3 cij = c1 - c2;
  Vec3 ckj = c3 - c2;
  Vec3 cr = cross(cij, ckj);
  Vec3 n = cr / norm(cr);
  Vec3 u = (c3 - c2) - dot(c3 - c2, c1 - c2) * (c1 - c2) / pysq(norm(c1 - c2));
  double t = dot(B - A, c3 - c2) / dot(u, c3 - c2);
  Vec3 X1 = A + t * u;
  double c = -pysq(norm(X1 - c1)) + pysq(r);
  if (c > 0) {
    double sq = std::sqrt(c);
    x1 = X1 - sq * n;
    x2 = X1 + sq * n;
    return true;
  }
  x1 = Vec3{};
  x2 = Vec3{};
  return false;
}

// Port of the local ``intersection_circle`` (two coplanar circles).
inline void intersection_circle(const Vec3& c1, double r1, const Vec3& c2,
                                double r2, const Vec3& n, Vec3& x1, Vec3& x2) {
  double dd = norm(c1 - c2);
  double t = (pysq(r1) - pysq(r2) + pysq(dd)) / (2.0 * dd);
  Vec3 v = (c2 - c1) / norm(c2 - c1);
  Vec3 O = c1 + t * v;
  Vec3 u = cross(n, v);
  u = u / norm(u);
  double h = std::sqrt(std::max(pysq(r1) - pysq(t), 0.0));
  x1 = O + h * u;
  x2 = O - h * u;
}

// ----- TreeNode (port of Classes/treenode.m) --------------------------------
struct TreeNode {
  int activenode = 0;
  std::vector<int> set;  // python list
  int activeelement = 0;
  int n1 = 0;
  int n2 = 0;
  int leftnode = 0;
  int rightnode = 0;
};

// ----- running concave state (av deferred) ----------------------------------
struct ConcaveState {
  double Rp;
  double d;
  int arg_eSAS = 1;
  ConcaveState(double Rp_, double d_) : Rp(Rp_), d(d_) {}
};

// stable_sort of indices 1..N by key, ties keeping input order. `key` is
// 1-based-with-dummy (key[0] unused). Returns a 1-based-with-dummy order array
// (order[0] dummy, order[1..N] the sorted old indices).
inline std::vector<int> stable_order(const std::vector<double>& key, int N) {
  std::vector<int> order(static_cast<std::size_t>(N));
  std::iota(order.begin(), order.end(), 1);  // 1..N
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
    return key[static_cast<std::size_t>(a)] < key[static_cast<std::size_t>(b)];
  });
  std::vector<int> out(static_cast<std::size_t>(N) + 1, 0);
  for (int i = 1; i <= N; ++i) out[static_cast<std::size_t>(i)] = order[static_cast<std::size_t>(i - 1)];
  return out;
}

// ----- sort_boundary --------------------------------------------------------
// boundary: rows 1..N, columns [_, x(1:3)=1:3, j=4]. Reorder by increasing angle
// about axis n, keeping the dummy row 0 intact.
std::vector<std::array<double, 5>> sort_boundary(
    const Vec3& c, const std::vector<std::array<double, 5>>& boundary, int N,
    const Vec3& n) {
  std::vector<double> alpha1(static_cast<std::size_t>(N) + 1, 0.0);
  Vec3 x0 = rget3(boundary[1], 1);
  for (int k = 2; k <= N; ++k) {
    alpha1[static_cast<std::size_t>(k)] =
        alpha_(x0 - c, rget3(boundary[static_cast<std::size_t>(k)], 1) - c, n);
  }
  std::vector<int> order = stable_order(alpha1, N);
  std::vector<std::array<double, 5>> out = boundary;
  for (int newpos = 1; newpos <= N; ++newpos) {
    out[static_cast<std::size_t>(newpos)] = boundary[static_cast<std::size_t>(order[static_cast<std::size_t>(newpos)])];
  }
  return out;
}

// ----- sort_segment ---------------------------------------------------------
// Returns (A, K1, K2, N_new). A/K1/K2 are 1-based-with-dummy, length N_A+1.
// The av-only V_eSES/V_cSES volume bookkeeping is SKIPPED (no-op when av==None);
// the A/K1/K2/N_new geometry it returns is reproduced faithfully.
struct SortSegOut {
  std::vector<double> A;
  std::vector<int> K1;
  std::vector<int> K2;
  int N_new = 0;
};
SortSegOut sort_segment(int j, const std::vector<int>& atom_row, int N,
                        const Vec3& c, const Vec3& n, double r,
                        const std::vector<Vec3>& I,
                        const std::vector<int>& hight_set,
                        const std::vector<int>& K, int Kn,
                        const std::vector<Vec3>& I_probe, double probe,
                        const Vec3& x, const Vec3& xi, const Vec3& xj,
                        const Vec3& xk) {
  // point(k,:) for k=1..N == I_probe(atom(k),:)
  std::vector<Vec3> point(static_cast<std::size_t>(N) + 1);
  for (int k = 1; k <= N; ++k) {
    point[static_cast<std::size_t>(k)] = I_probe[static_cast<std::size_t>(atom_row[static_cast<std::size_t>(k)])];
  }

  std::vector<double> alpha1(static_cast<std::size_t>(N) + 1, 0.0);
  Vec3 x0 = point[1];  // start point
  for (int k = 2; k <= N; ++k) {
    alpha1[static_cast<std::size_t>(k)] = alpha_(x0 - c, point[static_cast<std::size_t>(k)] - c, n);
  }

  std::vector<int> order = stable_order(alpha1, N);  // 1-based-with-dummy
  std::vector<double> alpha2(static_cast<std::size_t>(N) + 1, 0.0);
  for (int newpos = 1; newpos <= N; ++newpos) {
    alpha2[static_cast<std::size_t>(newpos)] = alpha1[static_cast<std::size_t>(order[static_cast<std::size_t>(newpos)])];
  }

  int nhalf = N / 2;
  SortSegOut out;
  out.A.assign(static_cast<std::size_t>(nhalf) + 1, 0.0);
  out.K1.assign(static_cast<std::size_t>(nhalf) + 1, 0);
  out.K2.assign(static_cast<std::size_t>(nhalf) + 1, 0);
  int N_A = 0;

  // order is 0-based in python (order[k-1]); here order[] is 1-based-with-dummy,
  // so python order[k-1] == order[k].
  for (int k = 1; k <= N; ++k) {
    int knext = (k % N) + 1;  // mod(k,N)+1
    const Vec3& pk = point[static_cast<std::size_t>(order[static_cast<std::size_t>(k)])];
    const Vec3& pkn = point[static_cast<std::size_t>(order[static_cast<std::size_t>(knext)])];
    if (norm(pk - pkn) > EPSILON) {
      int true_ = 1;
      Vec3 x_middle = 0.5 * (pk + pkn);
      if (norm(x_middle - c) != 0) {
        x_middle = c + (x_middle - c) / norm(x_middle - c) * r;
      } else {
        Vec3 v1 = pk - c;
        Vec3 u = -1.0 * cross(v1, n);
        x_middle = c + r * u / norm(u);
      }
      if (alpha_(pk - c, x_middle - c, n) > std::numbers::pi) {
        x_middle = 2.0 * c - x_middle;
      }
      for (int i = 1; i <= Kn; ++i) {
        Vec3 center = I[static_cast<std::size_t>(hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(i)])])];
        if (i != j && norm(x_middle - center) < probe) {
          true_ = 0;
          break;
        }
      }
      if (true_ == 1 && test2(x_middle, x, xi, xj, xk)) {
        N_A += 1;
        double theta = alpha2[static_cast<std::size_t>(knext)] - alpha2[static_cast<std::size_t>(k)];
        if (theta < 0) {
          out.A[static_cast<std::size_t>(N_A)] = theta + TWO_PI;
        } else {
          out.A[static_cast<std::size_t>(N_A)] = theta;
        }
        out.K1[static_cast<std::size_t>(N_A)] = atom_row[static_cast<std::size_t>(order[static_cast<std::size_t>(k)])];
        out.K2[static_cast<std::size_t>(N_A)] = atom_row[static_cast<std::size_t>(order[static_cast<std::size_t>(knext)])];
      }
    }
  }

  out.N_new = 2 * N_A;
  out.A.resize(static_cast<std::size_t>(N_A) + 1);
  out.K1.resize(static_cast<std::size_t>(N_A) + 1);
  out.K2.resize(static_cast<std::size_t>(N_A) + 1);
  // av-only area/volume bookkeeping skipped.
  return out;
}

// ----- interiorloop_concave -------------------------------------------------
// loop: 1-based-with-dummy seg indices; segment record matrix
// [_, c(1:3)=1:3, n(4:6)=4:6, r=7, angle=8, k1=9, k2=10].
bool interiorloop_concave(const Vec3& point, const Vec3& x, double probe,
                          const std::vector<int>& loop, int loopsize,
                          const std::vector<Vec3>& I_probe,
                          const std::vector<std::array<double, 11>>& segment,
                          int direct) {
  double theta0 = 0.0;
  int nearest = 0;
  for (int k = 1; k <= loopsize; ++k) {
    const auto& seg = segment[static_cast<std::size_t>(loop[static_cast<std::size_t>(k)])];
    Vec3 v1 = rget3(seg, 4) * probe;
    if (direct == -1) v1 = -1.0 * v1;
    Vec3 v2 = I_probe[static_cast<std::size_t>(static_cast<int>(seg[9]))] - x;
    Vec3 v3 = point - x;
    double theta = acos_clamped(dot(v3, v1) / pysq(probe)) -
                   acos_clamped(dot(v2, v1) / pysq(probe));
    if (k == 1 || theta < theta0) {
      theta0 = theta;
      nearest = k;
    }
  }
  if (theta0 < 0) return false;

  std::vector<int> K(static_cast<std::size_t>(loopsize) + 1, 0);
  int nK = 0;
  const auto& segn = segment[static_cast<std::size_t>(loop[static_cast<std::size_t>(nearest)])];
  for (int k = 1; k <= loopsize; ++k) {
    const auto& seg = segment[static_cast<std::size_t>(loop[static_cast<std::size_t>(k)])];
    if (norm(rget3(seg, 1) - rget3(segn, 1)) == 0 &&
        std::fabs(seg[7] - segn[7]) == 0) {
      nK += 1;
      K[static_cast<std::size_t>(nK)] = k;
    }
  }

  for (int k0 = 1; k0 <= nK; ++k0) {
    int k = K[static_cast<std::size_t>(k0)];
    const auto& seg = segment[static_cast<std::size_t>(loop[static_cast<std::size_t>(k)])];
    Vec3 v1 = rget3(seg, 4) * probe;
    if (direct == -1) v1 = -1.0 * v1;
    Vec3 v3 = point - x;
    Vec3 v = v3 - (dot(v1, v3) * v1) / pysq(probe);
    Vec3 n = rget3(seg, 4);
    Vec3 A = rget3(seg, 1);
    double alpha1 = alpha_(I_probe[static_cast<std::size_t>(static_cast<int>(seg[9]))] - A, v, n);
    double alpha2 = alpha_(I_probe[static_cast<std::size_t>(static_cast<int>(seg[9]))] - A,
                           I_probe[static_cast<std::size_t>(static_cast<int>(seg[10]))] - A, n);
    if (alpha1 < alpha2) return true;
  }
  return false;
}

// ----- loopconstruct_concave ------------------------------------------------
// segment columns: k1 at col 9, k2 at col 10. map_same is 1-based-with-dummy.
// Returns loops (1-based-with-dummy rows & cols of seg indices), loopsize, nloops.
struct LoopsOut {
  std::vector<std::vector<int>> loops;  // [k] is 1-based-with-dummy
  std::vector<int> loopsize;            // 1-based-with-dummy
  int nloops = 0;
};
LoopsOut loopconstruct_concave(const std::vector<std::array<double, 11>>& segment,
                               int N_segment, const std::vector<int>& map_same) {
  // s(k,1)=segment(k,9), s(k,2)=segment(k,10)
  auto s0 = [&](int k) { return static_cast<int>(segment[static_cast<std::size_t>(k)][9]); };
  auto s1 = [&](int k) { return static_cast<int>(segment[static_cast<std::size_t>(k)][10]); };

  std::vector<int> true_(static_cast<std::size_t>(N_segment) + 1, 0);
  int ntrue = 0;

  int nhalf = N_segment / 2;
  LoopsOut out;
  out.loops.assign(static_cast<std::size_t>(nhalf) + 1,
                   std::vector<int>(static_cast<std::size_t>(N_segment) + 1, 0));
  out.loopsize.assign(static_cast<std::size_t>(nhalf) + 1, 0);
  out.nloops = 0;

  for (int k1 = 1; k1 <= N_segment; ++k1) {
    std::vector<int> loop(static_cast<std::size_t>(N_segment) + 1, 0);
    int loopsi = 0;
    if (true_[static_cast<std::size_t>(k1)] == 0 && ntrue < N_segment) {
      out.nloops += 1;
      int pointer = k1;
      ntrue += 1;
      true_[static_cast<std::size_t>(k1)] = 1;
      loopsi += 1;
      loop[static_cast<std::size_t>(loopsi)] = pointer;

      int times = 0;
      while (map_same[static_cast<std::size_t>(s0(k1))] !=
             map_same[static_cast<std::size_t>(s1(pointer))]) {
        for (int k2 = k1 + 1; k2 <= N_segment; ++k2) {
          if (true_[static_cast<std::size_t>(k2)] == 0 &&
              map_same[static_cast<std::size_t>(s0(k2))] ==
                  map_same[static_cast<std::size_t>(s1(pointer))]) {
            pointer = k2;
            loopsi += 1;
            loop[static_cast<std::size_t>(loopsi)] = pointer;
            true_[static_cast<std::size_t>(pointer)] = 1;
            ntrue += 1;
          }
        }
        times += 1;
        if (times > 100) break;
      }

      // Trim the zero-padding so loops[k].size()-1 == loopsize[k]: mesh_sphpat
      // derives each loop's size from loops[k].size()-1 (no loopsize param), so a
      // padded row would mesh phantom dummy (segment0[0]) segments on multi-loop /
      // loop+circle concave patches. (single-loop patches are unaffected.)
      loop.resize(static_cast<std::size_t>(loopsi) + 1);
      out.loops[static_cast<std::size_t>(out.nloops)] = loop;
      out.loopsize[static_cast<std::size_t>(out.nloops)] = loopsi;
    }
  }

  out.loops.resize(static_cast<std::size_t>(out.nloops) + 1);
  out.loopsize.resize(static_cast<std::size_t>(out.nloops) + 1);
  return out;
}

// ----- patchesconstruct_concave ---------------------------------------------
// circle is the 1-based circle_interior record matrix
// [_, c(1:3)=1:3, n(4:6)=4:6, r=7]. Returns patches (1-based-with-dummy rows &
// cols of signed indices), patchesize, npatches.
struct PatchesOut {
  std::vector<std::vector<int>> patches;  // [k] 1-based-with-dummy
  std::vector<int> patchesize;            // 1-based-with-dummy
  int npatches = 0;
};
PatchesOut patchesconstruct_concave(
    int i_probe_center_idx, const std::vector<Vec3>& I,
    const std::vector<Vec3>& I_probe,
    const std::vector<std::array<double, 11>>& segment,
    const std::vector<std::array<double, 8>>& circle, int N_circle,
    const std::vector<std::vector<int>>& loops, const std::vector<int>& loopsize,
    int nloops, int direct, double Rp) {
  Vec3 ci_sphere = I[static_cast<std::size_t>(i_probe_center_idx)];

  int dim = nloops + N_circle;
  PatchesOut out;
  out.patches.assign(static_cast<std::size_t>(dim) + 1,
                     std::vector<int>(static_cast<std::size_t>(dim) + 1, 0));
  out.patchesize.assign(static_cast<std::size_t>(dim) + 1, 0);
  out.npatches = 0;

  // S0 = [1..nloops, -1..-N_circle]
  std::vector<int> S0;
  for (int t = 1; t <= nloops; ++t) S0.push_back(t);
  for (int t = 1; t <= N_circle; ++t) S0.push_back(-t);

  std::vector<TreeNode> tree(1);  // dummy index 0
  {
    TreeNode root;
    root.activenode = 1;
    root.set = S0;
    root.activeelement = S0.empty() ? 0 : S0[0];
    root.n1 = nloops;
    root.n2 = N_circle;
    tree.push_back(root);
  }
  int ntree = 1;
  int j = 1;

  // Python: range(1, 2*(nloops+N_circle)+1 + 1) -> _s = 1..2N+1 (exclusive upper).
  for (int _s = 1; _s <= 2 * (nloops + N_circle) + 1; ++_s) {
    if (j > ntree) break;
    if (j <= ntree && tree[static_cast<std::size_t>(j)].activenode == 1) {
      std::vector<int> S1, S2;
      int k1 = 0, k2 = 0;
      int t1 = 1, t2 = 1;
      int left_n1 = 0, left_n2 = 0, right_n1 = 0, right_n2 = 0;

      std::vector<int> S = tree[static_cast<std::size_t>(j)].set;
      int k = tree[static_cast<std::size_t>(j)].activeelement;

      if (S.size() == 1) {
        tree[static_cast<std::size_t>(j)].activenode = 0;
        j += 1;
        continue;
      }

      if (k > 0) {
        int n1 = tree[static_cast<std::size_t>(j)].n1;
        int n2 = tree[static_cast<std::size_t>(j)].n2;
        // loops part: S(1..n1)
        for (int s1 = 1; s1 <= n1; ++s1) {
          int Ssi = S[static_cast<std::size_t>(s1 - 1)];
          Vec3 point = I_probe[static_cast<std::size_t>(
              static_cast<int>(segment[static_cast<std::size_t>(loops[static_cast<std::size_t>(Ssi)][1])][9]))];
          if (Ssi == k ||
              interiorloop_concave(point, ci_sphere, Rp,
                                   loops[static_cast<std::size_t>(k)],
                                   loopsize[static_cast<std::size_t>(k)], I_probe,
                                   segment, direct)) {
            S1.push_back(Ssi);
            left_n1 += 1;
            if (t1 == 1 && Ssi > k) {
              t1 = 0;
              k1 = Ssi;
            }
          } else {
            S2.push_back(Ssi);
            right_n1 += 1;
            if (t2 == 1 && Ssi > k) {
              t2 = 0;
              k2 = Ssi;
            }
          }
        }
        // circles part: S(n1+1 .. n1+n2)
        for (int s2 = n1 + 1; s2 <= n1 + n2; ++s2) {
          int Ssi = S[static_cast<std::size_t>(s2 - 1)];
          auto [vector1, vector2] = orthogonalvectors(rget3(circle[static_cast<std::size_t>(-Ssi)], 4));
          (void)vector2;
          Vec3 point = rget3(circle[static_cast<std::size_t>(-Ssi)], 1) +
                       circle[static_cast<std::size_t>(-Ssi)][7] * vector1;
          if (interiorloop_concave(point, ci_sphere, Rp,
                                   loops[static_cast<std::size_t>(k)],
                                   loopsize[static_cast<std::size_t>(k)], I_probe,
                                   segment, direct)) {
            S1.push_back(Ssi);
            left_n2 += 1;
            if (t1 == 1) {
              t1 = 0;
              k1 = Ssi;
            }
          } else {
            S2.push_back(Ssi);
            right_n2 += 1;
            if (t2 == 1) {
              t2 = 0;
              k2 = Ssi;
            }
          }
        }

        if (S2.empty()) {
          if (t1 == 1) {
            tree[static_cast<std::size_t>(j)].activenode = 0;
            j += 1;
          } else {
            tree[static_cast<std::size_t>(j)].activeelement = k1;
          }
        } else {
          tree[static_cast<std::size_t>(j)].activenode = 0;
          tree[static_cast<std::size_t>(j)].leftnode = ntree + 1;
          tree[static_cast<std::size_t>(j)].rightnode = ntree + 2;

          ntree += 1;
          TreeNode nodeL;
          nodeL.activenode = 1;
          nodeL.set = S1;
          nodeL.activeelement = k1;
          nodeL.n1 = left_n1;
          nodeL.n2 = left_n2;
          if (t1 == 1) nodeL.activenode = 0;
          tree.push_back(nodeL);

          ntree += 1;
          TreeNode nodeR;
          nodeR.activenode = 1;
          nodeR.set = S2;
          nodeR.activeelement = k2;
          nodeR.n1 = right_n1;
          nodeR.n2 = right_n2;
          if (t2 == 1) nodeR.activenode = 0;
          tree.push_back(nodeR);

          j += 1;
        }
      } else {
        tree[static_cast<std::size_t>(j)].activenode = 0;
        j += 1;
      }
    } else {
      j += 1;
    }
  }

  for (int jj = 1; jj <= ntree; ++jj) {
    if (tree[static_cast<std::size_t>(jj)].activenode == 0 &&
        tree[static_cast<std::size_t>(jj)].leftnode == 0) {
      out.npatches += 1;
      out.patchesize[static_cast<std::size_t>(out.npatches)] =
          tree[static_cast<std::size_t>(jj)].n1 + tree[static_cast<std::size_t>(jj)].n2;
      const std::vector<int>& S = tree[static_cast<std::size_t>(jj)].set;
      for (int col = 1; col <= static_cast<int>(S.size()); ++col) {
        out.patches[static_cast<std::size_t>(out.npatches)][static_cast<std::size_t>(col)] =
            S[static_cast<std::size_t>(col - 1)];
      }
    }
  }

  out.patches.resize(static_cast<std::size_t>(out.npatches) + 1);
  out.patchesize.resize(static_cast<std::size_t>(out.npatches) + 1);
  return out;
}

// ----- fill_vatom_nearest3 --------------------------------------------------
// Attribute every vertex of a meshed concave patch to the nearer of the patch's
// three defining atoms (the "neighbour side" the spec allows for concave/torus).
// Pure post-pass over lm.P (1-based, P[0] dummy) -- topology untouched,
// deterministic and thread-safe. No-op when the patch is not emitted.
void fill_vatom_nearest3(LocalMesh& lm, int32_t ai, int32_t aj, int32_t ak,
                         const Vec3& ci, const Vec3& cj, const Vec3& ck) {
  if (!lm.emit) return;
  const int32_t atoms[3] = {ai, aj, ak};
  const Vec3 cs[3] = {ci, cj, ck};
  const std::size_t Np = lm.P.empty() ? 0 : lm.P.size() - 1;
  lm.vatom.assign(Np, 0);
  for (std::size_t k = 1; k <= Np; ++k)
    lm.vatom[k - 1] = nearest_atom(lm.P[k], atoms, cs, 3);
}

// ----- construct_concavepat -------------------------------------------------
// Emits its meshed patches (in order) into `out` instead of state.add_patch, so
// the per-probe loop can run in parallel into per-iteration LocalMesh lists.
void decomp_construct_concavepat(ProbePatchSet& pd, double Rp, int k1,
                          const std::vector<Vec3>& I,
                          const std::vector<Vec3>& I_probe, int N_probe,
                          const std::vector<std::array<double, 11>>& segment,
                          int N_segment,
                          const std::vector<std::array<double, 8>>& circle_interior,
                          int N_circle, int direct) {
  // --- group identical I_probe points (map_same) ---
  std::vector<int> map_same(static_cast<std::size_t>(N_probe) + 1, 0);
  std::vector<int> samepoint_rep(static_cast<std::size_t>(N_probe) + 1, 0);
  int nsamepoint = 0;
  for (int i = 1; i <= N_probe; ++i) {
    if (i == 1) {
      nsamepoint = 1;
      samepoint_rep[1] = i;
      map_same[static_cast<std::size_t>(i)] = nsamepoint;
    } else {
      int flag = 1;
      for (int jj = 1; jj <= nsamepoint; ++jj) {
        if (norm(I_probe[static_cast<std::size_t>(samepoint_rep[static_cast<std::size_t>(jj)])] -
                 I_probe[static_cast<std::size_t>(i)]) < EPSILON) {
          map_same[static_cast<std::size_t>(i)] = jj;
          flag = 0;
          break;
        }
      }
      if (flag == 1) {
        nsamepoint += 1;
        samepoint_rep[static_cast<std::size_t>(nsamepoint)] = i;
        map_same[static_cast<std::size_t>(i)] = nsamepoint;
      }
    }
  }

  // --- build loops ---
  LoopsOut lo = loopconstruct_concave(segment, N_segment, map_same);

  // --- build spherical patches ---
  PatchesOut po = patchesconstruct_concave(k1, I, I_probe, segment, circle_interior,
                                           N_circle, lo.loops, lo.loopsize,
                                           lo.nloops, direct, Rp);

  // --- build segment0 ([c,n,r,spoint,angle]) ---
  std::vector<std::array<double, 12>> segment0(static_cast<std::size_t>(N_segment) + 1,
                                               std::array<double, 12>{});
  for (int i = 1; i <= N_segment; ++i) {
    const auto& seg = segment[static_cast<std::size_t>(i)];
    auto& s0row = segment0[static_cast<std::size_t>(i)];
    if (direct == 1) {
      // segment0[i,1:8] = segment[i,1:8] (cols 1..7)
      for (int c = 1; c <= 7; ++c) s0row[static_cast<std::size_t>(c)] = seg[static_cast<std::size_t>(c)];
      rset3(s0row, 8, I_probe[static_cast<std::size_t>(static_cast<int>(seg[9]))]);
      s0row[11] = seg[8];
    } else {
      rset3(s0row, 1, rget3(seg, 1));
      rset3(s0row, 4, -1.0 * rget3(seg, 4));
      s0row[7] = seg[7];
      rset3(s0row, 8, I_probe[static_cast<std::size_t>(static_cast<int>(seg[10]))]);
      s0row[11] = seg[8];
    }
  }

  // --- build circle0 ([c,n,r]) ---
  // circle0 record matrix (*,9): [_, c(1:3)=1:3, n(4:6)=4:6, r=7, torus=8].
  std::vector<std::array<double, 9>> circle0(static_cast<std::size_t>(N_circle) + 1,
                                             std::array<double, 9>{});
  for (int i = 1; i <= N_circle; ++i) {
    const auto& ci = circle_interior[static_cast<std::size_t>(i)];
    auto& c0 = circle0[static_cast<std::size_t>(i)];
    for (int c = 1; c <= 7; ++c) c0[static_cast<std::size_t>(c)] = ci[static_cast<std::size_t>(c)];
    c0[8] = 0.0;  // no torus radius
  }
  if (direct == -1) {
    for (int i = 1; i <= N_circle; ++i) {
      rset3(circle0[static_cast<std::size_t>(i)], 4,
            -1.0 * rget3(circle_interior[static_cast<std::size_t>(i)], 4));
    }
    for (int i = 1; i <= lo.nloops; ++i) {
      int ls = lo.loopsize[static_cast<std::size_t>(i)];
      // flipud of loops(i,1:loopsize(i))
      std::vector<int>& loop = lo.loops[static_cast<std::size_t>(i)];
      std::reverse(loop.begin() + 1, loop.begin() + 1 + ls);
    }
  }

  // --- store the decomposition; meshing happens in mesh_probe_patchset -----
  pd.k1 = k1;
  pd.loops = std::move(lo.loops);
  pd.segment0 = std::move(segment0);
  pd.circle0 = std::move(circle0);
  pd.patches = std::move(po.patches);
  pd.patchesize = std::move(po.patchesize);
  pd.npatches = po.npatches;
}

// Mesh every patch of one probe's decomposition (the cSAS full-trim pass; the
// eSAS pass is area-only and never meshes). Bit-identical inputs and call order
// to the old construct_concavepat mesh loop + the fill_vatom_nearest3 pass.
void mesh_probe_patchset(std::vector<LocalMesh>& out, const ProbePatchSet& pd,
                         const std::vector<Vec3>& I, double Rp, double d) {
  const std::size_t vstart = out.size();
  for (int i = 1; i <= pd.npatches; ++i) {
    Tag btag{2, pd.k1, 0};  // ("probe", k1)
    out.push_back(mesh_sphpat(
        I[static_cast<std::size_t>(pd.k1)], Rp, pd.loops, pd.segment0, pd.circle0,
        pd.patches[static_cast<std::size_t>(i)],
        pd.patchesize[static_cast<std::size_t>(i)], Rp, d, nullptr, btag));
  }
  // Attribute the just-meshed concave patches to the nearest of {a_i,a_j,a_k}.
  for (std::size_t t = vstart; t < out.size(); ++t) {
    fill_vatom_nearest3(out[t], pd.a_i, pd.a_j, pd.a_k, pd.ci, pd.cj, pd.ck);
  }
}

// ----- data_concavepat ------------------------------------------------------
// Computes probe i's DENSITY-INDEPENDENT decomposition into `pd`; all scratch is
// call-local (thread-safe). Meshing happens separately (mesh_probe_patchset).
void data_concavepat(ProbePatchSet& pd, int i,
                     const std::vector<int>& hight_set,
                     const std::vector<int>& K_in, int Kn_in,
                     const std::vector<Vec3>& I,
                     const std::vector<std::array<int32_t, 3>>& Iijk,
                     const std::vector<Vec3>& C, double probe,
                     const std::vector<std::array<int32_t, 3>>& direction) {
  int i0 = hight_set[static_cast<std::size_t>(i)];
  // Iijk[i0,1:4] -> Iijk[i0][0..2]; direction[i0,1] -> direction[i0][0].
  int a_i = Iijk[static_cast<std::size_t>(i0)][0];
  int a_j = Iijk[static_cast<std::size_t>(i0)][1];
  int a_k = Iijk[static_cast<std::size_t>(i0)][2];
  int dir_i0 = direction[static_cast<std::size_t>(i0)][0];

  Vec3 x = I[static_cast<std::size_t>(i0)];
  Vec3 ci = C[static_cast<std::size_t>(a_i)];
  Vec3 cj = C[static_cast<std::size_t>(a_j)];
  Vec3 ck = C[static_cast<std::size_t>(a_k)];
  Vec3 ni = (ci - x) / norm(ci - x);
  Vec3 nj = (cj - x) / norm(cj - x);
  Vec3 nk = (ck - x) / norm(ck - x);
  Vec3 nij = cross(ni, nj);
  Vec3 njk = cross(nj, nk);
  Vec3 nki = cross(nk, ni);
  nij = nij / norm(nij);
  njk = njk / norm(njk);
  nki = nki / norm(nki);
  Vec3 xi = x + probe * ni;
  Vec3 xj = x + probe * nj;
  Vec3 xk = x + probe * nk;

  // Iijk membership helper: np.any(a - Iijk[j0,1:4] == 0).
  auto any_eq = [&](int a, int j0) {
    const auto& row = Iijk[static_cast<std::size_t>(j0)];
    return (a == row[0]) || (a == row[1]) || (a == row[2]);
  };

  // local copy of K (1-based-with-dummy) and Kn
  int Kn = Kn_in;
  std::vector<int> K(static_cast<std::size_t>(Kn) + 1, 0);
  for (int jj = 1; jj <= Kn; ++jj) K[static_cast<std::size_t>(jj)] = K_in[static_cast<std::size_t>(jj)];

  // --- modify K, Kn by removing some unuseful nearby probe ---
  std::vector<int> remove(static_cast<std::size_t>(Kn) + 1, 0);
  int Num_ij = 0, Num_jk = 0, Num_ki = 0;
  int Kij = 0, Kjk = 0, Kki = 0;
  double alphaij = 0.0, alphajk = 0.0, alphaki = 0.0;
  for (int j = 1; j <= Kn; ++j) {
    int j0 = hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(j)])];
    Vec3 x0 = I[static_cast<std::size_t>(j0)];
    // ij edge
    if (any_eq(a_i, j0) && any_eq(a_j, j0)) {
      remove[static_cast<std::size_t>(j)] = 1;
      Num_ij += 1;
      Vec3 cj_ci = cj - ci;
      Vec3 u = -1.0 * ni - (dot(-1.0 * ni, cj_ci) / norm(cj_ci)) * cj_ci / norm(cj_ci);
      Vec3 v = (x0 - ci) - (dot(x0 - ci, cj_ci) / norm(cj_ci)) * cj_ci / norm(cj_ci);
      double alpha0 = alpha_(u, v, cj_ci);
      if (dir_i0 == -1) alpha0 = TWO_PI - alpha0;
      if (Num_ij == 1) {
        Kij = j;
        alphaij = alpha0;
      } else if (alpha0 > alphaij) {
        Kij = j;
        alphaij = alpha0;
      }
    }
    // jk edge
    if (any_eq(a_j, j0) && any_eq(a_k, j0)) {
      remove[static_cast<std::size_t>(j)] = 1;
      Num_jk += 1;
      Vec3 ck_cj = ck - cj;
      Vec3 u = -1.0 * nj - (dot(-1.0 * nj, ck_cj) / norm(ck_cj)) * ck_cj / norm(ck_cj);
      Vec3 v = (x0 - cj) - (dot(x0 - cj, ck_cj) / norm(ck_cj)) * ck_cj / norm(ck_cj);
      double alpha0 = alpha_(u, v, ck_cj);
      if (dir_i0 == -1) alpha0 = TWO_PI - alpha0;
      if (Num_jk == 1) {
        Kjk = j;
        alphajk = alpha0;
      } else if (alpha0 > alphajk) {
        Kjk = j;
        alphajk = alpha0;
      }
    }
    // ki edge
    if (any_eq(a_k, j0) && any_eq(a_i, j0)) {
      remove[static_cast<std::size_t>(j)] = 1;
      Num_ki += 1;
      Vec3 ci_ck = ci - ck;
      Vec3 u = -1.0 * nk - (dot(-1.0 * nk, ci_ck) / norm(ci_ck)) * ci_ck / norm(ci_ck);
      Vec3 v = (x0 - ck) - (dot(x0 - ck, ci_ck) / norm(ci_ck)) * ci_ck / norm(ci_ck);
      double alpha0 = alpha_(u, v, ci_ck);
      if (dir_i0 == -1) alpha0 = TWO_PI - alpha0;
      if (Num_ki == 1) {
        Kki = j;
        alphaki = alpha0;
      } else if (alpha0 > alphaki) {
        Kki = j;
        alphaki = alpha0;
      }
    }
  }
  if (Kij > 0) remove[static_cast<std::size_t>(Kij)] = 0;
  if (Kjk > 0) remove[static_cast<std::size_t>(Kjk)] = 0;
  if (Kki > 0) remove[static_cast<std::size_t>(Kki)] = 0;

  int Kn0 = 0;
  std::vector<int> K0(static_cast<std::size_t>(Kn) + 1, 0);
  for (int j = 1; j <= Kn; ++j) {
    if (remove[static_cast<std::size_t>(j)] == 0) {
      Kn0 += 1;
      K0[static_cast<std::size_t>(Kn0)] = K[static_cast<std::size_t>(j)];
    }
  }
  K = K0;
  K.resize(static_cast<std::size_t>(Kn0) + 1);
  Kn = Kn0;

  // --- compute circles on the three planes (ij, jk, ki) ---
  std::vector<std::array<double, 6>> circle_ij(static_cast<std::size_t>(Kn) + 1,
                                               std::array<double, 6>{});
  std::vector<std::array<double, 6>> circle_jk(static_cast<std::size_t>(Kn) + 1,
                                               std::array<double, 6>{});
  std::vector<std::array<double, 6>> circle_ki(static_cast<std::size_t>(Kn) + 1,
                                               std::array<double, 6>{});
  int ncircle_ij = 0, ncircle_jk = 0, ncircle_ki = 0;

  for (int j = 1; j <= Kn; ++j) {
    int j0 = hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(j)])];
    Vec3 x0 = I[static_cast<std::size_t>(j0)];
    double tij = dot(x - x0, nij);
    if (std::fabs(tij) < probe) {
      Vec3 cij = x0 + tij * nij;
      double rij = std::sqrt(std::max(pysq(probe) - pysq(tij), 0.0));
      ncircle_ij += 1;
      rset3(circle_ij[static_cast<std::size_t>(ncircle_ij)], 1, cij);
      circle_ij[static_cast<std::size_t>(ncircle_ij)][4] = rij;
      circle_ij[static_cast<std::size_t>(ncircle_ij)][5] = j;
    }
    double tjk = dot(x - x0, njk);
    if (std::fabs(tjk) < probe) {
      Vec3 cjk = x0 + tjk * njk;
      double rjk = std::sqrt(std::max(pysq(probe) - pysq(tjk), 0.0));
      ncircle_jk += 1;
      rset3(circle_jk[static_cast<std::size_t>(ncircle_jk)], 1, cjk);
      circle_jk[static_cast<std::size_t>(ncircle_jk)][4] = rjk;
      circle_jk[static_cast<std::size_t>(ncircle_jk)][5] = j;
    }
    double tki = dot(x - x0, nki);
    if (std::fabs(tki) < probe) {
      Vec3 cki = x0 + tki * nki;
      double rki = std::sqrt(std::max(pysq(probe) - pysq(tki), 0.0));
      ncircle_ki += 1;
      rset3(circle_ki[static_cast<std::size_t>(ncircle_ki)], 1, cki);
      circle_ki[static_cast<std::size_t>(ncircle_ki)][4] = rki;
      circle_ki[static_cast<std::size_t>(ncircle_ki)][5] = j;
    }
  }

  // --- intersection points on the boundaries ---
  int bcap = 2 * Kn + 3;
  std::vector<std::array<double, 5>> boundary_ij(static_cast<std::size_t>(bcap) + 1,
                                                 std::array<double, 5>{});
  std::vector<std::array<double, 5>> boundary_jk(static_cast<std::size_t>(bcap) + 1,
                                                 std::array<double, 5>{});
  std::vector<std::array<double, 5>> boundary_ki(static_cast<std::size_t>(bcap) + 1,
                                                 std::array<double, 5>{});
  int N_ij = 2, N_jk = 2, N_ki = 2;
  rset3(boundary_ij[1], 1, xi);
  boundary_ij[1][4] = 0;
  rset3(boundary_ij[2], 1, xj);
  boundary_ij[2][4] = 0;
  rset3(boundary_jk[1], 1, xj);
  boundary_jk[1][4] = 0;
  rset3(boundary_jk[2], 1, xk);
  boundary_jk[2][4] = 0;
  rset3(boundary_ki[1], 1, xk);
  boundary_ki[1][4] = 0;
  rset3(boundary_ki[2], 1, xi);
  boundary_ki[2][4] = 0;

  auto min_abs_dist = [&](const std::vector<std::array<double, 5>>& boundary, int N,
                          const Vec3& xpt) {
    double best = std::numeric_limits<double>::infinity();
    for (int r = 1; r <= N; ++r) {
      Vec3 d3 = rget3(boundary[static_cast<std::size_t>(r)], 1) - xpt;
      double s = std::fabs(d3.x) + std::fabs(d3.y) + std::fabs(d3.z);
      best = std::min(best, s);
    }
    return best;
  };

  for (int k = 1; k <= ncircle_ij; ++k) {
    Vec3 cij = rget3(circle_ij[static_cast<std::size_t>(k)], 1);
    double rij = circle_ij[static_cast<std::size_t>(k)][4];
    int jcol = static_cast<int>(circle_ij[static_cast<std::size_t>(k)][5]);
    double nrm = norm(cij - x);
    if (probe - rij < nrm && nrm < probe + rij) {
      Vec3 x1, x2;
      intersection_circle(x, probe, cij, rij, nij, x1, x2);
      if (test1(x1, x, xi, xj, xk, k, circle_ij, ncircle_ij) &&
          min_abs_dist(boundary_ij, N_ij, x1) > EPSILON) {
        N_ij += 1;
        rset3(boundary_ij[static_cast<std::size_t>(N_ij)], 1, x1);
        boundary_ij[static_cast<std::size_t>(N_ij)][4] = jcol;
      }
      if (test1(x2, x, xi, xj, xk, k, circle_ij, ncircle_ij) &&
          min_abs_dist(boundary_ij, N_ij, x2) > EPSILON) {
        N_ij += 1;
        rset3(boundary_ij[static_cast<std::size_t>(N_ij)], 1, x2);
        boundary_ij[static_cast<std::size_t>(N_ij)][4] = jcol;
      }
    }
  }
  for (int k = 1; k <= ncircle_jk; ++k) {
    Vec3 cjk = rget3(circle_jk[static_cast<std::size_t>(k)], 1);
    double rjk = circle_jk[static_cast<std::size_t>(k)][4];
    int jcol = static_cast<int>(circle_jk[static_cast<std::size_t>(k)][5]);
    double nrm = norm(cjk - x);
    if (probe - rjk < nrm && nrm < probe + rjk) {
      Vec3 x1, x2;
      intersection_circle(x, probe, cjk, rjk, njk, x1, x2);
      if (test1(x1, x, xj, xk, xi, k, circle_jk, ncircle_jk) &&
          min_abs_dist(boundary_jk, N_jk, x1) > EPSILON) {
        N_jk += 1;
        rset3(boundary_jk[static_cast<std::size_t>(N_jk)], 1, x1);
        boundary_jk[static_cast<std::size_t>(N_jk)][4] = jcol;
      }
      if (test1(x2, x, xj, xk, xi, k, circle_jk, ncircle_jk) &&
          min_abs_dist(boundary_jk, N_jk, x2) > EPSILON) {
        N_jk += 1;
        rset3(boundary_jk[static_cast<std::size_t>(N_jk)], 1, x2);
        boundary_jk[static_cast<std::size_t>(N_jk)][4] = jcol;
      }
    }
  }
  for (int k = 1; k <= ncircle_ki; ++k) {
    Vec3 cki = rget3(circle_ki[static_cast<std::size_t>(k)], 1);
    double rki = circle_ki[static_cast<std::size_t>(k)][4];
    int jcol = static_cast<int>(circle_ki[static_cast<std::size_t>(k)][5]);
    double nrm = norm(cki - x);
    if (probe - rki < nrm && nrm < probe + rki) {
      Vec3 x1, x2;
      intersection_circle(x, probe, cki, rki, nki, x1, x2);
      if (test1(x1, x, xk, xi, xj, k, circle_ki, ncircle_ki) &&
          min_abs_dist(boundary_ki, N_ki, x1) > EPSILON) {
        N_ki += 1;
        rset3(boundary_ki[static_cast<std::size_t>(N_ki)], 1, x1);
        boundary_ki[static_cast<std::size_t>(N_ki)][4] = jcol;
      }
      if (test1(x2, x, xk, xi, xj, k, circle_ki, ncircle_ki) &&
          min_abs_dist(boundary_ki, N_ki, x2) > EPSILON) {
        N_ki += 1;
        rset3(boundary_ki[static_cast<std::size_t>(N_ki)], 1, x2);
        boundary_ki[static_cast<std::size_t>(N_ki)][4] = jcol;
      }
    }
  }

  // --- I_interior: interior intersection points (three probes) ---
  int icap = Kn * (Kn - 1) + 1;
  if (icap < 1) icap = 1;
  std::vector<std::array<double, 6>> I_interior(static_cast<std::size_t>(icap) + 1,
                                                std::array<double, 6>{});
  int N_interior = 0;
  for (int j = 1; j <= Kn; ++j) {
    int j0 = hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(j)])];
    Vec3 x0 = I[static_cast<std::size_t>(j0)];
    for (int k = j + 1; k <= Kn; ++k) {
      int k0 = hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(k)])];
      Vec3 y0 = I[static_cast<std::size_t>(k0)];
      Vec3 x1, x2;
      bool meet = intersection_three(x, x0, y0, probe, x1, x2);
      if (meet) {
        if (test2(x1, x, xi, xj, xk)) {
          int true_test2 = 1;
          for (int t = 1; t <= Kn; ++t) {
            if (t != j && t != k) {
              int t0 = hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(t)])];
              Vec3 z0 = I[static_cast<std::size_t>(t0)];
              if (norm(x1 - z0) < probe) {
                true_test2 = 0;
                break;
              }
            }
          }
          if (true_test2 == 1) {
            N_interior += 1;
            rset3(I_interior[static_cast<std::size_t>(N_interior)], 1, x1);
            I_interior[static_cast<std::size_t>(N_interior)][4] = j;
            I_interior[static_cast<std::size_t>(N_interior)][5] = k;
          }
        }
        if (test2(x2, x, xi, xj, xk)) {
          int true_test2 = 1;
          for (int t = 1; t <= Kn; ++t) {
            if (t != j && t != k) {
              int t0 = hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(t)])];
              Vec3 z0 = I[static_cast<std::size_t>(t0)];
              if (norm(x2 - z0) < probe) {
                true_test2 = 0;
                break;
              }
            }
          }
          if (true_test2 == 1) {
            N_interior += 1;
            rset3(I_interior[static_cast<std::size_t>(N_interior)], 1, x2);
            I_interior[static_cast<std::size_t>(N_interior)][4] = j;
            I_interior[static_cast<std::size_t>(N_interior)][5] = k;
          }
        }
      }
    }
  }

  // --- segments on the boundaries: sort each boundary by angle ---
  boundary_ij = sort_boundary(x, boundary_ij, N_ij, nij);
  boundary_jk = sort_boundary(x, boundary_jk, N_jk, njk);
  boundary_ki = sort_boundary(x, boundary_ki, N_ki, nki);

  // --- assemble I_probe (coordinate array, dummy row 0) ---
  int N_boundary = N_ij + N_jk + N_ki;
  int N_probe = N_ij + N_jk + N_ki + N_interior;
  std::vector<Vec3> I_probe(static_cast<std::size_t>(N_probe) + 1);
  for (int r = 1; r <= N_ij; ++r) {
    I_probe[static_cast<std::size_t>(r)] = rget3(boundary_ij[static_cast<std::size_t>(r)], 1);
  }
  for (int r = 1; r <= N_jk; ++r) {
    I_probe[static_cast<std::size_t>(N_ij + r)] = rget3(boundary_jk[static_cast<std::size_t>(r)], 1);
  }
  for (int r = 1; r <= N_ki; ++r) {
    I_probe[static_cast<std::size_t>(N_ij + N_jk + r)] = rget3(boundary_ki[static_cast<std::size_t>(r)], 1);
  }
  for (int r = 1; r <= N_interior; ++r) {
    I_probe[static_cast<std::size_t>(N_boundary + r)] = rget3(I_interior[static_cast<std::size_t>(r)], 1);
  }

  // --- atom: for each probe, which I_probe points lie on it ---
  int acap = std::max(11, N_probe + 2);
  std::vector<std::vector<int>> atom(static_cast<std::size_t>(Kn) + 1,
                                     std::vector<int>(static_cast<std::size_t>(acap), 0));
  std::vector<int> Natom(static_cast<std::size_t>(Kn) + 1, 0);
  for (int j = 1; j <= Kn; ++j) {
    for (int k = 1; k <= N_probe; ++k) {
      if (k <= N_ij) {
        int k0 = k;
        if (static_cast<int>(boundary_ij[static_cast<std::size_t>(k0)][4]) == j) {
          Natom[static_cast<std::size_t>(j)] += 1;
          atom[static_cast<std::size_t>(j)][static_cast<std::size_t>(Natom[static_cast<std::size_t>(j)])] = k;
        }
      } else if (k <= N_ij + N_jk) {
        int k0 = k - N_ij;
        if (static_cast<int>(boundary_jk[static_cast<std::size_t>(k0)][4]) == j) {
          Natom[static_cast<std::size_t>(j)] += 1;
          atom[static_cast<std::size_t>(j)][static_cast<std::size_t>(Natom[static_cast<std::size_t>(j)])] = k;
        }
      } else if (k <= N_boundary) {
        int k0 = k - N_ij - N_jk;
        if (static_cast<int>(boundary_ki[static_cast<std::size_t>(k0)][4]) == j) {
          Natom[static_cast<std::size_t>(j)] += 1;
          atom[static_cast<std::size_t>(j)][static_cast<std::size_t>(Natom[static_cast<std::size_t>(j)])] = k;
        }
      } else {
        int k0 = k - N_boundary;
        if (static_cast<int>(I_interior[static_cast<std::size_t>(k0)][4]) == j) {
          Natom[static_cast<std::size_t>(j)] += 1;
          atom[static_cast<std::size_t>(j)][static_cast<std::size_t>(Natom[static_cast<std::size_t>(j)])] = k;
        }
        if (static_cast<int>(I_interior[static_cast<std::size_t>(k0)][5]) == j) {
          Natom[static_cast<std::size_t>(j)] += 1;
          atom[static_cast<std::size_t>(j)][static_cast<std::size_t>(Natom[static_cast<std::size_t>(j)])] = k;
        }
      }
    }
  }

  // --- interior circles ---
  std::vector<std::array<double, 8>> circle_interior(static_cast<std::size_t>(Kn) + 1,
                                                     std::array<double, 8>{});
  int N_circle = 0;
  for (int j = 1; j <= Kn; ++j) {
    if (Natom[static_cast<std::size_t>(j)] == 0) {
      Vec3 x1 = 0.5 * (I[static_cast<std::size_t>(hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(j)])])] + x);
      if (test2(x1, x, xi, xj, xk)) {
        double r1 = std::sqrt(std::max(pysq(probe) - pysq(norm(x - x1)), 0.0));
        Vec3 n1 = static_cast<double>(dir_i0) * (x1 - x) / norm(x1 - x);
        int true1 = 1;
        for (int k = 1; k <= Kn; ++k) {
          if (j != k && norm(x1 - I[static_cast<std::size_t>(hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(k)])])]) < probe) {
            true1 = 0;
            break;
          }
        }
        if (true1 == 1) {
          N_circle += 1;
          rset3(circle_interior[static_cast<std::size_t>(N_circle)], 1, x1);
          rset3(circle_interior[static_cast<std::size_t>(N_circle)], 4, n1);
          circle_interior[static_cast<std::size_t>(N_circle)][7] = r1;
          // av-only V_eSES/V_cSES bookkeeping skipped.
        }
      }
    }
  }
  circle_interior.resize(static_cast<std::size_t>(N_circle) + 1);

  // --- segments: boundary then interior ---
  int segcap = 2 * N_probe + 2;
  std::vector<std::array<double, 11>> segment(static_cast<std::size_t>(segcap) + 1,
                                              std::array<double, 11>{});
  // boundary segments: j = 1,3,5,...,N_boundary-1
  for (int j = 1; j <= N_boundary - 1; j += 2) {
    Vec3 uj = I_probe[static_cast<std::size_t>(j)] - x;
    Vec3 vj = I_probe[static_cast<std::size_t>(j + 1)] - x;
    double aj = acos_clamped(dot(uj, vj) / pysq(probe));
    Vec3 nvec;
    int j1;
    if (j <= N_ij) {
      nvec = nij;
      if (j < N_ij - 1) j1 = j + 1; else j1 = j + 2;
    } else if (j <= N_ij + N_jk) {
      nvec = njk;
      if (j < N_ij + N_jk - 1) j1 = j + 1; else j1 = j + 2;
    } else {
      nvec = nki;
      if (j < N_boundary - 1) j1 = j + 1; else j1 = 1;
    }
    int row = (j + 1) / 2;
    rset3(segment[static_cast<std::size_t>(row)], 1, x);
    rset3(segment[static_cast<std::size_t>(row)], 4, nvec);
    segment[static_cast<std::size_t>(row)][7] = probe;
    segment[static_cast<std::size_t>(row)][8] = aj;
    segment[static_cast<std::size_t>(row)][9] = j;
    segment[static_cast<std::size_t>(row)][10] = j1;
  }

  int N_segment = N_boundary / 2;

  // interior segments
  for (int j = 1; j <= Kn; ++j) {
    if (Natom[static_cast<std::size_t>(j)] > 0) {
      Vec3 x1 = 0.5 * (I[static_cast<std::size_t>(hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(j)])])] + x);
      double r1 = std::sqrt(std::max(pysq(probe) - pysq(norm(x - x1)), 0.0));
      Vec3 hk = I[static_cast<std::size_t>(hight_set[static_cast<std::size_t>(K[static_cast<std::size_t>(j)])])];
      Vec3 n1 = static_cast<double>(dir_i0) * (hk - x) / norm(hk - x);

      std::vector<int> atom_row(static_cast<std::size_t>(Natom[static_cast<std::size_t>(j)]) + 1, 0);
      for (int t = 1; t <= Natom[static_cast<std::size_t>(j)]; ++t) {
        atom_row[static_cast<std::size_t>(t)] = atom[static_cast<std::size_t>(j)][static_cast<std::size_t>(t)];
      }
      SortSegOut ss = sort_segment(j, atom_row, Natom[static_cast<std::size_t>(j)], x1, n1,
                                   r1, I, hight_set, K, Kn, I_probe, probe, x, xi, xj, xk);

      N_segment = N_segment + ss.N_new / 2;
      int half = ss.N_new / 2;
      for (int t = 1; t <= half; ++t) {
        int row = N_segment - half + t;
        rset3(segment[static_cast<std::size_t>(row)], 1, x1);
        rset3(segment[static_cast<std::size_t>(row)], 4, n1);
        segment[static_cast<std::size_t>(row)][7] = r1;
        segment[static_cast<std::size_t>(row)][8] = ss.A[static_cast<std::size_t>(t)];
        segment[static_cast<std::size_t>(row)][9] = ss.K1[static_cast<std::size_t>(t)];
        segment[static_cast<std::size_t>(row)][10] = ss.K2[static_cast<std::size_t>(t)];
      }
    }
  }
  segment.resize(static_cast<std::size_t>(N_segment) + 1);

  // --- construct the concave-patch decomposition (no meshing here) ---
  decomp_construct_concavepat(pd, probe, i0, I, I_probe, N_probe, segment,
                              N_segment, circle_interior, N_circle, dir_i0);
  pd.a_i = static_cast<int32_t>(a_i);
  pd.a_j = static_cast<int32_t>(a_j);
  pd.a_k = static_cast<int32_t>(a_k);
  pd.ci = ci;
  pd.cj = cj;
  pd.ck = ck;
}

// ----- build_neighbors ------------------------------------------------------
struct NeighborsOut {
  std::vector<std::vector<int>> neighbor_I;  // [i] ragged list (0-based vector)
  std::vector<int> N_neighbor;               // 1-based-with-dummy
};
NeighborsOut build_neighbors(int nhight, const std::vector<int>& hight_set,
                             const std::vector<int>& inverse_hight,
                             const std::vector<Vec3>& I,
                             const std::vector<std::array<int32_t, 3>>& Iijk,
                             const DataI& di, const Neighbors& inter, double Rp) {
  NeighborsOut out;
  out.neighbor_I.assign(static_cast<std::size_t>(nhight) + 1, {});
  out.N_neighbor.assign(static_cast<std::size_t>(nhight) + 1, 0);
  // Each probe i writes only its own fixed slots (neighbor_I[i], N_neighbor[i])
  // and reads const inputs -> the loop is parallel with no merge step at all;
  // per-slot contents (first-seen order included) are unchanged.
  meshms::parallel_for(1, nhight + 1, [&](int i) {
    Vec3 x1 = I[static_cast<std::size_t>(hight_set[static_cast<std::size_t>(i)])];
    int i1 = Iijk[static_cast<std::size_t>(hight_set[static_cast<std::size_t>(i)])][0];
    // seen set seeded with 0 (first-seen order preserved).
    std::vector<int> seen{0};
    auto contains = [&](int v) {
      for (int s : seen) if (s == v) return true;
      return false;
    };
    std::vector<int>& row = out.neighbor_I[static_cast<std::size_t>(i)];
    int num_int_i1 = inter.count(i1);
    for (int j = 1; j <= num_int_i1; ++j) {
      int j1 = static_cast<int>(inter.of(i1)[static_cast<std::size_t>(j - 1)]);
      int In_j1 = static_cast<int>(di.Ii[static_cast<std::size_t>(j1)].size());
      for (int k = 1; k <= In_j1; ++k) {
        int k1 = static_cast<int>(di.Ii[static_cast<std::size_t>(j1)][static_cast<std::size_t>(k - 1)]);
        Vec3 x2 = I[static_cast<std::size_t>(k1)];
        int vid = inverse_hight[static_cast<std::size_t>(k1)];
        if (k1 != hight_set[static_cast<std::size_t>(i)] && norm(x1 - x2) < 2 * Rp &&
            !contains(vid)) {
          row.push_back(vid);
          seen.push_back(vid);
        }
      }
    }
    out.N_neighbor[static_cast<std::size_t>(i)] = static_cast<int>(row.size());
  });
  return out;
}

}  // namespace

// ----- precompute_concave (density-independent half) -------------------------
ConcaveDecomp precompute_concave(const Geom& geom, const DataI& di, double Rp,
                                 const Neighbors& inter) {
  int s = di.nI;
  const std::vector<Vec3>& I = di.I;
  const std::vector<std::array<int32_t, 3>>& Iijk = di.Iijk;
  const std::vector<Vec3>& C = geom.centers;
  const std::vector<int32_t>& hight = di.high_I;
  const std::vector<std::array<int32_t, 3>>& direction = di.direction;

  // ======================== cSAS case (arg_eSAS=0) ======================== //
  std::vector<int> hight_set(static_cast<std::size_t>(s) + 1, 0);
  int nhight = 0;
  std::vector<int> inverse_hight(static_cast<std::size_t>(s) + 1, 0);
  for (int i = 1; i <= s; ++i) {
    if (hight[static_cast<std::size_t>(i)] == 0) {
      nhight += 1;
      hight_set[static_cast<std::size_t>(nhight)] = i;
      inverse_hight[static_cast<std::size_t>(i)] = nhight;
    }
  }
  hight_set.resize(static_cast<std::size_t>(nhight) + 1);

  NeighborsOut nbo = build_neighbors(nhight, hight_set, inverse_hight, I, Iijk, di,
                                     inter, Rp);

  // PARALLEL: each probe's decomposition is independent (call-local scratch)
  // and lands in its own fixed slot -- no merge step, order-independent.
  ConcaveDecomp dec;
  dec.nhight = nhight;
  dec.probes.assign(static_cast<std::size_t>(nhight) + 1, ProbePatchSet{});
  meshms::parallel_for(1, nhight + 1, [&](int i) {
    int Kn = nbo.N_neighbor[static_cast<std::size_t>(i)];
    std::vector<int> K_in(static_cast<std::size_t>(Kn) + 1, 0);
    for (int t = 1; t <= Kn; ++t) {
      K_in[static_cast<std::size_t>(t)] = nbo.neighbor_I[static_cast<std::size_t>(i)][static_cast<std::size_t>(t - 1)];
    }
    data_concavepat(dec.probes[static_cast<std::size_t>(i)], i, hight_set, K_in,
                    Kn, I, Iijk, C, Rp, direction);
  });
  return dec;
}

// ----- SESconcavepat_mesh (density-dependent half) ---------------------------
void SESconcavepat_mesh(MeshState& state, const Geom& geom, const DataI& di,
                        const ConcaveDecomp& decomp, double Rp, double d) {
  int s = di.nI;
  const std::vector<Vec3>& I = di.I;
  const std::vector<std::array<int32_t, 3>>& Iijk = di.Iijk;
  const std::vector<Vec3>& C = geom.centers;
  const std::vector<double>& R = geom.R;
  const std::vector<int32_t>& hight = di.high_I;
  const std::vector<std::array<int32_t, 3>>& direction = di.direction;

  const int nhight = decomp.nhight;

  // PARALLEL S7 (cSAS pass): each probe meshes its precomputed patch set
  // independently -> a per-iteration vector<LocalMesh>, then a SERIAL ordered
  // merge add_patch's every patch in probe-then-patch order.
  std::vector<std::vector<LocalMesh>> probe_lm(static_cast<std::size_t>(nhight) + 1);
  meshms::parallel_for(1, nhight + 1, [&](int i) {
    mesh_probe_patchset(probe_lm[static_cast<std::size_t>(i)],
                        decomp.probes[static_cast<std::size_t>(i)], I, Rp, d);
  });
  {
    std::size_t add_v = 0, add_f = 0;
    for (int i = 1; i <= nhight; ++i) {
      for (const LocalMesh& lm : probe_lm[static_cast<std::size_t>(i)]) {
        if (lm.emit) {
          add_v += lm.P.empty() ? 0 : lm.P.size() - 1;
          add_f += lm.T.size();
        }
      }
    }
    state.reserve_extra(add_v, add_f);
  }
  for (int i = 1; i <= nhight; ++i) {
    for (LocalMesh& lm : probe_lm[static_cast<std::size_t>(i)]) {
      if (lm.emit) state.add_patch(lm.P, lm.T, lm.NV, std::move(lm.vids), lm.vatom);
    }
  }

  // ===================== simple concave triangles ======================== //

  // PARALLEL S7 (simple-triangle pass): each i with hight[i]==1 emits exactly
  // ONE LocalMesh; all scratch (loops/segment0/...) is call-local. A fixed
  // 1-slot-per-i layout (emit=false where hight[i]!=1) keeps the merge order
  // deterministic (ascending i), then a SERIAL ordered add_patch merge.
  std::vector<LocalMesh> tri_lm(static_cast<std::size_t>(s) + 1);
  meshms::parallel_for(1, s + 1, [&](int i) {
    if (hight[static_cast<std::size_t>(i)] == 1) {
      Vec3 c = I[static_cast<std::size_t>(i)];
      int a_i = Iijk[static_cast<std::size_t>(i)][0];
      int a_j = Iijk[static_cast<std::size_t>(i)][1];
      int a_k = Iijk[static_cast<std::size_t>(i)][2];
      Vec3 ci = C[static_cast<std::size_t>(a_i)];
      Vec3 cj = C[static_cast<std::size_t>(a_j)];
      Vec3 ck = C[static_cast<std::size_t>(a_k)];
      Vec3 ni = (ci - c) / norm(ci - c);
      Vec3 nj = (cj - c) / norm(cj - c);
      Vec3 nk = (ck - c) / norm(ck - c);
      Vec3 nij = cross(ni, nj);
      Vec3 njk = cross(nj, nk);
      Vec3 nki = cross(nk, ni);
      nij = nij / norm(nij);
      njk = njk / norm(njk);
      nki = nki / norm(nki);
      Vec3 spoint_i = c + (ci - c) * Rp / (R[static_cast<std::size_t>(a_i)] + Rp);
      Vec3 spoint_j = c + (cj - c) * Rp / (R[static_cast<std::size_t>(a_j)] + Rp);
      Vec3 spoint_k = c + (ck - c) * Rp / (R[static_cast<std::size_t>(a_k)] + Rp);
      double angle_i = acos_clamped(dot(spoint_i - c, spoint_j - c) / pysq(Rp));
      double angle_j = acos_clamped(dot(spoint_j - c, spoint_k - c) / pysq(Rp));
      double angle_k = acos_clamped(dot(spoint_k - c, spoint_i - c) / pysq(Rp));

      // loops as 1-based-with-dummy matrix: single loop of size 3.
      std::vector<Loop> loops(2);
      loops[0] = Loop(1);  // dummy row 0
      loops[1] = Loop{0, 0, 0, 0};  // [0] dummy + 3 entries
      std::vector<int> loopsize(2, 0);
      loopsize[1] = 3;
      if (direction[static_cast<std::size_t>(i)][0] == -1) {
        nij = -1.0 * nij;
        njk = -1.0 * njk;
        nki = -1.0 * nki;
        Vec3 si = spoint_i;
        spoint_i = spoint_j;
        spoint_j = spoint_k;
        spoint_k = si;
        loops[1][1] = 1;
        loops[1][2] = 3;
        loops[1][3] = 2;
      } else {
        loops[1][1] = 1;
        loops[1][2] = 2;
        loops[1][3] = 3;
      }

      // segment0 record matrix [_, c=1:3, n=4:6, r=7, spoint=8:10, angle=11].
      std::vector<std::array<double, 12>> segment0(4, std::array<double, 12>{});
      rset3(segment0[1], 1, c);
      rset3(segment0[1], 4, nij);
      segment0[1][7] = Rp;
      rset3(segment0[1], 8, spoint_i);
      segment0[1][11] = angle_i;
      rset3(segment0[2], 1, c);
      rset3(segment0[2], 4, njk);
      segment0[2][7] = Rp;
      rset3(segment0[2], 8, spoint_j);
      segment0[2][11] = angle_j;
      rset3(segment0[3], 1, c);
      rset3(segment0[3], 4, nki);
      segment0[3][7] = Rp;
      rset3(segment0[3], 8, spoint_k);
      segment0[3][11] = angle_k;

      // patches: single patch = loop 1, patchesize 1.
      std::vector<int> patches_row{0, 1};  // [0] dummy, [1]=1
      int patchesize_i = 1;

      // circle0 empty (Python passes None).
      std::vector<std::array<double, 9>> circle0_empty;
      Tag btag{2, i, 0};  // ("probe", i)
      tri_lm[static_cast<std::size_t>(i)] =
          mesh_sphpat(c, Rp, loops, segment0, circle0_empty, patches_row,
                      patchesize_i, Rp, d, nullptr, btag);
      // Attribute each triangle vertex to the nearest of the 3 atoms a_i/a_j/a_k.
      fill_vatom_nearest3(tri_lm[static_cast<std::size_t>(i)],
                          static_cast<int32_t>(a_i), static_cast<int32_t>(a_j),
                          static_cast<int32_t>(a_k), ci, cj, ck);
    }
  });
  {
    std::size_t add_v = 0, add_f = 0;
    for (int i = 1; i <= s; ++i) {
      const LocalMesh& lm = tri_lm[static_cast<std::size_t>(i)];
      if (lm.emit) {
        add_v += lm.P.empty() ? 0 : lm.P.size() - 1;
        add_f += lm.T.size();
      }
    }
    state.reserve_extra(add_v, add_f);
  }
  for (int i = 1; i <= s; ++i) {
    LocalMesh& lm = tri_lm[static_cast<std::size_t>(i)];
    if (lm.emit) state.add_patch(lm.P, lm.T, lm.NV, std::move(lm.vids), lm.vatom);
  }
}

// ----- SESconcavepat (driver) -----------------------------------------------
void SESconcavepat(MeshState& state, const Geom& geom, const DataI& di,
                   const Ext& ext, double Rp, double d, const Neighbors& inter) {
  (void)ext;  // ext_I (eSAS flags) is read only in the deferred av-only passes.
  // Volume of cSAS/eSAS and the eSAS interior-trim pass are av-only -> skipped.
  ConcaveDecomp dec = precompute_concave(geom, di, Rp, inter);
  SESconcavepat_mesh(state, geom, di, dec, Rp, d);
}

}  // namespace meshms
