// Faithful port of the sas module's data_I_Cir(). Float term/eval order is kept
// identical to the numpy expressions (via the vec3 helpers) so the triple-point
// topology and circle set match the Python to the last branch. Storage follows
// the append / per-atom CSR / jagged layout -- the OOM-wall fix.
//
// PARALLELISM (S6, ParallelSAS): the per-atom i-loop is split into two phases.
// PHASE 1 is '#pragma omp parallel for' over i=1..M and is fully READ-ONLY on
// shared state -- it reads only the immutable {C, R, nb, Rp} and writes ONLY its
// own per-atom event buffer ev[i] (no shared write -> race-free). For each atom
// it runs the EXACT current logic but, instead of mutating I/Ii/I_circle/s, it
// APPENDS events to ev[i] in PRODUCTION ORDER (Y1 before Y2 within the same
// (row1,k); the free CircleEvent at the end). PHASE 2 is SERIAL over i=1..M and
// does ZERO float math: it replays ev[i] in order, performing EXACTLY the same
// push_back / ++s / In/Ii / I_circle append sequence as the original serial
// loop. Because phase 2 visits atoms in the same 1..M order and replays each
// atom's events in recorded order, the s-numbering and every append order are
// reproduced BIT-FOR-BIT.
#include "meshms/sas.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace meshms {

namespace {

// A recorded triple point for phase-2 replay. Captures everything phase 2 needs
// so it performs no float math: the geometric data (Y, dir, hightvalue,
// true_hight) and the three append targets (row1, row3, row2; j and k for the
// Ii lists). One TripleEvent per accepted Y (true1 -> Y1, true2 -> Y2), pushed
// in the same order as the original push_back sequence.
struct TripleEvent {
  int j;
  int k;
  int row1;
  int row2;
  int row3;
  Vec3 Y;
  std::array<int32_t, 3> dir;
  double hightvalue;
  int true_hight;  // 0 -> high_I=0, else high_I=1
};

// A recorded free circle (circletest==1 after the row2 loop) for phase-2 replay.
struct CircleEvent {
  int j;
  Vec3 A;
  Vec3 nij;
  double rij;
};

// Per-atom event buffer written ONLY by the owning thread in phase 1.
struct AtomEvents {
  std::vector<TripleEvent> triples;
  std::vector<CircleEvent> circles;
};

}  // namespace

std::pair<DataI, DataCir> data_I_Cir(const Geom& geom, const Neighbors& nb, double Rp) {
  const int M = geom.M;
  const std::vector<Vec3>& C = geom.centers;  // size M+1, [0] dummy
  const std::vector<double>& R = geom.R;      // size M+1, [0] dummy

  // Row[i] -> nb.count(i); M_int[i,row] -> nb.of(i)[row-1] (1-based row in
  // 1..Row[i]; the CSR row is sorted ascending-unique, matching the sort at the
  // top of the Python data_I_Cir).
  auto Row = [&](int i) { return nb.count(i); };
  auto Mint = [&](int i, int row) { return nb.of(i)[static_cast<std::size_t>(row - 1)]; };

  DataI di;
  DataCir dc;

  // maxs is ONLY a reserve() hint (the maxs-family are append vectors): the
  // provable triple-point upper bound sum_i Row[i]*(Row[i]-1)/3 + 1.
  long long sum = 0;
  for (int i = 1; i <= M; ++i) {
    const long long r = Row(i);
    sum += r * (r - 1);
  }
  const std::size_t maxs = static_cast<std::size_t>(sum / 3) + 1;

  // ---------------------------------------------------------------------------
  // PHASE 1: parallel, READ-ONLY. Each thread fills its own ev[i] (no shared
  // write). The body is the exact serial logic with mutation replaced by
  // appends to ev[i] in production order.
  // ---------------------------------------------------------------------------
  std::vector<AtomEvents> ev(static_cast<std::size_t>(M) + 1);

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 1; i <= M; ++i) {
    AtomEvents& evi = ev[static_cast<std::size_t>(i)];
    const int Ri = Row(i);
    for (int row1 = 1; row1 <= Ri; ++row1) {
      const int j = Mint(i, row1);
      if (j <= i) continue;

      int circletest = 1;
      Vec3 A = circlecenter(C[i], C[j], R[i] + Rp, R[j] + Rp);
      // discriminant of the intersection-circle radius; <= 0 means one SAS ball
      // contains the other (or tangent) -> no real circle.
      const double disc = pysq(R[i] + Rp) - pysq(norm(C[i] - A));
      if (disc <= 0) continue;

      Vec3 normal = C[j] - C[i];
      const double dist_normal = norm(normal);
      Vec3 nij = normal / dist_normal;
      const double rij = std::sqrt(disc);

      // Test if the circle is entirely covered by a sphere Sk.
      for (int row = 1; row <= Ri; ++row) {
        const int k = Mint(i, row);
        if (k != j) {
          Vec3 a = A - C[k];
          const double dist_a = norm(a);
          const double ctheta = dot(a, normal) / (dist_a * dist_normal);
          if (pysq(std::sqrt(1 - pysq(ctheta)) * dist_a + rij) +
                  pysq(ctheta * dist_a) - pysq(R[k] + Rp) <
              0) {
            circletest = 0;
            break;
          }
        }
      }

      if (circletest == 1) {
        for (int row2 = 1; row2 <= Ri; ++row2) {
          const int k = Mint(i, row2);

          if (k < j && circletest == 1) {  // test if circle is a whole circle
            Vec3 a = A - C[k];
            const double dist_a = norm(a);
            const double ctheta = dot(a, normal) / (dist_a * dist_normal);
            if (pysq(-std::sqrt(1 - pysq(ctheta)) * dist_a + rij) +
                    pysq(ctheta * dist_a) - pysq(R[k] + Rp) <
                0) {
              circletest = 0;
            }
          }

          // 1-based column of k in atom j's neighbor row (row3) via binary
          // search: the CSR row nb.of(j) is sorted ascending-unique (csr.hpp),
          // so lower_bound is an exact membership test. row3 == 0 means k <= j
          // or k not a neighbor of j -- both cases the original condition
          // skipped. The integer column index is bit-identical to the map.
          int row3 = 0;
          if (k > j) {
            const std::span<const int32_t> jrow = nb.of(j);
            auto it = std::lower_bound(jrow.begin(), jrow.end(), k);
            if (it != jrow.end() && *it == k) {
              row3 = static_cast<int>(it - jrow.begin()) + 1;
            }
          }
          if (row3) {
            Vec3 B = circlecenter(C[k], C[j], R[k] + Rp, R[j] + Rp);
            Vec3 cij = C[i] - C[j];
            Vec3 ckj = C[k] - C[j];
            Vec3 cr = cross(cij, ckj);
            Vec3 n = cr / norm(cr);
            // u = ckj - np.dot(ckj, cij) * cij / np.dot(cij, cij)
            Vec3 u = ckj - dot(ckj, cij) * cij / dot(cij, cij);
            const double t = dot(B - A, ckj) / dot(u, ckj);
            Vec3 X1 = A + t * u;
            const double c = -(pysq(norm(X1 - C[i]))) + pysq(R[i] + Rp);

            if (c > 0) {
              circletest = 0;

              const double s1 = -std::sqrt(c);
              const double s2 = std::sqrt(c);
              Vec3 Y1 = X1 + s1 * n;
              Vec3 Y2 = X1 + s2 * n;

              int true_hight = 1;
              if (s2 < Rp) true_hight = 0;

              int true1 = 1;
              for (int row = 1; row <= Ri; ++row) {
                const int ll = Mint(i, row);
                const double dist1 = norm(C[ll] - Y1);
                if (dist1 < R[ll] + Rp && ll != j && ll != k) {
                  true1 = 0;
                  break;
                }
              }
              int true2 = 1;
              for (int row = 1; row <= Ri; ++row) {
                const int ll = Mint(i, row);
                const double dist2 = norm(C[ll] - Y2);
                if (dist2 < R[ll] + Rp && ll != j && ll != k) {
                  true2 = 0;
                  break;
                }
              }

              // Record Y1 then Y2 (same order as the original push_back blocks).
              if (true1 == 1) {
                evi.triples.push_back(TripleEvent{
                    j, k, row1, row2, row3, Y1, {1, 1, -1}, s2, true_hight});
              }
              if (true2 == 1) {
                evi.triples.push_back(TripleEvent{
                    j, k, row1, row2, row3, Y2, {-1, -1, 1}, s2, true_hight});
              }
            }
          }
        }
      }

      if (circletest == 1) {
        evi.circles.push_back(CircleEvent{j, A, nij, rij});
      }
    }
  }

  // ---------------------------------------------------------------------------
  // PHASE 2: SERIAL, atoms i=1..M in order, ZERO float math. Replay ev[i] in
  // order, performing the exact original mutation sequence. This reproduces the
  // serial s-numbering and every append order BIT-FOR-BIT.
  // ---------------------------------------------------------------------------

  // s is 1-based: index 0 of every maxs-family vector is the dummy row.
  di.I.assign(1, Vec3{});
  di.Iijk.assign(1, std::array<int32_t, 3>{0, 0, 0});
  di.direction.assign(1, std::array<int32_t, 3>{0, 0, 0});
  di.high_I.assign(1, 1);  // dummy row 0 == 1 (Python high_I = np.ones(...))
  di.hightvalue.assign(1, 0.0);
  di.I.reserve(maxs + 1);
  di.Iijk.reserve(maxs + 1);
  di.direction.reserve(maxs + 1);
  di.high_I.reserve(maxs + 1);
  di.hightvalue.reserve(maxs + 1);

  // Ii[a] = point ids on atom a (per-atom growable; In[a] == Ii[a].size()).
  di.Ii.assign(static_cast<std::size_t>(M) + 1, {});

  // I_circle[i] sized Row[i]+1 ([0] dummy); I_circle[i][row] grows as points
  // land on the (i,row)-th intersection circle. I_circle_num == .size().
  di.I_circle.assign(static_cast<std::size_t>(M) + 1, {});
  for (int i = 1; i <= M; ++i) {
    di.I_circle[static_cast<std::size_t>(i)].assign(
        static_cast<std::size_t>(Row(i)) + 1, {});
  }

  int s = 0;

  // circle rows accumulate as {_, i, j, A(3), nij(3), rij}; circle[0] dummy.
  // reserve() is only a capacity hint -- it never changes contents/order. A
  // free circle exists for at most one (i,j) directed neighbor pair with j>i,
  // so half the total CSR neighbour count is a safe (loose) upper bound.
  std::vector<std::array<double, 10>> circle_rows;
  {
    long long nb_total = 0;
    for (int i = 1; i <= M; ++i) nb_total += Row(i);
    circle_rows.reserve(static_cast<std::size_t>(nb_total / 2) + 1);
  }

  for (int i = 1; i <= M; ++i) {
    const AtomEvents& evi = ev[static_cast<std::size_t>(i)];

    for (const TripleEvent& e : evi.triples) {
      const int j = e.j;
      const int k = e.k;
      ++s;
      di.I.push_back(e.Y);
      di.Iijk.push_back({i, j, k});
      di.direction.push_back(e.dir);
      di.Ii[static_cast<std::size_t>(i)].push_back(s);
      di.Ii[static_cast<std::size_t>(j)].push_back(s);
      di.Ii[static_cast<std::size_t>(k)].push_back(s);
      di.I_circle[static_cast<std::size_t>(i)][static_cast<std::size_t>(e.row1)].push_back(s);
      di.I_circle[static_cast<std::size_t>(j)][static_cast<std::size_t>(e.row3)].push_back(s);
      di.I_circle[static_cast<std::size_t>(i)][static_cast<std::size_t>(e.row2)].push_back(s);
      di.hightvalue.push_back(e.hightvalue);
      di.high_I.push_back(e.true_hight == 0 ? 0 : 1);
    }

    for (const CircleEvent& e : evi.circles) {
      circle_rows.push_back({0.0, static_cast<double>(i), static_cast<double>(e.j),
                             e.A.x, e.A.y, e.A.z, e.nij.x, e.nij.y, e.nij.z, e.rij});
    }
  }

  di.nI = s;

  // circle[0] dummy, then circle_rows in order.
  const int ncircle = static_cast<int>(circle_rows.size());
  dc.ncircle = ncircle;
  dc.circle.assign(1, std::array<double, 10>{});
  dc.circle.reserve(static_cast<std::size_t>(ncircle) + 1);
  for (const auto& rrow : circle_rows) dc.circle.push_back(rrow);

  // circleindex[a] = free-circle ids touching atom a (ci then cj per circle).
  // Per-atom reserve hint (Row(a)): a free circle touching atom a corresponds
  // to one of a's neighbour pairs, so its count is bounded by Row(a). reserve()
  // is capacity-only and leaves the appended order untouched.
  dc.circleindex.assign(static_cast<std::size_t>(M) + 1, {});
  for (int a = 1; a <= M; ++a) {
    dc.circleindex[static_cast<std::size_t>(a)].reserve(
        static_cast<std::size_t>(Row(a)));
  }
  for (int k = 1; k <= ncircle; ++k) {
    const int ci = static_cast<int>(dc.circle[static_cast<std::size_t>(k)][1]);
    const int cj = static_cast<int>(dc.circle[static_cast<std::size_t>(k)][2]);
    dc.circleindex[static_cast<std::size_t>(ci)].push_back(k);
    dc.circleindex[static_cast<std::size_t>(cj)].push_back(k);
  }

  return {std::move(di), std::move(dc)};
}

}  // namespace meshms
