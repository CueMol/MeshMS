// Faithful port of the intersection module's interstructure().
//
// Differences vs Python (semantics identical):
//   - The result is emitted DIRECTLY as a CSR (Neighbors): off is a prefix sum
//     of per-atom counts, val is the concatenation of the ascending neighbour
//     lists. We never materialise the dense (M+1, kmax+1) M_int rectangle
//     (the OOM-wall fix).
//   - Cell binning uses std::floor then cast to int64 (NEVER int() truncation),
//     matching np.floor(...).astype(np.int64) for negative coordinates.
#include "meshms/intersection.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace meshms {

namespace {

// Hash a 3D integer cell key (i, j, k) for the spatial-hash grid.
struct CellKey {
  std::int64_t i, j, k;
  bool operator==(const CellKey& o) const { return i == o.i && j == o.j && k == o.k; }
};

struct CellKeyHash {
  std::size_t operator()(const CellKey& c) const noexcept {
    // Mix three int64 with a 64-bit splitmix-style combine.
    std::uint64_t h = static_cast<std::uint64_t>(c.i) * 0x9E3779B97F4A7C15ull;
    h ^= static_cast<std::uint64_t>(c.j) + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
    h ^= static_cast<std::uint64_t>(c.k) + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
    return static_cast<std::size_t>(h);
  }
};

}  // namespace

Neighbors interstructure(const Geom& geom, double Rp) {
  const int M = geom.M;
  const std::vector<Vec3>& C = geom.centers;  // size M+1, [0] dummy
  const std::vector<double>& R = geom.R;      // size M+1, [0] dummy

  Neighbors nb;
  nb.M = M;
  nb.off.assign(static_cast<std::size_t>(M) + 2, 0);  // size M+2, all zero
  if (M <= 0) {
    return nb;  // empty (off all zero, val empty)
  }

  // box side = max SAS radius over real atoms 1..M, plus probe.
  // float(R[1:M+1].max()) + Rp
  double rmaxR = R[1];
  for (int a = 2; a <= M; ++a) {
    if (R[a] > rmaxR) rmaxR = R[a];
  }
  const double rmax = rmaxR + Rp;

  // Bin each atom into EXACTLY ONE cell: floor(center / rmax) cast to int64.
  // (np.floor(...).astype(np.int64))
  std::vector<std::int64_t> cix(static_cast<std::size_t>(M) + 1);
  std::vector<std::int64_t> ciy(static_cast<std::size_t>(M) + 1);
  std::vector<std::int64_t> ciz(static_cast<std::size_t>(M) + 1);

  std::unordered_map<CellKey, std::vector<int>, CellKeyHash> cells;
  cells.reserve(static_cast<std::size_t>(M) * 2);
  for (int a = 1; a <= M; ++a) {
    const std::int64_t ix = static_cast<std::int64_t>(std::floor(C[a].x / rmax));
    const std::int64_t iy = static_cast<std::int64_t>(std::floor(C[a].y / rmax));
    const std::int64_t iz = static_cast<std::int64_t>(std::floor(C[a].z / rmax));
    cix[a] = ix;
    ciy[a] = iy;
    ciz[a] = iz;
    cells[CellKey{ix, iy, iz}].push_back(a);
  }

  // First pass: per-atom neighbour lists (sorted ascending), and counts.
  // We build the rows so we can prefix-sum off, then fill val.
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(M) + 1);

  // PARALLEL (output-preserving): each atom a builds its OWN rows[a] reading only
  // the immutable cells grid + geom (const find on the unordered_map is a safe
  // concurrent read); the row is sorted, so its content/order is deterministic.
  // The serial prefix-sum + concatenation below emits the CSR in atom order, so
  // the result is bit-identical to the serial build.
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for (int a = 1; a <= M; ++a) {
    const Vec3 ca = C[a];
    const double ra = R[a];
    const std::int64_t ix = cix[a];
    const std::int64_t iy = ciy[a];
    const std::int64_t iz = ciz[a];

    std::vector<int>& found = rows[static_cast<std::size_t>(a)];
    // +/-2 cell stencil; di,dj,dk in -2..2. The 125 stencil cells are distinct,
    // and each atom is binned into exactly one cell, so no candidate b is ever
    // visited twice -> the row is naturally unique once sorted.
    for (int di = -2; di <= 2; ++di) {
      for (int dj = -2; dj <= 2; ++dj) {
        for (int dk = -2; dk <= 2; ++dk) {
          auto it = cells.find(CellKey{ix + di, iy + dj, iz + dk});
          if (it == cells.end()) continue;
          for (int b : it->second) {
            if (b == a) continue;
            // dist(c_a, c_b): same term order as the numpy port (x,y,z).
            const double dx = C[b].x - ca.x;
            const double dy = C[b].y - ca.y;
            const double dz = C[b].z - ca.z;
            const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (ra + R[b] + 2.0 * Rp - dist > 0.0) {
              found.push_back(b);
            }
          }
        }
      }
    }
    std::sort(found.begin(), found.end());
  }

  // Prefix-sum off: off[i] = start of atom i's row, off[M+1] = total.
  // off[0] = 0 (dummy atom 0 has an empty range).
  std::int32_t total = 0;
  for (int a = 0; a <= M; ++a) {
    nb.off[static_cast<std::size_t>(a)] = total;
    total += static_cast<std::int32_t>(rows[static_cast<std::size_t>(a)].size());
  }
  nb.off[static_cast<std::size_t>(M) + 1] = total;

  // Concatenate ascending neighbour lists into val, asserting STRICT ascent
  // (no duplicates) --- the downstream col_of reverse index depends on it.
  nb.val.resize(static_cast<std::size_t>(total));
  std::int32_t w = 0;
  for (int a = 1; a <= M; ++a) {
    const std::vector<int>& r = rows[static_cast<std::size_t>(a)];
    for (std::size_t t = 0; t < r.size(); ++t) {
      if (t > 0) {
        assert(r[t] > r[t - 1] && "interstructure row must be strictly ascending");
      }
      nb.val[static_cast<std::size_t>(w++)] = static_cast<std::int32_t>(r[t]);
    }
  }

  return nb;
}

}  // namespace meshms
