// Faithful port of the fusion module's fuse_by_id (warn=False pipeline path).
#include "meshms/fusion.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

#include "meshms/parallel.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

namespace {

// Spatial cell index (floor(v * inv) per axis). 64-bit components match the
// Python arbitrary-precision int(np.floor(...)) for any coordinate we meet.
struct Cell {
  std::int64_t x, y, z;
};

// cell = (int(floor(v[0]*inv)), int(floor(v[1]*inv)), int(floor(v[2]*inv))).
inline Cell cell_of(const Vec3& v, double inv) {
  return Cell{static_cast<std::int64_t>(std::floor(v.x * inv)),
              static_cast<std::int64_t>(std::floor(v.y * inv)),
              static_cast<std::int64_t>(std::floor(v.z * inv))};
}

}  // namespace

std::tuple<std::vector<Vec3>, std::vector<Tri>, std::vector<int32_t>> fuse_by_id(
    const std::vector<Vec3>& V, const std::vector<Tri>& F,
    const std::vector<TagList>& tags, const std::vector<int32_t>& atom_id,
    double eps, std::vector<std::uint32_t>* kept_faces) {
  const std::size_t n = V.size();
  const bool have_aid = !atom_id.empty();

  // --- union-find (parent[max(ra,rb)] = min(ra,rb); find with full path
  //     compression to root, matching the Python `parent[x], x = root, parent[x]`) ---
  std::vector<std::int64_t> parent(n);
  for (std::size_t i = 0; i < n; ++i) parent[i] = static_cast<std::int64_t>(i);

  auto find = [&](std::int64_t x) -> std::int64_t {
    std::int64_t root = x;
    while (parent[static_cast<std::size_t>(root)] != root)
      root = parent[static_cast<std::size_t>(root)];
    // full path compression to root: parent[x], x = root, parent[x]
    while (parent[static_cast<std::size_t>(x)] != root) {
      std::int64_t next = parent[static_cast<std::size_t>(x)];
      parent[static_cast<std::size_t>(x)] = root;
      x = next;
    }
    return root;
  };
  auto unite = [&](std::int64_t a, std::int64_t b) {
    std::int64_t ra = find(a), rb = find(b);
    if (ra != rb) {
      std::int64_t hi = ra > rb ? ra : rb;
      std::int64_t lo = ra < rb ? ra : rb;
      parent[static_cast<std::size_t>(hi)] = lo;
    }
  };

  const double eps2 = eps * eps;
  const double inv_eps = 1.0 / eps;

  // SORT-AND-SWEEP match (the hash-grid probe it replaces was most of the mesh
  // time on a large molecule). Every (tag, cell, vertex) entry is sorted, so a
  // vertex only meets the later entries of the same tag whose x cell is at most
  // one above its own; a pair outside that window cannot share a tag with cells
  // within +-1 per axis. The match predicate is the old one -- shared tag, cell
  // within +-1 per axis, dot(dv, dv) < eps^2 with dv = (later - earlier vertex)
  // -- and the final union-find partition is order-independent:
  //   * the predicate does not depend on union-find state;
  //   * uniting the two ends of every matching pair in any order yields the
  //     connected components of the match-pair graph;
  //   * parent[max_root] = min_root keeps every component's root at its minimum
  //     member regardless of union order, so V2's first-seen order and the
  //     root's atom id are identical.
  struct Entry {
    Tag tag;
    Cell cell;
    std::int64_t idx;
  };
  // tags[i] is None (interior) when the TagList is empty; tags shorter than V
  // also reads as None past its end (Python: ti = tags[i] if i < len(tags)).
  std::vector<std::size_t> eoff(n + 1, 0);
  for (std::size_t i = 0; i < n; ++i)
    eoff[i + 1] = eoff[i] + ((i < tags.size()) ? tags[i].size() : 0);
  std::vector<Entry> ent(eoff[n]);
  meshms::parallel_for(0, static_cast<int>(n), [&](int ii) {
    const std::size_t i = static_cast<std::size_t>(ii);
    if (i >= tags.size() || tags[i].empty()) return;
    const Cell ce = cell_of(V[i], inv_eps);
    std::size_t k = eoff[i];
    for (const Tag& tag : tags[i]) ent[k++] = Entry{tag, ce, static_cast<std::int64_t>(i)};
  });
  // A total order (idx last), so the sorted array is the same on every run.
  meshms::parallel_sort(ent.begin(), ent.end(), [](const Entry& a, const Entry& b) {
    if (a.tag.kind != b.tag.kind) return a.tag.kind < b.tag.kind;
    if (a.tag.i != b.tag.i) return a.tag.i < b.tag.i;
    if (a.tag.j != b.tag.j) return a.tag.j < b.tag.j;
    if (a.cell.x != b.cell.x) return a.cell.x < b.cell.x;
    if (a.cell.y != b.cell.y) return a.cell.y < b.cell.y;
    if (a.cell.z != b.cell.z) return a.cell.z < b.cell.z;
    return a.idx < b.idx;
  });

  // Sweep in fixed blocks, each collecting its own pairs (a block's scan may
  // read past its end; the array is read-only here).
  const std::size_t ne = ent.size();
  const int nblk = static_cast<int>(std::min<std::size_t>(1024, ne / 4096 + 1));
  std::vector<std::vector<std::pair<std::int64_t, std::int64_t>>> blk_pairs(
      static_cast<std::size_t>(nblk));
  meshms::parallel_for(0, nblk, [&](int bk) {
    const std::size_t b0 = ne * static_cast<std::size_t>(bk) / static_cast<std::size_t>(nblk);
    const std::size_t b1 =
        ne * static_cast<std::size_t>(bk + 1) / static_cast<std::size_t>(nblk);
    auto& pairs = blk_pairs[static_cast<std::size_t>(bk)];
    for (std::size_t p = b0; p < b1; ++p) {
      const Entry& ep = ent[p];
      for (std::size_t q = p + 1; q < ne; ++q) {
        const Entry& eq = ent[q];
        if (eq.tag != ep.tag || eq.cell.x > ep.cell.x + 1) break;
        if (eq.cell.y < ep.cell.y - 1 || eq.cell.y > ep.cell.y + 1 ||
            eq.cell.z < ep.cell.z - 1 || eq.cell.z > ep.cell.z + 1)
          continue;
        if (eq.idx == ep.idx) continue;  // one vertex listing a tag twice
        const std::int64_t i = std::max(ep.idx, eq.idx);
        const std::int64_t j = std::min(ep.idx, eq.idx);
        Vec3 dv = V[static_cast<std::size_t>(i)] - V[static_cast<std::size_t>(j)];
        if (dot(dv, dv) < eps2) pairs.emplace_back(i, j);
      }
    }
  });
  for (const auto& pairs : blk_pairs)
    for (const auto& pr : pairs) unite(pr.first, pr.second);

  // --- deduplicated vertex list from union-find roots (first-seen order) ------
  // Roots are vertex indices, so a flat array maps each root to its new index.
  std::vector<std::int64_t> root_new(n, -1);
  std::vector<std::int64_t> o2n(n);
  std::vector<Vec3> V2;
  std::vector<std::int32_t> atom_id2;  // empty unless have_aid (root's atom id)
  for (std::size_t i = 0; i < n; ++i) {
    const std::int64_t r = find(static_cast<std::int64_t>(i));
    std::int64_t& g = root_new[static_cast<std::size_t>(r)];
    if (g < 0) {
      g = static_cast<std::int64_t>(V2.size());
      V2.push_back(V[static_cast<std::size_t>(r)]);
      // The kept vertex is the root r; carry its atom id (a neighbour atom).
      if (have_aid) atom_id2.push_back(atom_id[static_cast<std::size_t>(r)]);
    }
    o2n[i] = g;
  }

  // --- remap faces, dropping degenerate (two equal indices) -------------------
  std::vector<Tri> F2;
  F2.reserve(F.size());
  for (std::size_t fi = 0; fi < F.size(); ++fi) {
    const Tri& f = F[fi];
    std::int64_t na = o2n[static_cast<std::size_t>(f[0])];
    std::int64_t nb = o2n[static_cast<std::size_t>(f[1])];
    std::int64_t nc = o2n[static_cast<std::size_t>(f[2])];
    if (na == nb || nb == nc || nc == na) continue;
    F2.push_back(Tri{static_cast<int32_t>(na), static_cast<int32_t>(nb),
                     static_cast<int32_t>(nc)});
    if (kept_faces) kept_faces->push_back(static_cast<std::uint32_t>(fi));
  }

  return {std::move(V2), std::move(F2), std::move(atom_id2)};
}

// 2-output overload: delegate to the 3-output core with no atom ids and drop the
// (empty) atom_id2. Keeps the existing pipeline / test_fusion call sites working.
std::pair<std::vector<Vec3>, std::vector<Tri>> fuse_by_id(
    const std::vector<Vec3>& V, const std::vector<Tri>& F,
    const std::vector<TagList>& tags, double eps) {
  auto [V2, F2, aid2] = fuse_by_id(V, F, tags, std::vector<int32_t>{}, eps);
  (void)aid2;
  return {std::move(V2), std::move(F2)};
}

}  // namespace meshms
