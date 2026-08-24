// Faithful port of the fusion module's fuse_by_id (warn=False pipeline path).
#include "meshms/fusion.hpp"

#include <cmath>
#include <cstdint>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "meshms/parallel.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

namespace {

// Spatial-hash cell index (floor(v * inv) per axis). 64-bit components match the
// Python arbitrary-precision int(np.floor(...)) for any coordinate we meet.
struct Cell {
  std::int64_t x, y, z;
};
inline bool operator==(const Cell& a, const Cell& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

struct CellHash {
  std::size_t operator()(const Cell& c) const noexcept {
    // 64-bit mix (splitmix-style) folded over the three components.
    auto mix = [](std::uint64_t h) {
      h ^= h >> 30;
      h *= 0xbf58476d1ce4e5b9ULL;
      h ^= h >> 27;
      h *= 0x94d049bb133111ebULL;
      h ^= h >> 31;
      return h;
    };
    std::uint64_t h = mix(static_cast<std::uint64_t>(c.x));
    h ^= mix(static_cast<std::uint64_t>(c.y)) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h ^= mix(static_cast<std::uint64_t>(c.z)) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return static_cast<std::size_t>(h);
  }
};

struct TagHash {
  std::size_t operator()(const Tag& t) const noexcept {
    std::uint64_t h = static_cast<std::uint64_t>(static_cast<std::uint32_t>(t.kind));
    h = h * 0x100000001b3ULL ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(t.i));
    h = h * 0x100000001b3ULL ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(t.j));
    return static_cast<std::size_t>(h);
  }
};

// cell = (int(floor(v[0]*inv)), int(floor(v[1]*inv)), int(floor(v[2]*inv))).
inline Cell cell_of(const Vec3& v, double inv) {
  return Cell{static_cast<std::int64_t>(std::floor(v.x * inv)),
              static_cast<std::int64_t>(std::floor(v.y * inv)),
              static_cast<std::int64_t>(std::floor(v.z * inv))};
}

// Flat (tag, cell) key replacing the old nested tag -> {cell -> [idx]} maps:
// one hash probe per neighbour cell instead of a tag lookup plus a nested cell
// lookup. Exact key equality; the hash only affects speed.
struct TagCell {
  Tag tag;
  Cell cell;
};
inline bool operator==(const TagCell& a, const TagCell& b) {
  return a.tag == b.tag && a.cell == b.cell;
}
struct TagCellHash {
  std::size_t operator()(const TagCell& k) const noexcept {
    std::uint64_t h = TagHash{}(k.tag);
    h ^= CellHash{}(k.cell) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return static_cast<std::size_t>(h);
  }
};

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

  // TWO-PHASE match (replaces the old serial insert-as-you-go walk). The final
  // union-find partition is provably order-independent, which makes the split
  // safe:
  //   * the match predicate (shared tag, cell within +-1, dot(dv,dv) < eps^2)
  //     does not depend on union-find state, and dot(vi-vj, vi-vj) is evaluated
  //     with i as the LATER index -- exactly the pair and expression the old
  //     walk evaluated when the later of the two vertices was processed;
  //   * uniting i with the ROOTS of its matches is equivalent to uniting i with
  //     each matched j (unite re-finds both sides), so the final partition is
  //     the connected components of the match-pair graph either way;
  //   * parent[max_root] = min_root keeps every component's root at its minimum
  //     member regardless of union order, so V2's first-seen order and the
  //     root's atom id are identical.
  //
  // Phase 1 (serial): build the COMPLETE flat (tag, cell) -> [idx] grid.
  std::unordered_map<TagCell, std::vector<std::int64_t>, TagCellHash> grid_eps;
  std::vector<Cell> cell_cache(n);  // cell_of(V[i]) for tagged i, reused below
  {
    std::size_t ntagged = 0;
    for (std::size_t i = 0; i < n && i < tags.size(); ++i)
      if (!tags[i].empty()) ++ntagged;
    grid_eps.reserve(2 * ntagged + 1);
  }
  for (std::size_t i = 0; i < n; ++i) {
    // tags[i] is None (interior) when the TagList is empty; tags shorter than V
    // also reads as None past its end (Python: ti = tags[i] if i < len(tags)).
    if (i >= tags.size() || tags[i].empty()) continue;
    const Cell ce = cell_of(V[i], inv_eps);
    cell_cache[i] = ce;
    for (const Tag& tag : tags[i])
      grid_eps[TagCell{tag, ce}].push_back(static_cast<std::int64_t>(i));
  }

  // NB offsets in EXACT Python order: dx,dy,dz each over (-1, 0, 1).
  static const int NB[27][3] = {
      {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1},
      {-1, 1, -1},  {-1, 1, 0},  {-1, 1, 1},  {0, -1, -1}, {0, -1, 0}, {0, -1, 1},
      {0, 0, -1},   {0, 0, 0},   {0, 0, 1},   {0, 1, -1},  {0, 1, 0},  {0, 1, 1},
      {1, -1, -1},  {1, -1, 0},  {1, -1, 1},  {1, 0, -1},  {1, 0, 0},  {1, 0, 1},
      {1, 1, -1},   {1, 1, 0},   {1, 1, 1}};

  // Phase 2 (parallel per vertex): collect the eps-coincident same-tag partners
  // j < i of every tagged vertex into its own fixed slot (no shared writes; the
  // grid is read-only). Keeping only j < i visits each unordered pair exactly
  // once, from its later endpoint -- the old walk's pair set.
  std::vector<std::vector<std::int64_t>> match(n);
  meshms::parallel_for(0, static_cast<int>(n), [&](int ii) {
    const std::size_t i = static_cast<std::size_t>(ii);
    if (i >= tags.size() || tags[i].empty()) return;
    const Vec3& vi = V[i];
    const Cell ce = cell_cache[i];
    std::vector<std::int64_t>& mi = match[i];
    for (const Tag& tag : tags[i]) {
      for (const auto& off : NB) {
        auto cit = grid_eps.find(
            TagCell{tag, Cell{ce.x + off[0], ce.y + off[1], ce.z + off[2]}});
        if (cit == grid_eps.end()) continue;
        for (std::int64_t j : cit->second) {
          if (j >= static_cast<std::int64_t>(i)) continue;
          Vec3 dv = vi - V[static_cast<std::size_t>(j)];
          double dd = dot(dv, dv);
          if (dd < eps2) mi.push_back(j);
        }
      }
    }
  });

  // Phase 3 (serial): replay the unions in ascending i (any order would yield
  // the same partition -- see above).
  for (std::size_t i = 0; i < n; ++i) {
    for (std::int64_t j : match[i]) unite(static_cast<std::int64_t>(i), j);
  }

  // --- deduplicated vertex list from union-find roots (first-seen order) ------
  std::unordered_map<std::int64_t, std::int64_t> roots;
  std::vector<std::int64_t> o2n(n);
  std::vector<Vec3> V2;
  std::vector<std::int32_t> atom_id2;  // empty unless have_aid (root's atom id)
  for (std::size_t i = 0; i < n; ++i) {
    std::int64_t r = find(static_cast<std::int64_t>(i));
    auto it = roots.find(r);
    std::int64_t g;
    if (it == roots.end()) {
      g = static_cast<std::int64_t>(V2.size());
      roots.emplace(r, g);
      V2.push_back(V[static_cast<std::size_t>(r)]);
      // The kept vertex is the root r; carry its atom id (a neighbour atom).
      if (have_aid) atom_id2.push_back(atom_id[static_cast<std::size_t>(r)]);
    } else {
      g = it->second;
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
