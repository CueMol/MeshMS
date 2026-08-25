// Faithful port of the weld module (minus manifold_report, which lives in
// mesh_check.hpp). Python dicts / defaultdicts preserve INSERTION order, and the
// component/cycle/flap walks depend on that order, so every map here is backed by
// an insertion-ordered structure (a vector of entries + a lookup index) and every
// per-vertex / per-edge adjacency list keeps first-seen order. This reproduces
// the Python topology and (via %.17g coords) the exact welded coordinates.
#include "meshms/weld.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace meshms {

namespace {

// 64-bit hash for an integer key triple (used as the weld bucket key).
struct Int3 {
  int64_t a, b, c;
  bool operator==(const Int3& o) const { return a == o.a && b == o.b && c == o.c; }
};
struct Int3Hash {
  std::size_t operator()(const Int3& k) const {
    std::size_t h = std::hash<int64_t>{}(k.a);
    h ^= std::hash<int64_t>{}(k.b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h ^= std::hash<int64_t>{}(k.c) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
  }
};

// Sorted undirected edge key.
inline std::pair<int, int> ekey(int u, int v) {
  return (u < v) ? std::make_pair(u, v) : std::make_pair(v, u);
}
struct EdgeHash {
  std::size_t operator()(const std::pair<int, int>& e) const {
    return std::hash<int64_t>{}((static_cast<int64_t>(e.first) << 32) ^
                                static_cast<uint32_t>(e.second));
  }
};

// np.round on a scalar == banker's rounding (round-half-to-even). std::nearbyint
// honours the current rounding mode, which is FE_TONEAREST by default == ties to
// even. We compute V/tol first (the Python writes V / tol), then round, exactly
// matching np.round(V / tol).astype(np.int64).
inline int64_t round_key(double coord, double tol) {
  return static_cast<int64_t>(std::nearbyint(coord / tol));
}

}  // namespace

WeldResult weld(const std::vector<Vec3>& V, const std::vector<Tri>& F, double tol) {
  WeldResult out;
  const std::size_t n = V.size();

  // round(V/tol) integer-key spatial hash; first vertex into a bucket is the
  // representative, new indices follow first-seen order (Python range(n)).
  std::unordered_map<Int3, int, Int3Hash> bucket;
  bucket.reserve(n * 2 + 1);
  std::vector<int> old_to_new(n, 0);
  for (std::size_t old = 0; old < n; ++old) {
    const Vec3& v = V[old];
    Int3 key{round_key(v.x, tol), round_key(v.y, tol), round_key(v.z, tol)};
    auto it = bucket.find(key);
    int nw;
    if (it == bucket.end()) {
      nw = static_cast<int>(out.V.size());
      bucket.emplace(key, nw);
      out.V.push_back(v);  // representative = first vertex hashing here
    } else {
      nw = it->second;
    }
    old_to_new[old] = nw;
  }

  // remap faces; drop degenerate (any two welded indices coincide).
  for (const Tri& f : F) {
    const int na = old_to_new[static_cast<std::size_t>(f[0])];
    const int nb = old_to_new[static_cast<std::size_t>(f[1])];
    const int nc = old_to_new[static_cast<std::size_t>(f[2])];
    if (na == nb || nb == nc || nc == na) continue;
    out.F.push_back(Tri{na, nb, nc});
  }
  return out;
}

namespace {

// Element adapters so the facade-layout overload (facade MeshResult arrays) runs
// the identical loop without a whole-mesh Vec3/Tri conversion.
inline Vec3 blv_get(const std::vector<Vec3>& V, std::size_t i) { return V[i]; }
inline Vec3 blv_get(const std::vector<std::array<double, 3>>& V, std::size_t i) {
  return Vec3{V[i][0], V[i][1], V[i][2]};
}
inline int bli_get(const Tri& f, int k) { return f[static_cast<std::size_t>(k)]; }
inline int bli_get(const std::array<std::uint32_t, 3>& f, int k) {
  return static_cast<int>(f[static_cast<std::size_t>(k)]);
}

template <class VT, class FT>
BoundaryLoopsResult boundary_loops_impl(const std::vector<VT>& V,
                                        const std::vector<FT>& F) {
  BoundaryLoopsResult res;

  // Insertion-ordered undirected edge counts (key -> count), keyed by first-seen.
  std::vector<std::pair<int, int>> edge_order;
  std::vector<int> edge_cnt;
  std::unordered_map<std::pair<int, int>, int, EdgeHash> edge_idx;
  edge_idx.reserve(F.size() * 6 + 1);
  for (const FT& f : F) {
    const int a = bli_get(f, 0), b = bli_get(f, 1), c = bli_get(f, 2);
    const std::pair<int, int> es[3] = {ekey(a, b), ekey(b, c), ekey(c, a)};
    for (const auto& e : es) {
      auto it = edge_idx.find(e);
      if (it == edge_idx.end()) {
        edge_idx.emplace(e, static_cast<int>(edge_order.size()));
        edge_order.push_back(e);
        edge_cnt.push_back(1);
      } else {
        ++edge_cnt[static_cast<std::size_t>(it->second)];
      }
    }
  }

  // boundary edges (count==1) in first-seen order; nonmanifold (count>=3).
  std::vector<std::pair<int, int>> boundary;
  for (std::size_t i = 0; i < edge_order.size(); ++i) {
    const int n = edge_cnt[i];
    if (n == 1) {
      boundary.push_back(edge_order[i]);
    } else if (n >= 3) {
      res.nonmanifold.push_back(
          NonmanifoldEdge{edge_order[i].first, edge_order[i].second, n});
    }
  }

  // undirected adjacency over the boundary subgraph, insertion-ordered per vertex
  // (defaultdict(set): a vertex first appears in the order its first boundary
  // edge is added; the set membership only matters for the closed test / count).
  std::vector<int> vert_order;                  // vertices in first-seen order
  std::unordered_map<int, int> vert_idx;        // vertex -> slot in adj
  std::vector<std::vector<int>> adj;            // neighbour sets (insertion order)
  auto slot_of = [&](int x) -> int {
    auto it = vert_idx.find(x);
    if (it != vert_idx.end()) return it->second;
    int s = static_cast<int>(vert_order.size());
    vert_idx.emplace(x, s);
    vert_order.push_back(x);
    adj.emplace_back();
    return s;
  };
  auto add_adj = [&](int x, int y) {
    int s = slot_of(x);
    auto& nb = adj[static_cast<std::size_t>(s)];
    if (std::find(nb.begin(), nb.end(), y) == nb.end()) nb.push_back(y);
  };
  for (const auto& e : boundary) {
    add_adj(e.first, e.second);
    add_adj(e.second, e.first);
  }

  // connected components of the boundary-edge subgraph, iterating vertices in
  // first-seen order (Python `for s0 in adj`), DFS via a LIFO stack.
  std::unordered_map<int, int> comp_id;
  std::vector<std::vector<int>> comps;
  for (int s0 : vert_order) {
    if (comp_id.count(s0)) continue;
    const int cid = static_cast<int>(comps.size());
    std::vector<int> members;
    std::vector<int> stack{s0};
    while (!stack.empty()) {
      int x = stack.back();
      stack.pop_back();
      if (comp_id.count(x)) continue;
      comp_id.emplace(x, cid);
      members.push_back(x);
      // stack.extend(y for y in adj[x] if y not in comp_id) -- in adj order.
      for (int y : adj[static_cast<std::size_t>(vert_idx[x])]) {
        if (!comp_id.count(y)) stack.push_back(y);
      }
    }
    comps.push_back(std::move(members));
  }

  // boundary edges per component (each boundary edge counted via its u endpoint).
  std::vector<int> edges_per(comps.size(), 0);
  for (const auto& e : boundary) {
    ++edges_per[static_cast<std::size_t>(comp_id[e.first])];
  }

  for (std::size_t cid = 0; cid < comps.size(); ++cid) {
    const auto& members = comps[cid];
    bool closed = true;
    for (int x : members) {
      if (adj[static_cast<std::size_t>(vert_idx[x])].size() != 2) {
        closed = false;
        break;
      }
    }
    BoundaryLoop L;
    L.n_verts = static_cast<int>(members.size());
    L.n_edges = edges_per[cid];
    Vec3 sum{};
    for (int x : members) sum = sum + blv_get(V, static_cast<std::size_t>(x));
    L.centroid = sum / static_cast<double>(members.size());
    L.closed = closed;
    L.vids = members;
    std::sort(L.vids.begin(), L.vids.end());
    res.loops.push_back(std::move(L));
  }

  // sort largest-first by n_edges (stable: ties keep first-seen component order,
  // matching Python list.sort(key=lambda L: -L["n_edges"])).
  std::stable_sort(res.loops.begin(), res.loops.end(),
                   [](const BoundaryLoop& x, const BoundaryLoop& y) {
                     return x.n_edges > y.n_edges;
                   });
  return res;
}

}  // namespace

BoundaryLoopsResult boundary_loops(const std::vector<Vec3>& V, const std::vector<Tri>& F) {
  return boundary_loops_impl(V, F);
}

BoundaryLoopsResult boundary_loops(const std::vector<std::array<double, 3>>& V,
                                   const std::vector<std::array<std::uint32_t, 3>>& F) {
  return boundary_loops_impl(V, F);
}

FillResult fill_small_holes(std::vector<Vec3> V, const std::vector<Tri>& F,
                            int max_loop) {
  FillResult out;
  out.F = F;
  if (F.empty()) {
    out.V = std::move(V);
    return out;
  }

  // undirected boundary edges + the (unnormalised) normal of each edge's adjacent
  // face (last writer wins, like the Python dict assignment).
  std::vector<std::pair<int, int>> edge_order;
  std::vector<int> edge_cnt;
  std::unordered_map<std::pair<int, int>, int, EdgeHash> edge_idx;
  std::unordered_map<std::pair<int, int>, Vec3, EdgeHash> edge_face_normal;
  edge_idx.reserve(F.size() * 6 + 1);
  for (const Tri& f : F) {
    const int a = f[0], b = f[1], c = f[2];
    const Vec3 fn = cross(V[static_cast<std::size_t>(b)] - V[static_cast<std::size_t>(a)],
                          V[static_cast<std::size_t>(c)] - V[static_cast<std::size_t>(a)]);
    const std::pair<int, int> es[3] = {ekey(a, b), ekey(b, c), ekey(c, a)};
    for (const auto& e : es) {
      auto it = edge_idx.find(e);
      if (it == edge_idx.end()) {
        edge_idx.emplace(e, static_cast<int>(edge_order.size()));
        edge_order.push_back(e);
        edge_cnt.push_back(1);
      } else {
        ++edge_cnt[static_cast<std::size_t>(it->second)];
      }
      edge_face_normal[e] = fn;  // last writer wins
    }
  }
  std::unordered_set<std::pair<int, int>, EdgeHash> boundary;
  for (std::size_t i = 0; i < edge_order.size(); ++i) {
    if (edge_cnt[i] == 1) boundary.insert(edge_order[i]);
  }
  if (boundary.empty()) {
    out.V = std::move(V);
    return out;
  }

  // undirected adjacency (for components / simple cycles), insertion-ordered;
  // Python uses list (NOT set) here -> duplicates allowed, so do not dedup.
  std::vector<int> vert_order;
  std::unordered_map<int, int> vert_idx;
  std::vector<std::vector<int>> adj;
  auto slot_of = [&](int x) -> int {
    auto it = vert_idx.find(x);
    if (it != vert_idx.end()) return it->second;
    int s = static_cast<int>(vert_order.size());
    vert_idx.emplace(x, s);
    vert_order.push_back(x);
    adj.emplace_back();
    return s;
  };
  // NB: Python iterates `for u, v in boundary` over a SET, so adj insertion order
  // is the set's iteration order. adj is only used for component membership, the
  // `closed` degree test, and simple-cycle walking; for a simple cycle every
  // vertex has exactly two neighbours so the walk is order-independent. We build
  // adjacency in first-seen (edge_order) order -- deterministic and equivalent.
  for (std::size_t i = 0; i < edge_order.size(); ++i) {
    if (edge_cnt[i] != 1) continue;
    const auto& e = edge_order[i];
    adj[static_cast<std::size_t>(slot_of(e.first))].push_back(e.second);
    adj[static_cast<std::size_t>(slot_of(e.second))].push_back(e.first);
  }

  // directed half-edges along boundary edges, in F order (Python `out`).
  std::unordered_map<int, std::vector<int>> dout;
  for (const Tri& f : F) {
    const int a = f[0], b = f[1], c = f[2];
    const std::pair<int, int> dirs[3] = {{a, b}, {b, c}, {c, a}};
    for (const auto& d : dirs) {
      if (boundary.count(ekey(d.first, d.second))) {
        dout[d.first].push_back(d.second);
      }
    }
  }

  auto fan_fill = [&](const std::vector<int>& cyc) {
    const int m = static_cast<int>(cyc.size());
    if (!(3 <= m && m <= max_loop)) return;
    Vec3 nref{0.0, 0.0, 0.0};
    for (int i = 0; i < m; ++i) {
      const int a = cyc[static_cast<std::size_t>(i)];
      const int b = cyc[static_cast<std::size_t>((i + 1) % m)];
      auto it = edge_face_normal.find(ekey(a, b));
      if (it != edge_face_normal.end()) {
        const Vec3& fn = it->second;
        const double ln = std::sqrt(dot(fn, fn));
        if (ln > 0.0) nref = nref + fn / ln;
      }
    }
    if (nref.x == 0.0 && nref.y == 0.0 && nref.z == 0.0) {
      nref = Vec3{0.0, 0.0, 1.0};
    }
    const int v0 = cyc[0];
    for (int i = 1; i < m - 1; ++i) {
      const int a = cyc[static_cast<std::size_t>(i)];
      const int b = cyc[static_cast<std::size_t>(i + 1)];
      const Vec3 tn = cross(V[static_cast<std::size_t>(a)] - V[static_cast<std::size_t>(v0)],
                            V[static_cast<std::size_t>(b)] - V[static_cast<std::size_t>(v0)]);
      if (dot(tn, nref) >= 0.0) {
        out.F.push_back(Tri{v0, a, b});
      } else {
        out.F.push_back(Tri{v0, b, a});
      }
    }
  };

  // process each connected boundary component (vertices in first-seen order).
  std::unordered_set<int> seen;
  for (int s : vert_order) {
    if (seen.count(s)) continue;
    std::vector<int> comp;
    std::vector<int> stack{s};
    while (!stack.empty()) {
      int x = stack.back();
      stack.pop_back();
      if (seen.count(x)) continue;
      seen.insert(x);
      comp.push_back(x);
      for (int y : adj[static_cast<std::size_t>(vert_idx[x])]) {
        if (!seen.count(y)) stack.push_back(y);
      }
    }
    // n_edges = sum(len(adj[x]) for x in comp) // 2
    int deg_sum = 0;
    for (int x : comp) deg_sum += static_cast<int>(adj[static_cast<std::size_t>(vert_idx[x])].size());
    const int n_edges = deg_sum / 2;
    if (!(3 <= n_edges && n_edges <= max_loop)) continue;

    bool all_deg2 = true;
    for (int x : comp) {
      if (adj[static_cast<std::size_t>(vert_idx[x])].size() != 2) {
        all_deg2 = false;
        break;
      }
    }

    if (all_deg2) {
      // simple cycle: walk it undirected (orientation handled by nref).
      std::vector<int> cyc{comp[0]};
      int prev = -1, cur = comp[0];
      bool ok = true;
      bool has_prev = false;
      for (int step = 0; step < static_cast<int>(comp.size()) - 1; ++step) {
        int next = -1;
        for (int y : adj[static_cast<std::size_t>(vert_idx[cur])]) {
          if (!has_prev || y != prev) {
            next = y;
            break;
          }
        }
        if (next == -1) {
          ok = false;
          break;
        }
        prev = cur;
        has_prev = true;
        cur = next;
        cyc.push_back(cur);
      }
      if (ok && static_cast<int>(cyc.size()) == static_cast<int>(comp.size())) {
        fan_fill(cyc);
      }
    } else {
      // pinched component: split into directed cycles via the half-edges.
      std::unordered_set<int> comp_set(comp.begin(), comp.end());
      std::unordered_set<std::pair<int, int>, EdgeHash> used;
      for (int u0 : comp) {
        auto outit = dout.find(u0);
        if (outit == dout.end()) continue;
        for (int v0 : outit->second) {
          if (used.count({u0, v0})) continue;
          std::vector<int> cyc{u0};
          used.insert({u0, v0});
          int cur = v0;
          bool ok = true;
          while (cur != u0) {
            cyc.push_back(cur);
            if (static_cast<int>(cyc.size()) > max_loop) {
              ok = false;
              break;
            }
            // nbrs = [v for v in out[cur] if (cur,v) not in used and v in comp_set]
            int nv = -1;
            auto cit = dout.find(cur);
            if (cit != dout.end()) {
              for (int v : cit->second) {
                if (!used.count({cur, v}) && comp_set.count(v)) {
                  nv = v;
                  break;
                }
              }
            }
            if (nv == -1) {
              ok = false;
              break;
            }
            used.insert({cur, nv});
            cur = nv;
          }
          if (ok && cur == u0) fan_fill(cyc);
        }
      }
    }
  }
  out.V = std::move(V);
  return out;
}

FlapResult remove_nonmanifold_flaps(std::vector<Vec3> V, const std::vector<Tri>& F,
                                    int passes, std::vector<std::uint32_t>* kept_orig) {
  FlapResult out;
  // V is purely topological here: moved straight through to the result.
  out.V = std::move(V);
  std::vector<Tri> cur = F;
  // Original-face indices carried alongside `cur` so a parallel per-face array can
  // be filtered to the survivors. Identity until a pass actually drops faces.
  std::vector<std::uint32_t> idx(cur.size());
  for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = static_cast<std::uint32_t>(i);
  if (F.empty()) {
    out.F = cur;
    if (kept_orig) kept_orig->clear();
    return out;
  }

  for (int pass = 0; pass < passes; ++pass) {
    // vfaces[v] = face indices touching v (insertion order); edge_faces[e] = same.
    std::unordered_map<int, std::vector<int>> vfaces;
    std::vector<std::pair<int, int>> edge_order;
    std::unordered_map<std::pair<int, int>, std::vector<int>, EdgeHash> edge_faces;
    for (std::size_t fi = 0; fi < cur.size(); ++fi) {
      const int a = cur[fi][0], b = cur[fi][1], c = cur[fi][2];
      for (int x : {a, b, c}) vfaces[x].push_back(static_cast<int>(fi));
      const std::pair<int, int> es[3] = {ekey(a, b), ekey(b, c), ekey(c, a)};
      for (const auto& e : es) {
        auto it = edge_faces.find(e);
        if (it == edge_faces.end()) {
          edge_faces.emplace(e, std::vector<int>{static_cast<int>(fi)});
          edge_order.push_back(e);
        } else {
          it->second.push_back(static_cast<int>(fi));
        }
      }
    }

    // non-manifold edges in first-seen order (Python dict iteration order).
    std::vector<std::pair<int, int>> nm_edges;
    for (const auto& e : edge_order) {
      if (edge_faces[e].size() >= 3) nm_edges.push_back(e);
    }
    if (nm_edges.empty()) break;

    std::unordered_set<int> drop;
    for (const auto& e : nm_edges) {
      for (int fi : edge_faces[e]) {
        // third = the apex vertex of this face not on the edge.
        int third = -1;
        for (int x : {cur[static_cast<std::size_t>(fi)][0],
                      cur[static_cast<std::size_t>(fi)][1],
                      cur[static_cast<std::size_t>(fi)][2]}) {
          if (x != e.first && x != e.second) {
            third = x;
            break;
          }
        }
        // apex hangs ONLY off this edge -> the whole fan at `third` is a flap.
        bool all_on_edge = true;
        for (int gi : vfaces[third]) {
          const Tri& g = cur[static_cast<std::size_t>(gi)];
          const bool has_u = (e.first == g[0] || e.first == g[1] || e.first == g[2]);
          const bool has_v = (e.second == g[0] || e.second == g[1] || e.second == g[2]);
          if (!(has_u && has_v)) {
            all_on_edge = false;
            break;
          }
        }
        if (all_on_edge) {
          for (int gi : vfaces[third]) drop.insert(gi);
        }
      }
    }
    if (drop.empty()) break;

    std::vector<Tri> kept;
    std::vector<std::uint32_t> kept_idx;
    kept.reserve(cur.size());
    kept_idx.reserve(cur.size());
    for (std::size_t i = 0; i < cur.size(); ++i) {
      if (!drop.count(static_cast<int>(i))) {
        kept.push_back(cur[i]);
        kept_idx.push_back(idx[i]);
      }
    }
    cur = std::move(kept);
    idx = std::move(kept_idx);
  }

  out.F = std::move(cur);
  if (kept_orig) *kept_orig = std::move(idx);
  return out;
}

}  // namespace meshms
