// Cross-check ID-based boundary fusion (build_surface, fuse=True) against the
// Python golden <mol>.fused.txt.
//
// The golden is the full build_surface mesh AFTER orient_faces AND fuse_by_id
// (warn=False): fused vertices V2 and remapped 0-based faces F2 (no normals --
// N is dropped when fused). The golden floats are %.17g (full IEEE-double
// round-trip), so each parsed golden double equals the Python double EXACTLY.
//
// PASS criterion: nV/nF EXACT, F2 exact 0-based topology, V2 within 1e-9.
// ADDITIONALLY (reported): manifold_report boundary_edges BEFORE (fuse=False)
// and AFTER (fuse=True) fusion -- the fused mesh must have FEWER boundary edges
// (the C1 seams close), checked here as a regression guard.
//
// check_fuse_reference also compares fuse_by_id with a copy of its original
// hash-grid port on random tagged point sets; it needs no golden, so it runs
// under every FP policy.
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "meshms/fusion.hpp"
#include "meshms/mesh_check.hpp"
#include "meshms/pipeline.hpp"
#include "meshms/vec3.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

struct Golden {
  int nV = 0, nF = 0;
  std::vector<Vec3> V;                // 0-based store (V2[k-1])
  std::vector<std::array<int, 3>> F;  // 0-based vertex indices (post-fusion)
};

Golden read_golden(const std::string& path) {
  std::ifstream in(path);
  CHECK(static_cast<bool>(in));
  Golden g;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream ls(line);
    std::string tag;
    ls >> tag;
    if (tag == "nV") {
      ls >> g.nV;
    } else if (tag == "nF") {
      ls >> g.nF;
    } else if (tag == "V") {
      int k = 0;
      Vec3 v{};
      ls >> k >> v.x >> v.y >> v.z;
      g.V.push_back(v);
    } else if (tag == "F") {
      int k = 0;
      std::array<int, 3> t{};
      ls >> k >> t[0] >> t[1] >> t[2];
      g.F.push_back(t);
    }
  }
  return g;
}

void check_mol(const std::string& mol) {
  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
  const std::string refpath =
      std::string(MESHMS_REF_DIR) + "/" + mol + ".fused.txt";

  // Unfused mesh (for the boundary-edge before/after regression guard).
  Surface raw = build_surface(xyzr, 1.4, 0.5, /*fuse=*/false);
  ManifoldReport mr_before = manifold_report(raw.V, raw.F);

  // Fused mesh -- the subject under test.
  Surface s = build_surface(xyzr, 1.4, 0.5, /*fuse=*/true);
  Golden gold = read_golden(refpath);

  const int gotNV = static_cast<int>(s.V.size());
  const int gotNF = static_cast<int>(s.F.size());
  CHECK(gotNV == gold.nV);
  CHECK(gotNF == gold.nF);
  if (gotNV != gold.nV)
    std::fprintf(stderr, "  %s nV got=%d want=%d\n", mol.c_str(), gotNV, gold.nV);
  if (gotNF != gold.nF)
    std::fprintf(stderr, "  %s nF got=%d want=%d\n", mol.c_str(), gotNF, gold.nF);

  // fuse=True drops normals (Python returns N=None).
  CHECK(s.N.empty());

  // --- vertices: 1e-9 --------------------------------------------------------
  const int nV = std::min(gotNV, gold.nV);
  bool vfirst = true;
  for (int k = 0; k < nV; ++k) {
    const Vec3& got = s.V[static_cast<std::size_t>(k)];
    const Vec3& want = gold.V[static_cast<std::size_t>(k)];
    const bool near = std::fabs(got.x - want.x) <= 1e-9 &&
                      std::fabs(got.y - want.y) <= 1e-9 &&
                      std::fabs(got.z - want.z) <= 1e-9;
    CHECK(near);
    if (!near && vfirst) {
      std::fprintf(stderr,
                   "  %s V[%d] got=(%.17g,%.17g,%.17g) want=(%.17g,%.17g,%.17g)\n",
                   mol.c_str(), k + 1, got.x, got.y, got.z, want.x, want.y, want.z);
      vfirst = false;
    }
  }

  // --- faces: EXACT 0-based topology (post-fusion) ---------------------------
  const int nF = std::min(gotNF, gold.nF);
  bool ffirst = true;
  for (int t = 0; t < nF; ++t) {
    const Tri& f = s.F[static_cast<std::size_t>(t)];
    const std::array<int, 3>& w = gold.F[static_cast<std::size_t>(t)];
    const bool ok = (f[0] == w[0]) && (f[1] == w[1]) && (f[2] == w[2]);
    CHECK(ok);
    if (!ok && ffirst) {
      std::fprintf(stderr, "  %s F[%d] got=(%d,%d,%d) want=(%d,%d,%d)\n",
                   mol.c_str(), t + 1, f[0], f[1], f[2], w[0], w[1], w[2]);
      ffirst = false;
    }
  }

  // --- boundary-edge regression guard: fusion must CLOSE seams ---------------
  ManifoldReport mr_after = manifold_report(s.V, s.F);
  std::printf(
      "  %-9s nV %d->%d nF %d->%d  boundary_edges %d->%d  area=%.6f\n",
      mol.c_str(), static_cast<int>(raw.V.size()), gotNV,
      static_cast<int>(raw.F.size()), gotNF, mr_before.boundary_edges,
      mr_after.boundary_edges, mr_after.area);
  CHECK(mr_after.boundary_edges < mr_before.boundary_edges);
}

// The original fuse_by_id port: a (tag, cell) -> [vertex] hash grid probed over
// the 27 neighbour cells of every tagged vertex, unions replayed in vertex
// order, roots mapped through a hash map. Kept as the reference the rewrite
// must reproduce exactly.
struct RefKey {
  int kind, i, j;
  long long cx, cy, cz;
  bool operator==(const RefKey& o) const {
    return kind == o.kind && i == o.i && j == o.j && cx == o.cx && cy == o.cy && cz == o.cz;
  }
};
struct RefKeyHash {
  std::size_t operator()(const RefKey& k) const {
    std::size_t h = 1469598103934665603ULL;
    for (long long v : {static_cast<long long>(k.kind), static_cast<long long>(k.i),
                        static_cast<long long>(k.j), k.cx, k.cy, k.cz})
      h = (h ^ static_cast<std::size_t>(v)) * 1099511628211ULL;
    return h;
  }
};

std::tuple<std::vector<Vec3>, std::vector<Tri>, std::vector<int32_t>> ref_fuse(
    const std::vector<Vec3>& V, const std::vector<Tri>& F, const std::vector<TagList>& tags,
    const std::vector<int32_t>& atom_id, double eps, std::vector<std::uint32_t>& kept) {
  const std::size_t n = V.size();
  std::vector<long long> parent(n);
  for (std::size_t i = 0; i < n; ++i) parent[i] = static_cast<long long>(i);
  auto find = [&](long long x) {
    while (parent[static_cast<std::size_t>(x)] != x) x = parent[static_cast<std::size_t>(x)];
    return x;
  };
  auto unite = [&](long long a, long long b) {
    long long ra = find(a), rb = find(b);
    if (ra != rb) parent[static_cast<std::size_t>(std::max(ra, rb))] = std::min(ra, rb);
  };
  const double inv = 1.0 / eps;
  auto cell = [&](const Vec3& v) {
    return std::array<long long, 3>{static_cast<long long>(std::floor(v.x * inv)),
                                    static_cast<long long>(std::floor(v.y * inv)),
                                    static_cast<long long>(std::floor(v.z * inv))};
  };
  std::unordered_map<RefKey, std::vector<long long>, RefKeyHash> grid;
  for (std::size_t i = 0; i < n; ++i) {
    if (i >= tags.size()) continue;
    const auto c = cell(V[i]);
    for (const Tag& t : tags[i])
      grid[RefKey{t.kind, t.i, t.j, c[0], c[1], c[2]}].push_back(static_cast<long long>(i));
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (i >= tags.size()) continue;
    const auto c = cell(V[i]);
    for (const Tag& t : tags[i])
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
          for (int dz = -1; dz <= 1; ++dz) {
            auto it = grid.find(RefKey{t.kind, t.i, t.j, c[0] + dx, c[1] + dy, c[2] + dz});
            if (it == grid.end()) continue;
            for (long long j : it->second) {
              if (j >= static_cast<long long>(i)) continue;
              Vec3 dv = V[i] - V[static_cast<std::size_t>(j)];
              if (dot(dv, dv) < eps * eps) unite(static_cast<long long>(i), j);
            }
          }
  }
  std::unordered_map<long long, long long> roots;
  std::vector<long long> o2n(n);
  std::vector<Vec3> V2;
  std::vector<int32_t> A2;
  for (std::size_t i = 0; i < n; ++i) {
    const long long r = find(static_cast<long long>(i));
    auto it = roots.find(r);
    if (it == roots.end()) {
      it = roots.emplace(r, static_cast<long long>(V2.size())).first;
      V2.push_back(V[static_cast<std::size_t>(r)]);
      A2.push_back(atom_id[static_cast<std::size_t>(r)]);
    }
    o2n[i] = it->second;
  }
  std::vector<Tri> F2;
  kept.clear();
  for (std::size_t fi = 0; fi < F.size(); ++fi) {
    const long long a = o2n[static_cast<std::size_t>(F[fi][0])];
    const long long b = o2n[static_cast<std::size_t>(F[fi][1])];
    const long long c = o2n[static_cast<std::size_t>(F[fi][2])];
    if (a == b || b == c || c == a) continue;
    F2.push_back(Tri{static_cast<int32_t>(a), static_cast<int32_t>(b), static_cast<int32_t>(c)});
    kept.push_back(static_cast<std::uint32_t>(fi));
  }
  return {std::move(V2), std::move(F2), std::move(A2)};
}

// Random tagged point clouds built to stress the candidate search: clusters of
// copies within ~eps of each other (some just inside, some just outside eps),
// seeded on cell boundaries so copies straddle cells, vertices carrying two
// tags or one tag twice, untagged vertices, and a tag list shorter than V.
void check_fuse_reference() {
  const double eps = 1e-6;
  int total_pairs_fused = 0;
  for (unsigned seed = 1; seed <= 20; ++seed) {
    std::mt19937 rng(seed);
    auto unit = [&]() { return static_cast<double>(rng() % 1000000) / 1000000.0; };
    std::vector<Vec3> V;
    std::vector<TagList> tags;
    std::vector<int32_t> aid;
    for (int c = 0; c < 400; ++c) {
      // On a cell boundary half of the time, so the copies straddle cells.
      const double bx = static_cast<double>(rng() % 2000) * ((rng() % 2) ? eps : 1.37e-3);
      const Vec3 base{bx, unit() * 3.0, unit() * 3.0};
      const int ncopy = 1 + static_cast<int>(rng() % 4);
      for (int k = 0; k < ncopy; ++k) {
        const double s = eps * ((rng() % 3 == 0) ? 0.9 : 0.4);
        V.push_back(base + Vec3{s * (unit() - 0.5), s * (unit() - 0.5), s * (unit() - 0.5)});
        TagList tl;
        const int r = static_cast<int>(rng() % 10);
        if (r > 0) tl.push_back(Tag{static_cast<int>(rng() % 3), c % 7, 0});
        if (r > 6) tl.push_back(Tag{static_cast<int>(rng() % 3), c % 5, 1});
        if (r == 9) tl.push_back(tl.front());
        tags.push_back(tl);
        aid.push_back(static_cast<int32_t>(rng() % 50));
      }
    }
    tags.resize(tags.size() - 3);  // tags shorter than V: the tail reads as None
    std::vector<Tri> F;
    const int nv = static_cast<int>(V.size());
    for (int f = 0; f < 3 * nv; ++f)
      F.push_back(Tri{static_cast<int>(rng() % nv), static_cast<int>(rng() % nv),
                      static_cast<int>(rng() % nv)});

    std::vector<std::uint32_t> ref_kept, got_kept;
    auto [rV, rF, rA] = ref_fuse(V, F, tags, aid, eps, ref_kept);
    auto [gV, gF, gA] = fuse_by_id(V, F, tags, aid, eps, &got_kept);
    bool vsame = gV.size() == rV.size();
    for (std::size_t k = 0; vsame && k < gV.size(); ++k)
      vsame = gV[k].x == rV[k].x && gV[k].y == rV[k].y && gV[k].z == rV[k].z;
    CHECK(vsame);
    CHECK(gF == rF);
    CHECK(gA == rA);
    CHECK(got_kept == ref_kept);
    total_pairs_fused += static_cast<int>(V.size() - rV.size());
  }
  // The clouds must actually merge vertices, or the comparison proves nothing.
  CHECK(total_pairs_fused > 1000);
  std::printf("  fuse-ref 20 seeds match, %d vertices merged\n", total_pairs_fused);
}

}  // namespace

int main() {
  check_fuse_reference();
  // The golden part below needs strict FP; a fast build still reports the
  // reference check above instead of skipping the whole test.
  if (kFpFast) {
    std::printf("SKIP: golden part requires MESHMS_FP=strict\n");
    TEST_MAIN_RETURN();
  }
  check_mol("2spheres");
  check_mol("3spheres");
  check_mol("tetra");
  check_mol("ArgArg");
  TEST_MAIN_RETURN();
}
