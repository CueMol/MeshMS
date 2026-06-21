// Cross-check the weld.py post-processing routines (weld, boundary_loops,
// remove_nonmanifold_flaps, fill_small_holes) against the Python golden
// <case>.weld.txt.
//
// Each golden file holds, for one case, the raw mesh (RAW), boundary_loops on the
// raw mesh (BL), and the weld -> flaps -> fill chain (WELD, FLAP, FILL) with
// %.17g coords (full IEEE-double round-trip).  Cases:
//   tetra, ArgArg  -- real build_surface(fuse=False) meshes (weld exercised fully);
//   annulus        -- two 6-edge simple cycles fan-filled by fill_small_holes;
//   pinched        -- a degree-4 figure-8 boundary split into two directed cycles;
//   flap           -- a doubled flap dropped by remove_nonmanifold_flaps (then the
//                     exposed perimeter is fan-filled).
// The real meshes come from the C++ build_surface; the synthetic meshes are
// reconstructed here byte-for-byte from dump_weld.py.
//
// PASS criterion: at every stage nV/nF EXACT, F exact 0-based topology, V within
// 1e-9; boundary_loops loop count / per-loop (n_edges,n_verts,closed,vids) and
// the nonmanifold list EXACT, centroid within 1e-9.
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include "meshms/pipeline.hpp"
#include "meshms/vec3.hpp"
#include "meshms/weld.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

struct MeshBlock {
  int nV = 0, nF = 0;
  std::vector<Vec3> V;
  std::vector<Tri> F;
};

struct GoldenLoop {
  int n_edges = 0, n_verts = 0, closed = 0;
  Vec3 centroid{};
  std::vector<int> vids;
};

struct LoopSet {
  int nloops = 0, nnm = 0;
  std::vector<GoldenLoop> loops;
  std::vector<NonmanifoldEdge> nm;
};

struct Golden {
  std::string case_name;
  MeshBlock raw, weld_, flap, fill;
  int fill_nin = 0;   // #faces fill_small_holes received (carried-through prefix)
  LoopSet bl;         // boundary_loops on the RAW mesh
  LoopSet fbl;        // boundary_loops on the FINAL (post-fill) mesh
};

MeshBlock& block_for(Golden& g, const std::string& tag) {
  if (tag == "RAW") return g.raw;
  if (tag == "WELD") return g.weld_;
  if (tag == "FLAP") return g.flap;
  return g.fill;  // "FILL"
}

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
    auto loopset_for = [&](const std::string& pfx) -> LoopSet& {
      return (pfx == "FBL") ? g.fbl : g.bl;
    };
    if (tag == "case") {
      ls >> g.case_name;
    } else if (tag == "FILL_nin") {
      ls >> g.fill_nin;
    } else if (tag == "BL_nloops" || tag == "FBL_nloops") {
      ls >> loopset_for(tag.substr(0, tag.rfind('_'))).nloops;
    } else if (tag == "BL_nnm" || tag == "FBL_nnm") {
      ls >> loopset_for(tag.substr(0, tag.rfind('_'))).nnm;
    } else if (tag == "BL_loop" || tag == "FBL_loop") {
      int i;
      GoldenLoop L;
      ls >> i >> L.n_edges >> L.n_verts >> L.closed >> L.centroid.x >>
          L.centroid.y >> L.centroid.z;
      loopset_for(tag.substr(0, tag.rfind('_'))).loops.push_back(L);
    } else if (tag == "BL_vids" || tag == "FBL_vids") {
      int i, n;
      ls >> i >> n;
      auto& set = loopset_for(tag.substr(0, tag.rfind('_')));
      for (int k = 0; k < n; ++k) {
        int v;
        ls >> v;
        set.loops[static_cast<std::size_t>(i)].vids.push_back(v);
      }
    } else if (tag == "BL_nm" || tag == "FBL_nm") {
      int i;
      NonmanifoldEdge e;
      ls >> i >> e.u >> e.v >> e.count;
      loopset_for(tag.substr(0, tag.rfind('_'))).nm.push_back(e);
    } else {
      // mesh blocks: <TAG>_nV / <TAG>_nF / <TAG>_V / <TAG>_F
      const auto us = tag.rfind('_');
      if (us == std::string::npos) continue;
      const std::string base = tag.substr(0, us);
      const std::string kind = tag.substr(us + 1);
      if (base != "RAW" && base != "WELD" && base != "FLAP" && base != "FILL")
        continue;
      MeshBlock& b = block_for(g, base);
      if (kind == "nV") {
        ls >> b.nV;
      } else if (kind == "nF") {
        ls >> b.nF;
      } else if (kind == "V") {
        int k;
        Vec3 v{};
        ls >> k >> v.x >> v.y >> v.z;
        b.V.push_back(v);
      } else if (kind == "F") {
        int k;
        Tri t{};
        ls >> k >> t[0] >> t[1] >> t[2];
        b.F.push_back(t);
      }
    }
  }
  return g;
}

// Compare a got (V, F) against a golden MeshBlock: nV/nF exact, F exact, V 1e-9.
void check_block(const std::string& cs, const char* stage, const std::vector<Vec3>& V,
                 const std::vector<Tri>& F, const MeshBlock& want) {
  const int gotNV = static_cast<int>(V.size());
  const int gotNF = static_cast<int>(F.size());
  CHECK(gotNV == want.nV);
  CHECK(gotNF == want.nF);
  if (gotNV != want.nV || gotNF != want.nF) {
    std::fprintf(stderr, "  %s %s nV got=%d want=%d  nF got=%d want=%d\n",
                 cs.c_str(), stage, gotNV, want.nV, gotNF, want.nF);
    return;  // shapes differ -> per-element compare would be noise
  }
  bool vfirst = true;
  for (int k = 0; k < gotNV; ++k) {
    const Vec3& got = V[static_cast<std::size_t>(k)];
    const Vec3& w = want.V[static_cast<std::size_t>(k)];
    const bool near = std::fabs(got.x - w.x) <= 1e-9 &&
                      std::fabs(got.y - w.y) <= 1e-9 &&
                      std::fabs(got.z - w.z) <= 1e-9;
    CHECK(near);
    if (!near && vfirst) {
      std::fprintf(stderr,
                   "  %s %s V[%d] got=(%.17g,%.17g,%.17g) want=(%.17g,%.17g,%.17g)\n",
                   cs.c_str(), stage, k, got.x, got.y, got.z, w.x, w.y, w.z);
      vfirst = false;
    }
  }
  bool ffirst = true;
  for (int k = 0; k < gotNF; ++k) {
    const Tri& f = F[static_cast<std::size_t>(k)];
    const Tri& w = want.F[static_cast<std::size_t>(k)];
    const bool ok = (f[0] == w[0]) && (f[1] == w[1]) && (f[2] == w[2]);
    CHECK(ok);
    if (!ok && ffirst) {
      std::fprintf(stderr, "  %s %s F[%d] got=(%d,%d,%d) want=(%d,%d,%d)\n",
                   cs.c_str(), stage, k, f[0], f[1], f[2], w[0], w[1], w[2]);
      ffirst = false;
    }
  }
}

void check_loops(const std::string& cs, const char* which,
                 const std::vector<Vec3>& V, const std::vector<Tri>& F,
                 const LoopSet& g) {
  BoundaryLoopsResult r = boundary_loops(V, F);
  CHECK(static_cast<int>(r.loops.size()) == g.nloops);
  CHECK(static_cast<int>(r.nonmanifold.size()) == g.nnm);
  if (static_cast<int>(r.loops.size()) != g.nloops)
    std::fprintf(stderr, "  %s %s nloops got=%d want=%d\n", cs.c_str(), which,
                 static_cast<int>(r.loops.size()), g.nloops);
  const int nL = std::min(static_cast<int>(r.loops.size()), g.nloops);
  for (int i = 0; i < nL; ++i) {
    const BoundaryLoop& got = r.loops[static_cast<std::size_t>(i)];
    const GoldenLoop& w = g.loops[static_cast<std::size_t>(i)];
    const bool fields = got.n_edges == w.n_edges && got.n_verts == w.n_verts &&
                        (got.closed ? 1 : 0) == w.closed;
    CHECK(fields);
    if (!fields) {
      std::fprintf(stderr,
                   "  %s loop[%d] got(ne=%d nv=%d cl=%d) want(ne=%d nv=%d cl=%d)\n",
                   cs.c_str(), i, got.n_edges, got.n_verts, got.closed ? 1 : 0,
                   w.n_edges, w.n_verts, w.closed);
    }
    const bool vids_ok = got.vids == w.vids;
    CHECK(vids_ok);
    const bool cen = std::fabs(got.centroid.x - w.centroid.x) <= 1e-9 &&
                     std::fabs(got.centroid.y - w.centroid.y) <= 1e-9 &&
                     std::fabs(got.centroid.z - w.centroid.z) <= 1e-9;
    CHECK(cen);
  }
  const int nNM = std::min(static_cast<int>(r.nonmanifold.size()), g.nnm);
  for (int i = 0; i < nNM; ++i) {
    const NonmanifoldEdge& got = r.nonmanifold[static_cast<std::size_t>(i)];
    const NonmanifoldEdge& w = g.nm[static_cast<std::size_t>(i)];
    const bool ok = got.u == w.u && got.v == w.v && got.count == w.count;
    CHECK(ok);
    if (!ok)
      std::fprintf(stderr, "  %s nm[%d] got(%d,%d,%d) want(%d,%d,%d)\n", cs.c_str(),
                   i, got.u, got.v, got.count, w.u, w.v, w.count);
  }
}

// --- synthetic meshes (must match dump_weld.py byte-for-byte) ----------------

void make_annulus(std::vector<Vec3>& V, std::vector<Tri>& F) {
  V.clear();
  F.clear();
  for (int i = 0; i < 6; ++i) {
    const double a = 2.0 * std::numbers::pi * i / 6.0;
    V.push_back(Vec3{std::cos(a), std::sin(a), 0.0});
  }
  for (int i = 0; i < 6; ++i) {
    const double a = 2.0 * std::numbers::pi * i / 6.0;
    V.push_back(Vec3{2.0 * std::cos(a), 2.0 * std::sin(a), 0.3});
  }
  for (int i = 0; i < 6; ++i) {
    const int j = (i + 1) % 6;
    F.push_back(Tri{i, 6 + i, 6 + j});
    F.push_back(Tri{i, 6 + j, j});
  }
}

void make_pinched(std::vector<Vec3>& V, std::vector<Tri>& F) {
  V = {Vec3{0, 0, 0},   Vec3{1, 0, 0},  Vec3{1, 1, 0},
       Vec3{0, 1, 0},   Vec3{-1, 0, 0}, Vec3{-1, -1, 0},
       Vec3{0, -1, 0},  Vec3{2, 0.5, 0.2}, Vec3{-2, -0.5, 0.2}};
  F = {Tri{1, 2, 7}, Tri{2, 3, 7}, Tri{0, 1, 7}, Tri{3, 0, 7},
       Tri{4, 5, 8}, Tri{5, 6, 8}, Tri{0, 4, 8}, Tri{6, 0, 8}};
}

void make_flap(std::vector<Vec3>& V, std::vector<Tri>& F) {
  V = {Vec3{0, 0, 0}, Vec3{2, 0, 0}, Vec3{1, 1, 0},  Vec3{1, -1, 0},
       Vec3{1, 0, 1}, Vec3{-1, 0, 0}, Vec3{3, 0, 0}};
  F = {Tri{0, 1, 2}, Tri{1, 0, 3}, Tri{0, 1, 4}, Tri{5, 0, 2},
       Tri{0, 5, 3}, Tri{1, 6, 2}, Tri{6, 1, 3}};
}

void get_case(const std::string& cs, std::vector<Vec3>& V, std::vector<Tri>& F) {
  if (cs == "annulus") {
    make_annulus(V, F);
  } else if (cs == "pinched") {
    make_pinched(V, F);
  } else if (cs == "flap") {
    make_flap(V, F);
  } else {  // real build_surface mesh
    const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/" + cs + ".xyzr";
    Surface s = build_surface(xyzr, 1.4, 0.5, /*fuse=*/false);
    V = s.V;
    F = s.F;
  }
}

void check_case(const std::string& cs) {
  const std::string refpath =
      std::string(MESHMS_REF_DIR) + "/" + cs + ".weld.txt";
  Golden g = read_golden(refpath);

  std::vector<Vec3> V0;
  std::vector<Tri> F0;
  get_case(cs, V0, F0);

  // raw mesh must match the golden RAW (for real meshes this re-checks the C++
  // build_surface; for synthetic ones it pins the reconstruction).
  check_block(cs, "RAW", V0, F0, g.raw);

  // boundary_loops on the raw mesh -- EXACT (loops + nonmanifold).
  check_loops(cs, "BL", V0, F0, g.bl);

  // weld -> flaps : both deterministic -> EXACT mesh comparison.
  WeldResult w = weld(V0, F0);
  check_block(cs, "WELD", w.V, w.F, g.weld_);

  FlapResult fl = remove_nonmanifold_flaps(w.V, w.F);
  check_block(cs, "FLAP", fl.V, fl.F, g.flap);

  // fill : nV/nF EXACT, and the carried-through original-face prefix (the first
  // FILL_nin faces) is preserved verbatim -> compared EXACT to the FLAP mesh /
  // golden.  The APPENDED fan triangles' apex is a CPython set-iteration artifact
  // (a different-but-valid fan of the same hole), so instead of comparing their
  // bytes we check the routine's actual guarantee via boundary_loops on the final
  // mesh (FBL): the same small holes are closed (here: every hole -> watertight).
  FillResult fi = fill_small_holes(fl.V, fl.F);
  CHECK(static_cast<int>(fi.V.size()) == g.fill.nV);
  CHECK(static_cast<int>(fi.F.size()) == g.fill.nF);
  // carried-through prefix: first FILL_nin faces == the FLAP faces, verbatim.
  const int prefix = std::min(g.fill_nin, static_cast<int>(fi.F.size()));
  bool pfirst = true;
  for (int k = 0; k < prefix; ++k) {
    const Tri& f = fi.F[static_cast<std::size_t>(k)];
    const Tri& wf = g.fill.F[static_cast<std::size_t>(k)];
    const bool ok = (f[0] == wf[0]) && (f[1] == wf[1]) && (f[2] == wf[2]);
    CHECK(ok);
    if (!ok && pfirst) {
      std::fprintf(stderr, "  %s FILL prefix F[%d] got=(%d,%d,%d) want=(%d,%d,%d)\n",
                   cs.c_str(), k, f[0], f[1], f[2], wf[0], wf[1], wf[2]);
      pfirst = false;
    }
  }
  // appended fan triangles must be non-degenerate (real fill, not junk).
  for (std::size_t k = static_cast<std::size_t>(prefix); k < fi.F.size(); ++k) {
    const Tri& f = fi.F[k];
    CHECK(f[0] != f[1] && f[1] != f[2] && f[0] != f[2]);
  }
  // fill guarantee: the final boundary_loops matches the golden (FBL) -- same
  // small holes closed.
  check_loops(cs, "FBL", fi.V, fi.F, g.fbl);

  std::printf(
      "  %-9s RAW nV=%d nF=%d  WELD nF=%d  FLAP nF=%d  FILL nF=%d  "
      "rawloops=%d nm=%d  finalloops=%d\n",
      cs.c_str(), static_cast<int>(V0.size()), static_cast<int>(F0.size()),
      static_cast<int>(w.F.size()), static_cast<int>(fl.F.size()),
      static_cast<int>(fi.F.size()), g.bl.nloops, g.bl.nnm, g.fbl.nloops);
}

}  // namespace

int main() {
  check_case("tetra");
  check_case("ArgArg");
  check_case("annulus");
  check_case("pinched");
  check_case("flap");
  TEST_MAIN_RETURN();
}
