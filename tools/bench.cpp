// Per-stage + end-to-end timing benchmark for the C++ SES pipeline.
// Times the full build_surface (run(geom, para, fuse=false)) AND the per-stage
// breakdown (so data_I_Cir + mesher speedups can be reported). Each timing is a
// best-of-`reps` to suppress scheduler/cache noise; reps defaults to 5 and can
// be overridden by env BENCH_REPS.
//
// Built as the `meshms_bench` CMake target (top-level builds; MESHMS_BUILD_TOOLS),
// or stand-alone:
//   g++ -std=c++20 -O3 -march=native -ffp-contract=off \
//       -I include tools/bench.cpp build/libMeshMS.a -ltbb -o /tmp/bench
// Run:  meshms_bench <data_dir> <mol>:<mesh_size> ... (oneTBB auto-detects cores)
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "meshms/meshms.hpp"
#include "meshms/concave.hpp"
#include "meshms/convex.hpp"
#include "meshms/exterior.hpp"
#include "meshms/geom.hpp"
#include "meshms/intersection.hpp"
#include "meshms/mesh_check.hpp"
#include "meshms/mesh_state.hpp"
#include "meshms/params.hpp"
#include "meshms/pipeline.hpp"
#include "meshms/sas.hpp"
#include "meshms/sas_patches.hpp"
#include "meshms/toroidal.hpp"

using namespace meshms;
using clk = std::chrono::steady_clock;
static double ms(clk::time_point a, clk::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

struct Stage {
  double total = 0, inter = 0, dICir = 0, segpat = 0, ext = 0;
  double convex = 0, concave = 0, torus = 0, orient = 0;
  std::size_t nV = 0;
  double area = 0.0;  // computed outside the timed region
};

// One full build_surface (fuse=false) with a per-stage breakdown.
static Stage run_once(const Geom& geom, double Rp, double d) {
  Stage s;
  auto t0 = clk::now();
  Neighbors nb = interstructure(geom, Rp);
  auto t1 = clk::now();
  auto [di, dc] = data_I_Cir(geom, nb, Rp);
  auto t2 = clk::now();
  auto [ds, dl, dp] = data_Seg_Pat(geom, nb, di, dc, Rp);
  auto t3 = clk::now();
  Ext ext = data_ext(geom, nb, di, dc, ds, dl, dp, Rp);
  auto t4 = clk::now();
  MeshState st;
  data_SESsphpat_convex(st, geom, di, dc, ds, dl, dp, &ext, Rp, d);
  auto t5 = clk::now();
  SESconcavepat(st, geom, di, ext, Rp, d, nb);
  auto t6 = clk::now();
  data_SEStorpat(st, geom, di, ds, dc, &ext, Rp, d);
  // Move, don't copy: the real pipeline (pipeline.cpp) moves V/F out of the
  // MeshState, so timing a ~10 MB deep copy here inflated the orient column by
  // roughly 6x on the large molecules.
  std::vector<Vec3> V = std::move(st.V);
  std::vector<Tri> F = std::move(st.F);
  auto t7 = clk::now();
  orient_faces(V, F, st.N);
  auto t8 = clk::now();

  s.total = ms(t0, t8);
  s.inter = ms(t0, t1);
  s.dICir = ms(t1, t2);
  s.segpat = ms(t2, t3);
  s.ext = ms(t3, t4);
  s.convex = ms(t4, t5);
  s.concave = ms(t5, t6);
  s.torus = ms(t6, t7);
  s.orient = ms(t7, t8);
  s.nV = V.size();
  // Outside every timed span: shows how far the mesh itself moved between two
  // builds, in the same table as the timings.
  s.area = mesh_area(V, F);
  return s;
}

// The facade path a host application (cuemol2) runs: RS, fused mesh, then the
// flap clean-up. Each stage is timed on its own so remove_flaps can be compared
// across builds; nF is the face count left after remove_flaps.
struct FacadeStage {
  double total = 0, rs = 0, mesh = 0, flaps = 0;
  std::size_t nF = 0;
};

static FacadeStage run_facade_once(const std::vector<std::array<double, 4>>& xyzr,
                                   double Rp, double d) {
  FacadeStage s;
  auto t0 = clk::now();
  std::shared_ptr<RSCache> rs = compute_rs_from_array(xyzr, Rp);
  auto t1 = clk::now();
  MeshResult mesh = build_mesh_from_cache(rs, d, /*fuse=*/true);
  auto t2 = clk::now();
  MeshResult clean = remove_flaps(mesh);
  auto t3 = clk::now();
  s.total = ms(t0, t3);
  s.rs = ms(t0, t1);
  s.mesh = ms(t1, t2);
  s.flaps = ms(t2, t3);
  s.nF = clean.faces.size();
  return s;
}

// BENCH_MODE=facade: best-of-`reps` per stage (each column is its own minimum,
// since the stages are independent of each other's timing noise).
static void bench_facade(const std::string& dir, int argc, char** argv, int reps,
                         double Rp) {
  std::printf("%-9s %6s %10s | %9s %10s %9s %9s\n", "mol", "atoms", "total_ms", "rs",
              "mesh", "flaps", "nF");
  for (int ai = 2; ai < argc; ++ai) {
    std::string arg = argv[ai];
    auto colon = arg.find(':');
    std::string mol = arg.substr(0, colon);
    double d = std::stod(arg.substr(colon + 1));
    Geom geom = read_xyzr(dir + "/" + mol + ".xyzr");
    std::vector<std::array<double, 4>> xyzr;
    xyzr.reserve(static_cast<std::size_t>(geom.M));
    for (int i = 1; i <= geom.M; ++i) {
      const Vec3& c = geom.centers[static_cast<std::size_t>(i)];
      xyzr.push_back({c.x, c.y, c.z, geom.R[static_cast<std::size_t>(i)]});
    }

    FacadeStage best;
    best.total = best.rs = best.mesh = best.flaps = std::numeric_limits<double>::max();
    for (int r = 0; r < reps; ++r) {
      FacadeStage s = run_facade_once(xyzr, Rp, d);
      best.total = std::min(best.total, s.total);
      best.rs = std::min(best.rs, s.rs);
      best.mesh = std::min(best.mesh, s.mesh);
      best.flaps = std::min(best.flaps, s.flaps);
      best.nF = s.nF;
    }
    std::printf("%-9s %6d %10.3f | %9.3f %10.3f %9.3f %9zu\n", mol.c_str(), geom.M,
                best.total, best.rs, best.mesh, best.flaps, best.nF);
  }
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: bench <data_dir> <mol>:<mesh> ...\n");
    return 2;
  }
  std::string dir = argv[1];
  const double Rp = 1.4;
  int reps = 5;
  if (const char* e = std::getenv("BENCH_REPS")) reps = std::atoi(e);
  if (reps < 1) reps = 1;

  // Banner names the build, so a pasted table is self-describing about which FP
  // policy produced it.
  std::printf("# %s\n", build_info());
  if (const char* m = std::getenv("BENCH_MODE")) {
    if (std::string(m) == "facade") {
      bench_facade(dir, argc, argv, reps, Rp);
      return 0;
    }
  }
  std::printf("%-9s %6s %9s | %8s %9s %8s %7s %7s %7s %7s %6s %7s %10s\n", "mol", "atoms",
              "total_ms", "inter", "dICir", "segpat", "ext", "convex", "concav", "torus",
              "orient", "nV", "area");
  for (int ai = 2; ai < argc; ++ai) {
    std::string arg = argv[ai];
    auto colon = arg.find(':');
    std::string mol = arg.substr(0, colon);
    double d = std::stod(arg.substr(colon + 1));
    std::string path = dir + "/" + mol + ".xyzr";

    Geom geom = read_xyzr(path);

    // Best-of-`reps`: take the minimum per-stage time across runs (each run is a
    // complete independent build_surface). Minimum is the least-noisy estimate.
    Stage best;
    best.total = std::numeric_limits<double>::max();
    for (int r = 0; r < reps; ++r) {
      Stage s = run_once(geom, Rp, d);
      if (s.total < best.total) best = s;  // keep the fastest complete run
    }

    std::printf("%-9s %6d %9.2f | %8.2f %9.2f %8.2f %7.2f %7.2f %7.2f %7.2f %6.2f %7zu %10.3f\n",
                mol.c_str(), geom.M, best.total, best.inter, best.dICir, best.segpat,
                best.ext, best.convex, best.concave, best.torus, best.orient, best.nV,
                best.area);
  }
  return 0;
}
