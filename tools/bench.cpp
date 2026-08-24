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
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

#include "meshms/concave.hpp"
#include "meshms/convex.hpp"
#include "meshms/exterior.hpp"
#include "meshms/geom.hpp"
#include "meshms/intersection.hpp"
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
  auto t7 = clk::now();
  std::vector<Vec3> V = st.V;
  std::vector<Tri> F = st.F;
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
  return s;
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

  std::printf("%-9s %6s %9s | %8s %9s %8s %7s %7s %7s %7s %6s %6s\n", "mol", "atoms",
              "total_ms", "inter", "dICir", "segpat", "ext", "convex", "concav", "torus",
              "orient", "nV");
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

    std::printf("%-9s %6d %9.2f | %8.2f %9.2f %8.2f %7.2f %7.2f %7.2f %7.2f %6.2f %6zu\n",
                mol.c_str(), geom.M, best.total, best.inter, best.dICir, best.segpat,
                best.ext, best.convex, best.concave, best.torus, best.orient, best.nV);
  }
  return 0;
}
