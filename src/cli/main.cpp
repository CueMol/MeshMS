// Command-line entry point --- faithful port of meshms/__main__.py.
//
//   meshms_cli INPUT.xyzr -o OUTPUT.ply [--probe 1.4] [--mesh-size 0.5]
//               [--no-fuse] [--fuse-cusps] [--weld-tol 1e-4]
//               [--vertex-normals] [-v|--verbose]
//
// Default: ID-based boundary fusion connects the C1-smooth (non-cusp) surface by
// topological tag; cusp/singular boundaries stay UNFUSED (sharp shading).
// --fuse-cusps additionally welds+fills those for a closed manifold.
//
// NOTE: the Python --jitter / --no-jitter symmetry-breaking fallback is out of
// scope here (build_surface uses faithful coords); --no-jitter is accepted as a
// no-op and --jitter is rejected so the difference is explicit.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "meshms/geom.hpp"
#include "meshms/mesh_check.hpp"
#include "meshms/pipeline.hpp"
#include "meshms/ply.hpp"
#include "meshms/vec3.hpp"
#include "meshms/weld.hpp"

using namespace meshms;

// Area-weighted per-vertex normals (port of ply.vertex_normals_from_faces).
static std::vector<Vec3> vertex_normals_from_faces(const std::vector<Vec3>& V,
                                                   const std::vector<Tri>& F) {
  std::vector<Vec3> Nv(V.size(), Vec3{});
  for (const Tri& f : F) {
    const Vec3& a = V[static_cast<std::size_t>(f[0])];
    const Vec3& b = V[static_cast<std::size_t>(f[1])];
    const Vec3& c = V[static_cast<std::size_t>(f[2])];
    Vec3 fn = cross(b - a, c - a);  // magnitude == 2 * triangle area
    Nv[static_cast<std::size_t>(f[0])] = Nv[static_cast<std::size_t>(f[0])] + fn;
    Nv[static_cast<std::size_t>(f[1])] = Nv[static_cast<std::size_t>(f[1])] + fn;
    Nv[static_cast<std::size_t>(f[2])] = Nv[static_cast<std::size_t>(f[2])] + fn;
  }
  for (Vec3& nv : Nv) {
    double len = std::sqrt(pysq(nv.x) + pysq(nv.y) + pysq(nv.z));
    if (len == 0.0) len = 1.0;
    nv = nv / len;
  }
  return Nv;
}

static void usage() {
  std::fprintf(stderr,
               "usage: meshms_cli INPUT.xyzr -o OUTPUT.ply [--probe 1.4] "
               "[--mesh-size 0.5]\n"
               "                   [--no-fuse] [--fuse-cusps] [--weld-tol 1e-4] "
               "[--vertex-normals] [-v|--verbose]\n");
}

// Nearest-3 atom ids (1-based) to a point, by Euclidean distance (mirrors the
// np.argsort(norm(atoms - ctr))[:3] used in _log_non_watertight).
static std::string nearest_atoms(const Geom& g, const Vec3& ctr) {
  std::vector<std::pair<double, int>> d;
  d.reserve(static_cast<std::size_t>(g.M));
  for (int a = 1; a <= g.M; ++a) d.emplace_back(norm(g.centers[static_cast<std::size_t>(a)] - ctr), a);
  int k = std::min<int>(3, static_cast<int>(d.size()));
  std::partial_sort(d.begin(), d.begin() + k, d.end());
  std::string out;
  for (int i = 0; i < k; ++i) {
    if (i) out += ", ";
    out += "a" + std::to_string(d[static_cast<std::size_t>(i)].second);
  }
  return out;
}

int main(int argc, char** argv) {
  std::string input, output;
  double probe = 1.4, mesh_size = 0.5, weld_tol = 1e-4;
  bool no_fuse = false, fuse_cusps = false, vertex_normals = false, verbose = false;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto need = [&](const char* name) -> std::string {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "error: %s needs a value\n", name);
        std::exit(2);
      }
      return argv[++i];
    };
    if (a == "-o" || a == "--output") output = need("--output");
    else if (a == "--probe") probe = std::stod(need("--probe"));
    else if (a == "--mesh-size") mesh_size = std::stod(need("--mesh-size"));
    else if (a == "--weld-tol") weld_tol = std::stod(need("--weld-tol"));
    else if (a == "--no-fuse") no_fuse = true;
    else if (a == "--fuse-cusps") fuse_cusps = true;
    else if (a == "--vertex-normals") vertex_normals = true;
    else if (a == "-v" || a == "--verbose") verbose = true;
    else if (a == "--no-jitter") { /* C++ already uses faithful coords */ }
    else if (a == "-h" || a == "--help") { usage(); return 0; }
    else if (!a.empty() && a[0] == '-') {
      std::fprintf(stderr, "error: unknown option %s\n", a.c_str());
      usage();
      return 2;
    } else if (input.empty()) {
      input = a;
    } else {
      std::fprintf(stderr, "error: unexpected argument %s\n", a.c_str());
      usage();
      return 2;
    }
  }
  if (input.empty() || output.empty()) {
    usage();
    return 2;
  }

  Surface s = build_surface(input, probe, mesh_size, /*fuse=*/!no_fuse);
  std::vector<Vec3> V = std::move(s.V);
  std::vector<Tri> F = std::move(s.F);

  if (fuse_cusps) {
    WeldResult w = weld(V, F, weld_tol);
    V = std::move(w.V);
    F = std::move(w.F);
    FillResult fh = fill_small_holes(V, F);
    V = std::move(fh.V);
    F = std::move(fh.F);
  }

  // drop spurious doubled-flap triangles (a no-op on a clean mesh).
  FlapResult fl = remove_nonmanifold_flaps(V, F);
  V = std::move(fl.V);
  F = std::move(fl.F);

  std::vector<Vec3> Nv;
  const std::vector<Vec3>* nptr = nullptr;
  if (vertex_normals) {
    Nv = vertex_normals_from_faces(V, F);
    nptr = &Nv;
  }
  write_ply(output, V, F, nptr, vertex_normals);

  ManifoldReport rep = manifold_report(V, F);
  std::printf("wrote %s\n", output.c_str());
  std::printf("  vertices : %d\n", static_cast<int>(V.size()));
  std::printf("  faces    : %d\n", static_cast<int>(F.size()));
  std::printf("  area     : %.3f\n", rep.area);
  std::printf("  watertight: %s  (boundary edges: %d)\n",
              rep.watertight ? "True" : "False", rep.boundary_edges);

  if (!rep.watertight && verbose) {
    BoundaryLoopsResult diag = boundary_loops(V, F);
    Geom g = read_xyzr(input);
    std::printf("  non-watertight detail: %d open boundary loop(s), %d non-manifold edge(s)\n",
                static_cast<int>(diag.loops.size()), static_cast<int>(diag.nonmanifold.size()));
    const int max_loops = 20;
    for (int li = 0; li < static_cast<int>(diag.loops.size()) && li < max_loops; ++li) {
      const BoundaryLoop& L = diag.loops[static_cast<std::size_t>(li)];
      std::printf("    loop#%d: %3d edges (%s), center=[%.2f, %.2f, %.2f], near=[%s]\n",
                  li, L.n_edges, L.closed ? "closed" : "open chain",
                  L.centroid.x, L.centroid.y, L.centroid.z,
                  nearest_atoms(g, L.centroid).c_str());
    }
    if (static_cast<int>(diag.loops.size()) > max_loops)
      std::printf("    ... and %d more loop(s)\n", static_cast<int>(diag.loops.size()) - max_loops);
    if (!diag.nonmanifold.empty()) {
      std::printf("  non-manifold edges (shared by 3+ faces): %d\n",
                  static_cast<int>(diag.nonmanifold.size()));
      for (int ni = 0; ni < static_cast<int>(diag.nonmanifold.size()) && ni < 8; ++ni)
        std::printf("    edge (%d, %d) shared by %d faces\n", diag.nonmanifold[static_cast<std::size_t>(ni)].u,
                    diag.nonmanifold[static_cast<std::size_t>(ni)].v, diag.nonmanifold[static_cast<std::size_t>(ni)].count);
      if (static_cast<int>(diag.nonmanifold.size()) > 8)
        std::printf("    ... and %d more\n", static_cast<int>(diag.nonmanifold.size()) - 8);
    }
  }
  return 0;
}
