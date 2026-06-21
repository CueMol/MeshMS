// Verify the post-processing + diagnostics public API (capi.hpp) that the CLI and
// cuemol2/3 build on. This test deliberately includes ONLY the public facade
// (plus a local xyzr loader), exercising the library exactly as an external
// consumer would -- no internal MeshMS headers.
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "meshms/capi.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

// Minimal xyzr loader (consumer side): skip blank/'#' rows, take the first 4
// columns. Mirrors read_xyzr so the array feeds the same surface.
std::vector<std::array<double, 4>> load_array(const std::string& path) {
  std::ifstream in(path);
  std::vector<std::array<double, 4>> out;
  std::string line;
  while (std::getline(in, line)) {
    std::size_t i = 0;
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    if (i == line.size() || line[i] == '#') continue;
    std::istringstream ls(line);
    std::vector<double> c;
    double v;
    while (ls >> v) c.push_back(v);
    if (c.size() >= 4) out.push_back({c[0], c[1], c[2], c[3]});
  }
  return out;
}

double vlen(const std::array<double, 3>& n) {
  return std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
}

}  // namespace

int main() {
  // version() is a usable, non-empty string.
  CHECK(version() != nullptr);
  CHECK(std::string(version()).size() > 0);

  const std::string xyzr = std::string(MESHMS_DATA_DIR) + "/ArgArg.xyzr";
  std::vector<std::array<double, 4>> arr = load_array(xyzr);
  CHECK(!arr.empty());

  // --- fused build is a clean, watertight 2-manifold -------------------------
  MeshResult fused = build_surface_from_array(arr, 1.4, 0.5, /*fuse=*/true);
  MeshReport rf = analyze_mesh(fused);
  CHECK(rf.n_vertices == fused.verts.size());
  CHECK(rf.n_faces == fused.faces.size());
  CHECK(rf.watertight);
  CHECK(rf.boundary_edges == 0);
  CHECK(rf.nonmanifold_edges == 0);
  CHECK(rf.area > 0.0);
  std::printf("  fused: V=%zu F=%zu area=%.3f watertight=%d\n", fused.verts.size(),
              fused.faces.size(), rf.area, rf.watertight ? 1 : 0);

  // --- unfused build has open patch seams; diagnostics localize them ---------
  MeshResult open = build_surface_from_array(arr, 1.4, 0.5, /*fuse=*/false);
  MeshReport ro = analyze_mesh(open);
  CHECK(!ro.watertight);
  CHECK(ro.boundary_edges > 0);

  BoundaryDiagnostics diag = boundary_diagnostics(open);
  CHECK(!diag.loops.empty());
  // loops are sorted largest-first by edge count.
  bool sorted = true;
  for (std::size_t i = 1; i < diag.loops.size(); ++i)
    if (diag.loops[i].n_edges > diag.loops[i - 1].n_edges) sorted = false;
  CHECK(sorted);
  std::printf("  open : boundary_edges=%u loops=%zu (largest=%u edges)\n",
              ro.boundary_edges, diag.loops.size(),
              diag.loops.empty() ? 0u : diag.loops.front().n_edges);

  // --- remove_flaps is a no-op on the clean fused mesh, and keeps atom_id -----
  MeshResult flapped = remove_flaps(fused);
  CHECK(flapped.verts.size() == fused.verts.size());
  CHECK(flapped.faces.size() == fused.faces.size());
  CHECK(flapped.atom_id.size() == flapped.verts.size());
  // its vnormals are the area-weighted ones vertex_normals() recomputes.
  std::vector<std::array<double, 3>> nv = vertex_normals(flapped.verts, flapped.faces);
  CHECK(nv.size() == flapped.verts.size());
  CHECK(nv == flapped.vnormals);
  // every vertex of a closed mesh has faces -> unit normals.
  bool unit = true;
  for (const std::array<double, 3>& n : nv)
    if (std::fabs(vlen(n) - 1.0) > 1e-9) unit = false;
  CHECK(unit);

  // --- close_cusps stays watertight and drops atom_id (documented) -----------
  MeshResult closed = close_cusps(fused);
  MeshReport rc = analyze_mesh(closed);
  CHECK(rc.watertight);
  CHECK(closed.atom_id.empty());
  CHECK(closed.vnormals.size() == closed.verts.size());

  TEST_MAIN_RETURN();
}
