// Command-line front-end for libMeshMS.
//
//   meshms_cli INPUT.xyzr -o OUTPUT.ply [--probe 1.4] [--mesh-size 0.5]
//               [--no-fuse] [--fuse-cusps] [--weld-tol 1e-4]
//               [--vertex-normals] [-v|--verbose]
//
// This is the reference example of consuming libMeshMS the way cuemol2/3 does:
// EVERY surface/topology operation goes through the public facade
// (meshms/capi.hpp) and nothing else. The CLI itself only parses arguments,
// reads the xyzr file, writes the PLY, and reports -- it never reaches into the
// library's internal C++20 headers.
//
// Default: ID-based boundary fusion connects the C1-smooth (non-cusp) surface by
// topological tag; cusp/singular boundaries stay UNFUSED (sharp shading).
// --fuse-cusps additionally welds+fills those for a closed manifold.
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "meshms/capi.hpp"

namespace {

using Atom = std::array<double, 4>;  // x, y, z, radius
using meshms::MeshResult;

// Parse an xyzr file like numpy.loadtxt (mirrors the library's read_xyzr): skip
// blank and '#'-comment lines, split each row on whitespace, keep the first 4
// columns. Atoms are returned in file order, so a MeshResult.atom_id value i
// refers to the (i-1)-th atom here.
std::vector<Atom> read_xyzr_file(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    std::fprintf(stderr, "error: cannot open input file: %s\n", path.c_str());
    std::exit(2);
  }
  std::vector<Atom> atoms;
  std::string line;
  int row = 0;
  while (std::getline(in, line)) {
    std::size_t i = 0;
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    if (i == line.size() || line[i] == '#') continue;  // blank / comment line
    std::istringstream ls(line);
    std::vector<double> cols;
    double v;
    while (ls >> v) cols.push_back(v);
    if (cols.empty()) continue;
    ++row;
    if (cols.size() < 4) {
      std::fprintf(stderr, "error: xyzr row %d has <4 columns\n", row);
      std::exit(2);
    }
    atoms.push_back({cols[0], cols[1], cols[2], cols[3]});
  }
  if (atoms.empty()) {
    std::fprintf(stderr, "error: no data rows in %s\n", path.c_str());
    std::exit(2);
  }
  return atoms;
}

// Write the mesh as ASCII PLY: positions (%.6f), optional per-vertex normals,
// then faces ("3 i j k"). Byte-for-byte the same format the library uses.
void write_ply(const std::string& path, const MeshResult& m, bool with_normals) {
  std::ofstream fh(path);
  if (!fh) {
    std::fprintf(stderr, "error: cannot open PLY for writing: %s\n", path.c_str());
    std::exit(2);
  }
  const std::size_t nv = m.verts.size();
  const std::size_t nf = m.faces.size();
  fh << "ply\n" << "format ascii 1.0\n" << "element vertex " << nv << "\n"
     << "property float x\n" << "property float y\n" << "property float z\n";
  if (with_normals)
    fh << "property float nx\n" << "property float ny\n" << "property float nz\n";
  fh << "element face " << nf << "\n"
     << "property list uchar int vertex_indices\n" << "end_header\n";

  char buf[256];
  for (std::size_t i = 0; i < nv; ++i) {
    const std::array<double, 3>& v = m.verts[i];
    if (with_normals) {
      const std::array<double, 3>& n = m.vnormals[i];
      std::snprintf(buf, sizeof(buf), "%.6f %.6f %.6f %.6f %.6f %.6f", v[0], v[1],
                    v[2], n[0], n[1], n[2]);
    } else {
      std::snprintf(buf, sizeof(buf), "%.6f %.6f %.6f", v[0], v[1], v[2]);
    }
    fh << buf << "\n";
  }
  for (const std::array<std::uint32_t, 3>& f : m.faces) {
    std::snprintf(buf, sizeof(buf), "3 %d %d %d", static_cast<int>(f[0]),
                  static_cast<int>(f[1]), static_cast<int>(f[2]));
    fh << buf << "\n";
  }
}

// Nearest-3 atom labels (1-based "a<k>") to a point, by Euclidean distance --
// the same selection the diagnostic used (argsort of |atom - ctr|, first 3).
std::string nearest_atoms(const std::vector<Atom>& atoms,
                          const std::array<double, 3>& ctr) {
  std::vector<std::pair<double, int>> d;
  d.reserve(atoms.size());
  for (std::size_t a = 0; a < atoms.size(); ++a) {
    const double dx = atoms[a][0] - ctr[0];
    const double dy = atoms[a][1] - ctr[1];
    const double dz = atoms[a][2] - ctr[2];
    d.emplace_back(std::sqrt(dx * dx + dy * dy + dz * dz), static_cast<int>(a) + 1);
  }
  const int k = std::min<int>(3, static_cast<int>(d.size()));
  std::partial_sort(d.begin(), d.begin() + k, d.end());
  std::string out;
  for (int i = 0; i < k; ++i) {
    if (i) out += ", ";
    out += "a" + std::to_string(d[static_cast<std::size_t>(i)].second);
  }
  return out;
}

void usage() {
  std::fprintf(stderr,
               "usage: meshms_cli INPUT.xyzr -o OUTPUT.ply [--probe 1.4] "
               "[--mesh-size 0.5]\n"
               "                   [--no-fuse] [--fuse-cusps] [--weld-tol 1e-4] "
               "[--vertex-normals] [-v|--verbose]\n");
}

}  // namespace

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
    else if (a == "--no-jitter") { /* libMeshMS always uses faithful coords */ }
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

  // --- everything below goes through the public libMeshMS facade only --------
  std::vector<Atom> atoms = read_xyzr_file(input);

  MeshResult m =
      meshms::build_surface_from_array(atoms, probe, mesh_size, /*fuse=*/!no_fuse);
  if (fuse_cusps) m = meshms::close_cusps(m, weld_tol);  // weld + fill small holes
  m = meshms::remove_flaps(m);  // drop doubled flaps (no-op on a clean mesh)

  write_ply(output, m, vertex_normals);

  const meshms::MeshReport rep = meshms::analyze_mesh(m);
  std::printf("wrote %s\n", output.c_str());
  std::printf("  vertices : %d\n", static_cast<int>(m.verts.size()));
  std::printf("  faces    : %d\n", static_cast<int>(m.faces.size()));
  std::printf("  area     : %.3f\n", rep.area);
  std::printf("  watertight: %s  (boundary edges: %d)\n",
              rep.watertight ? "True" : "False", static_cast<int>(rep.boundary_edges));

  if (!rep.watertight && verbose) {
    const meshms::BoundaryDiagnostics diag = meshms::boundary_diagnostics(m);
    std::printf("  non-watertight detail: %d open boundary loop(s), %d non-manifold edge(s)\n",
                static_cast<int>(diag.loops.size()),
                static_cast<int>(diag.nonmanifold.size()));
    const int max_loops = 20;
    for (int li = 0; li < static_cast<int>(diag.loops.size()) && li < max_loops; ++li) {
      const meshms::BoundaryLoopInfo& L = diag.loops[static_cast<std::size_t>(li)];
      std::printf("    loop#%d: %3d edges (%s), center=[%.2f, %.2f, %.2f], near=[%s]\n",
                  li, static_cast<int>(L.n_edges), L.closed ? "closed" : "open chain",
                  L.centroid[0], L.centroid[1], L.centroid[2],
                  nearest_atoms(atoms, L.centroid).c_str());
    }
    if (static_cast<int>(diag.loops.size()) > max_loops)
      std::printf("    ... and %d more loop(s)\n",
                  static_cast<int>(diag.loops.size()) - max_loops);
    if (!diag.nonmanifold.empty()) {
      std::printf("  non-manifold edges (shared by 3+ faces): %d\n",
                  static_cast<int>(diag.nonmanifold.size()));
      for (int ni = 0; ni < static_cast<int>(diag.nonmanifold.size()) && ni < 8; ++ni)
        std::printf("    edge (%d, %d) shared by %d faces\n",
                    static_cast<int>(diag.nonmanifold[static_cast<std::size_t>(ni)].u),
                    static_cast<int>(diag.nonmanifold[static_cast<std::size_t>(ni)].v),
                    static_cast<int>(diag.nonmanifold[static_cast<std::size_t>(ni)].count));
      if (static_cast<int>(diag.nonmanifold.size()) > 8)
        std::printf("    ... and %d more\n", static_cast<int>(diag.nonmanifold.size()) - 8);
    }
  }
  return 0;
}
