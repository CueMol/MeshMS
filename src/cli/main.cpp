// Command-line front-end for libMeshMS.
//
//   meshms_cli INPUT.xyzr -o OUT [--probe 1.4] [--mesh-size 0.5]
//               [--no-fuse] [--fuse-cusps] [--weld-tol 1e-4]
//               [--ply] [--vertex-normals] [--no-jitter] [-v|--verbose]
//
// Output format: MSMS .vert/.face by DEFAULT -- writes <OUT>.vert and <OUT>.face
// (any .vert/.face/.ply suffix on OUT is stripped to form the base name). --ply
// (or an OUT ending in .ply) writes a single ASCII PLY instead.
//
// This is the reference example of consuming libMeshMS the way cuemol2/3 does:
// EVERY surface/topology operation goes through the public facade
// (meshms/meshms.hpp) and nothing else. The CLI itself only parses arguments,
// reads the xyzr file (io_xyzr), writes the output (io_ply / io_msms), and
// reports (report) -- it never reaches into the library's internal C++20 headers.
//
// Default: ID-based boundary fusion connects the C1-smooth (non-cusp) surface by
// topological tag; cusp/singular boundaries stay UNFUSED (sharp shading).
// --fuse-cusps additionally welds+fills those for a closed manifold.
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "meshms/meshms.hpp"

#include "io_msms.hpp"
#include "io_ply.hpp"
#include "io_xyzr.hpp"
#include "report.hpp"

int main(int argc, char** argv) {
  using meshms_cli::Atom;
  std::string input, output;
  double probe = 1.4, mesh_size = 0.5, weld_tol = 1e-4;
  bool no_fuse = false, fuse_cusps = false, vertex_normals = false, verbose = false;
  bool force_ply = false, force_msms = false, no_jitter = false;

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
    else if (a == "--ply") force_ply = true;
    else if (a == "--msms") force_msms = true;
    else if (a == "--vertex-normals") vertex_normals = true;
    else if (a == "-v" || a == "--verbose") verbose = true;
    // Skip the auto symmetry-jitter NaN-retry fallback and mesh the faithful
    // coordinates exactly (Jitter::None); the default Jitter::Auto only perturbs
    // a degenerate (NaN-poisoned) molecule, so this is a no-op on clean inputs.
    else if (a == "--no-jitter") no_jitter = true;
    else if (a == "-h" || a == "--help") { meshms_cli::usage(); return 0; }
    else if (!a.empty() && a[0] == '-') {
      std::fprintf(stderr, "error: unknown option %s\n", a.c_str());
      meshms_cli::usage();
      return 2;
    } else if (input.empty()) {
      input = a;
    } else {
      std::fprintf(stderr, "error: unexpected argument %s\n", a.c_str());
      meshms_cli::usage();
      return 2;
    }
  }
  if (input.empty() || output.empty()) {
    meshms_cli::usage();
    return 2;
  }

  // --- everything below goes through the public libMeshMS facade only --------
  std::vector<Atom> atoms = meshms_cli::read_xyzr_file(input);

  const meshms::Jitter jitter =
      no_jitter ? meshms::Jitter::None : meshms::Jitter::Auto;
  meshms::MeshResult m = meshms::build_surface_from_array(
      atoms, probe, mesh_size, /*fuse=*/!no_fuse, jitter);
  if (fuse_cusps) m = meshms::close_cusps(m, weld_tol);  // weld + fill small holes
  m = meshms::remove_flaps(m);  // drop doubled flaps (no-op on a clean mesh)

  // Format: explicit flag wins; otherwise MSMS by default, with a .ply suffix on
  // OUT auto-selecting PLY for convenience.
  bool ply = force_ply || (!force_msms && meshms_cli::ends_with(output, ".ply"));
  if (ply) {
    meshms_cli::write_ply(output, m, vertex_normals);
    std::printf("wrote %s\n", output.c_str());
  } else {
    std::string base = output;  // strip a known suffix to form the <base> name
    for (const char* e : {".vert", ".face", ".ply"}) {
      if (meshms_cli::ends_with(base, e)) {
        base.resize(base.size() - std::string(e).size());
        break;
      }
    }
    meshms_cli::write_msms(base, m, static_cast<int>(atoms.size()), mesh_size, probe);
    std::printf("wrote %s.vert + %s.face\n", base.c_str(), base.c_str());
  }

  const meshms::MeshReport rep = meshms::analyze_mesh(m);
  std::printf("  vertices : %d\n", static_cast<int>(m.verts.size()));
  std::printf("  faces    : %d\n", static_cast<int>(m.faces.size()));
  std::printf("  area     : %.3f\n", rep.area);
  std::printf("  watertight: %s  (boundary edges: %d)\n",
              rep.watertight ? "True" : "False", static_cast<int>(rep.boundary_edges));
  // A relaxed-FP build has no bit-exact reference to diff against, so surface
  // this unconditionally: anything but 0 means the mesh is unusable, not merely
  // different from a strict build.
  if (rep.nonfinite_vertices != 0) {
    std::printf("  WARNING  : %d non-finite vertex/vertices (mesh is unusable)\n",
                static_cast<int>(rep.nonfinite_vertices));
  }
  if (verbose) std::printf("  build    : %s\n", meshms::build_info());

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
                  meshms_cli::nearest_atoms(atoms, L.centroid).c_str());
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
