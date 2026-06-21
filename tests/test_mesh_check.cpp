// §5 equivalence gate: read the Python golden mesh (.golden.ply), run the C++
// manifold_report, and compare against the Python-dumped gate (.gate.txt).
//   nV, nF, boundary_edges, nonmanifold_edges : EXACT
//   watertight                                : matches
//   area, signed_volume                       : RELATIVE 1e-6 (|c-p| <= 1e-6*|p|)
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "meshms/mesh_check.hpp"
#include "meshms/ply.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

// Parsed contents of <mol>.gate.txt.
struct Gate {
  int nV{0};
  int nF{0};
  double area{0.0};
  double signed_volume{0.0};
  int watertight{0};
  int boundary_edges{0};
  int nonmanifold_edges{0};
};

Gate read_gate(const std::string& path) {
  std::ifstream in(path);
  CHECK(static_cast<bool>(in));
  Gate g;
  std::string key;
  while (in >> key) {
    if (key == "nV") in >> g.nV;
    else if (key == "nF") in >> g.nF;
    else if (key == "area") in >> g.area;
    else if (key == "signed_volume") in >> g.signed_volume;
    else if (key == "watertight") in >> g.watertight;
    else if (key == "boundary_edges") in >> g.boundary_edges;
    else if (key == "nonmanifold_edges") in >> g.nonmanifold_edges;
  }
  return g;
}

// |c - p| <= 1e-6 * |p|  (relative tolerance about the Python value p).
void check_rel(double c, double p, const char* what, const std::string& mol) {
  double d = std::fabs(c - p);
  double tol = 1e-6 * std::fabs(p);
  CHECK(d <= tol);
  if (!(d <= tol)) {
    std::fprintf(stderr, "  mol=%s %s got=%.17g want=%.17g |d|=%.3g tol=%.3g\n",
                 mol.c_str(), what, c, p, d, tol);
  }
}

void check_mol(const std::string& mol) {
  const std::string plypath = std::string(MESHMS_REF_DIR) + "/" + mol + ".golden.ply";
  const std::string gatepath = std::string(MESHMS_REF_DIR) + "/" + mol + ".gate.txt";

  PlyMesh m = read_ply(plypath);
  ManifoldReport rep = manifold_report(m.V, m.F);
  Gate g = read_gate(gatepath);

  // Exact integer counts.
  CHECK(rep.n_vertices == g.nV);
  CHECK(rep.n_faces == g.nF);
  CHECK(rep.boundary_edges == g.boundary_edges);
  CHECK(rep.nonmanifold_edges == g.nonmanifold_edges);
  CHECK(rep.watertight == (g.watertight != 0));

  if (rep.n_vertices != g.nV || rep.n_faces != g.nF ||
      rep.boundary_edges != g.boundary_edges ||
      rep.nonmanifold_edges != g.nonmanifold_edges) {
    std::fprintf(stderr,
                 "  mol=%s nV got=%d want=%d nF got=%d want=%d bnd got=%d want=%d nm got=%d want=%d\n",
                 mol.c_str(), rep.n_vertices, g.nV, rep.n_faces, g.nF,
                 rep.boundary_edges, g.boundary_edges, rep.nonmanifold_edges,
                 g.nonmanifold_edges);
  }

  // Relative-tolerance reals.
  check_rel(rep.area, g.area, "area", mol);
  check_rel(rep.signed_volume, g.signed_volume, "signed_volume", mol);
}

}  // namespace

int main() {
  check_mol("tetra");
  check_mol("ArgArg");
  TEST_MAIN_RETURN();
}
