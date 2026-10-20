// Regression: a concave probe patch whose boundary has single-segment loops.
//
// loopconstruct_concave used to size its loop table to N_segment/2 rows, on the
// assumption that every loop has at least two segments. A segment whose two
// ends map_same merges into one point closes a loop by itself, so a probe with
// such loops wrote past the end of the table and the build crashed. The input
// is the 38 atoms around the probe that hit it in human 80S ribosome 4V6X
// (CueMol vdW radii): 9 segments in loops of 3, 1, 3, 1 and 1 segments.
#include <array>
#include <string>
#include <vector>

#include "meshms/geom.hpp"
#include "meshms/meshms.hpp"
#include "test_util.hpp"

using namespace meshms;

int main() {
  Geom g = read_xyzr(std::string(MESHMS_DATA_DIR) + "/loop1seg.xyzr");
  std::vector<std::array<double, 4>> xyzr;
  for (int a = 1; a <= g.M; ++a) {
    const Vec3& c = g.centers[static_cast<std::size_t>(a)];
    xyzr.push_back({c.x, c.y, c.z, g.R[static_cast<std::size_t>(a)]});
  }

  MeshResult m = build_mesh_from_cache(compute_rs_from_array(xyzr, 1.4), 0.834, true);
  CHECK(!m.faces.empty());
  CHECK(analyze_mesh(m).nonfinite_vertices == 0);
  TEST_MAIN_RETURN();
}
