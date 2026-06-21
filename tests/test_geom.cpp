// Smoke test for read_xyzr / Geom (real cross-check is test_intersection later).
#include <cmath>
#include <string>

#include "meshms/geom.hpp"
#include "meshms/params.hpp"
#include "test_util.hpp"

using namespace meshms;

int main() {
  const std::string path = std::string(MESHMS_DATA_DIR) + "/ArgArg.xyzr";
  Geom g = read_xyzr(path);

  CHECK(g.M > 0);
  CHECK(g.centers.size() == static_cast<std::size_t>(g.M + 1));
  CHECK(g.R.size() == static_cast<std::size_t>(g.M + 1));

  // centers[1] is a real atom and must be finite.
  CHECK(std::isfinite(g.centers[1].x));
  CHECK(std::isfinite(g.centers[1].y));
  CHECK(std::isfinite(g.centers[1].z));
  CHECK(std::isfinite(g.R[1]));

  // Row 0 is dummy (left zero).
  CHECK(g.centers[0].x == 0.0 && g.centers[0].y == 0.0 && g.centers[0].z == 0.0);
  CHECK(g.R[0] == 0.0);

  // Para defaults match the params module.
  Para p;
  CHECK_NEAR(p.radius_probe, 1.4, 0.0);
  CHECK_NEAR(p.mesh_size, 0.5, 0.0);

  TEST_MAIN_RETURN();
}
