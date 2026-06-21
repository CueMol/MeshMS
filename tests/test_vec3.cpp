// Smoke test for the vec3 numerical core (faithful to the mathutil module).
#include "meshms/vec3.hpp"

#include "test_util.hpp"

using namespace meshms;

int main() {
  const double EPS = 1e-15;

  // dot / norm
  CHECK_NEAR(dot(Vec3{1, 2, 3}, Vec3{4, 5, 6}), 32.0, EPS);
  CHECK_NEAR(norm(Vec3{3, 4, 0}), 5.0, EPS);

  // cross: e_x x e_y == e_z
  Vec3 cz = cross(Vec3{1, 0, 0}, Vec3{0, 1, 0});
  CHECK_NEAR(cz.x, 0.0, EPS);
  CHECK_NEAR(cz.y, 0.0, EPS);
  CHECK_NEAR(cz.z, 1.0, EPS);

  // triple of the basis == det(I) == 1
  CHECK_NEAR(triple(Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}), 1.0, EPS);

  // unit
  Vec3 u = unit(Vec3{0, 0, 5});
  CHECK_NEAR(u.z, 1.0, EPS);

  // sign / acos_clamped
  CHECK_NEAR(sign(-3.0), -1.0, EPS);
  CHECK_NEAR(sign(0.0), 0.0, EPS);
  CHECK_NEAR(acos_clamped(2.0), 0.0, EPS);    // clamped to 1 -> acos(1) == 0
  CHECK_NEAR(acos_clamped(-2.0), std::numbers::pi, EPS);

  // arc_angle: 90 degrees in-plane (axis +z), direct = +1
  CHECK_NEAR(arc_angle(Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}, 1),
             std::numbers::pi / 2.0, 1e-14);
  // opposite winding -> 2*pi - pi/2
  CHECK_NEAR(arc_angle(Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, -1}, 1),
             TWO_PI - std::numbers::pi / 2.0, 1e-14);

  // orthogonalvectors: results are unit and mutually orthogonal with n
  auto [v1, v2] = orthogonalvectors(Vec3{0, 0, 1});
  CHECK_NEAR(norm(v1), 1.0, EPS);
  CHECK_NEAR(norm(v2), 1.0, EPS);
  CHECK_NEAR(dot(v1, Vec3{0, 0, 1}), 0.0, EPS);
  CHECK_NEAR(dot(v2, Vec3{0, 0, 1}), 0.0, EPS);
  CHECK_NEAR(dot(v1, v2), 0.0, EPS);

  // circlecenter: two equal unit spheres at +/-1 on x -> midpoint at origin
  Vec3 cc = circlecenter(Vec3{-1, 0, 0}, Vec3{1, 0, 0}, 1.0, 1.0);
  CHECK_NEAR(cc.x, 0.0, EPS);
  CHECK_NEAR(cc.y, 0.0, EPS);
  CHECK_NEAR(cc.z, 0.0, EPS);

  TEST_MAIN_RETURN();
}
