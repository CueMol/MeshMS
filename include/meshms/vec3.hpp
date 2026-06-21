#pragma once
// Fused scalar length-3 vector math --- faithful port of the mathutil module.
//
// The float evaluation order is kept identical to the numpy helpers (same term
// order in dot/cross/circlecenter) so the C++ port reproduces the Python result
// to the last ULP wherever the algorithm is the same. This is the deliberate
// "scalar fused len-3, no Eigen" choice: it dissolves
// the Python numpy-dispatch self-time at zero cost AND keeps the evaluation order
// under our control (a hard requirement of the advancing-front faithfulness trap).
#include <cmath>
#include <numbers>
#include <utility>

namespace meshms {

// 2*np.pi in IEEE-754 double (std::numbers::pi == numpy's np.pi).
inline constexpr double TWO_PI = 2.0 * std::numbers::pi;

struct Vec3 {
  double x{0.0}, y{0.0}, z{0.0};
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(const Vec3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline Vec3 operator*(double s, const Vec3& a) { return {a.x * s, a.y * s, a.z * s}; }
inline Vec3 operator/(const Vec3& a, double s) { return {a.x / s, a.y / s, a.z / s}; }

// np.dot(a, b) for length-3 vectors (x, then y, then z).
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

// norm(v) == float(np.sqrt(np.dot(v, v))).
inline double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }

// unit(v) == v / norm(v).
inline Vec3 unit(const Vec3& v) { return v / norm(v); }

// cross(a, b): same component order as mathutil.cross.
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

// triple(u, v, n) == u . (v x n) == det([u; v; n]).
inline double triple(const Vec3& u, const Vec3& v, const Vec3& n) { return dot(u, cross(v, n)); }

// np.sign: 0 for 0, +1/-1 otherwise.
inline double sign(double x) { return static_cast<double>((x > 0.0) - (x < 0.0)); }

// Python scalar `x ** 2` == C pow(x, 2.0), which is NOT bit-identical to x*x:
// they differ by up to 1 ULP on ~0.08% of inputs -- enough to flip a < / <= / >
// boundary test and so change the SAS triple-point set. Use pysq() wherever the
// Python source wrote `x ** 2` on a SCALAR, so the C++ matches the golden to the
// last ULP. (numpy ARRAY x**2 *is* x*x; and mathutil.circlecenter writes r1*r1
// literally -> those squares stay '*' in both Python and C++.)
inline double pysq(double x) { return std::pow(x, 2.0); }

// Clamp to [-1, 1] then arccos (mathutil._acos_clamped).
inline double acos_clamped(double x) {
  if (x > 1.0) x = 1.0;
  else if (x < -1.0) x = -1.0;
  return std::acos(x);
}

// arc_angle(u, v, n, direct): angle from u to v around axis n (mathutil.arc_angle).
inline double arc_angle(const Vec3& u, const Vec3& v, const Vec3& n, int direct = 1) {
  double t = sign(triple(u, v, n));
  double base = acos_clamped(dot(u, v) / (norm(u) * norm(v)));
  if (direct * t > 0.0) return base;
  return TWO_PI - base;
}

// orthogonalvectors(n) -> (v1, v2) with (v1, v2, n) right-handed
// (faithful port of mathutil.orthogonalvectors / orthogonalvectors.m).
inline std::pair<Vec3, Vec3> orthogonalvectors(const Vec3& n) {
  const double a = n.x, b = n.y, c = n.z;
  double x, y, z;
  if (a != 0.0 && b != 0.0) {
    x = -sign(a) * b; y = sign(a) * a; z = 0.0;
  } else if (b != 0.0 && c != 0.0) {
    x = 0.0; y = -sign(b) * c; z = sign(b) * b;
  } else if (c != 0.0 && a != 0.0) {
    x = sign(c) * c; y = 0.0; z = -sign(c) * a;
  } else if (a != 0.0) {
    x = 0.0; y = 1.0; z = 0.0;
  } else if (b != 0.0) {
    x = 0.0; y = 0.0; z = 1.0;
  } else {  // c != 0
    x = 1.0; y = 0.0; z = 0.0;
  }
  Vec3 v1 = unit(Vec3{x, y, z});
  Vec3 v2 = unit(Vec3{b * z - c * y, c * x - a * z, a * y - b * x});
  return {v1, v2};
}

// circlecenter(c1, c2, r1, r2): center of the intersection circle of two spheres
// (faithful port of mathutil.circlecenter / circlecenter.m). The term order
// (c2 - c1) * t / d matches the numpy elementwise (*t then /d) rounding.
inline Vec3 circlecenter(const Vec3& c1, const Vec3& c2, double r1, double r2) {
  double d = norm(c1 - c2);
  double t = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
  return c1 + (c2 - c1) * t / d;
}

}  // namespace meshms
