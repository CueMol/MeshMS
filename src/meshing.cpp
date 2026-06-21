// Advancing-front spherical patch mesher --- faithful port of the meshing module.
//
// BYTE-FRAGILE port: every scalar `x ** 2` -> pysq(); every float term / eval
// order reproduced with the vec3 helpers; the EXACT comparison order and tie-
// breaks, the index_np selection, the edgedist accumulation, the angle conditions,
// the m loop bounds (range 3..Nae, mod(m,Nae)+1) and the recursion are matched
// line-by-line against the Python. The Ae 1-based list (dummy Ae[0]) and the
// _ae_from_rows re-wrapping are preserved.
#include "meshms/meshing.hpp"

#include <cmath>
#include <deque>

#include "meshms/parallel.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

namespace {

using Edge = std::array<int, 2>;        // [tail, head], 1-based
using AeList = std::vector<Edge>;       // 1-based list, Ae[0] dummy

// ---------------------------------------------------------------------------
// FrontRing: the advancing front as a 1-based ring buffer. E(m), m in 1..n,
// lives at buf[(head + m - 1) & mask]. The per-triangle front rewrites
// (addnewpoint, both collapse_neighbor cases) are head/tail moves plus appends
// -- O(1) instead of rebuilding the whole front, which cost O(Nae) copies and a
// reallocation ladder per triangle (O(Nt x Nae) overall). Split and merge
// (collapse_nonneighbor{1,2}) still copy O(Nae) rows, but once per split rather
// than once per triangle. The logical row order m = 1..n is exactly the old
// AeList order, so every predicate and tie-break sees the same sequence.
// ---------------------------------------------------------------------------
struct FrontRing {
  std::vector<Edge> buf;  // capacity is always a power of two (empty until used)
  std::size_t head = 0;
  int n = 0;  // the live row count (the old Nae)

  std::size_t mask() const { return buf.size() - 1; }
  Edge& at(int m) {
    return buf[(head + static_cast<std::size_t>(m) - 1) & mask()];
  }
  const Edge& at(int m) const {
    return buf[(head + static_cast<std::size_t>(m) - 1) & mask()];
  }
  void clear() {
    head = 0;
    n = 0;
  }
  void push_back(const Edge& e) {
    if (buf.empty()) {
      buf.resize(64);  // lazy: pool slots for depths never reached stay free
    } else if (static_cast<std::size_t>(n) == buf.size()) {
      grow();
    }
    buf[(head + static_cast<std::size_t>(n)) & mask()] = e;
    ++n;
  }
  void pop_front() {
    head = (head + 1) & mask();
    --n;
  }
  void pop_back() { --n; }

 private:
  void grow() {
    std::vector<Edge> nb(buf.size() * 2);
    for (int m = 1; m <= n; ++m) nb[static_cast<std::size_t>(m - 1)] = at(m);
    buf = std::move(nb);
    head = 0;
  }
};
using Tri3 = std::array<int, 3>;        // 1-based [a,b,c]

// ---------------------------------------------------------------------------
// Mutable context standing in for the MATLAB globals threaded through the
// advancing front (active list, nactive, rj, Rj).
// ---------------------------------------------------------------------------
struct Ctx {
  std::vector<ActiveFront> active;  // 1-based: active[1..nactive] (active[0] dummy)
  int nactive = 0;
  double rj = 0.0;
  const std::vector<double>* Rj = nullptr;  // GLOBAL segment index -> neighbour VdW radius
  int depth = 0;                            // advancing-front recursion depth
  // Shared grow-only Dist scratch for sweep 2. Safe across the recursion: every
  // recursive call site breaks out of the while loop immediately after, so a
  // parent frame never reads Dist again once a child may have overwritten it,
  // and each iteration fully rewrites Dist[1..Nae] before any read.
  std::vector<double> dist;
  // The patch's main front, a scratch ring for building the inactive fronts,
  // and per-recursion-depth spare rings for the collapse_nonneighbor splits
  // (Ae1 stays alive across the first recursive call and Ae2 across both, so
  // they cannot share one slot; the depth is unique per live frame, which makes
  // ring_pool[depth] collision-free). All grow-only, so front (re)builds stop
  // allocating in the steady state.
  FrontRing root;
  FrontRing scratch;
  // deque, NOT vector: a deeper frame growing the pool appends new slots, and a
  // deque keeps references to existing elements valid -- the parent frames hold
  // live references (including their own Ae) into earlier slots.
  std::deque<std::array<FrontRing, 2>> ring_pool;
};

// Fetch the depth-indexed spare-ring pair for the current live frame.
inline std::array<FrontRing, 2>& split_rings(Ctx& ctx) {
  while (ctx.ring_pool.size() <= static_cast<std::size_t>(ctx.depth)) {
    ctx.ring_pool.emplace_back();
  }
  return ctx.ring_pool[static_cast<std::size_t>(ctx.depth)];
}

// The frame loop is already bounded -- `while (k <= N)` increments k on every
// iteration -- so the only unbounded path through the mesher is the mutual
// recursion advancing_front_approach <-> collapse_nonneighbor{1,2}_sphere, where
// collapse_nonneighbor2_sphere can hand the callee a LARGER front. A relaxed-FP
// build (MESHMS_FP=fast) may flip the boundary tests that decide how the front
// splits, so cap the depth rather than risk a stack overflow. A frame that hits
// the cap returns without closing its front, leaving a hole the equivalence gate
// sees as extra boundary edges -- a detectable degradation instead of a crash.
// A strict build never reaches it (the golden suite is the proof) and the guard
// costs one compare per frame, never per triangle.
constexpr int kMaxFrontDepth = 4096;

struct FrontDepthGuard {
  Ctx& ctx;
  const bool ok;
  explicit FrontDepthGuard(Ctx& c) : ctx(c), ok(c.depth < kMaxFrontDepth) { ++ctx.depth; }
  ~FrontDepthGuard() { --ctx.depth; }
  FrontDepthGuard(const FrontDepthGuard&) = delete;
  FrontDepthGuard& operator=(const FrontDepthGuard&) = delete;
};

// ---------------------------------------------------------------------------
// small geometric helpers (port of the local functions in
// advancing_front_approach.m). x ** 2 on a SCALAR -> pysq().
// ---------------------------------------------------------------------------
// Outward unit normal of the sphere at point x (normal_sphere.m).
Vec3 normal_sphere(const Vec3& c_sphere, const Vec3& x) {
  Vec3 d = x - c_sphere;
  return d / std::sqrt(pysq(d.x) + pysq(d.y) + pysq(d.z));
}

#ifdef MESHMS_FP_FAST
// A deploy build compares angles without taking acos. The pair (grp, c) --
// grp 0 when the orientation t < 0 (angle = acos(c), in [0,pi]) and grp 1
// otherwise (angle = 2pi - acos(c), in [pi,2pi]) -- orders exactly like the
// angle it stands for: the groups are separated at pi, the angle is decreasing
// in c inside grp 0 and increasing in c inside grp 1. Every angle produced in
// this file is only ever COMPARED (never consumed as a radian value), so the
// whole mesher loses its acos calls -- ~10% of serial self time on Apple M2.
// A tie at exactly pi can resolve differently than the double comparison did;
// that is inside the fast contract and covered by test_fp_gate.
struct AngleT {
  int grp;
  double c;  // clamped cosine of the in-group acos argument
};
inline bool operator<(const AngleT& x, const AngleT& y) {
  if (x.grp != y.grp) return x.grp < y.grp;
  return x.grp == 0 ? x.c > y.c : x.c < y.c;
}
inline bool operator>(const AngleT& x, const AngleT& y) { return y < x; }
inline double clamp_pm1(double x) {
  return x > 1.0 ? 1.0 : (x < -1.0 ? -1.0 : x);
}
// Thresholds live in [pi,2pi] (grp 1), where angle = 2pi - acos(c) gives
// c = cos(angle) directly: cos(5pi/4) = -sqrt(2)/2, cos(5pi/3) = 1/2 (exact).
inline constexpr AngleT kFivePiOver4{1, -std::numbers::sqrt2 / 2.0};
inline constexpr AngleT kFivePiOver3{1, 0.5};
#else
using AngleT = double;
inline constexpr double kFivePiOver4 = 5 * std::numbers::pi / 4;
inline constexpr double kFivePiOver3 = 5.0 / 3.0 * std::numbers::pi;
#endif

// Angle (counterclockwise) between two neighbour edges e and f; n the sphere
// normal (angle_sphere.m). e,f are [tail,head] 1-based index pairs.
AngleT angle_sphere(const Edge& e, const Edge& f, const Vec3& n,
                    const std::vector<Vec3>& P) {
  Vec3 u = P[static_cast<std::size_t>(e[0])] - P[static_cast<std::size_t>(e[1])];
  Vec3 v = P[static_cast<std::size_t>(f[1])] - P[static_cast<std::size_t>(f[0])];
  u = u - dot(u, n) * n;  // project u onto the tangent plane
  v = v - dot(v, n) * n;
  double t = sign(dot(u, cross(v, n)));  // sign(det([u;v;n]))
  double nu = std::sqrt(pysq(u.x) + pysq(u.y) + pysq(u.z));
  double nv = std::sqrt(pysq(v.x) + pysq(v.y) + pysq(v.z));
#ifdef MESHMS_FP_FAST
  if (nu == 0.0 || nv == 0.0) return AngleT{0, 1.0};  // angle 0
  return AngleT{t < 0 ? 0 : 1, clamp_pm1(dot(u, v) / (nu * nv))};
#else
  if (nu == 0.0 || nv == 0.0) return 0.0;  // degenerate; avoid div-by-zero
  double base = acos_clamped(dot(u, v) / (nu * nv));
  if (t < 0) return base;
  return TWO_PI - base;
#endif
}

// Angle (counterclockwise) between two vectors around axis n (angle_vectors.m).
AngleT angle_vectors(Vec3 u, Vec3 v, const Vec3& n) {
  u = u - dot(u, n) * n;
  v = v - dot(v, n) * n;
  double t = sign(dot(u, cross(v, n)));
  double nu = norm(u);
  double nv = norm(v);
#ifdef MESHMS_FP_FAST
  if (nu == 0.0 || nv == 0.0) return AngleT{0, 1.0};  // angle 0
  return AngleT{t < 0 ? 0 : 1, clamp_pm1(dot(u, v) / (nu * nv))};
#else
  if (nu == 0.0 || nv == 0.0) return 0.0;  // degenerate; avoid div-by-zero
  double base = acos_clamped(dot(u, v) / (nu * nv));
  if (t < 0) return base;
  return TWO_PI - base;
#endif
}

// Tangent vector pointing outwards of the edge e (ne_sphere.m).
Vec3 ne_sphere(const Edge& e, const Vec3& n1, const Vec3& n2,
               const std::vector<Vec3>& P) {
  Vec3 u = P[static_cast<std::size_t>(e[1])] - P[static_cast<std::size_t>(e[0])];
  Vec3 v = 0.5 * n1 + 0.5 * n2;
  u = u - dot(u, v) * v;  // project u onto the tangent plane
  Vec3 ne{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
  ne = ne / std::sqrt(pysq(ne.x) + pysq(ne.y) + pysq(ne.z));
  return ne;
}

// Test point off the edge e by height h (testpoint_sphere.m).
Vec3 testpoint_sphere(const Edge& e, const Vec3& ne,
                      const std::vector<Vec3>& P, double h) {
  return 0.5 * (P[static_cast<std::size_t>(e[0])] + P[static_cast<std::size_t>(e[1])]) + h * ne;
}

// Map any point p radially onto the sphere (map_sphere.m).
Vec3 map_sphere(const Vec3& c_sphere, double r_sphere, const Vec3& p, double Rp) {
  double r;
  if (r_sphere > Rp)  // mesh the spherical patch on the SES
    r = r_sphere - Rp;
  else
    r = r_sphere;
  Vec3 u = (p - c_sphere) / norm(p - c_sphere);
  return c_sphere + r * u;
}

// ---------------------------------------------------------------------------
// loop / arc / circle division (local to mesh_sphpat.m). Appends to P, Ae.
// ---------------------------------------------------------------------------
// Divide an arc into several points / edges (arc_division.m). c centre, r radius,
// P1/P2 start/end, angle from P1 to P2 (clockwise), n points to k1.
void arc_division(const Vec3& c, double r, const Vec3& P1, const Vec3& P2,
                  double angle, const Vec3& n, FrontRing& Ae,
                  std::vector<Vec3>& P, int& Np, double r_sphere, int flag,
                  double d, double Rp, Ctx& ctx) {
  double theta0 = std::numbers::pi / 3;  // maximum angle change
  double N_division = std::max(std::floor(angle / theta0) + 1.0,
                               std::floor(r * angle / d) + 1.0);
  if (flag == 1) {
    double N_arc = std::floor(angle / theta0) + 1.0;
    N_division = N_arc * (std::floor(r * angle / d / N_arc) + 1.0);
  }

  // angle <= TWO_PI (vs MATLAB's strict <) lets full circles be aligned too.
  if (flag == 0 && angle <= TWO_PI && r_sphere > Rp && r_sphere < ctx.rj + Rp) {
    double r_new = r * (ctx.rj * r_sphere) / ((r_sphere - Rp) * (ctx.rj + Rp));
    N_division = std::max(std::floor(angle / theta0) + 1.0,
                          std::floor(r_new * angle / d) + 1.0);
  }

  int Ndiv = static_cast<int>(N_division);

  // The degenerate-arc early return depends only on Ndiv, the endpoints and the
  // angle -- not on the generated points -- so test it BEFORE generating and push
  // straight into P/Ae. Same values in the same order; the per-arc points/edges
  // scratch vectors (two heap allocations per arc segment) disappear.
  if (Ndiv == 1 && norm(P2 - P1) < 1e-10 && angle < 0.1) {
    return;
  }

  Vec3 u = (P1 - c) / r;
  Vec3 v{n.y * u.z - n.z * u.y, n.z * u.x - n.x * u.z, n.x * u.y - n.y * u.x};

  for (int j = 0; j < Ndiv; ++j) {
    double angle_j = static_cast<double>(j) / Ndiv * angle;
    Vec3 P_j = r * std::cos(angle_j) * u + r * std::sin(angle_j) * v + c;
    P.push_back(P_j);
    Ae.push_back(Edge{Np + j + 1, Np + j + 2});  // ring push tracks the count
  }
  Np = Np + Ndiv;
}

// Divide a loop into several edges (loop_division.m). loop[1..loopsize] are
// indices into segment0; segment0 is the (*,12) record matrix.
void loop_division(const Vec3& c_sphere, double r_sphere, const Loop& loop,
                   int loopsize,
                   const std::vector<std::array<double, 12>>& segment0,
                   FrontRing& Ae, std::vector<Vec3>& P, int& Np, double d,
                   double Rp, Ctx& ctx) {
  Ae.clear();
  int Np0 = Np;

  for (int i = 1; i <= loopsize; ++i) {
    int s = loop[static_cast<std::size_t>(i)];  // the segment number
    const std::array<double, 12>& seg = segment0[static_cast<std::size_t>(s)];
    Vec3 c{seg[1], seg[2], seg[3]};
    Vec3 n{seg[4], seg[5], seg[6]};
    double r = seg[7];
    Vec3 P1{seg[8], seg[9], seg[10]};
    // P2 = start point of the next segment (MATLAB mod(i,loopsize)+1).
    int s_next = loop[static_cast<std::size_t>(i % loopsize + 1)];
    const std::array<double, 12>& seg_next = segment0[static_cast<std::size_t>(s_next)];
    Vec3 P2{seg_next[8], seg_next[9], seg_next[10]};
    double angle = seg[11];

    if (r_sphere == Rp)
      ctx.rj = Rp;
    else
      ctx.rj = (*ctx.Rj)[static_cast<std::size_t>(s)];

    int flag = 0;
    if (r < Rp && r_sphere > Rp) flag = 1;
    if (r_sphere > Rp) {  // mesh the spherical patch on the SES
      c = c_sphere + (c - c_sphere) * (r_sphere - Rp) / r_sphere;
      r = r * (r_sphere - Rp) / r_sphere;
      P1 = c_sphere + (P1 - c_sphere) * (r_sphere - Rp) / r_sphere;
    }

    arc_division(c, r, P1, P2, angle, n, Ae, P, Np, r_sphere, flag, d, Rp, ctx);
  }

  if (Ae.n == 0) return;

  // close the loop: last edge head wraps to the first added point.
  Ae.at(Ae.n)[1] = Np0 + 1;
}

// Divide a full circle into edges (circle_division.m). circle is the 1-based row
// [_, c(3), n(3), r, torusR].
void circle_division(const Vec3& c_sphere, double r_sphere,
                     const std::array<double, 9>& circle, FrontRing& Ae,
                     std::vector<Vec3>& P, int& Np, double d, double Rp, Ctx& ctx) {
  Ae.clear();

  Vec3 c{circle[1], circle[2], circle[3]};
  Vec3 n{circle[4], circle[5], circle[6]};
  double r = circle[7];

  int flag = 0;
  if (r < Rp && r_sphere > Rp) flag = 1;

  if (r_sphere > Rp) {  // mesh the spherical patch on the SES
    c = c_sphere + (c - c_sphere) * (r_sphere - Rp) / r_sphere;
    r = r * (r_sphere - Rp) / r_sphere;
  }

  auto [v1, v2] = orthogonalvectors(n);  // n pointing outside
  (void)v2;
  Vec3 P1 = c + r * v1;

  // rj for arc_division's near-cusp refinement; mirror loop_division.
  if (r_sphere == Rp) {
    ctx.rj = Rp;
  } else if (circle[8] > 0.0) {
    // Full circle on a convex patch: align to the toroidal side (circle0 col 8).
    ctx.rj = circle[8];
  }

  int Np0 = Np;
  arc_division(c, r, P1, P1, TWO_PI, n, Ae, P, Np, r_sphere, flag, d, Rp, ctx);

  // modify the active edge set: close the circle.
  Ae.at(Ae.n)[1] = Np0 + 1;
}

// ---------------------------------------------------------------------------
// forward decls for the recursion.
// ---------------------------------------------------------------------------
void advancing_front_approach(const Vec3& c_sphere, double r_sphere, int N,
                              std::vector<Tri3>& T, int& Nt, FrontRing& Ae,
                              std::vector<Vec3>& P, int& Np, double d,
                              double tolerance, double Rp, Ctx& ctx);

void collapse_neighbor_sphere(int index_case, const Vec3& c_sphere,
                              double r_sphere, std::vector<Tri3>& T, int& Nt,
                              FrontRing& Ae, std::vector<Vec3>& P, int& Np,
                              double tolerance, double Rp);

void collapse_nonneighbor1_sphere(int index_np, const Vec3& c_sphere,
                                  double r_sphere, int N, std::vector<Tri3>& T,
                                  int& Nt, FrontRing& Ae,
                                  std::vector<Vec3>& P, int& Np, double d,
                                  double tolerance, double Rp, Ctx& ctx);

void collapse_nonneighbor2_sphere(int index_np, const Vec3& c_sphere,
                                  double r_sphere, int N, int index_nactive,
                                  std::vector<Tri3>& T, int& Nt, FrontRing& Ae, std::vector<Vec3>& P, int& Np,
                                  double d, double tolerance, double Rp, Ctx& ctx);

void addnewpoint_sphere(const Vec3& x, std::vector<Tri3>& T, int& Nt, FrontRing& Ae, std::vector<Vec3>& P, int& Np);

// ---------------------------------------------------------------------------
// advancing front (advancing_front_approach.m).
// ---------------------------------------------------------------------------
void advancing_front_approach(const Vec3& c_sphere, double r_sphere, int N,
                              std::vector<Tri3>& T, int& Nt, FrontRing& Ae,
                              std::vector<Vec3>& P, int& Np, double d,
                              double tolerance, double Rp, Ctx& ctx) {
  const FrontDepthGuard depth_guard(ctx);
  if (!depth_guard.ok) return;

  const double theta = 0.25;
  const AngleT alpha1 = kFivePiOver3;  // the angle condition param (5pi/3)
  const double h = d * std::sqrt(3.0) / 2.0;

  auto Pt = [&](int idx) -> const Vec3& { return P[static_cast<std::size_t>(idx)]; };
  int& Nae = Ae.n;
  auto E = [&](int idx) -> Edge& { return Ae.at(idx); };

  // Dist lives in ctx (see Ctx::dist): one allocation per patch instead of one
  // per afa frame, and grow-only so the reallocation count is O(log max Nae).
  std::vector<double>& Dist = ctx.dist;

  int k = 1;
  while (k <= N) {
    k += 1;

    if (Nae < 3) {
      Ae.clear();
      break;
    }

    if (Nae == 3) {  // the end condition
      T.push_back(Tri3{E(1)[0], E(1)[1], E(2)[1]});
      Nt = Nt + 1;
      Ae.clear();
      break;
    }

    if (Nae >= 4 || Nt == 1) {
      Vec3 n1 = normal_sphere(c_sphere, Pt(E(1)[0]));
      Vec3 n2 = normal_sphere(c_sphere, Pt(E(1)[1]));
      AngleT angle1 = angle_sphere(E(Nae), E(1), n1, P);  // the left angle
      AngleT angle2 = angle_sphere(E(1), E(2), n2, P);    // the right angle

      Vec3 ne = ne_sphere(E(1), n1, n2, P);  // tangent vector outwards
      Vec3 p = testpoint_sphere(E(1), ne, P, h);  // the test point
      Vec3 x = map_sphere(c_sphere, r_sphere, p, Rp);  // x is on the sphere

      // ---------------------------------------------------------------
      // dist / edgedist condition (including neighbours)
      // ---------------------------------------------------------------
      int index_np = 0;
      int index_nactive = 0;
      Vec3 PA1 = Pt(E(1)[0]);
      Vec3 PA2 = Pt(E(1)[1]);

      double edgedist = 0.0;
      for (int m = 3; m <= Nae; ++m) {
        Vec3 Pm1 = Pt(E(m)[0]);
        double Dist1_m = norm(PA1 - Pm1);
        double Dist2_m = norm(PA2 - Pm1);
        if (Dist1_m < tolerance || Dist2_m < tolerance ||
            Dist1_m + Dist2_m < tolerance * theta + norm(PA1 - PA2)) {
          if (dot(ne, Pm1 - PA1) > 0 &&
              (E(1)[0] != E(m)[0] && E(1)[1] != E(m)[0])) {
            int flag = 1;
            if (m == 3 && angle2 < kFivePiOver4) {
              flag = 0;
            } else if (m == Nae && angle1 < kFivePiOver4) {
              flag = 0;
            } else if (m > 3 && m < Nae) {
              AngleT angle1_m = angle_vectors(Pt(E(Nae)[0]) - Pt(E(1)[0]),
                                              Pt(E(m)[0]) - Pt(E(1)[0]), n1);
              AngleT angle2_m = angle_vectors(Pt(E(m)[0]) - Pt(E(2)[0]),
                                              Pt(E(2)[1]) - Pt(E(2)[0]), n2);
              if (angle1_m < angle1 || angle2_m < angle2) flag = 0;
            }

            if (flag == 1) {
              if (index_np == 0) {
                index_np = m;
                edgedist = Dist1_m + Dist2_m;
              } else if (Dist1_m + Dist2_m < edgedist) {
                index_np = m;
                edgedist = Dist1_m + Dist2_m;
              } else if (E(m)[0] == E(index_np)[0]) {
                Vec3 n_m = (Pt(E(m)[0]) - c_sphere) / norm(Pt(E(m)[0]) - c_sphere);
                AngleT a1 = angle_vectors(Pt(E(1)[0]) - Pt(E(m)[0]),
                                          Pt(E(m)[1]) - Pt(E(m)[0]), n_m);
                AngleT a2 = angle_vectors(Pt(E(1)[0]) - Pt(E(index_np)[0]),
                                          Pt(E(index_np)[1]) - Pt(E(index_np)[0]),
                                          n_m);
                if (m == Nae || a1 > a2) index_np = m;
              }
            }
          }
        }
      }

      for (int i = 1; i <= ctx.nactive; ++i) {  // nonneighbour point on other loop
        if (ctx.active[static_cast<std::size_t>(i)].meshed == 0) {
          AeList& Ae0 = ctx.active[static_cast<std::size_t>(i)].Ae;
          int Nae0 = ctx.active[static_cast<std::size_t>(i)].Nae;
          for (int m = 1; m <= Nae0; ++m) {
            double dist1 = norm(Pt(E(1)[0]) - Pt(Ae0[static_cast<std::size_t>(m)][0]));
            double dist2 = norm(Pt(E(1)[1]) - Pt(Ae0[static_cast<std::size_t>(m)][0]));
            if ((dist1 < tolerance || dist2 < tolerance ||
                 dist1 + dist2 < tolerance * theta + norm(Pt(E(1)[0]) - Pt(E(1)[1]))) &&
                dot(ne, Pt(Ae0[static_cast<std::size_t>(m)][0]) - Pt(E(1)[0])) > 0) {
              if (index_np == 0) {
                index_np = m;
                edgedist = dist1 + dist2;
                index_nactive = i;
              } else if (dist1 + dist2 < edgedist) {
                index_np = m;
                edgedist = dist1 + dist2;
                index_nactive = i;
              }
            }
          }
        }
      }

      if (index_np > 1) {
        if (index_nactive == 0) {
          if (index_np > 3 && index_np < Nae) {
            collapse_nonneighbor1_sphere(index_np, c_sphere, r_sphere, N, T, Nt,
                                         Ae, P, Np, d, tolerance, Rp, ctx);
            break;
          } else if (index_np == Nae) {
            collapse_neighbor_sphere(1, c_sphere, r_sphere, T, Nt, Ae, P, Np,
                                     tolerance, Rp);
          } else if (index_np == 3) {
            collapse_neighbor_sphere(2, c_sphere, r_sphere, T, Nt, Ae, P, Np,
                                     tolerance, Rp);
          }
          continue;
        } else if (index_nactive > 0) {
          collapse_nonneighbor2_sphere(index_np, c_sphere, r_sphere, N,
                                       index_nactive, T, Nt, Ae, P, Np, d,
                                       tolerance, Rp, ctx);
          break;
        }
      }

      // ---------------------------------------------------------------
      // Test x: general dist condition with respect to x
      // ---------------------------------------------------------------
      index_np = 0;
      index_nactive = 0;

      PA1 = Pt(E(1)[0]);
      // Grow-only: Dist[1..Nae] is fully overwritten below before any read and
      // Dist[0] is never read, so the assign()'s zero-fill (an O(Nae) memset
      // every iteration) was pure waste.
      if (Dist.size() < static_cast<std::size_t>(Nae) + 1)
        Dist.resize(static_cast<std::size_t>(Nae) + 1);
      for (int m = 1; m <= Nae; ++m) {
        Dist[static_cast<std::size_t>(m)] = norm(x - Pt(E(m)[0]));
      }

      edgedist = 0.0;
      for (int m = 3; m <= Nae; ++m) {
        int m_next = m % Nae + 1;  // MATLAB mod(m,Nae)+1
        bool cond =
            Dist[static_cast<std::size_t>(m)] < tolerance ||
            Dist[static_cast<std::size_t>(m - 1)] + Dist[static_cast<std::size_t>(m)] <
                tolerance * theta + norm(Pt(E(m - 1)[0]) - Pt(E(m - 1)[1])) ||
            Dist[static_cast<std::size_t>(m)] + Dist[static_cast<std::size_t>(m_next)] <
                tolerance * theta + norm(Pt(E(m)[0]) - Pt(E(m)[1]));
        if (cond && dot(ne, Pt(E(m)[0]) - PA1) > 0 &&
            (E(1)[0] != E(m)[0] && E(1)[1] != E(m)[0])) {
          double dist1 = norm(Pt(E(1)[0]) - Pt(E(m)[0]));
          double dist2 = norm(Pt(E(1)[1]) - Pt(E(m)[0]));

          int flag = 1;
          if (m == 3 && angle2 < kFivePiOver4) {
            flag = 0;
          } else if (m == Nae && angle1 < kFivePiOver4) {
            flag = 0;
          } else if (m > 3 && m < Nae) {
            AngleT angle1_m = angle_vectors(Pt(E(Nae)[0]) - Pt(E(1)[0]),
                                            Pt(E(m)[0]) - Pt(E(1)[0]), n1);
            AngleT angle2_m = angle_vectors(Pt(E(m)[0]) - Pt(E(2)[0]),
                                            Pt(E(2)[1]) - Pt(E(2)[0]), n2);
            if (angle1_m < angle1 || angle2_m < angle2) flag = 0;
          }

          if (flag == 1) {
            if (index_np == 0) {
              index_np = m;
              edgedist = dist1 + dist2;
            } else if (dist1 + dist2 < edgedist) {
              index_np = m;
              edgedist = dist1 + dist2;
            } else if (E(m)[0] == E(index_np)[0]) {
              Vec3 n_m = (Pt(E(m)[0]) - c_sphere) / norm(Pt(E(m)[0]) - c_sphere);
              AngleT a1 = angle_vectors(Pt(E(1)[0]) - Pt(E(m)[0]),
                                        Pt(E(m)[1]) - Pt(E(m)[0]), n_m);
              AngleT a2 = angle_vectors(Pt(E(1)[0]) - Pt(E(index_np)[0]),
                                        Pt(E(index_np)[1]) - Pt(E(index_np)[0]),
                                        n_m);
              if (m == Nae || a1 > a2) index_np = m;
            }
          }
        }
      }

      for (int i = 1; i <= ctx.nactive; ++i) {
        if (ctx.active[static_cast<std::size_t>(i)].meshed == 0) {
          AeList& Ae0 = ctx.active[static_cast<std::size_t>(i)].Ae;
          int Nae0 = ctx.active[static_cast<std::size_t>(i)].Nae;
          for (int m = 1; m <= Nae0; ++m) {
            double dist = norm(x - Pt(Ae0[static_cast<std::size_t>(m)][0]));
            if (dist < tolerance &&
                dot(ne, Pt(Ae0[static_cast<std::size_t>(m)][0]) - Pt(E(1)[0])) > 0) {
              double dist1 = norm(Pt(E(1)[0]) - Pt(Ae0[static_cast<std::size_t>(m)][0]));
              double dist2 = norm(Pt(E(1)[1]) - Pt(Ae0[static_cast<std::size_t>(m)][0]));
              if (index_np == 0) {
                index_np = m;
                edgedist = dist1 + dist2;
                index_nactive = i;
              } else if (dist1 + dist2 < edgedist) {
                index_np = m;
                edgedist = dist1 + dist2;
                index_nactive = i;
              }
            }
          }
        }
      }

      if (index_np > 0) {
        if (index_nactive == 0) {
          if (index_np > 3 && index_np < Nae) {
            collapse_nonneighbor1_sphere(index_np, c_sphere, r_sphere, N, T, Nt,
                                         Ae, P, Np, d, tolerance, Rp, ctx);
            break;
          } else if (index_np == Nae) {
            collapse_neighbor_sphere(1, c_sphere, r_sphere, T, Nt, Ae, P, Np,
                                     tolerance, Rp);
          } else if (index_np == 3) {
            collapse_neighbor_sphere(2, c_sphere, r_sphere, T, Nt, Ae, P, Np,
                                     tolerance, Rp);
          }
          continue;
        } else if (index_nactive > 0) {
          collapse_nonneighbor2_sphere(index_np, c_sphere, r_sphere, N,
                                       index_nactive, T, Nt, Ae, P, Np, d,
                                       tolerance, Rp, ctx);
          break;
        }
      }

      // ---------------------------------------------------------------
      // angle condition
      // ---------------------------------------------------------------
      if (angle1 > alpha1 && angle2 > alpha1) {
        double edgedist_Nae =
            norm(Pt(E(Nae)[0]) - Pt(E(1)[0])) + norm(Pt(E(Nae)[0]) - Pt(E(1)[1]));
        double edgedist_2 =
            norm(Pt(E(2)[0]) - Pt(E(1)[0])) + norm(Pt(E(2)[0]) - Pt(E(1)[1]));
        if (edgedist_Nae < edgedist_2) {
          collapse_neighbor_sphere(1, c_sphere, r_sphere, T, Nt, Ae, P, Np,
                                   tolerance, Rp);
        } else {
          collapse_neighbor_sphere(2, c_sphere, r_sphere, T, Nt, Ae, P, Np,
                                   tolerance, Rp);
        }
        continue;
      } else if (angle1 > alpha1) {
        AngleT angle1_2 = angle_vectors(Pt(E(Nae)[0]) - Pt(E(2)[0]),
                                        Pt(E(2)[1]) - Pt(E(2)[0]), n2);
        if (angle1_2 > angle2) {
          collapse_neighbor_sphere(1, c_sphere, r_sphere, T, Nt, Ae, P, Np,
                                   tolerance, Rp);
          continue;
        }
      } else if (angle2 > alpha1) {
        AngleT angle1_1 = angle_vectors(Pt(E(Nae)[0]) - Pt(E(1)[0]),
                                        Pt(E(2)[1]) - Pt(E(1)[0]), n1);
        if (angle1_1 > angle1) {
          collapse_neighbor_sphere(2, c_sphere, r_sphere, T, Nt, Ae, P, Np,
                                   tolerance, Rp);
          continue;
        }
      }

      // ---------------------------------------------------------------
      // Add a new point
      // ---------------------------------------------------------------
      addnewpoint_sphere(x, T, Nt, Ae, P, Np);
    }
  }
}

// ---------------------------------------------------------------------------
// collapse / add-point helpers
// ---------------------------------------------------------------------------
void collapse_neighbor_sphere(int index_case, const Vec3& c_sphere,
                              double r_sphere, std::vector<Tri3>& T, int& Nt,
                              FrontRing& Ae, std::vector<Vec3>& P, int& Np,
                              double tolerance, double Rp) {
  const double gamma = 1.8;  // controls edge length
  auto Pt = [&](int idx) -> const Vec3& { return P[static_cast<std::size_t>(idx)]; };
  int& Nae = Ae.n;
  auto E = [&](int idx) -> Edge& { return Ae.at(idx); };

  if (index_case == 1) {  // connect the left edge
    if (norm(Pt(E(Nae)[0]) - Pt(E(1)[1])) < gamma * tolerance) {
      T.push_back(Tri3{E(Nae)[0], E(1)[1], E(1)[0]});
      Nt = Nt + 1;
      // Ae=[Ae(2:Nae-1,:);Ae(Nae,1),Ae(1,2)] -- drop edges 1 and Nae, append.
      Edge tail{E(Nae)[0], E(1)[1]};
      Ae.pop_back();
      Ae.pop_front();
      Ae.push_back(tail);
    } else {
      Vec3 p = 0.5 * (Pt(E(Nae)[0]) + Pt(E(1)[1]));
      Vec3 x = map_sphere(c_sphere, r_sphere, p, Rp);
      if (Nae == 4 && norm(x - Pt(E(3)[0])) < 1e-10) {
        T.push_back(Tri3{E(Nae)[0], E(3)[0], E(1)[0]});
        T.push_back(Tri3{E(3)[0], E(1)[1], E(1)[0]});
        Nt = Nt + 2;
        Ae.clear();
      } else {
        P.push_back(x);
        Np = Np + 1;
        T.push_back(Tri3{E(Nae)[0], Np, E(1)[0]});
        T.push_back(Tri3{Np, E(1)[1], E(1)[0]});
        Nt = Nt + 2;
        Edge t1{E(Nae)[0], Np};
        Edge t2{Np, E(1)[1]};
        Ae.pop_back();
        Ae.pop_front();
        Ae.push_back(t1);
        Ae.push_back(t2);
      }
    }
  } else {  // connect the right edge
    if (norm(Pt(E(1)[0]) - Pt(E(2)[1])) < gamma * tolerance) {
      T.push_back(Tri3{E(1)[0], E(2)[1], E(1)[1]});
      Nt = Nt + 1;
      // Ae=[Ae(3:Nae,:);Ae(1,1),Ae(2,2)] -- drop edges 1 and 2, append.
      Edge tail{E(1)[0], E(2)[1]};
      Ae.pop_front();
      Ae.pop_front();
      Ae.push_back(tail);
    } else {
      Vec3 p = 0.5 * (Pt(E(1)[0]) + Pt(E(2)[1]));
      Vec3 x = map_sphere(c_sphere, r_sphere, p, Rp);
      if (Nae == 4 && norm(x - Pt(E(Nae)[0])) < 1e-10) {
        T.push_back(Tri3{E(1)[0], E(Nae)[0], E(2)[0]});
        T.push_back(Tri3{E(Nae)[0], E(2)[1], E(2)[0]});
        Nt = Nt + 2;
        Ae.clear();
      } else {
        P.push_back(x);
        Np = Np + 1;
        T.push_back(Tri3{E(1)[0], Np, E(2)[0]});
        T.push_back(Tri3{Np, E(2)[1], E(2)[0]});
        Nt = Nt + 2;
        Edge t1{E(1)[0], Np};
        Edge t2{Np, E(2)[1]};
        Ae.pop_front();
        Ae.pop_front();
        Ae.push_back(t1);
        Ae.push_back(t2);
      }
    }
  }
}

void collapse_nonneighbor1_sphere(int index_np, const Vec3& c_sphere,
                                  double r_sphere, int N, std::vector<Tri3>& T,
                                  int& Nt, FrontRing& Ae,
                                  std::vector<Vec3>& P, int& Np, double d,
                                  double tolerance, double Rp, Ctx& ctx) {
  const double gamma = 1.8;
  int m = index_np;
  auto Pt = [&](int idx) -> const Vec3& { return P[static_cast<std::size_t>(idx)]; };
  int& Nae = Ae.n;
  auto E = [&](int idx) -> Edge& { return Ae.at(idx); };

  std::array<FrontRing, 2>& spare = split_rings(ctx);
  FrontRing& Ae1 = spare[0];
  FrontRing& Ae2 = spare[1];
  Ae1.clear();
  Ae2.clear();

  double dm1 = norm(Pt(E(m)[0]) - Pt(E(1)[0]));
  double dm2 = norm(Pt(E(m)[0]) - Pt(E(1)[1]));

  if (dm1 <= gamma * tolerance && dm2 <= gamma * tolerance) {
    T.push_back(Tri3{E(m)[0], E(1)[1], E(1)[0]});
    Nt = Nt + 1;
    for (int kk = 2; kk <= m - 1; ++kk) Ae1.push_back(E(kk));
    Ae1.push_back(Edge{E(m)[0], E(1)[1]});   // Nae1 = m - 1
    for (int kk = m; kk <= Nae; ++kk) Ae2.push_back(E(kk));
    Ae2.push_back(Edge{E(1)[0], E(m)[0]});   // Nae2 = Nae - m + 2
  } else if (dm1 > gamma * tolerance && dm2 > gamma * tolerance) {
    Vec3 p = 0.5 * (Pt(E(1)[0]) + Pt(E(m)[0]));
    Vec3 x = map_sphere(c_sphere, r_sphere, p, Rp);
    P.push_back(x);
    Np = Np + 1;
    p = 0.5 * (Pt(E(2)[0]) + Pt(E(m)[0]));
    x = map_sphere(c_sphere, r_sphere, p, Rp);
    P.push_back(x);
    Np = Np + 1;

    if (norm(Pt(E(m)[0]) - Pt(E(1)[0])) > norm(Pt(E(m)[0]) - Pt(E(1)[1]))) {
      T.push_back(Tri3{E(1)[0], Np - 1, E(2)[0]});
      T.push_back(Tri3{Np - 1, E(m)[0], Np});
      T.push_back(Tri3{Np, E(2)[0], Np - 1});
      Nt = Nt + 3;
    } else {
      T.push_back(Tri3{E(1)[0], Np - 1, Np});
      T.push_back(Tri3{Np - 1, E(m)[0], Np});
      T.push_back(Tri3{Np, E(2)[0], E(1)[0]});
      Nt = Nt + 3;
    }

    for (int kk = 2; kk <= m - 1; ++kk) Ae1.push_back(E(kk));
    Ae1.push_back(Edge{E(m)[0], Np});
    Ae1.push_back(Edge{Np, E(1)[1]});        // Nae1 = m

    for (int kk = m; kk <= Nae; ++kk) Ae2.push_back(E(kk));
    Ae2.push_back(Edge{E(1)[0], Np - 1});
    Ae2.push_back(Edge{Np - 1, E(m)[0]});    // Nae2 = Nae - m + 3
  } else if (dm1 > gamma * tolerance) {
    Vec3 p = 0.5 * (Pt(E(1)[0]) + Pt(E(m)[0]));
    Vec3 x = map_sphere(c_sphere, r_sphere, p, Rp);
    P.push_back(x);
    Np = Np + 1;

    T.push_back(Tri3{E(1)[0], Np, E(2)[0]});
    T.push_back(Tri3{Np, E(m)[0], E(1)[1]});
    Nt = Nt + 2;

    for (int kk = 2; kk <= m - 1; ++kk) Ae1.push_back(E(kk));
    Ae1.push_back(Edge{E(m)[0], E(1)[1]});   // Nae1 = m - 1
    for (int kk = m; kk <= Nae; ++kk) Ae2.push_back(E(kk));
    Ae2.push_back(Edge{E(1)[0], Np});
    Ae2.push_back(Edge{Np, E(m)[0]});        // Nae2 = Nae - m + 3
  } else if (dm2 > gamma * tolerance) {
    Vec3 p = 0.5 * (Pt(E(2)[0]) + Pt(E(m)[0]));
    Vec3 x = map_sphere(c_sphere, r_sphere, p, Rp);
    P.push_back(x);
    Np = Np + 1;

    T.push_back(Tri3{E(m)[0], Np, E(1)[0]});
    T.push_back(Tri3{Np, E(2)[0], E(1)[0]});
    Nt = Nt + 2;

    for (int kk = 2; kk <= m - 1; ++kk) Ae1.push_back(E(kk));
    Ae1.push_back(Edge{E(m)[0], Np});
    Ae1.push_back(Edge{Np, E(1)[1]});        // Nae1 = m
    for (int kk = m; kk <= Nae; ++kk) Ae2.push_back(E(kk));
    Ae2.push_back(Edge{E(1)[0], E(m)[0]});   // Nae2 = Nae - m + 2
  }

  advancing_front_approach(c_sphere, r_sphere, N, T, Nt, Ae1, P, Np, d,
                           tolerance, Rp, ctx);
  advancing_front_approach(c_sphere, r_sphere, N, T, Nt, Ae2, P, Np, d,
                           tolerance, Rp, ctx);
  // mirror Python: Ae is rebound to the (post-second-call) front. swap rather
  // than move, so the parent's old buffer stays in the pool slot for reuse.
  std::swap(Ae, Ae2);
}

void collapse_nonneighbor2_sphere(int index_np, const Vec3& c_sphere,
                                  double r_sphere, int N, int index_nactive,
                                  std::vector<Tri3>& T, int& Nt, FrontRing& Ae, std::vector<Vec3>& P, int& Np,
                                  double d, double tolerance, double Rp, Ctx& ctx) {
  ctx.active[static_cast<std::size_t>(index_nactive)].meshed = 1;
  AeList& Ae0 = ctx.active[static_cast<std::size_t>(index_nactive)].Ae;
  int Nae0 = ctx.active[static_cast<std::size_t>(index_nactive)].Nae;
  int& Nae = Ae.n;
  auto E = [&](int idx) -> Edge& { return Ae.at(idx); };
  auto E0 = [&](int idx) -> Edge& { return Ae0[static_cast<std::size_t>(idx)]; };

  int m = index_np;
  T.push_back(Tri3{E0(m)[0], E(1)[1], E(1)[0]});
  Nt = Nt + 1;

  // Ae1=[Ae(2:Nae,:); Ae(1,1),Ae0(m,1); Ae0(m:Nae0,:); Ae0(1:m-1,:); Ae0(m,1),Ae(1,2)]
  // Built in a depth-pool ring; the row count lands at Nae + Nae0 + 1 by
  // construction.
  FrontRing& Ae1 = split_rings(ctx)[0];
  Ae1.clear();
  for (int kk = 2; kk <= Nae; ++kk) Ae1.push_back(E(kk));
  Ae1.push_back(Edge{E(1)[0], E0(m)[0]});
  for (int kk = m; kk <= Nae0; ++kk) Ae1.push_back(E0(kk));
  for (int kk = 1; kk <= m - 1; ++kk) Ae1.push_back(E0(kk));
  Ae1.push_back(Edge{E0(m)[0], E(1)[1]});

  advancing_front_approach(c_sphere, r_sphere, N, T, Nt, Ae1, P, Np, d,
                           tolerance, Rp, ctx);
  std::swap(Ae, Ae1);
}

void addnewpoint_sphere(const Vec3& x, std::vector<Tri3>& T, int& Nt,
                        FrontRing& Ae, std::vector<Vec3>& P, int& Np) {
  auto E = [&](int idx) -> Edge& { return Ae.at(idx); };
  P.push_back(x);
  Np = Np + 1;
  T.push_back(Tri3{E(1)[0], Np, E(1)[1]});
  Nt = Nt + 1;
  // Ae=[Ae(2:Nae,:);Ae(1,1),Np;Np,Ae(1,2)] -- drop edge 1, append two.
  Edge t1{E(1)[0], Np};
  Edge t2{Np, E(1)[1]};
  Ae.pop_front();
  Ae.push_back(t1);
  Ae.push_back(t2);
}

}  // namespace

// ---------------------------------------------------------------------------
// compute_NV (port of Output/compute_NV.m)
// ---------------------------------------------------------------------------
std::vector<Vec3> compute_NV(const std::vector<std::array<int, 3>>& T,
                             const std::vector<Vec3>& P, const Vec3& c_sphere,
                             double arg_NV) {
  std::vector<Vec3> NV;
  NV.reserve(T.size());
  for (const auto& tri : T) {
    const Vec3& a = P[static_cast<std::size_t>(tri[0])];
    const Vec3& b = P[static_cast<std::size_t>(tri[1])];
    const Vec3& c = P[static_cast<std::size_t>(tri[2])];

    Vec3 V1 = b - a;
    Vec3 V2 = c - a;

    Vec3 nv = cross(V1, V2);
    nv = nv / norm(nv);

    double s = sign(dot(nv, a - c_sphere));

    if (s != arg_NV) nv = nv * -1.0;

    NV.push_back(nv);
  }
  return NV;
}

// ---------------------------------------------------------------------------
// mesh_sphpat (entry point, port of mesh_sphpat.m)
// ---------------------------------------------------------------------------
LocalMesh mesh_sphpat(const Vec3& c_sphere, double r_sphere,
                      const std::vector<Loop>& loops,
                      const std::vector<std::array<double, 12>>& segment0,
                      const std::vector<std::array<double, 9>>& circle0,
                      const std::vector<int>& patches, int patchesize, double Rp,
                      double d, const std::vector<double>* Rj, Tag boundary_tag) {
  LocalMesh out;  // emit=false until the old add_patch path is reached.
  double tolerance = 0.8 * std::min(d, r_sphere);

  Ctx ctx;
  ctx.nactive = patchesize - 1;
  ctx.Rj = Rj;

  // N denotes the max number of triangles on the spherical patch.
  int N = static_cast<int>(std::floor(
      4 * std::numbers::pi * pysq(r_sphere) /
      (std::sqrt(3.0) * pysq(tolerance) / 4)));
  std::vector<Tri3> T;
  int Nt = 0;
  std::vector<Vec3> P(1);  // 1-based, P[0] dummy
  int Np = 0;

  // The patch's main front lives in ctx (grow-only across patches);
  // loop_division/circle_division clear and rebuild it before any use, and the
  // r_sphere == 0 path never touches it.
  FrontRing& Ae = ctx.root;
  Ae.clear();

  if (r_sphere != 0) {
    int k = patches[1];  // the k-th loop or the -k-th circle
    if (k > 0) {
      loop_division(c_sphere, r_sphere, loops[static_cast<std::size_t>(k)],
                    static_cast<int>(loops[static_cast<std::size_t>(k)].size()) - 1,
                    segment0, Ae, P, Np, d, Rp, ctx);
    } else {
      circle_division(c_sphere, r_sphere, circle0[static_cast<std::size_t>(-k)], Ae,
                      P, Np, d, Rp, ctx);
    }

    // build the inactive fronts (other loops/circles of this patch).
    ctx.active.assign(static_cast<std::size_t>(ctx.nactive) + 1, ActiveFront{});

    for (int i = 1; i <= ctx.nactive; ++i) {
      int kk = patches[static_cast<std::size_t>(i + 1)];
      // Build in the scratch ring, then materialise as the plain 1-based AeList
      // ActiveFront stores (a public-header type the sweeps read directly).
      FrontRing& ring0 = ctx.scratch;
      if (kk > 0) {
        loop_division(c_sphere, r_sphere, loops[static_cast<std::size_t>(kk)],
                      static_cast<int>(loops[static_cast<std::size_t>(kk)].size()) - 1,
                      segment0, ring0, P, Np, d, Rp, ctx);
      } else {
        circle_division(c_sphere, r_sphere, circle0[static_cast<std::size_t>(-kk)],
                        ring0, P, Np, d, Rp, ctx);
      }
      AeList Ae0(1);
      Ae0.reserve(static_cast<std::size_t>(ring0.n) + 1);
      for (int mm = 1; mm <= ring0.n; ++mm) Ae0.push_back(ring0.at(mm));
      ctx.active[static_cast<std::size_t>(i)] = ActiveFront{0, std::move(Ae0), ring0.n};
    }

    // P[1..np_boundary] are the loop/circle division points (the patch boundary).
    int np_boundary = Np;

    advancing_front_approach(c_sphere, r_sphere, N, T, Nt, Ae, P, Np, d,
                             tolerance, Rp, ctx);

    if (Nt > 0) {  // there might exist an isolated point as a patch
      double arg_NV;
      if (r_sphere == Rp)
        arg_NV = -1;  // concave
      else
        arg_NV = 1;  // convex
      std::vector<Vec3> NV = compute_NV(T, P, c_sphere, arg_NV);

      std::vector<TagList> vids;
      if (boundary_tag.kind != 0) {
        // {boundary_tag} for the np_boundary boundary vertices, empty TagList
        // (interior) for the advancing-front interior vertices.
        vids.assign(static_cast<std::size_t>(Np), TagList{});
        for (int i = 0; i < np_boundary; ++i)
          vids[static_cast<std::size_t>(i)] = TagList{boundary_tag};
      }

      // RETURN the local mesh; the driver calls add_patch when emit==true.
      // T is std::vector<std::array<int,3>> (Tri3 is exactly that).
      out.P = std::move(P);
      out.T = std::move(T);
      out.NV = std::move(NV);
      out.vids = std::move(vids);
      out.emit = true;
    }
  }
  return out;
}


// ---------------------------------------------------------------------------
// Parallel ordered merge (see meshing.hpp).
// ---------------------------------------------------------------------------
void merge_local_meshes(MeshState& state, const std::vector<LocalMesh*>& lms) {
  const std::size_t n = lms.size();

  // Serial prefix sum of the write bases. These are exactly the
  // V.size()/F.size()/N.size() values the sequential add_patch calls would have
  // seen, so every element lands at the same index with the same value.
  std::vector<std::size_t> bv(n), bf(n), bn(n);
  std::size_t v = state.V.size(), f = state.F.size(), nn = state.N.size();
  for (std::size_t i = 0; i < n; ++i) {
    const LocalMesh& lm = *lms[i];
    bv[i] = v;
    bf[i] = f;
    bn[i] = nn;
    v += lm.P.empty() ? 0 : lm.P.size() - 1;
    f += lm.T.size();
    nn += lm.NV.size();  // add_patch appends N only when normals are present
  }
  state.V.resize(v);
  state.tags.resize(v);      // default TagList{} == the empty-vids case
  state.vatom.resize(v, 0);  // 0 == the empty-patch_vatom case
  state.F.resize(f);
  state.N.resize(nn);

  meshms::parallel_for(0, static_cast<int>(n), [&](int i) {
    LocalMesh& lm = *lms[static_cast<std::size_t>(i)];
    const std::size_t base = bv[static_cast<std::size_t>(i)];
    const std::size_t fbase = bf[static_cast<std::size_t>(i)];
    const std::size_t nbase = bn[static_cast<std::size_t>(i)];
    const std::size_t Np = lm.P.empty() ? 0 : lm.P.size() - 1;
    for (std::size_t k = 1; k <= Np; ++k) state.V[base + k - 1] = lm.P[k];
    if (!lm.vids.empty()) {
      for (std::size_t k = 0; k < Np; ++k)
        state.tags[base + k] = std::move(lm.vids[k]);
    }
    if (!lm.vatom.empty()) {
      for (std::size_t k = 0; k < Np; ++k) state.vatom[base + k] = lm.vatom[k];
    }
    for (std::size_t t = 0; t < lm.T.size(); ++t) {
      const auto& tt = lm.T[t];
      state.F[fbase + t] =
          Tri{static_cast<int32_t>(tt[0] - 1 + static_cast<int>(base)),
              static_cast<int32_t>(tt[1] - 1 + static_cast<int>(base)),
              static_cast<int32_t>(tt[2] - 1 + static_cast<int>(base))};
    }
    for (std::size_t t = 0; t < lm.NV.size(); ++t) state.N[nbase + t] = lm.NV[t];
    // Free the dead patch here, on the worker: the serial merge previously also
    // serialised every cross-thread free.
    lm = LocalMesh{};
  });
}

}  // namespace meshms
