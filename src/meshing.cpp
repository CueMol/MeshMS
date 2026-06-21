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

#include "meshms/vec3.hpp"

namespace meshms {

namespace {

using Edge = std::array<int, 2>;        // [tail, head], 1-based
using AeList = std::vector<Edge>;       // 1-based list, Ae[0] dummy
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
};

// ---------------------------------------------------------------------------
// Helpers for the 1-based Ae list.
// ---------------------------------------------------------------------------
// An empty 1-based Ae list (only the dummy slot).
AeList new_ae() { return AeList(1); }

// Build a fresh 1-based Ae list from a [first,last) range of [tail,head] rows.
AeList ae_from_rows(const AeList& src, int lo, int hi) {
  // src[lo..hi] inclusive (1-based, may be empty if lo > hi).
  AeList out(1);
  for (int k = lo; k <= hi; ++k) out.push_back(src[static_cast<std::size_t>(k)]);
  return out;
}

// ---------------------------------------------------------------------------
// small geometric helpers (port of the local functions in
// advancing_front_approach.m). x ** 2 on a SCALAR -> pysq().
// ---------------------------------------------------------------------------
// Outward unit normal of the sphere at point x (normal_sphere.m).
Vec3 normal_sphere(const Vec3& c_sphere, const Vec3& x) {
  Vec3 d = x - c_sphere;
  return d / std::sqrt(pysq(d.x) + pysq(d.y) + pysq(d.z));
}

double acos_clamped_local(double x) {
  if (x > 1.0) x = 1.0;
  else if (x < -1.0) x = -1.0;
  return std::acos(x);
}

// Angle (counterclockwise) between two neighbour edges e and f; n the sphere
// normal (angle_sphere.m). e,f are [tail,head] 1-based index pairs.
double angle_sphere(const Edge& e, const Edge& f, const Vec3& n,
                    const std::vector<Vec3>& P) {
  Vec3 u = P[static_cast<std::size_t>(e[0])] - P[static_cast<std::size_t>(e[1])];
  Vec3 v = P[static_cast<std::size_t>(f[1])] - P[static_cast<std::size_t>(f[0])];
  u = u - dot(u, n) * n;  // project u onto the tangent plane
  v = v - dot(v, n) * n;
  double t = sign(dot(u, cross(v, n)));  // sign(det([u;v;n]))
  double nu = std::sqrt(pysq(u.x) + pysq(u.y) + pysq(u.z));
  double nv = std::sqrt(pysq(v.x) + pysq(v.y) + pysq(v.z));
  if (nu == 0.0 || nv == 0.0) return 0.0;  // degenerate; avoid div-by-zero
  double base = acos_clamped_local(dot(u, v) / (nu * nv));
  if (t < 0) return base;
  return TWO_PI - base;
}

// Angle (counterclockwise) between two vectors around axis n (angle_vectors.m).
double angle_vectors(Vec3 u, Vec3 v, const Vec3& n) {
  u = u - dot(u, n) * n;
  v = v - dot(v, n) * n;
  double t = sign(dot(u, cross(v, n)));
  double nu = norm(u);
  double nv = norm(v);
  if (nu == 0.0 || nv == 0.0) return 0.0;  // degenerate; avoid div-by-zero
  double base = acos_clamped_local(dot(u, v) / (nu * nv));
  if (t < 0) return base;
  return TWO_PI - base;
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
                  double angle, const Vec3& n, AeList& Ae, int& Nae,
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

  Vec3 u = (P1 - c) / r;
  Vec3 v{n.y * u.z - n.z * u.y, n.z * u.x - n.x * u.z, n.x * u.y - n.y * u.x};

  std::vector<Vec3> points;
  std::vector<Edge> edges;
  for (int j = 0; j < Ndiv; ++j) {
    double angle_j = static_cast<double>(j) / Ndiv * angle;
    Vec3 P_j = r * std::cos(angle_j) * u + r * std::sin(angle_j) * v + c;
    points.push_back(P_j);
    edges.push_back(Edge{Np + j + 1, Np + j + 2});
  }

  if (Ndiv == 1 && norm(P2 - P1) < 1e-10 && angle < 0.1) {
    return;
  }

  for (const Vec3& pt : points) P.push_back(pt);
  Np = Np + Ndiv;

  for (const Edge& e : edges) Ae.push_back(e);
  Nae = Nae + Ndiv;
}

// Divide a loop into several edges (loop_division.m). loop[1..loopsize] are
// indices into segment0; segment0 is the (*,12) record matrix.
void loop_division(const Vec3& c_sphere, double r_sphere, const Loop& loop,
                   int loopsize,
                   const std::vector<std::array<double, 12>>& segment0,
                   AeList& Ae, int& Nae, std::vector<Vec3>& P, int& Np, double d,
                   double Rp, Ctx& ctx) {
  Ae = new_ae();
  Nae = 0;
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

    arc_division(c, r, P1, P2, angle, n, Ae, Nae, P, Np, r_sphere, flag, d, Rp, ctx);
  }

  if (Nae == 0) return;

  // close the loop: last edge head wraps to the first added point.
  Ae[static_cast<std::size_t>(Nae)][1] = Np0 + 1;
}

// Divide a full circle into edges (circle_division.m). circle is the 1-based row
// [_, c(3), n(3), r, torusR].
void circle_division(const Vec3& c_sphere, double r_sphere,
                     const std::array<double, 9>& circle, AeList& Ae, int& Nae,
                     std::vector<Vec3>& P, int& Np, double d, double Rp, Ctx& ctx) {
  Ae = new_ae();
  Nae = 0;

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
  arc_division(c, r, P1, P1, TWO_PI, n, Ae, Nae, P, Np, r_sphere, flag, d, Rp, ctx);

  // modify the active edge set: close the circle.
  Ae[static_cast<std::size_t>(Nae)][1] = Np0 + 1;
}

// ---------------------------------------------------------------------------
// forward decls for the recursion.
// ---------------------------------------------------------------------------
void advancing_front_approach(const Vec3& c_sphere, double r_sphere, int N,
                              std::vector<Tri3>& T, int& Nt, AeList& Ae, int& Nae,
                              std::vector<Vec3>& P, int& Np, double d,
                              double tolerance, double Rp, Ctx& ctx);

void collapse_neighbor_sphere(int index_case, const Vec3& c_sphere,
                              double r_sphere, std::vector<Tri3>& T, int& Nt,
                              AeList& Ae, int& Nae, std::vector<Vec3>& P, int& Np,
                              double tolerance, double Rp);

void collapse_nonneighbor1_sphere(int index_np, const Vec3& c_sphere,
                                  double r_sphere, int N, std::vector<Tri3>& T,
                                  int& Nt, AeList& Ae, int& Nae,
                                  std::vector<Vec3>& P, int& Np, double d,
                                  double tolerance, double Rp, Ctx& ctx);

void collapse_nonneighbor2_sphere(int index_np, const Vec3& c_sphere,
                                  double r_sphere, int N, int index_nactive,
                                  std::vector<Tri3>& T, int& Nt, AeList& Ae,
                                  int& Nae, std::vector<Vec3>& P, int& Np,
                                  double d, double tolerance, double Rp, Ctx& ctx);

void addnewpoint_sphere(const Vec3& x, std::vector<Tri3>& T, int& Nt, AeList& Ae,
                        int& Nae, std::vector<Vec3>& P, int& Np);

// ---------------------------------------------------------------------------
// advancing front (advancing_front_approach.m).
// ---------------------------------------------------------------------------
void advancing_front_approach(const Vec3& c_sphere, double r_sphere, int N,
                              std::vector<Tri3>& T, int& Nt, AeList& Ae, int& Nae,
                              std::vector<Vec3>& P, int& Np, double d,
                              double tolerance, double Rp, Ctx& ctx) {
  const double theta = 0.25;
  const double alpha1 = 5.0 / 3.0 * std::numbers::pi;  // the angle condition param
  const double h = d * std::sqrt(3.0) / 2.0;

  auto Pt = [&](int idx) -> const Vec3& { return P[static_cast<std::size_t>(idx)]; };
  auto E = [&](int idx) -> Edge& { return Ae[static_cast<std::size_t>(idx)]; };

  // One Dist scratch per afa frame, reused across while iterations (S3): it is
  // re-assigned to Nae+1 and fully overwritten by the m=1..Nae fill before any
  // read each iteration, so contents are identical -- only allocations drop.
  std::vector<double> Dist;

  int k = 1;
  while (k <= N) {
    k += 1;

    if (Nae < 3) {
      Ae = new_ae();
      Nae = 0;
      break;
    }

    if (Nae == 3) {  // the end condition
      T.push_back(Tri3{E(1)[0], E(1)[1], E(2)[1]});
      Nt = Nt + 1;
      Nae = 0;
      Ae = new_ae();
      break;
    }

    if (Nae >= 4 || Nt == 1) {
      Vec3 n1 = normal_sphere(c_sphere, Pt(E(1)[0]));
      Vec3 n2 = normal_sphere(c_sphere, Pt(E(1)[1]));
      double angle1 = angle_sphere(E(Nae), E(1), n1, P);  // the left angle
      double angle2 = angle_sphere(E(1), E(2), n2, P);    // the right angle

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
            if (m == 3 && angle2 < 5 * std::numbers::pi / 4) {
              flag = 0;
            } else if (m == Nae && angle1 < 5 * std::numbers::pi / 4) {
              flag = 0;
            } else if (m > 3 && m < Nae) {
              double angle1_m = angle_vectors(Pt(E(Nae)[0]) - Pt(E(1)[0]),
                                              Pt(E(m)[0]) - Pt(E(1)[0]), n1);
              double angle2_m = angle_vectors(Pt(E(m)[0]) - Pt(E(2)[0]),
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
                double a1 = angle_vectors(Pt(E(1)[0]) - Pt(E(m)[0]),
                                          Pt(E(m)[1]) - Pt(E(m)[0]), n_m);
                double a2 = angle_vectors(Pt(E(1)[0]) - Pt(E(index_np)[0]),
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
                                         Ae, Nae, P, Np, d, tolerance, Rp, ctx);
            break;
          } else if (index_np == Nae) {
            collapse_neighbor_sphere(1, c_sphere, r_sphere, T, Nt, Ae, Nae, P, Np,
                                     tolerance, Rp);
          } else if (index_np == 3) {
            collapse_neighbor_sphere(2, c_sphere, r_sphere, T, Nt, Ae, Nae, P, Np,
                                     tolerance, Rp);
          }
          continue;
        } else if (index_nactive > 0) {
          collapse_nonneighbor2_sphere(index_np, c_sphere, r_sphere, N,
                                       index_nactive, T, Nt, Ae, Nae, P, Np, d,
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
      Dist.assign(static_cast<std::size_t>(Nae) + 1, 0.0);  // Dist[0] dummy
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
          if (m == 3 && angle2 < 5 * std::numbers::pi / 4) {
            flag = 0;
          } else if (m == Nae && angle1 < 5 * std::numbers::pi / 4) {
            flag = 0;
          } else if (m > 3 && m < Nae) {
            double angle1_m = angle_vectors(Pt(E(Nae)[0]) - Pt(E(1)[0]),
                                            Pt(E(m)[0]) - Pt(E(1)[0]), n1);
            double angle2_m = angle_vectors(Pt(E(m)[0]) - Pt(E(2)[0]),
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
              double a1 = angle_vectors(Pt(E(1)[0]) - Pt(E(m)[0]),
                                        Pt(E(m)[1]) - Pt(E(m)[0]), n_m);
              double a2 = angle_vectors(Pt(E(1)[0]) - Pt(E(index_np)[0]),
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
                                         Ae, Nae, P, Np, d, tolerance, Rp, ctx);
            break;
          } else if (index_np == Nae) {
            collapse_neighbor_sphere(1, c_sphere, r_sphere, T, Nt, Ae, Nae, P, Np,
                                     tolerance, Rp);
          } else if (index_np == 3) {
            collapse_neighbor_sphere(2, c_sphere, r_sphere, T, Nt, Ae, Nae, P, Np,
                                     tolerance, Rp);
          }
          continue;
        } else if (index_nactive > 0) {
          collapse_nonneighbor2_sphere(index_np, c_sphere, r_sphere, N,
                                       index_nactive, T, Nt, Ae, Nae, P, Np, d,
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
          collapse_neighbor_sphere(1, c_sphere, r_sphere, T, Nt, Ae, Nae, P, Np,
                                   tolerance, Rp);
        } else {
          collapse_neighbor_sphere(2, c_sphere, r_sphere, T, Nt, Ae, Nae, P, Np,
                                   tolerance, Rp);
        }
        continue;
      } else if (angle1 > alpha1) {
        double angle1_2 = angle_vectors(Pt(E(Nae)[0]) - Pt(E(2)[0]),
                                        Pt(E(2)[1]) - Pt(E(2)[0]), n2);
        if (angle1_2 > angle2) {
          collapse_neighbor_sphere(1, c_sphere, r_sphere, T, Nt, Ae, Nae, P, Np,
                                   tolerance, Rp);
          continue;
        }
      } else if (angle2 > alpha1) {
        double angle1_1 = angle_vectors(Pt(E(Nae)[0]) - Pt(E(1)[0]),
                                        Pt(E(2)[1]) - Pt(E(1)[0]), n1);
        if (angle1_1 > angle1) {
          collapse_neighbor_sphere(2, c_sphere, r_sphere, T, Nt, Ae, Nae, P, Np,
                                   tolerance, Rp);
          continue;
        }
      }

      // ---------------------------------------------------------------
      // Add a new point
      // ---------------------------------------------------------------
      addnewpoint_sphere(x, T, Nt, Ae, Nae, P, Np);
    }
  }
}

// ---------------------------------------------------------------------------
// collapse / add-point helpers
// ---------------------------------------------------------------------------
void collapse_neighbor_sphere(int index_case, const Vec3& c_sphere,
                              double r_sphere, std::vector<Tri3>& T, int& Nt,
                              AeList& Ae, int& Nae, std::vector<Vec3>& P, int& Np,
                              double tolerance, double Rp) {
  const double gamma = 1.8;  // controls edge length
  auto Pt = [&](int idx) -> const Vec3& { return P[static_cast<std::size_t>(idx)]; };
  auto E = [&](int idx) -> Edge& { return Ae[static_cast<std::size_t>(idx)]; };

  if (index_case == 1) {  // connect the left edge
    if (norm(Pt(E(Nae)[0]) - Pt(E(1)[1])) < gamma * tolerance) {
      T.push_back(Tri3{E(Nae)[0], E(1)[1], E(1)[0]});
      Nt = Nt + 1;
      // Ae=[Ae(2:Nae-1,:);Ae(Nae,1),Ae(1,2)]
      Edge tail{E(Nae)[0], E(1)[1]};
      AeList nw = ae_from_rows(Ae, 2, Nae - 1);
      nw.push_back(tail);
      Ae = std::move(nw);
      Nae = Nae - 1;
    } else {
      Vec3 p = 0.5 * (Pt(E(Nae)[0]) + Pt(E(1)[1]));
      Vec3 x = map_sphere(c_sphere, r_sphere, p, Rp);
      if (Nae == 4 && norm(x - Pt(E(3)[0])) < 1e-10) {
        T.push_back(Tri3{E(Nae)[0], E(3)[0], E(1)[0]});
        T.push_back(Tri3{E(3)[0], E(1)[1], E(1)[0]});
        Nt = Nt + 2;
        Ae = new_ae();
        Nae = 0;
      } else {
        P.push_back(x);
        Np = Np + 1;
        T.push_back(Tri3{E(Nae)[0], Np, E(1)[0]});
        T.push_back(Tri3{Np, E(1)[1], E(1)[0]});
        Nt = Nt + 2;
        Edge t1{E(Nae)[0], Np};
        Edge t2{Np, E(1)[1]};
        AeList nw = ae_from_rows(Ae, 2, Nae - 1);
        nw.push_back(t1);
        nw.push_back(t2);
        Ae = std::move(nw);
        Nae = static_cast<int>(Ae.size()) - 1;
      }
    }
  } else {  // connect the right edge
    if (norm(Pt(E(1)[0]) - Pt(E(2)[1])) < gamma * tolerance) {
      T.push_back(Tri3{E(1)[0], E(2)[1], E(1)[1]});
      Nt = Nt + 1;
      // Ae=[Ae(3:Nae,:);Ae(1,1),Ae(2,2)]
      Edge tail{E(1)[0], E(2)[1]};
      AeList nw = ae_from_rows(Ae, 3, Nae);
      nw.push_back(tail);
      Ae = std::move(nw);
      Nae = Nae - 1;
    } else {
      Vec3 p = 0.5 * (Pt(E(1)[0]) + Pt(E(2)[1]));
      Vec3 x = map_sphere(c_sphere, r_sphere, p, Rp);
      if (Nae == 4 && norm(x - Pt(E(Nae)[0])) < 1e-10) {
        T.push_back(Tri3{E(1)[0], E(Nae)[0], E(2)[0]});
        T.push_back(Tri3{E(Nae)[0], E(2)[1], E(2)[0]});
        Nt = Nt + 2;
        Ae = new_ae();
        Nae = 0;
      } else {
        P.push_back(x);
        Np = Np + 1;
        T.push_back(Tri3{E(1)[0], Np, E(2)[0]});
        T.push_back(Tri3{Np, E(2)[1], E(2)[0]});
        Nt = Nt + 2;
        Edge t1{E(1)[0], Np};
        Edge t2{Np, E(2)[1]};
        AeList nw = ae_from_rows(Ae, 3, Nae);
        nw.push_back(t1);
        nw.push_back(t2);
        Ae = std::move(nw);
        Nae = static_cast<int>(Ae.size()) - 1;
      }
    }
  }
}

void collapse_nonneighbor1_sphere(int index_np, const Vec3& c_sphere,
                                  double r_sphere, int N, std::vector<Tri3>& T,
                                  int& Nt, AeList& Ae, int& Nae,
                                  std::vector<Vec3>& P, int& Np, double d,
                                  double tolerance, double Rp, Ctx& ctx) {
  const double gamma = 1.8;
  int m = index_np;
  auto Pt = [&](int idx) -> const Vec3& { return P[static_cast<std::size_t>(idx)]; };
  auto E = [&](int idx) -> Edge& { return Ae[static_cast<std::size_t>(idx)]; };

  AeList Ae1;
  AeList Ae2;
  int Nae1 = 0;
  int Nae2 = 0;

  double dm1 = norm(Pt(E(m)[0]) - Pt(E(1)[0]));
  double dm2 = norm(Pt(E(m)[0]) - Pt(E(1)[1]));

  if (dm1 <= gamma * tolerance && dm2 <= gamma * tolerance) {
    T.push_back(Tri3{E(m)[0], E(1)[1], E(1)[0]});
    Nt = Nt + 1;
    Ae1 = ae_from_rows(Ae, 2, m - 1);
    Ae1.push_back(Edge{E(m)[0], E(1)[1]});
    Nae1 = m - 1;
    Ae2 = ae_from_rows(Ae, m, Nae);
    Ae2.push_back(Edge{E(1)[0], E(m)[0]});
    Nae2 = Nae - m + 2;
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

    Ae1 = ae_from_rows(Ae, 2, m - 1);
    Ae1.push_back(Edge{E(m)[0], Np});
    Ae1.push_back(Edge{Np, E(1)[1]});
    Nae1 = m;

    Ae2 = ae_from_rows(Ae, m, Nae);
    Ae2.push_back(Edge{E(1)[0], Np - 1});
    Ae2.push_back(Edge{Np - 1, E(m)[0]});
    Nae2 = Nae - m + 3;
  } else if (dm1 > gamma * tolerance) {
    Vec3 p = 0.5 * (Pt(E(1)[0]) + Pt(E(m)[0]));
    Vec3 x = map_sphere(c_sphere, r_sphere, p, Rp);
    P.push_back(x);
    Np = Np + 1;

    T.push_back(Tri3{E(1)[0], Np, E(2)[0]});
    T.push_back(Tri3{Np, E(m)[0], E(1)[1]});
    Nt = Nt + 2;

    Ae1 = ae_from_rows(Ae, 2, m - 1);
    Ae1.push_back(Edge{E(m)[0], E(1)[1]});
    Nae1 = m - 1;
    Ae2 = ae_from_rows(Ae, m, Nae);
    Ae2.push_back(Edge{E(1)[0], Np});
    Ae2.push_back(Edge{Np, E(m)[0]});
    Nae2 = Nae - m + 3;
  } else if (dm2 > gamma * tolerance) {
    Vec3 p = 0.5 * (Pt(E(2)[0]) + Pt(E(m)[0]));
    Vec3 x = map_sphere(c_sphere, r_sphere, p, Rp);
    P.push_back(x);
    Np = Np + 1;

    T.push_back(Tri3{E(m)[0], Np, E(1)[0]});
    T.push_back(Tri3{Np, E(2)[0], E(1)[0]});
    Nt = Nt + 2;

    Ae1 = ae_from_rows(Ae, 2, m - 1);
    Ae1.push_back(Edge{E(m)[0], Np});
    Ae1.push_back(Edge{Np, E(1)[1]});
    Nae1 = m;
    Ae2 = ae_from_rows(Ae, m, Nae);
    Ae2.push_back(Edge{E(1)[0], E(m)[0]});
    Nae2 = Nae - m + 2;
  }

  advancing_front_approach(c_sphere, r_sphere, N, T, Nt, Ae1, Nae1, P, Np, d,
                           tolerance, Rp, ctx);
  advancing_front_approach(c_sphere, r_sphere, N, T, Nt, Ae2, Nae2, P, Np, d,
                           tolerance, Rp, ctx);
  // mirror Python: Ae,Nae are rebound to the (post-second-call) values.
  Ae = std::move(Ae2);
  Nae = Nae2;
}

void collapse_nonneighbor2_sphere(int index_np, const Vec3& c_sphere,
                                  double r_sphere, int N, int index_nactive,
                                  std::vector<Tri3>& T, int& Nt, AeList& Ae,
                                  int& Nae, std::vector<Vec3>& P, int& Np,
                                  double d, double tolerance, double Rp, Ctx& ctx) {
  ctx.active[static_cast<std::size_t>(index_nactive)].meshed = 1;
  AeList& Ae0 = ctx.active[static_cast<std::size_t>(index_nactive)].Ae;
  int Nae0 = ctx.active[static_cast<std::size_t>(index_nactive)].Nae;
  auto E = [&](int idx) -> Edge& { return Ae[static_cast<std::size_t>(idx)]; };
  auto E0 = [&](int idx) -> Edge& { return Ae0[static_cast<std::size_t>(idx)]; };

  int m = index_np;
  T.push_back(Tri3{E0(m)[0], E(1)[1], E(1)[0]});
  Nt = Nt + 1;

  // Ae1=[Ae(2:Nae,:); Ae(1,1),Ae0(m,1); Ae0(m:Nae0,:); Ae0(1:m-1,:); Ae0(m,1),Ae(1,2)]
  AeList Ae1(1);
  for (int kk = 2; kk <= Nae; ++kk) Ae1.push_back(E(kk));
  Ae1.push_back(Edge{E(1)[0], E0(m)[0]});
  for (int kk = m; kk <= Nae0; ++kk) Ae1.push_back(E0(kk));
  for (int kk = 1; kk <= m - 1; ++kk) Ae1.push_back(E0(kk));
  Ae1.push_back(Edge{E0(m)[0], E(1)[1]});
  int Nae1 = Nae + Nae0 + 1;

  advancing_front_approach(c_sphere, r_sphere, N, T, Nt, Ae1, Nae1, P, Np, d,
                           tolerance, Rp, ctx);
  Ae = std::move(Ae1);
  Nae = Nae1;
}

void addnewpoint_sphere(const Vec3& x, std::vector<Tri3>& T, int& Nt, AeList& Ae,
                        int& Nae, std::vector<Vec3>& P, int& Np) {
  auto E = [&](int idx) -> Edge& { return Ae[static_cast<std::size_t>(idx)]; };
  P.push_back(x);
  Np = Np + 1;
  T.push_back(Tri3{E(1)[0], Np, E(1)[1]});
  Nt = Nt + 1;
  // Ae=[Ae(2:Nae,:);Ae(1,1),Np;Np,Ae(1,2)]
  Edge t1{E(1)[0], Np};
  Edge t2{Np, E(1)[1]};
  AeList nw = ae_from_rows(Ae, 2, Nae);
  nw.push_back(t1);
  nw.push_back(t2);
  Ae = std::move(nw);
  Nae = Nae + 1;
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

  AeList Ae = new_ae();
  int Nae = 0;

  if (r_sphere != 0) {
    int k = patches[1];  // the k-th loop or the -k-th circle
    if (k > 0) {
      loop_division(c_sphere, r_sphere, loops[static_cast<std::size_t>(k)],
                    static_cast<int>(loops[static_cast<std::size_t>(k)].size()) - 1,
                    segment0, Ae, Nae, P, Np, d, Rp, ctx);
    } else {
      circle_division(c_sphere, r_sphere, circle0[static_cast<std::size_t>(-k)], Ae,
                      Nae, P, Np, d, Rp, ctx);
    }

    // build the inactive fronts (other loops/circles of this patch).
    ctx.active.assign(static_cast<std::size_t>(ctx.nactive) + 1, ActiveFront{});

    for (int i = 1; i <= ctx.nactive; ++i) {
      int kk = patches[static_cast<std::size_t>(i + 1)];
      AeList Ae0;
      int Nae0 = 0;
      if (kk > 0) {
        loop_division(c_sphere, r_sphere, loops[static_cast<std::size_t>(kk)],
                      static_cast<int>(loops[static_cast<std::size_t>(kk)].size()) - 1,
                      segment0, Ae0, Nae0, P, Np, d, Rp, ctx);
      } else {
        circle_division(c_sphere, r_sphere, circle0[static_cast<std::size_t>(-kk)],
                        Ae0, Nae0, P, Np, d, Rp, ctx);
      }
      ctx.active[static_cast<std::size_t>(i)] = ActiveFront{0, std::move(Ae0), Nae0};
    }

    // P[1..np_boundary] are the loop/circle division points (the patch boundary).
    int np_boundary = Np;

    advancing_front_approach(c_sphere, r_sphere, N, T, Nt, Ae, Nae, P, Np, d,
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

}  // namespace meshms
