// Faithful port of the toroidal module's data_SEStorpat() (geometry + mesh path).
//
// See toroidal.hpp for the faithfulness contract. The local-function structure
// of the Python module is preserved: _sqrt_pos, coord_toroide, index_P_toroide,
// coord_cusp, index_P, orient_face_toroidal, mesh_toroide, mesh_cusp,
// data_SEStorpat.
#include "meshms/toroidal.hpp"
#include "meshms/parallel.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

#include "meshms/meshing.hpp"  // LocalMesh
#include "meshms/vec3.hpp"

namespace meshms {

namespace {

// sqrt guarding against tiny negatives (MATLAB silently took the real part).
inline double sqrt_pos(double x) { return x > 0.0 ? std::sqrt(x) : 0.0; }

// math.acos(max(-1, min(1, x))) -- the Python clamp before acos.
inline double acos_clamp(double x) {
  return std::acos(std::max(-1.0, std::min(1.0, x)));
}

// ---------------------------------------------------------------------------
// coord helpers (literal ports of the MATLAB local functions)
// ---------------------------------------------------------------------------

// Coordinate of the (i,j) point on a regular thin toroide.
Vec3 coord_toroide(int i, int j, int N_probe, int N_arc, const Vec3& P_k1,
                   double r, const Vec3& A, double Rp, double theta1,
                   double theta2, double angle, int direct, const Vec3& n,
                   const Vec3& ci, int flag) {
  double theta_i;
  if (flag == 1) {
    theta_i = i * (theta1 + theta2) / N_probe;
  } else {
    theta_i = i * (theta1 - theta2) / N_probe;
  }
  double d = std::abs(norm(ci - A) - std::tan(theta1 - theta_i) * r);
  Vec3 Ai = ci + d * n;

  Vec3 u = (P_k1 - A) / r;
  Vec3 v = direct * Vec3{n.y * u.z - n.z * u.y, n.z * u.x - n.x * u.z,
                         n.x * u.y - n.y * u.x};
  double angle_j = static_cast<double>(j) / N_arc * angle;
  Vec3 P_j = r * std::cos(angle_j) * u + r * std::sin(angle_j) * v + A;

  return P_j + Rp * (Ai - P_j) / norm(Ai - P_j);
}

// 1-based vertex indices for a structured (N_probe+1)x(N_arc+1) grid.
inline std::array<int, 3> index_P_toroide(int i1, int j1, int i2, int j2, int i3,
                                          int j3, int N_arc) {
  int k1 = i1 * (N_arc + 1) + j1 + 1;
  int k2 = i2 * (N_arc + 1) + j2 + 1;
  int k3 = i3 * (N_arc + 1) + j3 + 1;
  return {k1, k2, k3};
}

// Coordinate of the (i,j) point on a singular cusp toroide. `t` is 1-based
// (dummy t[0]) of per-row arc subdivision multipliers, exactly as in MATLAB.
Vec3 coord_cusp(int i, int j, int N_probe, int N_arc, const Vec3& P_k1,
                const Vec3& A1, double r, const Vec3& A, double Rp,
                double theta1, double angle, int direct, const Vec3& n,
                const std::vector<long long>& t) {
  if (i == 0) {
    return A1;
  }

  double theta = acos_clamp(r / Rp);
  double cos_arg = r / std::cos(theta + theta1 * i / N_probe);
  double radius_i = (cos_arg - Rp) / cos_arg * r;
  Vec3 Ai = A + sqrt_pos(pysq(cos_arg) - pysq(r)) * (A1 - A) / norm(A1 - A);

  double division = static_cast<double>(t[i + 1] * N_arc);

  Vec3 u = (P_k1 - A) / r;
  Vec3 v = direct * Vec3{n.y * u.z - n.z * u.y, n.z * u.x - n.x * u.z,
                         n.x * u.y - n.y * u.x};
  double angle_j = static_cast<double>(j) / division * angle;
  Vec3 P_j = r * std::cos(angle_j) * u + r * std::sin(angle_j) * v + A;

  return radius_i / r * P_j + (r - radius_i) / r * Ai;
}

// 1-based vertex indices for the variable-row cusp grid (row_off prefix sums).
inline std::array<int, 3> index_P(int i1, int j1, int i2, int j2, int i3, int j3,
                                  const std::vector<long long>& row_off) {
  return {static_cast<int>(row_off[i1]) + j1 + 1,
          static_cast<int>(row_off[i2]) + j2 + 1,
          static_cast<int>(row_off[i3]) + j3 + 1};
}

// Orient ONE face by the analytic torus normal `nrm` (the compute_NV_toroidal
// loop body, fused into the per-face normal computation at both call sites so
// the intermediate per-face `normal` array is never materialized; every float
// expression and its evaluation order are unchanged).
inline Vec3 orient_face_toroidal(const std::array<int, 3>& tri,
                                 const std::vector<Vec3>& P, const Vec3& nrm) {
  const Vec3& a = P[static_cast<std::size_t>(tri[0])];
  const Vec3& b = P[static_cast<std::size_t>(tri[1])];
  const Vec3& c = P[static_cast<std::size_t>(tri[2])];
  Vec3 V1 = b - a;
  Vec3 V2 = c - a;
  Vec3 nv = cross(V1, V2);
  nv = nv / norm(nv);
  if (dot(nv, nrm) < 0.0) {
    nv = nv * -1.0;
  }
  return nv;
}

// ---------------------------------------------------------------------------
// mesh_toroide : regular thin toroide
// ---------------------------------------------------------------------------
LocalMesh mesh_toroide(const Vec3& ci, const Vec3& cj, double ri,
                       double rj, int direct, const Vec3& n, double angle, double r,
                       const Vec3& A, const Vec3& P_k1, const Vec3& /*P_k2*/,
                       double theta1, double theta2, double Rp, double d, int k,
                       int atom_i, int atom_j, int e1, int e2) {
  int flag;
  double theta;
  if (norm(ci - cj) >
      std::max(sqrt_pos(pysq(ri + Rp) - pysq(r)),
               sqrt_pos(pysq(rj + Rp) - pysq(r)))) {
    flag = 1;
    theta = theta1 + theta2;
  } else {
    flag = 0;
    theta = std::abs(theta1 - theta2);
  }

  double theta0 = std::numbers::pi / 3.0;
  // to guarantee the length of each edge is greater than tolerance
  int N_probe = std::max(static_cast<int>(std::floor(Rp * theta / d)) + 1,
                         static_cast<int>(std::floor(theta / theta0)) + 1);
  int N_arc = std::max(
      static_cast<int>(std::floor(angle / theta0)) + 1,
      static_cast<int>(std::floor(r * angle * std::max(ri, rj) /
                                  ((Rp + std::max(ri, rj)) * d))) +
          1);

  int N_point = (N_probe + 1) * (N_arc + 1);
  // 1-based point storage (dummy row 0), real points 1..N_point
  std::vector<Vec3> P(static_cast<std::size_t>(N_point) + 1, Vec3{});
  int Np = 0;

  for (int i = 0; i <= N_probe; ++i) {
    for (int j = 0; j <= N_arc; ++j) {
      ++Np;
      P[static_cast<std::size_t>(Np)] =
          coord_toroide(i, j, N_probe, N_arc, P_k1, r, A, Rp, theta1, theta2,
                        angle, direct, n, ci, flag);
    }
  }

  std::vector<std::array<int, 3>> T;
  T.reserve(static_cast<std::size_t>(2) * N_probe * N_arc);  // exact face count
  for (int i = 1; i <= N_probe; ++i) {
    for (int j = 1; j <= N_arc; ++j) {
      T.push_back(index_P_toroide(i - 1, j, i, j - 1, i, j, N_arc));
      T.push_back(index_P_toroide(i - 1, j, i, j - 1, i - 1, j - 1, N_arc));
    }
  }

  // analytic torus normal at each face centroid, fused with the face orientation
  // (identical per-face expressions; the intermediate array is gone)
  std::vector<Vec3> NV;
  NV.reserve(T.size());
  for (const auto& tri : T) {
    Vec3 p = (P[static_cast<std::size_t>(tri[0])] +
              P[static_cast<std::size_t>(tri[1])] +
              P[static_cast<std::size_t>(tri[2])]) /
             3.0;
    Vec3 nt = (p - A) - dot(p - A, n) * n;
    Vec3 xp = A + r * nt / norm(nt);
    NV.push_back(orient_face_toroidal(tri, P, (xp - p) / norm(xp - p)));
  }

  // ID-fusion tags: row i=0 lies on atom_i's sphere, row i=N_probe on atom_j;
  // column j=0 is the probe arc at endpoint e1, j=N_arc at e2. A dual-tagged
  // corner stores the FULL Python tag list (e.g. {("atom",j),("tseam",k,gi)}),
  // so it bridges the convex and concave groups during fuse_by_id.
  // The Python gates the whole vids construction on `atom_i is not None`; here
  // atom_i is always supplied by data_SEStorpat (>=1 for atoms, never 0), so we
  // build vids whenever this overload is used with a real atom_i.
  std::vector<TagList> vids;
  // Per-vertex owning atom: the toroide spans atom_i (row gi=0) to atom_j (row
  // gi=N_probe); attribute each grid row to the nearer rim atom (the "neighbour
  // side" the plan allows). Aligned with P[1..Np] like vids; never affects V/F/N.
  std::vector<int32_t> vatom;
  if (atom_i != 0) {
    bool full_circle = angle >= TWO_PI - 1e-9;
    vids.assign(static_cast<std::size_t>((N_probe + 1) * (N_arc + 1)), TagList{});
    vatom.assign(static_cast<std::size_t>((N_probe + 1) * (N_arc + 1)), 0);
    for (int gi = 0; gi <= N_probe; ++gi) {
      for (int gj = 0; gj <= N_arc; ++gj) {
        // Collect tags in the SAME order as Python (the FULL list per vertex).
        TagList tg;
        if (gi == 0) tg.push_back(Tag{1, atom_i, 0});           // ("atom", atom_i)
        if (gi == N_probe) tg.push_back(Tag{1, atom_j, 0});     // ("atom", atom_j)
        if (e1 != 0 && gj == 0) tg.push_back(Tag{2, e1, 0});    // ("probe", e1)
        if (e2 != 0 && gj == N_arc) tg.push_back(Tag{2, e2, 0}); // ("probe", e2)
        if (full_circle && (gj == 0 || gj == N_arc))
          tg.push_back(Tag{3, k, gi});                          // ("tseam", k, gi)
        vids[static_cast<std::size_t>(gi * (N_arc + 1) + gj)] = std::move(tg);
        vatom[static_cast<std::size_t>(gi * (N_arc + 1) + gj)] =
            (gi <= N_probe / 2) ? atom_i : atom_j;
      }
    }
  }

  // RETURN the local mesh; the driver adds it (this path was unconditional, so
  // emit is always true here -- behaviour-preserving).
  LocalMesh out;
  out.P = std::move(P);
  out.T = std::move(T);
  out.NV = std::move(NV);
  out.vids = std::move(vids);
  out.vatom = std::move(vatom);
  out.emit = true;
  return out;
}

// ---------------------------------------------------------------------------
// mesh_cusp : singular self-intersecting toroide
// ---------------------------------------------------------------------------
LocalMesh mesh_cusp(const Vec3& /*c1*/, double /*r1*/,
                    const Vec3& A1, int direct, const Vec3& n, double angle,
                    double r, const Vec3& A, const Vec3& P_k1, const Vec3& /*P_k2*/,
                    double theta1, double Rp, double d, int k, int atom_rim, int e1,
                    int e2, bool tagged) {
  // corresponds to the arc on the probe at fixed point
  int N_probe = static_cast<int>(std::floor(Rp * theta1 / d)) + 1;
  double theta0 = std::numbers::pi / 3.0;
  int N_arc = static_cast<int>(std::floor(angle / theta0)) + 1;

  long long N_point = 1;

  // t is 1-based (dummy t[0]); MATLAB sets t(i+1) for i=1..N_probe, t(1)=0.
  std::vector<long long> t(static_cast<std::size_t>(N_probe) + 2, 0);
  for (int i = 1; i <= N_probe; ++i) {
    double theta_i = acos_clamp(r / Rp) + theta1 / N_probe * i;
    double L_i = r / std::cos(theta_i);
    double r_i = (L_i - Rp) / L_i * r;
    double angle_division = angle / N_arc;
    t[static_cast<std::size_t>(i + 1)] =
        static_cast<long long>(std::floor(r_i * angle_division / d)) + 1;
    N_point = N_point + t[static_cast<std::size_t>(i + 1)] * N_arc + 1;
  }

  // 1-based point storage (dummy row 0), real points 1..N_point
  std::vector<Vec3> P(static_cast<std::size_t>(N_point) + 1, Vec3{});
  int Np = 0;

  // ID-fusion tags aligned with P[1:]; tagged==false -> no tagging. A vertex on
  // both the atom rim and a probe/seam meridian carries the FULL Python tag list.
  std::vector<TagList> vids;
  vids.reserve(static_cast<std::size_t>(N_point));  // capacity hint only
  bool full_circle = angle >= TWO_PI - 1e-9;

  for (int i = 0; i <= N_probe; ++i) {
    long long div_i = t[static_cast<std::size_t>(i + 1)] * N_arc;
    for (long long j = 0; j <= div_i; ++j) {
      ++Np;
      P[static_cast<std::size_t>(Np)] =
          coord_cusp(i, static_cast<int>(j), N_probe, N_arc, P_k1, A1, r, A, Rp,
                     theta1, angle, direct, n, t);
      if (tagged) {
        // Collect tags in Python order (the FULL list per vertex).
        TagList tg;
        if (i == N_probe) tg.push_back(Tag{1, atom_rim, 0});  // ("atom", atom_rim)
        if (full_circle) {
          if (j == 0 || j == div_i) tg.push_back(Tag{3, k, i});  // ("tseam", k, i)
        } else {
          if (e1 != 0 && j == 0) tg.push_back(Tag{2, e1, 0});       // ("probe", e1)
          if (e2 != 0 && j == div_i) tg.push_back(Tag{2, e2, 0});   // ("probe", e2)
        }
        vids.push_back(std::move(tg));
      }
    }
  }

  // Per-row vertex-index offsets: row_off[i] = sum_{ii<i}(t[ii+1]*N_arc+1).
  std::vector<long long> row_off(static_cast<std::size_t>(N_probe) + 1, 0);
  for (int i = 1; i <= N_probe; ++i) {
    row_off[static_cast<std::size_t>(i)] =
        row_off[static_cast<std::size_t>(i - 1)] +
        t[static_cast<std::size_t>(i)] * N_arc + 1;
  }

  std::vector<std::array<int, 3>> T;
  {
    // Exact face count: rows i < N_probe emit 2 triangles per j, the last row 1.
    long long nt = 0;
    for (int i = 1; i <= N_probe; ++i)
      nt += (i < N_probe ? 2 : 1) * t[static_cast<std::size_t>(i + 1)] * N_arc;
    T.reserve(static_cast<std::size_t>(nt));
  }
  for (int i = 1; i <= N_probe; ++i) {
    long long tip1 = t[static_cast<std::size_t>(i + 1)];
    long long ti = t[static_cast<std::size_t>(i)];
    for (long long j = 1; j <= tip1 * N_arc; ++j) {
      if (i < N_probe) {
        // t[i+2] is valid (and Python only reads it) when i < N_probe; reading it
        // at i==N_probe would be one past the end of t (size N_probe+2).
        long long tip2 = t[static_cast<std::size_t>(i + 2)];
        if (tip1 == ti + 1) {
          int j1 = static_cast<int>(j - (j - 1) / tip1 - 1);
          T.push_back(index_P(i - 1, j1, i, static_cast<int>(j - 1), i,
                              static_cast<int>(j), row_off));
        } else {
          T.push_back(index_P(i - 1, static_cast<int>(j), i,
                              static_cast<int>(j - 1), i, static_cast<int>(j),
                              row_off));
        }

        if (tip1 == tip2 - 1) {
          int j1 = static_cast<int>(j + (j - 1) / tip1);
          T.push_back(index_P(i + 1, j1, i, static_cast<int>(j - 1), i,
                              static_cast<int>(j), row_off));
        } else {
          T.push_back(index_P(i + 1, static_cast<int>(j - 1), i,
                              static_cast<int>(j - 1), i, static_cast<int>(j),
                              row_off));
        }
      } else {
        if (tip1 == ti + 1) {
          int j1 = static_cast<int>(j - (j - 1) / tip1 - 1);
          T.push_back(index_P(i - 1, j1, i, static_cast<int>(j - 1), i,
                              static_cast<int>(j), row_off));
        } else {
          T.push_back(index_P(i - 1, static_cast<int>(j), i,
                              static_cast<int>(j - 1), i, static_cast<int>(j),
                              row_off));
        }
      }
    }
  }

  // analytic torus normal at each face centroid, fused with the face orientation
  // (identical per-face expressions; the intermediate array is gone)
  std::vector<Vec3> NV;
  NV.reserve(T.size());
  for (const auto& tri : T) {
    Vec3 p = (P[static_cast<std::size_t>(tri[0])] +
              P[static_cast<std::size_t>(tri[1])] +
              P[static_cast<std::size_t>(tri[2])]) /
             3.0;
    Vec3 nt = (p - A) - (dot(p - A, n) / norm(n)) * n / norm(n);
    Vec3 xp = A + r * nt / norm(nt);
    NV.push_back(orient_face_toroidal(tri, P, (xp - p) / norm(xp - p)));
  }
  // The whole cusp sheet lies on its single rim atom: attribute every vertex to
  // atom_rim. Aligned with P[1..Np]; never affects V/F/N.
  std::vector<int32_t> vatom(static_cast<std::size_t>(Np), atom_rim);
  // RETURN the local mesh; this path was unconditional, so emit is always true.
  LocalMesh out;
  out.P = std::move(P);
  out.T = std::move(T);
  out.NV = std::move(NV);
  out.vids = std::move(vids);
  out.vatom = std::move(vatom);
  out.emit = true;
  return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// driver
// ---------------------------------------------------------------------------
void data_SEStorpat(MeshState& state, const Geom& geom, const DataI& di,
                    const DataSeg& ds, const DataCir& dc, const Ext* ext,
                    double Rp, double d) {
  (void)ext;  // both ext/int branches mesh -> result identical; ext unused.

  const std::vector<Vec3>& I = di.I;
  const auto& segment = ds.segment;
  const int nsegment = ds.nsegment;
  const auto& ncrasegment = ds.ncrasegment;
  const auto& circle = dc.circle;
  const int ncircle = dc.ncircle;

  const std::vector<Vec3>& C = geom.centers;
  const std::vector<double>& R = geom.R;

  // PARALLEL S9: each entity (segment / circle) meshes independently into a
  // fixed 2-slot layout (a regular toroide fills slot 0 only; a singular cusp
  // pair fills slots 0 and 1 in their original order). All inputs are immutable
  // and the mesh_* functions allocate only call-local state, so the per-entity
  // build is race-free. A SERIAL ordered merge then add_patch every emitted
  // LocalMesh in the EXACT original order (segments ascending, then circles
  // ascending; within a singular entity, the two cusp meshes in order).
  std::vector<std::array<LocalMesh, 2>> seg_lm(static_cast<std::size_t>(nsegment) + 1);
  std::vector<std::array<LocalMesh, 2>> crc_lm(static_cast<std::size_t>(ncircle) + 1);

  // --- Arc case ----------------------------------------------------------
  meshms::parallel_for(1, nsegment + 1, [&](int i) {
    const auto& seg = segment[static_cast<std::size_t>(i)];
    const auto& nc = ncrasegment[static_cast<std::size_t>(i)];
    double r = nc[6];  // radius of the i-th arc
    Vec3 A{nc[3], nc[4], nc[5]};
    Vec3 n{nc[0], nc[1], nc[2]};  // from ci to cj
    Vec3 P_k1 = I[static_cast<std::size_t>(seg[2])];
    Vec3 P_k2 = I[static_cast<std::size_t>(seg[3])];
    int direct = seg[4];
    double angle = nc[7];
    double ri = R[static_cast<std::size_t>(seg[0])];
    double rj = R[static_cast<std::size_t>(seg[1])];
    Vec3 ci = C[static_cast<std::size_t>(seg[0])];
    Vec3 cj = C[static_cast<std::size_t>(seg[1])];

    auto& slots = seg_lm[static_cast<std::size_t>(i)];
    // singular self-intersecting branch: r < Rp and dot(ci-A, cj-A) < 0
    if (r < Rp && dot(ci - A, cj - A) < 0.0) {
      Vec3 A1 = A - sqrt_pos(pysq(Rp) - pysq(r)) * n;
      double theta1 = acos_clamp(r / (ri + Rp)) - acos_clamp(r / Rp);
      Vec3 c1 = ci * Rp / (Rp + ri) + A * ri / (Rp + ri);
      double r1 = ri / (Rp + ri) * r;
      slots[0] = mesh_cusp(c1, r1, A1, direct, n, angle, r, A, P_k1, P_k2,
                           theta1, Rp, d, i, /*atom_rim=*/seg[0],
                           /*e1=*/seg[2], /*e2=*/seg[3], /*tagged=*/true);

      Vec3 A2 = A + sqrt_pos(pysq(Rp) - pysq(r)) * n;
      double theta2 = acos_clamp(r / (rj + Rp)) - acos_clamp(r / Rp);
      c1 = cj * Rp / (Rp + rj) + A * rj / (Rp + rj);
      r1 = rj / (Rp + rj) * r;
      slots[1] = mesh_cusp(c1, r1, A2, direct, n, angle, r, A, P_k1, P_k2, theta2,
                           Rp, d, i, /*atom_rim=*/seg[1], /*e1=*/seg[2],
                           /*e2=*/seg[3], /*tagged=*/true);
    } else {
      double theta1 = acos_clamp(r / (ri + Rp));
      double theta2 = acos_clamp(r / (rj + Rp));
      slots[0] = mesh_toroide(ci, cj, ri, rj, direct, n, angle, r, A, P_k1,
                              P_k2, theta1, theta2, Rp, d, i,
                              /*atom_i=*/seg[0], /*atom_j=*/seg[1],
                              /*e1=*/seg[2], /*e2=*/seg[3]);
    }
  });

  // --- Circle case -------------------------------------------------------
  meshms::parallel_for(1, ncircle + 1, [&](int i) {
    const auto& crc = circle[static_cast<std::size_t>(i)];
    double r = crc[9];  // radius of the i-th circle
    Vec3 A{crc[3], crc[4], crc[5]};
    Vec3 n{crc[6], crc[7], crc[8]};
    int direct = 1;
    auto [v1, v2] = orthogonalvectors(n);
    (void)v2;
    Vec3 P_k1 = A + r * v1;

    int ci_idx = static_cast<int>(crc[1]);
    int cj_idx = static_cast<int>(crc[2]);
    double ri = R[static_cast<std::size_t>(ci_idx)];
    double rj = R[static_cast<std::size_t>(cj_idx)];
    Vec3 ci = C[static_cast<std::size_t>(ci_idx)];
    Vec3 cj = C[static_cast<std::size_t>(cj_idx)];

    auto& slots = crc_lm[static_cast<std::size_t>(i)];
    if (r < Rp && dot(ci - A, cj - A) < 0.0) {
      // GEOMETRY FIX: each apex picked by the side its rim atom lies on.
      double h = sqrt_pos(pysq(Rp) - pysq(r));
      Vec3 A_i = A + std::copysign(h, dot(ci - A, n)) * n;
      Vec3 A_j = A + std::copysign(h, dot(cj - A, n)) * n;

      double theta1 = acos_clamp(r / (ri + Rp)) - acos_clamp(r / Rp);
      double theta2 = acos_clamp(r / (rj + Rp)) - acos_clamp(r / Rp);

      Vec3 c1 = ci * Rp / (Rp + ri) + A * ri / (Rp + ri);
      double r1 = ri / (Rp + ri) * r;
      slots[0] = mesh_cusp(c1, r1, A_i, direct, n, TWO_PI, r, A, P_k1, P_k1,
                           theta1, Rp, d, -i, /*atom_rim=*/ci_idx, /*e1=*/0,
                           /*e2=*/0, /*tagged=*/true);

      c1 = cj * Rp / (Rp + rj) + A * rj / (Rp + rj);
      r1 = rj / (Rp + rj) * r;
      slots[1] = mesh_cusp(c1, r1, A_j, direct, n, TWO_PI, r, A, P_k1, P_k1,
                           theta2, Rp, d, -i, /*atom_rim=*/cj_idx, /*e1=*/0,
                           /*e2=*/0, /*tagged=*/true);
    } else {
      double theta1 = acos_clamp(r / (ri + Rp));
      double theta2 = acos_clamp(r / (rj + Rp));
      slots[0] = mesh_toroide(ci, cj, ri, rj, direct, n, TWO_PI, r, A, P_k1,
                              P_k1, theta1, theta2, Rp, d, -i,
                              /*atom_i=*/ci_idx, /*atom_j=*/cj_idx, /*e1=*/0,
                              /*e2=*/0);
    }
  });

  // --- SERIAL ordered merge (segments ascending, then circles ascending) -
  // Pre-sum the emitted sizes so the accumulator reserves once (capacity-only),
  // and move the per-vertex tag lists out of the dead LocalMesh.
  std::size_t add_v = 0, add_f = 0;
  auto count_lm = [&](const std::array<LocalMesh, 2>& pair) {
    for (const LocalMesh& lm : pair) {
      if (lm.emit) {
        add_v += lm.P.empty() ? 0 : lm.P.size() - 1;
        add_f += lm.T.size();
      }
    }
  };
  for (int i = 1; i <= nsegment; ++i) count_lm(seg_lm[static_cast<std::size_t>(i)]);
  for (int i = 1; i <= ncircle; ++i) count_lm(crc_lm[static_cast<std::size_t>(i)]);
  state.reserve_extra(add_v, add_f);
  for (int i = 1; i <= nsegment; ++i) {
    for (LocalMesh& lm : seg_lm[static_cast<std::size_t>(i)]) {
      if (lm.emit) state.add_patch(lm.P, lm.T, lm.NV, std::move(lm.vids), lm.vatom);
    }
  }
  for (int i = 1; i <= ncircle; ++i) {
    for (LocalMesh& lm : crc_lm[static_cast<std::size_t>(i)]) {
      if (lm.emit) state.add_patch(lm.P, lm.T, lm.NV, std::move(lm.vids), lm.vatom);
    }
  }
}

}  // namespace meshms
