// Faithful port of the pipeline module (geometry path).
#include "meshms/pipeline.hpp"

#include <cmath>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "meshms/concave.hpp"
#include "meshms/convex.hpp"
#include "meshms/exterior.hpp"
#include "meshms/fusion.hpp"
#include "meshms/geom.hpp"
#include "meshms/intersection.hpp"
#include "meshms/mesh_state.hpp"
#include "meshms/parallel.hpp"
#include "meshms/sas.hpp"
#include "meshms/sas_patches.hpp"
#include "meshms/toroidal.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

void orient_faces(const std::vector<Vec3>& V, std::vector<Tri>& F,
                  const std::vector<Vec3>& N) {
  // Mirror pipeline.orient_faces: a = V[F[:,0]], b = V[F[:,1]], c = V[F[:,2]];
  // wn = cross(b-a, c-a); flip = (wn*N).sum(axis=1) < 0.0; swap F[flip,1] and
  // F[flip,2]. len(F)==0 short-circuit is naturally a no-op here.
  // Each face is independent (reads V/N const, swaps only F[f] in place), so the
  // loop is parallel; per-face float math and results are unchanged.
  parallel_for(0, static_cast<int>(F.size()), [&](int fi) {
    const std::size_t f = static_cast<std::size_t>(fi);
    const Vec3& a = V[static_cast<std::size_t>(F[f][0])];
    const Vec3& b = V[static_cast<std::size_t>(F[f][1])];
    const Vec3& c = V[static_cast<std::size_t>(F[f][2])];
    const Vec3 wn = cross(b - a, c - a);
    if (dot(wn, N[f]) < 0.0) {
      std::swap(F[f][1], F[f][2]);
    }
  });
}

namespace {

// Per-vertex outward normals (Surface.NV): accumulate each face normal onto its
// three vertices, then normalise. Non-fuse meshes still carry the per-face N, so
// those (already-oriented, unit) normals are summed directly; the fuse path drops
// N, so the face normal is recomputed from the oriented winding as
// normalize(cross(b-a, c-a)). NV is NOT a golden-checked field (the golden tests
// compare V/F/N only), so this never affects bit-equivalence.
void compute_vertex_normals(Surface& s) {
  s.NV.assign(s.V.size(), Vec3{0.0, 0.0, 0.0});
  const bool have_face_n = !s.N.empty();
  for (std::size_t f = 0; f < s.F.size(); ++f) {
    const Tri& t = s.F[f];
    Vec3 fn;
    if (have_face_n) {
      fn = s.N[f];
    } else {
      const Vec3& a = s.V[static_cast<std::size_t>(t[0])];
      const Vec3& b = s.V[static_cast<std::size_t>(t[1])];
      const Vec3& c = s.V[static_cast<std::size_t>(t[2])];
      Vec3 cn = cross(b - a, c - a);
      double cl = norm(cn);
      fn = (cl > 0.0) ? cn / cl : Vec3{0.0, 0.0, 0.0};
    }
    for (int k = 0; k < 3; ++k) {
      Vec3& acc = s.NV[static_cast<std::size_t>(t[k])];
      acc = acc + fn;
    }
  }
  for (Vec3& nv : s.NV) {
    double l = norm(nv);
    if (l > 0.0) nv = nv / l;
  }
}

}  // namespace

RSComponents compute_rs(const Geom& geom, double radius_probe) {
  const double Rp = radius_probe;

  RSComponents rs;
  rs.geom = geom;  // density-independent; kept so build_mesh needs only the cache
  rs.Rp = Rp;

  // --- SAS arrangement (data_I_Cir / data_Seg_Pat / data_ext) ---------------
  // Identical call order/arguments to the old run() L40-44; data_ext reads the
  // local structured-binding results BEFORE they are moved into `rs`, so the
  // computation is byte-identical (only ownership transfers afterwards).
  rs.inter = interstructure(geom, Rp);
  auto [di, dc] = data_I_Cir(geom, rs.inter, Rp);
  auto [ds, dl, dp] = data_Seg_Pat(geom, rs.inter, di, dc, Rp);
  rs.ext = data_ext(geom, rs.inter, di, dc, ds, dl, dp, Rp);
  rs.di = std::move(di);
  rs.dc = std::move(dc);
  rs.ds = std::move(ds);
  rs.dl = std::move(dl);
  rs.dp = std::move(dp);
  // Density-independent concave decomposition (was recomputed inside every
  // SESconcavepat call): pure memoization -- same function, same inputs,
  // computed once. build_mesh() then only re-meshes it per density.
  rs.concave = precompute_concave(rs.geom, rs.di, Rp, rs.inter);
  return rs;
}

Surface build_mesh(const RSComponents& rs, double mesh_size, bool fuse) {
  const Geom& geom = rs.geom;
  const double Rp = rs.Rp;
  const double d = mesh_size;

  // --- SES structure: mesh into the shared accumulator (convex, then concave,
  // then toroidal -- the exact Python order so indices line up) --------------
  // Snapshot the face count after each mesher: the three appended blocks are the
  // SES component types (MSMS codes 3 convex / 2 concave / 1 toroidal).
  MeshState state;
  data_SESsphpat_convex(state, geom, rs.di, rs.dc, rs.ds, rs.dl, rs.dp, &rs.ext, Rp, d);
  const std::size_t n_convex = state.F.size();
  SESconcavepat_mesh(state, geom, rs.di, rs.concave, Rp, d);
  const std::size_t n_concave = state.F.size();
  data_SEStorpat(state, geom, rs.di, rs.ds, rs.dc, &rs.ext, Rp, d);
  const std::size_t n_total = state.F.size();

  Surface out;
  // state is dead after this point except state.tags (read by the fuse branch),
  // so the four big arrays can be moved out instead of copied.
  out.V = std::move(state.V);
  out.F = std::move(state.F);
  out.N = std::move(state.N);
  out.atom_id = std::move(state.vatom);  // per-vertex owning atom (aligned with V)
  out.ftype.resize(n_total);
  for (std::size_t f = 0; f < n_total; ++f)
    out.ftype[f] = static_cast<std::uint8_t>(f < n_convex  ? 3   // convex / contact
                                             : f < n_concave ? 2  // concave / reentrant
                                                             : 1);  // toroidal
  orient_faces(out.V, out.F, out.N);  // in place per face: never reorders F

  if (fuse) {
    // ID-based boundary fusion (topological tags, not coordinate quantization):
    // structurally connects the non-cusp patches. The face count may drop
    // (degenerate triangles), so the per-face normals N are no longer aligned ->
    // N is dropped, mirroring Python's `return V, F, None`. warn=False path. The
    // per-vertex atom_id is remapped to the surviving union-find root's atom id;
    // ftype is filtered to the surviving faces via kept_faces (same drop).
    std::vector<std::uint32_t> kept_faces;
    auto [V2, F2, A2] = fuse_by_id(out.V, out.F, state.tags, out.atom_id, 1e-6, &kept_faces);
    out.V = std::move(V2);
    out.F = std::move(F2);
    out.atom_id = std::move(A2);
    out.N.clear();
    std::vector<std::uint8_t> ftype2;
    ftype2.reserve(kept_faces.size());
    for (std::uint32_t k : kept_faces) ftype2.push_back(out.ftype[k]);
    out.ftype = std::move(ftype2);
  }

  // Per-vertex normals for external consumers (cuemol2/BALL expect per-vertex).
  compute_vertex_normals(out);
  return out;
}

Surface run(const Geom& geom, const Para& para, bool fuse) {
  // Split into the density-independent RS components and the density-dependent
  // mesher; the only added cost vs the old monolithic run() is one Geom copy.
  RSComponents rs = compute_rs(geom, para.radius_probe);
  return build_mesh(rs, para.mesh_size, fuse);
}

Surface build_surface(const std::string& xyzr_path, double radius_probe,
                      double mesh_size, bool fuse) {
  Geom geom = read_xyzr(xyzr_path);
  Para para;
  para.radius_probe = radius_probe;
  para.mesh_size = mesh_size;

  // jitter=None semantics: build from faithful coordinates with no perturbation.
  // The jitter="auto" NaN-retry fallback lives in run_auto (call it on the geom
  // read from the file when faithful coordinates yield a degenerate mesh); this
  // file-path helper deliberately stays faithful-only.
  return run(geom, para, fuse);
}

bool surface_has_nan(const Surface& s) {
  // A single non-finite vertex component poisons the whole surface (areas /
  // normals propagate the NaN), so any non-finite component flags the mesh.
  for (const Vec3& v : s.V) {
    if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z))
      return true;
  }
  return false;
}

Surface run_auto(const Geom& geom, const Para& para, bool fuse,
                 double jitter_mag, int max_retries) {
  // Deterministic (non-numpy-matching) degenerate-symmetry fallback. Attempt 0
  // is the faithful run() so clean molecules stay bit-identical; only a NaN
  // result triggers the jittered retries. Build failures are caught so a single
  // unstable perturbation cannot abort the whole search.
  Surface best;
  bool have_best = false;
  std::exception_ptr last_error;

  for (int attempt = 0; attempt <= max_retries; ++attempt) {
    try {
      // attempt 0 == faithful centers; attempt k>=1 perturbs with seed = k.
      const Geom& g = (attempt == 0) ? geom : jitter_centers(geom, jitter_mag, attempt);
      // Re-run the full geom -> compute_rs -> build_mesh cycle: jittering the
      // centers invalidates the RS arrangement, so build_mesh_from_cache reuse
      // is not applicable here.
      Surface s = run(g, para, fuse);
      if (!surface_has_nan(s)) return s;  // first NaN-free mesh wins
      best = std::move(s);                // keep the last finite-sized attempt
      have_best = true;
      last_error = nullptr;
    } catch (...) {
      last_error = std::current_exception();  // failed build -> keep trying
    }
  }

  // Every attempt was degenerate: return the last successful (still-NaN) build
  // as a best effort, or rethrow if even the last attempt threw.
  if (have_best) return best;
  if (last_error) std::rethrow_exception(last_error);
  return best;  // unreachable for max_retries >= 0, but keeps the return total
}

}  // namespace meshms
