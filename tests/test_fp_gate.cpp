// Equivalence gate: the end-to-end contract a relaxed-FP deploy build must meet.
//
// MESHMS_FP=fast trades the bit-for-bit golden gate for speed, so the nine
// golden bit-regression tests skip themselves there (MESHMS_SKIP_IF_FAST). This
// test is what replaces them: it asserts the properties a consumer actually
// depends on -- the build does not throw, the mesh contains no NaN, the vertices
// stay inside the molecule's bounding box, the aggregate geometry matches a
// frozen strict-build baseline, and close_cusps still yields a watertight mesh.
//
// The SAME baseline serves both policies at different tolerances: tight in a
// strict build (where it is an ordinary regression lock) and loose in a fast one
// (where it is a "did the surface break" envelope). Regenerate it from a strict
// build with:
//
//   ./build/test_fp_gate --dump > tests/ref/fp_gate.txt
//
// Everything here goes through <meshms/capi.hpp> only -- the consumer's view --
// so the facade paths (multi-component split, isolated atoms, post-processing)
// are exercised alongside the pipeline.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <string>
#include <vector>

#include "meshms/capi.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

struct Case {
  const char* mol;
  double d;
};

// Small shapes first (the cusp cases are the degenerate-geometry stressors),
// then the proteins none of the other tests cover, then a density sweep on 101M
// where the advancing front runs longest.
//
// fullerene.xyzr is deliberately NOT here: a strict build already produces a
// broken mesh for it (about half the vertices come out NaN), which is the
// deferred symmetry-jitter degenerate fallback noted in README "Out of scope",
// not anything this gate introduces. Add it back once that is implemented -- it
// is exactly the kind of case this gate should cover.
const Case kCases[] = {
    {"tetra", 0.5},   {"cusp2", 0.5},   {"cusp3", 0.5}, {"ArgArg", 0.5},
    {"1YJO", 0.5},    {"1ETN", 0.5},    {"1crn", 0.5},  {"1B17", 0.5},
    {"barstar", 0.5}, {"101M", 0.5},    {"101M", 0.25},
};

// One measured (or frozen) row. Aggregates only: a fast build may legitimately
// land on a different vertex count, so per-vertex comparison is not the gate.
struct Row {
  std::string mol;
  double d = 0.0;
  long nV = 0, nF = 0;
  long nonfinite = 0, degenerate = 0, duplicate = 0;
  long boundary = 0, nonmanifold = 0;
  int watertight = 0;
  double area = 0.0, volume = 0.0;
  int closed_watertight = 0;
  double closed_area = 0.0, closed_volume = 0.0;
};

std::string data_path(const std::string& mol) {
  return std::string(MESHMS_DATA_DIR) + "/" + mol + ".xyzr";
}

// Every vertex must sit inside the atom bounding box grown by (max radius +
// probe + slack). Catches an exploded or sign-flipped vertex without needing a
// baseline or a reference build, so it holds on every platform.
long count_outside_bbox(const std::vector<std::array<double, 4>>& atoms,
                        const MeshResult& m, double probe) {
  if (atoms.empty() || m.verts.empty()) return 0;
  double lo[3] = {atoms[0][0], atoms[0][1], atoms[0][2]};
  double hi[3] = {atoms[0][0], atoms[0][1], atoms[0][2]};
  double rmax = 0.0;
  for (const auto& a : atoms) {
    for (int k = 0; k < 3; ++k) {
      if (a[k] < lo[k]) lo[k] = a[k];
      if (a[k] > hi[k]) hi[k] = a[k];
    }
    if (a[3] > rmax) rmax = a[3];
  }
  const double slack = rmax + probe + 1.0;
  long bad = 0;
  for (const auto& v : m.verts) {
    for (int k = 0; k < 3; ++k) {
      if (v[k] < lo[k] - slack || v[k] > hi[k] + slack) {
        ++bad;
        break;
      }
    }
  }
  return bad;
}

// Longest deviation of a vertex normal from unit length. NaN-poisoned normals
// collapse to 0 or NaN, so this catches them even when the positions survive.
double worst_normal_error(const MeshResult& m) {
  double worst = 0.0;
  for (const auto& n : m.vnormals) {
    const double len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    const double e = std::fabs(len - 1.0);
    if (!(e <= worst)) worst = e;  // NaN-safe: propagates a NaN into `worst`
  }
  return worst;
}

const double kProbe = 1.4;

// Build + post-process one case and collect the aggregates. Throws only if the
// pipeline does, which the caller reports as a gate failure.
Row measure(const Case& c) {
  Row r;
  r.mol = c.mol;
  r.d = c.d;
  const std::vector<std::array<double, 4>> atoms = load_xyzr_array(data_path(c.mol));
  CHECK(!atoms.empty());

  const MeshResult fused = build_surface_from_array(atoms, kProbe, c.d, /*fuse=*/true);
  const MeshReport rep = analyze_mesh(fused);
  r.nV = rep.n_vertices;
  r.nF = rep.n_faces;
  r.nonfinite = rep.nonfinite_vertices;
  r.degenerate = rep.degenerate_faces;
  r.duplicate = rep.duplicate_faces;
  r.boundary = rep.boundary_edges;
  r.nonmanifold = rep.nonmanifold_edges;
  r.watertight = rep.watertight ? 1 : 0;
  r.area = rep.area;
  r.volume = rep.signed_volume;

  const MeshResult closed = close_cusps(fused);
  const MeshReport rc = analyze_mesh(closed);
  r.closed_watertight = rc.watertight ? 1 : 0;
  r.closed_area = rc.area;
  r.closed_volume = rc.signed_volume;

  // Property checks that need no baseline at all.
  CHECK(rep.nonfinite_vertices == 0);
  CHECK(rc.nonfinite_vertices == 0);
  CHECK(count_outside_bbox(atoms, fused, kProbe) == 0);
  CHECK(worst_normal_error(fused) < 1e-6);
  CHECK(r.nV > 0 && r.nF > 0);
  CHECK(r.area > 0.0);
  CHECK(r.volume > 0.0);  // orient_faces produced outward-facing triangles
  // remove_flaps never adds faces and never makes the mesh less manifold. It is
  // a no-op on a clean mesh, but barstar genuinely carries flaps, so "unchanged"
  // is the wrong assertion -- "not worse" is the invariant.
  const MeshResult flapped = remove_flaps(fused);
  const MeshReport rflap = analyze_mesh(flapped);
  CHECK(flapped.faces.size() <= fused.faces.size());
  CHECK(flapped.face_type.size() == flapped.faces.size());
  CHECK(rflap.nonmanifold_edges <= rep.nonmanifold_edges);
  CHECK(rflap.duplicate_faces <= rep.duplicate_faces);
  // Same inputs, same process, twice: catches nondeterminism the relaxed flags
  // could introduce via reassociation or thread scheduling.
  const MeshResult again = build_surface_from_array(atoms, kProbe, c.d, /*fuse=*/true);
  CHECK(again.verts == fused.verts);
  CHECK(again.faces == fused.faces);
  return r;
}

void print_row(const Row& r) {
  std::printf("%s %g %ld %ld %ld %ld %ld %ld %ld %d %.17g %.17g %d %.17g %.17g\n",
              r.mol.c_str(), r.d, r.nV, r.nF, r.nonfinite, r.degenerate, r.duplicate,
              r.boundary, r.nonmanifold, r.watertight, r.area, r.volume,
              r.closed_watertight, r.closed_area, r.closed_volume);
}

bool read_baseline(std::vector<Row>& out) {
  const std::string path = std::string(MESHMS_REF_DIR) + "/fp_gate.txt";
  std::FILE* f = std::fopen(path.c_str(), "r");
  if (!f) {
    std::fprintf(stderr, "cannot open baseline %s\n", path.c_str());
    return false;
  }
  char mol[64];
  Row r;
  while (std::fscanf(f, "%63s %lf %ld %ld %ld %ld %ld %ld %ld %d %lf %lf %d %lf %lf",
                     mol, &r.d, &r.nV, &r.nF, &r.nonfinite, &r.degenerate, &r.duplicate,
                     &r.boundary, &r.nonmanifold, &r.watertight, &r.area, &r.volume,
                     &r.closed_watertight, &r.closed_area, &r.closed_volume) == 15) {
    r.mol = mol;
    out.push_back(r);
  }
  std::fclose(f);
  return !out.empty();
}

// Tolerances. Strict keeps the baseline as an exact-ish regression lock; fast
// widens it to the gap between "the last bits moved" and "the surface broke".
// A discretised area/volume does not shift by 0.1% from a retriangulation, but
// a dropped or exploded patch shifts it by whole percent.
#ifdef MESHMS_FP_FAST
const double kCountRel = 0.02;    // vertex/face count drift
const double kGeomRel = 1e-3;     // area / signed volume
const long kBoundarySlack = 64;   // absolute slack on open boundary edges
const double kBoundaryRel = 0.5;
#else
const double kCountRel = 0.001;
const double kGeomRel = 1e-6;
const long kBoundarySlack = 0;
const double kBoundaryRel = 0.001;
#endif

void compare(const Row& got, const Row& want) {
  std::printf("  %-9s d=%.2f nV=%ld (base %ld) area=%.6g (base %.6g)\n", got.mol.c_str(),
              got.d, got.nV, want.nV, got.area, want.area);
  CHECK_REL(static_cast<double>(got.nV), static_cast<double>(want.nV), kCountRel);
  CHECK_REL(static_cast<double>(got.nF), static_cast<double>(want.nF), kCountRel);
  // Degenerate faces may not increase beyond a per-mille of the mesh; duplicate
  // faces and non-manifold edges may not increase at all -- those are real
  // breakage, not a retriangulation artefact.
  CHECK(got.degenerate <= want.degenerate + (kFpFast ? got.nF / 1000 : 0));
  CHECK(got.duplicate <= want.duplicate);
  CHECK(got.nonmanifold <= want.nonmanifold);
  const long bmax = want.boundary + kBoundarySlack +
                    static_cast<long>(kBoundaryRel * static_cast<double>(want.boundary));
  CHECK(got.boundary <= bmax);
  CHECK(got.watertight == want.watertight);
  CHECK_REL(got.area, want.area, kGeomRel);
  CHECK_REL(got.volume, want.volume, kGeomRel);
  // The closed surface is what a consumer renders, so it is held to the same
  // standard in both policies: if the baseline closed watertight, so must this.
  CHECK(got.closed_watertight >= want.closed_watertight);
  CHECK_REL(got.closed_area, want.closed_area, kGeomRel);
  CHECK_REL(got.closed_volume, want.closed_volume, kGeomRel);
}

// Malformed input is rejected at the facade boundary instead of propagating
// through the pipeline and surfacing as a NaN mesh. Radius 0 stays valid (the
// xyzr convention uses it for non-contributing atoms) and is covered by the
// barstar case above, which has 503 of them.
void check_input_validation() {
  const std::vector<std::array<double, 4>> bad = {{0.0, 0.0, 0.0, 1.5},
                                                  {1.0, std::nan(""), 0.0, 1.5}};
  bool threw = false;
  try {
    build_surface_from_array(bad, kProbe, 0.5, /*fuse=*/true);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
}

}  // namespace

int main(int argc, char** argv) {
  const bool dump = (argc > 1 && std::strcmp(argv[1], "--dump") == 0);
  if (!dump) {
    std::printf("fp gate: %s (fp_mode=%d)\n", build_info(), fp_mode());
    // The library and these headers must agree on the policy, or the inline math
    // in the internal headers was compiled differently than the library was.
    CHECK(fp_mode() == (kFpFast ? 1 : 0));
  }

  std::vector<Row> baseline;
  if (!dump && !read_baseline(baseline)) return 1;

  for (const Case& c : kCases) {
    Row got;
    try {
      got = measure(c);
    } catch (const std::exception& e) {
      std::fprintf(stderr, "FAIL %s d=%g threw: %s\n", c.mol, c.d, e.what());
      ++g_fail;
      continue;
    }
    if (dump) {
      print_row(got);
      continue;
    }
    const Row* want = nullptr;
    for (const Row& b : baseline) {
      if (b.mol == got.mol && std::fabs(b.d - got.d) < 1e-12) want = &b;
    }
    if (!want) {
      std::fprintf(stderr, "FAIL no baseline row for %s d=%g\n", c.mol, c.d);
      ++g_fail;
      continue;
    }
    compare(got, *want);
  }

  if (dump) return 0;
  check_input_validation();
  TEST_MAIN_RETURN();
}
