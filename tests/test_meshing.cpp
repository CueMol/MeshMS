// Faithful-port cross-check of mesh_sphpat against the per-call golden
// <mol>.meshcalls.txt (every mesh_sphpat call during a Python build of
// {2spheres, 3spheres, tetra}, captured with radius_probe=1.4, mesh_size=0.5,
// fuse=False).
//
// Each block records the full mesh_sphpat inputs (csphere/rsphere/Rp/d, the
// patches, the loops it references, segment0, circle0, Rj, boundary tag) and the
// vertices/triangles/per-face-normals it appended to the MeshState. The port must
// reproduce: appended vertices == golden P (1e-9), faces == golden T (EXACT after
// the 1-based<->0-based rebase add_patch performs), per-face normals == golden NV
// (1e-9), and Np/Nt EXACT. A topology mismatch means a float term-order / pysq /
// comparison-tie deviation flipped an advancing-front decision.
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "meshms/mesh_state.hpp"
#include "meshms/meshing.hpp"
#include "meshms/vec3.hpp"
#include "test_util.hpp"

using namespace meshms;

namespace {

struct Call {
  int id = 0;
  Vec3 csphere{};
  double rsphere = 0.0;
  double Rp = 0.0;
  double d = 0.0;
  int patchesize = 0;
  std::vector<int> patches;                       // 1-based, [0] dummy
  std::vector<Loop> loops;                         // 1-based by loop id, [0] dummy; each loop 1-based
  std::vector<std::array<double, 12>> segment0;    // [0] dummy row
  std::vector<std::array<double, 9>> circle0;      // [0] dummy row
  int Rj_flag = 0;
  std::vector<double> Rj;                           // 1-based by local seg index, [0] dummy
  Tag btag{};
  // golden output
  int outNp = 0;
  std::vector<Vec3> outP;                           // outP[k] for k=1..outNp (0-based store: outP[k-1])
  int outNt = 0;
  std::vector<std::array<int, 3>> outT;            // 1-based into this call's P
  std::vector<Vec3> outNV;
};

std::vector<Call> read_calls(const std::string& path) {
  std::ifstream in(path);
  CHECK(static_cast<bool>(in));
  std::vector<Call> calls;
  Call cur;
  bool in_out = false;
  int max_loop_id = 0;
  std::string line;

  auto flush = [&]() {
    if (cur.id != 0) calls.push_back(cur);
  };

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream ls(line);
    std::string tag;
    ls >> tag;
    if (tag == "CALL") {
      flush();
      cur = Call{};
      max_loop_id = 0;
      in_out = false;
      ls >> cur.id;
    } else if (tag == "csphere") {
      ls >> cur.csphere.x >> cur.csphere.y >> cur.csphere.z;
    } else if (tag == "rsphere") {
      ls >> cur.rsphere;
    } else if (tag == "Rp") {
      ls >> cur.Rp;
    } else if (tag == "d") {
      ls >> cur.d;
    } else if (tag == "patchesize") {
      ls >> cur.patchesize;
      cur.patches.assign(static_cast<std::size_t>(cur.patchesize) + 1, 0);
    } else if (tag == "patches") {
      // signed entries 1..patchesize
      for (int i = 1; i <= cur.patchesize; ++i) ls >> cur.patches[static_cast<std::size_t>(i)];
    } else if (tag == "nloops") {
      int nl = 0;
      ls >> nl;
      (void)nl;  // loops are addressed by id; size lazily on LOOP
    } else if (tag == "LOOP") {
      int kk = 0, lsz = 0;
      ls >> kk >> lsz;
      if (kk > max_loop_id) {
        cur.loops.resize(static_cast<std::size_t>(kk) + 1);
        max_loop_id = kk;
      }
      Loop lp(1);  // [0] dummy
      for (int t = 0; t < lsz; ++t) {
        int s = 0;
        ls >> s;
        lp.push_back(s);
      }
      cur.loops[static_cast<std::size_t>(kk)] = std::move(lp);
    } else if (tag == "nseg0") {
      int ns = 0;
      ls >> ns;
      cur.segment0.assign(static_cast<std::size_t>(ns) + 1, std::array<double, 12>{});
    } else if (tag == "SEG0") {
      int s = 0;
      ls >> s;
      std::array<double, 12>& row = cur.segment0[static_cast<std::size_t>(s)];
      for (int c = 1; c <= 11; ++c) ls >> row[static_cast<std::size_t>(c)];
    } else if (tag == "ncirc0") {
      int nc = 0;
      ls >> nc;
      cur.circle0.assign(static_cast<std::size_t>(nc) + 1, std::array<double, 9>{});
    } else if (tag == "CIRC0") {
      int kk = 0;
      ls >> kk;
      std::array<double, 9>& row = cur.circle0[static_cast<std::size_t>(kk)];
      // c(3), n(3), r -> cols 1..7; col 8 (torus radius) absent here -> 0.
      for (int c = 1; c <= 7; ++c) ls >> row[static_cast<std::size_t>(c)];
      row[8] = 0.0;
    } else if (tag == "Rj") {
      ls >> cur.Rj_flag;
    } else if (tag == "Rjvals") {
      // aligned with segment0: Rj[1..S]. Read all remaining as 1..S.
      cur.Rj.assign(cur.segment0.size(), 0.0);  // size S+1, [0] dummy
      int idx = 1;
      double v;
      while (ls >> v) {
        if (static_cast<std::size_t>(idx) < cur.Rj.size())
          cur.Rj[static_cast<std::size_t>(idx)] = v;
        ++idx;
      }
    } else if (tag == "btag") {
      ls >> cur.btag.kind >> cur.btag.i >> cur.btag.j;
    } else if (tag == "--") {
      in_out = true;  // "-- out --"
    } else if (tag == "outNp") {
      ls >> cur.outNp;
    } else if (tag == "P" && in_out) {
      int kk = 0;
      Vec3 v{};
      ls >> kk >> v.x >> v.y >> v.z;
      cur.outP.push_back(v);
    } else if (tag == "outNt") {
      ls >> cur.outNt;
    } else if (tag == "T" && in_out) {
      int kk = 0;
      std::array<int, 3> t{};
      ls >> kk >> t[0] >> t[1] >> t[2];
      cur.outT.push_back(t);
    } else if (tag == "NV") {
      int kk = 0;
      Vec3 v{};
      ls >> kk >> v.x >> v.y >> v.z;
      cur.outNV.push_back(v);
    }
  }
  flush();
  return calls;
}

void check_vec(const std::string& mol, int call, const char* what, int idx,
               const Vec3& got, const Vec3& want, double tol) {
  if (!(std::fabs(got.x - want.x) <= tol && std::fabs(got.y - want.y) <= tol &&
        std::fabs(got.z - want.z) <= tol)) {
    std::fprintf(stderr,
                 "FAIL %s CALL %d %s[%d] got=(%.17g,%.17g,%.17g) "
                 "want=(%.17g,%.17g,%.17g)\n",
                 mol.c_str(), call, what, idx, got.x, got.y, got.z, want.x,
                 want.y, want.z);
    ++g_fail;
  }
}

void check_call(const std::string& mol, const Call& c) {
  MeshState st;
  const std::vector<double>* rj = (c.Rj_flag == 1) ? &c.Rj : nullptr;
  LocalMesh lm =
      mesh_sphpat(c.csphere, c.rsphere, c.loops, c.segment0, c.circle0, c.patches,
                  c.patchesize, c.Rp, c.d, rj, c.btag);
  if (lm.emit) st.add_patch(lm.P, lm.T, lm.NV, lm.vids);

  // Vertices appended this call (fresh state, base 0).
  const int gotNp = static_cast<int>(st.V.size());
  CHECK(gotNp == c.outNp);
  if (gotNp != c.outNp) {
    std::fprintf(stderr, "  %s CALL %d Np got=%d want=%d\n", mol.c_str(), c.id,
                 gotNp, c.outNp);
  }
  const int nP = std::min(gotNp, c.outNp);
  for (int i = 0; i < nP; ++i) {
    check_vec(mol, c.id, "P", i + 1, st.V[static_cast<std::size_t>(i)],
              c.outP[static_cast<std::size_t>(i)], 1e-9);
  }

  // Faces: golden T is 1-based into this call's P; add_patch rebased to 0-based
  // (base 0 for a fresh state). So st.F[t][k] == golden T[t][k] - 1.
  const int gotNt = static_cast<int>(st.F.size());
  CHECK(gotNt == c.outNt);
  if (gotNt != c.outNt) {
    std::fprintf(stderr, "  %s CALL %d Nt got=%d want=%d\n", mol.c_str(), c.id,
                 gotNt, c.outNt);
  }
  const int nT = std::min(gotNt, c.outNt);
  for (int t = 0; t < nT; ++t) {
    const Tri& f = st.F[static_cast<std::size_t>(t)];
    const std::array<int, 3>& w = c.outT[static_cast<std::size_t>(t)];
    bool ok = (f[0] == w[0] - 1) && (f[1] == w[1] - 1) && (f[2] == w[2] - 1);
    CHECK(ok);
    if (!ok) {
      std::fprintf(stderr,
                   "  %s CALL %d T[%d] got=(%d,%d,%d) want 0-based=(%d,%d,%d)\n",
                   mol.c_str(), c.id, t + 1, f[0], f[1], f[2], w[0] - 1, w[1] - 1,
                   w[2] - 1);
    }
  }

  // Per-face normals.
  const int gotNV = static_cast<int>(st.N.size());
  CHECK(gotNV == static_cast<int>(c.outNV.size()));
  const int nNV = std::min(gotNV, static_cast<int>(c.outNV.size()));
  for (int t = 0; t < nNV; ++t) {
    check_vec(mol, c.id, "NV", t + 1, st.N[static_cast<std::size_t>(t)],
              c.outNV[static_cast<std::size_t>(t)], 1e-9);
  }
}

void check_mol(const std::string& mol) {
  const std::string path =
      std::string(MESHMS_REF_DIR) + "/" + mol + ".meshcalls.txt";
  std::vector<Call> calls = read_calls(path);
  CHECK(!calls.empty());
  for (const Call& c : calls) check_call(mol, c);
}

}  // namespace

int main() {
  check_mol("2spheres");
  check_mol("3spheres");
  check_mol("tetra");
  TEST_MAIN_RETURN();
}
