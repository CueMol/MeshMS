#pragma once
// Minimal header-only check harness for the C++ port tests (no GoogleTest dep).
// One test_*.cpp == one executable; CHECK macros count failures into g_fail and
// TEST_MAIN_RETURN() makes main() exit non-zero (CTest failure) if any failed.
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

static int g_fail = 0;

// FP policy the tests were compiled under. A constant rather than a bare #ifdef
// so the golden-regression code still COMPILES in a fast build (no bit rot) and
// only the run is skipped.
#ifdef MESHMS_FP_FAST
inline constexpr bool kFpFast = true;
#else
inline constexpr bool kFpFast = false;
#endif

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::fprintf(stderr, "FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);     \
      ++g_fail;                                                                \
    }                                                                          \
  } while (0)

#define CHECK_NEAR(a, b, tol)                                                  \
  do {                                                                         \
    double _va = (a), _vb = (b), _d = std::fabs(_va - _vb);                    \
    if (!(_d <= (tol))) {                                                      \
      std::fprintf(stderr, "FAIL %s:%d  |%.17g - %.17g| = %.3g > %.3g\n",      \
                   __FILE__, __LINE__, _va, _vb, _d, (double)(tol));           \
      ++g_fail;                                                                \
    }                                                                          \
  } while (0)

// Relative comparison, for quantities whose absolute scale varies by molecule.
// Guards the |b| < 1 case so the tolerance never collapses to zero.
#define CHECK_REL(a, b, rel)                                                   \
  do {                                                                         \
    double _va = (a), _vb = (b);                                               \
    double _lim = (rel) * std::fmax(std::fabs(_vb), 1.0);                      \
    double _d = std::fabs(_va - _vb);                                          \
    if (!(_d <= _lim)) {                                                       \
      std::fprintf(stderr, "FAIL %s:%d  |%.17g - %.17g| = %.3g > %.3g (rel)\n", \
                   __FILE__, __LINE__, _va, _vb, _d, _lim);                    \
      ++g_fail;                                                                \
    }                                                                          \
  } while (0)

// Golden bit-regression only holds under MESHMS_FP=strict. In a fast build the
// integer topology itself may differ, so there is nothing to loosen -- the test
// exits with CTest's SKIP_RETURN_CODE and shows up as "Skipped", never silently
// as a pass. End-to-end coverage of a fast build lives in test_fp_gate.
#define MESHMS_SKIP_IF_FAST()                                                  \
  do {                                                                         \
    if (kFpFast) {                                                             \
      std::printf("SKIP: golden bit-regression requires MESHMS_FP=strict\n");  \
      return 77;                                                               \
    }                                                                          \
  } while (0)

// Consumer-side xyzr loader: skip blank/'#' rows, take the first 4 columns.
// Mirrors read_xyzr without using the internal API, so a facade test feeds the
// same atoms a consumer would pass to build_surface_from_array.
inline std::vector<std::array<double, 4>> load_xyzr_array(const std::string& path) {
  std::ifstream in(path);
  std::vector<std::array<double, 4>> out;
  std::string line;
  while (std::getline(in, line)) {
    std::size_t i = 0;
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    if (i == line.size() || line[i] == '#') continue;
    std::istringstream ls(line);
    std::vector<double> c;
    double v;
    while (ls >> v) c.push_back(v);
    if (c.size() >= 4) out.push_back({c[0], c[1], c[2], c[3]});
  }
  return out;
}

#define TEST_MAIN_RETURN()                                                     \
  do {                                                                         \
    if (g_fail) {                                                              \
      std::fprintf(stderr, "%d check(s) FAILED\n", g_fail);                    \
      return 1;                                                                \
    }                                                                          \
    std::printf("OK\n");                                                       \
    return 0;                                                                  \
  } while (0)
