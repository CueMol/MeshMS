#pragma once
// Minimal header-only check harness for the C++ port tests (no GoogleTest dep).
// One test_*.cpp == one executable; CHECK macros count failures into g_fail and
// TEST_MAIN_RETURN() makes main() exit non-zero (CTest failure) if any failed.
#include <cmath>
#include <cstdio>

static int g_fail = 0;

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

#define TEST_MAIN_RETURN()                                                     \
  do {                                                                         \
    if (g_fail) {                                                              \
      std::fprintf(stderr, "%d check(s) FAILED\n", g_fail);                    \
      return 1;                                                                \
    }                                                                          \
    std::printf("OK\n");                                                       \
    return 0;                                                                  \
  } while (0)
