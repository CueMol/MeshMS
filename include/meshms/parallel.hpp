#ifndef MESHMS_PARALLEL_HPP
#define MESHMS_PARALLEL_HPP
// Internal parallel-for over MeshMS's output-preserving per-entity stages. Each
// stage writes ONLY to its own per-entity slot and merges serially in ascending
// index order afterwards, so the iteration order never changes the bytes emitted
// --- only the wall-clock. The backend is selected at build time:
//   * MESHMS_WITH_TBB defined -> Intel oneTBB tbb::parallel_for
//   * otherwise               -> a plain serial loop (identical output, no dep)
// Because iteration order is irrelevant to the result, the stages' former
// dynamic scheduling needs no emulation; TBB's auto_partitioner over a
// blocked_range is an equivalent dynamic split. This header is INTERNAL: never
// included by the public facade (meshms.hpp), so no TBB type crosses the boundary.
#if defined(MESHMS_WITH_TBB)
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#endif

namespace meshms {

// Run body(i) for every i in [begin, end). Iterations are independent; the caller
// guarantees each one writes only its own slot. parallel_for is blocking, so the
// by-reference capture of body stays valid for the whole call.
template <typename Body>
inline void parallel_for(int begin, int end, Body&& body) {
#if defined(MESHMS_WITH_TBB)
  tbb::parallel_for(tbb::blocked_range<int>(begin, end),
                    [&body](const tbb::blocked_range<int>& r) {
                      for (int i = r.begin(); i != r.end(); ++i) body(i);
                    });
#else
  for (int i = begin; i < end; ++i) body(i);
#endif
}

}  // namespace meshms

#endif  // MESHMS_PARALLEL_HPP
