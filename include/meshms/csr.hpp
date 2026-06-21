#pragma once
// CSR (compressed-sparse-row) neighbour lists --- the OOM-wall fix:
// per-atom ragged storage instead of a global-kmax dense rectangle
// M_int[(M+1, kmax+1)].
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace meshms {

// Atom i's (1-based) intersecting SAS-balls are val[off[i] .. off[i+1]),
// ascending-unique. off has size M+2 (valid indices 0..M+1); val has size
// sum_i Row[i]. The dummy atom 0 has an empty range.
struct Neighbors {
  int M{0};
  std::vector<int32_t> off;  // size M+2
  std::vector<int32_t> val;  // size sum_i Row[i]

  std::span<const int32_t> of(int i) const {
    return std::span<const int32_t>(val.data() + off[i],
                                    static_cast<std::size_t>(off[i + 1] - off[i]));
  }
  int count(int i) const { return off[i + 1] - off[i]; }
};

}  // namespace meshms
