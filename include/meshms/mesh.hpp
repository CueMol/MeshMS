#pragma once
// Contiguous in-memory triangle mesh (CPP_PORT_DESIGN.md §5): the OpenGL-ready
// vertex/index buffers that MeshState accumulates into. Faces are 0-based, like
// the Python MeshState.to_arrays() output.
#include <array>
#include <cstdint>
#include <vector>

#include "meshms/vec3.hpp"

namespace meshms {

using Tri = std::array<int32_t, 3>;  // 0-based vertex indices

struct Mesh {
  std::vector<Vec3> V;
  std::vector<Tri> F;        // 0-based triangle indices
  std::vector<Vec3> N;       // per-face outward normals (or per-vertex)
  std::vector<int32_t> tag;  // per-vertex fusion tag id (-1 = interior/none)
};

}  // namespace meshms
