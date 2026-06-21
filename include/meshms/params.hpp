#pragma once
// Run parameters --- faithful port of the params module (Para dataclass).
namespace meshms {

struct Para {
  double radius_probe = 1.4;  // solvent probe radius Rp
  double mesh_size = 0.5;     // target triangle edge length d (= arg_meshing(2))
};

}  // namespace meshms
