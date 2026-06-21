#pragma once
// ASCII PLY writer/reader for the triangulated SES mesh --- faithful port of
// meshms/ply.py. Vertices are TRUE coordinates (no positivity shift); face
// indices are 0-based; normals are an optional per-vertex property.
#include <string>
#include <vector>

#include "meshms/mesh.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

// Write an ASCII PLY. Coordinates are printed with %.6f. When vertex_normals is
// true, nx/ny/nz properties are emitted from `normals` (must have V.size() rows).
void write_ply(const std::string& path, const std::vector<Vec3>& V,
               const std::vector<Tri>& F, const std::vector<Vec3>* normals = nullptr,
               bool vertex_normals = false);

struct PlyMesh {
  std::vector<Vec3> V;
  std::vector<Tri> F;
};

// Parse the ASCII PLY that write_ply emits: 'element vertex N', optional
// nx/ny/nz props, 'element face M', 'property list uchar int vertex_indices',
// 'end_header', then N vertex lines (x y z [nx ny nz]) and M face lines
// ('3 i j k'). Faces are 0-based.
PlyMesh read_ply(const std::string& path);

}  // namespace meshms
