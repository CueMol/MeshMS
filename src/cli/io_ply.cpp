#include "io_ply.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

namespace meshms_cli {

void write_ply(const std::string& path, const meshms::MeshResult& m,
               bool with_normals) {
  std::ofstream fh(path);
  if (!fh) {
    std::fprintf(stderr, "error: cannot open PLY for writing: %s\n", path.c_str());
    std::exit(2);
  }
  const std::size_t nv = m.verts.size();
  const std::size_t nf = m.faces.size();
  fh << "ply\n" << "format ascii 1.0\n" << "element vertex " << nv << "\n"
     << "property float x\n" << "property float y\n" << "property float z\n";
  if (with_normals)
    fh << "property float nx\n" << "property float ny\n" << "property float nz\n";
  fh << "element face " << nf << "\n"
     << "property list uchar int vertex_indices\n" << "end_header\n";

  char buf[256];
  for (std::size_t i = 0; i < nv; ++i) {
    const std::array<double, 3>& v = m.verts[i];
    if (with_normals) {
      const std::array<double, 3>& n = m.vnormals[i];
      std::snprintf(buf, sizeof(buf), "%.6f %.6f %.6f %.6f %.6f %.6f", v[0], v[1],
                    v[2], n[0], n[1], n[2]);
    } else {
      std::snprintf(buf, sizeof(buf), "%.6f %.6f %.6f", v[0], v[1], v[2]);
    }
    fh << buf << "\n";
  }
  for (const std::array<std::uint32_t, 3>& f : m.faces) {
    std::snprintf(buf, sizeof(buf), "3 %d %d %d", static_cast<int>(f[0]),
                  static_cast<int>(f[1]), static_cast<int>(f[2]));
    fh << buf << "\n";
  }
}

}  // namespace meshms_cli
