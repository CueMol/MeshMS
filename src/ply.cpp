// Faithful port of the ply module's ASCII PLY writer + matching reader.
#include "meshms/ply.hpp"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace meshms {

void write_ply(const std::string& path, const std::vector<Vec3>& V,
               const std::vector<Tri>& F, const std::vector<Vec3>* normals,
               bool vertex_normals) {
  const std::size_t n_verts = V.size();
  const std::size_t n_faces = F.size();

  const bool write_normals = vertex_normals;
  if (write_normals) {
    if (normals == nullptr)
      throw std::invalid_argument("vertex_normals=true requires normals");
    if (normals->size() != n_verts)
      throw std::invalid_argument("normals row count != number of vertices");
  }

  std::ofstream fh(path);
  if (!fh) throw std::runtime_error("cannot open PLY for writing: " + path);

  fh << "ply\n";
  fh << "format ascii 1.0\n";
  fh << "element vertex " << n_verts << "\n";
  fh << "property float x\n";
  fh << "property float y\n";
  fh << "property float z\n";
  if (write_normals) {
    fh << "property float nx\n";
    fh << "property float ny\n";
    fh << "property float nz\n";
  }
  fh << "element face " << n_faces << "\n";
  fh << "property list uchar int vertex_indices\n";
  fh << "end_header\n";

  char buf[256];
  if (write_normals) {
    for (std::size_t i = 0; i < n_verts; ++i) {
      const Vec3& v = V[i];
      const Vec3& nrm = (*normals)[i];
      std::snprintf(buf, sizeof(buf), "%.6f %.6f %.6f %.6f %.6f %.6f", v.x, v.y, v.z,
                    nrm.x, nrm.y, nrm.z);
      fh << buf << "\n";
    }
  } else {
    for (std::size_t i = 0; i < n_verts; ++i) {
      const Vec3& v = V[i];
      std::snprintf(buf, sizeof(buf), "%.6f %.6f %.6f", v.x, v.y, v.z);
      fh << buf << "\n";
    }
  }

  for (const Tri& f : F) {
    std::snprintf(buf, sizeof(buf), "3 %d %d %d", static_cast<int>(f[0]),
                  static_cast<int>(f[1]), static_cast<int>(f[2]));
    fh << buf << "\n";
  }
}

PlyMesh read_ply(const std::string& path) {
  std::ifstream fh(path);
  if (!fh) throw std::runtime_error("cannot open PLY for reading: " + path);

  std::string line;
  std::size_t n_verts = 0;
  std::size_t n_faces = 0;
  bool has_normals = false;
  // Count vertex properties so we know how many tokens precede normals.
  bool in_vertex_props = false;

  while (std::getline(fh, line)) {
    std::istringstream iss(line);
    std::string tok;
    iss >> tok;
    if (tok == "element") {
      std::string kind;
      std::size_t count = 0;
      iss >> kind >> count;
      if (kind == "vertex") {
        n_verts = count;
        in_vertex_props = true;
      } else if (kind == "face") {
        n_faces = count;
        in_vertex_props = false;
      }
    } else if (tok == "property" && in_vertex_props) {
      std::string ptype, pname;
      iss >> ptype >> pname;
      if (pname == "nx" || pname == "ny" || pname == "nz") has_normals = true;
    } else if (tok == "end_header") {
      break;
    }
  }

  PlyMesh mesh;
  mesh.V.reserve(n_verts);
  mesh.F.reserve(n_faces);

  for (std::size_t i = 0; i < n_verts; ++i) {
    if (!std::getline(fh, line))
      throw std::runtime_error("PLY truncated reading vertices: " + path);
    std::istringstream iss(line);
    Vec3 v{};
    iss >> v.x >> v.y >> v.z;
    if (has_normals) {
      double nx, ny, nz;
      iss >> nx >> ny >> nz;
    }
    mesh.V.push_back(v);
  }

  for (std::size_t i = 0; i < n_faces; ++i) {
    if (!std::getline(fh, line))
      throw std::runtime_error("PLY truncated reading faces: " + path);
    std::istringstream iss(line);
    int count = 0;
    iss >> count;
    int a = 0, b = 0, c = 0;
    iss >> a >> b >> c;
    mesh.F.push_back(Tri{a, b, c});
  }

  return mesh;
}

}  // namespace meshms
