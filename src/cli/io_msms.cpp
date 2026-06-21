#include "io_msms.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace meshms_cli {

void write_msms(const std::string& base, const meshms::MeshResult& m, int natom,
                double density, double probe) {
  const std::size_t nv = m.verts.size();
  const std::size_t nf = m.faces.size();
  const bool have_aid = m.atom_id.size() == nv;
  const bool have_nrm = m.vnormals.size() == nv;
  const bool have_ftype = m.face_type.size() == nf;

  // Per-vertex type = smallest incident face_type (0 if none / unavailable).
  std::vector<int> vtype(nv, 0);
  if (have_ftype) {
    for (std::size_t f = 0; f < nf; ++f) {
      const int t = static_cast<int>(m.face_type[f]);
      for (int k = 0; k < 3; ++k) {
        const std::size_t v = m.faces[f][static_cast<std::size_t>(k)];
        if (vtype[v] == 0 || t < vtype[v]) vtype[v] = t;
      }
    }
  }

  auto open_out = [](const std::string& path) -> std::ofstream {
    std::ofstream fh(path);
    if (!fh) {
      std::fprintf(stderr, "error: cannot open for writing: %s\n", path.c_str());
      std::exit(2);
    }
    return fh;
  };

  char hdr[128], buf[256];
  std::ofstream vf = open_out(base + ".vert");
  vf << "# MSMS solvent excluded surface vertices, libMeshMS " << meshms::build_info() << "\n"
     << "#vertex #sphere density probe_r\n";
  std::snprintf(hdr, sizeof(hdr), "%7zu %7d %8.3f %8.3f", nv, natom, density, probe);
  vf << hdr << "\n";
  for (std::size_t i = 0; i < nv; ++i) {
    const std::array<double, 3>& v = m.verts[i];
    const std::array<double, 3> n = have_nrm ? m.vnormals[i] : std::array<double, 3>{0, 0, 0};
    const int sphere = have_aid ? static_cast<int>(m.atom_id[i]) : 0;
    std::snprintf(buf, sizeof(buf), "%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d",
                  v[0], v[1], v[2], n[0], n[1], n[2], 0, sphere, vtype[i]);
    vf << buf << "\n";
  }

  std::ofstream ff = open_out(base + ".face");
  ff << "# MSMS solvent excluded surface faces, libMeshMS " << meshms::build_info() << "\n"
     << "#faces #sphere density probe_r\n";
  std::snprintf(hdr, sizeof(hdr), "%7zu %7d %8.3f %8.3f", nf, natom, density, probe);
  ff << hdr << "\n";
  for (std::size_t f = 0; f < nf; ++f) {
    const std::array<std::uint32_t, 3>& t = m.faces[f];
    const int ft = have_ftype ? static_cast<int>(m.face_type[f]) : 0;
    std::snprintf(buf, sizeof(buf), "%6d %6d %6d %2d %6d", static_cast<int>(t[0]) + 1,
                  static_cast<int>(t[1]) + 1, static_cast<int>(t[2]) + 1, ft, 0);
    ff << buf << "\n";
  }
}

}  // namespace meshms_cli
