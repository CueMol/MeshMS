#pragma once
// CLI ASCII-PLY writer.
//
// INTENTIONAL DUPLICATION: write_ply reproduces the library's ply.cpp format
// byte-for-byte but is kept CLI-local on purpose -- the library's writer trades
// in C++20 internal types (Vec3/Tri), which cannot cross the C++17-clean facade
// boundary. The format itself is documented once in docs/API.md.
#include <string>

#include "meshms/meshms.hpp"

namespace meshms_cli {

// Write the mesh as ASCII PLY: positions (%.6f), optional per-vertex normals,
// then faces ("3 i j k"). Byte-for-byte the same format the library uses.
void write_ply(const std::string& path, const meshms::MeshResult& m,
               bool with_normals);

}  // namespace meshms_cli
