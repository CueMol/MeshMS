#pragma once
// Standalone mesh validation oracle --- faithful port of the mesh-check oracle.
// Edge-incidence statistics, watertightness, degenerate/duplicate faces, total
// surface area and signed volume. Used to gate SES mesh quality (the §5
// equivalence gate). Indexing here is the 0-based PLY/mesh convention.
#include <array>
#include <cstdint>
#include <vector>

#include "meshms/mesh.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

// 0.5 * sum |(b-a) x (c-a)|  (mesh_area in the mesh-check oracle).
double mesh_area(const std::vector<Vec3>& V, const std::vector<Tri>& F);

// sum a . (b x c) / 6  (signed_volume by the divergence theorem).
double signed_volume(const std::vector<Vec3>& V, const std::vector<Tri>& F);

struct ManifoldReport {
  int n_vertices;
  int n_faces;
  int nonfinite_vertices;  // vertices with a NaN or Inf component
  int degenerate_faces;
  int duplicate_faces;
  int boundary_edges;
  int nonmanifold_edges;
  bool watertight;
  double area;
  double signed_volume;
};

// Undirected edge incidence via a map keyed by sorted (u,v):
//   boundary_edges     = #edges used exactly once,
//   nonmanifold_edges  = #edges used >= 3 times,
//   watertight         = boundary==0 && nonmanifold==0 && nF>0,
//   degenerate_faces   = repeated-index OR |cross| < 1e-12,
//   nonfinite_vertices = #vertices with a non-finite component (a relaxed-FP
//                        build has no golden to compare against, so this is the
//                        primary NaN tripwire; see tests/test_fp_gate.cpp),
//   duplicate_faces    = #sorted-triples seen more than once.
ManifoldReport manifold_report(const std::vector<Vec3>& V, const std::vector<Tri>& F);

// Facade-layout overload: identical math over the facade MeshResult arrays
// (x,y,z triples and 0-based uint32 index triples) without converting the whole
// mesh to Vec3/Tri first. Bit-identical results.
ManifoldReport manifold_report(const std::vector<std::array<double, 3>>& V,
                               const std::vector<std::array<std::uint32_t, 3>>& F);

}  // namespace meshms
