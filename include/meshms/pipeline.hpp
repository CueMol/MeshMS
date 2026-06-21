#pragma once
// End-to-end SES construction driver --- faithful port of the pipeline module.
//
// Wires the already-ported geometry/meshing stages together in the MATLAB order:
//   interstructure -> data_I_Cir -> data_Seg_Pat -> data_ext
//                  -> data_SESsphpat_convex (convex) + SESconcavepat (concave)
//                  + data_SEStorpat (toroidal)
//                  -> accumulate into ONE MeshState -> orient_faces -> (V, F, N).
//
// The three meshers accumulate into the SAME MeshState in the order
// convex -> concave -> toroidal, so the vertex/face indices line up with the
// Python pipeline exactly. orient_faces is faithful to pipeline.orient_faces.
//
// fuse=true runs the ID-boundary-fusion path (fuse_by_id, warn=False): after
// orient_faces the boundary vertices sharing a tag are welded by union-find, and
// the per-face normals N are dropped (Python `return V, F, None`).
//
// OUT OF SCOPE (stubbed, see the pipeline module): the build_surface "auto" NaN-retry
// jitter fallback (numpy-RNG-dependent, never triggered by the golden molecules;
// build_surface uses faithful coords == jitter=None semantics, noted as a TODO).
#include <cstdint>
#include <string>
#include <vector>

#include "meshms/exterior.hpp"
#include "meshms/geom.hpp"
#include "meshms/intersection.hpp"
#include "meshms/mesh.hpp"
#include "meshms/params.hpp"
#include "meshms/sas.hpp"
#include "meshms/sas_patches.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

// Make every triangle wind consistently with its outward normal: for each face,
// if dot(cross(b-a, c-a), N[f]) < 0 swap the last two vertex indices (in place).
// Faithful to pipeline.orient_faces (wn = cross(b-a,c-a); flip = (wn*N).sum<0).
void orient_faces(const std::vector<Vec3>& V, std::vector<Tri>& F,
                  const std::vector<Vec3>& N);

// Full SES mesh: vertices V, 0-based faces F (after orient_faces), per-face
// outward normals N, per-vertex outward normals NV (aligned with V), and
// per-vertex owning atom atom_id (1-based, 0 = unknown; aligned with V).
struct Surface {
  std::vector<Vec3> V;
  std::vector<Tri> F;
  std::vector<Vec3> N;             // per-face normals (empty after fuse)
  std::vector<Vec3> NV;            // per-vertex normals (always populated)
  std::vector<int32_t> atom_id;    // per-vertex owning atom (1-based, 0=unknown)
  // Per-face SES component type, aligned with F (MSMS face-type codes):
  //   3 = convex (contact), 2 = concave (spherical reentrant),
  //   1 = toroidal (toric reentrant). Stays aligned through orient_faces and the
  //   fuse degenerate-drop. Always populated by build_mesh.
  std::vector<uint8_t> ftype;
};

// Density-independent "RS components": everything in the SAS arrangement that
// depends ONLY on the geometry and the probe radius Rp, NOT on the mesh density
// (mesh_size). Compute once via compute_rs(), then call build_mesh() for any
// number of densities without recomputing the (expensive) arrangement. `geom`
// and `Rp` are stored so build_mesh() needs only the cache.
struct RSComponents {
  Geom geom;
  double Rp = 1.4;
  Neighbors inter;
  DataI di;
  DataCir dc;
  DataSeg ds;
  DataLoop dl;
  DataPat dp;
  Ext ext;
};

// Compute the density-independent RS components for `geom` at probe radius
// `radius_probe`. This is the expensive SAS-arrangement half of run().
RSComponents compute_rs(const Geom& geom, double radius_probe);

// Mesh the SES from precomputed RS components at the given `mesh_size` density.
// fuse=true welds tagged boundary vertices (fuse_by_id) and drops per-face N.
// Output (V/F/N) is bit-identical to run(rs.geom, {rs.Rp, mesh_size}, fuse).
Surface build_mesh(const RSComponents& rs, double mesh_size, bool fuse = false);

// Build the full SES mesh for `geom`. fuse=true welds tagged boundary vertices
// (fuse_by_id) and drops the per-face normals N. Equivalent to
// build_mesh(compute_rs(geom, para.radius_probe), para.mesh_size, fuse).
Surface run(const Geom& geom, const Para& para, bool fuse = false);

// Read an xyzr file and return the SES mesh. Uses faithful coordinates
// (jitter=None semantics); the auto NaN-retry jitter fallback is out of scope.
Surface build_surface(const std::string& xyzr_path, double radius_probe = 1.4,
                      double mesh_size = 0.5, bool fuse = false);

}  // namespace meshms
