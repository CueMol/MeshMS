#pragma once
// SAS arrangement: triple intersection points (I) and free circles --- faithful
// port of meshms/sas.py data_I_Cir().
//
// Storage follows CPP_PORT_DESIGN.md §2 (the OOM-wall fix): the maxs-family
// (I/Iijk/direction/high_I/hightvalue) are APPEND vectors (push_back as the
// triple-point counter s increments; maxs is only a reserve() hint). Ii is a
// per-atom growable list (CSR by In, §2.3 --- no kmax*(kmax-1) rectangle).
// I_circle is jagged by (atom, neighbour-row), §2.4 (no
// (M+1,kmax+1,nptmax) rectangle). The first real triple point is s=1; index 0
// is the reserved dummy row, preserving the 1-based-with-dummy-row-0 contract.
#include <array>
#include <cstdint>
#include <utility>
#include <vector>

#include "meshms/csr.hpp"
#include "meshms/geom.hpp"
#include "meshms/vec3.hpp"

namespace meshms {

// Triple intersection points and per-atom membership (data_I_Cir DataI half).
struct DataI {
  std::vector<Vec3> I;                           // [0] dummy, 1..nI real
  int nI{0};
  std::vector<std::array<int32_t, 3>> Iijk;      // [s] = (i,j,k); [0] dummy
  std::vector<std::array<int32_t, 3>> direction; // [s] = (dij,djk,dki); [0] dummy
  std::vector<int32_t> high_I;                   // [s]; [0] dummy
  std::vector<double> hightvalue;                // [s]; [0] dummy
  std::vector<std::vector<int32_t>> Ii;          // Ii[a] = point ids on atom a (In[a]==Ii[a].size())
  // I_circle[i] sized Row[i]+1 ([0] dummy); I_circle[i][row] = point ids on the
  // (i, row)-th intersection circle (I_circle_num == .size()).
  std::vector<std::vector<std::vector<int32_t>>> I_circle;
};

// Free (uncovered) intersection circles (data_I_Cir DataCir half).
struct DataCir {
  std::vector<std::array<double, 10>> circle;    // [k] = {_, ci,cj, Ax,Ay,Az, nx,ny,nz, r}; [0] dummy
  int ncircle{0};
  std::vector<std::vector<int32_t>> circleindex; // circleindex[a] = circle ids (ncircleindex == .size())
};

// Build I (triple points) and the free-circle list from the SAS-ball geometry
// and its intersection CSR (nb = interstructure(geom, Rp)).
std::pair<DataI, DataCir> data_I_Cir(const Geom& geom, const Neighbors& nb, double Rp);

}  // namespace meshms
