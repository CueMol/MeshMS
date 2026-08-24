// Library version string (declared in the public facade meshms/capi.hpp). This
// TU also force-compiles the header-only shared contracts so a syntax error in
// them surfaces at build time even before a consumer module exists.
#include "meshms/csr.hpp"
#include "meshms/mesh.hpp"
#include "meshms/vec3.hpp"

namespace meshms {
const char* version() { return "0.1.0"; }

const char* build_info() {
  return "0.1.0"
#ifdef MESHMS_FP_MODE
         " fp=" MESHMS_FP_MODE
#else
         " fp=strict"
#endif
#ifdef MESHMS_WITH_TBB
         " tbb=on"
#else
         " tbb=off"
#endif
      ;
}

int fp_mode() {
#ifdef MESHMS_FP_FAST
  return 1;
#else
  return 0;
#endif
}
}  // namespace meshms
