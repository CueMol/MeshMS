// Library version string (declared in the public facade meshms/meshms.hpp). The
// shared header-only contracts (csr.hpp/mesh.hpp/vec3.hpp) are already compiled
// by other translation units, so this TU no longer force-includes them.

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
