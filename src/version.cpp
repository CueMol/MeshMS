// Library version string (declared in the public facade meshms/capi.hpp). This
// TU also force-compiles the header-only shared contracts so a syntax error in
// them surfaces at build time even before a consumer module exists.
#include "meshms/csr.hpp"
#include "meshms/mesh.hpp"
#include "meshms/vec3.hpp"

namespace meshms {
const char* version() { return "0.1.0"; }
}  // namespace meshms
