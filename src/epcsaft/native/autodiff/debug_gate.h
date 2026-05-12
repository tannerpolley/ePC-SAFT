#pragma once

#include <cstdlib>
#include <string>

namespace epcsaft::autodiff {

inline bool debug_env_enabled(const char* name) {
    const char* value = std::getenv(name);
    return value != nullptr && std::string(value) == "1";
}

inline bool unsupported_derivative_debug_enabled() {
    return debug_env_enabled("EPCSAFT_ALLOW_DERIVATIVE_BACKEND_DEBUG");
}

inline std::string unsupported_derivative_debug_only_message(const std::string& surface) {
    return surface
        + " unsupported_derivative jacobian_backend is debug-only; "
          "set EPCSAFT_ALLOW_DERIVATIVE_BACKEND_DEBUG=1 to use it for explicit diagnostics.";
}

}  // namespace epcsaft::autodiff



