#include "carbonaraConfig.h"

#include <stdexcept>
#include <cstdlib>

std::string CARBONARA_ROOT;

std::string interpolationPath(const std::string& filename) {
    // The predictor explicitly initialises CARBONARA_ROOT from argv[20].
    // The trial generator does not. For that executable use the environment
    // passed by CarbonaraDataTools, without modifying global state here.
    std::string root = CARBONARA_ROOT;
    if (root.empty()) {
        const char* env_root = std::getenv("CARBONARA_ROOT");
        if (!env_root || !*env_root) {
            env_root = std::getenv("CARBONARADIR");
        }
        if (env_root && *env_root) {
            root = env_root;
        }
    }
    if (root.empty()) {
        throw std::runtime_error(
            "Carbonara installation directory not configured before interpolation loading"
        );
    }

    const std::string separator = root.back() == '/' ? "" : "/";
    return root + separator + "probabilityInterpolation/" + filename;
}
