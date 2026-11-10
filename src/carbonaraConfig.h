#pragma once

#include <string>

// Set by each executable before constructing a randomMol instance.
extern std::string CARBONARA_ROOT;

// Returns <CARBONARA_ROOT>/probabilityInterpolation/<filename>.
std::string interpolationPath(const std::string& filename);
