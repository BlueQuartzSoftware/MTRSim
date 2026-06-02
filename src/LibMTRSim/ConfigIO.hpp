#pragma once

#include "SimulationParams.hpp"
#include "libmtrsim_export.h"

#include <filesystem>

namespace mtrsim {

/**
 * @brief Parse an MTRSim config JSON into a SimulationParams.
 *
 * Recognized keys: xLen,yLen,zLen, dx,dy,dz, volumeFractions, thetaList,
 * nuggetVariance (parsed but unused by the simulation), odfInputPath, seed.
 * Any comments / unknown keys are ignored. Fields absent from the JSON keep
 * their SimulationParams defaults.
 *
 * @throws std::runtime_error if the file cannot be opened or the JSON is invalid.
 */
LIBMTRSIM_EXPORT SimulationParams parseConfigJson(const std::filesystem::path& path);

} // namespace mtrsim
