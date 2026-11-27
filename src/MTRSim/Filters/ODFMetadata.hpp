#pragma once

#include "simplnx/Common/StringLiteral.hpp"
#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"

#include <fmt/format.h>

#include <string>
#include <string_view>

/**
 * @namespace nx::core::mtrsim::odf_metadata
 * @brief Defines persistent units metadata for ODF components.
 */
namespace nx::core::mtrsim::odf_metadata
{
/**
 * @brief Names the metadata group directly under the ODF geometry.
 */
inline constexpr StringLiteralType k_GroupName = "ODF Metadata";
/**
 * @brief Identifies a component's units array.
 */
inline constexpr StringLiteralType k_UnitsSuffix = " Units";

/**
 * @brief Locates the one-tuple units StringArray for an ODF component.
 * @param geometryPath ODF geometry that owns the metadata group.
 * @param componentName Name of the ODF cell array.
 * @return Path to the persistent units array.
 */
inline DataPath CreateUnitsPath(const DataPath& geometryPath, std::string_view componentName)
{
  return geometryPath.createChildPath(k_GroupName.str()).createChildPath(fmt::format("{}{}", componentName, k_UnitsSuffix.view()));
}

/**
 * @brief Converts the validated units choice to its persistent label.
 * @param units Zero selects Count-Density; one selects MUD.
 * @return Exact units label for serialization and presentation.
 */
inline std::string UnitsToString(ChoicesParameter::ValueType units)
{
  return units == 0 ? "Count-Density" : "MUD";
}
} // namespace nx::core::mtrsim::odf_metadata
