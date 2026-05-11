#include "MTRSimPlugin.hpp"

#include "MTRSim/MTRSim_filter_registration.hpp"

using namespace nx::core;

namespace
{
// Plugin Uuid
constexpr AbstractPlugin::IdType k_ID = *Uuid::FromString("f6bacee6-310a-4853-80f2-8092f4333560");
} // namespace

MTRSimPlugin::MTRSimPlugin()
: AbstractPlugin(k_ID, "MTRSim", "Plugin to hold highly experimental filters", "BlueQuartz Software, LLC")
{
  std::vector<::FilterCreationFunc> filterFuncs = ::GetPluginFilterList();
  for(const auto& filterFunc : filterFuncs)
  {
    addFilter(filterFunc);
  }
}

MTRSimPlugin::~MTRSimPlugin() = default;

AbstractPlugin::SIMPLMapType MTRSimPlugin::getSimplToSimplnxMap() const
{
  return {};
}

SIMPLNX_DEF_PLUGIN(MTRSimPlugin)
