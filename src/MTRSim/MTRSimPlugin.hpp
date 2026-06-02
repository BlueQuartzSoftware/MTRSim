#pragma once

#include "MTRSim/MTRSim_export.hpp"

#include "simplnx/Plugin/AbstractPlugin.hpp"

class MTRSIM_EXPORT MTRSimPlugin : public nx::core::AbstractPlugin
{
public:
  MTRSimPlugin();
  ~MTRSimPlugin() override;

  MTRSimPlugin(const MTRSimPlugin&) = delete;
  MTRSimPlugin(MTRSimPlugin&&) = delete;

  MTRSimPlugin& operator=(const MTRSimPlugin&) = delete;
  MTRSimPlugin& operator=(MTRSimPlugin&&) = delete;

  /**
   * @brief Returns a map of UUIDs as strings, where SIMPL UUIDs are keys to
   * their simplnx counterpart
   * @return std::map<nx::core::Uuid, nx::core::Uuid>
   */
  SIMPLMapType getSimplToSimplnxMap() const override;
};
