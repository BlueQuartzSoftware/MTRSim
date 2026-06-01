#include "MTRSim.hpp"

#include "simplnx/DataStructure/DataArray.hpp"

using namespace nx::core;

// -----------------------------------------------------------------------------
MTRSim::MTRSim(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, MTRSimInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(mesgHandler)
{
}

// -----------------------------------------------------------------------------
MTRSim::~MTRSim() noexcept = default;

// -----------------------------------------------------------------------------
Result<> MTRSim::operator()()
{
  // NOTE: The algorithm body is implemented in a later task. For now this is a
  // no-op stub so the filter scaffolding (preflight + parameter validation) can
  // be built and tested independently.
  return {};
}
