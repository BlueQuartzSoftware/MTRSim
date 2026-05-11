#pragma once

#include "MTRSim/MTRSim_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"

#include <cstdint>
#include <string>

namespace nx::core
{

struct MTRSIM_EXPORT ComputeODFInputValues
{
  bool applySmoothing;
  double binSizeDeg;
  bool useMask;
  DataPath eulerAnglesPath;
  DataPath phasesPath;
  DataPath crystalStructuresPath;
  DataPath maskPath;
  DataPath outputImageGeometry;
  std::string cellAttrMatName;
  std::string componentName;
  int32_t nphi1;
  int32_t nPHI;
  int32_t nphi2;
  // 0 = Count-Density (raw normalized histogram, sums to 1 over Bunge bins;
  //                    matches MATLAB calc_ODF.m).
  // 1 = MUD          (per-row Jacobian sin(Phi) correction applied so values
  //                    are density on SO(3) divided by uniform density;
  //                    matches MTEX's plot(odf) MUD convention).
  uint64_t outputUnits;
};

/**
 * @class ComputeODF
 * @brief This algorithm implements support code for the ComputeODFFilter
 */

class MTRSIM_EXPORT ComputeODF
{
public:
  ComputeODF(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, ComputeODFInputValues* inputValues);
  ~ComputeODF() noexcept;

  ComputeODF(const ComputeODF&) = delete;
  ComputeODF(ComputeODF&&) noexcept = delete;
  ComputeODF& operator=(const ComputeODF&) = delete;
  ComputeODF& operator=(ComputeODF&&) noexcept = delete;

  Result<> operator()();

private:
  DataStructure& m_DataStructure;
  const ComputeODFInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;
};

} // namespace nx::core
