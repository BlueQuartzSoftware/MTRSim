#pragma once

#include "MTRSim/MTRSim_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"

#include <filesystem>
#include <vector>

namespace nx::core {

struct MTRSIM_EXPORT WriteMTRSimODFInputValues {
  std::filesystem::path outputFile;
  std::string hdf5PathPrefix;
  DataPath inputImageGeometry;
  std::vector<DataPath> odfComponents;
};

/**
 * @class WriteMTRSimODF
 * @brief Algorithm that writes the selected Float64 cell-data arrays on an
 * ODF ImageGeom out to a MATLAB-format MTRSim ODF HDF5 file via
 * mtrsim::writeODFFile. Reverses the (phi2=X, PHI=Y, phi1=Z) ImageGeom
 * convention back to the on-disk phi1-first axis order.
 */

class MTRSIM_EXPORT WriteMTRSimODF {
public:
  WriteMTRSimODF(DataStructure &dataStructure,
                 const IFilter::MessageHandler &mesgHandler,
                 const std::atomic_bool &shouldCancel,
                 WriteMTRSimODFInputValues *inputValues);
  ~WriteMTRSimODF() noexcept;

  WriteMTRSimODF(const WriteMTRSimODF &) = delete;
  WriteMTRSimODF(WriteMTRSimODF &&) noexcept = delete;
  WriteMTRSimODF &operator=(const WriteMTRSimODF &) = delete;
  WriteMTRSimODF &operator=(WriteMTRSimODF &&) noexcept = delete;

  Result<> operator()();

private:
  DataStructure &m_DataStructure;
  const WriteMTRSimODFInputValues *m_InputValues = nullptr;
  const std::atomic_bool &m_ShouldCancel;
  const IFilter::MessageHandler &m_MessageHandler;
};

} // namespace nx::core
