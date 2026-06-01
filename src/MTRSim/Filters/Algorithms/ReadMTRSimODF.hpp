#pragma once

#include "MTRSim/MTRSim_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"

#include <filesystem>
#include <string>

namespace nx::core {

struct MTRSIM_EXPORT ReadMTRSimODFInputValues {
  std::filesystem::path inputFile;
  std::string hdf5PathPrefix;
  DataPath outputImageGeometryPath;
  std::string cellAttrMatName;
};

/**
 * @class ReadMTRSimODF
 * @brief Algorithm that copies ODF component values out of a MATLAB-format
 * MTRSim ODF HDF5 file into the pre-created Float64 DataArrays under the
 * output ImageGeom's cell AttributeMatrix.
 */

class MTRSIM_EXPORT ReadMTRSimODF {
public:
  ReadMTRSimODF(DataStructure &dataStructure,
                const IFilter::MessageHandler &mesgHandler,
                const std::atomic_bool &shouldCancel,
                ReadMTRSimODFInputValues *inputValues);
  ~ReadMTRSimODF() noexcept;

  ReadMTRSimODF(const ReadMTRSimODF &) = delete;
  ReadMTRSimODF(ReadMTRSimODF &&) noexcept = delete;
  ReadMTRSimODF &operator=(const ReadMTRSimODF &) = delete;
  ReadMTRSimODF &operator=(ReadMTRSimODF &&) noexcept = delete;

  Result<> operator()();

private:
  DataStructure &m_DataStructure;
  const ReadMTRSimODFInputValues *m_InputValues = nullptr;
  const std::atomic_bool &m_ShouldCancel;
  const IFilter::MessageHandler &m_MessageHandler;
};

} // namespace nx::core
