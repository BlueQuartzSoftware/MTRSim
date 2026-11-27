#include "ReadMTRSimODF.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/DataGroup.hpp"

#include "LibMTRSim/ODFFileIO.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <exception>

using namespace nx::core;

// -----------------------------------------------------------------------------
ReadMTRSimODF::ReadMTRSimODF(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, ReadMTRSimODFInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(mesgHandler)
{
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
Result<> ReadMTRSimODF::operator()()
{
  m_MessageHandler.sendInfoMessage(fmt::format("Reading MTRSim ODF file: {}", m_InputValues->input_file.string()));

  std::vector<mtrsim::ODFFileComponent> components;
  try
  {
    components = mtrsim::ReadOdfComponents(m_InputValues->input_file, m_InputValues->hdf5_path_prefix);
  } catch(const std::exception& e)
  {
    return MakeErrorResult(-12010, fmt::format("ODF file read failed: {}", e.what()));
  }

  // The preflight/execute contract guarantees the DataStructure has one
  // Float64 array per component (created via CreateArrayAction in preflight).
  // Per-array size consistency is still verified below against store.getSize().
  const DataPath cellAttrMatPath = m_InputValues->output_image_geometry_path.createChildPath(m_InputValues->cell_attr_mat_name);

  for(size_t c = 0; c < components.size(); ++c)
  {
    if(m_ShouldCancel)
    {
      return {};
    }

    const DataPath componentPath = cellAttrMatPath.createChildPath(fmt::format("component_{}", c));
    auto& componentArray = m_DataStructure.getDataRefAs<Float64ArrayType>(componentPath);
    auto& store = componentArray.getDataStoreRef();

    const auto& srcValues = components[c].values;
    if(store.getSize() != srcValues.size())
    {
      return MakeErrorResult(-12012, fmt::format("Component {} size mismatch: array has {} tuples "
                                                 "but file provided {} values",
                                                 c, store.getSize(), srcValues.size()));
    }

    m_MessageHandler.sendInfoMessage(fmt::format("Copying component_{} ({} values)", c, srcValues.size()));

    std::copy(srcValues.begin(), srcValues.end(), store.begin());
  }

  return {};
}
