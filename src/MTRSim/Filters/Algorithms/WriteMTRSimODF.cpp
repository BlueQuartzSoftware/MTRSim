#include "WriteMTRSimODF.hpp"

#include "simplnx/DataStructure/DataArray.hpp"
#include "simplnx/DataStructure/Geometry/ImageGeom.hpp"

#include "LibMTRSim/ODFFileIO.hpp"

#include <fmt/format.h>

#include <array>
#include <cstdint>
#include <exception>
#include <vector>

using namespace nx::core;

// -----------------------------------------------------------------------------
WriteMTRSimODF::WriteMTRSimODF(DataStructure &dataStructure,
                               const IFilter::MessageHandler &mesgHandler,
                               const std::atomic_bool &shouldCancel,
                               WriteMTRSimODFInputValues *inputValues)
    : m_DataStructure(dataStructure), m_InputValues(inputValues),
      m_ShouldCancel(shouldCancel), m_MessageHandler(mesgHandler) {}

// -----------------------------------------------------------------------------
WriteMTRSimODF::~WriteMTRSimODF() noexcept = default;

// -----------------------------------------------------------------------------
Result<> WriteMTRSimODF::operator()() {
  m_MessageHandler(IFilter::Message::Type::Info,
                   fmt::format("Writing MTRSim ODF file: {}",
                               m_InputValues->outputFile.string()));

  const auto &geom = m_DataStructure.getDataRefAs<ImageGeom>(
      m_InputValues->inputImageGeometry);

  // Reverse the (phi2=X, PHI=Y, phi1=Z) ImageGeom convention back to phi1-first
  // on-disk order.
  const std::array<int64_t, 3> dimsPhi1PHIPhi2 = {
      static_cast<int64_t>(geom.getNumZCells()),
      static_cast<int64_t>(geom.getNumYCells()),
      static_cast<int64_t>(geom.getNumXCells())};

  const auto spacingXYZ = geom.getSpacing();
  const std::array<double, 3> spacingDegPhi1PHIPhi2 = {
      static_cast<double>(spacingXYZ[2]), static_cast<double>(spacingXYZ[1]),
      static_cast<double>(spacingXYZ[0])};

  // Per-component size validation against the geometry happens in preflight;
  // writeODFFile re-validates values.size() == dims[0]*dims[1]*dims[2] as a
  // final defence in the library layer.
  std::vector<mtrsim::ODFFileComponent> components;
  components.reserve(m_InputValues->odfComponents.size());

  for (std::size_t c = 0; c < m_InputValues->odfComponents.size(); ++c) {
    if (m_ShouldCancel) {
      return {};
    }

    const auto &componentArray = m_DataStructure.getDataRefAs<Float64Array>(
        m_InputValues->odfComponents[c]);
    const auto &store = componentArray.getDataStoreRef();

    m_MessageHandler(IFilter::Message::Type::Info,
                     fmt::format("Packing component_{} ({} values) from '{}'",
                                 c, store.getSize(),
                                 m_InputValues->odfComponents[c].toString()));

    mtrsim::ODFFileComponent comp;
    comp.values.assign(store.begin(), store.end());
    components.push_back(std::move(comp));
  }

  try {
    mtrsim::writeODFFile(m_InputValues->outputFile, dimsPhi1PHIPhi2,
                         spacingDegPhi1PHIPhi2, components,
                         m_InputValues->hdf5PathPrefix);
  } catch (const std::exception &e) {
    return MakeErrorResult(-12110,
                           fmt::format("ODF write failed: {}", e.what()));
  }

  return {};
}
