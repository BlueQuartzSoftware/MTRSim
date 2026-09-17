#include "WriteODFSections.hpp"

using namespace nx::core;

WriteODFSections::WriteODFSections(DataStructure& dataStructure, const IFilter::MessageHandler& messageHandler, const std::atomic_bool& shouldCancel, const WriteODFSectionsInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_MessageHandler(messageHandler)
, m_ShouldCancel(shouldCancel)
, m_InputValues(inputValues)
{
}

WriteODFSections::~WriteODFSections() noexcept = default;

Result<> WriteODFSections::operator()()
{
  return MakeErrorResult(-12320, "Write ODF Sections execution is not implemented.");
}
