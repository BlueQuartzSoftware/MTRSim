#include "ODFFileIO.hpp"

#include <catch2/catch.hpp>

#include <cstdint>
#include <filesystem>
#include <numeric>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace fs = std::filesystem;
using namespace mtrsim;

namespace
{
const fs::path k_DataDir{MTRSIM_TEST_DATA_DIR};
const fs::path k_ExemplarFile = k_DataDir / "simulation_ODF.h5";
} // namespace

TEST_CASE("ODFFileIO::readODFMetadata against exemplar", "[odffileio]")
{
  REQUIRE(fs::exists(k_ExemplarFile));

  const ODFFileMetadata md = readODFMetadata(k_ExemplarFile);

  CHECK(md.numComponents == 3);
  CHECK(md.dimsPhi1PHIPhi2[0] == 72);
  CHECK(md.dimsPhi1PHIPhi2[1] == 36);
  CHECK(md.dimsPhi1PHIPhi2[2] == 72);
  CHECK(md.spacingDegPhi1PHIPhi2[0] == Approx(5.0).margin(1e-9));
  CHECK(md.spacingDegPhi1PHIPhi2[1] == Approx(5.0).margin(1e-9));
  CHECK(md.spacingDegPhi1PHIPhi2[2] == Approx(5.0).margin(1e-9));
}

TEST_CASE("ODFFileIO::readODFMetadata throws on missing file", "[odffileio]")
{
  const fs::path bogus = k_DataDir / "this_file_does_not_exist_XYZ.h5";
  REQUIRE_FALSE(fs::exists(bogus));

  CHECK_THROWS_AS(readODFMetadata(bogus), std::runtime_error);
}

TEST_CASE("ODFFileIO::readODFComponents against exemplar", "[odffileio]")
{
  REQUIRE(fs::exists(k_ExemplarFile));

  const std::vector<ODFFileComponent> comps = readODFComponents(k_ExemplarFile);

  REQUIRE(comps.size() == 3);

  const std::size_t expectedSize = static_cast<std::size_t>(72 * 36 * 72); // 186624
  for(std::size_t ci = 0; ci < comps.size(); ++ci)
  {
    CHECK(comps[ci].values.size() == expectedSize);

    bool allNonNegative = true;
    for(double v : comps[ci].values)
    {
      if(v < 0.0)
      {
        allNonNegative = false;
        break;
      }
    }
    CHECK(allNonNegative);

    // The exemplar ODFs are produced by MATLAB and stored un-normalised; each
    // component sums to approximately 1.0 within ~1% (observed: 1.0009,
    // 0.9999, 0.9973). Downstream code (see src/app/main.cpp) explicitly
    // re-normalises after reading.
    const double sum = std::accumulate(comps[ci].values.begin(), comps[ci].values.end(), 0.0);
    CHECK(sum == Approx(1.0).margin(1e-2));
  }
}

TEST_CASE("ODFFileIO round-trip read/write/read", "[odffileio]")
{
  REQUIRE(fs::exists(k_ExemplarFile));

  const ODFFileMetadata srcMd = readODFMetadata(k_ExemplarFile);
  const std::vector<ODFFileComponent> srcComps = readODFComponents(k_ExemplarFile);

  const fs::path outPath = fs::temp_directory_path() / "mtrsim_roundtrip.h5";

  std::error_code ec;
  fs::remove(outPath, ec); // ignore error if file doesn't exist

  writeODFFile(outPath, srcMd.dimsPhi1PHIPhi2, srcMd.spacingDegPhi1PHIPhi2, srcComps);

  REQUIRE(fs::exists(outPath));

  const ODFFileMetadata rtMd = readODFMetadata(outPath);
  const std::vector<ODFFileComponent> rtComps = readODFComponents(outPath);

  CHECK(rtMd.numComponents == srcMd.numComponents);
  CHECK(rtMd.dimsPhi1PHIPhi2 == srcMd.dimsPhi1PHIPhi2);
  CHECK(rtMd.spacingDegPhi1PHIPhi2[0] == Approx(srcMd.spacingDegPhi1PHIPhi2[0]).margin(1e-9));
  CHECK(rtMd.spacingDegPhi1PHIPhi2[1] == Approx(srcMd.spacingDegPhi1PHIPhi2[1]).margin(1e-9));
  CHECK(rtMd.spacingDegPhi1PHIPhi2[2] == Approx(srcMd.spacingDegPhi1PHIPhi2[2]).margin(1e-9));

  REQUIRE(rtComps.size() == srcComps.size());
  for(std::size_t ci = 0; ci < srcComps.size(); ++ci)
  {
    CHECK(rtComps[ci].values == srcComps[ci].values);
  }

  // Cleanup
  fs::remove(outPath, ec);
}
