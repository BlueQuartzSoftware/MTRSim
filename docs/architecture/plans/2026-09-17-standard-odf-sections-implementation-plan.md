# Standard ODF Sections Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reusable EbsdLib renderer for conventional constant-phi2 ODF montages and expose it through MTRSim's `WriteODFSectionsFilter` with PNG and RGB Image Geometry outputs.

**Architecture:** EbsdLib owns grid validation, Count-Density-to-MUD conversion, Laue-class plotting limits, exact section interpolation, colorization, and RGBA composition. MTRSim owns SIMPLNX parameters, Data Structure validation/adaptation, color-preset resolution, optional RGB geometry creation, and atomic PNG output.

**Tech Stack:** C++20, EbsdLib `EbsdDataArray`/`canvas_ity`, Catch2, SIMPLNX filter/action APIs, `CreateColorMapParameter`, `ColorTableUtilities`, `AtomicFile`, `PngWriter`, CMake/Ninja.

**Spec:** `docs/architecture/2026-09-17-odf-section-presentation-design.md`

## Global Constraints

- Repository roots:
  - EbsdLib: `/Users/mjackson/Workspace7/EbsdLib`
  - MTRSim: `/Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim`
  - DREAM3D-NX: `/Users/mjackson/Workspace7/DREAM3DNX`
- Implement EbsdLib first on its `feature/odf_generation` branch; MTRSim consumes that public API through the local-EbsdLib DREAM3D-NX preset.
- Preserve all unrelated untracked content in MTRSim, including `.superpowers/`, `docs/Images/mtrsim_generated_microstructure.png`, and `output/avatar_slides/`.
- Do not modify simplnx OrientationAnalysis for v1.
- Input grids must be full, zero-origin, uniformly spaced Bunge grids with `(X,Y,Z) = (phi2,PHI,phi1)`.
- Support scalar Float32 and Float64 ODF arrays; render only MUD, converting Count-Density with the approved Jacobian.
- Require Crystal Structures, Material Names, and a positive phase index.
- Default to six exact phi2 sections, three sections per row, 512-pixel section width, and the SIMPLNX `Rainbow` color preset.
- Support Automatic `[0, displayed maximum]` and Manual `[minimum, maximum]` MUD scaling.
- Require PNG output, RGB Image Geometry output, or both.
- Use error range `-7500` through `-7599` in EbsdLib result objects and `-12300` through `-12399` in the MTRSim filter/algorithm.
- Before editing C++ in each repository, invoke `bluequartz-skills:clang-format`, `bluequartz-skills:cpp-code-documentation`, and `bluequartz-skills:cross-platform-cpp`; before MTRSim filter work also invoke `bluequartz-skills:filter-implementation`, `bluequartz-skills:error-messages`, and `bluequartz-skills:type-dispatch`.
- Follow test-driven development: add a failing focused test, run it and observe the expected failure, implement the minimum behavior, then rerun it green.
- Format every touched C/C++ file with `/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format --style="file:<repo-root>/.clang-format" --fallback-style=none -i`, then verify with `--dry-run --Werror`.
- Run all CTest commands serially with `-j 1`.
- Use signed-off BlueQuartz commit messages and never stage unrelated files.

---

## File map

### EbsdLib repository

- Create `Source/EbsdLib/Utilities/ODFSectionUtilities.h`: public grid/unit/result types and scientific preparation functions.
- Create `Source/EbsdLib/Utilities/ODFSectionUtilities.cpp`: validation, Laue limits, MUD conversion, section-angle generation, and phi2 interpolation.
- Create `Source/EbsdLib/Utilities/ODFSectionChrome.h`: canvas drawing declarations for axes, labels, title, information block, and color bar.
- Create `Source/EbsdLib/Utilities/ODFSectionChrome.cpp`: deterministic canvas chrome implementation.
- Create `Source/EbsdLib/Utilities/ODFSectionCompositor.h`: public render configuration/result/layout API.
- Create `Source/EbsdLib/Utilities/ODFSectionCompositor.cpp`: colorization, layout, canvas composition, and RGBA extraction.
- Create `Source/Test/ODFSectionUtilitiesTest.cpp`: analytical scientific-transform tests and invalid-input coverage.
- Create `Source/Test/ODFSectionCompositorTest.cpp`: layout, scale, color, metadata, and selected-pixel rendering tests.
- Modify `Source/EbsdLib/Utilities/SourceList.cmake`: compile and install the six new public/source files.
- Modify `Source/Test/CMakeLists.txt`: register both new Catch2 test files.

### MTRSim repository

- Create `src/MTRSim/Filters/WriteODFSectionsFilter.hpp`: filter traits, parameter keys, and filter declaration.
- Create `src/MTRSim/Filters/WriteODFSectionsFilter.cpp`: parameters, preflight validation/actions/preview, and algorithm dispatch.
- Create `src/MTRSim/Filters/Algorithms/WriteODFSections.hpp`: immutable execution input values and algorithm declaration.
- Create `src/MTRSim/Filters/Algorithms/WriteODFSections.cpp`: type adaptation, EbsdLib invocation, RGB copy, and atomic PNG writing.
- Create `test/WriteODFSectionsTest.cpp`: synthetic preflight, execution, equivalence, file, and error tests.
- Create `docs/WriteODFSectionsFilter.md`: user-facing filter documentation.
- Create `docs/Images/odf_phi2_sections.png`: approved canonical output from the external validation dataset.
- Modify `MTRSimPlugin.cmake`: register the new filter and algorithm.
- Modify `test/CMakeLists.txt`: compile the new unit test.
- Modify `docs/ComputeODFFilter.md`: cross-link the presentation filter and explain phase-isolated ODF input.

### External validation workspace

- Create `/Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/compute_and_write_odf_sections.d3dpipeline`: Read ANG, construct indexed/phase mask, compute MUD, and write the montage.
- Create `/Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/validation_notes.md`: commands, counts, scale, section angles, MTEX comparison, and artifact inventory.
- Generate all pipelines, logs, DREAM3D data, PNGs, and comparison scripts for the 324 MB dataset only in that validation directory.

---

### Task 1: Implement EbsdLib ODF scientific preparation

**Files:**
- Create: `/Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionUtilities.h`
- Create: `/Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionUtilities.cpp`
- Create: `/Users/mjackson/Workspace7/EbsdLib/Source/Test/ODFSectionUtilitiesTest.cpp`
- Modify: `/Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/SourceList.cmake`
- Modify: `/Users/mjackson/Workspace7/EbsdLib/Source/Test/CMakeLists.txt`

**Interfaces:**
- Produces:
  - `enum class ODFValueUnits : uint8_t { MUD = 0, CountDensity = 1 };`
  - `struct ODFGridView`
  - `struct ODFEulerPlotLimits`
  - `struct PreparedODFSections`
  - `struct ODFSectionPreparationResult`
  - `ODFEulerPlotLimits GetODFEulerPlotLimits(uint32_t laueOpsIndex);`
  - `std::vector<double> GenerateODFSectionAngles(uint32_t laueOpsIndex, size_t sectionCount);`
  - `ODFSectionPreparationResult PrepareODFSections(const ODFGridView& grid, size_t sectionCount);`
- Consumed by: Task 2's compositor and Task 4's MTRSim adapter.

- [ ] **Step 1: Add the public declarations and failing analytical tests**

Define the public contract in `ODFSectionUtilities.h`:

```cpp
enum class ODFValueUnits : uint8_t
{
  MUD = 0,
  CountDensity = 1
};

struct EbsdLib_EXPORT ODFGridView
{
  DoubleArrayType* values = nullptr;
  std::array<size_t, 3> dimensions = {0, 0, 0}; // phi1, PHI, phi2
  std::array<double, 3> originDeg = {0.0, 0.0, 0.0};
  std::array<double, 3> spacingDeg = {0.0, 0.0, 0.0};
  ODFValueUnits units = ODFValueUnits::MUD;
  uint32_t laueOpsIndex = CrystalStructure::UnknownCrystalStructure;
};

struct EbsdLib_EXPORT ODFEulerPlotLimits
{
  double phi1MaxDeg = 0.0;
  double phiMaxDeg = 0.0;
  double phi2MaxDeg = 0.0;
};

struct EbsdLib_EXPORT PreparedODFSections
{
  std::vector<double> mudValues; // section, PHI, phi1; phi1 fastest
  std::vector<double> sectionAnglesDeg;
  size_t phi1Count = 0;
  size_t phiCount = 0;
  double maximumDisplayedMUD = 0.0;
  ODFEulerPlotLimits limits;
};

struct EbsdLib_EXPORT ODFSectionPreparationResult
{
  PreparedODFSections sections;
  int32_t errorCode = 0;
  std::string errorMessage;
  explicit operator bool() const noexcept { return errorCode >= 0; }
};
```

Add `ODFSectionUtilitiesTest.cpp` to the test source list and add both new utility files to `Utilities/SourceList.cmake`. Write these first tests:

```cpp
TEST_CASE("ebsdlib::ODFSectionUtilities::LaueLimitsAndAngles", "[EbsdLib][ODFSectionUtilities]")
{
  const auto limits = GetODFEulerPlotLimits(CrystalStructure::Hexagonal_High);
  REQUIRE(limits.phi1MaxDeg == Approx(360.0));
  REQUIRE(limits.phiMaxDeg == Approx(90.0));
  REQUIRE(limits.phi2MaxDeg == Approx(60.0));
  REQUIRE(GenerateODFSectionAngles(CrystalStructure::Hexagonal_High, 6) ==
          std::vector<double>{0.0, 10.0, 20.0, 30.0, 40.0, 50.0});
}

TEST_CASE("ebsdlib::ODFSectionUtilities::CountDensityToMUD", "[EbsdLib][ODFSectionUtilities]")
{
  auto values = DoubleArrayType::CreateArray(32, "CountDensity", true);
  values->initializeWithValue(1.0 / 32.0);
  ODFGridView grid{values.get(), {4, 2, 4}, {0.0, 0.0, 0.0}, {90.0, 90.0, 90.0}, ODFValueUnits::CountDensity, CrystalStructure::Triclinic};
  const auto result = PrepareODFSections(grid, 4);
  REQUIRE(result);
  const double step = constants::k_PiD / 2.0;
  const double expected = (1.0 / 32.0) * 8.0 * constants::k_PiD * constants::k_PiD /
                          (step * step * step * std::sin(constants::k_PiD / 4.0));
  REQUIRE(result.sections.mudValues.front() == Approx(expected));
}
```

- [ ] **Step 2: Run the focused EbsdLib test and confirm the red state**

Run from `/Users/mjackson/Workspace7/EbsdLib`:

```bash
cmake --preset ebsdlib-release
cmake --build --preset "ebsdlib-release [BUILD]" --target EbsdLibUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/ebsdlib-release \
  -R 'EbsdLib::ebsdlib::ODFSectionUtilities' -j 1 --output-on-failure
```

Expected: compilation or link failure because the three declared functions are not implemented.

- [ ] **Step 3: Implement Laue limits, section angles, grid validation, and MUD conversion**

Implement the approved Laue table exactly:

```cpp
constexpr std::array<ODFEulerPlotLimits, 11> k_Limits = {{
    {360.0, 90.0, 60.0},   // Hexagonal High
    {360.0, 90.0, 90.0},   // Cubic High
    {360.0, 180.0, 60.0},  // Hexagonal Low
    {360.0, 90.0, 180.0},  // Cubic Low
    {360.0, 180.0, 360.0}, // Triclinic
    {360.0, 90.0, 360.0},  // Monoclinic, EbsdLib unique-axis convention
    {360.0, 90.0, 180.0},  // Orthorhombic
    {360.0, 180.0, 90.0},  // Tetragonal Low
    {360.0, 90.0, 90.0},   // Tetragonal High
    {360.0, 180.0, 120.0}, // Trigonal Low
    {360.0, 90.0, 120.0}   // Trigonal High
}};
```

Use error codes and messages:

```cpp
constexpr int32_t k_NullValues = -7500;
constexpr int32_t k_InvalidLaueClass = -7501;
constexpr int32_t k_InvalidGrid = -7502;
constexpr int32_t k_ValueCountMismatch = -7503;
constexpr int32_t k_InvalidValue = -7504;
constexpr int32_t k_InvalidSectionCount = -7505;
```

Validate zero origin, equal positive spacing, full-cube dimensions, tuple count, supported Laue index, section count of at least 2, and every source value being finite and nonnegative. Convert Count-Density with:

```cpp
const double stepRad = grid.spacingDeg[0] * constants::k_PiOver180D;
const double phiCenterRad = (static_cast<double>(phiIndex) + 0.5) * stepRad;
const double mud = countDensity * 8.0 * constants::k_PiD * constants::k_PiD /
                   (stepRad * stepRad * stepRad * std::sin(phiCenterRad));
```

Add complete Doxygen comments for every public enum, struct, field, and function,
including array order, degrees/radians, ownership, units, and failure semantics.

- [ ] **Step 4: Add failing interpolation and invalid-input tests**

Use a four-plane MUD grid with phi2-center values `{1,2,3,4}` repeated over phi1 and PHI. Assert that the exact zero-degree section equals `2.5`, the periodic midpoint between the 315-degree and 45-degree planes. Add table-driven assertions for all 11 Laue limits and errors for null values, transposed dimensions, nonuniform spacing, negative values, NaN, an unsupported crystal code, and section count 1.

```cpp
REQUIRE(result.sections.sectionAnglesDeg == std::vector<double>{0.0, 90.0, 180.0, 270.0});
REQUIRE(result.sections.mudValues.front() == Approx(2.5));
REQUIRE_FALSE(PrepareODFSections(negativeGrid, 4));
REQUIRE(PrepareODFSections(negativeGrid, 4).errorCode == -7504);
```

- [ ] **Step 5: Run the focused test and confirm the interpolation assertions fail**

Run the focused CTest command from Step 2. Expected: limit/conversion tests pass and interpolation/invalid-input cases fail until the remaining logic is implemented.

- [ ] **Step 6: Implement periodic phi2 interpolation and cropped section buffers**

For each exact section angle, compute its fractional position relative to the first cell center, wrap both neighboring phi2 indices modulo `nphi2`, and write section-major output with phi1 fastest:

```cpp
const double coordinate = (angleDeg - firstCenterDeg) / spacingDeg;
const auto floorCoordinate = std::floor(coordinate);
const auto lower = WrapIndex(static_cast<int64_t>(floorCoordinate), nphi2);
const auto upper = (lower + 1) % nphi2;
const double fraction = coordinate - floorCoordinate;
const double value = (1.0 - fraction) * fullMud[FullIndex(phi1, phi, lower)] +
                     fraction * fullMud[FullIndex(phi1, phi, upper)];
prepared.mudValues[SectionIndex(section, phi, phi1)] = value;
```

Crop phi1 and PHI by selecting every cell center strictly below the approved maximum. Record the maximum over the displayed interpolated values.

- [ ] **Step 7: Format, verify, and run focused tests**

```bash
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/EbsdLib/.clang-format" --fallback-style=none -i \
  Source/EbsdLib/Utilities/ODFSectionUtilities.h \
  Source/EbsdLib/Utilities/ODFSectionUtilities.cpp \
  Source/Test/ODFSectionUtilitiesTest.cpp
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/EbsdLib/.clang-format" --fallback-style=none \
  --dry-run --Werror \
  Source/EbsdLib/Utilities/ODFSectionUtilities.h \
  Source/EbsdLib/Utilities/ODFSectionUtilities.cpp \
  Source/Test/ODFSectionUtilitiesTest.cpp
cmake --build --preset "ebsdlib-release [BUILD]" --target EbsdLibUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/ebsdlib-release \
  -R 'EbsdLib::ebsdlib::ODFSectionUtilities' -j 1 --output-on-failure
```

Expected: all ODFSectionUtilities tests pass.

- [ ] **Step 8: Commit the scientific preparation layer**

```bash
git add Source/EbsdLib/Utilities/ODFSectionUtilities.h \
  Source/EbsdLib/Utilities/ODFSectionUtilities.cpp \
  Source/EbsdLib/Utilities/SourceList.cmake \
  Source/Test/ODFSectionUtilitiesTest.cpp Source/Test/CMakeLists.txt
git commit -s -m "ENH: Add ODF section preparation utilities"
```

---

### Task 2: Implement the EbsdLib ODF compositor

**Files:**
- Create: `/Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionChrome.h`
- Create: `/Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionChrome.cpp`
- Create: `/Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionCompositor.h`
- Create: `/Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionCompositor.cpp`
- Create: `/Users/mjackson/Workspace7/EbsdLib/Source/Test/ODFSectionCompositorTest.cpp`
- Modify: `/Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/SourceList.cmake`
- Modify: `/Users/mjackson/Workspace7/EbsdLib/Source/Test/CMakeLists.txt`

**Interfaces:**
- Consumes: Task 1's `ODFGridView`, `PreparedODFSections`, and `PrepareODFSections()`.
- Produces:
  - `enum class ODFScaleMode : uint8_t { Automatic = 0, Manual = 1 };`
  - `struct ODFSectionConfiguration`
  - `struct ODFSectionLayoutMetrics`
  - `struct ODFSectionResult`
  - `ODFSectionLayoutMetrics ComputeODFSectionLayout(const ODFSectionConfiguration& config);`
  - `ODFSectionResult ODFSectionCompositor::generateCompositeImage(const ODFSectionConfiguration& config) const;`
- Consumed by: Tasks 3 and 4.

- [ ] **Step 1: Add compositor declarations and failing configuration/layout tests**

Define this public shape in `ODFSectionCompositor.h`:

```cpp
enum class ODFScaleMode : uint8_t
{
  Automatic = 0,
  Manual = 1
};

struct EbsdLib_EXPORT ODFSectionConfiguration
{
  ODFGridView grid;
  size_t sectionCount = 6;
  size_t sectionsPerRow = 3;
  int32_t sectionWidth = 512;
  ODFScaleMode scaleMode = ODFScaleMode::Automatic;
  double manualMinimumMUD = 0.0;
  double manualMaximumMUD = 1.0;
  std::vector<float> colorControlPoints;
  std::string title;
  std::string materialName;
  int32_t phaseNumber = 1;
};

struct EbsdLib_EXPORT ODFSectionLayoutMetrics
{
  int32_t panelWidth = 0;
  int32_t panelHeight = 0;
  int32_t columns = 0;
  int32_t rows = 0;
  float fontPtSize = 0.0f;
  float margin = 0.0f;
  float panelSlotWidth = 0.0f;
  float panelSlotHeight = 0.0f;
  float titleHeight = 0.0f;
  float legendWidth = 0.0f;
  int32_t pageWidth = 0;
  int32_t pageHeight = 0;
};

struct EbsdLib_EXPORT ODFSectionResult
{
  UInt8ArrayType::Pointer image;
  int32_t width = 0;
  int32_t height = 0;
  std::vector<double> sectionAnglesDeg;
  double appliedMinimumMUD = 0.0;
  double appliedMaximumMUD = 0.0;
  int32_t errorCode = 0;
  std::string errorMessage;
  explicit operator bool() const noexcept { return errorCode >= 0 && image != nullptr; }
};
```

Write tests asserting default values and, for Hexagonal High with width 512, a data-panel height of 128, three columns, and two rows. Assert that invalid section width, row count, manual range, and fewer than two color controls return `-7510` through `-7513`.

- [ ] **Step 2: Run focused compositor tests and confirm the red state**

```bash
cmake --build --preset "ebsdlib-release [BUILD]" --target EbsdLibUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/ebsdlib-release \
  -R 'EbsdLib::ebsdlib::ODFSectionCompositor' -j 1 --output-on-failure
```

Expected: compilation or link failure for the undeclared compositor implementation.

- [ ] **Step 3: Implement deterministic layout and color mapping**

Compute panel height from the Laue limit and arrange the requested grid:

```cpp
metrics.panelWidth = config.sectionWidth;
metrics.panelHeight = static_cast<int32_t>(std::lround(
    static_cast<double>(config.sectionWidth) * limits.phiMaxDeg / limits.phi1MaxDeg));
metrics.columns = static_cast<int32_t>(config.sectionsPerRow);
metrics.rows = static_cast<int32_t>((config.sectionCount + config.sectionsPerRow - 1) /
                                    config.sectionsPerRow);
```

Use these deterministic page metrics:

```cpp
metrics.fontPtSize = std::max(10.0f, static_cast<float>(config.sectionWidth) / 24.0f);
metrics.margin = std::max(8.0f, static_cast<float>(config.sectionWidth) / 32.0f);
metrics.panelSlotWidth = static_cast<float>(metrics.panelWidth) + 2.0f * metrics.margin;
metrics.panelSlotHeight = static_cast<float>(metrics.panelHeight) +
                          3.0f * metrics.fontPtSize + 3.0f * metrics.margin;
metrics.titleHeight = metrics.fontPtSize + 2.0f * metrics.margin;
metrics.legendWidth = std::max(static_cast<float>(config.sectionWidth) / 3.0f, 180.0f);
metrics.pageWidth = static_cast<int32_t>(std::ceil(
    metrics.columns * metrics.panelSlotWidth + metrics.legendWidth + 2.0f * metrics.margin));
metrics.pageHeight = static_cast<int32_t>(std::ceil(
    metrics.titleHeight + metrics.rows * metrics.panelSlotHeight + metrics.margin));
```

Use these fields directly in tests and chrome so both share one source of truth.

Normalize each MUD value against the applied range, clamp to `[0,1]`, locate the surrounding `[position,r,g,b]` controls, and linearly interpolate RGB. Automatic range is exactly `[0, prepared.maximumDisplayedMUD]`; reject an automatic maximum of zero. Manual range requires `0 <= minimum < maximum` and clips endpoint values.

- [ ] **Step 4: Add failing chrome and RGBA composition tests**

Use a 4x2x4 Triclinic grid, four sections, two per row, a 128-pixel panel width, manual scale `[0,4]`, and two controls `[0,0,0,1]` and `[1,1,0,0]`. Assert:

```cpp
ODFSectionResult result = ODFSectionCompositor{}.generateCompositeImage(config);
REQUIRE(result);
REQUIRE(result.image->getNumberOfComponents() == 4);
REQUIRE(result.image->getNumberOfTuples() == static_cast<size_t>(result.width * result.height));
REQUIRE(result.sectionAnglesDeg == std::vector<double>{0.0, 90.0, 180.0, 270.0});
REQUIRE(result.appliedMinimumMUD == Approx(0.0));
REQUIRE(result.appliedMaximumMUD == Approx(4.0));
```

Check selected pixels inside a low-value panel, a high-value panel, the white background, a black axis, and both color-bar endpoints. Permit a per-channel tolerance of 1 for `canvas_ity` floating-point differences; do not use a broad whole-image tolerance to hide layout regressions.

- [ ] **Step 5: Implement `ODFSectionChrome` and RGBA composition**

Use embedded Lato/Fira fonts and `canvas_ity` to draw:

- White page background.
- Figure title.
- Each RGB panel scaled to its data rectangle.
- `phi1` and `PHI` axes with degree ticks.
- Per-panel `phi2 = N degrees` title.
- One shared MUD color bar and numeric endpoints.
- Phase number, material name, Laue group, input/rendered units, bin spacing,
  and section count.

Allocate the final image as:

```cpp
auto image = UInt8ArrayType::CreateArray(
    static_cast<size_t>(layout.pageWidth) * layout.pageHeight,
    {4ULL}, "ODFSectionsComposite", true);
context.get_image_data(image->getPointer(0), layout.pageWidth, layout.pageHeight,
                       layout.pageWidth * 4, 0, 0);
```

Add complete Doxygen comments for the public configuration, result, layout,
free layout function, compositor class, and `generateCompositeImage()` method.

- [ ] **Step 6: Format, build, and run both ODF section suites**

```bash
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/EbsdLib/.clang-format" --fallback-style=none -i \
  Source/EbsdLib/Utilities/ODFSectionChrome.h \
  Source/EbsdLib/Utilities/ODFSectionChrome.cpp \
  Source/EbsdLib/Utilities/ODFSectionCompositor.h \
  Source/EbsdLib/Utilities/ODFSectionCompositor.cpp \
  Source/Test/ODFSectionCompositorTest.cpp
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/EbsdLib/.clang-format" --fallback-style=none \
  --dry-run --Werror \
  Source/EbsdLib/Utilities/ODFSectionChrome.h \
  Source/EbsdLib/Utilities/ODFSectionChrome.cpp \
  Source/EbsdLib/Utilities/ODFSectionCompositor.h \
  Source/EbsdLib/Utilities/ODFSectionCompositor.cpp \
  Source/Test/ODFSectionCompositorTest.cpp
cmake --build --preset "ebsdlib-release [BUILD]" --target EbsdLibUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/ebsdlib-release \
  -R 'EbsdLib::ebsdlib::ODFSection' -j 1 --output-on-failure
```

Expected: utility and compositor suites pass.

- [ ] **Step 7: Commit the compositor**

```bash
git add Source/EbsdLib/Utilities/ODFSectionChrome.h \
  Source/EbsdLib/Utilities/ODFSectionChrome.cpp \
  Source/EbsdLib/Utilities/ODFSectionCompositor.h \
  Source/EbsdLib/Utilities/ODFSectionCompositor.cpp \
  Source/EbsdLib/Utilities/SourceList.cmake \
  Source/Test/ODFSectionCompositorTest.cpp Source/Test/CMakeLists.txt
git commit -s -m "ENH: Add standard ODF section compositor"
```

---

### Task 3: Add the MTRSim filter surface and preflight

**Files:**
- Create: `src/MTRSim/Filters/WriteODFSectionsFilter.hpp`
- Create: `src/MTRSim/Filters/WriteODFSectionsFilter.cpp`
- Create: `src/MTRSim/Filters/Algorithms/WriteODFSections.hpp`
- Create: `src/MTRSim/Filters/Algorithms/WriteODFSections.cpp`
- Create: `test/WriteODFSectionsTest.cpp`
- Modify: `MTRSimPlugin.cmake`
- Modify: `test/CMakeLists.txt`

**Interfaces:**
- Consumes: Task 2's `ODFSectionConfiguration`, `ComputeODFSectionLayout()`, and `ODFSectionCompositor`.
- Produces:
  - Filter UUID `5ec65380-b1a8-4aff-8bcd-310b5067b236`.
  - `WriteODFSectionsFilter` parameter contract.
  - `WriteODFSectionsInputValues` consumed by Task 4.
  - Preflight-created `Output ODF Sections/Cell Data/ODF Sections RGB` UInt8 array with three components.

- [ ] **Step 1: Register skeleton files and write failing preflight tests**

Add `WriteODFSectionsFilter` and `WriteODFSections` to `FilterList` and `AlgorithmList`, and add `WriteODFSectionsTest.cpp` to the plugin unit-test list.

Define these keys in the filter header:

```cpp
static constexpr StringLiteral k_InputImageGeometryPath_Key = "input_image_geometry_path";
static constexpr StringLiteral k_ODFArrayPath_Key = "odf_array_path";
static constexpr StringLiteral k_InputUnits_Key = "input_units_index";
static constexpr StringLiteral k_CrystalStructuresPath_Key = "crystal_structures_path";
static constexpr StringLiteral k_MaterialNamesPath_Key = "material_names_path";
static constexpr StringLiteral k_PhaseIndex_Key = "phase_index";
static constexpr StringLiteral k_Title_Key = "title";
static constexpr StringLiteral k_NumSections_Key = "number_of_sections";
static constexpr StringLiteral k_SectionsPerRow_Key = "sections_per_row";
static constexpr StringLiteral k_SectionWidth_Key = "section_width";
static constexpr StringLiteral k_ColorTablePreset_Key = "color_table_preset";
static constexpr StringLiteral k_ScaleMode_Key = "scale_mode_index";
static constexpr StringLiteral k_ManualMinimum_Key = "manual_minimum_mud";
static constexpr StringLiteral k_ManualMaximum_Key = "manual_maximum_mud";
static constexpr StringLiteral k_WritePng_Key = "write_png";
static constexpr StringLiteral k_OutputPngPath_Key = "output_png_file";
static constexpr StringLiteral k_SaveAsImageGeometry_Key = "save_as_image_geometry";
static constexpr StringLiteral k_OutputImageGeometryPath_Key = "output_image_geometry_path";
```

Use:

```cpp
SIMPLNX_DEF_FILTER_TRAITS(nx::core, WriteODFSectionsFilter,
                          "5ec65380-b1a8-4aff-8bcd-310b5067b236");
```

Build test data with a 72x36x72, five-degree `ImageGeom`, a scalar Float64 ODF array, a two-tuple UInt32 Crystal Structures array whose phase-1 value is `Hexagonal_High`, and `StringArray::CreateWithValues(..., {"Invalid", "Titanium (Alpha)"})`.

Write a valid-preflight test asserting:

- No errors.
- Exact section-angle preflight value `0, 10, 20, 30, 40, 50 degrees`.
- Phase/material/Laue preview values.
- One output Image Geometry action and one UInt8 RGB array action.
- Composite dimensions equal `ComputeODFSectionLayout()`.

- [ ] **Step 2: Run the focused MTRSim test and confirm the red state**

From `/Users/mjackson/Workspace7/DREAM3DNX`:

```bash
cmake --preset NX-Com-Qt69-Vtk96-Rel-EbsdLib
env LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 \
  cmake --build --preset "NX-Com-Qt69-Vtk96-Rel-EbsdLib [BUILD]" \
  --target MTRSimUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R 'MTRSim::WriteODFSectionsFilter' -j 1 --output-on-failure --timeout 1500
```

Expected: compilation failure because the new filter surface is not implemented.

- [ ] **Step 3: Implement parameters and links**

Use these parameter types and defaults:

```cpp
GeometrySelectionParameter  // ODF Image Geometry
ArraySelectionParameter     // scalar float32/float64 ODF array
ChoicesParameter            // MUD, Count-Density; default MUD
ArraySelectionParameter     // UInt32 Crystal Structures
DataPathSelectionParameter  // StringArray Material Names
Int32Parameter              // Phase Index = 1
StringParameter             // Figure Title
Int32Parameter              // Number of Sections = 6
Int32Parameter              // Sections Per Row = 3
Int32Parameter              // Section Width = 512
CreateColorMapParameter     // default "Rainbow"
ChoicesParameter            // Automatic, Manual; default Automatic
Float64Parameter            // Manual Minimum MUD = 0.0
Float64Parameter            // Manual Maximum MUD = 1.0
BoolParameter               // Write PNG = true
FileSystemPathParameter     // Output PNG, PathType::OutputFile, extension .png
BoolParameter               // Save as Image Geometry = true
DataGroupCreationParameter  // Output Image Geometry = "ODF Sections"
```

Link manual limits to Manual mode, PNG path to Write PNG, and output geometry path to Save as Image Geometry.

- [ ] **Step 4: Add failing table-driven preflight error tests**

Cover codes in `-12300..-12319` for:

- Both outputs disabled.
- ODF array not scalar or not Float32/Float64.
- ODF array not under the selected geometry's cell Attribute Matrix.
- Nonzero origin.
- Nonuniform, zero, or negative spacing.
- Dimensions not matching the full Bunge cube.
- Ensemble tuple-count mismatch.
- Phase index 0 or out of range.
- Unknown crystal structure.
- Section count below 2.
- Sections per row outside `[1, sectionCount]`.
- Section width below 128.
- Invalid preset.
- Manual minimum negative, non-finite, or greater than/equal to maximum.
- Non-PNG output extension.

Assert each message includes the selected path or actual invalid value.

- [ ] **Step 5: Implement preflight validation, actions, and preview**

Reorder Image Geometry XYZ values to EbsdLib `(phi1,PHI,phi2)` values:

```cpp
config.grid.dimensions = {imageGeom.getNumZCells(), imageGeom.getNumYCells(),
                          imageGeom.getNumXCells()};
config.grid.originDeg = {origin[2], origin[1], origin[0]};
config.grid.spacingDeg = {spacing[2], spacing[1], spacing[0]};
```

Resolve the phase tuple, call `GetODFEulerPlotLimits()` and `GenerateODFSectionAngles()`, call `ComputeODFSectionLayout()`, and append `CreateImageGeometryAction` plus `CreateArrayAction` only when geometry output is enabled. Use cell Attribute Matrix name `Cell Data` and array name `ODF Sections RGB`.

Create Image Geometry dimensions as `{pageWidth, pageHeight, 1}` and the RGB
array tuple shape as `{1, pageHeight, pageWidth}` with component shape `{3}`.

- [ ] **Step 6: Format, build, and run preflight tests**

```bash
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/.clang-format" \
  --fallback-style=none -i \
  src/MTRSim/Filters/WriteODFSectionsFilter.hpp \
  src/MTRSim/Filters/WriteODFSectionsFilter.cpp \
  src/MTRSim/Filters/Algorithms/WriteODFSections.hpp \
  src/MTRSim/Filters/Algorithms/WriteODFSections.cpp \
  test/WriteODFSectionsTest.cpp
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/.clang-format" \
  --fallback-style=none --dry-run --Werror \
  src/MTRSim/Filters/WriteODFSectionsFilter.hpp \
  src/MTRSim/Filters/WriteODFSectionsFilter.cpp \
  src/MTRSim/Filters/Algorithms/WriteODFSections.hpp \
  src/MTRSim/Filters/Algorithms/WriteODFSections.cpp \
  test/WriteODFSectionsTest.cpp
env LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 \
  cmake --build --preset "NX-Com-Qt69-Vtk96-Rel-EbsdLib [BUILD]" \
  --target MTRSimUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R 'MTRSim::WriteODFSectionsFilter.*Preflight' -j 1 \
  --output-on-failure --timeout 1500
```

Expected: all preflight and parameter tests pass. Execution tests remain absent until Task 4.

- [ ] **Step 7: Commit the filter surface**

```bash
git add MTRSimPlugin.cmake test/CMakeLists.txt test/WriteODFSectionsTest.cpp \
  src/MTRSim/Filters/WriteODFSectionsFilter.hpp \
  src/MTRSim/Filters/WriteODFSectionsFilter.cpp \
  src/MTRSim/Filters/Algorithms/WriteODFSections.hpp \
  src/MTRSim/Filters/Algorithms/WriteODFSections.cpp
git commit -s -m "FILT: Add ODF section image filter preflight"
```

---

### Task 4: Implement MTRSim execution and outputs

**Files:**
- Modify: `src/MTRSim/Filters/WriteODFSectionsFilter.cpp`
- Modify: `src/MTRSim/Filters/Algorithms/WriteODFSections.hpp`
- Modify: `src/MTRSim/Filters/Algorithms/WriteODFSections.cpp`
- Modify: `test/WriteODFSectionsTest.cpp`

**Interfaces:**
- Consumes: Task 3's filter arguments and Task 2's EbsdLib compositor.
- Produces: complete `WriteODFSections::operator()()` behavior, optional RGB geometry, and optional atomic PNG output.

- [ ] **Step 1: Add failing MUD/Count-Density equivalence tests**

Create two Data Structures containing the same analytical distribution:

```cpp
const double phiCenter = (static_cast<double>(phiIndex) + 0.5) * stepRad;
countDensity[index] = mud[index] * stepRad * stepRad * stepRad *
                      std::sin(phiCenter) /
                      (8.0 * std::numbers::pi * std::numbers::pi);
```

Execute once with MUD and once with Count-Density. Save geometry only, then assert equal dimensions and byte-identical `ODF Sections RGB` arrays.

- [ ] **Step 2: Run the focused execution test and confirm the red state**

Run:

```bash
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R 'MTRSim::WriteODFSectionsFilter.*Units' -j 1 \
  --output-on-failure --timeout 1500
```

Expected: execution returns the algorithm's not-implemented error or produces an empty RGB array.

- [ ] **Step 3: Implement typed ODF adaptation and EbsdLib invocation**

Populate `WriteODFSectionsInputValues` with every resolved argument. Use runtime type dispatch to copy Float32 or Float64 values to an EbsdLib `DoubleArrayType`:

```cpp
struct CopyODFToDouble
{
  template <class T>
  void operator()(const IDataArray& input, ebsdlib::DoubleArrayType& output) const
  {
    const auto& array = dynamic_cast<const DataArray<T>&>(input);
    for(usize index = 0; index < array.getSize(); ++index)
    {
      output[index] = static_cast<double>(array[index]);
    }
  }
};

ExecuteDataFunction(CopyODFToDouble{}, inputArray.getDataType(), inputArray,
                    *ebsdValues);
```

Resolve the color preset with `ColorTableUtilities::ExtractControlPoints()`, build the EbsdLib configuration, call `generateCompositeImage()`, and convert any EbsdLib error to a SIMPLNX error containing the ODF path and phase index.

- [ ] **Step 4: Add failing output-mode and invalid-value tests**

Add test sections for:

- Geometry only: output RGB array exists and is nonzero.
- PNG only: no output geometry is created and PNG exists.
- Both: decoded PNG dimensions and RGB pixels match the geometry array.
- Float32 and Float64 inputs render identically.
- Negative, NaN, and infinity source values fail with the actual flat index and value.
- All-zero automatic input fails.
- All-zero manual input succeeds and uses the first color.
- Cancellation before composition produces no destination file.
- A writer failure leaves an existing destination byte-identical.

Use `SimplnxCore::ReadImageFilter` in the test to decode the PNG and compare its RGB array to the generated geometry; MTRSimUnitTest already links `SimplnxCore`.

- [ ] **Step 5: Implement RGB geometry copying and atomic PNG writing**

Copy only RGB from the RGBA result:

```cpp
for(size_t tuple = 0; tuple < pixelCount; ++tuple)
{
  outputRgb[tuple * 3] = rgba[tuple * 4];
  outputRgb[tuple * 3 + 1] = rgba[tuple * 4 + 1];
  outputRgb[tuple * 3 + 2] = rgba[tuple * 4 + 2];
}
```

Write atomically:

```cpp
auto atomicResult = AtomicFile::Create(outputPath);
if(atomicResult.invalid())
{
  return ConvertResult(std::move(atomicResult));
}
AtomicFile atomicFile = std::move(atomicResult.value());
auto pngResult = PngWriter::WriteColorImage(atomicFile.tempFilePath().string(),
                                            result.width, result.height, 4,
                                            result.image->getPointer(0));
if(pngResult.first < 0)
{
  return MakeErrorResult(-12330, fmt::format(
      "Could not write ODF section PNG '{}': {} (writer code {}).",
      outputPath.string(), pngResult.second, pngResult.first));
}
return atomicFile.commit();
```

Check cancellation before conversion, before composition, before RGB copy, and before file creation.

- [ ] **Step 6: Add and pass the end-to-end Compute ODF test**

Build a small deterministic HCP Euler dataset using the existing `ComputeODFTest.cpp` helper pattern. Execute `ComputeODFFilter` with five-degree MUD output, create matching phase metadata, then execute `WriteODFSectionsFilter`. Assert six section angles in preflight, a valid RGB geometry, a finite positive automatic maximum, and a readable PNG.

- [ ] **Step 7: Format, build, and run the complete filter suite**

Run the MTRSim formatting commands from Task 3, then:

```bash
env LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 \
  cmake --build --preset "NX-Com-Qt69-Vtk96-Rel-EbsdLib [BUILD]" \
  --target MTRSimUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R 'MTRSim::WriteODFSectionsFilter' -j 1 \
  --output-on-failure --timeout 1500
```

Expected: all WriteODFSectionsFilter tests pass.

- [ ] **Step 8: Commit filter execution**

```bash
git add src/MTRSim/Filters/WriteODFSectionsFilter.cpp \
  src/MTRSim/Filters/Algorithms/WriteODFSections.hpp \
  src/MTRSim/Filters/Algorithms/WriteODFSections.cpp \
  test/WriteODFSectionsTest.cpp
git commit -s -m "ENH: Render standard ODF section images"
```

---

### Task 5: Document and validate the complete workflow

**Files:**
- Create: `docs/WriteODFSectionsFilter.md`
- Create: `docs/Images/odf_phi2_sections.png`
- Modify: `docs/ComputeODFFilter.md`
- Create externally: `/Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/compute_and_write_odf_sections.d3dpipeline`
- Create externally: `/Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/validation_notes.md`

**Interfaces:**
- Consumes: complete EbsdLib and MTRSim implementation.
- Produces: reproducible real-data evidence, canonical documentation image, and user-facing workflow documentation.

- [ ] **Step 1: Write the new filter documentation before generating its image**

Invoke `bluequartz-skills:filter-documentation` before editing the filter page.

Document these sections in `docs/WriteODFSectionsFilter.md`:

- Description and standard phi2-section convention.
- Input ODF axis mapping and full-grid requirement.
- MUD and Count-Density conversion formula.
- Crystal Structures, Material Names, phase index, and phase-isolated ODF requirement.
- Laue plotting-limit table from the spec.
- Exact section interpolation and default angles.
- Color-table and automatic/manual scale behavior.
- PNG and RGB Image Geometry outputs.
- Error conditions and references to MTEX/Euler-section conventions.

Add a `See Also` link from `docs/ComputeODFFilter.md` that directs users to compute a phase-isolated MUD or Count-Density array and then run `WriteODFSectionsFilter`.

- [ ] **Step 2: Create and preflight the external validation pipeline**

Invoke `bluequartz-skills:write-d3dpipeline`. Create the validation directory and pipeline with this sequence:

```text
ReadAngDataFilter(1.ang)
  -> MultiThresholdObjectsFilter(Confidence Index > 0, Bool "Indexed Mask")
  -> ComputeODFFilter(5 degrees, smoothing on, MUD, Indexed Mask)
  -> WriteODFSectionsFilter(phase 1, 6 sections, 3 per row, Rainbow, automatic scale)
  -> WriteDREAM3DFilter(validation.dream3d)
```

Use absolute input/output paths under `/Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation`. Preflight with:

```bash
/Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib/Bin/nxrunner \
  -p /Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/compute_and_write_odf_sections.d3dpipeline \
  -l /Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/preflight.log
```

Expected: preflight succeeds and reports phase 1, Titanium (Alpha), Hexagonal High, angles 0/10/20/30/40/50 degrees, and a three-by-two layout.

- [ ] **Step 3: Execute the pipeline and verify quantitative metadata**

```bash
/Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib/Bin/nxrunner \
  -e /Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/compute_and_write_odf_sections.d3dpipeline \
  -l /Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/execute.log
```

Verify the log reports exactly 2,493,822 contributing indexed pixels. Read the PNG and DREAM3D output independently and verify:

- Six sections.
- Three columns by two rows.
- Shared automatic scale starts at 0 MUD.
- RGB geometry dimensions equal PNG dimensions.
- RGB bytes equal decoded PNG RGB bytes.
- All ODF and output image values are finite.

- [ ] **Step 4: Compare the montage with the retained MTEX reference**

Use the existing reference files in `/Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_comparison`. Record in `validation_notes.md`:

- MTEX version and kernel halfwidth.
- EbsdLib/MTRSim grid spacing and section angles.
- Matching phase population.
- Applied MUD maximum.
- Visual correspondence of every major peak/valley across the six panels.
- Expected differences caused by MTEX's continuous kernel versus the 27-bin MTRSim smoothing stencil.
- The known MATLAB R2025b post-output teardown fault, clearly separated from EbsdLib/MTRSim behavior.

- [ ] **Step 5: Promote the approved montage into filter documentation**

Copy only the approved generated montage to `docs/Images/odf_phi2_sections.png`, reference it from `docs/WriteODFSectionsFilter.md`, and include a caption stating phase, input dataset, bin size, units, section angles, and color-table preset.

- [ ] **Step 6: Verify documentation and focused tests**

```bash
git diff --check -- docs/WriteODFSectionsFilter.md docs/ComputeODFFilter.md
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R 'MTRSim::(ComputeODFFilter|WriteODFSectionsFilter)' -j 1 \
  --output-on-failure --timeout 1500
```

Expected: no whitespace errors and all ComputeODF/WriteODFSections tests pass.

- [ ] **Step 7: Commit documentation and the approved image**

```bash
git add docs/WriteODFSectionsFilter.md docs/ComputeODFFilter.md \
  docs/Images/odf_phi2_sections.png
git commit -s -m "DOC: Document standard ODF section output"
```

---

### Task 6: Run cross-repository completion verification

**Files:**
- Verify only; no source mutation is expected.

**Interfaces:**
- Consumes: all prior tasks.
- Produces: final build/test/validation evidence and dependency-release readiness notes.

- [ ] **Step 1: Verify both worktrees and commit boundaries**

```bash
git -C /Users/mjackson/Workspace7/EbsdLib status --short
git -C /Users/mjackson/Workspace7/EbsdLib log --oneline --decorate -5
git -C /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim status --short
git -C /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim log --oneline --decorate -7
```

Expected: EbsdLib contains only the two planned implementation commits; MTRSim contains the preflight, execution, and documentation commits plus the existing design commit. Existing unrelated untracked MTRSim paths remain uncommitted.

- [ ] **Step 2: Run the complete EbsdLib build and test suite**

```bash
cd /Users/mjackson/Workspace7/EbsdLib
cmake --preset ebsdlib-release
cmake --build --preset "ebsdlib-release [BUILD]" -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/ebsdlib-release \
  -j 1 --output-on-failure
```

Expected: configure/build succeed and the full EbsdLib suite, including both ODF section suites, has zero failures.

- [ ] **Step 3: Run the full local-EbsdLib DREAM3D-NX build**

```bash
cd /Users/mjackson/Workspace7/DREAM3DNX
cmake --preset NX-Com-Qt69-Vtk96-Rel-EbsdLib
env LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 \
  cmake --build --preset "NX-Com-Qt69-Vtk96-Rel-EbsdLib [BUILD]" -j
```

Expected: full build succeeds, including documentation generation, MTRSim, and `nxrunner`.

- [ ] **Step 4: Run all MTRSim tests serially**

```bash
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R '^MTRSim::' -j 1 --output-on-failure --timeout 1500
```

Expected: every registered MTRSim test passes.

- [ ] **Step 5: Re-run the shipped validation pipeline from scratch**

Run the Task 5 preflight and execute commands again. Confirm a zero-error exit from both, compare the regenerated PNG and RGB geometry, and update `validation_notes.md` only if the measured scale or dimensions changed for an explained implementation reason.

- [ ] **Step 6: Record the normal-preset dependency boundary**

The normal `NX-Com-Qt69-Vtk96-Rel` preset consumes the released vcpkg EbsdLib package. Do not claim it builds the new filter until an EbsdLib release containing `ODFSectionCompositor` is available and the vcpkg dependency is updated. Record the exact EbsdLib commit tested by the local-EbsdLib preset and the future released version required by the normal preset.

- [ ] **Step 7: Run final diff and format gates**

```bash
git -C /Users/mjackson/Workspace7/EbsdLib diff --check upstream/develop...HEAD
git -C /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim diff --check
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/EbsdLib/.clang-format" --fallback-style=none \
  --dry-run --Werror \
  /Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionUtilities.h \
  /Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionUtilities.cpp \
  /Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionChrome.h \
  /Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionChrome.cpp \
  /Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionCompositor.h \
  /Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionCompositor.cpp \
  /Users/mjackson/Workspace7/EbsdLib/Source/Test/ODFSectionUtilitiesTest.cpp \
  /Users/mjackson/Workspace7/EbsdLib/Source/Test/ODFSectionCompositorTest.cpp
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/.clang-format" \
  --fallback-style=none --dry-run --Werror \
  /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/src/MTRSim/Filters/WriteODFSectionsFilter.hpp \
  /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/src/MTRSim/Filters/WriteODFSectionsFilter.cpp \
  /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/src/MTRSim/Filters/Algorithms/WriteODFSections.hpp \
  /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/src/MTRSim/Filters/Algorithms/WriteODFSections.cpp \
  /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/test/WriteODFSectionsTest.cpp
```

Expected: no diff-check or formatting errors.
