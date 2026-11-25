# Combined ODF Computation and Presentation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor MTRSim to one public `ComputeODFFilter` that computes a selected-phase quantitative ODF, persists its units, and optionally produces the standard annotated section montage.

**Architecture:** EbsdLib remains the reusable scientific renderer and receives targeted chrome improvements. MTRSim persists units as DataGroup/StringArray objects, extends `ComputeODFFilter` with phase selection and presentation parameters, invokes an internal presentation helper after numeric computation, and removes the unreleased standalone `WriteODFSectionsFilter`.

**Tech Stack:** C++20, EbsdLib `EbsdDataArray`/`canvas_ity`, Catch2, SIMPLNX filter/action APIs, `CreateColorMapParameter`, `CreateDataGroupAction`, `CreateStringArrayAction`, `AtomicFile`, `PngWriter`, CMake/Ninja.

**Spec:** `docs/architecture/2026-09-17-odf-section-presentation-design.md`

## Global Constraints

- Repository roots:
  - EbsdLib: `/Users/mjackson/Workspace7/EbsdLib`
  - MTRSim: `/Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim`
  - DREAM3D-NX: `/Users/mjackson/Workspace7/DREAM3DNX`
- Starting heads:
  - EbsdLib `feature/odf_generation`: `c67e78c986acbc9addde7a94c2447861f18bdaa9`
  - MTRSim `topic/create_mtr_sim_filter`: `75be1a43497cf9383417463690f720982a87a2d3`
- Preserve all unrelated untracked MTRSim content: `.superpowers/`, `docs/Images/mtrsim_generated_microstructure.png`, and `output/avatar_slides/`.
- Do not modify simplnx OrientationAnalysis.
- `ComputeODFFilter` is the only public ODF computation/presentation filter after this plan.
- Compute exactly one selected phase per execution; the optional mask further restricts that phase.
- Preserve the full numeric ODF output in selected units and persist those units at `<ODF Geometry>/ODF Metadata/<Component Name> Units`.
- Accepted unit strings are exactly `MUD` and `Count-Density`; there is no manual units override.
- Presentation is optional and always rendered in MUD; Count-Density presentation is labeled `MUD (from Count-Density)`.
- Default presentation: six sections, three per row, width 512, `Rainbow Blended White`, automatic MUD scaling.
- Use `φ₂`, `φ₁`, `Φ`, and `°`; axis ticks are 10 degrees when pixel spacing is at least 24 pixels, otherwise 20 degrees, without unit suffixes on numeric labels.
- Omit phase-dependent preflight updated values when ensemble values are unavailable; never emit `Deferred until execution`.
- Keep quotient-based full-grid tolerance `1.0e-6` and phi2-fastest full-source indexing.
- Before C++ edits invoke the applicable BlueQuartz clang-format, C++ documentation, cross-platform C++, filter-implementation, error-message, type-dispatch, and progress/cancellation skills; use superpowers:test-driven-development and systematic-debugging.
- Format every touched C/C++ file with `/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format --style="file:<repo-root>/.clang-format" --fallback-style=none -i`, then verify with `--dry-run --Werror`.
- Run all CTest commands serially with `-j 1`.
- Use signed-off BlueQuartz commit messages; do not push.

---

## File map

### EbsdLib

- Modify `Source/EbsdLib/Utilities/ODFSectionChrome.h`: deterministic title/color-bar/tick helper declarations.
- Modify `Source/EbsdLib/Utilities/ODFSectionChrome.cpp`: Greek typography, degree symbol, adaptive ticks, and source-unit-aware MUD title.
- Modify `Source/Test/ODFSectionCompositorTest.cpp`: helper and rendered-chrome regressions.

### MTRSim

- Create `src/MTRSim/Filters/ODFMetadata.hpp`: persistent units naming, paths, and exact string conversion.
- Modify `src/MTRSim/Filters/ComputeODFFilter.hpp/.cpp`: phase/material parameters, units actions, presentation parameters/actions/previews, and combined execution.
- Modify `src/MTRSim/Filters/Algorithms/ComputeODF.hpp/.cpp`: selected-phase accumulation.
- Create `src/MTRSim/Filters/Algorithms/GenerateODFSections.hpp/.cpp`: internal presentation helper refactored from the standalone algorithm.
- Delete `src/MTRSim/Filters/WriteODFSectionsFilter.hpp/.cpp`.
- Delete `src/MTRSim/Filters/Algorithms/WriteODFSections.hpp/.cpp` after the helper refactor lands.
- Modify `MTRSimPlugin.cmake`: remove standalone filter/algorithm registration and register the internal helper source through the plugin's algorithm list.
- Modify `test/ComputeODFTest.cpp`: selected-phase, units metadata, and numeric regressions.
- Create `test/ComputeODFPresentationTest.cpp`: combined-filter presentation/preflight/output tests.
- Delete `test/WriteODFSectionsTest.cpp` after its useful cases are migrated.
- Modify `test/CMakeLists.txt`: replace standalone test registration.
- Modify `docs/ComputeODFFilter.md`: combined workflow documentation and final image.
- Delete `docs/WriteODFSectionsFilter.md`.
- Replace `docs/Images/odf_phi2_sections.png` with the final combined-filter output.
- Update the external pipeline/scripts under `/Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/`.

---

### Task 1: Refine EbsdLib ODF section chrome

**Files:**
- Modify: `/Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionChrome.h`
- Modify: `/Users/mjackson/Workspace7/EbsdLib/Source/EbsdLib/Utilities/ODFSectionChrome.cpp`
- Modify: `/Users/mjackson/Workspace7/EbsdLib/Source/Test/ODFSectionCompositorTest.cpp`

**Interfaces:**
- Produces:
  - `std::string FormatODFSectionTitle(double angleDeg);`
  - `double SelectODFAxisTickInterval(double maximumDeg, int32_t pixelLength);`
  - `std::vector<double> GenerateODFAxisTicks(double maximumDeg, int32_t pixelLength);`
  - `std::string GetODFColorBarTitle(ODFValueUnits sourceUnits);`
- Consumed by: `DrawODFSectionChrome()` and Task 3's Count-Density presentation tests.

- [ ] **Step 1: Write failing text/tick helper tests**

Add tests before implementations:

```cpp
REQUIRE(FormatODFSectionTitle(10.0) == "φ₂ = 10°");
REQUIRE(SelectODFAxisTickInterval(360.0, 900) == Approx(10.0));
REQUIRE(SelectODFAxisTickInterval(360.0, 512) == Approx(20.0));
REQUIRE(GenerateODFAxisTicks(90.0, 128) == std::vector<double>{0.0, 20.0, 40.0, 60.0, 80.0, 90.0});
REQUIRE(GetODFColorBarTitle(ODFValueUnits::MUD) == "MUD");
REQUIRE(GetODFColorBarTitle(ODFValueUnits::CountDensity) == "MUD (from Count-Density)");
```

Assert every generated numeric tick string contains neither `deg` nor `°`.

- [ ] **Step 2: Run the focused tests and capture RED**

```bash
cmake --build --preset "ebsdlib-release [BUILD]" --target EbsdLibUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/ebsdlib-release \
  -R 'EbsdLib::ebsdlib::ODFSectionCompositor' -j 1 --output-on-failure
```

Expected: compile/link failure because the four helper functions do not exist.

- [ ] **Step 3: Implement the chrome helpers**

Use this exact tick rule:

```cpp
double SelectODFAxisTickInterval(double maximumDeg, int32_t pixelLength)
{
  const double tenDegreePixelSpacing = static_cast<double>(pixelLength) /
                                       (maximumDeg / 10.0);
  return tenDegreePixelSpacing >= 24.0 ? 10.0 : 20.0;
}
```

`GenerateODFAxisTicks()` starts at zero, advances by the selected interval, and
always appends `maximumDeg` when the interval does not land on it exactly.

Use UTF-8 literals `φ₂`, `φ₁`, `Φ`, and `°`. Keep PHI increasing downward.
Draw short hash marks at every generated tick and numeric labels without units.
Render `φ₁` centered below and `Φ` vertically centered to the left.

- [ ] **Step 4: Update rendered-pixel regressions**

Retain existing panel/background/color-bar assertions. Add separate renders for
MUD and Count-Density sources and assert:

- Dimensions/layout remain identical.
- Panel data remains identical for mathematically equivalent values.
- Chrome differs in the expected title/info region.
- Selected non-background pixels exist where the Greek titles/axes and dense
  tick marks are expected.

- [ ] **Step 5: Format, test, and commit**

```bash
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/EbsdLib/.clang-format" --fallback-style=none -i \
  Source/EbsdLib/Utilities/ODFSectionChrome.h \
  Source/EbsdLib/Utilities/ODFSectionChrome.cpp \
  Source/Test/ODFSectionCompositorTest.cpp
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/EbsdLib/.clang-format" --fallback-style=none \
  --dry-run --Werror \
  Source/EbsdLib/Utilities/ODFSectionChrome.h \
  Source/EbsdLib/Utilities/ODFSectionChrome.cpp \
  Source/Test/ODFSectionCompositorTest.cpp
cmake --build --preset "ebsdlib-release [BUILD]" --target EbsdLibUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/ebsdlib-release \
  -R 'EbsdLib::ebsdlib::ODFSection' -j 1 --output-on-failure
git add Source/EbsdLib/Utilities/ODFSectionChrome.h \
  Source/EbsdLib/Utilities/ODFSectionChrome.cpp \
  Source/Test/ODFSectionCompositorTest.cpp
git commit -s -m "ENH: Refine ODF section annotations"
```

Expected: all focused ODF section tests pass and the EbsdLib worktree is clean.

---

### Task 2: Add selected-phase computation and persistent units metadata

**Files:**
- Create: `src/MTRSim/Filters/ODFMetadata.hpp`
- Modify: `src/MTRSim/Filters/ComputeODFFilter.hpp`
- Modify: `src/MTRSim/Filters/ComputeODFFilter.cpp`
- Modify: `src/MTRSim/Filters/Algorithms/ComputeODF.hpp`
- Modify: `src/MTRSim/Filters/Algorithms/ComputeODF.cpp`
- Modify: `test/ComputeODFTest.cpp`

**Interfaces:**
- Produces:
  - `ComputeODFFilter::k_MaterialNames_Key = "material_names_path"`
  - `ComputeODFFilter::k_PhaseIndex_Key = "phase_index"`
  - `mtrsim::odf_metadata::k_GroupName = "ODF Metadata"`
  - `CreateUnitsPath(const DataPath& geometryPath, std::string_view componentName)`
  - `UnitsToString(ChoicesParameter::ValueType units)` returning exactly `MUD` or `Count-Density`
  - Selected-phase field in `ComputeODFInputValues`.
- Consumed by: Task 3's combined presentation.

- [ ] **Step 1: Add failing selected-phase tests**

Extend the synthetic helper with a two-tuple Material Names array and set phase
labels containing phases 1 and 2. Add tests proving:

```cpp
args.insertOrAssign(ComputeODFFilter::k_PhaseIndex_Key, int32{1});
// Only phase-1 voxels contribute.

args.insertOrAssign(ComputeODFFilter::k_PhaseIndex_Key, int32{2});
// Only phase-2 voxels contribute and that phase's Laue operators are used.
```

Use a mixed HCP/Cubic fixture whose nonzero-bin count or exact sentinel bins
distinguish the selections. Add a mask case that removes one selected-phase
voxel while non-selected voxels remain irrelevant.

- [ ] **Step 2: Add failing units-metadata action tests**

For Create-New and Append modes, assert preflight proposes:

```text
ODF/ODF Metadata/Component 1 Units
```

as a one-tuple StringArray initialized to the selected exact string. Apply
actions in Preflight mode and assert the value is immediately readable. Cover
MUD and Count-Density and reject units-array collisions or a non-DataGroup
`ODF Metadata` object.

- [ ] **Step 3: Run RED tests**

```bash
env LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 \
  cmake --build --preset "NX-Com-Qt69-Vtk96-Rel-EbsdLib [BUILD]" \
  --target MTRSimUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R 'MTRSim::ComputeODFFilter' -j 1 --output-on-failure --timeout 1500
```

Expected: missing parameter/action behavior and mixed-phase accumulation failures.

- [ ] **Step 4: Implement metadata helpers and parameters**

In `ODFMetadata.hpp` define:

```cpp
namespace nx::core::mtrsim::odf_metadata
{
inline constexpr StringLiteral k_GroupName = "ODF Metadata";
inline constexpr StringLiteral k_UnitsSuffix = " Units";

inline DataPath CreateUnitsPath(const DataPath& geometryPath,
                                std::string_view componentName)
{
  return geometryPath.createChildPath(k_GroupName.str())
      .createChildPath(fmt::format("{}{}", componentName, k_UnitsSuffix.view()));
}

inline std::string UnitsToString(ChoicesParameter::ValueType units)
{
  return units == 0 ? "Count-Density" : "MUD";
}
}
```

Add required Material Names and Phase Index parameters. Increment the parameter
version. Validate equal ensemble tuple counts, positive/in-range phase index,
and the material object type whenever values are available.

Append `CreateDataGroupAction` only when the group does not exist, then append
`CreateStringArrayAction({1}, unitsPath, unitsString)`. Reject wrong-type group
objects and units-path collisions.

- [ ] **Step 5: Enforce selected-phase accumulation**

Carry `selectedPhase` into `ComputeODFInputValues` and `AccumulateWorker`. Replace
the former nonzero-phase condition with:

```cpp
const int32_t phase = m_Phases[i];
if(phase != m_SelectedPhase || (m_Mask != nullptr && !m_Mask->isTrue(i)))
{
  progressMessenger.sendProgressMessage(1);
  continue;
}
```

Validate the selected phase's crystal code at execution before launching the
parallel algorithm. Update progress/completion messages to name the phase.

- [ ] **Step 6: Verify DREAM3D serialization of units metadata**

Execute Compute ODF, write the Data Structure to a temporary DREAM3D file, read
it into a fresh Data Structure, and assert the exact metadata path/value
survives. Use existing `WriteDREAM3DFilter`/`ReadDREAM3DFilter` test patterns;
do not infer persistence from in-memory behavior alone.

- [ ] **Step 7: Format, test, and commit**

```bash
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/.clang-format" \
  --fallback-style=none -i \
  src/MTRSim/Filters/ODFMetadata.hpp \
  src/MTRSim/Filters/ComputeODFFilter.hpp \
  src/MTRSim/Filters/ComputeODFFilter.cpp \
  src/MTRSim/Filters/Algorithms/ComputeODF.hpp \
  src/MTRSim/Filters/Algorithms/ComputeODF.cpp \
  test/ComputeODFTest.cpp
/opt/local/LLVM-19.1.1-macOS-ARM64/bin/clang-format \
  --style="file:/Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/.clang-format" \
  --fallback-style=none --dry-run --Werror \
  src/MTRSim/Filters/ODFMetadata.hpp \
  src/MTRSim/Filters/ComputeODFFilter.hpp \
  src/MTRSim/Filters/ComputeODFFilter.cpp \
  src/MTRSim/Filters/Algorithms/ComputeODF.hpp \
  src/MTRSim/Filters/Algorithms/ComputeODF.cpp \
  test/ComputeODFTest.cpp
env LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 \
  cmake --build --preset "NX-Com-Qt69-Vtk96-Rel-EbsdLib [BUILD]" \
  --target MTRSimUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R 'MTRSim::ComputeODFFilter' -j 1 --output-on-failure --timeout 1500
git add src/MTRSim/Filters/ODFMetadata.hpp \
  src/MTRSim/Filters/ComputeODFFilter.hpp \
  src/MTRSim/Filters/ComputeODFFilter.cpp \
  src/MTRSim/Filters/Algorithms/ComputeODF.hpp \
  src/MTRSim/Filters/Algorithms/ComputeODF.cpp \
  test/ComputeODFTest.cpp
git commit -s -m "ENH: Add phase and units provenance to ODF output"
```

---

### Task 3: Fold presentation into Compute ODF and remove the standalone filter

**Files:**
- Create: `src/MTRSim/Filters/Algorithms/GenerateODFSections.hpp`
- Create: `src/MTRSim/Filters/Algorithms/GenerateODFSections.cpp`
- Modify: `src/MTRSim/Filters/ComputeODFFilter.hpp`
- Modify: `src/MTRSim/Filters/ComputeODFFilter.cpp`
- Delete: `src/MTRSim/Filters/WriteODFSectionsFilter.hpp`
- Delete: `src/MTRSim/Filters/WriteODFSectionsFilter.cpp`
- Delete after refactor: `src/MTRSim/Filters/Algorithms/WriteODFSections.hpp`
- Delete after refactor: `src/MTRSim/Filters/Algorithms/WriteODFSections.cpp`
- Modify: `MTRSimPlugin.cmake`
- Create: `test/ComputeODFPresentationTest.cpp`
- Delete: `test/WriteODFSectionsTest.cpp`
- Modify: `test/CMakeLists.txt`

**Interfaces:**
- Consumes: Task 1's final EbsdLib chrome and Task 2's selected phase/units metadata.
- Produces:
  - One public `ComputeODFFilter` with optional presentation parameters.
  - Internal `GenerateODFSectionsInputValues` and `GenerateODFSections` helper.
  - No `WriteODFSectionsFilter` registration, UUID, or public class.

- [ ] **Step 1: Migrate standalone execution tests before deleting code**

Create `ComputeODFPresentationTest.cpp` and port the useful standalone cases to
execute the public Compute filter from Euler/phase inputs:

- MUD and Count-Density numeric modes produce equivalent section/color-bar data.
- Count-Density chrome uses `MUD (from Count-Density)`.
- PNG-only, RGB-only, both, and numeric-only output modes.
- RGB geometry equals decoded PNG RGB channels.
- Float32/Float64 renderer adaptation remains covered through the internal
  helper where the public Compute output is Float64.
- Automatic/manual scale, invalid values, cancellation, and atomic destination
  preservation.
- Deferred metadata conservative sizing and execution shrink.
- Float32 fractional bin spacing.

First register the new test alongside the old one and demonstrate the new
combined-filter tests fail because presentation parameters/actions do not exist.

- [ ] **Step 2: Add presentation parameters to Compute ODF**

Move these unreleased keys/defaults to `ComputeODFFilter`:

```text
title                         = "ODF Sections"
number_of_sections            = 6
sections_per_row              = 3
section_width                 = 512
color_table_preset            = "Rainbow Blended White"
scale_mode_index              = Automatic
manual_minimum_mud            = 0.0
manual_maximum_mud            = 1.0
write_png                     = false
output_png_file
save_as_image_geometry        = true
output_presentation_geometry  = "ODF Sections"
```

Presentation is enabled when `write_png || save_as_image_geometry`. Both false
is valid numeric-only mode.

- [ ] **Step 3: Integrate presentation preflight**

After numeric ODF and units-metadata actions, validate presentation parameters.
When ensemble values exist, use the selected phase's exact Laue limits/angles.
When the Crystal Structures store is Empty/EmptyOutOfCore:

- Do not index Crystal Structures or Material Names.
- Use the Triclinic maximum layout for presentation actions.
- Omit Section Angles, Material, Laue Class, and Euler Plot Limits updated
  values entirely.
- Emit no string containing `Deferred until execution`.
- Retain known phase index, numeric/rendered units, full-grid dimensions,
  requested layout/scale, and conservative composite dimensions.

Add an assertion over every updated value:

```cpp
for(const auto& value : preflight.outputValues)
{
  REQUIRE(value.value.find("Deferred until execution") == std::string::npos);
}
```

- [ ] **Step 4: Refactor the internal helper**

Rename `WriteODFSectionsInputValues`/`WriteODFSections` to
`GenerateODFSectionsInputValues`/`GenerateODFSections`. Remove filter-argument
semantics. Pass the known units enum/string from Compute execution; do not read
or ask for a units override.

Keep reviewed behavior for:

- Typed input copy.
- Color-control normalization.
- EbsdLib invocation.
- Contextual errors.
- Conservative geometry resize.
- RGB copy.
- Atomic PNG write.
- Stage cancellation.

- [ ] **Step 5: Chain compute and presentation execution**

`ComputeODFFilter::executeImpl()` first runs `ComputeODF`. If valid and
presentation is enabled, construct `GenerateODFSectionsInputValues` using the
same selected phase, material, output units, numeric ODF path, and presentation
arguments, then invoke the helper.

Numeric ODF and units metadata remain valid if presentation is disabled. A
presentation failure returns an invalid filter result and does not commit a PNG.

- [ ] **Step 6: Remove the standalone public filter**

Remove filter registration, source/header files, standalone algorithm files,
and standalone test registration only after all combined tests are green. Search
the repository:

```bash
rg -n "WriteODFSections(Filter)?" src test MTRSimPlugin.cmake
```

Expected: no public implementation, registration, or test references remain.
Task 4 removes the standalone documentation and validation-pipeline references.

- [ ] **Step 7: Run focused tests, format, and commit**

Format every touched/created C++ file with the pinned command and verify
`--dry-run --Werror`. Then:

```bash
env LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 \
  cmake --build --preset "NX-Com-Qt69-Vtk96-Rel-EbsdLib [BUILD]" \
  --target MTRSimUnitTest -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R 'MTRSim::ComputeODFFilter' -j 1 --output-on-failure --timeout 1500
git add MTRSimPlugin.cmake test/CMakeLists.txt \
  src/MTRSim/Filters test/ComputeODFTest.cpp test/ComputeODFPresentationTest.cpp \
  test/WriteODFSectionsTest.cpp
git commit -s -m "REV/ENH: Combine ODF computation and presentation"
```

Use `git diff --cached --name-status` before committing to verify intended
creates/modifications/deletions and no protected paths.

---

### Task 4: Update documentation and regenerate real-data validation

**Files:**
- Modify: `docs/ComputeODFFilter.md`
- Delete: `docs/WriteODFSectionsFilter.md`
- Replace: `docs/Images/odf_phi2_sections.png`
- Modify externally: `/Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/compute_and_write_odf_sections.d3dpipeline`
- Modify externally: validation scripts, logs, results, and notes in the same directory.

**Interfaces:**
- Consumes: the completed combined filter.
- Produces: the final public documentation and corrected reference montage.

- [ ] **Step 1: Update the validation pipeline**

Invoke `bluequartz-skills:write-d3dpipeline`. The final sequence is:

```text
ReadAngDataFilter(1.ang)
  -> MultiThresholdObjectsFilter(Confidence Index > 0, Bool "Indexed Mask")
  -> ComputeODFFilter(
       phase 1,
       5 degrees,
       smoothing on,
       MUD,
       Indexed Mask,
       six sections,
       three per row,
       Rainbow Blended White,
       automatic scale,
       PNG + RGB geometry)
  -> WriteDREAM3DFilter(validation.dream3d)
```

There is no standalone presentation-filter step.

- [ ] **Step 2: Preflight and execute**

```bash
/Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib/Bin/nxrunner \
  -p /Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/compute_and_write_odf_sections.d3dpipeline \
  -l /Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/preflight.log
/Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib/Bin/nxrunner \
  -e /Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/compute_and_write_odf_sections.d3dpipeline \
  -l /Users/Shared/Data/MTR_Data/RR_MTR_Examples/1/odf_section_validation/execute.log
```

Assert preflight output contains no `Deferred until execution` string.

- [ ] **Step 3: Extend independent validation**

Verify:

- 2,493,822 contributing phase-1 indexed points and 1,506,178 excluded points.
- Numeric ODF path, finite values, full 72x36x72 grid, and five-degree spacing.
- Units metadata path exists, equals `MUD`, and survives the DREAM3D file.
- Six sections at 0/10/20/30/40/50 degrees.
- PNG/RGB dimensions and byte equality.
- `Rainbow Blended White` default.
- MUD maximum agrees with independent extraction within `1e-12`.
- MTEX panel correlations remain 0.9730–0.9783.
- The standalone filter UUID/name is absent from the saved pipeline.

Add a Count-Density validation variant or focused runner invocation that proves
the image reports `MUD (from Count-Density)` without manual unit input.

- [ ] **Step 4: Visually inspect the final montage**

Use image inspection and confirm:

- `φ₂ = N°` titles render with correct glyphs.
- `φ₁` and vertical `Φ` render correctly.
- Default-width axes use 20-degree ticks with no unit suffixes.
- Tick marks/labels do not overlap or clip.
- Rainbow Blended White matches the filter default.
- MUD/phase/material information remains readable.

- [ ] **Step 5: Consolidate documentation**

Invoke `bluequartz-skills:filter-documentation`. Expand
`docs/ComputeODFFilter.md` with selected-phase semantics, units metadata path,
combined presentation parameters, optional output modes, exact unit labels,
preflight behavior, and final figure. Remove the standalone filter document.

Caption the figure with dataset, phase, bin size, MUD rendering, section angles,
and `Rainbow Blended White`.

- [ ] **Step 6: Verify and commit docs**

```bash
git diff --check -- docs/ComputeODFFilter.md docs/WriteODFSectionsFilter.md
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R 'MTRSim::ComputeODFFilter' -j 1 --output-on-failure --timeout 1500
git add docs/ComputeODFFilter.md docs/WriteODFSectionsFilter.md \
  docs/Images/odf_phi2_sections.png
git commit -s -m "DOC: Document combined ODF presentation"
```

---

### Task 5: Cross-repository completion verification

**Files:**
- Verification only; no source mutation expected.

**Interfaces:**
- Consumes: all prior tasks.
- Produces: final evidence and release-gate handoff.

- [ ] **Step 1: Verify branch scope and protected paths**

```bash
git -C /Users/mjackson/Workspace7/EbsdLib status --short
git -C /Users/mjackson/Workspace7/EbsdLib log --oneline --decorate -8
git -C /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim status --short
git -C /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim log --oneline --decorate -12
```

Expected: EbsdLib contains the annotation commit; MTRSim contains the phase/unit,
combined-filter, and documentation commits. Only the known protected untracked
MTRSim paths remain.

- [ ] **Step 2: Run complete EbsdLib verification**

```bash
cd /Users/mjackson/Workspace7/EbsdLib
cmake --preset ebsdlib-release
cmake --build --preset "ebsdlib-release [BUILD]" -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/ebsdlib-release \
  -j 1 --output-on-failure
```

- [ ] **Step 3: Run full local-EbsdLib DREAM3D-NX build and MTRSim tests**

```bash
cd /Users/mjackson/Workspace7/DREAM3DNX
cmake --preset NX-Com-Qt69-Vtk96-Rel-EbsdLib
env LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 \
  cmake --build --preset "NX-Com-Qt69-Vtk96-Rel-EbsdLib [BUILD]" -j
ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk96-Rel-EbsdLib \
  -R '^MTRSim::' -j 1 --output-on-failure --timeout 1500
```

- [ ] **Step 4: Rerun external pipeline and validators from scratch**

Run Task 4 preflight, execute, structural validator, and MTEX comparison again.
Require zero-error exits, no stale standalone-filter step, correct unit metadata,
PNG/RGB equality, correct maximum, and expected correlations.

- [ ] **Step 5: Run final searches and format gates**

```bash
rg -n "WriteODFSections(Filter)?" /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim \
  -g '!docs/architecture/**' -g '!.superpowers/**'
rg -n "Deferred until execution" /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim \
  -g '!docs/architecture/**' -g '!.superpowers/**'
git -C /Users/mjackson/Workspace7/EbsdLib diff --check upstream/develop...HEAD
git -C /Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim diff --check
```

Both searches must return no product/docs/pipeline matches. Run pinned
clang-format `--dry-run --Werror` over every touched C/C++ file in both repos.

- [ ] **Step 6: Record the release boundary**

The normal `NX-Com-Qt69-Vtk96-Rel` preset remains gated until an EbsdLib release
contains the final annotation/compositor API and the dependency is updated.
Record the exact tested EbsdLib commit and do not change vcpkg in this plan.
