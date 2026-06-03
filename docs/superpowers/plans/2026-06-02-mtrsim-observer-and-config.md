# MTRSim Observer + Config-File Input — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a progress/cancellation observer (`mtrsim::ISimulationObserver`) threaded through `simulateMTR`/`PGRFSimulation`/`ODFSampler`, and an optional MTRSim config-JSON input on `MTRSimFilter` gated by a linked boolean.

**Architecture:** A header-only pure-virtual `ISimulationObserver` (progress + cancel) is passed (nullable) into the library simulation entry points; `nullptr` preserves today's behavior exactly. The standalone CLI and the filter each supply a concrete observer. A shared `mtrsim::parseConfigJson` parser (extracted from `main.cpp`) lets both the CLI and the filter build `SimulationParams` from a config file; the filter exposes it via a linked boolean that shows config-mode vs manual-mode parameters.

**Tech Stack:** C++17, Eigen, simplnx filter framework, Catch2 v2 (bare `Approx`, not `Catch::Approx`), nlohmann::json, spdlog, CMake.

---

## Background the engineer needs

**Reference design:** `docs/superpowers/specs/2026-06-02-mtrsim-observer-and-config-design.md`.

**Current signatures (verified):**
- `mtrsim::simulateMTR(const SimulationParams&, const std::vector<ODFComponent>&, std::mt19937_64&, int n1, int nPHI, int n2)` → `MTRSimResult` (`src/LibMTRSim/MTRSimDriver.hpp`). `MTRSimResult` has `int nx,ny,nz; std::vector<int32_t> mtrIndex; std::vector<double> phi1, phi, phi2;`.
- `mtrsim::ODFSampler::sampleN(int n, const ODFComponent& component, const ODFComponent& uniform)` → `Eigen::MatrixXd (n×3)`. Two `for i in [0,n)` loops inside (sampling at `ODFSampler.cpp:56-64`, assignment at `77-83`).
- `mtrsim::PGRFSimulation::run(const SimulationParams&)` → `PGRFResult` (`src/LibMTRSim/PGRFSimulation.hpp`); loops over `numGaussians = numComponents-1` latent fields.
- `SimulationParams` fields: `double xLen,yLen,zLen,dx,dy,dz; std::vector<double> volumeFractions; std::vector<std::vector<double>> thetaList; std::vector<double> nuggetVariance; std::string odfInputPath, outputDir; uint64_t seed;`.
- `main.cpp` parses config JSON inline at lines 153-179.

**Filter (`src/MTRSim/Filters/`):** `MTRSimFilter` keys end with `_path`/`_index` per `FilterValidationTest`. The algorithm `MTRSim` holds `const IFilter::MessageHandler& m_MessageHandler` and `const std::atomic_bool& m_ShouldCancel`. The standard simplnx seed pattern (`use_seed`/`seed_value`/`seed_array_name`) and `linkParameters` are already used.

**Build & test commands:**
- **Standalone library** (Phase A, B tests): configure once `cmake --preset mtrsim-Rel`; build `cmake --build /Users/mjackson/Workspace7/Build/mtrsim-Rel`; run tests `/Users/mjackson/Workspace7/Build/mtrsim-Rel/bin/mtrsim_tests "[tag]"` or `ctest --test-dir /Users/mjackson/Workspace7/Build/mtrsim-Rel --output-on-failure`. Catch2 uses **bare `Approx`**.
- **Plugin** (Phase C): build `cmake --build /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk95-Rel-EbsdLib --target MTRSimUnitTest`; run `ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk95-Rel-EbsdLib -R "MTRSim::" --output-on-failure`.
- **FilterValidationTest** (part of `simplnx_test`): build `cmake --build <plugin-build> --target simplnx_test`; run `ctest --test-dir <plugin-build> -R "FilterValidation" --output-on-failure` (or run `Bin/simplnx_test "[FilterValidation]"` — confirm the exact ctest name with `ctest -N | grep -i filtervalidation`).
- **Debug smoke** (Phase D): `cmake --build /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk95-Dbg-EbsdLib --target MTRSim nxrunner` then `<dbg>/Bin/nxrunner_d --execute pipelines/MTRSim_smoke_test.d3dpipeline`.

**No git worktree** (DREAM3DNX build expects MTRSim at a fixed path). Work on the current branch. Sync first: `git fetch upstream && git reset --soft upstream/topic/create_mtr_sim_filter` if behind (it changes no files when trees match).

---

## File Structure

**Create:**
- `src/LibMTRSim/ISimulationObserver.hpp` — the interface.
- `src/LibMTRSim/SimulationObservers.hpp` — `NullObserver`, `ConsoleObserver`.
- `src/LibMTRSim/ConfigIO.hpp`, `src/LibMTRSim/ConfigIO.cpp` — `parseConfigJson`.
- `tests/test_observer.cpp`, `tests/test_config_io.cpp`.

**Modify:**
- `src/LibMTRSim/MTRSimDriver.{hpp,cpp}` — observer arg, `MTRSimResult::cancelled`.
- `src/LibMTRSim/PGRFSimulation.{hpp,cpp}` — observer arg.
- `src/LibMTRSim/ODFSampler.{hpp,cpp}` — observer arg + in-loop cancel/progress.
- `src/app/main.cpp` — `ConsoleObserver`, `parseConfigJson`.
- `src/MTRSim/Filters/MTRSimFilter.{hpp,cpp}` — config params + linking + config-mode preflight/execute.
- `src/MTRSim/Filters/Algorithms/MTRSim.{hpp,cpp}` — `FilterObserver`, config-mode params, pass observer.
- `src/LibMTRSim/CMakeLists.txt`, `MTRSimPlugin.cmake`, `tests/CMakeLists.txt` — register.
- `docs/MTRSimFilter.md`.

---

## Phase A — Observer (LibMTRSim)

### Task A1: `ISimulationObserver` + default observers

**Files:** Create `src/LibMTRSim/ISimulationObserver.hpp`, `src/LibMTRSim/SimulationObservers.hpp`, `tests/test_observer.cpp`; Modify `src/LibMTRSim/CMakeLists.txt`, `MTRSimPlugin.cmake`, `tests/CMakeLists.txt`.

- [ ] **Step 1: Create `src/LibMTRSim/ISimulationObserver.hpp`**

```cpp
#pragma once

#include "libmtrsim_export.h"

#include <cstdint>
#include <string>

namespace mtrsim {

/**
 * @brief Observer interface for long-running MTR simulations.
 *
 * Implementations receive throttled progress updates and are polled for
 * cancellation from the thread running simulateMTR. They must be cheap and
 * must not throw.
 */
class LIBMTRSIM_EXPORT ISimulationObserver {
public:
  ISimulationObserver() = default;
  virtual ~ISimulationObserver() = default;

  ISimulationObserver(const ISimulationObserver&) = delete;
  ISimulationObserver& operator=(const ISimulationObserver&) = delete;
  ISimulationObserver(ISimulationObserver&&) = delete;
  ISimulationObserver& operator=(ISimulationObserver&&) = delete;

  /// Report progress. `done`/`total` describe the current phase; `message`
  /// names it. `total <= 0` means "indeterminate".
  virtual void updateProgress(int64_t done, int64_t total, const std::string& message) = 0;

  /// Polled at checkpoints; returning true stops the simulation early.
  [[nodiscard]] virtual bool shouldCancel() const = 0;
};

} // namespace mtrsim
```

- [ ] **Step 2: Create `src/LibMTRSim/SimulationObservers.hpp`**

```cpp
#pragma once

#include "ISimulationObserver.hpp"

#include <spdlog/spdlog.h>

#include <cstdint>
#include <string>

namespace mtrsim {

/// No-op observer; the implicit default when none is supplied.
class NullObserver : public ISimulationObserver {
public:
  void updateProgress(int64_t /*done*/, int64_t /*total*/, const std::string& /*message*/) override {}
  [[nodiscard]] bool shouldCancel() const override { return false; }
};

/// Logs progress via spdlog (throttled). Never cancels.
class ConsoleObserver : public ISimulationObserver {
public:
  void updateProgress(int64_t done, int64_t total, const std::string& message) override {
    const int pct = (total > 0) ? static_cast<int>(done * 100 / total) : -1;
    if (pct != m_LastPct) {
      m_LastPct = pct;
      if (pct >= 0) {
        spdlog::info("[{:3d}%] {}", pct, message);
      } else {
        spdlog::info("{}", message);
      }
    }
  }
  [[nodiscard]] bool shouldCancel() const override { return false; }

private:
  int m_LastPct = -2;
};

} // namespace mtrsim
```

- [ ] **Step 3: Register headers in CMake**

In `src/LibMTRSim/CMakeLists.txt`, add `ISimulationObserver.hpp` and `SimulationObservers.hpp` to the `MTRSIM_HEADERS` list (mirror the existing `ODFSampler.hpp` entry). In `MTRSimPlugin.cmake`, add both to the plugin header list (mirror the `MTRSimDriver.hpp` entry).

- [ ] **Step 4: Write the failing test**

Create `tests/test_observer.cpp`:

```cpp
#include "SimulationObservers.hpp"

#include <catch2/catch.hpp>

#include <vector>

namespace {
// Test double: records progress calls; can be told to cancel after K updates.
class RecordingObserver : public mtrsim::ISimulationObserver {
public:
  explicit RecordingObserver(int cancelAfter = -1) : m_CancelAfter(cancelAfter) {}
  void updateProgress(int64_t done, int64_t total, const std::string& message) override {
    calls.push_back({done, total, message});
  }
  bool shouldCancel() const override {
    return m_CancelAfter >= 0 && static_cast<int>(calls.size()) >= m_CancelAfter;
  }
  struct Call { int64_t done; int64_t total; std::string message; };
  std::vector<Call> calls;
private:
  int m_CancelAfter;
};
} // namespace

TEST_CASE("NullObserver never cancels and ignores progress", "[observer]") {
  mtrsim::NullObserver obs;
  obs.updateProgress(1, 10, "x");
  REQUIRE_FALSE(obs.shouldCancel());
}

TEST_CASE("RecordingObserver records and cancels after K", "[observer]") {
  RecordingObserver obs(2);
  REQUIRE_FALSE(obs.shouldCancel());
  obs.updateProgress(1, 10, "a");
  REQUIRE_FALSE(obs.shouldCancel());
  obs.updateProgress(2, 10, "b");
  REQUIRE(obs.shouldCancel());
  REQUIRE(obs.calls.size() == 2);
  REQUIRE(obs.calls[1].done == 2);
}
```

Add `test_observer.cpp` to the `mtrsim_tests` executable in `tests/CMakeLists.txt` (mirror `test_odf_sampler.cpp`). Add include dir if needed (the tests target already adds `${PROJECT_SOURCE_DIR}/src/LibMTRSim`, so `#include "SimulationObservers.hpp"` resolves).

- [ ] **Step 5: Build and run; verify pass**

Run: `cmake --build /Users/mjackson/Workspace7/Build/mtrsim-Rel && /Users/mjackson/Workspace7/Build/mtrsim-Rel/bin/mtrsim_tests "[observer]"`
Expected: PASS (3 assertions / 2 cases).

- [ ] **Step 6: Commit**

```bash
git add src/LibMTRSim/ISimulationObserver.hpp src/LibMTRSim/SimulationObservers.hpp \
        src/LibMTRSim/CMakeLists.txt MTRSimPlugin.cmake tests/test_observer.cpp tests/CMakeLists.txt
git commit -m "feat(lib): add ISimulationObserver + Null/Console observers"
```
(End every commit body with: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`)

---

### Task A2: thread observer through `simulateMTR` + `PGRFSimulation` (coarse: stages + per-component)

**Files:** Modify `src/LibMTRSim/MTRSimDriver.{hpp,cpp}`, `src/LibMTRSim/PGRFSimulation.{hpp,cpp}`, `tests/test_mtrsim_driver.cpp`.

- [ ] **Step 1: Add `cancelled` to `MTRSimResult` and the observer param to `simulateMTR` (header)**

In `src/LibMTRSim/MTRSimDriver.hpp`: add `#include "ISimulationObserver.hpp"`. In `struct MTRSimResult` add `bool cancelled = false;`. Change the declaration:

```cpp
LIBMTRSIM_EXPORT MTRSimResult simulateMTR(const SimulationParams& params, const std::vector<ODFComponent>& odfComponents, std::mt19937_64& rng, int n1, int nPHI, int n2, ISimulationObserver* observer = nullptr);
```

- [ ] **Step 2: Add observer param to `PGRFSimulation::run` (header)**

In `src/LibMTRSim/PGRFSimulation.hpp`: add `#include "ISimulationObserver.hpp"` and change:

```cpp
PGRFResult run(const SimulationParams& params, ISimulationObserver* observer = nullptr);
```

- [ ] **Step 3: Thread cancel/progress into `PGRFSimulation::run` (impl)**

In `src/LibMTRSim/PGRFSimulation.cpp`, update the signature and, inside the latent-field loop (`for (int h = 0; h < numGaussians; ++h)`), before generating each field add:

```cpp
    if (observer != nullptr) {
      if (observer->shouldCancel()) {
        return PGRFResult{};  // empty; caller checks observer->shouldCancel()
      }
      observer->updateProgress(h, numGaussians, fmt::format("Simulating latent Gaussian field {}/{}", h + 1, numGaussians));
    }
```

Add `#include <fmt/format.h>` if not present. (Leave the existing spdlog logs.)

- [ ] **Step 4: Thread observer into `simulateMTR` (impl)** — update `src/LibMTRSim/MTRSimDriver.cpp`'s `simulateMTR` to accept `ISimulationObserver* observer`, pass it to `pgrf.run(params, observer)`, and add checkpoints. Use a tiny local helper for null-safety:

```cpp
MTRSimResult simulateMTR(const SimulationParams& params, const std::vector<ODFComponent>& odfComponents, std::mt19937_64& rng, int n1, int nPHI, int n2, ISimulationObserver* observer) {
  const int nx = static_cast<int>(std::round(params.xLen / params.dx));
  const int ny = static_cast<int>(std::round(params.yLen / params.dy));
  const int nz = std::max(static_cast<int>(std::round(params.zLen / params.dz)), 1);
  const int N = nx * ny * nz;

  if (static_cast<int>(odfComponents.size()) != static_cast<int>(params.volumeFractions.size())) {
    throw std::invalid_argument("simulateMTR: odfComponents count must equal volumeFractions count");
  }

  auto cancelled = [&]() { return observer != nullptr && observer->shouldCancel(); };
  auto report = [&](int64_t done, int64_t total, const std::string& msg) {
    if (observer != nullptr) { observer->updateProgress(done, total, msg); }
  };

  // 1. PGRF assignment.
  report(0, 100, "Running plurigaussian field simulation");
  PGRFSimulation pgrf{rng};
  const PGRFResult pgrf_result = pgrf.run(params, observer);
  if (cancelled()) { MTRSimResult out; out.cancelled = true; return out; }
  if (static_cast<int>(pgrf_result.mtrIndex.size()) != N) {
    throw std::runtime_error("simulateMTR: PGRF result size does not match grid dimensions");
  }

  // 2. Sample N orientations per component.
  const ODFComponent uniformOdf = buildUniformODF(n1, nPHI, n2);
  const int numComponents = static_cast<int>(odfComponents.size());
  std::vector<Eigen::MatrixXd> orientSamples(static_cast<std::size_t>(numComponents));
  ODFSampler sampler{rng};
  for (int j = 0; j < numComponents; ++j) {
    report(j, numComponents, fmt::format("Sampling orientations (component {}/{})", j + 1, numComponents));
    orientSamples[static_cast<std::size_t>(j)] = sampler.sampleN(N, odfComponents[static_cast<std::size_t>(j)], uniformOdf, observer);
    if (cancelled()) { MTRSimResult out; out.cancelled = true; return out; }
  }

  // 3. Per-voxel assignment (sim order).
  std::vector<int32_t> mtrSim(N);
  std::vector<double> phi1Sim(N), phiSim(N), phi2Sim(N);
  for (int i = 0; i < N; ++i) {
    const int comp = pgrf_result.mtrIndex[i] - 1;
    mtrSim[i] = pgrf_result.mtrIndex[i];
    phi1Sim[i] = orientSamples[static_cast<std::size_t>(comp)](i, 0);
    phiSim[i] = orientSamples[static_cast<std::size_t>(comp)](i, 1);
    phi2Sim[i] = orientSamples[static_cast<std::size_t>(comp)](i, 2);
  }

  // 4. Remap to SIMPLNX z,y,x order.
  report(100, 100, "Finalizing microstructure");
  MTRSimResult out;
  out.nx = nx; out.ny = ny; out.nz = nz;
  out.mtrIndex = remapSimToZYX(mtrSim, nx, ny, nz);
  out.phi1 = remapSimToZYX(phi1Sim, nx, ny, nz);
  out.phi = remapSimToZYX(phiSim, nx, ny, nz);
  out.phi2 = remapSimToZYX(phi2Sim, nx, ny, nz);
  return out;
}
```

Ensure `#include "ISimulationObserver.hpp"` and `#include <fmt/format.h>` are present in the .cpp. The `sampleN(..., observer)` 4-arg form is added in Task A3; until then call it with the existing 3-arg form and add `observer` in A3 — OR implement A3 first. **Implement A3 before building A2 to completion** (the `sampleN` overload is needed).

- [ ] **Step 5: Add the cancel test** to `tests/test_mtrsim_driver.cpp` (reuse the existing `[mtrsim_driver][statistical]` setup; add a `RecordingObserver`-style local class or include from a shared header). Add:

```cpp
#include "ISimulationObserver.hpp"

namespace {
class CancelAfterObserver : public mtrsim::ISimulationObserver {
public:
  explicit CancelAfterObserver(int k) : m_K(k) {}
  void updateProgress(int64_t, int64_t, const std::string&) override { ++m_Count; }
  bool shouldCancel() const override { return m_Count >= m_K; }
  int count() const { return m_Count; }
private:
  int m_K; mutable int m_Count = 0;
};
}

TEST_CASE("simulateMTR cancels early when observer requests it", "[mtrsim_driver]") {
  mtrsim::SimulationParams params;
  params.xLen = 6.0; params.yLen = 6.0; params.zLen = 0.0;
  params.dx = 0.02; params.dy = 0.02; params.dz = 0.02;
  params.volumeFractions = {0.30, 0.35, 0.35};
  params.thetaList = {{0.10, 0.45, 0.10}, {0.08, 0.37, 0.08}};
  params.seed = 42;
  std::vector<mtrsim::ODFComponent> comps = {
      mtrsim::buildUniformODF(72, 36, 72), mtrsim::buildUniformODF(72, 36, 72), mtrsim::buildUniformODF(72, 36, 72)};
  std::mt19937_64 rng(params.seed);
  CancelAfterObserver obs(1); // cancel at the first progress checkpoint
  const mtrsim::MTRSimResult r = mtrsim::simulateMTR(params, comps, rng, 72, 36, 72, &obs);
  REQUIRE(r.cancelled);
  REQUIRE(r.mtrIndex.empty()); // no full output produced
}

TEST_CASE("simulateMTR with nullptr observer is unaffected", "[mtrsim_driver]") {
  mtrsim::SimulationParams params;
  params.xLen = 2.0; params.yLen = 2.0; params.zLen = 0.0;
  params.dx = 0.02; params.dy = 0.02; params.dz = 0.02;
  params.volumeFractions = {0.30, 0.35, 0.35};
  params.thetaList = {{0.10, 0.45, 0.10}, {0.08, 0.37, 0.08}};
  std::mt19937_64 rng(7);
  const mtrsim::MTRSimResult r = mtrsim::simulateMTR(params, comps_or_uniform(), rng, 72, 36, 72);
  // build comps inline:
  REQUIRE_FALSE(r.cancelled);
}
```

Replace `comps_or_uniform()` with three `mtrsim::buildUniformODF(72,36,72)` entries built locally (repeat the vector). Keep both tests self-contained.

- [ ] **Step 6: Build, run `[mtrsim_driver]` + `[observer]`, verify pass; commit**

```bash
git add src/LibMTRSim/MTRSimDriver.hpp src/LibMTRSim/MTRSimDriver.cpp \
        src/LibMTRSim/PGRFSimulation.hpp src/LibMTRSim/PGRFSimulation.cpp tests/test_mtrsim_driver.cpp
git commit -m "feat(lib): thread ISimulationObserver through simulateMTR + PGRFSimulation"
```

---

### Task A3: deep cancel/progress in `ODFSampler::sampleN`

**Files:** Modify `src/LibMTRSim/ODFSampler.{hpp,cpp}`, `tests/test_odf_sampler.cpp`.

- [ ] **Step 1: Add observer param (header)** — in `src/LibMTRSim/ODFSampler.hpp` add `#include "ISimulationObserver.hpp"` and change `sampleN`:

```cpp
Eigen::MatrixXd sampleN(int n, const ODFComponent& component, const ODFComponent& uniform, ISimulationObserver* observer = nullptr);
```

- [ ] **Step 2: Add in-loop cancel checks (impl)** — in `src/LibMTRSim/ODFSampler.cpp`, update the signature and add a throttled cancel check inside the sampling loop (`for i in [0,n)` at line 56) and the assignment loop (line 77). Check every 4096 iterations to keep it cheap; on cancel, return an empty matrix (the caller treats an empty/short result as cancellation via `observer->shouldCancel()`):

```cpp
  constexpr int kCheck = 4096;
  for (int i = 0; i < n; ++i) {
    if (observer != nullptr && (i % kCheck) == 0 && observer->shouldCancel()) {
      return Eigen::MatrixXd(0, 3);
    }
    // ... existing sampling body ...
  }
```
Apply the same `if (observer != nullptr && (i % kCheck) == 0 && observer->shouldCancel()) { return Eigen::MatrixXd(0, 3); }` guard at the top of the assignment loop body. Do NOT change the RNG draw order on the non-cancel path (the modulo check must not consume RNG), so the regression test stays bit-stable.

> Note: `simulateMTR` (Task A2) calls `observer->shouldCancel()` right after each `sampleN` and sets `cancelled`. A short/empty matrix from a cancelled `sampleN` is therefore never consumed.

- [ ] **Step 3: Add a test** to `tests/test_odf_sampler.cpp`:

```cpp
#include "ISimulationObserver.hpp"

namespace {
class ImmediateCancel : public mtrsim::ISimulationObserver {
public:
  void updateProgress(int64_t, int64_t, const std::string&) override {}
  bool shouldCancel() const override { return true; }
};
}

TEST_CASE("sampleN bails out promptly when observer cancels", "[odf_sampler]") {
  mtrsim::ODFComponent uni = mtrsim::buildUniformODF_or_local(); // see note
  // Build a uniform component directly if buildUniformODF isn't linked here.
  std::mt19937_64 rng(1);
  mtrsim::ODFSampler sampler{rng};
  ImmediateCancel cancel;
  Eigen::MatrixXd out = sampler.sampleN(100000, uni, uni, &cancel);
  REQUIRE(out.rows() == 0);
}
```

If `buildUniformODF` is not available to this test translation unit, construct a small valid `ODFComponent` inline: set `odfVal`, `phi1Bins`, `phiBins`, `phi2Bins` each to a `Eigen::VectorXd` of equal length (e.g. 10) with positive `odfVal` and increasing bin centres. Replace `buildUniformODF_or_local()` accordingly. (Check whether `test_odf_sampler.cpp` already builds a component it can reuse.)

- [ ] **Step 4: Build, run `[odf_sampler]` + `[mtrsim_driver]` (regression), verify pass**

Run: `cmake --build /Users/mjackson/Workspace7/Build/mtrsim-Rel && /Users/mjackson/Workspace7/Build/mtrsim-Rel/bin/mtrsim_tests "[odf_sampler],[mtrsim_driver]"`
Expected: PASS, including the pre-existing statistical test (proves the cancel-check modulo didn't perturb the RNG path).

- [ ] **Step 5: Commit**

```bash
git add src/LibMTRSim/ODFSampler.hpp src/LibMTRSim/ODFSampler.cpp tests/test_odf_sampler.cpp
git commit -m "feat(lib): in-loop cancel checks in ODFSampler::sampleN"
```

---

### Task A4: standalone CLI uses `ConsoleObserver`

**Files:** Modify `src/app/main.cpp`.

- [ ] **Step 1: Pass a ConsoleObserver into simulateMTR** — in `src/app/main.cpp`, add `#include "SimulationObservers.hpp"`, and change the `simulateMTR(...)` call to:

```cpp
  mtrsim::ConsoleObserver observer;
  mtrsim::MTRSimResult sim = mtrsim::simulateMTR(params, odfComponents, rng, k_OdfBinsPhi1, k_OdfBinsPHI, k_OdfBinsPhi2, &observer);
```
(Use the existing `k_OdfBins*` constants. If `sim.cancelled` is true, log a warning and skip the CSV/PNG; the CLI never cancels but handle it for completeness.)

- [ ] **Step 2: Build the app, run it on the smoke config, verify clean**

Run: `cmake --build /Users/mjackson/Workspace7/Build/mtrsim-Rel` then run the `MTRsim` binary: `/Users/mjackson/Workspace7/Build/mtrsim-Rel/bin/MTRsim -c configs/smoke_test.json -o /tmp/mtr_a4 --seed 42` (create `/tmp/mtr_a4` first).
Expected: exit 0, progress `[ NN%]` lines logged, CSV/PNG written.

- [ ] **Step 3: Commit**

```bash
git add src/app/main.cpp
git commit -m "feat(app): standalone driver reports progress via ConsoleObserver"
```

---

## Phase B — Shared config parser (LibMTRSim)

### Task B1: `parseConfigJson`

**Files:** Create `src/LibMTRSim/ConfigIO.hpp`, `src/LibMTRSim/ConfigIO.cpp`, `tests/test_config_io.cpp`; Modify `src/LibMTRSim/CMakeLists.txt`, `MTRSimPlugin.cmake`, `tests/CMakeLists.txt`.

- [ ] **Step 1: Create `src/LibMTRSim/ConfigIO.hpp`**

```cpp
#pragma once

#include "SimulationParams.hpp"
#include "libmtrsim_export.h"

#include <filesystem>

namespace mtrsim {

/**
 * @brief Parse an MTRSim config JSON into a SimulationParams.
 *
 * Recognized keys: xLen,yLen,zLen, dx,dy,dz, volumeFractions, thetaList, seed.
 * Unknown keys (odfInputPath, nuggetVariance, comments) are ignored. Fields
 * absent from the JSON keep their SimulationParams defaults.
 *
 * @throws std::runtime_error if the file cannot be opened or the JSON is invalid.
 */
LIBMTRSIM_EXPORT SimulationParams parseConfigJson(const std::filesystem::path& path);

} // namespace mtrsim
```

- [ ] **Step 2: Create `src/LibMTRSim/ConfigIO.cpp`** (logic lifted from `main.cpp:153-179`)

```cpp
#include "ConfigIO.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <stdexcept>

namespace mtrsim {

SimulationParams parseConfigJson(const std::filesystem::path& path) {
  std::ifstream f(path);
  if (!f.is_open()) {
    throw std::runtime_error("parseConfigJson: cannot open config file: " + path.string());
  }
  SimulationParams params;
  try {
    const nlohmann::json j = nlohmann::json::parse(f);
    if (j.contains("xLen")) params.xLen = j["xLen"].get<double>();
    if (j.contains("yLen")) params.yLen = j["yLen"].get<double>();
    if (j.contains("zLen")) params.zLen = j["zLen"].get<double>();
    if (j.contains("dx")) params.dx = j["dx"].get<double>();
    if (j.contains("dy")) params.dy = j["dy"].get<double>();
    if (j.contains("dz")) params.dz = j["dz"].get<double>();
    if (j.contains("volumeFractions")) params.volumeFractions = j["volumeFractions"].get<std::vector<double>>();
    if (j.contains("thetaList")) params.thetaList = j["thetaList"].get<std::vector<std::vector<double>>>();
    if (j.contains("nuggetVariance")) params.nuggetVariance = j["nuggetVariance"].get<std::vector<double>>();
    if (j.contains("seed")) params.seed = j["seed"].get<uint64_t>();
  } catch (const nlohmann::json::exception& e) {
    throw std::runtime_error(std::string("parseConfigJson: invalid JSON: ") + e.what());
  }
  return params;
}

} // namespace mtrsim
```

> `nuggetVariance` is parsed for fidelity but unused by the simulation; `odfInputPath` is intentionally NOT applied (the filter/caller supplies the ODF separately).

- [ ] **Step 3: Register in CMake** — add `ConfigIO.cpp`/`.hpp` to `src/LibMTRSim/CMakeLists.txt` (source + header lists) and `MTRSimPlugin.cmake` (the `${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/ConfigIO.cpp` source + `.hpp` header).

- [ ] **Step 4: Write the failing test** — create `tests/test_config_io.cpp`:

```cpp
#include "ConfigIO.hpp"

#include <catch2/catch.hpp>

#include <cstdio>
#include <fstream>
#include <string>

namespace {
std::string writeTemp(const std::string& contents) {
  const std::string path = std::string(MTRSIM_TEST_DATA_DIR) + "/_tmp_config_io.json";
  std::ofstream o(path); o << contents; o.close();
  return path;
}
}

TEST_CASE("parseConfigJson reads known fields and ignores extras", "[config_io]") {
  const std::string path = writeTemp(R"({
    "xLen": 38.1, "yLen": 12.7, "zLen": 0.0,
    "dx": 0.02, "dy": 0.02, "dz": 0.02,
    "volumeFractions": [0.30, 0.35, 0.35],
    "thetaList": [[0.10,0.45,0.10],[0.08,0.37,0.08]],
    "odfInputPath": "ignored.h5", "nuggetVariance": [0.6,0.7,0.7], "seed": 99
  })");
  const mtrsim::SimulationParams p = mtrsim::parseConfigJson(path);
  REQUIRE(p.xLen == Approx(38.1));
  REQUIRE(p.dx == Approx(0.02));
  REQUIRE(p.volumeFractions.size() == 3);
  REQUIRE(p.volumeFractions[1] == Approx(0.35));
  REQUIRE(p.thetaList.size() == 2);
  REQUIRE(p.thetaList[0][1] == Approx(0.45));
  REQUIRE(p.seed == 99);
  std::remove(path.c_str());
}

TEST_CASE("parseConfigJson throws on missing file", "[config_io]") {
  REQUIRE_THROWS_AS(mtrsim::parseConfigJson("/nonexistent/path/nope.json"), std::runtime_error);
}

TEST_CASE("parseConfigJson throws on malformed JSON", "[config_io]") {
  const std::string path = writeTemp("{ this is not json");
  REQUIRE_THROWS_AS(mtrsim::parseConfigJson(path), std::runtime_error);
  std::remove(path.c_str());
}
```

Add `test_config_io.cpp` to `tests/CMakeLists.txt`. (`MTRSIM_TEST_DATA_DIR` is already defined for the test target.)

- [ ] **Step 5: Build, run `[config_io]`, verify FAIL→PASS; commit**

```bash
git add src/LibMTRSim/ConfigIO.hpp src/LibMTRSim/ConfigIO.cpp \
        src/LibMTRSim/CMakeLists.txt MTRSimPlugin.cmake tests/test_config_io.cpp tests/CMakeLists.txt
git commit -m "feat(lib): add shared parseConfigJson config reader"
```

---

### Task B2: standalone CLI uses `parseConfigJson`

**Files:** Modify `src/app/main.cpp`.

- [ ] **Step 1: Replace the inline JSON parse** — in `src/app/main.cpp`, add `#include "ConfigIO.hpp"`. Replace the inline block (`main.cpp:153-179`) that parses the config into `params` with:

```cpp
    try {
      mtrsim::SimulationParams cfg = mtrsim::parseConfigJson(configPath);
      cfg.outputDir = params.outputDir;            // keep CLI-provided output dir
      if (seed != 0) { cfg.seed = seed; }          // CLI --seed overrides JSON seed
      params = cfg;
    } catch (const std::exception& e) {
      spdlog::error("{}", e.what());
      return 1;
    }
```
Preserve the existing CLI semantics: `--seed` (non-zero) overrides the JSON seed; otherwise the JSON seed is used. Remove the now-unused `nlohmann/json.hpp` include from `main.cpp` only if nothing else there needs it (the HDF5 ODF loader doesn't); otherwise leave it.

- [ ] **Step 2: Build + run on default and smoke configs; verify identical behavior**

Run the app on `configs/smoke_test.json` and confirm exit 0 and the same grid/output as before this change.

- [ ] **Step 3: Commit**

```bash
git add src/app/main.cpp
git commit -m "refactor(app): standalone driver uses shared parseConfigJson"
```

---

## Phase C — Filter (config-file mode + filter observer)

### Task C1: `FilterObserver` in the algorithm; pass observer to `simulateMTR`

**Files:** Modify `src/MTRSim/Filters/Algorithms/MTRSim.{hpp,cpp}`.

- [ ] **Step 1: Declare a private FilterObserver** — in `src/MTRSim/Filters/Algorithms/MTRSim.cpp` (anonymous namespace at top, after includes), add an observer adapting to the message handler + cancel flag. Add `#include "LibMTRSim/ISimulationObserver.hpp"`:

```cpp
namespace {
class FilterObserver : public mtrsim::ISimulationObserver {
public:
  FilterObserver(const nx::core::IFilter::MessageHandler& mh, const std::atomic_bool& cancel)
  : m_MessageHandler(mh), m_ShouldCancel(cancel) {}
  void updateProgress(int64_t done, int64_t total, const std::string& message) override {
    const int32_t pct = (total > 0) ? static_cast<int32_t>(done * 100 / total) : 0;
    m_MessageHandler(nx::core::IFilter::ProgressMessage{nx::core::IFilter::Message::Type::Progress, message, pct});
  }
  bool shouldCancel() const override { return m_ShouldCancel.load(); }
private:
  const nx::core::IFilter::MessageHandler& m_MessageHandler;
  const std::atomic_bool& m_ShouldCancel;
};
}
```
> Verify the exact `ProgressMessage` construction against a simplnx filter that emits progress (grep `ProgressMessage` in `/Users/mjackson/Workspace7/simplnx/src/Plugins`); adjust the struct/arg form to match the real API. See the `progress-messaging` skill conventions.

- [ ] **Step 2: Use the observer in `operator()`** — in `MTRSim::operator()`, construct `FilterObserver observer{m_MessageHandler, m_ShouldCancel};` and pass `&observer` to `mtrsim::simulateMTR(params, components, rng, n1, nPHI, n2, &observer);`. After the call, replace the existing `if (m_ShouldCancel) { return {}; }` so it also covers `sim.cancelled`:

```cpp
  if (m_ShouldCancel || sim.cancelled) { return {}; }
```

- [ ] **Step 3: Build plugin + run MTRSim tests; verify pass; commit**

Run: `cmake --build /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk95-Rel-EbsdLib --target MTRSimUnitTest && ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk95-Rel-EbsdLib -R "MTRSim::" --output-on-failure`
Expected: all MTRSim tests PASS.

```bash
git add src/MTRSim/Filters/Algorithms/MTRSim.hpp src/MTRSim/Filters/Algorithms/MTRSim.cpp
git commit -m "feat(filter): report progress + honor cancel via FilterObserver"
```

---

### Task C2: config-file parameters + linked show/hide

**Files:** Modify `src/MTRSim/Filters/MTRSimFilter.{hpp,cpp}`.

- [ ] **Step 1: Add parameter keys (header)** — in `src/MTRSim/Filters/MTRSimFilter.hpp` add:

```cpp
  static constexpr StringLiteral k_UseConfigFile_Key = "use_config_file";
  static constexpr StringLiteral k_ConfigFilePath_Key = "config_file_path";
```

- [ ] **Step 2: Add the parameters + linking (parameters())** — in `MTRSimFilter.cpp` `parameters()`, add includes `#include "simplnx/Parameters/FileSystemPathParameter.hpp"` (already used by Read filter) and, in a new separator before the simulation parameters:

```cpp
  params.insertSeparator(Parameters::Separator{"Configuration Source"});
  params.insertLinkableParameter(std::make_unique<BoolParameter>(k_UseConfigFile_Key, "Load Simulation Parameters from Config File",
      "When ON, read volume fractions, theta list, physical size/spacing, and seed from an MTRSim JSON config file instead of the fields below.", false));
  params.insert(std::make_unique<FileSystemPathParameter>(k_ConfigFilePath_Key, "MTRSim Config File (JSON)",
      "MTRSim configuration JSON (same schema as the standalone tool). odfInputPath and nuggetVariance are ignored.", fs::path(""),
      FileSystemPathParameter::ExtensionsType{".json"}, FileSystemPathParameter::PathType::InputFile));
```

Then add link calls (near the existing `linkParameters` block):

```cpp
  params.linkParameters(k_UseConfigFile_Key, k_ConfigFilePath_Key, true);
  params.linkParameters(k_UseConfigFile_Key, k_VolumeFractions_Key, false);
  params.linkParameters(k_UseConfigFile_Key, k_ThetaList_Key, false);
  params.linkParameters(k_UseConfigFile_Key, k_PhysicalSize_Key, false);
  params.linkParameters(k_UseConfigFile_Key, k_PhysicalSpacing_Key, false);
  params.linkParameters(k_UseConfigFile_Key, k_UseSeed_Key, false);
```
> Verify nested linking: `k_UseSeed_Key` already links `k_SeedValue_Key`. If a parameter may be the dependent of one linkable AND the linkable of another, confirm simplnx honors it (grep simplnx for a filter doing this; otherwise, link `k_SeedValue_Key` and `k_SeedArrayName_Key` to `k_UseConfigFile_Key=false` directly instead of relying on the `use_seed` chain, and accept that in config mode the seed group is simply hidden). Whichever works, the user-visible result must be: in config mode the seed group is hidden.

- [ ] **Step 3: Build plugin; run FilterValidationTest** — the new keys (`use_config_file` bool, `config_file_path` ends with `_path`) must satisfy `FilterValidationTest`.

Run: `cmake --build /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk95-Rel-EbsdLib --target simplnx_test && ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk95-Rel-EbsdLib -R "FilterValidation" --output-on-failure`
Expected: PASS (no MTRSim parameter-naming violations). If `simplnx_test`/`FilterValidation` target names differ, find them: `ctest --test-dir <build> -N | grep -i validation`.

- [ ] **Step 4: Commit**

```bash
git add src/MTRSim/Filters/MTRSimFilter.hpp src/MTRSim/Filters/MTRSimFilter.cpp
git commit -m "feat(filter): add config-file source toggle + linked show/hide"
```

---

### Task C3: config-mode preflight + execute

**Files:** Modify `src/MTRSim/Filters/MTRSimFilter.{hpp,cpp}`, `src/MTRSim/Filters/Algorithms/MTRSim.{hpp,cpp}`.

- [ ] **Step 1: Thread config values into the algorithm input** — in `src/MTRSim/Filters/Algorithms/MTRSim.hpp`, add to `MTRSimInputValues`:

```cpp
  bool useConfigFile = false;
  std::filesystem::path configFilePath;
```
(add `#include <filesystem>`).

- [ ] **Step 2: Centralize "resolve SimulationParams" in the algorithm** — in `MTRSim.cpp`, where `SimulationParams params;` is built, branch on config mode. Add `#include "LibMTRSim/ConfigIO.hpp"`:

```cpp
  mtrsim::SimulationParams params;
  if (m_InputValues->useConfigFile) {
    params = mtrsim::parseConfigJson(m_InputValues->configFilePath); // throws -> caught below
    // size/spacing/volumeFractions/thetaList/seed all come from the file
  } else {
    params.xLen = m_InputValues->physicalSize[0];
    params.yLen = m_InputValues->physicalSize[1];
    params.zLen = m_InputValues->physicalSize[2];
    params.dx = m_InputValues->physicalSpacing[0];
    params.dy = m_InputValues->physicalSpacing[1];
    params.dz = m_InputValues->physicalSpacing[2];
    params.volumeFractions = m_InputValues->volumeFractions[0];
    params.thetaList = m_InputValues->thetaList;
    params.seed = m_InputValues->seed;
  }
```
Wrap the parse in the existing try/catch (or add one) so a parse error returns `MakeErrorResult(-13520, ...)`. The seed in config mode comes from the JSON; the seed-array write in the FILTER's `executeImpl` should record `params.seed` actually used (see Step 4).

- [ ] **Step 3: Preflight reads the config to validate + size the geometry** — in `MTRSimFilter.cpp` `preflightImpl`, at the top, branch:

```cpp
  const bool useConfig = filterArgs.value<bool>(k_UseConfigFile_Key);
  std::vector<std::vector<double>> volumeFractions;
  std::vector<std::vector<double>> thetaList;
  std::vector<float32> size, spacing;
  if (useConfig) {
    mtrsim::SimulationParams cfg;
    try {
      cfg = mtrsim::parseConfigJson(filterArgs.value<FileSystemPathParameter::ValueType>(k_ConfigFilePath_Key));
    } catch (const std::exception& e) {
      return {MakeErrorResult<OutputActions>(-13520, fmt::format("MTRSim config file error: {}", e.what()))};
    }
    volumeFractions = {cfg.volumeFractions};
    thetaList = cfg.thetaList;
    size = {static_cast<float32>(cfg.xLen), static_cast<float32>(cfg.yLen), static_cast<float32>(cfg.zLen)};
    spacing = {static_cast<float32>(cfg.dx), static_cast<float32>(cfg.dy), static_cast<float32>(cfg.dz)};
  } else {
    volumeFractions = filterArgs.value<DynamicTableParameter::ValueType>(k_VolumeFractions_Key);
    thetaList = filterArgs.value<DynamicTableParameter::ValueType>(k_ThetaList_Key);
    size = filterArgs.value<std::vector<float32>>(k_PhysicalSize_Key);
    spacing = filterArgs.value<std::vector<float32>>(k_PhysicalSpacing_Key);
  }
```
Then the EXISTING validation (VF count == numComponents, sum ≈ 1, in [0,1], theta rows ≥ components−1 with 3 cols, spacing X/Y > 0) and grid-dim computation run against these local `volumeFractions/thetaList/size/spacing` variables instead of the direct `filterArgs` reads. Add a preflight info note when `useConfig`: `preflightUpdatedValues.push_back({"Parameter Source", "Config file: " + path.string()})`. Include `#include "LibMTRSim/ConfigIO.hpp"`.

- [ ] **Step 4: executeImpl passes config fields + records the right seed** — in `MTRSimFilter.cpp` `executeImpl`, set `inputValues.useConfigFile` and `inputValues.configFilePath`. For the seed-array record: in config mode the seed comes from the JSON, so read it for the array write:

```cpp
  inputValues.useConfigFile = filterArgs.value<bool>(k_UseConfigFile_Key);
  inputValues.configFilePath = filterArgs.value<FileSystemPathParameter::ValueType>(k_ConfigFilePath_Key);

  uint64 seed;
  if (inputValues.useConfigFile) {
    seed = mtrsim::parseConfigJson(inputValues.configFilePath).seed;  // file is validated in preflight
  } else {
    seed = filterArgs.value<uint64>(k_SeedValue_Key);
    if (!filterArgs.value<bool>(k_UseSeed_Key)) {
      seed = static_cast<uint64>(std::chrono::steady_clock::now().time_since_epoch().count());
    }
  }
  dataStructure.getDataRefAs<UInt64Array>(DataPath({filterArgs.value<std::string>(k_SeedArrayName_Key)}))[0] = seed;
  inputValues.seed = seed;
```
The algorithm builds `params` from the config in config mode (Step 2) but uses `m_InputValues->seed` for the RNG. **Make Step 2 use `m_InputValues->seed` for the rng regardless of mode** (the filter already resolved the correct seed): keep the line `std::mt19937_64 rng(m_InputValues->seed);` and, in config mode, after `params = parseConfigJson(...)`, set `params.seed = m_InputValues->seed;` for consistency (it isn't read by simulateMTR but keeps the struct truthful).

- [ ] **Step 5: Write filter tests** — add to `test/MTRSimTest.cpp`:
  - **Config-mode preflight VALID:** build the ODF DataStructure (3 components), set `k_UseConfigFile_Key=true`, `k_ConfigFilePath_Key` = an absolute path to `configs/default.json` (compose from a known repo path or a test-data macro — mirror how other MTRSim tests locate files; if none, write a temp 3-component config in the test). Assert preflight VALID and that the output `MTRIds`/`Eulers` arrays are created with the config-derived grid (1905×635 for default.json → check tuple count `1905*635`).
  - **Config-mode missing file → INVALID (−13520):** set `use_config_file=true`, `config_file_path` to a nonexistent path; assert preflight INVALID.
  - **Config-mode VF/component mismatch → INVALID:** point at a config whose `volumeFractions` length ≠ selected component count (write a temp 2-fraction config but select 3 components); assert INVALID.
  - Keep an existing manual-mode test to prove no regression.

```cpp
// sketch of the valid case
args.insertOrAssign(MTRSimFilter::k_UseConfigFile_Key, true);
args.insertOrAssign(MTRSimFilter::k_ConfigFilePath_Key, fs::path("/Users/mjackson/Workspace7/DREAM3D_Plugins/MTRSim/configs/default.json"));
auto pre = filter.preflight(ds, args);
SIMPLNX_RESULT_REQUIRE_VALID(pre.outputActions);
```
> Prefer writing a small temp config inside the test (with 3 fractions and a tiny size like 1.0×0.6) over depending on an absolute repo path, so the test is portable. Use a `std::filesystem::temp_directory_path()`-based temp file.

- [ ] **Step 6: Build, run MTRSim tests, verify pass; commit**

Run: `cmake --build /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk95-Rel-EbsdLib --target MTRSimUnitTest && ctest --test-dir /Users/mjackson/Workspace7/DREAM3D-Build/NX-Com-Qt69-Vtk95-Rel-EbsdLib -R "MTRSim::" --output-on-failure`
Expected: all PASS.

```bash
git add src/MTRSim/Filters/MTRSimFilter.hpp src/MTRSim/Filters/MTRSimFilter.cpp \
        src/MTRSim/Filters/Algorithms/MTRSim.hpp src/MTRSim/Filters/Algorithms/MTRSim.cpp test/MTRSimTest.cpp
git commit -m "feat(filter): config-file mode in preflight + execute, with tests"
```

---

## Phase D — Docs + integration

### Task D1: documentation

**Files:** Modify `docs/MTRSimFilter.md`.

- [ ] **Step 1: Document the two new behaviors** — add a "Configuration Source" subsection under Description explaining the `Load Simulation Parameters from Config File` toggle (which params it replaces, that `odfInputPath`/`nuggetVariance` are ignored, output names/ODF still come from the UI). Update the Performance section to note that progress is now reported continuously and the filter can be cancelled mid-run. Add error `-13520` (config file missing/invalid) to the error table. Commit.

```bash
git add docs/MTRSimFilter.md
git commit -m "docs: document config-file mode and progress/cancel in MTRSimFilter"
```

### Task D2: debug smoke + validation sweep

- [ ] **Step 1: Standalone dual check** — `cmake --build /Users/mjackson/Workspace7/Build/mtrsim-Rel && ctest --test-dir /Users/mjackson/Workspace7/Build/mtrsim-Rel --output-on-failure`. Expected: all library tests (including `[observer]`, `[config_io]`) PASS.
- [ ] **Step 2: Plugin tests + FilterValidationTest** — build `MTRSimUnitTest` and `simplnx_test`; run `ctest -R "MTRSim::"` and `ctest -R "FilterValidation"`. Expected: PASS.
- [ ] **Step 3: Debug smoke execution** — build the Debug `MTRSim` + `nxrunner`, then `<dbg>/Bin/nxrunner_d --execute pipelines/MTRSim_smoke_test.d3dpipeline`. Expected: clean run, no Eigen asserts, progress lines visible. Also add a tiny config-mode pipeline variant or temporarily flip `use_config_file` in a copy to exercise the config path once in Debug.
- [ ] **Step 4: clang-format** — the user runs a clang-format script before pushing; ensure new files are included. (Format C++ sources if a `clang-format` binary is available; otherwise note that the new files need formatting before push.)

---

## Self-Review Notes (for the executor)

- **Spec coverage:** interface (A1), simulateMTR/PGRF threading (A2), deep sampleN cancel (A3), default observers + CLI (A1/A4), `MTRSimResult::cancelled` (A2), filter observer (C1), config params + linking (C2), config-mode preflight/execute + parser (B1/B2/C3), tests (each task + D2), docs (D1). All spec sections map to a task.
- **APIs to verify against the live tree** (evolve over time): `IFilter::ProgressMessage` construction (C1), nested `linkParameters` behavior (C2), `FileSystemPathParameter::ValueType`/ctor (C2/C3), `SIMPLNX_RESULT_REQUIRE_VALID/INVALID` macros, the exact `simplnx_test`/`FilterValidation` ctest names (C2/D2). Grep the local simplnx at `/Users/mjackson/Workspace7/simplnx` and mirror a sibling.
- **Bit-stability guard:** the A3 modulo cancel check must not consume RNG on the non-cancel path; the pre-existing `[mtrsim_driver][statistical]` test is the regression guard (A3 Step 4).
- **Back-compat:** every new observer arg defaults to `nullptr`; existing callers/tests compile unchanged.
