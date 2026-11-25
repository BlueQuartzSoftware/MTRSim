# MTRSim: Simulation Observer + Config-File Input — Design

**Date:** 2026-06-02
**Status:** Design approved — ready for implementation plan
**Scope:** Two cohesive enhancements to the MTRSim library + `MTRSimFilter`:
1. A progress/cancellation observer threaded through `simulateMTR`.
2. An optional MTRSim config-JSON input on the filter (MATLAB-migration convenience).

---

## 1. Background

`mtrsim::simulateMTR(params, odfComponents, rng, n1, nPHI, n2)` (in
`src/LibMTRSim/MTRSimDriver.{hpp,cpp}`) runs the full pipeline: PGRF assignment →
per-component ODF sampling → per-voxel orientation assignment → remap to SIMPLNX
z,y,x order. It is currently a single opaque call with no progress reporting and
no way to cancel; for the default ~1.2M-voxel grid it can run for a long time.

It has two callers:
- `src/MTRSim/Filters/Algorithms/MTRSim.cpp` (the DREAM3D-NX filter algorithm),
  which holds a `MessageHandler` and an `std::atomic_bool& m_ShouldCancel`.
- `src/app/main.cpp` (the standalone CLI), which currently relies on
  `PGRFSimulation`'s internal `spdlog` logging.

The long pole is the per-component ODF-sampling loop inside `simulateMTR`
(`ODFSampler::sampleN` is called once per component, drawing N orientations
each) and the per-voxel assignment loop.

---

## 2. Part A — Progress & Cancellation Observer

### 2.1 The interface

New header `src/LibMTRSim/ISimulationObserver.hpp`:

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
 * cancellation. The library calls these from the thread running simulateMTR;
 * implementations must be cheap and must not throw.
 */
class LIBMTRSIM_EXPORT ISimulationObserver {
public:
  ISimulationObserver() = default;
  virtual ~ISimulationObserver() = default;

  // Polymorphic, non-copyable base (prevents slicing).
  ISimulationObserver(const ISimulationObserver&) = delete;
  ISimulationObserver& operator=(const ISimulationObserver&) = delete;
  ISimulationObserver(ISimulationObserver&&) = delete;
  ISimulationObserver& operator=(ISimulationObserver&&) = delete;

  /// Report progress. `done`/`total` describe the current phase; `message`
  /// names the phase (e.g. "Sampling orientations (component 2/3)").
  virtual void updateProgress(int64_t done, int64_t total, const std::string& message) = 0;

  /// Polled at checkpoints. Returning true stops the simulation early.
  [[nodiscard]] virtual bool shouldCancel() const = 0;
};

} // namespace mtrsim
```

C++ best-practice notes:
- Pure-virtual methods (`= 0`), defaulted virtual destructor, `[[nodiscard]]`
  on the predicate, `const` correctness on `shouldCancel`.
- Copy/move deleted on the polymorphic base (rule-of-five, anti-slicing).
- Header-only interface; no `.cpp`.

### 2.2 Default implementations

In `src/LibMTRSim/SimulationObservers.hpp` (header-only):
- `NullObserver` — both methods no-op / `return false`. Used as the implicit
  default so existing callers are unaffected.
- `ConsoleObserver` — `updateProgress` logs via `spdlog` (throttled to avoid
  log spam); `shouldCancel` returns `false`. Used by `main.cpp` to preserve its
  console feedback.

### 2.3 Signature changes

```cpp
// MTRSimDriver.hpp
LIBMTRSIM_EXPORT MTRSimResult simulateMTR(
    const SimulationParams& params,
    const std::vector<ODFComponent>& odfComponents,
    std::mt19937_64& rng,
    int n1, int nPHI, int n2,
    ISimulationObserver* observer = nullptr);   // NEW, nullable

// ODFSampler.hpp
Eigen::MatrixXd sampleN(int n, const ODFComponent& component,
                        const ODFComponent& uniform,
                        ISimulationObserver* observer = nullptr,  // NEW
                        int64_t progressBase = 0, int64_t progressTotal = 0); // for global progress
```

- `observer == nullptr` ⇒ run exactly as today (no progress, no cancel checks):
  full backward compatibility.
- `simulateMTR` reports at PGRF stage boundaries (threshold selection, each
  latent Gaussian field, assignment) and forwards the observer into each
  `sampleN` call and the per-voxel assignment loop.
- **Throttling:** a small internal helper updates progress at most every ~1% of
  `total` (or every K iterations), so the observer/UI isn't flooded. Cancel is
  polled at the same checkpoints — frequently enough to feel responsive but not
  per-iteration.

### 2.4 Cancellation semantics

- `MTRSimResult` gains `bool cancelled = false;`.
- On `observer->shouldCancel()` at any checkpoint, `simulateMTR` returns early
  with `cancelled = true` and whatever arrays are allocated (callers must not
  consume them). No exceptions are used for control flow.
- The filter algorithm checks `m_ShouldCancel` (its observer's source) after the
  call — as it already does — and returns `{}` (success, no output) without
  writing arrays, matching the simplnx cancel convention.

### 2.5 Filter observer

`MTRSim.cpp` defines a private `FilterObserver : mtrsim::ISimulationObserver`
that captures `const IFilter::MessageHandler&` and `const std::atomic_bool&`:
- `updateProgress` → emits an `IFilter::ProgressMessage` (percent =
  `done*100/total`) with the phase text (see `progress-messaging` conventions).
- `shouldCancel` → returns the captured atomic's value.

`main.cpp` passes a `ConsoleObserver` (or `nullptr`).

---

## 3. Part B — Optional MTRSim Config-File Input (filter only)

### 3.1 New parameters

On `MTRSimFilter`:
- `k_UseConfigFile_Key` = `"use_config_file"` — linkable `BoolParameter`
  "Load Simulation Parameters from Config File" (default `false`).
- `k_ConfigFilePath_Key` = `"config_file_path"` — `FileSystemPathParameter`
  (InputFile, extension `.json`).

### 3.2 Linked show/hide UX

Using `params.linkParameters(...)`:
- Shown when `use_config_file == true`: `config_file_path`.
- Shown when `use_config_file == false`: `volume_fractions`, `theta_list`,
  `physical_size`, `physical_spacing`, and the seed group (`use_seed`,
  `seed_value`, `seed_array_name`).
- Always shown: `input_odf_geometry_path`, `odf_component_arrays`, and all
  output names/paths.

> Implementation note: `use_seed` is itself a linkable that controls
> `seed_value`. Verify simplnx supports a parameter being both a dependent (of
> `use_config_file`) and a linkable (for `seed_value`). If nested linking is not
> supported, fall back to leaving the seed group visible in config mode but
> documenting that the config's `seed` takes precedence there.

### 3.3 Shared config parser

Extract `main.cpp`'s JSON-parsing block into a reusable LibMTRSim helper so the
CLI and the filter share one implementation:

```cpp
// src/LibMTRSim/ConfigIO.hpp
namespace mtrsim {
/// Parse an MTRSim config JSON into a SimulationParams. Throws
/// std::runtime_error on missing file / invalid JSON. Unknown keys
/// (odfInputPath, nuggetVariance, comments) are ignored. Fields not present
/// in the JSON keep their SimulationParams defaults.
LIBMTRSIM_EXPORT SimulationParams parseConfigJson(const std::filesystem::path& path);
}
```

`main.cpp` is refactored to call `parseConfigJson` (removing its inline parsing).

### 3.4 Filter behavior in config mode

When `use_config_file == true`, both `preflightImpl` and `executeImpl` obtain
the simulation parameters from `parseConfigJson(config_file_path)` instead of
the UI fields. Specifically the JSON supplies: `xLen/yLen/zLen`, `dx/dy/dz`,
`volumeFractions`, `thetaList`, `seed` (treated as a fixed seed; equivalent to
`use_seed == true`).

- **Preflight** parses the file; a missing file or parse error is a preflight
  error (`-13520`). The JSON-derived values are validated with the *same* rules
  as manual mode: `volumeFractions` count must equal the number of selected ODF
  component arrays; sum ≈ 1.0; each in [0,1]; theta rows ≥ components−1 with 3
  columns; spacing X/Y > 0. The output geometry dims are computed from the
  JSON's size/spacing.
- **Ignored JSON keys:** `odfInputPath` (the ODF comes from the selected
  DataStructure geometry, not a file) and `nuggetVariance` (unused by the
  simulation). A preflight info note lists the effective grid + that parameters
  came from the config file.
- ODF selection and all output names/paths always come from the UI.

When `use_config_file == false`, behavior is exactly as today.

---

## 4. Testing

### 4.1 Observer (Part A) — LibMTRSim Catch2
- A `RecordingObserver` test double that records each `updateProgress` call and
  can be configured to return `true` from `shouldCancel` after the Kth update.
- `simulateMTR` with the recorder: assert progress is reported, `done ≤ total`
  always, and the final update reaches `total` on a completed run.
- Cancel test: recorder cancels mid-run → `simulateMTR` returns
  `cancelled == true` and does not run to completion (bounded number of updates).
- `sampleN` directly with an observer that cancels → returns promptly.
- `NullObserver`/`nullptr` path produces identical results to the pre-change
  behavior (regression guard via a fixed seed).

### 4.2 Config parser (Part B) — LibMTRSim Catch2
- `parseConfigJson` on a known config (`configs/default.json`) → expected
  `SimulationParams` fields; confirms `odfInputPath`/`nuggetVariance` ignored.
- Missing file and malformed JSON → throws.

### 4.3 Filter (Part B) — simplnx FilterValidationTest + unit tests
- `FilterValidationTest` must pass (new keys end with `_path` where required;
  `use_config_file`/`config_file_path` follow conventions).
- Config-mode preflight: valid config → VALID with the info note; component-count
  mismatch vs ODF arrays → invalid; missing/garbled config → invalid (−13520).
- Manual-mode tests unchanged.
- A cancelled execute writes no output arrays.

### 4.4 Integration
- The Debug smoke pipeline (`pipelines/MTRSim_smoke_test.d3dpipeline`) executes
  clean in the Debug build after the changes (Eigen assertions on).
- Run `FilterValidationTest` (now part of `simplnx_test`) after the filter
  parameter changes.

---

## 5. Files

**Create:**
- `src/LibMTRSim/ISimulationObserver.hpp`
- `src/LibMTRSim/SimulationObservers.hpp` (`NullObserver`, `ConsoleObserver`)
- `src/LibMTRSim/ConfigIO.{hpp,cpp}` (`parseConfigJson`)
- Tests: extend `tests/test_mtrsim_driver.cpp`; add `tests/test_config_io.cpp`.

**Modify:**
- `src/LibMTRSim/MTRSimDriver.{hpp,cpp}` — `simulateMTR` observer arg + cancel;
  `MTRSimResult::cancelled`.
- `src/LibMTRSim/ODFSampler.{hpp,cpp}` — `sampleN` observer arg + cancel checks.
- `src/app/main.cpp` — use `ConsoleObserver`; use `parseConfigJson`.
- `src/MTRSim/Filters/MTRSimFilter.{hpp,cpp}` — config-file params + linking +
  config-mode preflight/execute.
- `src/MTRSim/Filters/Algorithms/MTRSim.{hpp,cpp}` — `FilterObserver`; pass it
  into `simulateMTR`; build `SimulationParams` from config when in config mode.
- `src/LibMTRSim/CMakeLists.txt`, `MTRSimPlugin.cmake`, `tests/CMakeLists.txt` —
  register new sources/tests.
- `docs/MTRSimFilter.md` — document config-file mode + progress/cancel.

---

## 6. Out of Scope
- Parallelizing the simulation (the observer enables cancel/progress, not speed).
- Changing the ODF on-disk format or the MATLAB code.
- Writing the config's `odfInputPath` ODF from within the filter (ODF still comes
  from upstream Read/Compute ODF filters).
