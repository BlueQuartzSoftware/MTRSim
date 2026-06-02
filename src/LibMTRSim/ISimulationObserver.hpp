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
