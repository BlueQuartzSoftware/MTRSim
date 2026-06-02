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
