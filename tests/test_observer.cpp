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
