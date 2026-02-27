#include "SimulationParams.hpp"

#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

#include <fstream>
#include <string>

int main(int argc, char** argv)
{
  CLI::App app{"MTRsim — Microtexture Region Simulator"};

  std::string configPath;
  app.add_option("-c,--config", configPath, "Path to JSON configuration file");

  std::string outputDir = ".";
  app.add_option("-o,--output", outputDir, "Output directory");

  uint64_t seed = 0;
  app.add_option("--seed", seed, "Random seed (0 = use random_device)");

  CLI11_PARSE(app, argc, argv);

  spdlog::info("MTRsim v{}.{}.{}", 1, 0, 0);

  mtrsim::SimulationParams params;
  params.outputDir = outputDir;
  params.seed = seed;

  if(!configPath.empty())
  {
    std::ifstream f(configPath);
    if(!f.is_open())
    {
      spdlog::error("Cannot open config file: {}", configPath);
      return 1;
    }
    // TODO: parse JSON into params
  }

  // TODO: run simulation
  spdlog::info("Simulation complete.");

  return 0;
}
