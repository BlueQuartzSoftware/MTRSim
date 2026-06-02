#include "ConfigIO.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <stdexcept>

namespace mtrsim {

SimulationParams parseConfigJson(const std::filesystem::path& path)
{
  std::ifstream f(path);
  if(!f.is_open())
  {
    throw std::runtime_error("parseConfigJson: cannot open config file: " + path.string());
  }
  SimulationParams params;
  try
  {
    const nlohmann::json j = nlohmann::json::parse(f);
    if(j.contains("xLen"))
      params.xLen = j["xLen"].get<double>();
    if(j.contains("yLen"))
      params.yLen = j["yLen"].get<double>();
    if(j.contains("zLen"))
      params.zLen = j["zLen"].get<double>();
    if(j.contains("dx"))
      params.dx = j["dx"].get<double>();
    if(j.contains("dy"))
      params.dy = j["dy"].get<double>();
    if(j.contains("dz"))
      params.dz = j["dz"].get<double>();
    if(j.contains("volumeFractions"))
      params.volumeFractions = j["volumeFractions"].get<std::vector<double>>();
    if(j.contains("thetaList"))
      params.thetaList = j["thetaList"].get<std::vector<std::vector<double>>>();
    if(j.contains("nuggetVariance"))
      params.nuggetVariance = j["nuggetVariance"].get<std::vector<double>>();
    if(j.contains("seed"))
      params.seed = j["seed"].get<uint64_t>();
  } catch(const nlohmann::json::exception& e)
  {
    throw std::runtime_error(std::string("parseConfigJson: invalid JSON: ") + e.what());
  }
  return params;
}

} // namespace mtrsim
