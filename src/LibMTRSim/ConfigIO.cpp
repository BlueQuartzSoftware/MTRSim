#include "ConfigIO.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <stdexcept>

namespace mtrsim
{

SimulationParams ParseConfigJson(const std::filesystem::path& path)
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
    {
      params.x_len = j["xLen"].get<double>();
    }
    if(j.contains("yLen"))
    {
      params.y_len = j["yLen"].get<double>();
    }
    if(j.contains("zLen"))
    {
      params.z_len = j["zLen"].get<double>();
    }
    if(j.contains("dx"))
    {
      params.dx = j["dx"].get<double>();
    }
    if(j.contains("dy"))
    {
      params.dy = j["dy"].get<double>();
    }
    if(j.contains("dz"))
    {
      params.dz = j["dz"].get<double>();
    }
    if(j.contains("volumeFractions"))
    {
      params.volume_fractions = j["volumeFractions"].get<std::vector<double>>();
    }
    if(j.contains("thetaList"))
    {
      params.theta_list = j["thetaList"].get<std::vector<std::vector<double>>>();
    }
    if(j.contains("nuggetVariance"))
    {
      params.nugget_variance = j["nuggetVariance"].get<std::vector<double>>();
    }
    if(j.contains("odfInputPath"))
    {
      params.odf_input_path = j["odfInputPath"].get<std::string>();
    }
    if(j.contains("seed"))
    {
      params.seed = j["seed"].get<uint64_t>();
    }
  } catch(const nlohmann::json::exception& e)
  {
    throw std::runtime_error(std::string("parseConfigJson: invalid JSON: ") + e.what());
  }
  return params;
}

} // namespace mtrsim
