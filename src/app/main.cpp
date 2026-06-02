// MTRsim simulation driver — CLI entry point.
// Orchestrates PGRF simulation, ODF loading, orientation sampling, and output.

#include "IPFMapper.hpp"
#include "MTRSimDriver.hpp"
#include "ODFSampler.hpp"
#include "SimulationParams.hpp"

#include <CLI/CLI.hpp>
#include <Eigen/Dense>
#include <hdf5.h>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

#include <cmath>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

// ─────────────────────────────────────────────────────────────────────────────
// HDF5 helper utilities

int64_t readH5Int64(hid_t fileId, const std::string& path)
{
  const hid_t dsId = H5Dopen2(fileId, path.c_str(), H5P_DEFAULT);
  if(dsId < 0)
  {
    throw std::runtime_error("HDF5: cannot open dataset: " + path);
  }
  int64_t val = 0;
  H5Dread(dsId, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
  H5Dclose(dsId);
  return val;
}

std::vector<double> readH5Vector(hid_t fileId, const std::string& path)
{
  const hid_t dsId = H5Dopen2(fileId, path.c_str(), H5P_DEFAULT);
  if(dsId < 0)
  {
    throw std::runtime_error("HDF5: cannot open dataset: " + path);
  }

  const hid_t spaceId = H5Dget_space(dsId);
  hsize_t dims[1] = {0};
  H5Sget_simple_extent_dims(spaceId, dims, nullptr);
  H5Sclose(spaceId);

  std::vector<double> buf(dims[0]);
  H5Dread(dsId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
  H5Dclose(dsId);
  return buf;
}

// ─────────────────────────────────────────────────────────────────────────────
// Load all ODF components from an HDF5 file produced by
// convert_simulation_ODF.py.
//
// Expected HDF5 layout:
//   /ODF_best/num_components               — int64 scalar
//   /ODF_best/component_N/ODFval           — (186624,) float64
//   /ODF_best/component_N/phi1_bins        — (73,)  float64  [rad]
//   /ODF_best/component_N/PHI_bins         — (37,)  float64  [rad]
//   /ODF_best/component_N/phi2_bins        — (73,)  float64  [rad]
std::vector<mtrsim::ODFComponent> loadODFComponents(const std::string& hdfPath)
{
  const hid_t fileId = H5Fopen(hdfPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if(fileId < 0)
  {
    throw std::runtime_error("Cannot open HDF5 file: " + hdfPath);
  }

  const int64_t numComponents = readH5Int64(fileId, "/ODF_best/num_components");
  spdlog::info("ODF file: {} component(s)", numComponents);

  std::vector<mtrsim::ODFComponent> components;
  components.reserve(static_cast<std::size_t>(numComponents));

  for(int64_t j = 0; j < numComponents; ++j)
  {
    const std::string prefix = "/ODF_best/component_" + std::to_string(j);

    auto odfVec = readH5Vector(fileId, prefix + "/ODFval");
    auto phi1Vec = readH5Vector(fileId, prefix + "/phi1_bins");
    auto phiVec = readH5Vector(fileId, prefix + "/PHI_bins");
    auto phi2Vec = readH5Vector(fileId, prefix + "/phi2_bins");

    mtrsim::ODFComponent comp;
    comp.odfVal = Eigen::Map<Eigen::VectorXd>(odfVec.data(), static_cast<Eigen::Index>(odfVec.size()));
    comp.phi1Bins = Eigen::Map<Eigen::VectorXd>(phi1Vec.data(), static_cast<Eigen::Index>(phi1Vec.size()));
    comp.phiBins = Eigen::Map<Eigen::VectorXd>(phiVec.data(), static_cast<Eigen::Index>(phiVec.size()));
    comp.phi2Bins = Eigen::Map<Eigen::VectorXd>(phi2Vec.data(), static_cast<Eigen::Index>(phi2Vec.size()));

    // Normalise so that odfVal sums to 1 (matches MATLAB pre-processing in
    // simulate_MTRs.m)
    const double total = comp.odfVal.sum();
    if(total > 0.0)
    {
      comp.odfVal /= total;
    }

    components.push_back(std::move(comp));
  }

  H5Fclose(fileId);
  return components;
}

} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────
// main

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

  // ── Parse JSON config
  // ────────────────────────────────────────────────────────
  mtrsim::SimulationParams params;
  params.outputDir = outputDir;
  params.seed = seed; // may be overridden by JSON (unless CLI --seed was set)

  if(!configPath.empty())
  {
    std::ifstream f(configPath);
    if(!f.is_open())
    {
      spdlog::error("Cannot open config file: {}", configPath);
      return 1;
    }

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
      if(j.contains("odfInputPath"))
        params.odfInputPath = j["odfInputPath"].get<std::string>();
      // JSON seed only applies when CLI --seed was not explicitly provided
      // (seed == 0)
      if(j.contains("seed") && params.seed == 0)
      {
        params.seed = j["seed"].get<uint64_t>();
      }
    } catch(const nlohmann::json::exception& e)
    {
      spdlog::error("JSON parse error: {}", e.what());
      return 1;
    }
  }

  // ── Setup RNG
  // ────────────────────────────────────────────────────────────────
  std::mt19937_64 rng;
  if(params.seed == 0)
  {
    std::random_device rd;
    rng.seed(rd());
    spdlog::info("Using random seed from std::random_device");
  }
  else
  {
    rng.seed(params.seed);
    spdlog::info("Using fixed seed: {}", params.seed);
  }

  // ── Compute grid dimensions
  // ──────────────────────────────────────────────────
  const int nx = static_cast<int>(std::round(params.xLen / params.dx));
  const int ny = static_cast<int>(std::round(params.yLen / params.dy));
  const int nz = std::max(static_cast<int>(std::round(params.zLen / params.dz)), 1);
  const int N = nx * ny * nz;

  spdlog::info("Grid: nx={} ny={} nz={} N={}", nx, ny, nz, N);
  spdlog::info("  xLen={:.3f} mm  yLen={:.3f} mm  zLen={:.3f} mm", params.xLen, params.yLen, params.zLen);
  spdlog::info("  dx={}  dy={}  dz={} [mm]", params.dx, params.dy, params.dz);

  // ── Build spatial coordinate matrix
  // ──────────────────────────────────────────
  // Row order matches simulateMTR's SIMPLNX z,y,x output: z slowest, x fastest.
  //   k = iz*(ny*nx) + iy*nx + ix
  //   spatialCoords(k) = [(ix+1)*dx, (iy+1)*dy, (iz+1)*dz]
  // 1-based coordinate values are preserved (matching the original convention).
  Eigen::MatrixXd spatialCoords(N, 3);
  {
    for(int iz = 0; iz < nz; ++iz)
    {
      for(int iy = 0; iy < ny; ++iy)
      {
        for(int ix = 0; ix < nx; ++ix)
        {
          const int k = iz * (ny * nx) + iy * nx + ix;
          spatialCoords(k, 0) = (ix + 1) * params.dx;
          spatialCoords(k, 1) = (iy + 1) * params.dy;
          spatialCoords(k, 2) = (iz + 1) * params.dz;
        }
      }
    }
  }

  // ── Load ODF components from HDF5
  // ───────────────────────────────────────────
  spdlog::info("Loading ODF from: {}", params.odfInputPath);
  std::vector<mtrsim::ODFComponent> odfComponents;
  try
  {
    odfComponents = loadODFComponents(params.odfInputPath);
  } catch(const std::exception& e)
  {
    spdlog::error("ODF load failed: {}", e.what());
    return 1;
  }

  if(odfComponents.empty())
  {
    spdlog::error("No ODF components loaded.");
    return 1;
  }

  // ── Run full MTR simulation (PGRF + ODF sampling + remap to z,y,x order)
  // ─────────────────────────────────────────────────────────────────────────
  // The MATLAB ODF HDF5 layout is a fixed 5-degree Bunge-Euler grid:
  // 72 (phi1) x 36 (PHI) x 72 (phi2) = 186624 bins.
  constexpr int k_OdfBinsPhi1 = 72;
  constexpr int k_OdfBinsPHI = 36;
  constexpr int k_OdfBinsPhi2 = 72;

  spdlog::info("Running MTR simulation...");
  mtrsim::MTRSimResult sim = mtrsim::simulateMTR(params, odfComponents, rng, k_OdfBinsPhi1, k_OdfBinsPHI, k_OdfBinsPhi2);
  spdlog::info("MTR simulation complete.");

  Eigen::VectorXd phi1Vec = Eigen::Map<Eigen::VectorXd>(sim.phi1.data(), static_cast<Eigen::Index>(sim.phi1.size()));
  Eigen::VectorXd phiVec = Eigen::Map<Eigen::VectorXd>(sim.phi.data(), static_cast<Eigen::Index>(sim.phi.size()));
  Eigen::VectorXd phi2Vec = Eigen::Map<Eigen::VectorXd>(sim.phi2.data(), static_cast<Eigen::Index>(sim.phi2.size()));

  // ── Write IPF map PNG
  // ────────────────────────────────────────────────────────
  const std::string ipfPath = params.outputDir + "/sim_IPF_map.png";
  spdlog::info("Writing IPF map: {}", ipfPath);
  try
  {
    mtrsim::IPFMapper mapper{mtrsim::CrystalSystem::HCP};
    const Eigen::MatrixXd coords2D = spatialCoords.leftCols(2);
    mapper.writePNG(coords2D, phi1Vec, phiVec, phi2Vec, ipfPath);
    spdlog::info("IPF map written.");
  } catch(const std::exception& e)
  {
    spdlog::warn("IPF map write failed: {}", e.what());
  }

  // ── Write results CSV
  // ────────────────────────────────────────────────────────
  const std::string csvPath = params.outputDir + "/sim_results.csv";
  spdlog::info("Writing results CSV: {}", csvPath);
  {
    std::ofstream csv(csvPath);
    if(!csv.is_open())
    {
      spdlog::error("Cannot open output file: {}", csvPath);
      return 1;
    }

    csv << "x,y,z,phi1,PHI,phi2,mtr_index\n";
    csv << std::fixed;
    csv.precision(6);

    // Use sim dimensions to tie loop bounds to the actual result.
    const int simN = sim.nx * sim.ny * sim.nz;
    for(int i = 0; i < simN; ++i)
    {
      csv << spatialCoords(i, 0) << ',' << spatialCoords(i, 1) << ',' << spatialCoords(i, 2) << ',' << phi1Vec[i] << ',' << phiVec[i] << ',' << phi2Vec[i] << ','
          << sim.mtrIndex[static_cast<std::size_t>(i)] << '\n';
    }
  }
  spdlog::info("Results CSV written.");

  spdlog::info("Simulation complete.");
  return 0;
}
