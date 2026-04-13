// generate_pole_figures.cpp
// Reads a sim_results.csv and generates HCP pole figure images for each
// MTR component using EbsdLib's PoleFigureCompositor.
//
// Usage:
//   generate_pole_figures <csv_path> <output_dir> [--label <prefix>] [--dim <pixels>]
//
// Outputs:
//   <output_dir>/<prefix>_component_N_pf.png   (per component)
//   <output_dir>/<prefix>_all_pf.png           (all orientations combined)

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include <EbsdLib/Core/EbsdDataArray.hpp>
#include <EbsdLib/LaueOps/HexagonalOps.h>
#include <EbsdLib/Utilities/PoleFigureCompositor.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
struct VoxelRow
{
  double phi1;
  double PHI;
  double phi2;
  int mtrIndex;
};

// ---------------------------------------------------------------------------
std::vector<VoxelRow> readCSV(const std::string& path)
{
  std::vector<VoxelRow> rows;
  std::ifstream in(path);
  if(!in.is_open())
  {
    std::cerr << "Error: cannot open " << path << "\n";
    std::exit(1);
  }

  std::string line;
  std::getline(in, line); // skip header

  while(std::getline(in, line))
  {
    // Schema: x,y,z,phi1,PHI,phi2,mtr_index
    std::istringstream ss(line);
    std::string tok;
    VoxelRow r;

    std::getline(ss, tok, ','); // x
    std::getline(ss, tok, ','); // y
    std::getline(ss, tok, ','); // z
    std::getline(ss, tok, ',');
    r.phi1 = std::stod(tok);
    std::getline(ss, tok, ',');
    r.PHI = std::stod(tok);
    std::getline(ss, tok, ',');
    r.phi2 = std::stod(tok);
    std::getline(ss, tok, ',');
    r.mtrIndex = std::stoi(tok);

    rows.push_back(r);
  }
  return rows;
}

// ---------------------------------------------------------------------------
ebsdlib::FloatArrayType::Pointer toEulerArray(const std::vector<VoxelRow>& rows)
{
  std::vector<size_t> cDims = {3};
  auto eulers = ebsdlib::FloatArrayType::CreateArray(rows.size(), cDims, "Eulers", true);
  for(size_t i = 0; i < rows.size(); i++)
  {
    float* ptr = eulers->getTuplePointer(i);
    ptr[0] = static_cast<float>(rows[i].phi1);
    ptr[1] = static_cast<float>(rows[i].PHI);
    ptr[2] = static_cast<float>(rows[i].phi2);
  }
  return eulers;
}

// ---------------------------------------------------------------------------
void writePNG(const ebsdlib::CompositePoleFigureResult& result, const std::string& path)
{
  // Convert RGBA to RGB for stb_image_write
  int w = result.width;
  int h = result.height;
  std::vector<uint8_t> rgb(w * h * 3);
  const uint8_t* rgba = result.image->getPointer(0);
  for(int i = 0; i < w * h; i++)
  {
    rgb[i * 3 + 0] = rgba[i * 4 + 0];
    rgb[i * 3 + 1] = rgba[i * 4 + 1];
    rgb[i * 3 + 2] = rgba[i * 4 + 2];
  }
  if(stbi_write_png(path.c_str(), w, h, 3, rgb.data(), w * 3) == 0)
  {
    std::cerr << "Error: failed to write " << path << "\n";
    std::exit(1);
  }
  std::cout << "  wrote " << path << " (" << w << "x" << h << ")\n";
}

// ---------------------------------------------------------------------------
void generatePoleFigure(const std::vector<VoxelRow>& rows, const std::string& title, const std::string& outPath, int imageDim)
{
  auto eulers = toEulerArray(rows);

  ebsdlib::CompositePoleFigureConfiguration_t config;
  config.eulers = eulers.get();
  config.imageDim = imageDim;
  config.lambertDim = imageDim / 2;
  config.numColors = 32;
  config.discrete = false;
  config.discreteHeatMap = false;
  config.flipFinalImage = true;
  config.laueOpsIndex = 0; // HexagonalOps (6/mmm)
  config.layoutType = ebsdlib::PoleFigureLayoutType::Square;
  config.phaseName = "Ti-6Al-4V (alpha)";
  config.phaseNumber = 1;
  config.title = title;

  ebsdlib::PoleFigureCompositor compositor;
  ebsdlib::CompositePoleFigureResult result = compositor.generateCompositeImage(config);

  writePNG(result, outPath);
}

// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  if(argc < 3)
  {
    std::cerr << "Usage: generate_pole_figures <csv_path> <output_dir> [--label <prefix>] [--dim <pixels>]\n";
    return 1;
  }

  std::string csvPath = argv[1];
  std::string outDir = argv[2];
  std::string label = "cpp";
  int imageDim = 1024;

  for(int i = 3; i < argc; i++)
  {
    if(std::string(argv[i]) == "--label" && i + 1 < argc)
    {
      label = argv[++i];
    }
    else if(std::string(argv[i]) == "--dim" && i + 1 < argc)
    {
      imageDim = std::stoi(argv[++i]);
    }
  }

  fs::create_directories(outDir);

  std::cout << "Reading " << csvPath << "...\n";
  auto allRows = readCSV(csvPath);
  std::cout << "  " << allRows.size() << " voxels\n";

  // Find unique components
  int maxComp = 0;
  for(const auto& r : allRows)
  {
    if(r.mtrIndex > maxComp)
    {
      maxComp = r.mtrIndex;
    }
  }

  // Split by component
  std::vector<std::vector<VoxelRow>> byComponent(maxComp);
  for(const auto& r : allRows)
  {
    byComponent[r.mtrIndex - 1].push_back(r);
  }

  // Generate per-component pole figures
  for(int c = 0; c < maxComp; c++)
  {
    std::string title = "MTR Component " + std::to_string(c + 1) + " - ODF " + std::to_string(c + 1) + " (n=" + std::to_string(byComponent[c].size()) + ")";
    std::string path = outDir + "/" + label + "_component_" + std::to_string(c + 1) + "_pf.png";
    std::cout << "Generating pole figure: " << title << "\n";
    generatePoleFigure(byComponent[c], title, path, imageDim);
  }

  // Generate combined pole figure
  {
    std::string title = "All MTR Components Combined (n=" + std::to_string(allRows.size()) + ")";
    std::string path = outDir + "/" + label + "_all_pf.png";
    std::cout << "Generating pole figure: " << title << "\n";
    generatePoleFigure(allRows, title, path, imageDim);
  }

  std::cout << "Done.\n";
  return 0;
}
