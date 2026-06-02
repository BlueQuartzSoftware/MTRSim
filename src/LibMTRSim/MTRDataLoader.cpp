// Port of load_MTR_data.m — load experimental EBSD data from CSV files.
// Daniel M. Sparkman, 04/10/2018.

#include "MTRDataLoader.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mtrsim
{

namespace
{

// Read a single-column CSV (one numeric value per non-empty line).
std::vector<double> readColumn(const std::string& path)
{
  std::ifstream f(path);
  if(!f.is_open())
  {
    throw std::runtime_error("MTRDataLoader: cannot open file: " + path);
  }

  std::vector<double> vals;
  std::string line;
  while(std::getline(f, line))
  {
    // Strip trailing commas and whitespace (MATLAB csvwrite sometimes adds
    // them)
    while(!line.empty() && (line.back() == ',' || line.back() == '\r' || std::isspace(static_cast<unsigned char>(line.back()))))
    {
      line.pop_back();
    }
    if(line.empty())
    {
      continue;
    }
    vals.push_back(std::stod(line));
  }
  return vals;
}

// Read a multi-column CSV (comma-separated, N rows × ncols columns).
// Returns a flat row-major vector: [row0_col0, row0_col1, ..., row1_col0, ...].
std::vector<double> readMultiColumn(const std::string& path, int& outNRows, int& outNCols)
{
  std::ifstream f(path);
  if(!f.is_open())
  {
    throw std::runtime_error("MTRDataLoader: cannot open file: " + path);
  }

  std::vector<std::vector<double>> rows;
  std::string line;
  while(std::getline(f, line))
  {
    // Strip trailing whitespace/CR
    while(!line.empty() && (line.back() == '\r' || std::isspace(static_cast<unsigned char>(line.back()))))
    {
      line.pop_back();
    }
    if(line.empty())
    {
      continue;
    }

    std::vector<double> row;
    std::stringstream ss(line);
    std::string token;
    while(std::getline(ss, token, ','))
    {
      // Trim token whitespace
      while(!token.empty() && std::isspace(static_cast<unsigned char>(token.front())))
      {
        token.erase(token.begin());
      }
      while(!token.empty() && std::isspace(static_cast<unsigned char>(token.back())))
      {
        token.pop_back();
      }
      if(!token.empty())
      {
        row.push_back(std::stod(token));
      }
    }
    if(!row.empty())
    {
      rows.push_back(std::move(row));
    }
  }

  outNRows = static_cast<int>(rows.size());
  outNCols = rows.empty() ? 0 : static_cast<int>(rows[0].size());

  std::vector<double> flat;
  flat.reserve(static_cast<std::size_t>(outNRows * outNCols));
  for(const auto& row : rows)
  {
    for(int c = 0; c < outNCols; ++c)
    {
      flat.push_back(c < static_cast<int>(row.size()) ? row[c] : 0.0);
    }
  }
  return flat;
}

} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────
// load — reads five CSV files from directoryPath:
//   X_Position.csv, Y_Position.csv, EulerAngles.csv, ParentIds.csv, BoolMTR.csv

EBSDData MTRDataLoader::load(const std::string& directoryPath)
{
  const std::string sep = directoryPath + "/";

  // ── Spatial coordinates
  // ──────────────────────────────────────────────────────
  const std::vector<double> xVals = readColumn(sep + "X_Position.csv");
  const std::vector<double> yVals = readColumn(sep + "Y_Position.csv");

  const int N = static_cast<int>(xVals.size());
  if(static_cast<int>(yVals.size()) != N)
  {
    throw std::runtime_error("MTRDataLoader: X_Position and Y_Position row counts differ");
  }

  // ── Euler angles (N × 3)
  // ─────────────────────────────────────────────────────
  int nRowsEuler = 0, nColsEuler = 0;
  const std::vector<double> eulerFlat = readMultiColumn(sep + "EulerAngles.csv", nRowsEuler, nColsEuler);
  if(nRowsEuler != N || nColsEuler < 3)
  {
    throw std::runtime_error("MTRDataLoader: EulerAngles.csv must have N rows × 3 columns");
  }

  // ── Parent IDs
  // ───────────────────────────────────────────────────────────────
  const std::vector<double> parentVals = readColumn(sep + "ParentIds.csv");
  if(static_cast<int>(parentVals.size()) != N)
  {
    throw std::runtime_error("MTRDataLoader: ParentIds.csv row count differs from N");
  }

  // ── MTR boolean mask
  // ─────────────────────────────────────────────────────────
  const std::vector<double> boolVals = readColumn(sep + "BoolMTR.csv");
  if(static_cast<int>(boolVals.size()) != N)
  {
    throw std::runtime_error("MTRDataLoader: BoolMTR.csv row count differs from N");
  }

  // ── Pack into EBSDData
  // ───────────────────────────────────────────────────────
  EBSDData data;

  data.spatialCoords.resize(N, 2);
  for(int i = 0; i < N; ++i)
  {
    data.spatialCoords(i, 0) = xVals[static_cast<std::size_t>(i)];
    data.spatialCoords(i, 1) = yVals[static_cast<std::size_t>(i)];
  }

  data.eulerAngles.resize(N, 3);
  for(int i = 0; i < N; ++i)
  {
    const std::size_t base = static_cast<std::size_t>(i * nColsEuler);
    data.eulerAngles(i, 0) = eulerFlat[base];
    data.eulerAngles(i, 1) = eulerFlat[base + 1];
    data.eulerAngles(i, 2) = eulerFlat[base + 2];
  }

  data.parentIds.resize(N);
  data.isMTR.resize(N);
  for(int i = 0; i < N; ++i)
  {
    const std::size_t si = static_cast<std::size_t>(i);
    data.parentIds[i] = static_cast<int>(parentVals[si]);
    data.isMTR[i] = (boolVals[si] != 0.0);
  }

  return data;
}

} // namespace mtrsim
