#include "MTRSimDriver.hpp"

#include <numbers>

namespace mtrsim {

ODFComponent buildUniformODF(int n1, int nPHI, int n2) {
  const int nTotal = n1 * nPHI * n2;
  const double twoPiOverN1 = 2.0 * std::numbers::pi / static_cast<double>(n1);
  const double piOverNPHI = std::numbers::pi / static_cast<double>(nPHI);
  const double twoPiOverN2 = 2.0 * std::numbers::pi / static_cast<double>(n2);

  Eigen::VectorXd phi1Bins(nTotal);
  Eigen::VectorXd phiBins(nTotal);
  Eigen::VectorXd phi2Bins(nTotal);

  for (int ix = 0; ix < nTotal; ++ix) {
    const int i1 = ix / (nPHI * n2);
    const int iPHI = (ix % (nPHI * n2)) / n2;
    const int i2 = ix % n2;
    phi1Bins[ix] = (i1 + 0.5) * twoPiOverN1;
    phiBins[ix] = (iPHI + 0.5) * piOverNPHI;
    phi2Bins[ix] = (i2 + 0.5) * twoPiOverN2;
  }

  ODFComponent uni;
  uni.odfVal = Eigen::VectorXd::Constant(nTotal, 1.0 / static_cast<double>(nTotal));
  uni.phi1Bins = std::move(phi1Bins);
  uni.phiBins = std::move(phiBins);
  uni.phi2Bins = std::move(phi2Bins);
  return uni;
}

ODFComponent gridToODFComponent(const std::vector<double>& values, int n1, int nPHI, int n2, double stepDeg1, double stepDegPHI, double stepDeg2) {
  const int nTotal = n1 * nPHI * n2;
  const double deg2rad = std::numbers::pi / 180.0;
  const double s1 = stepDeg1 * deg2rad;
  const double sP = stepDegPHI * deg2rad;
  const double s2 = stepDeg2 * deg2rad;

  Eigen::VectorXd phi1Bins(nTotal);
  Eigen::VectorXd phiBins(nTotal);
  Eigen::VectorXd phi2Bins(nTotal);
  for (int ix = 0; ix < nTotal; ++ix) {
    const int i1 = ix / (nPHI * n2);
    const int iPHI = (ix % (nPHI * n2)) / n2;
    const int i2 = ix % n2;
    phi1Bins[ix] = (i1 + 0.5) * s1;
    phiBins[ix] = (iPHI + 0.5) * sP;
    phi2Bins[ix] = (i2 + 0.5) * s2;
  }

  ODFComponent c;
  c.odfVal = Eigen::Map<const Eigen::VectorXd>(values.data(), nTotal);
  const double total = c.odfVal.sum();
  if (total > 0.0) {
    c.odfVal /= total;
  }
  c.phi1Bins = std::move(phi1Bins);
  c.phiBins = std::move(phiBins);
  c.phi2Bins = std::move(phi2Bins);
  return c;
}

} // namespace mtrsim
