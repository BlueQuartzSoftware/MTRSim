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

} // namespace mtrsim
