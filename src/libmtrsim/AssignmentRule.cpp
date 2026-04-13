// Port of select_AR.m (otherwise branch) and eval_AR.m.
// Daniel M. Sparkman, 08/27/2013 and 9/19/2013.

#include "AssignmentRule.hpp"

#include "QSimVN.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace mtrsim
{

namespace
{
constexpr double k_Inf = std::numeric_limits<double>::infinity();
constexpr double k_PTol = 1e-6;
constexpr double k_QMin0 = -6.0;
constexpr double k_QMax0 = 6.0;
constexpr int k_QSimVNSamples = 20;
constexpr int k_MaxBisectIter = 2000;

// Bisect q in (qMin, qMax) so that qsimvn(identity, sTmp, tTmp with tTmp[goi]=q) == poi.
// tTmp is passed by value — the caller receives the final tTmp back via the reference returned.
double bisectThreshold(QSimVN& qsimvn, const Eigen::MatrixXd& identity, const Eigen::VectorXd& sTmp, Eigen::VectorXd& tTmp, int goi, double poi)
{
  double q = 0.0;
  double qMin = k_QMin0;
  double qMax = k_QMax0;

  tTmp(goi) = q;
  double pTmp = qsimvn.compute(k_QSimVNSamples, identity, sTmp, tTmp).first;

  for(int iter = 0; iter < k_MaxBisectIter && std::abs(poi - pTmp) > k_PTol; ++iter)
  {
    if(pTmp < poi)
    {
      qMin = q;
    }
    else
    {
      qMax = q;
    }
    q = (qMax + qMin) / 2.0;
    tTmp(goi) = q;
    pTmp = qsimvn.compute(k_QSimVNSamples, identity, sTmp, tTmp).first;
  }

  return q;
}

} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────

AssignmentRule::AssignmentRule(std::mt19937_64& rng)
: m_Rng(rng)
{
}

// ─────────────────────────────────────────────────────────────────────────────
// selectThresholds — port of select_AR.m "otherwise" branch.
//
// The MATLAB switch selects a specific topology based on AR_index = num_components.
// Cases 1, 2, 5 handle exotic 2D topologies.  For the standard 3-component
// simulation (and any N != 1,2,5) the "otherwise" branch runs, implementing
// the general N-component all-touching topology with N-1 independent Gaussians.
//
// The thresholds are chosen by bisection so that each component's target volume
// fraction is matched.  The last component's thresholds are set to cover the
// remaining region — no bisection needed.

AssignmentRuleThresholds AssignmentRule::selectThresholds(const std::vector<double>& volumeFractions)
{
  const int numComponents = static_cast<int>(volumeFractions.size());
  if(numComponents < 1)
  {
    throw std::invalid_argument("AssignmentRule: volumeFractions must not be empty");
  }

  const int numGaussians = numComponents - 1;

  AssignmentRuleThresholds result;
  result.numGaussians = numGaussians;
  result.minThresholds = Eigen::MatrixXd::Zero(numComponents, numGaussians);
  result.maxThresholds = Eigen::MatrixXd::Zero(numComponents, numGaussians);

  // Degenerate case: single component gets everything
  if(numComponents == 1)
  {
    // numGaussians = 0, matrices already have shape [1 x 0], nothing to do
    return result;
  }

  QSimVN qsimvn(m_Rng);
  const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(numGaussians, numGaussians);

  // ── Component 0 (MATLAB coi=1, goi=1) ──────────────────────────────────────
  // s = [-inf, ..., -inf],  t = [inf, ..., inf] except t[0] is bisected.
  {
    const int goi = 0;
    Eigen::VectorXd s = Eigen::VectorXd::Constant(numGaussians, -k_Inf);
    Eigen::VectorXd t = Eigen::VectorXd::Constant(numGaussians, k_Inf);

    Eigen::VectorXd sTmp = s;
    Eigen::VectorXd tTmp = t;

    const double q = bisectThreshold(qsimvn, identity, sTmp, tTmp, goi, volumeFractions[0]);
    t(goi) = q;

    result.minThresholds.row(0) = s.transpose();
    result.maxThresholds.row(0) = t.transpose();
  }

  // ── Intermediate components (MATLAB loop aix = 2:num_components-1) ─────────
  for(int coi = 1; coi < numComponents - 1; ++coi)
  {
    const int goi = coi;

    // s starts all -inf, then the diagonal entries from previous components
    // set s[i] = max_thresholds(i, i) for i = 0..coi-1
    Eigen::VectorXd s = Eigen::VectorXd::Constant(numGaussians, -k_Inf);
    for(int i = 0; i < coi; ++i)
    {
      s(i) = result.maxThresholds(i, i);
    }
    Eigen::VectorXd t = Eigen::VectorXd::Constant(numGaussians, k_Inf);

    Eigen::VectorXd sTmp = s;
    Eigen::VectorXd tTmp = t;

    const double q = bisectThreshold(qsimvn, identity, sTmp, tTmp, goi, volumeFractions[coi]);
    t(goi) = q;

    result.minThresholds.row(coi) = s.transpose();
    result.maxThresholds.row(coi) = t.transpose();
  }

  // ── Last component (MATLAB coi = num_components) ────────────────────────────
  // No bisection: this component covers the remainder of the sample space.
  {
    const int coi = numComponents - 1;
    Eigen::VectorXd s = Eigen::VectorXd::Constant(numGaussians, -k_Inf);
    for(int i = 0; i < coi; ++i)
    {
      s(i) = result.maxThresholds(i, i);
    }
    const Eigen::VectorXd t = Eigen::VectorXd::Constant(numGaussians, k_Inf);

    result.minThresholds.row(coi) = s.transpose();
    result.maxThresholds.row(coi) = t.transpose();
  }

  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// evaluate — port of eval_AR.m.
//
// For each voxel (row of z), a component_list is initialised to {1..numComponents}.
// A component is eliminated (zeroed) if any of its Gaussian conditions is violated.
// The assignment is the maximum surviving component index (1-based, matching MATLAB).

Eigen::VectorXi AssignmentRule::evaluate(const Eigen::MatrixXd& z, const AssignmentRuleThresholds& thresholds) const
{
  const int N = static_cast<int>(z.rows());
  const int numGaussians = static_cast<int>(z.cols());
  const int numComponents = static_cast<int>(thresholds.minThresholds.rows());

  Eigen::VectorXi assignments = Eigen::VectorXi::Zero(N);

  std::vector<int> componentList(static_cast<std::size_t>(numComponents));

  for(int i = 0; i < N; ++i)
  {
    // Initialise component list to {1, 2, ..., numComponents} (1-based)
    for(int j = 0; j < numComponents; ++j)
    {
      componentList[static_cast<std::size_t>(j)] = j + 1;
    }

    // Zero out components that violate any Gaussian threshold
    for(int j = 0; j < numComponents; ++j)
    {
      for(int k = 0; k < numGaussians; ++k)
      {
        if(z(i, k) < thresholds.minThresholds(j, k) || z(i, k) > thresholds.maxThresholds(j, k))
        {
          componentList[static_cast<std::size_t>(j)] = 0;
          break; // once eliminated, no need to check remaining Gaussians for this component
        }
      }
    }

    // Assignment = maximum surviving component (matches MATLAB: max(component_list))
    assignments(i) = *std::max_element(componentList.begin(), componentList.end());
  }

  return assignments;
}

} // namespace mtrsim
