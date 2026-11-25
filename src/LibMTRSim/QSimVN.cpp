// Port of qsimvn.m — Alan Genz quasi-Monte Carlo MVN CDF estimator (1992).
// Original copyright: Alan Genz, All rights reserved (BSD-style, see
// matlab/qsimvn.m). Adaptation: Daniel M. Sparkman 10/01/2013.

#include "QSimVN.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace mtrsim
{

namespace
{

// MATLAB sign(): -1, 0, or +1
inline double matlabSign(double x)
{
  return (x > 0.0) ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
}

// Build quasi-random Halton base vector: sqrt of first (n-1) primes.
// Matches MATLAB: ps = sqrt(primes(5*n*log(n+1)/4)); q = ps(1:n-1)
Eigen::VectorXd primeBasesForDim(int n)
{
  const int need = n - 1;
  if(need <= 0)
  {
    return Eigen::VectorXd{};
  }
  // Upper bound matches MATLAB formula; +20 ensures enough room
  const int limit = std::max(static_cast<int>(5.0 * n * std::log(static_cast<double>(n) + 1.0) / 4.0) + 20, 20);

  std::vector<bool> isComposite(static_cast<std::size_t>(limit + 1), false);
  for(int i = 2; i * i <= limit; ++i)
  {
    if(!isComposite[static_cast<std::size_t>(i)])
    {
      for(int j = i * i; j <= limit; j += i)
      {
        isComposite[static_cast<std::size_t>(j)] = true;
      }
    }
  }

  Eigen::VectorXd result(need);
  int found = 0;
  for(int i = 2; i <= limit && found < need; ++i)
  {
    if(!isComposite[static_cast<std::size_t>(i)])
    {
      result(found++) = std::sqrt(static_cast<double>(i));
    }
  }
  return result;
}

// Inverse normal CDF (probit) — Peter Acklam's rational polynomial
// approximation. Replaces boost::math::erfc_inv with no external dependency.
// Max absolute error < 1.15e-9 over (0, 1).
// Reference: P. J. Acklam, "An algorithm for computing the inverse normal
// cumulative distribution function", 2003.
double normalCDFInverse(double p)
{
  // Coefficients for the central-region rational approximation
  static constexpr double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static constexpr double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
  // Coefficients for the tail rational approximation
  static constexpr double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static constexpr double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};

  constexpr double p_low = 0.02425;
  constexpr double p_high = 1.0 - p_low;

  if(p < p_low)
  {
    // Lower tail
    const double q = std::sqrt(-2.0 * std::log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }
  if(p <= p_high)
  {
    // Central region
    const double q = p - 0.5;
    const double r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
  }
  // Upper tail — mirror of lower tail
  const double q = std::sqrt(-2.0 * std::log(1.0 - p));
  return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
}

} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────

QSimVN::QSimVN(std::mt19937_64& rng)
: m_Rng(rng)
{
}

// ─────────────────────────────────────────────────────────────────────────────
// Static helpers

double QSimVN::phi(double z)
{
  // Phi(z) = P(Z <= z) for Z ~ N(0,1).  MATLAB: erfc(-z/sqrt(2))/2
  return std::erfc(-z / std::sqrt(2.0)) / 2.0;
}

double QSimVN::phiInv(double p)
{
  // Inverse normal CDF.  MATLAB: -sqrt(2)*erfcinv(2*p)  =  normalCDFInverse(p)
  p = std::clamp(p, 1e-15, 1.0 - 1e-15);
  return normalCDFInverse(p);
}

// ─────────────────────────────────────────────────────────────────────────────
// chlrdr — greedy pivoted incomplete Cholesky with bound rescaling.
// Port of the nested function chlrdr() in qsimvn.m.
//
// The MATLAB source at line 173 contains a typo: "bi = (bp(i) -2)/cii"
// The correct expression is "bi = (bp(i) - s)/cii" (uses the accumulated
// inner product s).  The C++ port uses the correct formula.

QSimVN::ChlrdrResult QSimVN::chlrdr(const Eigen::MatrixXd& r, const Eigen::VectorXd& a, const Eigen::VectorXd& b) const
{
  const double ep = 1e-10;
  const int n = static_cast<int>(r.rows());
  const double sqtp = std::sqrt(2.0 * std::numbers::pi);

  Eigen::MatrixXd c = r;
  Eigen::VectorXd ap = a;
  Eigen::VectorXd bp = b;

  // Rescale to correlation matrix: divide each row/col by sqrt of diagonal
  Eigen::VectorXd d = c.diagonal().cwiseMax(0.0).cwiseSqrt();
  for(int i = 0; i < n; ++i)
  {
    if(d(i) > 0.0)
    {
      c.col(i) /= d(i);
      c.row(i) /= d(i);
      ap(i) /= d(i);
      bp(i) /= d(i);
    }
  }

  // y holds the conditional mean used when searching for the optimal pivot
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n);

  for(int k = 0; k < n; ++k)
  {
    // ── Search for the pivot that minimises de = Phi(bi) - Phi(ai) ──────────
    int pivotIdx = k;
    double ckk = 0.0;
    double dem = 1.0;
    double am = 0.0;
    double bm = 0.0;

    for(int i = k; i < n; ++i)
    {
      if(c(i, i) > std::numeric_limits<double>::epsilon())
      {
        const double cii = std::sqrt(std::max(c(i, i), 0.0));
        // s = c(i, 0:k-1) . y(0:k-1)  (accumulated inner product)
        const double s = (k > 0) ? c.row(i).head(k).dot(y.head(k)) : 0.0;
        const double ai = (ap(i) - s) / cii;
        const double bi = (bp(i) - s) / cii; // corrected from MATLAB bug
        const double de = phi(bi) - phi(ai);
        if(de <= dem)
        {
          ckk = cii;
          dem = de;
          am = ai;
          bm = bi;
          pivotIdx = i;
        }
      }
    }

    // ── Symmetric row/col swap: move pivot to position k ────────────────────
    if(pivotIdx > k)
    {
      const int p = pivotIdx;

      std::swap(ap(p), ap(k));
      std::swap(bp(p), bp(k));

      // Preserve old c(k,k) in c(p,p); c(k,k) will be overwritten by ckk below
      c(p, p) = c(k, k);

      // Swap completed Cholesky columns (first k elements of rows k and p)
      if(k > 0)
      {
        Eigen::VectorXd tmp = c.row(p).head(k);
        c.row(p).head(k) = c.row(k).head(k);
        c.row(k).head(k) = tmp;
      }

      // Swap sub-columns strictly below row p (columns k and p, rows p+1..n-1)
      if(p + 1 < n)
      {
        Eigen::VectorXd tmp = c.col(p).tail(n - p - 1);
        c.col(p).tail(n - p - 1) = c.col(k).tail(n - p - 1);
        c.col(k).tail(n - p - 1) = tmp;
      }

      // Swap off-diagonal block between k+1 and p-1
      // c(k+1..p-1, k)  <->  c(p, k+1..p-1)
      if(p - 1 >= k + 1)
      {
        const int segLen = p - k - 1;
        Eigen::VectorXd tmp = c.col(k).segment(k + 1, segLen);
        c.col(k).segment(k + 1, segLen) = c.row(p).segment(k + 1, segLen).transpose();
        c.row(p).segment(k + 1, segLen) = tmp.transpose();
      }
    }

    // ── Cholesky factorisation step ──────────────────────────────────────────
    if(ckk > ep * static_cast<double>(k + 1))
    {
      c(k, k) = ckk;
      if(n - k - 1 > 0)
      {
        c.row(k).tail(n - k - 1).setZero();
      }

      for(int i = k + 1; i < n; ++i)
      {
        c(i, k) /= ckk;
        // Rank-1 Schur complement update: c(i, k+1..i) -= c(i,k) * c(k+1..i, k)
        const int len = i - k;
        c.row(i).segment(k + 1, len) -= c(i, k) * c.col(k).segment(k + 1, len).transpose();
      }

      // Conditional mean for use in future pivot search
      if(std::abs(dem) > ep)
      {
        y(k) = (std::exp(-am * am / 2.0) - std::exp(-bm * bm / 2.0)) / (sqtp * dem);
      }
      else
      {
        if(am < -10.0)
        {
          y(k) = bm;
        }
        else if(bm > 10.0)
        {
          y(k) = am;
        }
        else
        {
          y(k) = (am + bm) / 2.0;
        }
      }
    }
    else
    {
      // Degenerate column — zero out and continue
      c.col(k).segment(k, n - k).setZero();
      y(k) = 0.0;
    }
  }

  return {c, ap, bp};
}

// ─────────────────────────────────────────────────────────────────────────────
// mvndns — inner integrand density evaluation.
// Port of nested function mvndns() in qsimvn.m.

double QSimVN::mvndns(int n, const Eigen::MatrixXd& ch, double ci, double dci, const Eigen::VectorXd& x, const Eigen::VectorXd& a, const Eigen::VectorXd& b) const
{
  const double cn = 37.5;
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n - 1);
  double c = ci;
  double dc = dci;
  double p = dc;

  // MATLAB loop: for i = 2:n  (1-based).  C++ loop variable k = i-1, 0-based.
  for(int k = 1; k < n; ++k)
  {
    y(k - 1) = phiInv(c + x(k - 1) * dc);

    // MATLAB: s = ch(i, 1:i-1) * y(1:i-1)  →  C++: ch.row(k).head(k) .
    // y.head(k)
    const double s = ch.row(k).head(k).dot(y.head(k));
    const double ct = ch(k, k);
    const double ai = a(k) - s;
    const double bi = b(k) - s;

    double newC;
    if(std::abs(ai) < cn * ct)
    {
      newC = phi(ai / ct);
    }
    else
    {
      newC = (1.0 + matlabSign(ai)) / 2.0;
    }

    double newD;
    if(std::abs(bi) < cn * ct)
    {
      newD = phi(bi / ct);
    }
    else
    {
      newD = (1.0 + matlabSign(bi)) / 2.0;
    }

    dc = newD - newC;
    c = newC;
    p *= dc;
  }

  return p;
}

// ─────────────────────────────────────────────────────────────────────────────
// compute — randomised quasi-Monte Carlo outer loop.
// Port of the main body of qsimvn() in qsimvn.m.

std::pair<double, double> QSimVN::compute(int m, const Eigen::MatrixXd& r, const Eigen::VectorXd& a, const Eigen::VectorXd& b)
{
  const int n = static_cast<int>(r.rows());
  auto [ch, as, bs] = chlrdr(r, a, b);

  // Initialise ci, dci from the first dimension
  const double ct0 = ch(0, 0);
  const double a0 = as(0);
  const double b0 = bs(0);
  const double cn = 37.5;

  double ci;
  if(std::abs(a0) < cn * ct0)
  {
    ci = phi(a0 / ct0);
  }
  else
  {
    ci = (1.0 + matlabSign(a0)) / 2.0;
  }

  double di;
  if(std::abs(b0) < cn * ct0)
  {
    di = phi(b0 / ct0);
  }
  else
  {
    di = (1.0 + matlabSign(b0)) / 2.0;
  }

  const double dci = di - ci;
  double p = 0.0;
  double e = 0.0;

  const int ns = 12;
  const int nv = std::max(m / ns, 1); // integer division matches MATLAB floor(m/ns)

  // Quasi-random Halton base vector: sqrt of first (n-1) primes
  const Eigen::VectorXd q = primeBasesForDim(n);

  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  for(int i = 1; i <= ns; ++i)
  {
    double vi = 0.0;

    // Random shift vector for this scramble
    Eigen::VectorXd xr(n - 1);
    for(int idx = 0; idx < n - 1; ++idx)
    {
      xr(idx) = uniform(m_Rng);
    }

    for(int j = 1; j <= nv; ++j)
    {
      // MATLAB: x = abs(2*mod(j*q + xr, 1) - 1)
      Eigen::VectorXd x(n - 1);
      for(int idx = 0; idx < n - 1; ++idx)
      {
        double val = std::fmod(static_cast<double>(j) * q(idx) + xr(idx), 1.0);
        if(val < 0.0)
        {
          val += 1.0;
        }
        x(idx) = std::abs(2.0 * val - 1.0);
      }

      const double vp = mvndns(n, ch, ci, dci, x, as, bs);
      vi += (vp - vi) / static_cast<double>(j);
    }

    const double d = (vi - p) / static_cast<double>(i);
    p += d;

    // Running error estimate (MATLAB lines 104-108)
    if(std::abs(d) > 0.0)
    {
      e = std::abs(d) * std::sqrt(1.0 + (e / d) * (e / d) * static_cast<double>(i - 2) / static_cast<double>(i));
    }
    else if(i > 1)
    {
      e *= std::sqrt(static_cast<double>(i - 2) / static_cast<double>(i));
    }
  }

  e *= 3.0;
  return {p, e};
}

} // namespace mtrsim
