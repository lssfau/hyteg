#ifndef MINRESSOLVER_HPP
#define MINRESSOLVER_HPP

#include <fmt/format.h>

namespace hhg
{

template<class F>
class MinResSolver
{
public:

  MinResSolver(Mesh& mesh, size_t minLevel, size_t maxLevel)
    : v0("v0", mesh, minLevel, maxLevel),
      v1("v1", mesh, minLevel, maxLevel),
      v2("v2", mesh, minLevel, maxLevel),
      w0("w0", mesh, minLevel, maxLevel),
      w1("w1", mesh, minLevel, maxLevel),
      w2("w2", mesh, minLevel, maxLevel),
      ap("ap", mesh, minLevel, maxLevel)
  {
  }

  void init(size_t level, size_t flag)
  {
    std::function<double(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
    v0.interpolate(zero, level, flag);
    w0.interpolate(zero, level, flag);
    w1.interpolate(zero, level, flag);
  }

  void solve(Operator& A, F& x, F& b, F& r, size_t level, double tolerance, size_t maxiter, size_t flag = All, bool printInfo = false)
  {
    init(level, flag);

    x.apply(A, ap, level, flag);
    v1.assign({1.0, -1.0}, {&b, &ap}, level, flag);

    double gamma0 = 1.0;
    double gamma1 = std::sqrt(v1.dot(v1, level, flag));
    double gamma2;

    double eta = gamma1;
    double s0 = 0.0;
    double s1 = 0.0;
    double s2;

    double c0 = 1.0;
    double c1 = 1.0;
    double c2;

    if (std::sqrt(gamma1) < tolerance && printInfo)
    {
      WALBERLA_LOG_INFO_ON_ROOT("[MinRes] converged");
      return;
    }

    for(size_t i = 0; i < maxiter; ++i) {
      v1.assign({1.0 / gamma1}, {&v1}, level, flag);
      v1.apply(A, ap, level, flag);
      double delta = ap.dot(v1, level, flag);

      v2.assign({1.0, -delta / gamma1, -gamma1 / gamma0}, {&ap, &v1, &v0}, level, flag);

      gamma2 = std::sqrt(v2.dot(v2, level, flag));

      double alpha0 = c1 * delta - c0 * s1 * gamma1;
      double alpha1 = std::sqrt(alpha0 * alpha0 + gamma2 * gamma2);
      double alpha2 = s1 * delta + c0 * c1 * gamma1;
      double alpha3 = s0 * gamma1;

      c2 = alpha0 / alpha1;
      s2 = gamma2 / alpha1;

      w2.assign({1.0 / alpha1, -alpha3 / alpha1, -alpha2 / alpha1}, {&v1, &w0, &w1}, level, flag);
      x.add({c2 * eta}, {&w2}, level, flag);

      eta = -s2 * eta;

      s0 = s1;
      s1 = s2;
      c0 = c1;
      c1 = c2;
      gamma0 = gamma1;
      gamma1 = gamma2;
      v0.assign({1.0}, {&v1}, level, flag);
      v1.assign({1.0}, {&v2}, level, flag);
      w0.assign({1.0}, {&w1}, level, flag);
      w1.assign({1.0}, {&w2}, level, flag);

      x.apply(A, ap, level, flag);
      r.assign({1.0, -1.0}, {&b, &ap}, level, flag);

      double res = std::sqrt(r.dot(r, level, flag));
      fmt::printf("res = %e\n", res);
    }
  }

private:
  F v0;
  F v1;
  F v2;
  F w0;
  F w1;
  F w2;
  F ap;

};

}

#endif /* MINRESSOLVER_HPP */