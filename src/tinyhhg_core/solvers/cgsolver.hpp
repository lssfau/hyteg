#ifndef CGSOLVER_HPP
#define CGSOLVER_HPP

#include "tinyhhg_core/comm.hpp"

#include <fmt/format.h>

namespace hhg
{

template<class F>
class CGSolver
{
public:

  CGSolver(Mesh& mesh, size_t minLevel, size_t maxLevel)
    : p("p", mesh, minLevel, maxLevel), ap("ap", mesh, minLevel, maxLevel)
  {
  }

  void solve(Operator& A, F& x, F& b, F& r, size_t level, double tolerance, size_t maxiter, size_t flag = All, bool printInfo = false)
  {
    x.apply(A, p, level, flag);
    r.assign({1.0, -1.0}, {&b, &p}, level, flag);
    double res_start = std::sqrt(r.dot(r, level, flag));
    p.assign({1.0}, {&r}, level, flag);
    double rsold = r.dot(r, level, flag);

    if (std::sqrt(rsold) < tolerance && printInfo && Comm::get().rk == 0)
    {
      fmt::printf("[CG] converged\n");
      return;
    }

    for(size_t i = 0; i < maxiter; ++i)
    {
      p.apply(A, ap, level, flag);
      double pAp = p.dot(ap, level, flag);

      double alpha = rsold / pAp;
      x.assign({1.0, alpha}, {&x, &p}, level, flag);
      r.add({ -alpha }, { &ap }, level, flag);
      double rsnew = r.dot(r, level, flag);
      double sqrsnew = std::sqrt(rsnew);

      if (printInfo && Comm::get().rk == 0)
      {
        fmt::printf("[CG] residuum: %e\n", sqrsnew);
      }

      if (sqrsnew/res_start < tolerance)
      {
        if (printInfo && Comm::get().rk == 0)
        {
          fmt::printf("[CG] converged after %d iterations\n", i);
        }
        break;
      }

      p.assign({1.0, rsnew/rsold}, {&r, &p}, level, flag);
      rsold = rsnew;
    }
  }

private:
  F p;
  F ap;

};

}

#endif /* CGSOLVER_HPP */