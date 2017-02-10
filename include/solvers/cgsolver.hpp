#ifndef CGSOLVER_HPP
#define CGSOLVER_HPP

#include <fmt/format.h>
#include <comm.hpp>

namespace hhg
{

class CGSolver
{
public:

  CGSolver(P1FunctionSpace& fspace, size_t minLevel, size_t maxLevel)
    : p("p", fspace, minLevel, maxLevel), ap("ap", fspace, minLevel, maxLevel)
  {
  }

  void solve(P1LaplaceOperator& A, Function<P1FunctionSpace>& x, Function<P1FunctionSpace>& b, Function<P1FunctionSpace>& r, size_t level, double tolerance, size_t maxiter, size_t flag = All, bool printInfo = false)
  {
    x.apply(A, p, level, flag);
    r.assign<2>({1.0, -1.0}, {&b, &p}, level, flag);
    double res_start = std::sqrt(r.dot(r, level, flag));
    p.assign<1>({1.0}, {&r}, level, flag);
    double rsold = r.dot(r, level, flag);

    for(size_t i = 0; i < maxiter; ++i)
    {
      p.apply(A, ap, level, flag);
      double alpha = rsold / p.dot(ap, level, flag);
      x.assign<2>({1.0, alpha}, {&x, &p}, level, flag);
      r.add<1>({ -alpha }, { &ap }, level, flag);
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

      p.assign<2>({1.0, rsnew/rsold}, {&r, &p}, level, flag);
      rsold = rsnew;
    }
  }

private:
  Function<P1FunctionSpace> p;
  Function<P1FunctionSpace> ap;

};

}

#endif /* CGSOLVER_HPP */