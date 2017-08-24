#ifndef CGSOLVER_HPP
#define CGSOLVER_HPP

#include <fmt/format.h>

namespace hhg
{

template<class F, class O>
class CGSolver
{
public:

  CGSolver(const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel)
    : p("p", storage, minLevel, maxLevel), ap("ap", storage, minLevel, maxLevel)
  {
  }

  void solve(O& A, F& x, F& b, F& r, size_t level, real_t tolerance, size_t maxiter, DoFType flag = All, bool printInfo = false)
  {
    A.apply(x, p, level, flag, Replace);
    r.assign({1.0, -1.0}, {&b, &p}, level, flag);
    real_t res_start = std::sqrt(r.dot(r, level, flag));
    p.assign({1.0}, {&r}, level, flag);
    real_t rsold = r.dot(r, level, flag);

    if (std::sqrt(rsold) < tolerance && printInfo)
    {
      WALBERLA_LOG_INFO_ON_ROOT("[CG] converged");
      return;
    }

    for(size_t i = 0; i < maxiter; ++i)
    {
      A.apply(p, ap, level, flag, Replace);
      real_t pAp = p.dot(ap, level, flag);

      real_t alpha = rsold / pAp;
      x.add({alpha}, {&p}, level, flag);
      r.add({ -alpha }, { &ap }, level, flag);
      real_t rsnew = r.dot(r, level, flag);
      real_t sqrsnew = std::sqrt(rsnew);

      if (printInfo)
      {
        WALBERLA_LOG_INFO_ON_ROOT(fmt::format("[CG] residual: {}", sqrsnew));
      }

      if (sqrsnew/res_start < tolerance)
      {
        if (printInfo)
        {
          WALBERLA_LOG_INFO_ON_ROOT(fmt::format("[CG] converged after {} iterations", i));
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
