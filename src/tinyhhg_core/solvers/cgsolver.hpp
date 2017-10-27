#pragma once

#include <fmt/format.h>

#include "tinyhhg_core/solvers/preconditioners/IdentityPreconditioner.hpp"

namespace hhg
{

template<class F, class O, class Preconditioner = IdentityPreconditioner<F> >
class CGSolver
{
public:

  CGSolver(const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel)
    : p("p", storage, minLevel, maxLevel), z("z", storage, minLevel, maxLevel), ap("ap", storage, minLevel, maxLevel),
      prec(std::make_shared<Preconditioner>())
  {
  }

  CGSolver(const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel, std::shared_ptr<Preconditioner> prec_)
      : p("p", storage, minLevel, maxLevel), z("z", storage, minLevel, maxLevel), ap("ap", storage, minLevel, maxLevel),
        prec(prec_)
  {
  }

  void solve(O& A, F& x, F& b, F& r, size_t level, real_t tolerance, size_t maxiter, DoFType flag = All, bool printInfo = false)
  {
    A.apply(x, p, level, flag, Replace);
    r.assign({1.0, -1.0}, {&b, &p}, level, flag);
    real_t res_start = std::sqrt(r.dot(r, level, flag));
    prec->apply(r, z, level, flag);
    p.assign({1.0}, {&z}, level, flag);
    real_t rsold = r.dot(r, level, flag);
    real_t prsold = r.dot(z, level, flag);

    if (std::sqrt(rsold) < tolerance)
    {
      if (printInfo) {
        WALBERLA_LOG_INFO_ON_ROOT("[CG] converged");
      }
      return;
    }
    real_t pAp, alpha, rsnew, sqrsnew, prsnew, beta;

    for(size_t i = 0; i < maxiter; ++i)
    {
      A.apply(p, ap, level, flag, Replace);
      pAp = p.dot(ap, level, flag);

      alpha = prsold / pAp;
      x.add({alpha}, {&p}, level, flag);
      r.add({ -alpha }, { &ap }, level, flag);
      rsnew = r.dot(r, level, flag);
      sqrsnew = std::sqrt(rsnew);

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

      prec->apply(r, z, level, flag);
      prsnew = r.dot(z, level, flag);
      beta = prsnew / prsold;

      p.assign({1.0, beta}, {&z, &p}, level, flag);
      prsold = prsnew;
    }
  }

private:
  F p;
  F z;
  F ap;
  std::shared_ptr<Preconditioner> prec;
};

}