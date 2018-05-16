#pragma once

#include "tinyhhg_core/solvers/preconditioners/IdentityPreconditioner.hpp"

namespace hhg
{

template<class F, class O, class Preconditioner = IdentityPreconditioner<F> >
class CGSolver
{
public:

  CGSolver(const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel,
           uint_t restartFrequency = std::numeric_limits<uint_t>::max(),
           std::shared_ptr<Preconditioner> prec_ = std::make_shared<Preconditioner>())
    : p("p", storage, minLevel, maxLevel), z("z", storage, minLevel, maxLevel), ap("ap", storage, minLevel, maxLevel),
      restartFrequency_(restartFrequency), prec(prec_)
  {
  }

  void init(O& A, F& x, F& b, F& r, size_t level, DoFType flag) {
    A.apply(x, p, level, flag, Replace);
    r.assign({1.0, -1.0}, {&b, &p}, level, flag);
    prec->apply(r, z, level, flag);
    p.assign({1.0}, {&z}, level, flag);
    prsold = r.dotGlobal(z, level, flag);
  }

  void solve(O& A, F& x, F& b, F& r, size_t level, real_t tolerance, size_t maxiter, DoFType flag = All, bool printInfo = false)
  {
    init(A, x, b, r, level, flag);
    res_start = std::sqrt(r.dotGlobal(r, level, flag));

    if (res_start < tolerance)
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
      pAp = p.dotGlobal(ap, level, flag);

      alpha = prsold / pAp;
      x.add({alpha}, {&p}, level, flag);
      r.add({ -alpha }, { &ap }, level, flag);
      rsnew = r.dotGlobal(r, level, flag);
      sqrsnew = std::sqrt(rsnew);

      if (printInfo)
      {
        WALBERLA_LOG_INFO_ON_ROOT("[CG] residual: " << sqrsnew);
      }

      if (sqrsnew/res_start < tolerance)
      {
        if (printInfo)
        {
          WALBERLA_LOG_INFO_ON_ROOT("[CG] converged after " << i << " iterations");
        }
        break;
      }

      prec->apply(r, z, level, flag);
      prsnew = r.dotGlobal(z, level, flag);
      beta = prsnew / prsold;

      p.assign({1.0, beta}, {&z, &p}, level, flag);
      prsold = prsnew;

      if (i % restartFrequency_ == 0) {
        init(A, x, b, r, level, flag);
        if (printInfo)
        {
          WALBERLA_LOG_INFO_ON_ROOT("[CG] restarted");
        }
      }
    }
  }

private:
  F p;
  F z;
  F ap;
  std::shared_ptr<Preconditioner> prec;

  real_t prsold;
  real_t res_start;
  uint_t restartFrequency_;
};

}
