#ifndef CGSOLVER_HPP
#define CGSOLVER_HPP

#include <fmt/format.h>

namespace hhg
{

template<class F, class O>
class CGSolver
{
public:

  CGSolver(Mesh& mesh, size_t minLevel, size_t maxLevel)
    : p("p", mesh, minLevel, maxLevel), ap("ap", mesh, minLevel, maxLevel)
  {
  }

  void solve(O& A, F& x, F& b, F& r, size_t level, real_t tolerance, size_t maxiter, DoFType flag = All, bool printInfo = false)
  {
    A.apply(x, p, level, flag);
    r.assign({1.0, -1.0}, {&b, &p}, level, flag);
    real_t res_start = std::sqrt(r.dot(r, level, flag));
    p.assign({1.0}, {&r}, level, flag);
    real_t rsold = r.dot(r, level, flag);

    if (std::sqrt(rsold) < tolerance && printInfo && walberla::mpi::MPIManager::instance()->rank() == 0)
    {
      fmt::printf("[CG] converged\n");
      return;
    }

    for(size_t i = 0; i < maxiter; ++i)
    {
      A.apply(p, ap, level, flag);
      real_t pAp = p.dot(ap, level, flag);

      real_t alpha = rsold / pAp;
      x.add({alpha}, {&p}, level, flag);
      r.add({ -alpha }, { &ap }, level, flag);
      real_t rsnew = r.dot(r, level, flag);
      real_t sqrsnew = std::sqrt(rsnew);

      if (printInfo && walberla::mpi::MPIManager::instance()->rank() == 0)
      {
        fmt::printf("[CG] residuum: %e\n", sqrsnew);
      }

      if (sqrsnew/res_start < tolerance)
      {
        if (printInfo && walberla::mpi::MPIManager::instance()->rank() == 0)
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
