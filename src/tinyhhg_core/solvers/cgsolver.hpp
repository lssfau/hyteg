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

  template<size_t Level>
  void solve(O& A, F& x, F& b, F& r, walberla::real_t tolerance, size_t maxiter, DoFType flag = All, bool printInfo = false)
  {
    A.template apply<Level>(x, p, flag);
    r.template assign<Level>({1.0, -1.0}, {&b, &p}, flag);
    walberla::real_t res_start = std::sqrt(r.template dot<Level>(r, flag));
    p.template assign<Level>({1.0}, {&r}, flag);
    walberla::real_t rsold = r.template dot<Level>(r, flag);

    if (std::sqrt(rsold) < tolerance && printInfo && walberla::mpi::MPIManager::instance()->rank() == 0)
    {
      fmt::printf("[CG] converged\n");
      return;
    }

    for(size_t i = 0; i < maxiter; ++i)
    {
      A.template apply<Level>(p, ap, flag);
      walberla::real_t pAp = p.template dot<Level>(ap, flag);

      walberla::real_t alpha = rsold / pAp;
      x.template add<Level>({alpha}, {&p}, flag);
      r.template add<Level>({ -alpha }, { &ap }, flag);
      walberla::real_t rsnew = r.template dot<Level>(r, flag);
      walberla::real_t sqrsnew = std::sqrt(rsnew);

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

      p.template assign<Level>({1.0, rsnew/rsold}, {&r, &p}, flag);
      rsold = rsnew;
    }
  }

private:
  F p;
  F ap;

};

}

#endif /* CGSOLVER_HPP */
