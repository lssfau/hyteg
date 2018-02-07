#pragma once

namespace hhg
{

template<class F, class O, class CoarseSolver>
class GMultigridSolver
{
public:

  enum class CycleType
  {
    VCYCLE,
    WCYCLE
  };

  GMultigridSolver(const std::shared_ptr<PrimitiveStorage> & storage, const std::shared_ptr<CoarseSolver>& coarseSolver, uint_t minLevel, uint_t maxLevel)
    : minLevel_(minLevel), maxLevel_(maxLevel), coarseSolver_(coarseSolver),
      ax_("gmg_ax", storage, minLevel, maxLevel), tmp_("gmg_tmp", storage, minLevel, maxLevel)
  {
    // TODO: remove hardcoded parameters
    nuPre_ = 3;
    nuPost_ = 3;

    zero_ = [](const hhg::Point3D&) { return 0.0; };
  }

  ~GMultigridSolver()
  {
  }

  void solve(O& A, F& x, F& b, F& r, uint_t level, real_t tolerance, size_t maxiter, DoFType flag = All, CycleType cycleType = CycleType::VCYCLE, bool printInfo = false)
  {

    if (level == minLevel_)
    {
      coarseSolver_->solve(A, x, b, r, minLevel_, tolerance, maxiter, flag, printInfo);
    }
    else
    {
      // pre-smooth
      for (size_t i = 0; i < nuPre_; ++i)
      {
        A.smooth_gs(x, b, level, flag);
      }

      A.apply(x, ax_, level, flag);
      r.assign({1.0, -1.0}, { &b, &ax_ }, level, flag);

      // restrict
      r.restrict(level, flag);
      b.assign({1.0}, { &r }, level - 1, flag);

      x.interpolate(zero_, level-1);

      solve(A, x, b, r, level-1, tolerance, maxiter, flag, cycleType, printInfo);

      if (cycleType == CycleType::WCYCLE) {
        solve(A, x, b, r, level-1, tolerance, maxiter, flag, cycleType, printInfo);
      }

      // prolongate
      tmp_.assign({1.0}, { &x }, level, flag);
      x.prolongate(level-1, flag);
      x.add({1.0}, { &tmp_ }, level, flag);

      // post-smooth
      for (size_t i = 0; i < nuPost_; ++i)
      {
        A.smooth_gs(x, b, level, flag);
      }
    }

  }

private:

  uint_t nuPre_;
  uint_t nuPost_;

  uint_t minLevel_;
  uint_t maxLevel_;

  std::shared_ptr<CoarseSolver> coarseSolver_;
  F ax_;
  F tmp_;

  std::function<real_t(const hhg::Point3D&)> zero_;

};

}
