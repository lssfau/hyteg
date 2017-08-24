#ifndef MINRESSOLVER_HPP
#define MINRESSOLVER_HPP

#include <fmt/format.h>

namespace hhg
{

template<class F, class O>
class MinResSolver
{
public:

  MinResSolver(const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel)
  {
    p_vm = new F("vm", storage, minLevel, maxLevel);
    p_v = new F("v", storage, minLevel, maxLevel);
    p_vp = new F("vp", storage, minLevel, maxLevel);

    p_z = new F("z", storage, minLevel, maxLevel);
    p_zp = new F("zp", storage, minLevel, maxLevel);

    p_wm = new F("wm", storage, minLevel, maxLevel);
    p_w = new F("w", storage, minLevel, maxLevel);
    p_wp = new F("wp", storage, minLevel, maxLevel);
  }

  ~MinResSolver()
  {
    delete p_vm;
    delete p_v;
    delete p_vp;

    delete p_z;
    delete p_zp;

    delete p_wm;
    delete p_w;
    delete p_wp;
  }

  void init(size_t level, DoFType flag)
  {
    std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
    p_vm->interpolate(zero, level, flag);
    p_wm->interpolate(zero, level, flag);
    p_w->interpolate(zero, level, flag);
  }

  void solve(O& A, F& x, F& b, F& r, size_t level, real_t tolerance, size_t maxiter, DoFType flag = All, bool printInfo = false)
  {
    init(level, flag);

    A.apply(x, r, level, flag);
    p_v->assign({1.0, -1.0}, {&b, &r}, level, flag);

    // identity preconditioner
    p_z->assign({1.0}, {p_v}, level, flag);

    real_t gamma_old = 1.0;
    real_t gamma_new = std::sqrt(p_z->dot(*p_v, level, flag));

    real_t res_start = gamma_new;

    real_t eta = gamma_new;
    real_t s_old = 0.0;
    real_t s_new = 0.0;

    real_t c_old = 1.0;
    real_t c_new = 1.0;

    if (gamma_new < tolerance && printInfo)
    {
      WALBERLA_LOG_INFO_ON_ROOT("[MinRes] converged");
      return;
    }

    for(size_t i = 0; i < maxiter; ++i) {
      p_z->assign({real_t(1) / gamma_new}, {p_z}, level, flag);
      A.apply(*p_z, *p_vp, level, flag);
      real_t delta = p_vp->dot(*p_z, level, flag);

      p_vp->assign({1.0, -delta / gamma_new, -gamma_new / gamma_old}, {p_vp, p_v, p_vm}, level, flag);

      // identity preconditioner
      p_zp->assign({1.0}, {p_vp}, level, flag);

      gamma_old = gamma_new;
      gamma_new = std::sqrt(p_zp->dot(*p_vp, level, flag));

      real_t alpha0 = c_new * delta - c_old * s_new * gamma_old;
      real_t alpha1 = std::sqrt(alpha0 * alpha0 + gamma_new * gamma_new);
      real_t alpha2 = s_new * delta + c_old * c_new * gamma_old;
      real_t alpha3 = s_old * gamma_old;

      c_old = c_new;
      c_new = alpha0 / alpha1;
      s_old = s_new;
      s_new = gamma_new / alpha1;

      p_wp->assign({real_t(1)/alpha1, -alpha3/alpha1, -alpha2/alpha1}, {p_z, p_wm, p_w}, level, flag);
      x.add({c_new * eta}, {p_wp}, level, flag);

      eta = -s_new * eta;

      p_tmp = p_vp;
      p_vp = p_vm;
      p_vm = p_v;
      p_v = p_tmp;

      p_tmp = p_wp;
      p_wp = p_wm;
      p_wm = p_w;
      p_w = p_tmp;

      p_tmp = p_zp;
      p_zp = p_z;
      p_z = p_tmp;

      if (printInfo)
      {
        WALBERLA_LOG_INFO_ON_ROOT(fmt::format("[MinRes] residuum: {:e}", std::abs(eta)));
      }

      if (std::abs(eta)/res_start < tolerance)
      {
        if (printInfo)
        {
          WALBERLA_LOG_INFO_ON_ROOT(fmt::format("[MinRes] converged after {:d} iterations", i));
        }
        break;
      }
    }
  }

private:

  F* p_vm;
  F* p_v;
  F* p_vp;

  F* p_z;
  F* p_zp;

  F* p_wm;
  F* p_w;
  F* p_wp;

  F* p_tmp;
};

}

#endif /* MINRESSOLVER_HPP */
