#ifndef MINRESSOLVER_HPP
#define MINRESSOLVER_HPP

#include <fmt/format.h>

namespace hhg
{

template<class F, class O>
class MinResSolver
{
public:

  MinResSolver(Mesh& mesh, size_t minLevel, size_t maxLevel)
  {
    p_vm = new F("vm", mesh, minLevel, maxLevel);
    p_v = new F("v", mesh, minLevel, maxLevel);
    p_vp = new F("vp", mesh, minLevel, maxLevel);

    p_z = new F("z", mesh, minLevel, maxLevel);
    p_zp = new F("zp", mesh, minLevel, maxLevel);

    p_wm = new F("wm", mesh, minLevel, maxLevel);
    p_w = new F("w", mesh, minLevel, maxLevel);
    p_wp = new F("wp", mesh, minLevel, maxLevel);
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

  template<size_t Level>
  void init(DoFType flag)
  {
    std::function<walberla::real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
    p_vm->template interpolate<Level>(zero, flag);
    p_wm->template interpolate<Level>(zero, flag);
    p_w->template interpolate<Level>(zero, flag);
  }

  template<size_t Level>
  void solve(O& A, F& x, F& b, F& r, walberla::real_t tolerance, size_t maxiter, DoFType flag = All, bool printInfo = false)
  {
    init<Level>(flag);

    A.template apply<Level>(x, r, flag);
    p_v->template assign<Level>({1.0, -1.0}, {&b, &r}, flag);

    // identity preconditioner
    p_z->template assign<Level>({1.0}, {p_v}, flag);

    walberla::real_t gamma_old = 1.0;
    walberla::real_t gamma_new = std::sqrt(p_z->template dot<Level>(*p_v, flag));

    walberla::real_t res_start = gamma_new;

    walberla::real_t eta = gamma_new;
    walberla::real_t s_old = 0.0;
    walberla::real_t s_new = 0.0;

    walberla::real_t c_old = 1.0;
    walberla::real_t c_new = 1.0;

    if (gamma_new < tolerance && printInfo)
    {
      WALBERLA_LOG_INFO_ON_ROOT("[MinRes] converged");
      return;
    }

    for(size_t i = 0; i < maxiter; ++i) {
      p_z->template assign<Level>({walberla::real_t(1) / gamma_new}, {p_z}, flag);
      A.template apply<Level>(*p_z, *p_vp, flag);
      walberla::real_t delta = p_vp->template dot<Level>(*p_z, flag);

      p_vp->template assign<Level>({1.0, -delta / gamma_new, -gamma_new / gamma_old}, {p_vp, p_v, p_vm}, flag);

      // identity preconditioner
      p_zp->template assign<Level>({1.0}, {p_vp}, flag);

      gamma_old = gamma_new;
      gamma_new = std::sqrt(p_zp->template dot<Level>(*p_vp, flag));

      walberla::real_t alpha0 = c_new * delta - c_old * s_new * gamma_old;
      walberla::real_t alpha1 = std::sqrt(alpha0 * alpha0 + gamma_new * gamma_new);
      walberla::real_t alpha2 = s_new * delta + c_old * c_new * gamma_old;
      walberla::real_t alpha3 = s_old * gamma_old;

      c_old = c_new;
      c_new = alpha0 / alpha1;
      s_old = s_new;
      s_new = gamma_new / alpha1;

      p_wp->template assign<Level>({walberla::real_t(1)/alpha1, -alpha3/alpha1, -alpha2/alpha1}, {p_z, p_wm, p_w}, flag);
      x.template add<Level>({c_new * eta}, {p_wp}, flag);

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
