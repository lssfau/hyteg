#pragma once

#include "core/timing/TimingTree.h"

#include "tinyhhg_core/solvers/preconditioners/IdentityPreconditioner.hpp"
#include "tinyhhg_core/solvers/Solver.hpp"

namespace hhg
{

template<class OperatorType >
class MinResSolver : public Solver< OperatorType >
{
public:

  typedef typename OperatorType::srcType FunctionType;

  MinResSolver(
      const std::shared_ptr< PrimitiveStorage >& storage,
      size_t                                     minLevel,
      size_t                                     maxLevel,
      uint_t                                     maxIter        = std::numeric_limits< uint_t >::max(),
      real_t                                     tolerance      = 1e-16,
      std::shared_ptr< Solver< OperatorType > >  preconditioner = std::make_shared< IdentityPreconditioner< OperatorType > >() )
  : maxIter_( maxIter )
  , tolerance_( tolerance )
  , printInfo_( false )
  , flag_( hhg::Inner | hhg::NeumannBoundary )
  , preconditioner_( preconditioner )
  , p_vm( "vm", storage, minLevel, maxLevel )
  , p_v( "v", storage, minLevel, maxLevel )
  , p_vp( "vp", storage, minLevel, maxLevel )
  , p_z( "z", storage, minLevel, maxLevel )
  , p_zp( "zp", storage, minLevel, maxLevel )
  , p_wm( "wm", storage, minLevel, maxLevel )
  , p_w( "w", storage, minLevel, maxLevel )
  , p_wp( "wp", storage, minLevel, maxLevel )
  , p_tmp( "tmp", storage, minLevel, maxLevel )
  , r_( "r", storage, minLevel, maxLevel )
  , timingTree_( storage->getTimingTree() )
  {}

  void solve( const OperatorType& A,const FunctionType& x, const FunctionType& b, const uint_t level ) override
  {
    timingTree_->start( "MinRes Solver" );
    std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
    p_vm.interpolate(zero, level, flag_);
    p_wm.interpolate(zero, level, flag_);
    p_w.interpolate(zero, level, flag_);

    A.apply(x, r_, level, flag_);
    p_v.assign({1.0, -1.0}, {b, r_}, level, flag_);

    preconditioner_->solve(A, p_v, p_z, level );

    real_t gamma_old = 1.0;
    real_t gamma_new = std::sqrt(p_z.dotGlobal(p_v, level, flag_));

    real_t res_start = gamma_new;

    real_t eta = gamma_new;
    real_t s_old = 0.0;
    real_t s_new = 0.0;

    real_t c_old = 1.0;
    real_t c_new = 1.0;

    if (printInfo_) {
      WALBERLA_LOG_INFO_ON_ROOT("[MinRes] residuum: "<< std::scientific << std::abs(gamma_new));
    }

    if (gamma_new < tolerance_)
    {
      if (printInfo_) {
        WALBERLA_LOG_INFO_ON_ROOT("[MinRes] converged");
      }
      timingTree_->stop( "MinRes Solver" );
      return;
    }

    for(size_t i = 0; i < maxIter_; ++i) {
      p_z.assign({real_t(1) / gamma_new}, {p_z}, level, flag_);
      A.apply(p_z, p_vp, level, flag_);
      real_t delta = p_vp.dotGlobal(p_z, level, flag_);

      p_vp.assign({1.0, -delta / gamma_new, -gamma_new / gamma_old}, {p_vp, p_v, p_vm}, level, flag_);

      preconditioner_->solve(A, p_vp, p_zp, level );

      gamma_old = gamma_new;
      gamma_new = std::sqrt(p_zp.dotGlobal(p_vp, level, flag_));

      real_t alpha0 = c_new * delta - c_old * s_new * gamma_old;
      real_t alpha1 = std::sqrt(alpha0 * alpha0 + gamma_new * gamma_new);
      real_t alpha2 = s_new * delta + c_old * c_new * gamma_old;
      real_t alpha3 = s_old * gamma_old;

      c_old = c_new;
      c_new = alpha0 / alpha1;
      s_old = s_new;
      s_new = gamma_new / alpha1;

      p_wp.assign({real_t(1)/alpha1, -alpha3/alpha1, -alpha2/alpha1}, {p_z, p_wm, p_w}, level, flag_);
      x.add({c_new * eta}, {p_wp}, level, flag_);

      eta = -s_new * eta;

      p_tmp.assign( {1.0}, {p_vp}, level );
      p_vp.assign( {1.0}, {p_vm}, level );
      p_vm.assign( {1.0}, {p_v}, level );
      p_v.assign( {1.0}, {p_tmp}, level );

      p_tmp.assign( {1.0}, {p_wp}, level );
      p_wp.assign( {1.0}, {p_wm}, level );
      p_wm.assign( {1.0}, {p_w}, level );
      p_w.assign( {1.0}, {p_tmp}, level );

      p_tmp.assign( {1.0}, {p_zp}, level );
      p_zp.assign( {1.0}, {p_z}, level );
      p_z.assign( {1.0}, {p_tmp}, level );

      if (printInfo_)
      {
        WALBERLA_LOG_INFO_ON_ROOT("[MinRes] residuum: " << std::scientific << std::abs(eta));
      }

      if (std::abs(eta)/res_start < tolerance_)
      {
        if (printInfo_)
        {
          WALBERLA_LOG_INFO_ON_ROOT("[MinRes] converged after " << std::defaultfloat <<  i << " iterations");
        }
        break;
      }
    }
    timingTree_->stop( "MinRes Solver" );
  }

private:

  uint_t       maxIter_;
  real_t       tolerance_;
  bool         printInfo_;
  hhg::DoFType flag_;

  std::shared_ptr< Solver< OperatorType > > preconditioner_;

  FunctionType p_vm;
  FunctionType p_v;
  FunctionType p_vp;

  FunctionType p_z;
  FunctionType p_zp;

  FunctionType p_wm;
  FunctionType p_w;
  FunctionType p_wp;

  FunctionType p_tmp;

  FunctionType r_;

  std::shared_ptr< walberla::WcTimingTree > timingTree_;
};

}
