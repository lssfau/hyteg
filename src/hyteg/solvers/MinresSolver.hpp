/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "core/timing/TimingTree.h"
#include "core/Format.hpp"

#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"
#include "hyteg/solvers/Solver.hpp"

namespace hyteg {

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
  , relativeTolerance_( tolerance )
  , absoluteTolerance_( 1e-16 )
  , printInfo_( false )
  , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
  , preconditioner_( preconditioner )
  , p_vm( "minres_vm", storage, minLevel, maxLevel )
  , p_v( "minres_v", storage, minLevel, maxLevel )
  , p_vp( "minres_vp", storage, minLevel, maxLevel )
  , p_z( "minres_z", storage, minLevel, maxLevel )
  , p_zp( "minres_zp", storage, minLevel, maxLevel )
  , p_wm( "minres_wm", storage, minLevel, maxLevel )
  , p_w( "minres_w", storage, minLevel, maxLevel )
  , p_wp( "minres_wp", storage, minLevel, maxLevel )
  , p_tmp( "minres_tmp", storage, minLevel, maxLevel )
  , r_( "minres_r", storage, minLevel, maxLevel )
  , timingTree_( storage->getTimingTree() )
  , storage_(storage)
  {}

  void solve( const OperatorType& A,const FunctionType& x, const FunctionType& b, const uint_t level ) override
  {
    timingTree_->start( "MinRes Solver" );

    p_vm.copyBoundaryConditionFromFunction( x );
    p_v.copyBoundaryConditionFromFunction( x );
    p_vp.copyBoundaryConditionFromFunction( x );
    p_z.copyBoundaryConditionFromFunction( x );
    p_zp.copyBoundaryConditionFromFunction( x );
    p_wm.copyBoundaryConditionFromFunction( x );
    p_w.copyBoundaryConditionFromFunction( x );
    p_wp.copyBoundaryConditionFromFunction( x );
    p_tmp.copyBoundaryConditionFromFunction( x );
    r_.copyBoundaryConditionFromFunction( x );

    std::function<real_t(const hyteg::Point3D&)> zero = [](const hyteg::Point3D&) { return 0.0; };
    p_vm.interpolate(zero, level, flag_);
    p_wm.interpolate(zero, level, flag_);
    p_w.interpolate(zero, level, flag_);

    A.apply(x, r_, level, flag_);
    p_v.assign({1.0, -1.0}, {b, r_}, level, flag_);

    p_z.interpolate(0, level, flag_);
    preconditioner_->solve(A, p_z, p_v, level );
    
    /*
     VTKOutput vtkOutput(
          "../../output", "MINRES_tmp", storage_ );
       vtkOutput.add( p_z );
     
    vtkOutput.write( level, 0 );
*/
     
 

    real_t gamma_old = 1.0;
    real_t gamma_new = std::sqrt(p_z.dotGlobal(p_v, level, flag_));

    real_t res_start = gamma_new;

    real_t eta = gamma_new;
    real_t s_old = 0.0;
    real_t s_new = 0.0;

    real_t c_old = 1.0;
    real_t c_new = 1.0;

    if (printInfo_) {
      WALBERLA_LOG_INFO_ON_ROOT("[MinRes] residuum: "<< std::scientific << std::setprecision (15) <<  std::abs(gamma_new));
    }

    if (gamma_new < absoluteTolerance_ )
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

      p_zp.interpolate(0, level, flag_);
      preconditioner_->solve(A, p_zp, p_vp, level );

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

      p_tmp.swap( p_vp, level );
      p_vp.swap( p_vm, level );
      p_vm.swap( p_v, level );
      p_v.swap( p_tmp, level );

      p_tmp.swap( p_wp, level );
      p_wp.swap( p_wm, level );
      p_wm.swap( p_w, level );
      p_w.swap( p_tmp, level );

      p_tmp.swap( p_zp, level );
      p_zp.swap( p_z, level );
      p_z.swap( p_tmp, level );

      if (printInfo_)
      {
        WALBERLA_LOG_INFO_ON_ROOT(walberla::format("[MinRes] iter: %6d | residuum: %10.5e", i, std::abs(eta)));
      }

      if (std::abs(eta)/res_start < relativeTolerance_ || std::abs(eta) < absoluteTolerance_ )
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

  void setPrintInfo( bool printInfo ) { printInfo_ = printInfo; }
  void setAbsoluteTolerance( real_t absoluteTolerance ) { absoluteTolerance_ = absoluteTolerance; }
  void setRelativeTolerance( real_t relativeTolerance ) { relativeTolerance_ = relativeTolerance; }

private:

  uint_t       maxIter_;
  real_t         relativeTolerance_;
  real_t         absoluteTolerance_;
  bool         printInfo_;
  hyteg::DoFType flag_;

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
  
   const std::shared_ptr< PrimitiveStorage >      storage_;
};

}
