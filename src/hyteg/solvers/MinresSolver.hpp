/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Marcus Mohr, Nils Kohl, Andreas Burkhart.
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

#include "hyteg/Format.hpp"
#include "hyteg/functions/FunctionTools.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

namespace hyteg {

template < class OperatorType >
class MinResSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   MinResSolver(
       const std::shared_ptr< PrimitiveStorage >& storage,
       size_t                                     minLevel,
       size_t                                     maxLevel,
       uint_t                                     maxIter           = std::numeric_limits< uint_t >::max(),
       real_t                                     relativeTolerance = 1e-16,
       real_t                                     absoluteTolerance = 1e-16,
       std::shared_ptr< Solver< OperatorType > >  preconditioner = std::make_shared< IdentityPreconditioner< OperatorType > >(),
       bool                                       lowMemoryMode  = false,
       std::function<
           bool( uint_t, const OperatorType&, const FunctionType&, const FunctionType&, uint_t, real_t ) > iterationHook =
           []( uint_t it, const OperatorType& A, const FunctionType& x, const FunctionType& b, uint_t level, real_t initRes ) {
              WALBERLA_UNUSED( it );
              WALBERLA_UNUSED( A );
              WALBERLA_UNUSED( x );
              WALBERLA_UNUSED( b );
              WALBERLA_UNUSED( level );
              WALBERLA_UNUSED( initRes );
              return false;
           } )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , maxIter_( maxIter )
   , iterations_( maxIter_ )
   , relativeTolerance_( relativeTolerance )
   , absoluteTolerance_( absoluteTolerance )
   , printInfo_( false )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   , preconditioner_( preconditioner )
   , timingTree_( storage->getTimingTree() )
   , iterationHook_( iterationHook )
   , name_( "MinRes" )
   , lowMemoryMode_( lowMemoryMode )
   {
      if ( !lowMemoryMode_ )
      {
         p_vm  = std::make_shared< FunctionType >( "minres_vm", storage, minLevel, maxLevel );
         p_v   = std::make_shared< FunctionType >( "minres_v", storage, minLevel, maxLevel );
         p_vp  = std::make_shared< FunctionType >( "minres_vp", storage, minLevel, maxLevel );
         p_z   = std::make_shared< FunctionType >( "minres_z", storage, minLevel, maxLevel );
         p_zp  = std::make_shared< FunctionType >( "minres_zp", storage, minLevel, maxLevel );
         p_wm  = std::make_shared< FunctionType >( "minres_wm", storage, minLevel, maxLevel );
         p_w   = std::make_shared< FunctionType >( "minres_w", storage, minLevel, maxLevel );
         p_wp  = std::make_shared< FunctionType >( "minres_wp", storage, minLevel, maxLevel );
         p_tmp = std::make_shared< FunctionType >( "minres_tmp", storage, minLevel, maxLevel );
         r_    = std::make_shared< FunctionType >( "minres_r", storage, minLevel, maxLevel );
      }
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      timingTree_->start( "MinRes Solver" );

      std::shared_ptr< FunctionType > p_vmSolve;
      std::shared_ptr< FunctionType > p_vSolve;
      std::shared_ptr< FunctionType > p_vpSolve;
      std::shared_ptr< FunctionType > p_zSolve;
      std::shared_ptr< FunctionType > p_zpSolve;
      std::shared_ptr< FunctionType > p_wmSolve;
      std::shared_ptr< FunctionType > p_wSolve;
      std::shared_ptr< FunctionType > p_wpSolve;
      std::shared_ptr< FunctionType > p_tmpSolve;
      std::shared_ptr< FunctionType > r_Solve;

      if ( !lowMemoryMode_ )
      {
         p_vmSolve  = p_vm;
         p_vSolve   = p_v;
         p_vpSolve  = p_vp;
         p_zSolve   = p_z;
         p_zpSolve  = p_zp;
         p_wmSolve  = p_wm;
         p_wSolve   = p_w;
         p_wpSolve  = p_wp;
         p_tmpSolve = p_tmp;
         r_Solve    = r_;
      }
      else
      {
         p_vmSolve  = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         p_vSolve   = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         p_vpSolve  = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         p_zSolve   = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         p_zpSolve  = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         p_wmSolve  = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         p_wSolve   = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         p_wpSolve  = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         p_tmpSolve = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         r_Solve    = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
      }

      copyBCs( x, *p_vmSolve );
      copyBCs( x, *p_vSolve );
      copyBCs( x, *p_vpSolve );
      copyBCs( x, *p_zSolve );
      copyBCs( x, *p_zpSolve );
      copyBCs( x, *p_wmSolve );
      copyBCs( x, *p_wSolve );
      copyBCs( x, *p_wpSolve );
      copyBCs( x, *p_tmpSolve );
      copyBCs( x, *r_Solve );

      p_vmSolve->setToZero( level );
      p_vSolve->setToZero( level );
      p_vpSolve->setToZero( level );
      p_zSolve->setToZero( level );
      p_zpSolve->setToZero( level );
      p_wmSolve->setToZero( level );
      p_wSolve->setToZero( level );
      p_wpSolve->setToZero( level );
      p_tmpSolve->setToZero( level );
      r_Solve->setToZero( level );

      A.apply( x, *r_Solve, level, flag_ );
      p_vSolve->assign( { 1.0, -1.0 }, { b, *r_Solve }, level, flag_ );

      preconditioner_->solve( A, *p_zSolve, *p_vSolve, level );

      real_t gamma_old = 1.0;
      real_t gamma_new = std::sqrt( p_zSolve->dotGlobal( *p_vSolve, level, flag_ ) );

      real_t res_start = gamma_new;

      real_t eta   = gamma_new;
      real_t s_old = 0.0;
      real_t s_new = 0.0;

      real_t c_old = 1.0;
      real_t c_new = 1.0;

      if ( printInfo_ )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "[" << name_ << "] residuum: " << std::scientific << std::abs( gamma_new ) );
      }

      if ( gamma_new < absoluteTolerance_ )
      {
         if ( printInfo_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "[" << name_ << "] converged" );
         }
         timingTree_->stop( "MinRes Solver" );
         return;
      }

      iterations_ = maxIter_;
      for ( size_t i = 0; i < maxIter_; ++i )
      {
         p_zSolve->assign( { real_t( 1 ) / gamma_new }, { *p_zSolve }, level, flag_ );
         A.apply( *p_zSolve, *p_vpSolve, level, flag_ );
         real_t delta = p_vpSolve->dotGlobal( *p_zSolve, level, flag_ );

         p_vpSolve->assign(
             { 1.0, -delta / gamma_new, -gamma_new / gamma_old }, { *p_vpSolve, *p_vSolve, *p_vmSolve }, level, flag_ );

         p_zpSolve->interpolate( 0, level, flag_ );
         preconditioner_->solve( A, *p_zpSolve, *p_vpSolve, level );

         gamma_old = gamma_new;
         gamma_new = std::sqrt( p_zpSolve->dotGlobal( *p_vpSolve, level, flag_ ) );

         real_t alpha0 = c_new * delta - c_old * s_new * gamma_old;
         real_t alpha1 = std::sqrt( alpha0 * alpha0 + gamma_new * gamma_new );
         real_t alpha2 = s_new * delta + c_old * c_new * gamma_old;
         real_t alpha3 = s_old * gamma_old;

         c_old = c_new;
         c_new = alpha0 / alpha1;
         s_old = s_new;
         s_new = gamma_new / alpha1;

         p_wpSolve->assign(
             { real_t( 1 ) / alpha1, -alpha3 / alpha1, -alpha2 / alpha1 }, { *p_zSolve, *p_wmSolve, *p_wSolve }, level, flag_ );
         x.add( { c_new * eta }, { *p_wpSolve }, level, flag_ );

         eta = -s_new * eta;

         p_tmpSolve->swap( *p_vpSolve, level );
         p_vpSolve->swap( *p_vmSolve, level );
         p_vmSolve->swap( *p_vSolve, level );
         p_vSolve->swap( *p_tmpSolve, level );

         p_tmpSolve->swap( *p_wpSolve, level );
         p_wpSolve->swap( *p_wmSolve, level );
         p_wmSolve->swap( *p_wSolve, level );
         p_wSolve->swap( *p_tmpSolve, level );

         p_tmpSolve->swap( *p_zpSolve, level );
         p_zpSolve->swap( *p_zSolve, level );
         p_zSolve->swap( *p_tmpSolve, level );

         if ( iterationHook_( i, A, x, b, level, res_start ) )
         {
            if ( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[" << name_ << "] iteration hook stopped the solver after " << std::defaultfloat << i
                                              << " iterations" );
            }
            iterations_ = i;
            break;
         }

         if ( printInfo_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
                std::string( "[" ) + name_ + std::string( "] iter: %6d | residuum: %10.5e" ), i, std::abs( eta ) ) );
         }

         if ( std::abs( eta ) / res_start < relativeTolerance_ || std::abs( eta ) < absoluteTolerance_ )
         {
            if ( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[" << name_ << "] converged after " << std::defaultfloat << i << " iterations" );
            }
            iterations_ = i;
            break;
         }
      }
      timingTree_->stop( "MinRes Solver" );
   }

   uint_t getIterations() { return iterations_; }

   void setPrintInfo( bool printInfo ) { printInfo_ = printInfo; }
   void setName( std::string newName ) { name_ = newName; }
   void setAbsoluteTolerance( real_t absoluteTolerance ) { absoluteTolerance_ = absoluteTolerance; }
   void setRelativeTolerance( real_t relativeTolerance ) { relativeTolerance_ = relativeTolerance; }
   void setIterationHook(
       std::function< bool( uint_t, const OperatorType&, const FunctionType&, const FunctionType&, uint_t, real_t ) >
           iterationHook )
   {
      iterationHook_ = iterationHook;
   };

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   uint_t         maxIter_;
   uint_t         iterations_;
   real_t         relativeTolerance_;
   real_t         absoluteTolerance_;
   bool           printInfo_;
   hyteg::DoFType flag_;

   std::shared_ptr< Solver< OperatorType > > preconditioner_;

   std::shared_ptr< FunctionType > p_vm;
   std::shared_ptr< FunctionType > p_v;
   std::shared_ptr< FunctionType > p_vp;

   std::shared_ptr< FunctionType > p_z;
   std::shared_ptr< FunctionType > p_zp;

   std::shared_ptr< FunctionType > p_wm;
   std::shared_ptr< FunctionType > p_w;
   std::shared_ptr< FunctionType > p_wp;

   std::shared_ptr< FunctionType > p_tmp;

   std::shared_ptr< FunctionType > r_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;

   std::function< bool( uint_t, const OperatorType&, const FunctionType&, const FunctionType&, uint_t, real_t ) > iterationHook_;

   std::string name_;

   bool lowMemoryMode_;
};

} // namespace hyteg
