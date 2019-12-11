/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#define UZAWA_OLD_VARIANT 0

#include "hyteg/composites/StokesOperatorTraits.hpp"
#include "hyteg/composites/P2P1UzawaDampingFactorEstimationOperator.hpp"
#include "hyteg/numerictools/SpectrumEstimation.hpp"
#include "hyteg/solvers/Solver.hpp"

#include "core/math/Random.h"

namespace hyteg {

template < class OperatorType >
class UzawaSmoother : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   UzawaSmoother( const std::shared_ptr< PrimitiveStorage >& storage,
                  const uint_t                               minLevel,
                  const uint_t                               maxLevel,
                  real_t                                     relaxParam,
                  hyteg::DoFType                             flag                    = hyteg::Inner | hyteg::NeumannBoundary,
                  const bool                                 symmetricGSVelocity     = false,
                  const uint_t                               numGSIterationsVelocity = 2,
                  const bool                                 symmetricGSPressure     = false,
                  const uint_t                               numGSIterationsPressure = 1,
                  const real_t                               velocityRelaxParam      = 1.0 )
   : UzawaSmoother( storage,
                    FunctionType( "uzawa_smoother_r", storage, minLevel, maxLevel ),
                    minLevel,
                    maxLevel,
                    relaxParam,
                    flag,
                    symmetricGSVelocity,
                    numGSIterationsVelocity,
                    symmetricGSPressure,
                    numGSIterationsPressure )
   {}

   UzawaSmoother( const std::shared_ptr< PrimitiveStorage >& storage,
                  const FunctionType&                        tmpFunction,
                  const uint_t                               minLevel,
                  const uint_t                               maxLevel,
                  real_t                                     relaxParam,
                  hyteg::DoFType                               flag                    = hyteg::Inner | hyteg::NeumannBoundary,
                  const bool                                 symmetricGSVelocity     = false,
                  const uint_t                               numGSIterationsVelocity = 2,
                  const bool                                 symmetricGSPressure     = false,
                  const uint_t                               numGSIterationsPressure = 1,
                  const real_t                               velocityRelaxParam      = 1.0 )
   : storage_( storage )
   , flag_( flag )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , relaxParam_( relaxParam )
   , symmetricGSVelocity_( symmetricGSVelocity )
   , numGSIterationsVelocity_( numGSIterationsVelocity )
   , symmetricGSPressure_( symmetricGSPressure )
   , numGSIterationsPressure_( numGSIterationsPressure )
   , velocityRelaxParam_( velocityRelaxParam )
   , r_( tmpFunction )
#if UZAWA_OLD_VARIANT
   , tmp_( "uzawa_smoother_tmp", storage, minLevel, maxLevel )
#endif
   {}

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const uint_t                          level ) override
   {
      uzawaSmooth( A,
                   x,
                   b,
                   level,
                   std::integral_constant< bool, tensor_variant< OperatorType >::value >(),
                   std::integral_constant< bool, has_pspg_block< OperatorType >::value >() );
   }

   real_t estimateAndSetRelaxationParameter( const OperatorType& A, const uint_t & level, const uint_t & numPowerIterations )
   {
      const bool isStableDiscretization = std::is_same< OperatorType, P2P1TaylorHoodStokesOperator >::value;
      WALBERLA_CHECK( isStableDiscretization, "Relaxation parameter estimation only implemented for P2-P1" );
      P2P1UzawaDampingFactorEstimationOperator estimator( storage_, level, level, symmetricGSVelocity_, numGSIterationsVelocity_ );
      P1Function< real_t > iterationVector( "iterationVector", storage_, level, level );
      walberla::math::seedRandomGenerator( 42 );
      auto randFunction = []( const Point3D & ) {
         return walberla::math::realRandom();
      };
      iterationVector.interpolate( randFunction, level, All );
      const real_t estimatedRelaxationParameter = estimateSpectralRadiusWithPowerIteration( estimator, iterationVector, numPowerIterations, storage_, level );
      relaxParam_ = estimatedRelaxationParameter;
      return estimatedRelaxationParameter;
   }

   void setRelaxationParameter( const real_t & omega ) { relaxParam_ = omega; }

 private:
   // Block-Laplace variant
   void uzawaSmooth( const OperatorType& A,
                     const FunctionType& x,
                     const FunctionType& b,
                     const uint_t        level,
                     std::false_type /* tensor */,
                     std::true_type /* PSPG */ ) const
   {
      A.divT_x.apply( x.p, r_.u, level, flag_, Replace );
      r_.u.assign( {1.0, -1.0}, {b.u, r_.u}, level, flag_ );
      for ( uint_t i = 0; i < numGSIterationsVelocity_; i++ )
      {
         if ( symmetricGSVelocity_ )
         {
            A.A.smooth_gs( x.u, r_.u, level, flag_ );
            A.A.smooth_gs_backwards( x.u, r_.u, level, flag_ );
         }
         else
         {
            A.A.smooth_gs( x.u, r_.u, level, flag_ );
         }
      }

      A.divT_y.apply( x.p, r_.v, level, flag_, Replace );
      r_.v.assign( {1.0, -1.0}, {b.v, r_.v}, level, flag_ );
      for ( uint_t i = 0; i < numGSIterationsVelocity_; i++ )
      {
         if ( symmetricGSVelocity_ )
         {
            A.A.smooth_gs( x.v, r_.v, level, flag_ );
            A.A.smooth_gs_backwards( x.v, r_.v, level, flag_ );
         }
         else
         {
            A.A.smooth_gs( x.v, r_.v, level, flag_ );
         }
      }

      if ( hasGlobalCells_ )
      {
         A.divT_z.apply( x.p, r_.w, level, flag_, Replace );
         r_.w.assign( {1.0, -1.0}, {b.w, r_.w}, level, flag_ );
         for ( uint_t i = 0; i < numGSIterationsVelocity_; i++ )
         {
            if ( symmetricGSVelocity_ )
            {
               A.A.smooth_gs( x.w, r_.w, level, flag_ );
               A.A.smooth_gs_backwards( x.w, r_.w, level, flag_ );
            }
            else
            {
               A.A.smooth_gs( x.w, r_.w, level, flag_ );
            }
         }
      }

      A.div_x.apply( x.u, r_.p, level, flag_, Replace );
      A.div_y.apply( x.v, r_.p, level, flag_, Add );

      if ( hasGlobalCells_ )
      {
         A.div_z.apply( x.w, r_.p, level, flag_, Add );
      }

      r_.p.assign( {1.0, -1.0}, {b.p, r_.p}, level, flag_ );

      for ( uint_t i = 0; i < numGSIterationsPressure_; i++ )
      {
         if ( symmetricGSPressure_ )
         {
            A.pspg.smooth_sor( x.p, r_.p, relaxParam_, level, flag_ );
            A.pspg.smooth_sor_backwards( x.p, r_.p, relaxParam_, level, flag_ );
         }
         else
         {
            A.pspg.smooth_sor( x.p, r_.p, relaxParam_, level, flag_ );
         }
      }
   }

   // Tensor variant
   void uzawaSmooth( const OperatorType& A,
                     const FunctionType& x,
                     const FunctionType& b,
                     const uint_t        level,
                     std::true_type /* tensor */,
                     std::true_type /* PSPG */ ) const
   {
      A.divT_x.apply( x.p, r_.u, level, flag_, Replace );
      A.A_uv.apply( x.v, r_.u, level, flag_, Add );
      r_.u.assign( {1.0, -1.0}, {b.u, r_.u}, level, flag_ );
      A.A_uu.smooth_gs( x.u, r_.u, level, flag_ );

      A.divT_y.apply( x.p, r_.v, level, flag_, Replace );
      A.A_vu.apply( x.u, r_.v, level, flag_, Add );
      r_.v.assign( {1.0, -1.0}, {b.v, r_.v}, level, flag_ );
      A.A_vv.smooth_gs( x.v, r_.v, level, flag_ );

      A.div_x.apply( x.u, r_.p, level, flag_, Replace );
      A.div_y.apply( x.v, r_.p, level, flag_, Add );

      r_.p.assign( {1.0, -1.0}, {b.p, r_.p}, level, flag_ );

      A.pspg.smooth_sor( x.p, r_.p, relaxParam_, level, flag_ );
   }

   // Block-Laplace variant without stabilization
   void uzawaSmooth( const OperatorType& A,
                     const FunctionType& x,
                     const FunctionType& b,
                     const uint_t        level,
                     std::false_type /* tensor */,
                     std::false_type /* PSPG */ ) const
   {
      A.divT_x.apply( x.p, r_.u, level, flag_, Replace );
      r_.u.assign( {1.0, -1.0}, {b.u, r_.u}, level, flag_ );

      for ( uint_t i = 0; i < numGSIterationsVelocity_; i++ )
      {
         if ( symmetricGSVelocity_ )
         {
            A.A.smooth_sor( x.u, r_.u, velocityRelaxParam_, level, flag_ );
            A.A.smooth_sor_backwards( x.u, r_.u, velocityRelaxParam_, level, flag_ );
         }
         else
         {
            A.A.smooth_sor( x.u, r_.u, velocityRelaxParam_, level, flag_ );
         }
      }

      A.divT_y.apply( x.p, r_.v, level, flag_, Replace );
      r_.v.assign( {1.0, -1.0}, {b.v, r_.v}, level, flag_ );

      for ( uint_t i = 0; i < numGSIterationsVelocity_; i++ )
      {
         if ( symmetricGSVelocity_ )
         {
            A.A.smooth_sor( x.v, r_.v, velocityRelaxParam_, level, flag_ );
            A.A.smooth_sor_backwards( x.v, r_.v, velocityRelaxParam_, level, flag_ );
         }
         else
         {
            A.A.smooth_sor( x.v, r_.v, velocityRelaxParam_, level, flag_ );
         }
      }

      if ( hasGlobalCells_ )
      {
         A.divT_z.apply( x.p, r_.w, level, flag_, Replace );
         r_.w.assign( {1.0, -1.0}, {b.w, r_.w}, level, flag_ );

         for ( uint_t i = 0; i < numGSIterationsVelocity_; i++ )
         {
            if ( symmetricGSVelocity_ )
            {
               A.A.smooth_sor( x.w, r_.w, velocityRelaxParam_, level, flag_ );
               A.A.smooth_sor_backwards( x.w, r_.w, velocityRelaxParam_, level, flag_ );
            }
            else
            {
               A.A.smooth_sor( x.w, r_.w, velocityRelaxParam_, level, flag_ );
            }
         }
      }

      A.div_x.apply( x.u, r_.p, level, flag_, Replace );
      A.div_y.apply( x.v, r_.p, level, flag_, Add );

      if ( hasGlobalCells_ )
      {
         A.div_z.apply( x.w, r_.p, level, flag_, Add );
      }

      r_.p.assign( {1.0, -1.0}, {b.p, r_.p}, level, flag_ );

#if UZAWA_OLD_VARIANT
      tmp_.p.interpolate( 0.0, level );

      for ( uint_t i = 0; i < numGSIterationsPressure_; i++ )
      {
         if ( symmetricGSPressure_ )
         {
            A.pspg_.smooth_sor( tmp_.p, r_.p, relaxParam_, level, flag_ );
            A.pspg_.smooth_sor_backwards( tmp_.p, r_.p, relaxParam_, level, flag_ );
         }
         else
         {
            A.pspg_.smooth_sor( tmp_.p, r_.p, relaxParam_, level, flag_ );
         }
      }

      x.p.add( {1.0}, {tmp_.p}, level, flag_ | DirichletBoundary );
#else
     // This variant is similar to the one published in
     // Gaspar et al. (2014): A simple and efficient segregated smoother for the discrete stokes equations
     // however, we additionally scale Bu with the inverse of the diagonal of the PSPG operator.
     // This is similar to the old variant where we performed one SOR relaxation step on
     // the zero vector and added the result to the solution.
     r_.p.assign( {relaxParam_}, {r_.p}, level, flag_ );
     A.pspg_inv_diag_.apply( r_.p, x.p, level, flag_, Add );
#endif

   }

 private:
   std::shared_ptr< PrimitiveStorage > storage_;
   DoFType flag_;
   bool    hasGlobalCells_;
   real_t  relaxParam_;
   real_t  velocityRelaxParam_;
   bool    symmetricGSVelocity_;
   uint_t  numGSIterationsVelocity_;
   bool    symmetricGSPressure_;
   uint_t  numGSIterationsPressure_;

   FunctionType r_;

#if UZAWA_OLD_VARIANT
   FunctionType tmp_;
#endif
};
} // namespace hyteg
