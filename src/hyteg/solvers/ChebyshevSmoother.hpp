/*
 * Copyright (c) 2020-2025 Andreas Wagner, Marcus Mohr, Andreas Burkhart
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

#include "core/DataTypes.h"

#include "hyteg/functions/FunctionTools.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/numerictools/SpectrumEstimation.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/ApplyInverseDiagonalWrapper.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/solvers/SubstitutePreconditioner.hpp"

#include "Solver.hpp"

namespace hyteg {

template < typename OperatorType, typename ProjectionOperatorType = hyteg::NoOperator, bool solverInterpolateZero = true >
class ChebyshevSmoother : public Solver< OperatorType >
{
 public:
   using FunctionType = typename OperatorType::srcType;
   // using ValueType    = typename FunctionTrait< FunctionType >::ValueType;

   ChebyshevSmoother( const std::shared_ptr< PrimitiveStorage >& storage,
                      size_t                                     minLevel,
                      size_t                                     maxLevel,
                      bool                                       lowMemoryMode  = false,
                      std::shared_ptr< ProjectionOperatorType >  projection     = nullptr,
                      DoFType                                    projectionFlag = FreeslipBoundary )
   : coefficients{}
   , flag_( Inner | NeumannBoundary | FreeslipBoundary )
   , chebyPrec_( nullptr )
   , projection_( projection )
   , projectionFlag_( projectionFlag )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , lowMemoryMode_( lowMemoryMode )
   {
      coefficients.resize( maxLevel_ - minLevel_ + 1 );

      if ( !lowMemoryMode_ )
      {
         tmp1_ = std::make_shared< FunctionType >( "cheb_tmp1", storage_, minLevel_, maxLevel_ );
         tmp2_ = std::make_shared< FunctionType >( "cheb_tmp2", storage_, minLevel_, maxLevel_ );
      }
   }

   ChebyshevSmoother( const FunctionType&                       tmp1,
                      const FunctionType&                       tmp2,
                      bool                                      lowMemoryMode  = false,
                      std::shared_ptr< ProjectionOperatorType > projection     = nullptr,
                      DoFType                                   projectionFlag = FreeslipBoundary )
   : coefficients{}
   , flag_( Inner | NeumannBoundary | FreeslipBoundary )
   , chebyPrec_( nullptr )
   , projection_( projection )
   , projectionFlag_( projectionFlag )
   , storage_( tmp1.getStorage() )
   , minLevel_( tmp1.getMinLevel() )
   , maxLevel_( tmp1.getMaxLevel() )
   , lowMemoryMode_( lowMemoryMode )
   {
      coefficients.resize( maxLevel_ - minLevel_ + 1 );

      if ( !lowMemoryMode_ )
      {
         tmp1_ = std::make_shared< FunctionType >( "cheb_tmp1", storage_, minLevel_, maxLevel_ );
         tmp2_ = std::make_shared< FunctionType >( "cheb_tmp2", storage_, minLevel_, maxLevel_ );
      }
   }

   template < class SolverOperatorType >
   void setPreconditioner( std::shared_ptr< hyteg::Solver< SolverOperatorType > > preconditioner         = nullptr,
                           std::shared_ptr< SolverOperatorType >                  preconditionerOperator = nullptr )
   {
      if ( ( preconditioner == nullptr ) || ( preconditionerOperator == nullptr ) )
      {
         chebyPrec_ = nullptr;
         WALBERLA_LOG_WARNING( "Chebyshev preconditioner or preconditioner operator is nullptr." );
      }
      else
      {
         chebyPrec_ = std::make_shared< SubstitutePreconditioner< OperatorType, SolverOperatorType > >( preconditioner,
                                                                                                        preconditionerOperator );
      }
   }

   /// Executes an iteration step of the smoother.
   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      if ( chebyPrec_ == nullptr )
      {
         std::shared_ptr< typename OperatorType::srcType > inverseDiagonalValues = nullptr;

         if ( const auto* A_with_inv_diag =
                  dynamic_cast< const OperatorWithInverseDiagonal< typename OperatorType::srcType >* >( &A ) )
         {
            inverseDiagonalValues = A_with_inv_diag->getInverseDiagonalValues();
         }
         else
         {
            throw std::runtime_error( "The Chebyshev-Smoother requires the OperatorWithInverseDiagonal interface." );
         }

         WALBERLA_DEBUG_SECTION()
         {
            const real_t localNormSqr = inverseDiagonalValues->dotLocal( *inverseDiagonalValues, level, flag_ );
            WALBERLA_UNUSED( localNormSqr );
            WALBERLA_ASSERT_GREATER( localNormSqr, 0.0, "diagonal not set" );
         }

         auto ApplyDiag = std::make_shared< ApplyFunctionMultiplicationWrapper< OperatorType, FunctionType > >(
             A.getStorage(), A.getMinLevel(), A.getMaxLevel(), inverseDiagonalValues );

         auto applyPrec =
             std::make_shared< ApplyPreconditioner< ApplyFunctionMultiplicationWrapper< OperatorType, FunctionType > > >(
                 ApplyDiag, flag_ );

         chebyPrec_ = std::make_shared<
             SubstitutePreconditioner< OperatorType, ApplyFunctionMultiplicationWrapper< OperatorType, FunctionType > > >(
             applyPrec, ApplyDiag );
      }

      // handle low memory mode
      std::shared_ptr< FunctionType > tmp1Solve;
      std::shared_ptr< FunctionType > tmp2Solve;

      if ( !lowMemoryMode_ )
      {
         tmp1Solve = tmp1_;
         tmp2Solve = tmp2_;
      }
      else
      {
         tmp1Solve = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         tmp2Solve = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
      }

      copyBCs( x, *tmp1Solve );
      copyBCs( x, *tmp2Solve );

      WALBERLA_DEBUG_SECTION()
      {
         WALBERLA_ASSERT( coefficients.size() > 0, "coefficients must be setup" );
      }

      // ##########################
      // ########## Main ##########
      // ##########################
      WALBERLA_ASSERT_GREATER_EQUAL( level, minLevel_, "level needs to be within minLevel and maxLevel" );
      WALBERLA_ASSERT_LESS_EQUAL( level, maxLevel_, "level needs to be within minLevel and maxLevel" );

      uint_t coefficients_index = level - minLevel_;

      // tmp1_ := Ax
      A.apply( x, *tmp2Solve, level, flag_ );
      // tmp1_ := b-Ax
      tmp2Solve->assign( { real_c( 1 ), real_c( -1 ) }, { b, *tmp2Solve }, level, flag_ );
      // tmp1_ := D^{-1} (b-Ax)
      if constexpr ( solverInterpolateZero )
      {
         tmp1Solve->setToZero( level );
      }
      chebyPrec_->solve( A, *tmp1Solve, *tmp2Solve, level );

      if constexpr ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value )
      {
         if ( projection_ != nullptr )
         {
            projection_->project( *tmp1Solve, level, projectionFlag_ );
         }
      }
      // x := x + omega_0 D^{-1} (b-Ax)
      x.assign( { real_c( 1 ), coefficients[coefficients_index][0] }, { x, *tmp1Solve }, level, flag_ );

      // Loop preconditions before the kth iteration:
      // 1.)  x := x + sum_{i<k} omega_i (D^{-1}A)^{i} D^{-1} (b-Ax)
      // 2.) tmp1_ := (D^{-1}A)^{k-1} D^{-1} (b-Ax)
      for ( uint_t k = 1; k < coefficients[coefficients_index].size(); k += 1 )
      {
         // tmp2_ := A (D^{-1}A)^{k-1} D^{-1} (b-Ax)
         A.apply( *tmp1Solve, *tmp2Solve, level, flag_ );
         // tmp1_ := (D^{-1}A)^{k} D^{-1} (b-Ax)
         if constexpr ( solverInterpolateZero )
         {
            tmp1Solve->setToZero( level );
         }
         chebyPrec_->solve( A, *tmp1Solve, *tmp2Solve, level );

         if constexpr ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value )
         {
            if ( projection_ != nullptr )
            {
               projection_->project( *tmp1Solve, level, projectionFlag_ );
            }
         }

         // x := x + sum_{i<k+1} omega_i (D^{-1}A)^{i} D^{-1} (b-Ax)
         x.assign( { real_c( 1 ), coefficients[coefficients_index][k] }, { x, *tmp1Solve }, level, flag_ );
      }

      if constexpr ( !std::is_same< ProjectionOperatorType, hyteg::NoOperator >::value )
      {
         if ( projection_ != nullptr )
         {
            projection_->project( *tmp1Solve, level, projectionFlag_ );
         }
      }
   }

   /// Calculates the coefficients for our Chebyshev-Smoother.
   /// Has to be called prior to the first usage of the solver.
   ///
   /// This function is based on the mfem implementation at
   ///     https://github.com/mfem/mfem/blob/7cbb0d484863bf661e88378eac4ab0247f10545c/linalg/solvers.cpp#L185
   /// which references the article
   ///     Parallel multigrid smoothing: polynomial versus Gauss-Seidel by Adams et al.
   ///
   /// \param order The order of our polynomial smoother. Only 0 to 9 are supported.
   /// \param spectralRadii An estimate for the spectral radius on every level, given as a vector starting with index 0 representing minLevel.
   /// \param upperFactor The spectral radius is multiplied by this factor to define an upper bound.
   /// \param lowerFactor The spectral radius is multiplied by this factor to define a lower bound.
   void setupCoefficients( const uint_t                order,
                           const std::vector< real_t > spectralRadii,
                           real_t                      upperFactor = real_c( 1.2 ),
                           real_t                      lowerFactor = real_c( 0.3 ) )
   {
      WALBERLA_ASSERT_EQUAL(
          spectralRadii.size(), maxLevel_ - minLevel_ + 1, "The spectral radius vector must be of the correct length" );

      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         const uint_t coefficients_index = level - minLevel_;
         const real_t spectralRadius     = spectralRadii.at( coefficients_index );

         const real_t upperBound = upperFactor * spectralRadius;
         const real_t lowerBound = lowerFactor * spectralRadius;

         setupCoefficientsInternal( order, lowerBound, upperBound, level );
      }
   }

   /// Calculates the coefficients for our Chebyshev-Smoother.
   /// Has to be called prior to the first usage of the solver.
   ///
   /// This function is based on the mfem implementation at
   ///     https://github.com/mfem/mfem/blob/7cbb0d484863bf661e88378eac4ab0247f10545c/linalg/solvers.cpp#L185
   /// which references the article
   ///     Parallel multigrid smoothing: polynomial versus Gauss-Seidel by Adams et al.
   ///
   /// \param order The order of our polynomial smoother. Only 0 to 9 are supported.
   /// \param spectralRadius An estimate for the spectral radius, that will be used on every level.
   /// \param upperFactor The spectral radius is multiplied by this factor to define an upper bound.
   /// \param lowerFactor The spectral radius is multiplied by this factor to define a lower bound.
   void setupCoefficients( const uint_t order,
                           const real_t spectralRadius,
                           real_t       upperFactor = real_c( 1.2 ),
                           real_t       lowerFactor = real_c( 0.3 ) )
   {
      const real_t upperBound = upperFactor * spectralRadius;
      const real_t lowerBound = lowerFactor * spectralRadius;

      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         setupCoefficientsInternal( order, lowerBound, upperBound, level );
      }
   }

   /// Calculates the coefficients for our Chebyshev-Smoother.
   /// Has to be called prior to the first usage of the solver.
   ///
   /// This function is based on the mfem implementation at
   ///     https://github.com/mfem/mfem/blob/7cbb0d484863bf661e88378eac4ab0247f10545c/linalg/solvers.cpp#L185
   /// which references the article
   ///     Parallel multigrid smoothing: polynomial versus Gauss-Seidel by Adams et al.
   ///
   /// \param order The order of our polynomial smoother. Only 0 to 9 are supported.
   /// \param spectralRadius An estimate for the spectral radius, that will be used on every level.
   /// \param level The level on which the given spectral radius is set.
   /// \param upperFactor The spectral radius is multiplied by this factor to define an upper bound.
   /// \param lowerFactor The spectral radius is multiplied by this factor to define a lower bound.
   void setupCoefficientsOnLevel( const uint_t order,
                                  const real_t spectralRadius,
                                  const uint_t level,
                                  real_t       upperFactor = real_c( 1.2 ),
                                  real_t       lowerFactor = real_c( 0.3 ) )
   {
      const real_t upperBound = upperFactor * spectralRadius;
      const real_t lowerBound = lowerFactor * spectralRadius;

      setupCoefficientsInternal( order, lowerBound, upperBound, level );
   }

   /// Calculates the coefficients for our Chebyshev-Smoother.
   /// Has to be called prior to the first usage of the solver.
   ///
   /// This function is based on the mfem implementation at
   ///     https://github.com/mfem/mfem/blob/7cbb0d484863bf661e88378eac4ab0247f10545c/linalg/solvers.cpp#L185
   /// which references the article
   ///     Parallel multigrid smoothing: polynomial versus Gauss-Seidel by Adams et al.
   ///
   /// \param order The order of our polynomial smoother. Only 1 to 5 are supported.
   /// \param lowerBound The smallest eigenvalue whose corresponding eigenvector will be smoothed.
   /// \param upperBound The largest eigenvalue.
   /// \param level The level on which the coefficients are set.
   void setupCoefficientsInternal( const uint_t order, const real_t lowerBound, const real_t upperBound, uint_t level )
   {
      WALBERLA_ASSERT_GREATER( order, 0, "Order cannot be 0." );

      WALBERLA_ASSERT_GREATER_EQUAL( level, minLevel_, "level needs to be within minLevel and maxLevel" );
      WALBERLA_ASSERT_LESS_EQUAL( level, maxLevel_, "level needs to be within minLevel and maxLevel" );

      uint_t coefficients_index = level - minLevel_;

      const double theta = 0.5 * ( upperBound + lowerBound );
      const double delta = 0.5 * ( upperBound - lowerBound );

      coefficients[coefficients_index].resize( order, 0 );

      switch ( order - 1 )
      {
      case 0: {
         coefficients[coefficients_index][0] = 1.0 / theta;
         break;
      }
      case 1: {
         const real_t tmp0 = 1.0 / ( std::pow( delta, 2 ) - 2 * std::pow( theta, 2 ) );

         coefficients[coefficients_index][0] = -4 * theta * tmp0;
         coefficients[coefficients_index][1] = 2 * tmp0;
         break;
      }
      case 2: {
         const real_t tmp0 = 3 * std::pow( delta, 2 );
         const real_t tmp1 = std::pow( theta, 2 );
         const real_t tmp2 = 1.0 / ( -4 * std::pow( theta, 3 ) + theta * tmp0 );

         coefficients[coefficients_index][0] = tmp2 * ( tmp0 - 12 * tmp1 );
         coefficients[coefficients_index][1] = 12 / ( tmp0 - 4 * tmp1 );
         coefficients[coefficients_index][2] = -4 * tmp2;
         break;
      }
      case 3: {
         const real_t tmp0 = std::pow( delta, 2 );
         const real_t tmp1 = std::pow( theta, 2 );
         const real_t tmp2 = 8 * tmp0;
         const real_t tmp3 = 1.0 / ( std::pow( delta, 4 ) + 8 * std::pow( theta, 4 ) - tmp1 * tmp2 );

         coefficients[coefficients_index][0] = tmp3 * ( 32 * std::pow( theta, 3 ) - 16 * theta * tmp0 );
         coefficients[coefficients_index][1] = tmp3 * ( -48 * tmp1 + tmp2 );
         coefficients[coefficients_index][2] = 32 * theta * tmp3;
         coefficients[coefficients_index][3] = -8 * tmp3;
         break;
      }
      case 4: {
         const real_t tmp0 = 5 * std::pow( delta, 4 );
         const real_t tmp1 = std::pow( theta, 4 );
         const real_t tmp2 = std::pow( theta, 2 );
         const real_t tmp3 = std::pow( delta, 2 );
         const real_t tmp4 = 60 * tmp3;
         const real_t tmp5 = 20 * tmp3;
         const real_t tmp6 = 1.0 / ( 16 * std::pow( theta, 5 ) - std::pow( theta, 3 ) * tmp5 + theta * tmp0 );
         const real_t tmp7 = -160 * tmp2;
         const real_t tmp8 = 1.0 / ( tmp0 + 16 * tmp1 - tmp2 * tmp5 );

         coefficients[coefficients_index][0] = tmp6 * ( tmp0 + 80 * tmp1 - tmp2 * tmp4 );
         coefficients[coefficients_index][1] = tmp8 * ( tmp4 + tmp7 );
         coefficients[coefficients_index][2] = tmp6 * ( -tmp5 - tmp7 );
         coefficients[coefficients_index][3] = -80 * tmp8;
         coefficients[coefficients_index][4] = 16 * tmp6;
         break;
      }
      case 5: {
         const real_t tmp0 = std::pow( delta, 4 );
         const real_t tmp1 = std::pow( delta, 2 );
         const real_t tmp2 = std::pow( theta, 3 );
         const real_t tmp3 = std::pow( theta, 4 );
         const real_t tmp4 = 48 * tmp1;
         const real_t tmp5 = std::pow( theta, 2 );
         const real_t tmp6 = 18 * tmp0;
         const real_t tmp7 = 1.0 / ( std::pow( delta, 6 ) - 32 * std::pow( theta, 6 ) + tmp3 * tmp4 - tmp5 * tmp6 );

         coefficients[coefficients_index][0] = tmp7 * ( -192 * std::pow( theta, 5 ) - 36 * theta * tmp0 + 192 * tmp1 * tmp2 );
         coefficients[coefficients_index][1] = tmp7 * ( -288 * tmp1 * tmp5 + 480 * tmp3 + tmp6 );
         coefficients[coefficients_index][2] = tmp7 * ( 192 * theta * tmp1 - 640 * tmp2 );
         coefficients[coefficients_index][3] = tmp7 * ( -tmp4 + 480 * tmp5 );
         coefficients[coefficients_index][4] = -192 * theta * tmp7;
         coefficients[coefficients_index][5] = 32 * tmp7;
         break;
      }
      case 6: {
         const real_t tmp0 = 7 * std::pow( delta, 6 );
         const real_t tmp1 = std::pow( theta, 6 );
         const real_t tmp2 = std::pow( theta, 4 );
         const real_t tmp3 = std::pow( delta, 2 );
         const real_t tmp4 = 560 * tmp3;
         const real_t tmp5 = std::pow( theta, 2 );
         const real_t tmp6 = std::pow( delta, 4 );
         const real_t tmp7 = 168 * tmp6;
         const real_t tmp8 = 112 * tmp3;
         const real_t tmp9 = 56 * tmp6;
         const real_t tmp10 =
             1.0 / ( -64 * std::pow( theta, 7 ) + std::pow( theta, 5 ) * tmp8 - std::pow( theta, 3 ) * tmp9 + theta * tmp0 );
         const real_t tmp11 = -1120 * tmp3 * tmp5;
         const real_t tmp12 = 1.0 / ( tmp0 - 64 * tmp1 + tmp2 * tmp8 - tmp5 * tmp9 );

         coefficients[coefficients_index][0] = tmp10 * ( tmp0 - 448 * tmp1 + tmp2 * tmp4 - tmp5 * tmp7 );
         coefficients[coefficients_index][1] = tmp12 * ( tmp11 + 1344 * tmp2 + tmp7 );
         coefficients[coefficients_index][2] = tmp10 * ( -tmp11 - 2240 * tmp2 - tmp9 );
         coefficients[coefficients_index][3] = tmp12 * ( -tmp4 + 2240 * tmp5 );
         coefficients[coefficients_index][4] = tmp10 * ( -1344 * tmp5 + tmp8 );
         coefficients[coefficients_index][5] = 448 * tmp12;
         coefficients[coefficients_index][6] = -64 * tmp10;
         break;
      }
      case 7: {
         const real_t tmp0  = std::pow( delta, 6 );
         const real_t tmp1  = std::pow( theta, 5 );
         const real_t tmp2  = std::pow( delta, 2 );
         const real_t tmp3  = 1536 * tmp2;
         const real_t tmp4  = std::pow( delta, 4 );
         const real_t tmp5  = std::pow( theta, 3 );
         const real_t tmp6  = std::pow( theta, 6 );
         const real_t tmp7  = 256 * tmp2;
         const real_t tmp8  = std::pow( theta, 4 );
         const real_t tmp9  = 160 * tmp4;
         const real_t tmp10 = std::pow( theta, 2 );
         const real_t tmp11 = 32 * tmp0;
         const real_t tmp12 =
             1.0 / ( std::pow( delta, 8 ) + 128 * std::pow( theta, 8 ) - tmp10 * tmp11 - tmp6 * tmp7 + tmp8 * tmp9 );

         coefficients[coefficients_index][0] =
             tmp12 * ( 1024 * std::pow( theta, 7 ) - 64 * theta * tmp0 - tmp1 * tmp3 + 640 * tmp4 * tmp5 );
         coefficients[coefficients_index][1] = tmp12 * ( -960 * tmp10 * tmp4 + tmp11 + 3840 * tmp2 * tmp8 - 3584 * tmp6 );
         coefficients[coefficients_index][2] = tmp12 * ( 640 * theta * tmp4 + 7168 * tmp1 - 5120 * tmp2 * tmp5 );
         coefficients[coefficients_index][3] = tmp12 * ( 3840 * tmp10 * tmp2 - 8960 * tmp8 - tmp9 );
         coefficients[coefficients_index][4] = tmp12 * ( -theta * tmp3 + 7168 * tmp5 );
         coefficients[coefficients_index][5] = tmp12 * ( -3584 * tmp10 + tmp7 );
         coefficients[coefficients_index][6] = 1024 * theta * tmp12;
         coefficients[coefficients_index][7] = -128 * tmp12;
         break;
      }
      case 8: {
         const real_t tmp0  = 9 * std::pow( delta, 8 );
         const real_t tmp1  = std::pow( theta, 8 );
         const real_t tmp2  = std::pow( theta, 6 );
         const real_t tmp3  = std::pow( delta, 2 );
         const real_t tmp4  = 4032 * tmp3;
         const real_t tmp5  = std::pow( theta, 4 );
         const real_t tmp6  = std::pow( delta, 4 );
         const real_t tmp7  = 2160 * tmp6;
         const real_t tmp8  = std::pow( theta, 2 );
         const real_t tmp9  = std::pow( delta, 6 );
         const real_t tmp10 = 360 * tmp9;
         const real_t tmp11 = 576 * tmp3;
         const real_t tmp12 = 432 * tmp6;
         const real_t tmp13 = 120 * tmp9;
         const real_t tmp14 = 1.0 / ( 256 * std::pow( theta, 9 ) - std::pow( theta, 7 ) * tmp11 + std::pow( theta, 5 ) * tmp12 -
                                      std::pow( theta, 3 ) * tmp13 + theta * tmp0 );
         const real_t tmp15 = tmp3 * tmp5;
         const real_t tmp16 = -4320 * tmp6 * tmp8;
         const real_t tmp17 = 1.0 / ( tmp0 + 256 * tmp1 - tmp11 * tmp2 + tmp12 * tmp5 - tmp13 * tmp8 );
         const real_t tmp18 = 32256 * tmp5;

         coefficients[coefficients_index][0] = tmp14 * ( tmp0 + 2304 * tmp1 - tmp10 * tmp8 - tmp2 * tmp4 + tmp5 * tmp7 );
         coefficients[coefficients_index][1] = tmp17 * ( tmp10 + 12096 * tmp15 + tmp16 - 9216 * tmp2 );
         coefficients[coefficients_index][2] = tmp14 * ( -tmp13 - 20160 * tmp15 - tmp16 + 21504 * tmp2 );
         coefficients[coefficients_index][3] = tmp17 * ( -tmp18 + 20160 * tmp3 * tmp8 - tmp7 );
         coefficients[coefficients_index][4] = tmp14 * ( tmp12 + tmp18 - 12096 * tmp3 * tmp8 );
         coefficients[coefficients_index][5] = tmp17 * ( tmp4 - 21504 * tmp8 );
         coefficients[coefficients_index][6] = tmp14 * ( -tmp11 + 9216 * tmp8 );
         coefficients[coefficients_index][7] = -2304 * tmp17;
         coefficients[coefficients_index][8] = 256 * tmp14;
         break;
      }
      case 9: {
         const real_t tmp0  = std::pow( delta, 8 );
         const real_t tmp1  = std::pow( delta, 2 );
         const real_t tmp2  = std::pow( theta, 7 );
         const real_t tmp3  = std::pow( theta, 5 );
         const real_t tmp4  = std::pow( delta, 4 );
         const real_t tmp5  = 6720 * tmp4;
         const real_t tmp6  = std::pow( delta, 6 );
         const real_t tmp7  = std::pow( theta, 3 );
         const real_t tmp8  = std::pow( theta, 8 );
         const real_t tmp9  = 1280 * tmp1;
         const real_t tmp10 = std::pow( theta, 6 );
         const real_t tmp11 = 1120 * tmp4;
         const real_t tmp12 = std::pow( theta, 4 );
         const real_t tmp13 = 400 * tmp6;
         const real_t tmp14 = std::pow( theta, 2 );
         const real_t tmp15 = 50 * tmp0;
         const real_t tmp16 = 1.0 / ( std::pow( delta, 10 ) - 512 * std::pow( theta, 10 ) - tmp10 * tmp11 + tmp12 * tmp13 -
                                      tmp14 * tmp15 + tmp8 * tmp9 );
         const real_t tmp17 = 35840 * tmp1;

         coefficients[coefficients_index][0] = tmp16 * ( -5120 * std::pow( theta, 9 ) - 100 * theta * tmp0 + 10240 * tmp1 * tmp2 -
                                                         tmp3 * tmp5 + 1600 * tmp6 * tmp7 );
         coefficients[coefficients_index][1] =
             tmp16 * ( -tmp10 * tmp17 + 16800 * tmp12 * tmp4 - 2400 * tmp14 * tmp6 + tmp15 + 23040 * tmp8 );
         coefficients[coefficients_index][2] =
             tmp16 * ( 1600 * theta * tmp6 + 71680 * tmp1 * tmp3 - 61440 * tmp2 - 22400 * tmp4 * tmp7 );
         coefficients[coefficients_index][3] = tmp16 * ( -89600 * tmp1 * tmp12 + 107520 * tmp10 - tmp13 + 16800 * tmp14 * tmp4 );
         coefficients[coefficients_index][4] = tmp16 * ( -theta * tmp5 + 71680 * tmp1 * tmp7 - 129024 * tmp3 );
         coefficients[coefficients_index][5] = tmp16 * ( tmp11 + 107520 * tmp12 - tmp14 * tmp17 );
         coefficients[coefficients_index][6] = tmp16 * ( 10240 * theta * tmp1 - 61440 * tmp7 );
         coefficients[coefficients_index][7] = tmp16 * ( 23040 * tmp14 - tmp9 );
         coefficients[coefficients_index][8] = -5120 * theta * tmp16;
         coefficients[coefficients_index][9] = 512 * tmp16;
         break;
      }
      default: {
         WALBERLA_ABORT( "Chebyshev smoother not implemented for order " << order );
      }
      }
   }

 protected:
   std::vector< std::vector< real_t > >             coefficients;
   std::shared_ptr< FunctionType >                  tmp1_;
   std::shared_ptr< FunctionType >                  tmp2_;
   DoFType                                          flag_;
   std::shared_ptr< hyteg::Solver< OperatorType > > chebyPrec_;
   std::shared_ptr< ProjectionOperatorType >        projection_;
   DoFType                                          projectionFlag_;

   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   bool lowMemoryMode_;
};

/// Namespace for utility functions of the Chebyshev-Smoother
namespace chebyshev {

/// Wraps an operator A, such that the apply method evaluates D^{-1}A, where D is the diagonal part of A.
/// This can be used, e.g. to determine the spectral radius of D^{-1}A with the spectral tools.
template < typename WrappedOperatorType >
class InvDiagOperatorWrapper : public Operator< typename WrappedOperatorType::srcType, typename WrappedOperatorType::dstType >
{
 public:
   using SrcFunctionType = typename WrappedOperatorType::srcType;
   using DstFunctionType = typename WrappedOperatorType::dstType;

   explicit InvDiagOperatorWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                    const uint_t&                              minLevel,
                                    const uint_t&                              maxLevel,
                                    const WrappedOperatorType&                 wrappedOperator )
   : Operator< SrcFunctionType, DstFunctionType >( storage, minLevel, maxLevel )
   , wrappedOperator_( wrappedOperator ){};

   /// Evaluates D^{-1} A.
   inline void apply( const SrcFunctionType& src,
                      const DstFunctionType& dst,
                      size_t                 level,
                      DoFType                flag,
                      UpdateType             updateType = Replace ) const
   {
      auto A_with_inv_diag       = dynamic_cast< const OperatorWithInverseDiagonal< SrcFunctionType >* >( &wrappedOperator_ );
      auto inverseDiagonalValues = A_with_inv_diag->getInverseDiagonalValues();

      WALBERLA_DEBUG_SECTION()
      {
         WALBERLA_ASSERT_NOT_NULLPTR( inverseDiagonalValues, "diagonal not initialized" );

         const real_t localNormSqr = inverseDiagonalValues->dotLocal( *inverseDiagonalValues, level, flag );
         WALBERLA_UNUSED( localNormSqr );
         WALBERLA_ASSERT_GREATER( localNormSqr, 0.0, "diagonal not set" );
      }

      wrappedOperator_.apply( src, dst, level, flag, updateType );
      dst.multElementwise( { *inverseDiagonalValues, dst }, level, flag );
   }

 private:
   const WrappedOperatorType& wrappedOperator_;
};

/// Wraps an operator A, such that the apply method evaluates D^{-1}A, where D is a given preconditioner.
/// This can be used, e.g. to determine the spectral radius of D^{-1}A with the spectral tools.
template < typename WrappedOperatorType >
class PreconditionerOperatorWrapper
: public Operator< typename WrappedOperatorType::srcType, typename WrappedOperatorType::dstType >
{
 public:
   using SrcFunctionType = typename WrappedOperatorType::srcType;
   using DstFunctionType = typename WrappedOperatorType::dstType;

   /// The createAllInnerBC() is used since estimateSpectralRadiusWithPowerIteration uses only DoFType All.
   explicit PreconditionerOperatorWrapper( const std::shared_ptr< PrimitiveStorage >&              storage,
                                           const uint_t&                                           minLevel,
                                           const uint_t&                                           maxLevel,
                                           const WrappedOperatorType&                              wrappedOperator,
                                           const std::shared_ptr< Solver< WrappedOperatorType > >& preconditioner )
   : Operator< SrcFunctionType, DstFunctionType >( storage, minLevel, maxLevel )
   , wrappedOperator_( wrappedOperator )
   , preconditioner_( preconditioner )
   , tmp_( "PreconditionerOperatorWrapper tmp", storage, minLevel, maxLevel, hyteg::BoundaryCondition::createAllInnerBC() )
   {}

   /// Evaluates D^{-1} A.
   inline void apply( const SrcFunctionType& src,
                      const DstFunctionType& dst,
                      size_t                 level,
                      DoFType                flag,
                      UpdateType             updateType = Replace ) const
   {
      wrappedOperator_.apply( src, dst, level, flag, updateType );

      tmp_.interpolate( real_c( 0 ), level, flag );
      preconditioner_->solve( wrappedOperator_, tmp_, dst, level );

      dst.assign( { real_c( 1 ) }, { tmp_ }, level, flag );
   }

 private:
   const WrappedOperatorType&                       wrappedOperator_;
   std::shared_ptr< Solver< WrappedOperatorType > > preconditioner_;
   mutable DstFunctionType                          tmp_;
};

/// Estimates the spectral radius of the Operator D^{-1}A for the Chebyshev iteration.
///
/// \tparam OperatorType
/// \param A Our operator.
/// \param level The level on which we should estimate the eigenvalue.
/// \param maxIter How many iterations, should we spend?
/// \param storage Our primitive storage.
/// \param x An initial guess of the eigenvector, which MUST NOT be in the kernel of A. Contains after the iteration an estimate for the eigenvector.
/// \return An estimate of the spectral radius.
template < typename OperatorType >
inline real_t estimateRadius( const OperatorType&                        A,
                              const uint_t&                              level,
                              const uint_t&                              maxIter,
                              const std::shared_ptr< PrimitiveStorage >& storage,
                              typename OperatorType::srcType&            x,
                              typename OperatorType::srcType&            tmp )
{
   InvDiagOperatorWrapper< OperatorType > invDiagA( storage, level, level, A );
   return hyteg::estimateSpectralRadiusWithPowerIteration( invDiagA, x, tmp, maxIter, storage, level );
}

/// Estimates the spectral radius of the Operator D^{-1}A for the Chebyshev iteration.
///
/// \tparam OperatorType
/// \param A Our operator.
/// \param level The level on which we should estimate the eigenvalue.
/// \param maxIter How many iterations, should we spend?
/// \param storage Our primitive storage.
/// \param x An initial guess of the eigenvector, which MUST NOT be in the kernel of A. Contains after the iteration an estimate for the eigenvector.
/// \param preconditioner A preconditioner D^{-1}. Set to nullptr if you want to use the inverse diagonal.
/// \return An estimate of the spectral radius.
template < typename OperatorType >
inline real_t estimateRadius( const OperatorType&                              A,
                              const uint_t&                                    level,
                              const uint_t&                                    maxIter,
                              const std::shared_ptr< PrimitiveStorage >&       storage,
                              typename OperatorType::srcType&                  x,
                              typename OperatorType::srcType&                  tmp,
                              const std::shared_ptr< Solver< OperatorType > >& preconditioner )
{
   if ( preconditioner != nullptr )
   {
      PreconditionerOperatorWrapper< OperatorType > preconditionedA( storage, level, level, A, preconditioner );
      return hyteg::estimateSpectralRadiusWithPowerIteration( preconditionedA, x, tmp, maxIter, storage, level );
   }
   else
   {
      WALBERLA_ABORT( "preconditioner undefined." );
      return real_c( 0 );
   }
}

// see Lemma 1 from "Matrix-free Monolithic Multigrid Methods for Stokes and Generalized Stokes Problems" ; Jodlbauer, Langer, Wick and Zulehner 2022
// https://doi.org/10.48550/arXiv.2205.15770
//
// beta is the biggest eigenvalue of D^-1 M
// (You can use the function estimateRadius in this header)
// 0 < alpha < beta
inline real_t determineScalingInternal( uint_t k, real_t alpha, real_t beta )
{
   real_t sbar = ( beta + alpha ) / ( beta - alpha );

   real_t chebyshevEval;

   switch ( k + 1 )
   {
   case 1: {
      chebyshevEval = sbar;
      break;
   }
   case 2: {
      chebyshevEval = real_c( 2.0 ) * std::pow( sbar, real_c( 2.0 ) ) - real_c( 1.0 );
      break;
   }
   case 3: {
      chebyshevEval = real_c( 4.0 ) * std::pow( sbar, real_c( 3.0 ) ) - real_c( 3.0 ) * sbar;
      break;
   }
   case 4: {
      chebyshevEval =
          real_c( 8.0 ) * std::pow( sbar, real_c( 4.0 ) ) - real_c( 8.0 ) * std::pow( sbar, real_c( 2.0 ) ) + real_c( 1.0 );
      break;
   }
   case 5: {
      chebyshevEval = real_c( 16.0 ) * std::pow( sbar, real_c( 5.0 ) ) - real_c( 20.0 ) * std::pow( sbar, real_c( 3.0 ) ) +
                      real_c( 5.0 ) * sbar;
      break;
   }
   default:
      WALBERLA_ABORT( "determineScaling not implemented for order " << k << " Chebyshev smoother." );
   }

   return chebyshevEval / ( real_c( 1.0 ) + chebyshevEval );
}

inline real_t
    determineScaling( uint_t k, real_t spectralRadius, real_t upperFactor = real_c( 1.2 ), real_t lowerFactor = real_c( 0.3 ) )
{
   return determineScalingInternal( k, lowerFactor * spectralRadius, upperFactor * spectralRadius );
}

} // namespace chebyshev

} // namespace hyteg