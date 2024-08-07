/*
 * Copyright (c) 2020-2022 Andreas Wagner, Marcus Mohr
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

#include "hyteg/numerictools/SpectrumEstimation.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Smoothables.hpp"

#include "Solver.hpp"

namespace hyteg {

template < typename OperatorType >
class ChebyshevSmoother : public Solver< OperatorType >
{
 public:
   using FunctionType = typename OperatorType::srcType;
   using ValueType    = typename FunctionTrait< FunctionType >::ValueType;

   ChebyshevSmoother( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : coefficients{}
   , tmp1_( "cheb_tmp1", storage, minLevel, maxLevel )
   , tmp2_( "cheb_tmp2", storage, minLevel, maxLevel )
   , flag_( Inner | NeumannBoundary | FreeslipBoundary )
   , timingTree_( storage->getTimingTree() )
   {}

   ChebyshevSmoother( const FunctionType& tmp1, const FunctionType& tmp2 )
   : coefficients{}
   , tmp1_( tmp1 )
   , tmp2_( tmp2 )
   , flag_( Inner | NeumannBoundary | FreeslipBoundary )
   , timingTree_( tmp1.getStorage()->getTimingTree() )
   {}

   /// Executes an iteration step of the smoother.
   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
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

      tmp1_.copyBoundaryConditionFromFunction( x );
      tmp2_.copyBoundaryConditionFromFunction( x );

      WALBERLA_DEBUG_SECTION()
      {
         WALBERLA_ASSERT( coefficients.size() > 0, "coefficients must be setup" );
         WALBERLA_ASSERT_NOT_NULLPTR( inverseDiagonalValues, "diagonal not initialized" );

         const real_t localNormSqr = inverseDiagonalValues->dotLocal( *inverseDiagonalValues, level, flag_ );
         WALBERLA_UNUSED( localNormSqr );
         WALBERLA_ASSERT_GREATER( localNormSqr, 0.0, "diagonal not set" );
      }
      timingTree_->start( "Chebyshev Smoother" );

      // tmp1_ := Ax
      A.apply( x, tmp1_, level, flag_ );
      // tmp1_ := b-Ax
      tmp1_.assign(
          { walberla::numeric_cast< ValueType >( 1. ), walberla::numeric_cast< ValueType >( -1. ) }, { b, tmp1_ }, level, flag_ );
      // tmp1_ := D^{-1} (b-Ax)
      tmp1_.multElementwise( { *inverseDiagonalValues, tmp1_ }, level, flag_ );
      // x := x + omega_0 D^{-1} (b-Ax)
      x.assign( { walberla::numeric_cast< ValueType >( 1. ), walberla::numeric_cast< ValueType >( coefficients[0] ) },
                { x, tmp1_ },
                level,
                flag_ );

      // Loop preconditions before the kth iteration:
      // 1.)  x := x + sum_{i<k} omega_i (D^{-1}A)^{i} D^{-1} (b-Ax)
      // 2.) tmp1_ := (D^{-1}A)^{k-1} D^{-1} (b-Ax)
      for ( uint_t k = 1; k < coefficients.size(); k += 1 )
      {
         // tmp2_ := A (D^{-1}A)^{k-1} D^{-1} (b-Ax)
         A.apply( tmp1_, tmp2_, level, flag_ );
         // tmp1_ := (D^{-1}A)^{k} D^{-1} (b-Ax)
         tmp1_.multElementwise( { *inverseDiagonalValues, tmp2_ }, level, flag_ );

         // x := x + sum_{i<k+1} omega_i (D^{-1}A)^{i} D^{-1} (b-Ax)
         x.assign( { walberla::numeric_cast< ValueType >( 1 ), walberla::numeric_cast< ValueType >( coefficients[k] ) },
                   { x, tmp1_ },
                   level,
                   flag_ );
      }

      timingTree_->stop( "Chebyshev Smoother" );
   }

   /// Calculates the coefficients for our Chebyshev-Smoother.
   /// Has to be called prior to the first usage of the solver.
   ///
   /// This function is based on the mfem implementation at
   ///     https://github.com/mfem/mfem/blob/7cbb0d484863bf661e88378eac4ab0247f10545c/linalg/solvers.cpp#L185
   /// which references the article
   ///     Parallel multigrid smoothing: polynomial versus Gauss-Seidel by Adams et al.
   ///
   /// \param order The order of our polynomial smoother. Only 0 to 4 are supported.
   /// \param spectralRadius An estimate for our spectral radius.
   void setupCoefficients( const uint_t order,
                           const real_t spectralRadius,
                           real_t       upperFactor = real_c( 1.2 ),
                           real_t       lowerFactor = real_c( 0.3 ) )
   {
      const real_t upperBound = upperFactor * spectralRadius;
      const real_t lowerBound = lowerFactor * spectralRadius;

      setupCoefficientsInternal( order, lowerBound, upperBound );
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
   void setupCoefficientsInternal( const uint_t order, const real_t lowerBound, const real_t upperBound )
   {
      WALBERLA_ASSERT_GREATER( order, 0, "Order cannot be 0." );
      WALBERLA_ASSERT_LESS( order, 6, "Chebyshev-Smoother does not support polynomial orders larger than 5." );

      const double theta = 0.5 * ( upperBound + lowerBound );
      const double delta = 0.5 * ( upperBound - lowerBound );

      coefficients.resize( order, 0 );

      switch ( order - 1 )
      {
      case 0: {
         coefficients[0] = 1.0 / theta;
         break;
      }
      case 1: {
         const double tmp_0 = 1.0 / ( pow( delta, 2 ) - 2 * pow( theta, 2 ) );
         coefficients[0]    = -4 * theta * tmp_0;
         coefficients[1]    = 2 * tmp_0;
         break;
      }
      case 2: {
         const double tmp_0 = 3 * pow( delta, 2 );
         const double tmp_1 = pow( theta, 2 );
         const double tmp_2 = 1.0 / ( -4 * pow( theta, 3 ) + theta * tmp_0 );
         coefficients[0]    = tmp_2 * ( tmp_0 - 12 * tmp_1 );
         coefficients[1]    = 12 / ( tmp_0 - 4 * tmp_1 );
         coefficients[2]    = -4 * tmp_2;
         break;
      }
      case 3: {
         const double tmp_0 = pow( delta, 2 );
         const double tmp_1 = pow( theta, 2 );
         const double tmp_2 = 8 * tmp_0;
         const double tmp_3 = 1.0 / ( pow( delta, 4 ) + 8 * pow( theta, 4 ) - tmp_1 * tmp_2 );
         coefficients[0]    = tmp_3 * ( 32 * pow( theta, 3 ) - 16 * theta * tmp_0 );
         coefficients[1]    = tmp_3 * ( -48 * tmp_1 + tmp_2 );
         coefficients[2]    = 32 * theta * tmp_3;
         coefficients[3]    = -8 * tmp_3;
         break;
      }
      case 4: {
         const double tmp_0 = 5 * pow( delta, 4 );
         const double tmp_1 = pow( theta, 4 );
         const double tmp_2 = pow( theta, 2 );
         const double tmp_3 = pow( delta, 2 );
         const double tmp_4 = 60 * tmp_3;
         const double tmp_5 = 20 * tmp_3;
         const double tmp_6 = 1.0 / ( 16 * pow( theta, 5 ) - pow( theta, 3 ) * tmp_5 + theta * tmp_0 );
         const double tmp_7 = 160 * tmp_2;
         const double tmp_8 = 1.0 / ( tmp_0 + 16 * tmp_1 - tmp_2 * tmp_5 );
         coefficients[0]    = tmp_6 * ( tmp_0 + 80 * tmp_1 - tmp_2 * tmp_4 );
         coefficients[1]    = tmp_8 * ( tmp_4 - tmp_7 );
         coefficients[2]    = tmp_6 * ( -tmp_5 + tmp_7 );
         coefficients[3]    = -80 * tmp_8;
         coefficients[4]    = 16 * tmp_6;
         break;
      }
      default:
         WALBERLA_ABORT( "Chebyshev smoother not implemented for order " << order );
      }
   }

 protected:
   std::vector< real_t >                     coefficients;
   FunctionType                              tmp1_;
   FunctionType                              tmp2_;
   DoFType                                   flag_;
   std::shared_ptr< walberla::WcTimingTree > timingTree_;
};

/***************************************************************************
NOTE: This is similar to the Chebyshev smoother except that projection is
      applied to set normal components to zero at the FreeslipBoundary at
      every step of working
***************************************************************************/
template < typename OperatorType >
class ChebyshevSmootherWithFreeSlipProjection : public ChebyshevSmoother< OperatorType >
{
 public:
   using FunctionType = typename OperatorType::srcType;
   using ChebyshevSmoother< OperatorType >::tmp1_;
   using ChebyshevSmoother< OperatorType >::tmp2_;
   using ChebyshevSmoother< OperatorType >::flag_;
   using ChebyshevSmoother< OperatorType >::coefficients;

   ChebyshevSmootherWithFreeSlipProjection( const std::shared_ptr< PrimitiveStorage >& storage,
                                            size_t                                     minLevel,
                                            size_t                                     maxLevel,
                                            std::shared_ptr< P2ProjectNormalOperator > projection )
   : ChebyshevSmoother< OperatorType >( storage, minLevel, maxLevel )
   , projection_( projection )
   {}

   /// Executes an iteration step of the smoother.
   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      std::shared_ptr< typename OperatorType::srcType > inverseDiagonalValues = nullptr;

      if ( const auto* A_with_inv_diag =
               dynamic_cast< const OperatorWithInverseDiagonal< typename OperatorType::srcType >* >( &( A.viscousOperator ) ) )
      {
         inverseDiagonalValues = A_with_inv_diag->getInverseDiagonalValues();
      }
      else
      {
         throw std::runtime_error( "The Chebyshev-Smoother requires the OperatorWithInverseDiagonal interface." );
      }

      tmp1_.copyBoundaryConditionFromFunction( x );
      tmp2_.copyBoundaryConditionFromFunction( x );

      WALBERLA_DEBUG_SECTION()
      {
         WALBERLA_ASSERT( coefficients.size() > 0, "coefficients must be setup" );
         WALBERLA_ASSERT_NOT_NULLPTR( inverseDiagonalValues, "diagonal not initialized" );

         const real_t localNormSqr = inverseDiagonalValues->dotLocal( *inverseDiagonalValues, level, flag_ );
         WALBERLA_UNUSED( localNormSqr );
         WALBERLA_ASSERT_GREATER( localNormSqr, 0.0, "diagonal not set" );
      }

      // tmp1_ := Ax
      A.apply( x, tmp1_, level, flag_ );
      // tmp1_ := b-Ax
      tmp1_.assign( { real_t( 1. ), real_t( -1. ) }, { b, tmp1_ }, level, flag_ );
      // tmp1_ := D^{-1} (b-Ax)
      tmp1_.multElementwise( { *inverseDiagonalValues, tmp1_ }, level, flag_ );
      projection_->project( tmp1_, level, FreeslipBoundary );
      // x := x + omega_0 D^{-1} (b-Ax)
      x.assign( { 1., coefficients[0] }, { x, tmp1_ }, level, flag_ );

      // Loop preconditions before the kth iteration:
      // 1.)  x := x + sum_{i<k} omega_i (D^{-1}A)^{i} D^{-1} (b-Ax)
      // 2.) tmp1_ := (D^{-1}A)^{k-1} D^{-1} (b-Ax)
      for ( uint_t k = 1; k < coefficients.size(); k += 1 )
      {
         // tmp2_ := A (D^{-1}A)^{k-1} D^{-1} (b-Ax)
         A.apply( tmp1_, tmp2_, level, flag_ );
         // tmp1_ := (D^{-1}A)^{k} D^{-1} (b-Ax)
         tmp1_.multElementwise( { *inverseDiagonalValues, tmp2_ }, level, flag_ );
         projection_->project( tmp1_, level, FreeslipBoundary );
         // x := x + sum_{i<k+1} omega_i (D^{-1}A)^{i} D^{-1} (b-Ax)
         x.assign( { 1, coefficients[k] }, { x, tmp1_ }, level, flag_ );
      }

      projection_->project( tmp1_, level, FreeslipBoundary );
   }

 private:
   std::shared_ptr< P2ProjectNormalOperator > projection_;
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
