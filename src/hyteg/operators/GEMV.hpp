/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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
#include "core/config/Config.h"

#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/misc/SFINAE.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/types/types.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

// this type gets returned by the GEMV functions, signaling which operation was used
enum GEMVType
{
   GEMV,
   SCALED,
   MANUAL,
   NOSCALINGREQUIRED
};

// dst = alpha * A(src) + gamma * dst
template < class OperatorType >
GEMVType applyGEMV( const OperatorType&                       Op,
                    typename OperatorType::srcType::valueType alpha,
                    const typename OperatorType::srcType&     src,
                    typename OperatorType::dstType::valueType gamma,
                    const typename OperatorType::dstType&     dst,
                    const uint_t                              level,
                    const hyteg::DoFType                      flag )
{
   const bool gammaZero = std::fpclassify( gamma ) == FP_ZERO;
   const bool gammaOne  = std::fpclassify( gamma - real_c( 1.0 ) ) == FP_ZERO;
   const bool alphaOne  = std::fpclassify( alpha - real_c( 1.0 ) ) == FP_ZERO;

   if constexpr ( SFINAE::has_gemv< OperatorType >() )
   {
      Op.gemv( alpha, src, gamma, dst, level, flag );
      return GEMVType::GEMV;
   }
   else
   {
      if ( gammaZero )
      {
         if ( alphaOne )
         {
            // just apply the operator
            Op.apply( src, dst, level, flag, hyteg::UpdateType::Replace );
            return GEMVType::NOSCALINGREQUIRED;
         }
         else
         {
            // overwrite dst
            // scale the result afterwards
            if constexpr ( SFINAE::has_applyScaled< OperatorType >() )
            {
               Op.applyScaled( alpha, src, dst, level, flag, hyteg::UpdateType::Replace );
               return GEMVType::SCALED;
            }
            else
            {
               Op.apply( src, dst, level, flag, hyteg::UpdateType::Replace );
               dst.assign( { alpha }, { dst }, level, flag );
               return GEMVType::MANUAL;
            }
         }
      }
      else if ( gammaOne )
      {
         if ( alphaOne )
         {
            // apply the operator with Add
            Op.apply( src, dst, level, flag, hyteg::UpdateType::Add );
            return GEMVType::NOSCALINGREQUIRED;
         }
         else
         {
            if constexpr ( SFINAE::has_applyScaled< OperatorType >() )
            {
               // apply the scaled operator with Add
               Op.applyScaled( alpha, src, dst, level, flag, hyteg::UpdateType::Add );
               return GEMVType::SCALED;
            }
            else
            {
               // apply operator to tmp function
               // add scaled tmp to unscaled dst
               // the real_c(1) cannot be removed
               auto tmpFct = getTemporaryFunction< typename OperatorType::dstType >(
                   dst.getStorage(), dst.getMinLevel(), dst.getMaxLevel() );

               tmpFct->copyBoundaryConditionFromFunction( dst );

               Op.apply( src, *tmpFct, level, flag, hyteg::UpdateType::Replace );

               dst.assign( { real_c( 1 ), alpha }, { dst, *tmpFct }, level, flag );
               return GEMVType::MANUAL;
            }
         }
      }
      else
      {
         if ( alphaOne )
         {
            // scale dst
            // apply the operator afterwards with Add
            dst.assign( { gamma }, { dst }, level, flag );
            Op.apply( src, dst, level, flag, hyteg::UpdateType::Add );
            return GEMVType::NOSCALINGREQUIRED;
         }
         else
         {
            if constexpr ( SFINAE::has_applyScaled< OperatorType >() )
            {
               // scale dst
               // apply the scaled operator afterwards with Add
               dst.assign( { gamma }, { dst }, level, flag );
               Op.applyScaled( alpha, src, dst, level, flag, hyteg::UpdateType::Add );
               return GEMVType::SCALED;
            }
            else
            {
               // apply operator to tmp function
               // add scaled tmp to scaled dst
               auto tmpFct = getTemporaryFunction< typename OperatorType::dstType >(
                   dst.getStorage(), dst.getMinLevel(), dst.getMaxLevel() );

               tmpFct->copyBoundaryConditionFromFunction( dst );

               Op.apply( src, *tmpFct, level, flag, hyteg::UpdateType::Replace );

               dst.assign( { gamma, alpha }, { dst, *tmpFct }, level, flag );
               return GEMVType::MANUAL;
            }
         }
      }
   }
}

template < class OperatorType >
GEMVType applyToMatrixScaled( const OperatorType&                                                          Op,
                              const typename OperatorType::srcType::valueType&                             alpha,
                              const std::shared_ptr< hyteg::SparseMatrixProxy >&                           mat,
                              const typename OperatorType::srcType::template FunctionType< hyteg::idx_t >& src,
                              const typename OperatorType::dstType::template FunctionType< hyteg::idx_t >& dst,
                              uint_t                                                                       level,
                              hyteg::DoFType                                                               flag )
{
   const bool scalingZero = std::fpclassify( alpha ) == FP_ZERO;
   const bool scalingOne  = std::fpclassify( alpha - real_c( 1.0 ) ) == FP_ZERO;

   if ( !scalingZero )
   {
      if ( scalingOne )
      {
         // just call toMatrix
         Op.toMatrix( mat, src, dst, level, flag );
         return GEMVType::NOSCALINGREQUIRED;
      }
      else
      {
         if constexpr ( SFINAE::has_toMatrixScaled< OperatorType >() )
         {
            Op.toMatrixScaled( alpha, mat, src, dst, level, flag );
            return GEMVType::SCALED;
         }
         else
         {
            // toMatrix to temporary copyMat matrix
            // add scaled copyMat to mat
            std::shared_ptr< hyteg::SparseMatrixProxy > copyMat = mat->createEmptyCopy();

            Op.toMatrix( copyMat, src, dst, level, flag );

            mat->createFromMatrixLinComb( { real_c( 1.0 ), alpha }, { mat, copyMat } );
            return GEMVType::MANUAL;
         }
      }
   }

   return GEMVType::NOSCALINGREQUIRED;
}

template < class OperatorType >
GEMVType applyComputeInverseDiagonalOperatorValuesScaled( OperatorType&                                    Op,
                                                          const typename OperatorType::srcType::valueType& alpha )
{
   const bool scalingOne = std::fpclassify( alpha - real_c( 1.0 ) ) == FP_ZERO;

   if ( scalingOne )
   {
      // just call computeInverseDiagonalOperatorValues
      Op.computeInverseDiagonalOperatorValues();
      return GEMVType::NOSCALINGREQUIRED;
   }
   else
   {
      if constexpr ( SFINAE::has_computeInverseDiagonalOperatorValuesScaled< OperatorType >() )
      {
         Op.computeInverseDiagonalOperatorValuesScaled( alpha );
         return GEMVType::SCALED;
      }
      else
      {
         Op.computeInverseDiagonalOperatorValues();

         std::shared_ptr< typename OperatorType::srcType > diag = Op.getInverseDiagonalValues();

         for ( uint_t level = Op.getMinLevel(); level <= Op.getMaxLevel(); level++ )
         {
            diag->assign( { real_c( 1 ) / alpha }, { *diag }, level, hyteg::DoFType::All );
         }

         return GEMVType::MANUAL;
      }
   }
}

/// @brief Assumes that you have precalculated the scaled inverse diagonals, e.g., via applyComputeInverseDiagonalOperatorValuesScaled.
template < class OperatorType >
GEMVType applySmoothJacScaled( OperatorType&                                    Op,
                               const typename OperatorType::srcType::valueType& alpha,
                               const typename OperatorType::srcType&            dst,
                               const typename OperatorType::srcType&            rhs,
                               const typename OperatorType::srcType&            src,
                               typename OperatorType::srcType::valueType        omega,
                               uint_t                                           level,
                               DoFType                                          flag )
{
   const bool scalingOne = std::fpclassify( alpha - real_c( 1.0 ) ) == FP_ZERO;

   if ( scalingOne )
   {
      // just call smooth_jac
      Op.smooth_jac( dst, rhs, src, omega, level, flag );
      return GEMVType::NOSCALINGREQUIRED;
   }
   else
   {
      if constexpr ( SFINAE::has_smooth_jac_scaled< OperatorType >() )
      {
         Op.smooth_jac_scaled( alpha, dst, rhs, src, omega, level, flag );
         return GEMVType::SCALED;
      }
      else
      {
         if constexpr ( SFINAE::has_getInverseDiagonalValues< OperatorType >() )
         {
            Op.startTiming( "smooth_jac" );

            // compute the current residual
            Op.apply( src, dst, level, flag, hyteg::UpdateType::Replace );
            dst.assign( { real_c( 1 ), real_c( -alpha ) }, { rhs, dst }, level, flag );

            // perform Jacobi update step
            dst.multElementwise( { *Op.getInverseDiagonalValues(), dst }, level, flag );
            dst.assign( { real_c( 1 ), omega }, { src, dst }, level, flag );

            Op.stopTiming( "smooth_jac" );

            return GEMVType::MANUAL;
         }
         else
         {
            WALBERLA_ABORT( "Cannot perform applySmoothJacScaled on an operator without getInverseDiagonalValues." )
            return GEMVType::MANUAL;
         }
      }
   }
}
} // namespace hyteg