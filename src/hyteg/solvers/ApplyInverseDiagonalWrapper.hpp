/*
 * Copyright (c) 2023-2025 Andreas Burkhart.
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

#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

template < typename OperatorType >
class ApplyDiagonalWrapper : public Operator< typename OperatorType::srcType, typename OperatorType::dstType >
{
 public:
   ApplyDiagonalWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                         const uint_t&                              minLevel,
                         const uint_t&                              maxLevel,
                         const OperatorType&                        wrappedOperator )
   : Operator< typename OperatorType::srcType, typename OperatorType::dstType >( storage, minLevel, maxLevel )
   , Operator_( wrappedOperator )
   , tmp_( "DiagTmp", storage, minLevel, maxLevel )
   {}

   /// Evaluates the multiplication with the diagonal.
   inline void apply( const typename OperatorType::srcType& src,
                      const typename OperatorType::dstType& dst,
                      size_t                                level,
                      DoFType                               flag,
                      UpdateType                            updateType = Replace ) const
   {
      auto DiagonalValues = Operator_.getDiagonalValues();

      tmp_.copyBoundaryConditionFromFunction( src );

      if ( updateType == Replace )
      {
         dst.multElementwise( { *DiagonalValues, src }, level, flag );
      }
      else
      {
         tmp_.multElementwise( { *DiagonalValues, src }, level, flag );
         dst.assign( { 1.0, 1.0 }, { dst, tmp_ }, level, flag );
      }
   }

 private:
   const OperatorType&                    Operator_;
   mutable typename OperatorType::srcType tmp_;
};

template < typename OperatorType >
class ApplyInverseDiagonalWrapper : public Operator< typename OperatorType::srcType, typename OperatorType::dstType >
{
 public:
   ApplyInverseDiagonalWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                const uint_t&                              minLevel,
                                const uint_t&                              maxLevel,
                                const OperatorType&                        wrappedOperator )
   : Operator< typename OperatorType::srcType, typename OperatorType::dstType >( storage, minLevel, maxLevel )
   , Operator_( wrappedOperator )
   , tmp_( "inverseDiagTmp", storage, minLevel, maxLevel )
   {}

   /// Evaluates the multiplication with the inverse diagonal.
   inline void apply( const typename OperatorType::srcType& src,
                      const typename OperatorType::dstType& dst,
                      size_t                                level,
                      DoFType                               flag,
                      UpdateType                            updateType = Replace ) const
   {
      auto A_with_inv_diag = dynamic_cast< const OperatorWithInverseDiagonal< typename OperatorType::srcType >* >( &Operator_ );
      auto inverseDiagonalValues = A_with_inv_diag->getInverseDiagonalValues();

      tmp_.copyBoundaryConditionFromFunction( src );

      if ( updateType == Replace )
      {
         dst.multElementwise( { *inverseDiagonalValues, src }, level, flag );
      }
      else
      {
         tmp_.multElementwise( { *inverseDiagonalValues, src }, level, flag );
         dst.assign( { 1.0, 1.0 }, { dst, tmp_ }, level, flag );
      }
   }

 private:
   const OperatorType&                    Operator_;
   mutable typename OperatorType::srcType tmp_;
};

template < typename OperatorType >
class ApplyLumpedInverseDiagonalWrapper : public Operator< typename OperatorType::srcType, typename OperatorType::dstType >
{
 public:
   ApplyLumpedInverseDiagonalWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                      const uint_t&                              minLevel,
                                      const uint_t&                              maxLevel,
                                      const OperatorType&                        wrappedOperator )
   : Operator< typename OperatorType::srcType, typename OperatorType::dstType >( storage, minLevel, maxLevel )
   , Operator_( wrappedOperator )
   , tmp_( "inverseDiagTmp", storage, minLevel, maxLevel )
   {}

   /// Evaluates the multiplication with the lumped inverse diagonal.
   inline void apply( const typename OperatorType::srcType& src,
                      const typename OperatorType::dstType& dst,
                      size_t                                level,
                      DoFType                               flag,
                      UpdateType                            updateType = Replace ) const
   {
      auto lumpedInverseDiagonalValues = Operator_.getLumpedInverseDiagonalValues();

      tmp_.copyBoundaryConditionFromFunction( src );

      if ( updateType == Replace )
      {
         dst.multElementwise( { *lumpedInverseDiagonalValues, src }, level, flag );
      }
      else
      {
         tmp_.multElementwise( { *lumpedInverseDiagonalValues, src }, level, flag );
         dst.assign( { 1.0, 1.0 }, { dst, tmp_ }, level, flag );
      }
   }

 private:
   const OperatorType&                    Operator_;
   mutable typename OperatorType::srcType tmp_;
};

template < typename OperatorType, typename FunctionType >
class ApplyFunctionMultiplicationWrapper : public Operator< typename OperatorType::srcType, typename OperatorType::dstType >
{
 public:
   ApplyFunctionMultiplicationWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                       const uint_t&                              minLevel,
                                       const uint_t&                              maxLevel,
                                       const std::shared_ptr< FunctionType >&     inverseFunction )
   : Operator< typename OperatorType::srcType, typename OperatorType::dstType >( storage, minLevel, maxLevel )
   , inverseFunction_( inverseFunction )
   {}

   /// Evaluates the multiplication with the inverse diagonal.
   inline void apply( const typename OperatorType::srcType& src,
                      const typename OperatorType::dstType& dst,
                      size_t                                level,
                      DoFType                               flag,
                      UpdateType                            updateType = Replace ) const
   {
      if ( updateType == Replace )
      {
         dst.multElementwise( { *inverseFunction_, src }, level, flag );
      }
      else
      {
         std::shared_ptr< typename OperatorType::srcType > tmp =
             getTemporaryFunction< typename OperatorType::srcType >( src.getStorage(), src.getMinLevel(), src.getMaxLevel() );

         tmp->copyBoundaryConditionFromFunction( src );
         tmp->multElementwise( { *inverseFunction_, src }, level, flag );
         dst.assign( { 1.0, 1.0 }, { dst, *tmp }, level, flag );
      }
   }

 private:
   std::shared_ptr< FunctionType > inverseFunction_;
};

} // namespace hyteg