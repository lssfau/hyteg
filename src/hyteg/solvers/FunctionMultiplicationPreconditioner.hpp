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

#include "hyteg/operators/Operator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

template < typename OperatorType >
class FunctionMultiplicationPreconditioner : hyteg::Solver< OperatorType >
{
 public:
   FunctionMultiplicationPreconditioner( const std::shared_ptr< hyteg::PrimitiveStorage >&        storage,
                                         const uint_t&                                            minLevel,
                                         const uint_t&                                            maxLevel,
                                         const std::shared_ptr< typename OperatorType::dstType >& multFunction,
                                         hyteg::DoFType flag = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   : multFunction_( multFunction )
   , flag_( flag )
   {
      WALBERLA_UNUSED( storage );
      WALBERLA_UNUSED( minLevel );
      WALBERLA_UNUSED( maxLevel );
   }

   /// Evaluates the multiplication with the inverse diagonal.
   inline void apply( const typename OperatorType::dstType& src,
                      const typename OperatorType::srcType& dst,
                      uint_t                                level,
                      DoFType                               flag ) const
   {
      dst.multElementwise( { *multFunction_, src }, level, flag );
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      WALBERLA_UNUSED( A );
      apply( b, x, level, flag_ );
   };

 private:
   std::shared_ptr< typename OperatorType::dstType > multFunction_;
   hyteg::DoFType                                    flag_;
};

} // namespace hyteg