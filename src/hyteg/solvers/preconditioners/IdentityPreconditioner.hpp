/*
 * Copyright (c) 2017-2024 Daniel Drzisga, Dominik Thoennes, Nils Kohl, Andreas Burkhart.
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
#include <vector>

#include "core/DataTypes.h"

#include "hyteg/operators/Operator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

template < class OperatorType >
class IdentityPreconditioner : public Solver< OperatorType >
{
 public:
   using FunctionType = typename OperatorType::srcType;
   using ValueType    = typename FunctionTrait< FunctionType >::ValueType;

   IdentityPreconditioner()
   : updateType_( Replace )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   {}

   void solve( const OperatorType&,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      x.assign( { walberla::numeric_cast< ValueType >( 1.0 ) }, { b }, level, flag_ );
   }

 private:
   UpdateType updateType_;
   DoFType    flag_;
};

template < class SrcType, class DstType >
class IdentityOperator : public Operator< SrcType, DstType >
{
 public:
   IdentityOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< SrcType, DstType >( storage, minLevel, maxLevel )
   {}

   void apply( const SrcType& src,
               const DstType& dst,
               const uint_t   level,
               const DoFType  flag,
               UpdateType     updateType = Replace ) const
   {
      if ( updateType == Replace )
      {
         dst.assign( { 1.0 }, { src }, level, flag );
      }
      else
      {
         dst.assign( { 1.0, 1.0 }, { dst, src }, level, flag );
      }
   }
};

} // namespace hyteg