/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/solvers/Solver.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

template < class OperatorType >
class IdentityPreconditioner : public Solver< OperatorType >
{
 public:
   IdentityPreconditioner()
   : updateType_( Replace )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   {}

   void solve( const OperatorType&,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      x.assign( { 1.0 }, { b }, level, flag_ );
   }

 private:
   UpdateType updateType_;
   DoFType    flag_;
};

template < typename functionType >
class IdentityOperator : public Operator< functionType, functionType >
{
 public:
   IdentityOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t level )
   : Operator< functionType, functionType >( storage, level, level )
   {}

   void apply( const functionType&    src,
               const functionType&    dst,
               const walberla::uint_t level,
               DoFType                flag,
               UpdateType             updateType )
   {
      dst.assign( { 1.0 }, { src }, level, flag );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const functionType&                         src,
                  const functionType&                         dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }
};
} // namespace hyteg