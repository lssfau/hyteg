/*
 * Copyright (c) 2021 Marcus Mohr.
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

#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/operators/GenericOperator.hpp"

namespace hyteg {

using walberla::uint_t;

template < typename oper_t >
class OperatorWrapper final : public GenericOperator< typename FunctionTrait< typename oper_t::srcType >::ValueType >
{
 public:
   typedef typename FunctionTrait< typename oper_t::srcType >::ValueType value_t;

   /// No need for this one, if we do not implement a setter method for wrappedFunc_;
   OperatorWrapper() = delete;

   /// Constructor that constructs the operator which the class wraps itself around
   OperatorWrapper( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   {
      wrappedOper_ = std::make_unique< oper_t >( storage, minLevel, maxLevel );
   };

   ~OperatorWrapper() {}

   /// provide access to wrapped operator
   /// @{
   oper_t& unwrap() { return *wrappedOper_; }

   const oper_t& unwrap() const { return *wrappedOper_; }
   /// @}

   void apply( const GenericFunction< value_t >& src,
               const GenericFunction< value_t >& dst,
               size_t                            level,
               DoFType                           flag,
               UpdateType                        updateType = Replace ) const
   {
      wrappedOper_->apply( src.template unwrap< typename oper_t::srcType >(),
                           dst.template unwrap< typename oper_t::dstType >(),
                           level,
                           flag,
                           updateType );
   };

 private:
   std::unique_ptr< oper_t > wrappedOper_;
};

} // namespace hyteg
