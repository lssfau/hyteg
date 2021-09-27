/*
 * Copyright (c) 2020-2021 Marcus Mohr.
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

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

using walberla::uint_t;

// forward declaration of child for use in mother class
template < typename func_t >
class OperatorWrapper;

template < typename value_t >
class GenericOperator
{
 public:
   virtual ~GenericOperator(){};

   template < typename oper_t >
   oper_t& unwrap()
   {
      auto realMe = static_cast< OperatorWrapper< oper_t >* >( this );
      return realMe->unwrap();
   };

   template < typename oper_t >
   const oper_t& unwrap() const
   {
      auto realMe = static_cast< const OperatorWrapper< oper_t >* >( this );
      return realMe->unwrap();
   };

   virtual void apply( const GenericFunction< value_t >& src,
                       const GenericFunction< value_t >& dst,
                       size_t                            level,
                       DoFType                           flag,
                       UpdateType                        updateType = Replace ) const = 0;

   virtual void smooth_gs( const GenericFunction< value_t >& src,
                           const GenericFunction< value_t >& dst,
                           size_t                            level,
                           DoFType                           flag ) const = 0;
};

} // namespace hyteg
