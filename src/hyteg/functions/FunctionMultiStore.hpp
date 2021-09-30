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

#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/functions/FunctionWrapper.hpp"

namespace hyteg {

using walberla::uint_t;

/// Class for storing multiple functions from the same family but with (potentially) different value types
///
/// This class allows to store multiple functions from the same family, which can but need not differ in
/// their value types. The class is currently intended as a replacement for the std::vector< "func_t"< real_t >
/// containes in VTKOuput and thus tries to behave as much as possbile like an std::vector. Currently
/// supported value types are
/// - double
/// - float
/// - int32_t
/// - int64_t
/// Extension to other datatypes is straightforward.
template < template < class > class func_t >
class FunctionMultiStore
{
 public:
   template < typename value_t >
   inline void push_back( const func_t< value_t >& function )
   {
      if constexpr ( std::is_same< value_t, double >::value )
      {
         r64Funcs.push_back( function );
      }
      else if constexpr ( std::is_same< value_t, float >::value )
      {
         r32Funcs.push_back( function );
      }
      else if constexpr ( std::is_same< value_t, int32_t >::value )
      {
         i32Funcs.push_back( function );
      }
      else if constexpr ( std::is_same< value_t, int64_t >::value )
      {
         i64Funcs.push_back( function );
      }
      else if constexpr ( std::is_same< value_t, long long >::value )
      {
         i64Funcs.push_back( function );
      }
      else
      {
         WALBERLA_ABORT( "FunctionMultiStore::add() detected unsupported datatype '" << typeid( function ).name() << "'" );
      };
   };

   uint_t size() const { return r64Funcs.size() + r32Funcs.size() + i32Funcs.size() + i64Funcs.size(); };

   template < typename value_t >
   const std::vector< func_t< value_t > >& getFunctions() const
   {
      if constexpr ( std::is_same< value_t, double >::value )
      {
         return r64Funcs;
      }
      else if constexpr ( std::is_same< value_t, float >::value )
      {
         return r32Funcs;
      }
      else if constexpr ( std::is_same< value_t, int >::value )
      {
         return i32Funcs;
      }
      else if constexpr ( std::is_same< value_t, long >::value )
      {
         return i64Funcs;
      }
      else if constexpr ( std::is_same< value_t, long long >::value )
      {
         return i64Funcs;
      }
      else
      {
         WALBERLA_ABORT( "FunctionMultiStore::getFunctions() detected unsupported datatype '" << typeid( value_t ).name() );
      }
   }

 private:
   std::vector< func_t< double > >  r64Funcs;
   std::vector< func_t< float > >   r32Funcs;
   std::vector< func_t< int32_t > > i32Funcs;
   std::vector< func_t< int64_t > > i64Funcs;
};

} // namespace hyteg
