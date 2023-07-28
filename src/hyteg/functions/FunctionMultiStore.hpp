/*
 * Copyright (c) 2021-2023 Marcus Mohr.
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
/// their value types. The multistore ensures uniqueness of the functions, by looking at their names. An
/// attempt to add a function a second time will make the program abort. Currently supported value types are
///
/// - double
/// - float
/// - int32_t
/// - int64_t
///
/// Extension to other datatypes is straightforward.
template < template < class > class func_t >
class FunctionMultiStore
{
 public:
   /// Adds a function to the multistore
   template < typename value_t >
   inline void add( const func_t< value_t >& function )
   {
      bool        functionPresent{ false };
      bool        abort{ true }; // abort, if function is present?
      std::string functionName = function.getFunctionName();

      auto isFunctionPresent = [&functionName, &functionPresent, abort]( const auto& func ) {
         functionPresent = functionPresent || func.getFunctionName() == functionName;
         if ( functionPresent && abort )
         {
            WALBERLA_ABORT( "Attempt to add function '" << functionName << "' which is already present in FunctionMultiStore!" );
         }
      };

      if constexpr ( std::is_same< value_t, double >::value )
      {
         std::for_each( r64Funcs_.begin(), r64Funcs_.end(), isFunctionPresent );
         if ( !functionPresent )
            r64Funcs_.push_back( function );
      }
      else if constexpr ( std::is_same< value_t, float >::value )
      {
         std::for_each( r32Funcs_.begin(), r32Funcs_.end(), isFunctionPresent );
         if ( !functionPresent )
            r32Funcs_.push_back( function );
      }
      else if constexpr ( std::is_same< value_t, int32_t >::value )
      {
         std::for_each( i32Funcs_.begin(), i32Funcs_.end(), isFunctionPresent );
         if ( !functionPresent )
            i32Funcs_.push_back( function );
      }
      else if constexpr ( std::is_same< value_t, int64_t >::value )
      {
         std::for_each( i64Funcs_.begin(), i64Funcs_.end(), isFunctionPresent );
         if ( !functionPresent )
            i64Funcs_.push_back( function );
      }
      else if constexpr ( std::is_same< value_t, long long >::value )
      {
         std::for_each( i64Funcs_.begin(), i64Funcs_.end(), isFunctionPresent );
         if ( !functionPresent )
            i64Funcs_.push_back( function );
      }
      else
      {
         WALBERLA_ABORT( "FunctionMultiStore::add() detected unsupported datatype '" << typeid( function ).name() << "'" );
      };
   };

   /// Return the size of the FunctionMultiStore, i.e. the number of all functions stored (independent of their value type)
   uint_t size() const { return r64Funcs_.size() + r32Funcs_.size() + i32Funcs_.size() + i64Funcs_.size(); };

   /// Return a vector with all stored functions of a certain value type
   template < typename value_t >
   const std::vector< func_t< value_t > >& getFunctions() const
   {
      if constexpr ( std::is_same< value_t, double >::value )
      {
         return r64Funcs_;
      }
      else if constexpr ( std::is_same< value_t, float >::value )
      {
         return r32Funcs_;
      }
      else if constexpr ( std::is_same< value_t, int >::value )
      {
         return i32Funcs_;
      }
      else if constexpr ( std::is_same< value_t, long >::value )
      {
         return i64Funcs_;
      }
      else if constexpr ( std::is_same< value_t, long long >::value )
      {
         return i64Funcs_;
      }
      else
      {
         WALBERLA_ABORT( "FunctionMultiStore::getFunctions() detected unsupported datatype '" << typeid( value_t ).name() );
      }
   }

   /// Return a vector with names of all functions contained in the store
   std::vector< std::string > getFunctionNames() const
   {
      std::vector< std::string > names;
      names.reserve( size() );
      for ( const auto& func : r64Funcs_ )
      {
         names.push_back( func.getFunctionName() );
      }
      for ( const auto& func : r32Funcs_ )
      {
         names.push_back( func.getFunctionName() );
      }
      for ( const auto& func : i64Funcs_ )
      {
         names.push_back( func.getFunctionName() );
      }
      for ( const auto& func : i32Funcs_ )
      {
         names.push_back( func.getFunctionName() );
      }
      return names;
   }

   /// remove function by name
   template < typename value_t >
   inline void remove( const func_t< value_t >& function )
   {
      std::string funcName = function.getFunctionName();
      uint_t      oldSize  = this->size();

      std::function< bool( const func_t< value_t >& ) > predicate = [&funcName]( const func_t< value_t >& func ) {
         return func.getFunctionName() == funcName;
      };

      if constexpr ( std::is_same_v< value_t, double > )
      {
         for ( const auto& entry : r64Funcs_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( " *** before: " << entry.getFunctionName() );
         }
         r64Funcs_.erase( std::remove_if( r64Funcs_.begin(), r64Funcs_.end(), predicate ), r64Funcs_.end() );
         for ( const auto& entry : r64Funcs_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( " *** after:  " << entry.getFunctionName() );
         }
      }

      else if constexpr ( std::is_same_v< value_t, float > )
      {
         r32Funcs_.erase( std::remove_if( r32Funcs_.begin(), r32Funcs_.end(), predicate ), r32Funcs_.end() );
      }
      else if constexpr ( std::is_same_v< value_t, int > )
      {
         i32Funcs_.erase( std::remove_if( i32Funcs_.begin(), i32Funcs_.end(), predicate ), i32Funcs_.end() );
      }
      else if constexpr ( std::is_same_v< value_t, long > )
      {
         i64Funcs_.erase( std::remove_if( i64Funcs_.begin(), i64Funcs_.end(), predicate ), i64Funcs_.end() );
      }
      else if constexpr ( std::is_same_v< value_t, long long > )
      {
         i64Funcs_.erase( std::remove_if( i64Funcs_.begin(), i64Funcs_.end(), predicate ), i64Funcs_.end() );
      }
      else
      {
         WALBERLA_ABORT( "FunctionMultiStore::remove() detected unsupported datatype '" << typeid( value_t ).name() );
      }

      if ( this->size() == oldSize )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "Could not remove function '" << function.getFunctionName()
                                                                     << "' as it is not inside the FunctionMultiStore object!" );
      }
   }

 private:
   std::vector< func_t< double > >  r64Funcs_;
   std::vector< func_t< float > >   r32Funcs_;
   std::vector< func_t< int32_t > > i32Funcs_;
   std::vector< func_t< int64_t > > i64Funcs_;
};

} // namespace hyteg
