/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl, Andreas Wagner.
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

#include "hyteg/functions/CSFVectorFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/functions/VectorFunctionTools.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

namespace hyteg {

template < typename ValueType >
class P2VectorFunction final : public CSFVectorFunction< P2VectorFunction< ValueType > >
{
 public:
   using valueType = ValueType;

   template < typename VType >
   using FunctionType = P2VectorFunction< VType >;

   using VectorComponentType = P2Function< ValueType >;

   using Tag = typename FunctionTrait< P2VectorFunction< ValueType > >::Tag;

   P2VectorFunction( const std::string&                         _name,
                     const std::shared_ptr< PrimitiveStorage >& storage,
                     size_t                                     minLevel,
                     size_t                                     maxLevel )
   : P2VectorFunction( _name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
   {}

   P2VectorFunction( const std::string&                         _name,
                     const std::shared_ptr< PrimitiveStorage >& storage,
                     size_t                                     minLevel,
                     size_t                                     maxLevel,
                     BoundaryCondition                          bc )
   : CSFVectorFunction< P2VectorFunction< ValueType > >( _name )
   {
      this->compFunc_.clear();
      this->compFunc_.push_back( std::make_shared< VectorComponentType >( _name + "_u", storage, minLevel, maxLevel, bc ) );
      this->compFunc_.push_back( std::make_shared< VectorComponentType >( _name + "_v", storage, minLevel, maxLevel, bc ) );

      if ( storage->hasGlobalCells() )
      {
         this->compFunc_.push_back( std::make_shared< VectorComponentType >( _name + "_w", storage, minLevel, maxLevel, bc ) );
      }
   }

   P2VectorFunction( const std::string name, const std::vector< std::shared_ptr< P2Function< ValueType > > >& compFunc )
   : CSFVectorFunction< P2VectorFunction< ValueType > >( name, compFunc )
   {}
};

} // namespace hyteg
