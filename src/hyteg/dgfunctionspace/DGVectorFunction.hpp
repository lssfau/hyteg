/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/functions/CSFVectorFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/functions/VectorFunctionTools.hpp"

namespace hyteg {
namespace dg {

template < typename ValueType >
class DGVectorFunction final : public CSFVectorFunction< DGVectorFunction< ValueType > >
{
 public:
   using valueType = ValueType;

   template < typename VType >
   using FunctionType = DGVectorFunction< VType >;

   using VectorComponentType = DGFunction< ValueType >;

   using Tag = typename FunctionTrait< DGVectorFunction< ValueType > >::Tag;

   DGVectorFunction( const std::string&                         _name,
                     const std::shared_ptr< PrimitiveStorage >& storage,
                     size_t                                     minLevel,
                     size_t                                     maxLevel,
                     const std::shared_ptr< DGBasisInfo >&      basis,
                     uint_t                                     initialPolyDegree )
   : DGVectorFunction( _name, storage, minLevel, maxLevel, basis, initialPolyDegree, BoundaryCondition::create0123BC() )
   {}

   DGVectorFunction( const std::string&                         _name,
                     const std::shared_ptr< PrimitiveStorage >& storage,
                     size_t                                     minLevel,
                     size_t                                     maxLevel,
                     const std::shared_ptr< DGBasisInfo >&      basis,
                     uint_t                                     initialPolyDegree,
                     BoundaryCondition                          bc )
   : CSFVectorFunction< DGVectorFunction< ValueType > >( _name )
   {
      this->compFunc_.clear();
      this->compFunc_.push_back(
          std::make_shared< VectorComponentType >( _name + "_u", storage, minLevel, maxLevel, basis, initialPolyDegree, bc ) );
      this->compFunc_.push_back(
          std::make_shared< VectorComponentType >( _name + "_v", storage, minLevel, maxLevel, basis, initialPolyDegree, bc ) );

      if ( storage->hasGlobalCells() )
      {
         this->compFunc_.push_back(
             std::make_shared< VectorComponentType >( _name + "_w", storage, minLevel, maxLevel, basis, initialPolyDegree, bc ) );
      }
   }

   DGVectorFunction( const std::string name, const std::vector< std::shared_ptr< DGFunction< ValueType > > >& compFunc )
   : CSFVectorFunction< DGVectorFunction< ValueType > >( name, compFunc )
   {}
};

} // namespace dg
} // namespace hyteg
