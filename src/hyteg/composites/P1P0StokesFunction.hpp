/*
* Copyright (c) 2022 Nils Kohl.
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

#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

namespace hyteg {

template < typename ValueType >
class P1P0StokesFunction : public BlockFunction< ValueType >
{
 public:
   using valueType = ValueType;

   template < typename VType >
   using FunctionType = P1P0StokesFunction< VType >;
   
   P1P0StokesFunction( const std::string&                         name,
                       const std::shared_ptr< PrimitiveStorage >& storage,
                       size_t                                     minLevel,
                       size_t                                     maxLevel )
   : BlockFunction< ValueType >( name )
   {
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1VectorFunction< ValueType > > >( name + "_uvw", storage, minLevel, maxLevel ) );
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P0Function< ValueType > > >( name + "_p", storage, minLevel, maxLevel, BoundaryCondition::create0123BC() ) );
   };

   [[nodiscard]] const P1VectorFunction< ValueType >& uvw() const
   {
      return this->subFunc_[0]->template unwrap< P1VectorFunction< ValueType > >();
   }

   [[nodiscard]] P1VectorFunction< ValueType >& uvw()
   {
      return this->subFunc_[0]->template unwrap< P1VectorFunction< ValueType > >();
   }

   [[nodiscard]] const P0Function< ValueType >& p() const
   {
      return this->subFunc_[1]->template unwrap< P0Function< ValueType > >();
   }

   [[nodiscard]] P0Function< ValueType >& p() { return this->subFunc_[1]->template unwrap< P0Function< ValueType > >(); }
};

} // namespace hyteg
