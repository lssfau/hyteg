/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

template < typename ValueType >
class P2P2StokesFunction : public BlockFunction< ValueType >
{
 public:
   using valueType = ValueType;

   template < typename VType >
   using FunctionType = P2P2StokesFunction< VType >;

   using VelocityFunction_T = P2VectorFunction< ValueType >;
   using PressureFunction_T = P2Function< ValueType >;

   using Tag = typename FunctionTrait< P2P2StokesFunction< ValueType > >::Tag;

   P2P2StokesFunction( const std::string&                         _name,
                       const std::shared_ptr< PrimitiveStorage >& storage,
                       size_t                                     minLevel,
                       size_t                                     maxLevel )
   : P2P2StokesFunction( _name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
   {}

   P2P2StokesFunction( const std::string&                         _name,
                       const std::shared_ptr< PrimitiveStorage >& storage,
                       size_t                                     minLevel,
                       size_t                                     maxLevel,
                       BoundaryCondition                          velocityBC )
   : BlockFunction< ValueType >( _name )
   {
      this->subFunc_.push_back( std::make_shared< FunctionWrapper< P2VectorFunction< ValueType > > >(
          _name + "_uvw", storage, minLevel, maxLevel, velocityBC ) );
      this->subFunc_.push_back( std::make_shared< FunctionWrapper< P2Function< ValueType > > >(
          _name + "_p", storage, minLevel, maxLevel, BoundaryCondition::createAllInnerBC() ) );
   }

   [[nodiscard]] const P2VectorFunction< ValueType >& uvw() const {
      return this->subFunc_[0]->template unwrap< P2VectorFunction< ValueType > >();
   }

   [[nodiscard]] P2VectorFunction< ValueType >& uvw()
   {
      return this->subFunc_[0]->template unwrap< P2VectorFunction< ValueType > >();
   }

   [[nodiscard]] const P2Function< ValueType >& p() const {
      return this->subFunc_[1]->template unwrap< P2Function< ValueType > >();
   }

   [[nodiscard]] P2Function< ValueType >& p()
   {
      return this->subFunc_[1]->template unwrap< P2Function< ValueType > >();
   }
};

} // namespace hyteg
