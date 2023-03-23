/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"

namespace hyteg {

template < typename ValueType >
class EGP0StokesFunction : public BlockFunction< ValueType >
{
 public:
   using valueType = ValueType;

   template < typename VType >
   using FunctionType = EGP0StokesFunction< VType >;

   using VelocityFunction_T = EGFunction< ValueType >;
   using PressureFunction_T = P0Function< ValueType >;

   using Tag = typename FunctionTrait< EGP0StokesFunction< ValueType > >::Tag;

   EGP0StokesFunction( const std::string&                         _name,
                          const std::shared_ptr< PrimitiveStorage >& storage,
                          size_t                                     minLevel,
                          size_t                                     maxLevel )
   : EGP0StokesFunction( _name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
   {}

   EGP0StokesFunction( const std::string&                         _name,
                          const std::shared_ptr< PrimitiveStorage >& storage,
                          size_t                                     minLevel,
                          size_t                                     maxLevel,
                          BoundaryCondition                          velocityBC )
   : BlockFunction< ValueType >( _name )
   {
      this->subFunc_.push_back( std::make_shared< FunctionWrapper< EGFunction< ValueType > > >(
          _name + "_uvw", storage, minLevel, maxLevel, velocityBC ) );
      this->subFunc_.push_back( std::make_shared< FunctionWrapper< P0Function< ValueType > > >(
          _name + "_p", storage, minLevel, maxLevel, BoundaryCondition::createAllInnerBC() ) );
   }

   [[nodiscard]] const EGFunction< ValueType >& uvw() const
   {
      return this->subFunc_[0]->template unwrap< EGFunction< ValueType > >();
   }

   [[nodiscard]] EGFunction< ValueType >& uvw() { return this->subFunc_[0]->template unwrap< EGFunction< ValueType > >(); }

   [[nodiscard]] const P0Function< ValueType >& p() const
   {
      return this->subFunc_[1]->template unwrap< P0Function< ValueType > >();
   }

   [[nodiscard]] P0Function< ValueType >& p() { return this->subFunc_[1]->template unwrap< P0Function< ValueType > >(); }

   /// \todo Get rid of this
   bool isDummy() const { return false; }
};

void applyDirichletBC( const EGP0StokesFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level );

} // namespace hyteg
