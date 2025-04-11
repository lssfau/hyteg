/*
 * Copyright (c) 2017-2025 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/p2functionspace/P2PlusBubbleVectorFunction.hpp"

namespace hyteg {

/// Composite function for the Conforming Crouzeix-Raviart element
template < typename ValueType >
class CCRStokesFunction : public BlockFunction< ValueType >
{
 public:
   using valueType = ValueType;

   template < typename VType >
   using FunctionType = CCRStokesFunction< VType >;

   using VelocityFunction_T = P2PlusBubbleVectorFunction< ValueType >;
   using PressureFunction_T = DG1Function< ValueType >;

   using Tag = typename FunctionTrait< CCRStokesFunction< ValueType > >::Tag;

   CCRStokesFunction( const std::string&                         _name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      size_t                                     minLevel,
                      size_t                                     maxLevel )
   : CCRStokesFunction( _name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
   {}

   CCRStokesFunction( const std::string&                         _name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      size_t                                     minLevel,
                      size_t                                     maxLevel,
                      BoundaryCondition                          velocityBC )
   : BlockFunction< ValueType >( _name )
   {
      this->subFunc_.push_back( std::make_shared< FunctionWrapper< P2PlusBubbleVectorFunction< ValueType > > >(
          _name + "_uvw", storage, minLevel, maxLevel, velocityBC ) );
      this->subFunc_.push_back( std::make_shared< FunctionWrapper< DG1Function< ValueType > > >(
          _name + "_p", storage, minLevel, maxLevel, BoundaryCondition::createAllInnerBC() ) );
   }

   [[nodiscard]] const P2PlusBubbleVectorFunction< ValueType >& uvw() const
   {
      return this->subFunc_[0]->template unwrap< P2PlusBubbleVectorFunction< ValueType > >();
   }

   [[nodiscard]] P2PlusBubbleVectorFunction< ValueType >& uvw()
   {
      return this->subFunc_[0]->template unwrap< P2PlusBubbleVectorFunction< ValueType > >();
   }

   [[nodiscard]] const DG1Function< ValueType >& p() const
   {
      return this->subFunc_[1]->template unwrap< DG1Function< ValueType > >();
   }

   [[nodiscard]] DG1Function< ValueType >& p() { return this->subFunc_[1]->template unwrap< DG1Function< ValueType > >(); }
};

} // namespace hyteg
