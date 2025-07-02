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

#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

template < typename value_t >
class P2P1TaylorHoodBlockFunction : public BlockFunction< value_t >
{
 public:
   using valueType = value_t;

   template < typename VType >
   using FunctionType = P2P1TaylorHoodBlockFunction< VType >;

   P2P1TaylorHoodBlockFunction( const std::string&                         name,
                                const std::shared_ptr< PrimitiveStorage >& storage,
                                size_t                                     minLevel,
                                size_t                                     maxLevel )
   : BlockFunction< value_t >( name )
   {
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P2VectorFunction< value_t > > >( name + "_uvw", storage, minLevel, maxLevel ) );
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1Function< value_t > > >( name + "_p", storage, minLevel, maxLevel ) );
   };
};

} // namespace hyteg
