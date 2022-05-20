/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

using walberla::real_t;

/// dummy operator used for cases where scalar operators are only available for 2D settings
class P2ToP1DummyOperator : public Operator< P2Function< real_t >, P1Function< real_t > >
{
 public:
   P2ToP1DummyOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   {}

   void apply( const P2Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      WALBERLA_ABORT( "P2ToP1DummyOperator::apply() should never be called!" );
   }
};

class P1ToP2DummyOperator : public Operator< P1Function< real_t >, P2Function< real_t > >
{
 public:
   P1ToP2DummyOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   {}

   void apply( const P1Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      WALBERLA_ABORT( "P1ToP2DummyOperator::apply() should never be called!" );
   }
};

} // namespace hyteg
