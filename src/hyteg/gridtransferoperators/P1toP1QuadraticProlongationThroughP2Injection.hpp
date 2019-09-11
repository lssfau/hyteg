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

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/FunctionMemory.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"

namespace hyteg {

class P1toP1QuadraticProlongationThroughP2Injection : public ProlongationOperator< P1Function< real_t > >
{
public:
  P1toP1QuadraticProlongationThroughP2Injection( const std::shared_ptr< PrimitiveStorage >& storage,
                                                 const uint_t&                              p1MinLevel,
                                                 const uint_t&                              p1MaxLevel )
  : tmp_( "tmp", storage, p1MinLevel - 1, p1MaxLevel - 1)
  {
    WALBERLA_CHECK_GREATER( p1MinLevel, 2 );
  }

  inline void prolongate( const P1Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
  {
    tmp_.assign( function, sourceLevel - 1, All );
    quadraticProlongation_.prolongate( tmp_, sourceLevel - 1, flag );
    function.assign( tmp_, sourceLevel + 1, All );
  }

private:

    P2Function< real_t > tmp_;
    P2toP2QuadraticProlongation quadraticProlongation_;
};

}