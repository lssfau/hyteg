/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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

#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

using walberla::real_t;

template < class operX_t, class operY_t, class operZ_t >
class P2ScalarToP2VectorOperator : public Operator< P2Function< real_t >, P2VectorFunction< real_t > >
{
 public:
   P2ScalarToP2VectorOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , operX( storage, minLevel, maxLevel )
   , operY( storage, minLevel, maxLevel )
   , operZ( storage, minLevel, maxLevel )
   {}

   void apply( const P2Function< real_t >&       src,
               const P2VectorFunction< real_t >& dst,
               size_t                            level,
               DoFType                           flag,
               UpdateType                        updateType = Replace ) const override final
   {
      operX.apply( src, dst[0], level, flag, updateType );
      operY.apply( src, dst[1], level, flag, updateType );
      if ( dst.getDimension() == 3 )
         operZ.apply( src, dst[2], level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2Function< idx_t >&                  src,
                  const P2VectorFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      operX.toMatrix( mat, src, dst[0], level, flag );
      operY.toMatrix( mat, src, dst[1], level, flag );
      if ( dst.getDimension() == 3 )
         operZ.toMatrix( mat, src, dst[2], level, flag );
   }

 private:
   operX_t operX;
   operY_t operY;
   operZ_t operZ;
};

// Some operators we might use more often than others
typedef P2ScalarToP2VectorOperator< P2ConstantDivTxOperator, P2ConstantDivTyOperator, P2ConstantDivTzOperator >
    P2ConstantDivTOperator;

} // namespace hyteg
