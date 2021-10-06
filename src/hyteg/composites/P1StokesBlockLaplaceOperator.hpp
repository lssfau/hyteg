/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

class P1StokesBlockLaplaceOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   P1StokesBlockLaplaceOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   {}

   void apply( const P1StokesFunction< real_t >& src,
               const P1StokesFunction< real_t >& dst,
               const size_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType ) const
   {
      for ( uint_t k = 0; k < src.uvw.getDimension(); k++ )
      {
         A.apply( src.uvw[k], dst.uvw[k], level, flag, updateType );
      }
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1StokesFunction< idx_t >&            src,
                  const P1StokesFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      A.toMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
      A.toMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );

      if ( src.uvw[0].getStorage()->hasGlobalCells() )
      {
         A.toMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );
      }
   }

   P1ConstantLaplaceOperator A;
};

} // namespace hyteg
