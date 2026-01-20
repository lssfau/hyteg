/*
 * Copyright (c) 2020-2026 Daniel Drzisga, Marcus Mohr.
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

#include "hyteg/ccrfunctionspace/CCRStokesFunction.hpp"
#include "hyteg/ccrfunctionspace/P2PlusBubbleFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFProjectNormalOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1ProjectNormalOperator.hpp"

namespace hyteg {

using walberla::real_t;

class P2PlusBubbleProjectNormalOperator
: public Operator< P2PlusBubbleVectorFunction< real_t >, P2PlusBubbleVectorFunction< real_t > >
{
 public:
   P2PlusBubbleProjectNormalOperator( const std::shared_ptr< PrimitiveStorage >&               storage,
                                      size_t                                                   minLevel,
                                      size_t                                                   maxLevel,
                                      const std::function< void( const Point3D&, Point3D& ) >& normalFunction );

   ~P2PlusBubbleProjectNormalOperator() override = default;

   void project( const P2PlusBubbleVectorFunction< real_t >& dst, size_t level, DoFType flag ) const;

   void project( const CCRStokesFunction< real_t >& dst, size_t level, DoFType flag ) const;

   /// Assemble operator as sparse matrix
   ///
   /// \param mat     a sparse matrix proxy
   /// \param numSrc  P2PlusBubbleVectorFunction for determining row and column indices
   /// \param numDst  unsued
   /// \param level   level in mesh hierarchy for which local operator is to be assembled
   /// \param flag    determines on which primitives this operator is assembled
   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2PlusBubbleVectorFunction< idx_t >&  numSrc,
                  const P2PlusBubbleVectorFunction< idx_t >&  numDst,
                  uint_t                                      level,
                  DoFType                                     flag ) const override;

 private:
   P1ProjectNormalOperator      p1Operator_;
   EdgeDoFProjectNormalOperator edgeDoFOperator_;
};

} // namespace hyteg
